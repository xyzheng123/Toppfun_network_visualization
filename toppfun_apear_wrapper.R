
library(aPEAR)
library(ggplot2)
library(dplyr)


apear_find_clusters <- function(apear_input_list,
                                initial_minCS=3,
                                max_label_threshold = 50,
                                max_sign_pathway= 300,
                                min_leading_edge_threshold = 3,
                                output_dir = NULL,
                                save_cluster_parameters_and_results=TRUE, 
                                ...
                                
){
  
  
  ############### Check before find path calculation ########################### 
  if(save_cluster_parameters_and_results) {
    # Check if output_dir is NULL when save_cluster_parameters_and_results is TRUE
    if(is.null(output_dir)) {
      stop("save_cluster_parameters_and_results is TRUE but output_dir is not set. Please provide a valid output directory before proceeding.")
    } else{
      # Ensure the output directory exists or create it
      if (!dir.exists(output_dir)) {
        message(paste0("Output directory:",output_dir, " does NOT exist! Will create one!") )
        dir.create(output_dir, recursive = TRUE)
      }
    } 
  } else {
    # Handle case for save_cluster_parameters_and_results being FALSE
    # You can add logic here if there's anything specific to do when not saving results
    if(!is.null(output_dir)) {
      message("Note: Output directory is provided but save_cluster_parameters_and_results is FALSE, so no results will be saved.")
    }
    # Else, nothing needs to be saved, so you can optionally include logic for handling this scenario.
  }
  ############################################################################ 
  
  
  
  
  ########################Count Leading Edge genes########################
  #for each comparison in apear_input_list, count Leading Edge genes in core enrichment. number of "/" +1
  for (x in 1:length(apear_input_list)) {
    tmp_apear_input <- apear_input_list[[x]]
    for (y in 1:length(tmp_apear_input)) {
      tmp_str <- apear_input_list[[x]][[y]]$core_enrichment
      apear_input_list[[x]][[y]]$number_of_LE_genes <- sapply(tmp_str, function(s) sum(gregexpr("/", s)[[1]] >= 0)+1)
      apear_input_list[[x]][[y]] <- as.data.frame(apear_input_list[[x]][[y]])
    }
    
  }  
  
  ############################################################################
  
  
  
  
  
  
  ######################Create an output result list##########################
  
  #define unique_clusters
  unique_clusters <- lapply(apear_input_list, function(df) names(df))
  
  # create a list to store results for findPathCluster results of each cluster
  findPathClusterres <- vector("list",length = length(apear_input_list))
  names(findPathClusterres) <- names(apear_input_list)
  
  
  # loop through each comparison 
  for (x in 1:length(apear_input_list)) {
    # create a tmp cluster to work with in this iteration
    tmp_cluster <- unique_clusters[[names(apear_input_list)[x]]]
    
    # create a list to store result for each cluster
    findPathClusterres[[x]] <- vector("list",length = length(tmp_cluster))
    names(findPathClusterres[[x]]) <- tmp_cluster
    
    
    # loop through each cluster
    for (y in 1:length(tmp_cluster)) {
      
      findPathClusterres[[x]][[y]] <- vector("list",length = 2)
      names(findPathClusterres[[x]][[y]])<- c("findPathCluster_opt_inputs","findPathClusterres") #store both parameters and results
      
    }
  }
  ############################################################################
  
  
  
  
  
  
  ######################Calculate findPathCluster##############################
  # Loop through each comparison in apear_input_list
  for (x in 1:length(apear_input_list)) {
    tmp_apear_input <- apear_input_list[[x]]
    
    # Loop through each tmp_apear_input and pull out 1 cluster at a time
    for (y in 1:length(tmp_apear_input)) {
      tmp_cluster <- tmp_apear_input[[y]]
      
      # Filter out pathways with too few leading edge genes that cause troubles for markov clustering
      tmp_cluster <- tmp_cluster[which(tmp_cluster$number_of_LE_genes>min_leading_edge_threshold),]
      
      #filter number of pathways with too many pvals
      if(nrow(tmp_cluster) > max_sign_pathway) {
        tmp_cluster <- tmp_cluster %>% arrange(pval)
        tmp_cluster <- tmp_cluster %>% top_n(., max_sign_pathway, -pval)
      }
      
      msg_name <- paste0(names(apear_input_list)[x], "(", names(tmp_apear_input)[y], ")")
      
      
      # Function to create enrichment network; w/ dynamic minCS cutoff
      cal_minCS <- function(clustMethod) {
        tmp_minCS <- initial_minCS
        while (TRUE) {
          apear_clusters <- aPEAR::findPathClusters(tmp_cluster,cluster = clustMethod,minClusterSize = tmp_minCS) 
          num_clusters <- length(unique(apear_clusters$clusters$Cluster))
          if (num_clusters < max_label_threshold) {
            print(paste0(msg_name, ": ",clustMethod," clustering Done! ","minClusterSize=",tmp_minCS))
            return(tmp_minCS) 
          } else {
            tmp_minCS <- tmp_minCS + 1
          }
        }
      }
      
      # In the original dataframe indicate whether a pathway is used in the findPathClusters; In other words, it passes the min_leading_edge_threshold and max_sign_pathway
      tmp_df <- apear_input_list[[x]][[y]]
      Select_for_aPEAR_vec <- rep("No", times = nrow(tmp_df))
      Select_for_aPEAR_vec[tmp_df$Description %in% tmp_cluster$Description] <- "Yes"
      apear_input_list[[x]][[y]]$Select_for_aPEAR <- Select_for_aPEAR_vec
      
      
      # Attempt 1: Markov clustering
      result_markov <- tryCatch({
        cal_minCS("markov")
      }, error = function(e1) {
        cat(msg_name, ":", conditionMessage(e1), "Will try hierarchical clustering!\n")
        NULL  # Return NULL to trigger the second attempt
      })
      
      # Attempt 2: Hierarchical clustering
      if (is.null(result_markov)) {
        result_hier <- tryCatch({
          cal_minCS("hier")
        }, error = function(e2) {
          cat(msg_name, ":", conditionMessage(e2), msg_name, "cannot be processed via aPEAR!\n")
          NULL  # Return NULL if hierarchical clustering also fails
        })
        
        if (!is.null(result_hier)) {
          # store both parameters and results 
          findPathClusterres[[x]][[y]][["findPathCluster_opt_inputs"]] <- list("hier",result_hier) 
          names(findPathClusterres[[x]][[y]][["findPathCluster_opt_inputs"]]) <- c("clustMethod","minCS")
          findPathClusterres[[x]][[y]][["findPathClusterres"]] <- aPEAR::findPathClusters(tmp_cluster,cluster = "hier",minClusterSize = result_hier)
        }
      } else {
        # store both parameters and results
        findPathClusterres[[x]][[y]][["findPathCluster_opt_inputs"]] <- list("markov",result_markov)
        names(findPathClusterres[[x]][[y]][["findPathCluster_opt_inputs"]]) <- c("clustMethod","minCS")
        findPathClusterres[[x]][[y]][["findPathClusterres"]] <- aPEAR::findPathClusters(tmp_cluster,cluster = "markov",minClusterSize = result_markov)
      }
      
    }
  }
  ############################################################################
  
  
  ######################Packing findPathCluster and parameters##################
  # Package the results and the used parameters
  findPathClusterres <- list(
    clusterres = findPathClusterres,  
    params = list(
      max_sign_pathway = max_sign_pathway,
      min_leading_edge_threshold = min_leading_edge_threshold
    )
  )
  ############################################################################
  
  
  ######################Save findPathCluster##################################
  #save findPathCluster results and parameters and results to avoid recalculation
  if(save_cluster_parameters_and_results){
    if(!is.null(output_dir)){
      saveRDS(findPathClusterres,paste0(output_dir,"findPathClusterres.rds"))
    }
  }
  ############################################################################
  
  return(findPathClusterres)
}


apear_clean_failed_clusters <- function(findPathClusterres, apear_input_list) {
  
  clusterres <- findPathClusterres[["clusterres"]]
  
  # Normalize all names to avoid mismatch (spaces, tabs, unicode spaces)
  .norm <- function(x) {
    x <- trimws(x)                # remove leading/trailing spaces
    x <- gsub("\\s+", " ", x)     # collapse multiple spaces
    return(x)
  }
  
  for (goi in names(clusterres)) {
    names(clusterres[[goi]]) <- .norm(names(clusterres[[goi]]))
    names(apear_input_list[[goi]]) <- .norm(names(apear_input_list[[goi]]))
  }
  
  # Now perform cleaning
  for (goi in names(clusterres)) {
    for (cat in names(clusterres[[goi]])) {
      
      subitem <- clusterres[[goi]][[cat]]
      
      cond1 <- is.null(subitem$findPathCluster_opt_inputs)
      cond2 <- is.null(subitem$findPathClusterres)
      
      if (cond1 || cond2) {
        
        reason <- c()
        if (cond1) reason <- c(reason, "findPathCluster_opt_inputs == NULL")
        if (cond2) reason <- c(reason, "findPathClusterres == NULL")
        
        warning(sprintf(
          "Removing category '%s' from GOI '%s' due to: %s",
          cat, goi, paste(reason, collapse=" & ")
        ))
        
        # remove from findPathClusterres
        clusterres[[goi]][[cat]] <- NULL
        
        # remove from apear_input_list (only if present)
        if (cat %in% names(apear_input_list[[goi]])) {
          apear_input_list[[goi]][[cat]] <- NULL
        }
      }
    }
  }
  
  # Return cleaned objects
  findPathClusterres[["clusterres"]] <- clusterres
  
  return(list(
    findPathClusterres = findPathClusterres,
    apear_input_list = apear_input_list
  ))
}


apear_update_and_save_input <- function(apear_input_list,
                                        findPathClusterres,
                                        output_dir = NULL,
                                        save_apear_input = TRUE,
                                        ...){
  
  ############### Check before find path calculation ########################### 
  if(save_apear_input) {
    # Check if output_dir is NULL when save_apear_input is TRUE
    if(is.null(output_dir)) {
      stop("save_apear_input is TRUE but output_dir is not set. Please provide a valid output directory before proceeding.")
    }else{
      # Ensure the output directory exists or create it
      if (!dir.exists(output_dir)) {
        message(paste0("Output directory:",output_dir, " does NOT exist! Will create one!") )
        dir.create(output_dir, recursive = TRUE)
      }
    }  
  } else {
    # Handle case for save_cluster_parameters_and_results being FALSE
    # You can add logic here if there's anything specific to do when not saving results
    if(!is.null(output_dir)) {
      message("Note: Output directory is provided but save_cluster_parameters_and_results is FALSE, so no results will be saved.")
    }
    # Else, nothing needs to be saved, so you can optionally include logic for handling this scenario.
  }
  ############################################################################ 
  
  
  
  
  ########################Define and create output dir########################
  #for each comparison in apear_input_list, if the dir does not exist, create one
  if(save_apear_input) {
    for(x in 1:length(apear_input_list)){
      tmp_dir <- paste0(output_dir,names(apear_input_list)[x],"/")
      if(!dir.exists(paste0(tmp_dir,"inputs/"))){
        dir.create(paste0(tmp_dir,"inputs/"),recursive = T)
      }
    }
  }
  ############################################################################
  
  
  ######################Unpacking findPathCluster and parameters##################
  # Extract parameters from the findPathClusterres object
  max_sign_pathway <- findPathClusterres$params$max_sign_pathway
  min_leading_edge_threshold <- findPathClusterres$params$min_leading_edge_threshold
  
  findPathClusterres <- findPathClusterres$clusterres
  ############################################################################
  
  
  
  
  ########################Count Leading Edge genes########################
  #for each comparison in apear_input_list, count Leading Edge genes in core enrichment. number of "/" +1
  for (x in 1:length(apear_input_list)) {
    tmp_apear_input <- apear_input_list[[x]]
    for (y in 1:length(tmp_apear_input)) {
      tmp_str <- apear_input_list[[x]][[y]]$core_enrichment
      apear_input_list[[x]][[y]]$number_of_LE_genes <- sapply(tmp_str, function(s) sum(gregexpr("/", s)[[1]] >= 0)+1)
      apear_input_list[[x]][[y]] <- as.data.frame(apear_input_list[[x]][[y]])
    }
    
  }  
  
  ############################################################################
  
  
  ######################Find selected pathway##############################
  # Loop through each comparison in apear_input_list
  for (x in 1:length(apear_input_list)) {
    tmp_apear_input <- apear_input_list[[x]]
    
    # Loop through each tmp_apear_input and pull out 1 cluster at a time
    for (y in 1:length(tmp_apear_input)) {
      tmp_cluster <- tmp_apear_input[[y]]
      
      # Filter out pathways with too few leading edge genes that cause troubles for markov clustering
      tmp_cluster <- tmp_cluster[which(tmp_cluster$number_of_LE_genes>min_leading_edge_threshold),]
      
      #filter number of pathways with too many pvals
      if(nrow(tmp_cluster) > max_sign_pathway) {
        tmp_cluster <- tmp_cluster %>% arrange(pval)
        tmp_cluster <- tmp_cluster %>% top_n(., max_sign_pathway, -pval)
      }
      
      # In the original dataframe indicate whether a pathway is used in the findPathClusters; In other words, it passes the min_leading_edge_threshold and max_sign_pathway
      tmp_df <- apear_input_list[[x]][[y]]
      Select_for_aPEAR_vec <- rep("No", times = nrow(tmp_df))
      Select_for_aPEAR_vec[tmp_df$Description %in% tmp_cluster$Description] <- "Yes"
      apear_input_list[[x]][[y]]$Select_for_aPEAR <- Select_for_aPEAR_vec
      
    }
  }
  
  ############################################################################
  
  
  
  
  #######################Update apear_input obj #############################
  # Loop through each comparison in apear_input_list
  for (x in 1:length(apear_input_list)) {
    tmp_apear_input <- apear_input_list[[x]]
    
    # Loop through each tmp_apear_input and pull out 1 cluster at a time
    for (y in 1:length(tmp_apear_input)) {
      tmp_cluster <- tmp_apear_input[[y]]
      tmp_cluster$Pathway <- tmp_cluster$Description
      tmp_clusterres <- as.data.frame(findPathClusterres[[x]][[y]][["findPathClusterres"]][["clusters"]])
      
      
      
      # combine findPathClusterres and apear_input dataframes based on their pathway to pass the cluster info
      merged_df <- left_join(tmp_cluster, tmp_clusterres, by = "Pathway")
      merged_df$Pathway <- NULL
      merged_df$Cluster[is.na(merged_df$Cluster)] <- c("NA")
      
      #store the final cluster results 
      apear_input_list[[x]][[y]]$aPEAR_Cluster <- merged_df$Cluster
    }
  }
  ############################################################################
  
  
  
  
  #######################Save the outputs#####################################
  if(save_apear_input) {
    if(!is.null(output_dir)){
      # Loop through each comparison in apear_input_list
      for (x in 1:length(apear_input_list)) {
        tmp_output_dir <- paste0(output_dir, names(apear_input_list)[x], "/")
        tmp_apear_input <- apear_input_list[[x]]
        # Loop through each tmp_apear_input and pull out 1 cluster at a time
        for (y in 1:length(tmp_apear_input)) {
          
          apear_input_tmp <-  apear_input_list[[x]][[y]]
          
          # Remove NES column if present
          if ("NES" %in% colnames(apear_input_tmp)) {
            apear_input_tmp$NES <- NULL
          }
          
          
          #save the result inputs as csv for future usage 
          write.csv(apear_input_tmp, paste0(tmp_output_dir,"inputs/",names(apear_input_list)[x],"_",names(apear_input_list[[x]])[y],"_aPEAR_input.csv"),row.names = F,quote = F)
        }
      }
    }else{
      message("Please provide an output directory for the output_dir paramater")
    }  
  } 
  
  ############################################################################
  
  
  return(apear_input_list)
}



apear_plot_cluster_mod <- function(apear_input_list,
                                   findPathClusterres,
                                   repelLabels = FALSE,
                                   drawEllipses = FALSE,
                                   fontSize = 2.5,
                                   colorBy_pval = FALSE,
                                   default_nodecolor= "firebrick",
                                   show_parameters=TRUE,
                                   output_dir = NULL,
                                   output_width = 10,
                                   output_height = 10,
                                   save_network_plot=TRUE, 
                                   save_apearres_obj=TRUE, 
                                   ...){
  
  
  ############### Check before plot cluster ########################### 
  # Check for `output_dir` only if needed to save
  if((save_network_plot || save_apearres_obj) && is.null(output_dir)) {
    stop("Please provide an output directory for the output_dir parameter")
  }else{
    if(!is.null(output_dir)) {
      if(!(save_network_plot || save_apearres_obj)){
        message("Note: Output directory is provided but save_network_plot and/or save_apearres_obj is FALSE, some results will be saved.")
      }
    }else{
      # Ensure the output directory exists or create it
      if (!dir.exists(output_dir)) {
        message(paste0("Output directory:",output_dir, " does NOT exist! Will create one!") )
        dir.create(output_dir, recursive = TRUE)
      }
    } 
  }
  
  
  #########################################################################################
  
  
  ######################Unpacking findPathCluster and parameters##################
  # Extract parameters from the findPathClusterres object
  max_sign_pathway <- findPathClusterres$params$max_sign_pathway
  min_leading_edge_threshold <- findPathClusterres$params$min_leading_edge_threshold
  
  findPathClusterres <- findPathClusterres$clusterres
  ############################################################################
  
  
  
  ########################Create apearres##############################################
  #define unique_clusters
  unique_clusters <- lapply(apear_input_list, function(df) names(df))
  
  # create a list to store results for apear outputs of each cluster
  apearres <- vector("list",length = length(apear_input_list))
  names(apearres) <- names(apear_input_list)
  
  
  # loop through each comparison 
  for (x in 1:length(apear_input_list)) {
    # create a tmp cluster to work with in this iteration
    tmp_cluster <- unique_clusters[[names(apear_input_list)[x]]]
    
    # create a list to store result for each cluster
    apearres[[x]] <- vector("list",length = length(tmp_cluster))
    names(apearres[[x]]) <- tmp_cluster
    
    
    # loop through each cluster
    for (y in 1:length(tmp_cluster)) {
      apearres[[x]][[y]] <- vector("list",length = 2)
      names(apearres[[x]][[y]])<- c("apear_input","netp") #store both inputs and networkplot
      
      
      apearres[[x]][[y]][["apear_input"]] <- apear_input_list[[x]][[y]]
    }
  }
  #########################################################################################
  
  
  
  
  
  ########################Output and Plotting##############################################
  
  # Loop through each comparison in apear_input_list
  for (x in 1:length(apear_input_list)) {
    tmp_apear_input <- apear_input_list[[x]]
    tmp_output_dir <- paste0(output_dir, names(apear_input_list)[x], "/")
    
    if (!dir.exists(tmp_output_dir)) {
      dir.create(tmp_output_dir, recursive = TRUE)
    }
    
    # Loop through each tmp_apear_input and pull out 1 cluster at a time
    for (y in 1:length(tmp_apear_input)) {
      
      # Obtain the parameters and results from findPathCluster_obj
      tmp_cluster <- tmp_apear_input[[y]]  
      tmp_parameters <- findPathClusterres[[x]][[y]][["findPathCluster_opt_inputs"]]
      clustMethod <- tmp_parameters[[1]]
      ClusterSize <- tmp_parameters[[2]]
      
      
      tmp_res<- findPathClusterres[[x]][[y]][["findPathClusterres"]]
      
      msg_name <- paste0(names(apear_input_list), "(", names(tmp_apear_input)[y], ")")
      
      
      # whether to show the parameter in the final output
      if(show_parameters){
        parameters_string <- paste0(clustMethod," clustering; Minimum cluster size: ",ClusterSize, "\nMaxmium # of significant pathways: ", max_sign_pathway, "\nMinimum # of leading edge (LE) genes per pathway: ", min_leading_edge_threshold)
      }else{
        parameters_string <-""
      }
      
      if(colorBy_pval==TRUE){
        # Plot the network plot
        apearres_output <- aPEAR::plotPathClusters(tmp_cluster,
                                                   sim = tmp_res$similarity,
                                                   clusters = tmp_res$clusters,
                                                   fontSize = fontSize,
                                                   nodeSize = "number_of_LE_genes",
                                                   repelLabels = repelLabels, drawEllipses = drawEllipses,
                                                   colorBy = "pval")+
          ggtitle(paste0(msg_name),
                  subtitle = parameters_string)+
          labs(size = "# of LE genes")  + 
          theme(legend.position = "right") +
          scale_colour_gradient(low="firebrick", high="steelblue") + theme(legend.position = "right") }
      else{
        
        
        # Plot the network plot
        apearres_output <- aPEAR::plotPathClusters(tmp_cluster,
                                                   sim = tmp_res$similarity,
                                                   clusters = tmp_res$clusters,
                                                   fontSize = fontSize,
                                                   nodeSize = "number_of_LE_genes",
                                                   repelLabels = repelLabels, drawEllipses = drawEllipses)+
          ggtitle(paste0(msg_name),
                  subtitle = parameters_string)+
          labs(size = "# of LE genes")  + 
          theme(legend.position = "right") +
          scale_colour_gradient2(low = default_nodecolor,high = default_nodecolor, mid = default_nodecolor, midpoint = 0) + theme(legend.position = "right") + guides(color = guide_none())
      }
      
      
      if(save_network_plot){
        if(dir.exists(tmp_output_dir)){
          #save the network plot as pdfs 
          pdf(paste0(tmp_output_dir, names(tmp_apear_input)[y], "_aPEAR_", clustMethod, "_network.pdf"),width = output_width, height = output_height)
          print(apearres_output)
          dev.off()
        }
      }
      
      #update the network slots
      apearres[[x]][[y]][["netp"]] <- apearres_output 
      
    }
  }
  #########################################################################################
  
  
  
  
  
  ##################################save the apearres obj#######################################
  
  if(save_apearres_obj){
    if(!is.null(output_dir)){
      saveRDS(apearres,paste0(output_dir,"aPEARres.rds"))
    }else{
      message("Please provide an output directory for the output_dir paramater")
    }  
  }
  #########################################################################################
  
  
  
  return(apearres)
  
}



apear_move_saved_objects <- function(
    output_directory,
    result_of_interest,
    files = c("aPEARres.rds", "findPathClusterres.rds"),
    verbose = TRUE
) {
  
  # Normalize directory paths
  output_directory <- normalizePath(output_directory, mustWork = FALSE)
  
  # Define destination
  dest_dir <- file.path(output_directory, result_of_interest)
  
  # Create destination directory if needed
  if (!dir.exists(dest_dir)) {
    if (verbose) message("Creating destination folder: ", dest_dir)
    dir.create(dest_dir, recursive = TRUE)
  }
  
  # Full file paths
  src_paths <- file.path(output_directory, files)
  dest_paths <- file.path(dest_dir, files)
  
  # Check existence
  existing <- file.exists(src_paths)
  
  if (verbose) {
    for (i in seq_along(files)) {
      if (!existing[i]) {
        message("File not found, skipping: ", src_paths[i])
      }
    }
  }
  
  # Attempt to move only existing files
  moved <- rep(FALSE, length(files))
  moved[existing] <- file.rename(from = src_paths[existing], to = dest_paths[existing])
  
  # Messages
  if (verbose) {
    for (i in seq_along(files)) {
      if (moved[i]) {
        message("✅ Moved: ", files[i], " → ", dest_dir)
      } else if (existing[i]) {
        message("Failed to move: ", files[i])
      }
    }
  }
  
  # Return result
  return(invisible(moved))
}



run_toppfun_apear_pipeline <- function(
    result_of_interest,
    toppfun_result_directory,
    output_directory,
    excluded_toppfun_category = c("Cytoband","ToppCell Atlas",
                                  "MicroRNA","Drug","Coexpression Atlas","Computational"),
    colorBy_pval=FALSE,
    default_nodecolor = "firebrick",
    initial_minCS = 3,
    max_sign_pathway = 300,
    min_leading_edge_threshold = 3,
    fontSize = 2.5,
    output_width = 10,
    output_height = 10,
    ...
) {
  
  message("===== 1. Reading TOPPFUN input =====")
  
  ## 1. Read ToppFun Input
  fp <- file.path(toppfun_result_directory, paste0(result_of_interest, ".txt"))
  toppfun_res <- read.csv(fp, sep="\t", header=TRUE, quote = "", 
                          fill = TRUE, check.names = FALSE)
  
  if ("Pubmed" %in% toppfun_res$Category) {
    message("Removing 'Pubmed' rows — Pubmed enrichment results are messy and not suitable for aPEAR clustering.")
  }
  
  toppfun_res_noPubmed <- toppfun_res[!(toppfun_res$Category %in% "Pubmed"),]
  rm(toppfun_res)
  
  colnames(toppfun_res_noPubmed) <- gsub("\\.","_", colnames(toppfun_res_noPubmed))
  colnames(toppfun_res_noPubmed)[7]  <- "q_value_FDR_BH"
  colnames(toppfun_res_noPubmed)[8]  <- "q_value_FDR_BY"
  
  ## 2. Prepare apear_input
  apear_input <- toppfun_res_noPubmed
  colnames(apear_input)[3]  <- "Description"
  colnames(apear_input)[10] <- "setSize"
  colnames(apear_input)[5]  <- "pval"
  colnames(apear_input)[7]  <- "p.adjust"
  colnames(apear_input)[11] <- "core_enrichment"
  
  apear_input$ID <- NULL
  apear_input$Source <- NULL
  apear_input$q_value_Bonferroni <- NULL
  apear_input$q_value_FDR_BY <- NULL
  apear_input$Hit_Count_in_Query_List <- NULL
  
  apear_input$NES <- 1
  apear_input$core_enrichment <- gsub(",", "/", apear_input$core_enrichment)
  
  ## 3. Encode categories
  apear_input$Category <- gsub("GO: Molecular Function","GO_MF", apear_input$Category)
  apear_input$Category <- gsub("GO: Biological Process","GO_BP", apear_input$Category)
  apear_input$Category <- gsub("GO: Cellular Component","GO_CC", apear_input$Category)
  apear_input$Category <- gsub("Human Phenotype","Human_Phenotype", apear_input$Category)
  apear_input$Category <- gsub("Mouse Phenotype","Mouse_Phenotype", apear_input$Category)
  
  ## 4. Split into Category List
  uniq_cat <- unique(apear_input$Category)
  apear_input_list <- vector("list", length(uniq_cat))
  
  for (i in seq_along(uniq_cat)) {
    apear_input_list[[i]] <- apear_input[apear_input$Category == uniq_cat[i], ]
  }
  names(apear_input_list) <- uniq_cat
  
  ## 5. Deduplicate by minimum p-value
  message("===== 2. Deduplicating pathways =====")
  for (i in seq_along(apear_input_list)) {
    apear_input_list[[i]] <- apear_input_list[[i]] %>%
      dplyr::group_by(Description) %>%
      dplyr::slice_min(pval, with_ties = FALSE) %>%
      dplyr::ungroup()
  }
  
  ## 6. Remove unwanted categories
  message("===== 3. Filtering unwanted categories =====")
  apear_input_list_short <- apear_input_list[!(names(apear_input_list) %in% excluded_toppfun_category)]
  
  ## 7. Nest into a single GOI
  apear_input_list_nested <- list(apear_input_list_short)
  names(apear_input_list_nested) <- result_of_interest
  
  ## 8. Run aPEAR clustering
  message("===== 4. Running aPEAR clustering =====")
  findPathClusterres <- apear_find_clusters(
    apear_input_list_nested,
    initial_minCS = initial_minCS,
    max_sign_pathway = max_sign_pathway,
    min_leading_edge_threshold = min_leading_edge_threshold,
    output_dir = output_directory,
    save_cluster_parameters_and_results = TRUE
  )
  
  ## 9. Clean failed clusters
  message("===== 5. Cleaning failed clusters =====")
  cleaned <- apear_clean_failed_clusters(
    findPathClusterres,
    apear_input_list_nested
  )
  
  findPathClusterres <- cleaned$findPathClusterres
  apear_input_list_nested <- cleaned$apear_input_list
  
  ## 10. Update apear input + save
  message("===== 6. Updating aPEAR input =====")
  apear_input_list_updated <- apear_update_and_save_input(
    apear_input_list_nested,
    findPathClusterres,
    output_dir = output_directory,
    save_apear_input = TRUE
  )
  
  ## 11. Plot networks
  message("===== 7. Generating network plots =====")
  apear_network_list <- apear_plot_cluster_mod(
    apear_input_list_updated,
    findPathClusterres,
    fontSize = fontSize,
    colorBy_pval=colorBy_pval,
    default_nodecolor = default_nodecolor,
    output_dir = output_directory,
    output_width = output_width,
    output_height = output_height,
    save_network_plot = TRUE,
    save_apearres_obj = TRUE
    
  )
  
  ## 12. Move objects into subfolder
  message("===== 8. Moving saved objects =====")
  apear_move_saved_objects(
    output_directory = output_directory,
    result_of_interest = result_of_interest
  )
  
  message("===== PIPELINE FINISHED SUCCESSFULLY =====")
}



run_toppfun_apear_pipeline(
  result_of_interest = "Demo_ToppFun_Analysis",  # same name as the ToppFun txt file w/o .txt
  toppfun_result_directory = "example_data/toppfun_results/",
  output_directory = "example_output/aPEAR_clusters/",
  excluded_toppfun_category = c(
    "Cytoband", "ToppCell Atlas", "Disease", "Gene Family",
    "MicroRNA", "Drug", "Coexpression Atlas", "Computational"
  ),
  colorBy_pval = FALSE
)





