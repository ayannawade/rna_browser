function(input, output, session) { 
  
  ##########################################################################################################################
  # REACTIVE VALUES
  ##########################################################################################################################
  
  RV_data = reactiveValues(counts_matrix = data.frame(),
                           metadata = data.frame(),
                           selected_genes = data.frame(),
                           exp_min = 0,
                           exp_max = 100,
                           r_type = "l",
                           DGE_matrix = NULL,
                           all_mouse = T,
                           generated = F,
                           user_data_open = F)
  
  RV_UI = reactiveValues(gene_names = list(),
                         label_names = list(),
                         condition_names = list(),
                         fc_names = list(),
                         pc_choices = c("PC1", "PC2"))
  
  RV_pca_plot = reactiveValues(plot = ggplot(),
                               title = "PCA Scatter Plot",
                               clicks = NULL)
  
  RV_mds_plot = reactiveValues(plot = ggplot(),
                               title = "MDS Scatter Plot",
                               clicks = NULL)
  
  RV_box_plot = reactiveValues(plot = ggplot(),
                               title = "DGE Box Plot",
                               clicks = NULL)
  
  RV_fc_table = reactiveValues(selector = list(),
                               table_list = list(),
                               table = data.frame(),
                               display_table = data.frame(),
                               subset_table = data.frame(),
                               logfc_min_1 = 0,
                               logfc_max_1 = 10,
                               logfc_min_2 = 0,
                               logfc_max_2 = 10)
  
  RV_fc_plot = reactiveValues(plot = ggplot(),
                              title = "Fold Change",
                              clicks = NULL)
  
  RV_heatmap = reactiveValues(matrix = data.frame("x"=c(0,1),"y"=c(1,0)),
                              title = "Heatmap",
                              clicks = NULL)
  
  RV_dynamic = reactiveValues(val = 0)
  
  reset = function() {
    # Blank current data
    RV_data$counts_matrix = data.frame()
    RV_data$metadata = data.frame()
    RV_data$selected_genes = data.frame("Sample_ID"=character())
    RV_data$exp_min = 0
    RV_data$exp_max = 100
    RV_data$DGE_matrix = data.frame()
    RV_data$all_mouse = T
    RV_data$r_type = "l"
    
    RV_UI$gene_names = list()
    RV_UI$label_names = list()
    RV_UI$condition_names = list()
    RV_UI$fc_names = list()
    
    RV_pca_plot$plot = ggplot()
    RV_pca_plot$clicks = NULL
    
    RV_mds_plot$plot = ggplot()
    RV_mds_plot$clicks = NULL
    
    RV_box_plot$plot = ggplot()
    RV_box_plot$clicks = NULL
    
    RV_fc_plot$plot = ggplot()
    RV_fc_plot$clicks = NULL
    
    RV_heatmap$matrix = data.frame("x"=c(0,1),"y"=c(1,0))
    
    RV_fc_table$selector = list()
    RV_fc_table$table_list = list()
    RV_fc_table$table = data.frame()
    RV_fc_table$display_table = data.frame()
    RV_fc_table$subset_table = data.frame()
    RV_fc_table$logfc_min_1 = 0
    RV_fc_table$logfc_max_1 = 10
    RV_fc_table$logfc_min_2 = 0
    RV_fc_table$logfc_max_2 = 10
  }
  
  ##########################################################################################################################
  # INITIAL DATA FORMATTING / READINS
  ##########################################################################################################################
  
  # Update server dataset(s)
  observeEvent(input$generate_button, {
    
    # Reset certain inputs the default
    reset()
    
    # Check if any data is selected (preloaded + user)
    dataset_count = F
    if (length(input$dataset) > 0 | !is.null(input$user_matrix$datapath)) {
      dataset_count = T
    }
    
    # If there is data selected ...
    if (dataset_count) {
      
      # For progress bar
      n = 7
      if (input$normalize == "r") {
        RV_data$r_type = "r"
        n = 9
      }
      
      withProgress(message = 'Initializing', value = 0, {
        
        # 1
        incProgress(1/n, message = "Collecting metadata")
        
        # 1.1 Get user sample info, if present
        if (!is.null(input$user_metadata$datapath)) {
          user_metadata = read.table(input$user_metadata$datapath, header = T, sep = "\t", quote = "")
          annotations[[length(annotations)+1]] = user_metadata
        }
        
        # 1.2 Merge sample info tables
        merged_annotations = annotations[[1]]
        if (length(annotations) > 1) {
          for (i in 2:length(annotations)) {
            merged_annotations = merge(merged_annotations, annotations[[i]], all = T)
          }
        }
        merged_annotations[is.na(merged_annotations)] = "N/A"
        
        # 1.3 Set reactive values for metadata
        RV_UI$label_names = colnames(merged_annotations)
        RV_UI$condition_names = colnames(merged_annotations)
        
        # 2
        incProgress(1/n, message = "Checking user data")
        
        # 2.1 Get gene matrix
        working_matrix = data.frame("gene"=character())
        if (!is.null(input$user_matrix$datapath)) {
          user_matrix = read.table(input$user_matrix$datapath, header = T, sep = "\t", quote = "")
          
          # 2.1.1 Change user table gene names if not mouse
          if (input$user_organism == "Human") {
            RV_data$all_mouse = F
            working_matrix = convertGenes(user_matrix, "Human")
          } else {
            working_matrix = user_matrix
          }
        }
        working_matrix[is.na(working_matrix)] = 0
        colnames(working_matrix)[colnames(working_matrix) == "gene_id"] = "gene"
        
        # 2.2 Upload fold change tables
        fc_names = list()
        fc_frames = list()
        fc_index = 1
        
        if (!is.null(input$user_fold_change$datapath)) {
          for(i in 1:length(input$user_fold_change[, 1])) {
            fc_frames[[fc_index]] = read.table(input$user_fold_change[[i, 'datapath']], header = T, sep = "\t", quote = "")
            fc_names[[fc_index]] = input$user_fold_change$name
            fc_index = fc_index + 1
          }
        }
        
        # 3
        
        if (length(input$dataset) > 0) {
          print("PRELOADED DATA SELECTED")
          
          # 3.1 Check selected sets for human data
          selected_sets = input$dataset
          for (i in 1:length(selected_sets)) {
            organism = datasets$organism[datasets$Datasets == selected_sets[i]]
            if ("Human" == organism) {
              RV_data$all_mouse = F
            }
          }
          
          # 3.2 Loop: Import gene tables and fold change tables
          for (i in 1:length(selected_sets)) {
            x = 1/length(selected_sets)
            incProgress(x/n, message = paste("Loading Dataset: ", i, " of ", length(selected_sets), sep = ""))
            
            # 3.2.1 Import gene matrix
            dir = datasets$dir[datasets$Datasets == selected_sets[i]]
            organism = datasets$organism[datasets$Datasets == selected_sets[i]]
            path = NULL
            if (RV_data$r_type == "r") {
              path = paste("./data/", dir, "/GENE_MATRIX.txt", sep = "")
            }
            else {
              path = paste("./data/", dir, "/logRPKM_MATRIX.txt", sep = "")
            }
            print(path)
            data = read.table(path, sep = "\t", header = T, stringsAsFactors = F)
            colnames(data)[colnames(data) == "gene_id"] = "gene"
            
            # 3.2.2 Convert gene names if data isn't all mouse data
            if (RV_data$all_mouse == F) {
              data = convertGenes(data, organism)
            }
            
            # 3.2.3 Loop: Import all fold change tables for a dataset
            current_fc_names = list.files(paste("./data/", dir, "/FOLD_CHANGE", sep = ""))
            for (name in current_fc_names) {
              fc_names[[fc_index]] = name
              
              current_path = paste("./data/", dir, "/FOLD_CHANGE/", name, sep = "")
              current_frame = read.table(current_path, sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
              current_frame = cbind("gene"=rownames(current_frame), current_frame)
              row.names(current_frame) = seq(1:nrow(current_frame))
              
              fc_frames[[fc_index]] = current_frame
              fc_index = fc_index + 1
            }
            
            # 3.3 Merge current matrix with previously looped matrices
            merged_matrix = merge(working_matrix, data, by = "gene", all = T)
            merged_matrix[is.na(merged_matrix)] = 0
            working_matrix = merged_matrix
          }
        }
        
        # 3.4 Save fold change reactive values
        RV_fc_table$table_list = fc_frames
        RV_UI$fc_names = unlist(lapply(fc_names, function(i) unlist(strsplit(i, "[.]"))[1]))
        
        # 4
        incProgress(1/n, message = "Calculating MDS")
        
        # 4.1 Prep gene count matrix for MDS and DGE analysis
        working_matrix = working_matrix[!duplicated(working_matrix$gene), ]
        row.names(working_matrix) = working_matrix$gene
        working_matrix$gene = NULL
        transposed_matrix = data.frame(t(working_matrix))
        
        d = dist(transposed_matrix)
        fit = cmdscale(d, eig=TRUE, k=2)
        
        MDS_data = data.frame(fit$points)
        MDS_data$Sample_ID = row.names(MDS_data)
        MDS_data = merge(MDS_data, merged_annotations, by = "Sample_ID")
        
        # 4.2 Set reactive values
        RV_UI$gene_names = row.names(working_matrix)
        RV_data$counts_matrix = working_matrix
        RV_data$metadata = MDS_data
        
        # 4.3 Set plot names for saving
        title_list = unique(MDS_data$Dataset)
        title_name = ""
        for (i in 1:length(title_list)) {
          title_name = paste(title_name, title_list[i], ",", sep = "")
        }
        title_name = substr(title_name, 1, nchar(title_name) - 1)
        
        # Only run DGE on raw data, NOT RPKM (can only handle integers)
        if (RV_data$r_type == "r") {
          
          # 5
          incProgress(1/n, message = "Calculating DGE")
          
          # 5.1 Set parameters for DGE analysis
          annotations = RV_data$metadata
          matrix = RV_data$counts_matrix
          cpm.Sample_ID.cutoff = 2
          min.cpm = 1
          
          matrix = matrix[,which(!apply(matrix,2,FUN = function(x){all(x == 0)}))]
          
          # 5.2 Set up edgeR expression object
          DGE_data = DGEList(counts=matrix)
          
          keep = rowSums(cpm(DGE_data)>min.cpm) >= cpm.Sample_ID.cutoff
          DGE_data = DGE_data[keep, , keep.lib.sizes=FALSE]
          
          # 6
          
          # 6.1  If DGE is successfull, run filters
          if (is.na(DGE_data$counts[1,1]) | DGE_data$counts[1,1] < 0) {
            incProgress(1/n, message = "DGE Failed")
            RV_data$DGE_matrix = NULL
          }
          
          else {
            # 6.2 Perform simple exact test on genotype
            incProgress(1/(3*n), message = "Running DGE filter 1 of 3")
            DGE_data <- calcNormFactors(DGE_data)
            incProgress(1/(3*n), message = "Running DGE filter 2 of 3")
            DGE_data <- estimateCommonDisp(DGE_data)
            incProgress(1/(3*n), message = "Running DGE filter 3 of 3")
            DGE_data <- estimateTagwiseDisp(DGE_data)
            
            # 6.3 Extract pseudo count data
            pseudo_counts = DGE_data$pseudo.counts
            pseudo_counts[pseudo_counts < 0] = 0
            
            # 6.4 Save reactive values
            RV_data$DGE_matrix = pseudo_counts
          }
        }
        
        # 7
        incProgress(1/n, message = "Calculating PCA")
        
        # 7.1 Run PCA
        pca = prcomp(as.data.frame(t(RV_data$counts_matrix)), scale. = F, center = F)
        
        # 7.2 Format coordinates (only need )
        coords = as.data.frame(pca$x)
        coords = coords[, 1:10]
        RV_UI$pc_choices = colnames(coords)
        coords$Sample_ID = row.names(coords)
        
        # 7.3 Save PCA coordinates to reactive value
        RV_data$metadata = merge(RV_data$metadata, coords, by = "Sample_ID")
        
        # 8
        incProgress(1/n, message = "Creating heatmap")
        
        # 9
        incProgress(1/n, message = "Done")
      })
    }
  })
  
  # Heatmap calcultions could be done separately because it might take a while to
  observeEvent(input$heatmap_button, {
    if (nrow(RV_fc_table$subset_table) > 2) {
      gene_list = RV_fc_table$subset_table$gene
      full_matrix = RV_data$counts_matrix
      full_matrix$gene = rownames(full_matrix)
      
      subset_matrix = subset(full_matrix, gene %in% gene_list)
      subset_matrix$gene = NULL
      
      RV_heatmap$matrix = as.matrix(subset_matrix)
    }
  })
  
  ##########################################################################################################################
  # DYNAMIC DATA UPDATES
  ##########################################################################################################################
  
  # Activated by gene list entry: Update reactive values
  observeEvent(input$gene, ignoreNULL = F, {
    return_frame = data.frame("Sample_ID"=character())
    if (!is.null(input$gene)) {
      for (i in 1:length(input$gene)) {
        current_row = RV_data$counts_matrix[input$gene[[i]], ]
        transposed_gene = data.frame(t(current_row))
        transposed_gene$Sample_ID = row.names(transposed_gene)
        return_frame = merge(return_frame, transposed_gene, by = "Sample_ID", all = T)
      }
    }
    
    # Set coex threshold
    if (ncol(return_frame) == 2) {
      RV_data$exp_max = round(max(return_frame[2]), -1) + 20
      RV_data$exp_min = round(min(return_frame[2]), -1) + 20
    } else if (ncol(return_frame) == 3) {
      RV_data$exp_max = round(max(max(return_frame[2]), max(return_frame[3])), -1) + 20
      RV_data$exp_min = round(min(min(return_frame[2]), min(return_frame[3])), -1) - 20
    }
    
    RV_data$selected_genes = return_frame
  })
  
  # Activated for clicked points
  observeEvent(input$plot_click_1, {
    if (input$plot_tab == RV_pca_plot$title & !is.null(input$plot_click_1)) {
      RV_pca_plot$clicks = input$plot_click_1
    }
  })
  # Activated for clicked points
  observeEvent(input$plot_click_2, {
    if (input$plot_tab == RV_mds_plot$title & !is.null(input$plot_click_2)) {
      RV_mds_plot$clicks = input$plot_click_2
    }
  })
  # Activated for clicked points
  observeEvent(input$plot_click_3, {
    if (input$plot_tab == RV_box_plot$title & !is.null(input$plot_click_3)) {
      RV_box_plot$clicks = input$plot_click_3
    }
  })
  # Activated for clicked points
  observeEvent(input$plot_click_4, {
    if (input$plot_tab == RV_fc_plot$title & !is.null(input$plot_click_4)) {
      RV_fc_plot$clicks = input$plot_click_4
    }
  })
  
  # Select user data button toggle
  observeEvent(input$user_data_button, {
    if (input$user_data_button %% 2 == 0) {
      RV_data$user_data_open = F
    } else {
      RV_data$user_data_open = T
    }
  })
  
  # Select new data button toggle
  observeEvent(input$generate_button, {
    if (length(input$dataset) > 0 | !is.null(input$user_matrix$datapath)) {
      RV_data$generated = T
    }
  })
  observeEvent(input$select_new_data_button, {
    RV_data$generated = F
  })
  
  # Fold change table and plot
  observe({
    RV_fc_table$selector = input$fc_table_select
  })
  observeEvent(RV_fc_table$selector, ignoreNULL = F, {
    return_table = data.frame()
    if (length(RV_fc_table$selector) == 1) {
      index = match(RV_fc_table$selector[[1]], RV_UI$fc_names)
      return_table = RV_fc_table$table_list[[index]]
      return_table = return_table[order(-return_table$logFC), ]
    }
    else if (length(RV_fc_table$selector) == 2) {
      index_1 = match(RV_fc_table$selector[[1]], RV_UI$fc_names)
      index_2 = match(RV_fc_table$selector[[2]], RV_UI$fc_names)
      table_1 = RV_fc_table$table_list[[index_1]]
      table_2 = RV_fc_table$table_list[[index_2]]
      
      t1 = "Mouse"
      t2 = "Mouse"
      if (grepl("Human", RV_fc_table$selector[[1]])) {
        t1 = "Human"
      }
      
      if (grepl("Human", RV_fc_table$selector[[2]])) {
        t2 = "Human"
      }
      # If one set of human and one set of mouse data are present
      if ((t1 == "Mouse" & t2 == "Human") | (t1 == "Human" & t2 == "Mouse")) {
        
        master = data.frame("Master"=master_name_table$Master, "gene"=master_name_table$Mouse)
        
        mouse_data = NULL
        if (t1 == "Mouse") {
          mouse_data = table_1
          
          converted_data = merge(mouse_data, master, by = "gene")
          converted_data$gene = NULL
          colnames(converted_data)[colnames(converted_data) == "Master"] = "gene"
          
          table_1 = converted_data
        } 
        else if (t2 == "Mouse") {
          mouse_data = table_2
          
          converted_data = merge(mouse_data, master, by = "gene")
          converted_data$gene = NULL
          colnames(converted_data)[colnames(converted_data) == "Master"] = "gene"
          
          table_2 = converted_data
        }
      }
      
      colnames(table_1) = paste(colnames(table_1), "1", sep = "_")
      colnames(table_2) = paste(colnames(table_2), "2", sep = "_")
      colnames(table_1)[colnames(table_1) == "gene_1"] = "gene"
      colnames(table_2)[colnames(table_2) == "gene_2"] = "gene"
      
      return_table = merge(table_1, table_2, by = "gene")
      return_table = return_table[order(-return_table$logFC_1), ]
    }
    
    RV_fc_table$table = return_table
    
  })
  
  # Observe selected rows in fold change table
  observe({
    current_table = RV_fc_table$display_table
    return_table = current_table[input$fold_change_table_rows_selected, ]
    RV_fc_table$subset_table = return_table
  })
  
  ##########################################################################################################################
  # GENERATING PLOTS/TABLES
  ##########################################################################################################################
  
  # All possible reactives for generating plots are stored here
  observeEvent(c(RV_data$metadata,
                 RV_data$selected_genes,
                 RV_data$DGE_matrix,
                 RV_pca_plot$clicks,
                 RV_mds_plot$clicks,
                 RV_box_plot$clicks,
                 RV_fc_plot$clicks,
                 input$condition,
                 input$label_type,
                 input$pcs_chosen,
                 input$binary_button,
                 input$coex_threshold,
                 input$multi_label
  ), {
    RV_dynamic$val = RV_dynamic$val + 1
  })
  
  # Create/update pca plot
  observeEvent(RV_dynamic$val, {
    if (nrow(RV_data$metadata) > 0) {
      
      plot = NULL
      
      if (is.null(input$pcs_chosen)) {
        plot = ggplot() + xlab("PC.X") + ylab("PC.Y") + geom_text(aes(x=0, y=0, label="Two principal components must be selected in order to plot."), size = 8, colour = plot_text_color)
      } else if (length(input$pcs_chosen) < 2) {
        plot = ggplot() + xlab("PC.X") + ylab("PC.Y") + geom_text(aes(x=0, y=0, label="Two principal components must be selected in order to plot."), size = 8, colour = plot_text_color)
      } else {
        # Set parameters for graph
        plot_data = RV_data$metadata
        gene_data = RV_data$selected_genes
        P1_name = input$pcs_chosen[[1]]
        P2_name = input$pcs_chosen[[2]]
        
        # Adjust plot data to fit selected condition
        plot_data = findConditionCol(plot_data, input$condition, RV_data$selected_genes)
        # Adjust plot for showing gene expression
        plot = graphGenes(plot_data, gene_data, P1_name, P2_name, "PCA Plot", input$binary_button, input$coex_threshold)
        # Add clicked point
        plot = addLabel(plot, plot_data, RV_pca_plot$clicks, input$label_type, input$multi_label, 0.03, P1_name, P2_name)
      }
      
      RV_pca_plot$plot = plot
    }
  })
  
  # Create/update mds plot
  observeEvent(RV_dynamic$val, {
    if (nrow(RV_data$metadata) > 0) {
      
      # Set parameters for graph
      plot_data = RV_data$metadata
      gene_data = RV_data$selected_genes
      plot = NULL
      
      # Adjust plot data to fit selected condition
      plot_data = findConditionCol(plot_data, input$condition, RV_data$selected_genes)
      # Adjust plot for showing gene expression
      plot = graphGenes(plot_data, gene_data, "X1", "X2", "MDS Plot", input$binary_button, input$coex_threshold)
      # Add clicked point
      plot = addLabel(plot, plot_data, RV_mds_plot$clicks, input$label_type, input$multi_label, 0.03, "X1", "X2")
      
      RV_mds_plot$plot = plot
      
    }
  })
  
  # Create/update box plot
  observeEvent(RV_dynamic$val, {
    box_plot = NULL
    # Can only handle raw count expression
    if (RV_data$r_type != "r") {
      box_plot = ggplot() + xlab("X1") + ylab("X2") + geom_text(aes(x=0, y=0, label="DGE analysis can only be done using pseudo counts data."), size = 8, colour = plot_text_color)
    }
    else {
      pseudo_genes = RV_data$DGE_matrix
      
      # DGE might fail, in which case DGE object will remain null
      if (is.null(pseudo_genes)) {
        box_plot = ggplot() + xlab("X1") + ylab("X2") + geom_text(aes(x=0, y=0, label="DGE analysis failed. Likely insufficient counts in selected data."), size = 8, colour = plot_text_color)
      }
      else {
        # Gene must be selected in order to create box plot
        if (ncol(RV_data$selected_genes) <= 1) {
          box_plot = ggplot() + xlab("X1") + ylab("X2") + geom_text(aes(x=0, y=0, label="Enter a gene name to view it's expression here."), size = 8, colour = plot_text_color)
        }
        else {
          gene_name = colnames(RV_data$selected_genes[2])
          
          # Gene must exist in matrix (some are filtered out during DGE)
          if (!gene_name %in% row.names(pseudo_genes)) {
            box_plot = ggplot() + xlab("X1") + ylab("X2") + geom_text(aes(x=0, y=0, label="Sorry, the selected gene does not have sufficient counts to plot."), colour = plot_text_color, size = 8)
          }
          else {
            # All conditions passed
            current_data = data.frame("Sample_ID"=colnames(pseudo_genes), "count"=pseudo_genes[gene_name, ])
            plot_data = merge(current_data, RV_data$metadata, by = "Sample_ID")
            plot_data = findConditionCol(plot_data, input$condition, RV_data$selected_genes)
            
            box_plot = ggplot(plot_data, aes(x=factor(plot_cond_1), y = log(count), fill = factor(plot_cond_1))) +
              geom_boxplot(colour = plot_text_color, width = 0.5) +
              # geom_jitter(aes(fill = factor(plot_cond_1)), size = 4, shape = 21, colour = plot_text_color) +
              geom_point(aes(fill = factor(plot_cond_1)), size = dot_size, shape = 21, alpha = dot_alpha, colour = plot_text_color) +
              ggtitle(gene_name) +
              scale_fill_manual(values = cbPalette) +
              xlab("Group") +
              ylab("Log Count") +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              theme(legend.position="none")
            
            no_coords = plot_data
            no_coords$X1 = NULL
            no_coords$X2 = NULL
            df_coords = cbind(layer_data(box_plot, i=2), no_coords)
            colnames(df_coords)[colnames(df_coords) == "x" ] = "X1"
            colnames(df_coords)[colnames(df_coords) == "y" ] = "X2"
            
            # Add clicked point
            box_plot = addLabel(box_plot, df_coords, RV_box_plot$clicks, input$label_type, input$multi_label, 0.03, "X1", "X2")
          }
        }
      }
    }
    
    RV_box_plot$plot = box_plot
  })
  
  # Table generation
  observeEvent(c(RV_fc_table$table, input$logfc_threshold_1, input$logfc_threshold_2, input$pval_threshold_1, input$pval_threshold_2, input$fc_yes, input$pval_yes), {
    current_table = RV_fc_table$table
    
    if (nrow(current_table) > 0) {
      logFC_lower_1 = input$logfc_threshold_1[1]
      logFC_upper_1 = input$logfc_threshold_1[2]
      logFC_lower_2 = input$logfc_threshold_2[1]
      logFC_upper_2 = input$logfc_threshold_2[2]
      pval_1 = input$pval_threshold_1
      pval_2 = input$pval_threshold_2
      
      # Subsetting tables
      if (length(RV_fc_table$selector) == 1) {
        RV_fc_table$logfc_min_1 = round(min(current_table$logFC), 2)
        RV_fc_table$logfc_max_1 = round(max(current_table$logFC), 2)
        
        # Conditionals
        if (input$fc_yes == "Yes") {
          current_table = subset(current_table, logFC >= logFC_lower_1 & logFC < logFC_upper_1)
        }
        if (input$pval_yes == "Yes") {
          current_table = subset(current_table, PValue <= pval_1)
        }
      }
      
      else if (length(RV_fc_table$selector) == 2) {
        RV_fc_table$logfc_min_1 = round(min(current_table$logFC_1), 2)
        RV_fc_table$logfc_max_1 = round(max(current_table$logFC_1), 2)
        RV_fc_table$logfc_min_2 = round(min(current_table$logFC_2), 2)
        RV_fc_table$logfc_max_2 = round(max(current_table$logFC_2), 2)
        
        # Conditionals
        if (input$fc_yes == "Yes") {
          current_table = subset(current_table, logFC_1 >= logFC_lower_1 & logFC_1 < logFC_upper_1)
          current_table = subset(current_table, logFC_2 >= logFC_lower_2 & logFC_2 < logFC_upper_2)
        }
        if (input$pval_yes == "Yes") {
          current_table = subset(current_table, PValue_1 <= pval_1)
          current_table = subset(current_table, PValue_2 <= pval_2)
        }
      }
    }
    
    RV_fc_table$display_table = current_table
  })
  
  # Correlation fold change plot generation
  observeEvent(c(RV_fc_table$display_table, RV_fc_plot$clicks, RV_fc_table$subset_table), {
    plot_table = RV_fc_table$display_table
    
    return_plot = NULL
    if (length(RV_fc_table$selector) < 2) {
      return_plot = ggplot() + xlab("X1") + ylab("X2") + geom_text(aes(x=0, y=0, label="Two fold change tables must be selected in order to plot correlation."), size = 8, colour = plot_text_color)
    }
    else if (length(RV_fc_table$selector) == 2) {
      return_plot = ggplot(plot_table, aes(x = logFC_1, y = logFC_2)) +
        geom_point(size = 3, shape = 1, colour = plot_text_color) +
        ggtitle("logFC Correlation")
      
      df_coords = data.frame("X1"=plot_table$logFC_1, "X2"=plot_table$logFC_2, "gene"=plot_table$gene)
      highlights = data.frame("X1"=RV_fc_table$subset_table$logFC_1, "X2"=RV_fc_table$subset_table$logFC_2, "gene"=RV_fc_table$subset_table$gene)
      return_plot = addLabel(return_plot, df_coords, RV_fc_plot$clicks, "gene", input$multi_label, 0.03)
      return_plot = tableLabel(return_plot, df_coords, highlights, "gene", 0.03)
    }
    
    RV_fc_plot$plot = return_plot
  })
  
  ##########################################################################################################################
  # UPDATE UI
  ##########################################################################################################################
  
  # Update gene list based on dataset(s) chosen
  observeEvent(RV_UI$gene_names, {
    updateSelectizeInput(session, inputId = "gene",
                         choices = RV_UI$gene_names,
                         server = TRUE,
                         options = list(maxItems = 2, placeholder = "Click here and select up to 2 genes to plot"))
  })
  
  # Update pc list based on dataset(s) chosen
  observeEvent(RV_UI$pc_choices, {
    updateSelectizeInput(session, inputId = "pcs_chosen",
                         choices = RV_UI$pc_choices,
                         selected = c("PC1", "PC2"),
                         server = TRUE,
                         options = list(maxItems = 2, placeholder = "Click here and select 2 PCs to plot"))
  })
  
  # Update color category
  observeEvent(RV_UI$condition_names, {
    names = RV_UI$condition_names
    if (length(names) > 0) {
      updateSelectizeInput(session, inputId = "condition",
                           choices = names,
                           selected = "ExperimentalDesign",
                           options = list(maxItems = 2, placeholder = "Click here and select up to 2 conditions to plot"))
    }
    else {
      updateSelectizeInput(session, inputId = "condition",
                           choices = c("ExperimentalDesign"),
                           selected = "ExperimentalDesign")
    }
  })
  
  # Update label type
  observeEvent(RV_UI$label_names, {
    updateSelectizeInput(session, inputId = "label_type",
                         choices = c(RV_UI$label_names, "None"),
                         selected = "Sample_ID")
    
  })
  
  # Update co-expression threshold
  observeEvent(RV_data$exp_max, {
    if (RV_data$r_type == "l") {
      updateSliderInput(session,
                        inputId = "coex_threshold",
                        label = NULL,
                        value = 1,
                        min = RV_data$exp_min,
                        max = RV_data$exp_max,
                        step = 1)
    }
    else if (RV_data$r_type == "r") {
      updateSliderInput(session,
                        inputId = "coex_threshold",
                        label = NULL,
                        value = 10,
                        min = 0,
                        max = RV_data$exp_max,
                        step = 10)
    }
    
  })
  
  # Update fold change table options
  observeEvent(RV_UI$fc_names, {
    print("FC FRAMES UPDATED TO:")
    print(RV_UI$fc_names)
    updateSelectizeInput(session, inputId = "fc_table_select",
                         choices = RV_UI$fc_names,
                         server = TRUE,
                         options = list(maxItems = 2, placeholder = "Click here and select up to 2 fc tables"))
    
  })
  
  # Update p-value threshold
  observeEvent(RV_fc_table$logfc_max_1, {
    
    fc_label_1 = NULL
    fc_label_2 = NULL
    pval_label_1 = NULL
    pval_label_2 = NULL
    
    if (length(RV_fc_table$selector) == 2) {
      fc_label_1 = "Table 1"
      fc_label_2 = "Table 2"
      pval_label_1 = "Table 1"
      pval_label_2 = "Table 2"
    }
    
    updateSliderInput(session,
                      inputId = "logfc_threshold_1",
                      label = fc_label_1,
                      value = c(RV_fc_table$logfc_min_1, RV_fc_table$logfc_max_1),
                      min = RV_fc_table$logfc_min_1,
                      max = RV_fc_table$logfc_max_1,
                      step = 0.1)
    
    updateSliderInput(session,
                      inputId = "logfc_threshold_2",
                      label = fc_label_2,
                      value = c(RV_fc_table$logfc_min_2, RV_fc_table$logfc_max_2),
                      min = RV_fc_table$logfc_min_2,
                      max = RV_fc_table$logfc_max_2,
                      step = 0.1)
    
    updateSliderInput(session,
                      inputId = "pval_threshold_1",
                      label = pval_label_1,
                      value = 0.05,
                      min = 0,
                      max = 0.1,
                      step = 0.001)
    
    updateSliderInput(session,
                      inputId = "pval_threshold_2",
                      label = pval_label_2,
                      value = 0.05,
                      min = 0,
                      max = 0.1,
                      step = 0.001)
    
  })
  
  ##########################################################################################################################
  # RENDERS
  ##########################################################################################################################
  
  interactive_theme = theme(
    plot.background = element_blank(),
    # legend.background = element_blank(),
    legend.key = element_blank(),
    axis.line = element_line(colour = plot_text_color),
    panel.grid = element_line(colour = plot_text_color, size = 0.1),
    axis.text = element_text(colour = plot_text_color),
    legend.title = element_blank()
  )
  
  # PCA plot render
  output$pca_plot = renderPlot({
    RV_pca_plot$plot + interactive_theme + theme(text = element_text(colour=plot_text_color, size=22))
  }, height=700, width=925, bg = "transparent")
  
  # MDS plot render
  output$mds_plot = renderPlot({
    RV_mds_plot$plot + interactive_theme + theme(text = element_text(colour=plot_text_color, size=22))
  }, height=700, width=925, bg = "transparent")
  
  # Box plot render
  output$box_plot = renderPlot({
    RV_box_plot$plot + interactive_theme + theme(text = element_text(colour=plot_text_color, size=22))
  }, height=700, width=925, bg = "transparent")
  
  # Corr plot render
  output$fc_plot = renderPlot({
    RV_fc_plot$plot + interactive_theme + theme(text = element_text(colour=plot_text_color, size=22))
  }, height=700, width=925, bg = "transparent")
  
  # Fold change table render
  output$fold_change_table = DT::renderDataTable({
    datatable(RV_fc_table$display_table, rownames = F, options = list(searching = T, lengthChange = F, pageLength = 50, searchHighlight = T))
  })
  
  output$heatmap = renderPlotly({
    heatmaply(log(RV_heatmap$matrix+1))
              #Rowv = F,
              #Colv = F,
              #label_names = c("gene", "sample", "exp"),
              #xlab = "sample",
              #ylab = "gene",
              #showticklabels = T,
              #fontsize_row = 4,
              #fontsize_col = 4,
              #key.title = "Exp",
              #margins = c(100,100,100,100))
  })
  
  
  
  # Download Button
  output$downloadData = downloadHandler(
    filename = function() {
      if (input$plot_tab == RV_pca_plot$title) {
        paste("my_pcaPlot", ".png", sep="")
      }
      else if (input$plot_tab == RV_mds_plot$title) {
        paste("my_mdsPlot", ".png", sep="")
      }
      else if (input$plot_tab == RV_box_plot$title) {
        paste("my_boxPlot", ".png", sep="")
      }
      else if (input$plot_tab == RV_fc_plot$title) {
        if (input$fold_change_panel == "Table") {
          paste("my_fcTable", ".txt", sep="")
        }
        else if ((input$fold_change_panel == "Plot")) {
          paste("my_correlationPlot", ".png", sep="")
        }
      }
      
    },
    content = function(file) {
      new_theme = theme(text = element_text(colour=plot_text_color, size=20), 
                        panel.background = element_rect(fill = background_color),
                        plot.background = element_rect(fill = background_color),
                        panel.border = element_rect(colour=plot_text_color),
                        plot.margin=grid::unit(c(1,1,1,1),"cm"))
      
      output_plot = NULL
      if (input$plot_tab == RV_pca_plot$title) {
        output_plot = RV_pca_plot$plot + interactive_theme
        output_plot = output_plot + new_theme
        ggsave(file, plot = output_plot, device = "png", height = 9, width = 12)
      }
      else if (input$plot_tab == RV_mds_plot$title) {
        output_plot = RV_mds_plot$plot + interactive_theme
        output_plot = output_plot + new_theme
        ggsave(file, plot = output_plot, device = "png", height = 9, width = 12)
      }
      else if (input$plot_tab == RV_box_plot$title) {
        output_plot = RV_box_plot$plot + interactive_theme
        output_plot = output_plot + new_theme
        ggsave(file, plot = output_plot, device = "png", height = 9, width = 12)
      }
      else if (input$plot_tab == RV_fc_plot$title) {
        if (input$fold_change_panel == "Table") {
          write.table(RV_fc_table$display_table, file, sep = "\t", quote = F, col.names = T, row.names = F)
        }
        else if ((input$fold_change_panel == "Plot")) {
          output_plot = RV_fc_plot$plot + interactive_theme
          output_plot = output_plot + new_theme
          ggsave(file, plot = output_plot, device = "png", height = 9, width = 12)
        }
      }
    }
  )
  
  ##########################################################################################################################
  # CONDITIONAL BOOLEANS
  ##########################################################################################################################
  
  # Conditional panel for main panel
  output$plot_display = reactive({
    nrow(RV_data$metadata) > 0
  })
  outputOptions(output, "plot_display", suspendWhenHidden = FALSE)
  
  # Conditional panel for whether pc choices is displayed
  output$pc_choose = reactive({
    input$plot_tab == "PCA Scatter Plot"
  })
  outputOptions(output, "pc_choose", suspendWhenHidden = FALSE)
  
  # Conditional panel for co-expression threshold slider
  output$co_expression_check = reactive({
    ncol(RV_data$selected_genes) > 1 & input$binary_button == "On" & (input$plot_tab == "MDS Scatter Plot" | input$plot_tab == "PCA Scatter Plot")
  })
  outputOptions(output, "co_expression_check", suspendWhenHidden = FALSE)
  
  # Conditional panel for whether binary output is possible
  output$binary_possible = reactive({
    ncol(RV_data$selected_genes) == 2 & (input$plot_tab == "MDS Scatter Plot" | input$plot_tab == "PCA Scatter Plot")
  })
  outputOptions(output, "binary_possible", suspendWhenHidden = FALSE)
  
  # Condition panel for fold change data
  output$fc_tables_yes = reactive({
    length(RV_UI$fc_names) == 0
  })
  outputOptions(output, "fc_tables_yes", suspendWhenHidden = FALSE)
  
  
  # Condition panel for fold change side bar
  output$fc_check = reactive({
    input$plot_tab == "Fold Change"
  })
  outputOptions(output, "fc_check", suspendWhenHidden = FALSE)
  
  # Condition panel for not fold change side bar
  output$not_fc_check = reactive({
    input$plot_tab != "Fold Change"
  })
  outputOptions(output, "not_fc_check", suspendWhenHidden = FALSE)
  
  # Condition panel for pval/fold change thresholds
  output$fc_sliders = reactive({
    length(RV_fc_table$selector) > 0
  })
  outputOptions(output, "fc_sliders", suspendWhenHidden = FALSE)
  
  # Condition hide data selection
  output$data_check = reactive({
    RV_data$generated == F
  })
  outputOptions(output, "data_check", suspendWhenHidden = FALSE)
  
  # Data button was clicked: show data selection
  output$data_button = reactive({
    RV_data$generated == T
  })
  outputOptions(output, "data_button", suspendWhenHidden = FALSE)
  
  # Data button was clicked: show data selection
  output$multi_fold_data = reactive({
    length(RV_fc_table$selector) == 2
  })
  outputOptions(output, "multi_fold_data", suspendWhenHidden = FALSE)
  
  # Instruction conditionals
  output$genes_yes = reactive({
    length(RV_data$selected_genes) > 1
  })
  outputOptions(output, "genes_yes", suspendWhenHidden = FALSE)
  
  output$genes_no = reactive({
    length(RV_data$selected_genes) <= 1
  })
  outputOptions(output, "genes_no", suspendWhenHidden = FALSE)
  
  output$heatmap_yes = reactive({
    nrow(RV_heatmap$matrix) > 2
  })
  outputOptions(output, "heatmap_yes", suspendWhenHidden = FALSE)
  
  session$onSessionEnded(stopApp)
}