function(input, output, session) { 
  
  ##########################################################################################################################
  # REACTIVE VALUES
  ##########################################################################################################################
  
  RV_data = reactiveValues(matrix = data.frame(),
                           table = data.frame(),
                           selected_genes = data.frame(),
                           exp_min = 0,
                           exp_max = 100,
                           r_type = "l",
                           DGE_matrix = NULL,
                           all_mouse = T,
                           generated = F)
  
  RV_UI = reactiveValues(gene_names = list(),
                         label_names = list(),
                         condition_names = list(),
                         fc_names = list())
  
  RV_interactive_plot = reactiveValues(plot = ggplot(),
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
  
  RV_dynamic = reactiveValues(val = 0)
  
  reset = function() {
    # Blank current data
    RV_data$matrix = data.frame()
    RV_data$table = data.frame()
    RV_data$selected_genes = data.frame("Lib_ID"=character())
    RV_data$exp_min = 0
    RV_data$exp_max = 100
    RV_data$DGE_matrix = data.frame()
    RV_data$all_mouse = T
    
    RV_UI$gene_names = list()
    RV_UI$label_names = list()
    RV_UI$condition_names = list()
    RV_UI$fc_names = list()
    
    RV_interactive_plot$plot = ggplot()
    RV_interactive_plot$clicks = NULL
    
    RV_box_plot$plot = ggplot()
    RV_box_plot$clicks = NULL
    
    RV_fc_plot$plot = ggplot()
    RV_fc_plot$clicks = NULL
    
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
    
    reset()
    n = 0
    
    if (input$normalize == "l") {
      RV_data$r_type = "l"
    }
    else if (input$normalize == "r") {
      RV_data$r_type = "r"
      n = n + 6
    }
    
    if (!is.null(input$file1$datapath)) {
      user_table = read.table(input$file1$datapath,
                     header = T,
                     sep = "\t",
                     quote = "")
      
      print(head(user_table))
    }
    
    
    if (length(input$dataset) > 0) {
      withProgress(message = 'Initializing', value = 0, {
        
        selected_sets = input$dataset
        working_matrix = data.frame("gene"=character())
        lib_ID_names = data.frame()
        
        n = n + (2*length(selected_sets) + 3)
        
        # IMPORT DATA
        
        # Import sample info tables
        label_info = annotations[[1]]
        if (length(annotations) > 1) {
          for (i in 2:length(annotations)) {
            label_info = merge(label_info, annotations[[i]], all = T)
          }
        }
        label_info[is.na(label_info)] = "N/A"
        
        RV_UI$label_names = colnames(label_info)
        RV_UI$condition_names = colnames(label_info)
        
        # Check selected sets for human data
        for (i in 1:length(selected_sets)) {
          organism = datasets$organism[datasets$Datasets == selected_sets[i]]
          if ("Human" == organism) {
            RV_data$all_mouse = F
          }
        }
        
        # NEW ADDITION: FOLD CHANGE NAME SAVE
        fc_names = list()
        fc_frames = list()
        fc_index = 1
        
        # Import gene tables and fold change tables
        for (i in 1:length(selected_sets)) {
          incProgress(1/n, message = paste("Dataset ", i, " of ", length(selected_sets),  ": ", "Loading data", sep = ""))
          
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
          
          # NEW ADDITION: FOLD CHANGE TABLE IMPORT
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
          
          # merge matrices into one matrix
          incProgress(1/n, message = paste("Dataset ", i, " of ", length(selected_sets),  ": ", "Merging data", sep = ""))
          
          colnames(data)[1] = "gene"
          if (RV_data$all_mouse == F) {
            if (organism == "Human") {
              master = data.frame("Master"=master_name_table$Master, "gene"=master_name_table$Human)
              data = merge(data, master, by = "gene")
              data$gene = NULL
              colnames(data)[colnames(data) == "Master"] = "gene"
            }
            else if (organism == "Mouse") {
              master = data.frame("Master"=master_name_table$Master, "gene"=master_name_table$Mouse)
              data = merge(data, master, by = "gene")
              data$gene = NULL
              colnames(data)[colnames(data) == "Master"] = "gene"
            }
            else {
              # And disconnect...
              print("Incorrect organism input")
            }
          }
          
          merged_matrix = merge(working_matrix, data, by = "gene", all = T)
          merged_matrix[is.na(merged_matrix)] = 0
          
          working_matrix = merged_matrix
        }
        
        incProgress(1/n, message = "Merging all data")
        
        # Set matrix variable
        working_matrix = working_matrix[!duplicated(working_matrix$gene), ]
        row.names(working_matrix) = working_matrix$gene
        working_matrix$gene = NULL
        
        # MDS function
        incProgress(1/n, message = "Calculating MDS")
        transposed_matrix = data.frame(t(working_matrix))
        
        d = dist(transposed_matrix)
        fit = cmdscale(d,eig=TRUE, k=2)
        
        MDS_data = data.frame(fit$points)
        MDS_data$Lib_ID = row.names(MDS_data)
        MDS_data = merge(MDS_data, label_info, by = "Lib_ID")
 
        # Set gene names variable
        RV_UI$gene_names = row.names(working_matrix)
        RV_data$matrix = working_matrix
        RV_data$table = MDS_data
        
        # Set plot names for saving
        title_list = unique(MDS_data$Dataset)
        title_name = ""
        for (i in 1:length(title_list)) {
          title_name = paste(title_name, title_list[i], ",", sep = "")
        }
        title_name = substr(title_name, 1, nchar(title_name) - 1)
        
        # ALEX'S SCRIPT FOR RUNNING DGE
        if (RV_data$r_type == "r") {
          annotations = RV_data$table
          matrix = RV_data$matrix
          cpm.Lib_ID.cutoff = 2
          min.cpm = 1
          
          matrix = matrix[,which(!apply(matrix,2,FUN = function(x){all(x == 0)}))]
          
          # Set up edgeR expression object
          incProgress(1/n, message = "Calculating DGE")
          DGE_data = DGEList(counts=matrix)
          
          keep <- rowSums(cpm(DGE_data)>min.cpm) >= cpm.Lib_ID.cutoff
          DGE_data <- DGE_data[keep, , keep.lib.sizes=FALSE]
          
          if (is.na(DGE_data$counts[1,1]) | DGE_data$counts[1,1] < 0) {
            RV_data$DGE_matrix = NULL
          }
          else {
            # Perform simple exact test on genotype
            incProgress(1/n, message = "Running DGE filter 1 of 3")
            DGE_data <- calcNormFactors(DGE_data)
            incProgress(1/n, message = "Running DGE filter 2 of 3")
            DGE_data <- estimateCommonDisp(DGE_data)
            incProgress(1/n, message = "Running DGE filter 3 of 3")
            DGE_data <- estimateTagwiseDisp(DGE_data)
            
            pseudo_counts = DGE_data$pseudo.counts
            pseudo_counts[pseudo_counts < 0] = 0
            
            RV_data$DGE_matrix = pseudo_counts
            write.table(pseudo_counts, "FIND_ME.txt", sep = "\t", quote = F, col.names = T, row.names = T)
          }
        }
        
        RV_fc_table$table_list = fc_frames
        fc_names = unlist(lapply(fc_names, function(i) unlist(strsplit(i, "[.]"))[1]))
        RV_UI$fc_names = fc_names
        
        incProgress(1/n, message = "Done")
      })
    }
  })
  
  ##########################################################################################################################
  # DYNAMIC DATA UPDATES
  ##########################################################################################################################
  
  # Activated by gene list entry: Update reactive values
  observeEvent(input$gene, ignoreNULL = F, {
    return_frame = data.frame("Lib_ID"=character())
    if (!is.null(input$gene)) {
      for (i in 1:length(input$gene)) {
        current_row = RV_data$matrix[input$gene[[i]], ]
        transposed_gene = data.frame(t(current_row))
        transposed_gene$Lib_ID = row.names(transposed_gene)
        return_frame = merge(return_frame, transposed_gene, by = "Lib_ID", all = T)
      }
    }
    RV_data$selected_genes = return_frame
  })
  
  # Activated for clicked points
  observeEvent(input$plot_click, {
    if (!is.null(input$plot_click) & input$plot_tab == RV_interactive_plot$title) {
      RV_interactive_plot$clicks = input$plot_click
    }
  })
  observeEvent(input$plot_click_2, {
    if (!is.null(input$plot_click_2) & input$plot_tab == RV_box_plot$title) {
      RV_box_plot$clicks = input$plot_click_2
    }
  })
  observeEvent(input$plot_click_3, {
    if (!is.null(input$plot_click_3) & input$plot_tab == RV_fc_plot$title) {
      RV_fc_plot$clicks = input$plot_click_3
    }
  })
  
  # Select new data button toggle
  observeEvent(input$generate_button, {
    if (length(input$dataset) > 0) {
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
  observeEvent(c(RV_data$table,
                 RV_data$selected_genes,
                 RV_data$DGE_matrix,
                 RV_interactive_plot$clicks,
                 RV_box_plot$clicks,
                 RV_fc_plot$clicks,
                 input$condition,
                 input$label_type,
                 input$multi_label,
                 input$coex_threshold), {
    RV_dynamic$val = RV_dynamic$val + 1
  })
  
  # Create/update interactive plot
  observeEvent(RV_dynamic$val, {
    if (nrow(RV_data$table) > 0) {
      
      plot_data = RV_data$table
      gene_data = RV_data$selected_genes
      plot = NULL
      
      plot_data = findConditionCol(plot_data, input$condition)
      
      # No genes are present
      if (ncol(RV_data$selected_genes) == 1) {
        plot = ggplot(plot_data, aes(x = X1, y = X2, colour = factor(current_plot_group))) +
          ggtitle("MDS Plot") +
          geom_point(size = 4, alpha = 0.7) +
          scale_colour_manual(values=cbPalette)
      }
      
      # 1 gene is present
      else if (ncol(RV_data$selected_genes) == 2) {
        
        gene_name = colnames(gene_data[2])
        colnames(gene_data) = c("Lib_ID", "Exp")
        plot_data = merge(plot_data, gene_data, by = "Lib_ID", all.x = T)
        
        current_max = max(plot_data$Exp)
        current_min = min(plot_data$Exp)
        current_mid = mean(plot_data$Exp)
        
        plot = ggplot(plot_data, aes(x = X1, y = X2, shape = factor(current_plot_group), colour = Exp)) +
          ggtitle(gene_name) +
          geom_point(size = 3) +
          scale_shape_manual(values=shapes) +
          scale_colour_gradient2(low="blue", mid="grey", high=theme_orange, midpoint=current_mid, space="Lab",
                                 breaks=c(current_min, current_max),
                                 labels=c(round(current_min, 2), round(current_max, 2)),
                                 limits=c(current_min, current_max))
      }
      
      # 2 genes are present
      else if (ncol(RV_data$selected_genes) == 3) {
        
        gene_name_1 = colnames(gene_data[2])
        gene_name_2 = colnames(gene_data[3])
        colnames(gene_data) = c("Lib_ID", "Exp1", "Exp2")
        plot_data = merge(plot_data, gene_data, by = "Lib_ID", all.x = T)
        
        RV_data$exp_max = round(max(max(plot_data$Exp1), max(plot_data$Exp2)), -1) + 20
        RV_data$exp_min = round(min(min(plot_data$Exp1), min(plot_data$Exp2)), -1) - 20
        
        threshold = input$coex_threshold
        
        current_red = subset(plot_data, (Exp1 >= threshold & Exp2 < threshold))
        current_green = subset(plot_data, (Exp2 >= threshold & Exp1 < threshold))
        current_yellow = subset(plot_data, (Exp1 >= threshold & Exp2 >= threshold))
        current_grey = subset(plot_data, (Exp1 < threshold & Exp2 < threshold))
        
        colours = c("G1"="red", "G2"="green", "Both"="#dfdf02", "Neg"="grey")
        plot = ggplot() +
          ggtitle(paste("G1 = ", gene_name_1, ", G2 = ", gene_name_2, sep = "")) +
          scale_colour_manual(values = colours) +
          scale_shape_manual(values = shapes)
        
        if (nrow(current_red) > 0) { plot = plot + geom_point(data = current_red, aes(x = X1, y = X2, colour = "G1", shape = factor(current_plot_group)), size = 3) }
        if (nrow(current_green) > 0) { plot = plot + geom_point(data = current_green, aes(x = X1, y = X2, colour = "G2", shape = factor(current_plot_group)), size = 3) }
        if (nrow(current_yellow) > 0) { plot = plot + geom_point(data = current_yellow, aes(x = X1, y = X2, colour = "Both", shape = factor(current_plot_group)), size = 3) }
        if (nrow(current_grey) > 0) { plot = plot + geom_point(data = current_grey, aes(x = X1, y = X2, colour = "Neg", shape = factor(current_plot_group)), size = 3) }
      }
      
      # Add clicked point
      plot = addLabel(plot, plot_data, RV_interactive_plot$clicks, input$label_type, input$multi_label, 0.03)
      
      RV_interactive_plot$plot = plot
      
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
            current_data = data.frame("Lib_ID"=colnames(pseudo_genes), "count"=pseudo_genes[gene_name, ])
            plot_data = merge(current_data, RV_data$table, by = "Lib_ID")
            plot_data = findConditionCol(plot_data, input$condition)
            
            box_plot = ggplot(plot_data, aes(x=factor(current_plot_group), y = log(count), fill = factor(current_plot_group))) +
              geom_boxplot(colour = plot_text_color, width = 0.5) +
              # geom_jitter(aes(fill = factor(current_plot_group)), size = 4, shape = 21, colour = plot_text_color) +
              geom_point(aes(fill = factor(current_plot_group)), size = 4, shape = 21, colour = plot_text_color) +
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
            colnames(df_coords)[colnames(df_coords) == "x" ] <- "X1"
            colnames(df_coords)[colnames(df_coords) == "y" ] <- "X2"
            
            # Add clicked point
            box_plot = addLabel(box_plot, df_coords, RV_box_plot$clicks, input$label_type, input$multi_label, 0.03)
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
  
  # Update color category
  observeEvent(RV_UI$condition_names, {
    names = RV_UI$condition_names
    if (length(names) > 0) {
      updateSelectizeInput(session, inputId = "condition",
                           choices = names,
                           selected = "ExperimentalDesign")
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
                         choices = c(RV_UI$label_names, "none"),
                         selected = "Lib_ID")

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
      fc_label_1 = "Dataset 1"
      fc_label_2 = "Dataset 2"
      pval_label_1 = "Dataset 1"
      pval_label_2 = "Dataset 2"
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
  
  # Interactive plot render
  output$interactive_plot = renderPlot({
    RV_interactive_plot$plot + interactive_theme + theme(text = element_text(colour=plot_text_color, size=22))
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
  
  # Download Button
  output$downloadData = downloadHandler(
    filename = function() {
      if (input$plot_tab == RV_interactive_plot$title) {
        paste("my_mdsplot", ".png", sep="")
      }
      else if (input$plot_tab == RV_box_plot$title) {
        paste("my_boxplot", ".png", sep="")
      }
      else if (input$plot_tab == RV_fc_plot$title) {
        if (input$fold_change_panel == "Table") {
          paste("my_fctable", ".txt", sep="")
        }
        else if ((input$fold_change_panel == "Plot")) {
          paste("my_correlationplot", ".png", sep="")
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
      if (input$plot_tab == RV_interactive_plot$title) {
        output_plot = RV_interactive_plot$plot + interactive_theme
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
  output$MDS_display = reactive({
    nrow(RV_data$table) > 0
  })
  outputOptions(output, "MDS_display", suspendWhenHidden = FALSE)
  
  # Conditional panel for co-expression threshold slider
  output$co_expression_check = reactive({
    ncol(RV_data$selected_genes) == 3 & input$plot_tab == "MDS Scatter Plot"
  })
  outputOptions(output, "co_expression_check", suspendWhenHidden = FALSE)
  
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
  
  session$onSessionEnded(stopApp)
}