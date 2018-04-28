library(rsconnect)
library(shiny)
library(shinythemes)
library(plotly)
library(ggplot2)
library(devtools)
library(Matrix)
library(plyr)
library(edgeR)
library(GenomicFeatures)
library(colorRamps)
library(DT)
library(magrittr)
library(htmlwidgets)

theme_set(theme_bw(15))
plot_text_color = "#2b3e50"
# "#f0f0f0"
background_color = "white"
# "#2b3e50"
theme_orange = "#ffb62f"
cbPalette = c(theme_orange, "#56B4E9", "#66FF66", "#f63093", "#0000ff", "#ff0000", "#CC79A7", "#D55E00", "#009E73", "#33FFFF", "#FF3399", "#F0E442", "#FF3333", "#66FFCC", "#0000CC", "#999999", primary.colors(), blue2green2red(50), magenta2green(50))
shapes = c(16,15,17,18,6,3,4,8,10,7,9,13,12,1,0,2,5,97:122,65:90,1:25,97:122,65:90)
dot_size = 5
dot_alpha = 0.6

# Functions
# addLabel takes a plot, the data frame that was used to make the plot, clicked coordinates, and a label and creates a new plot with the
# label generated onto the plot
old = data.frame("x"="x", "y"="y", stringsAsFactors = F)
addLabel = function(p, df, clicks, label, sm, ratio, x_name, y_name) {
  
  return_plot = p
  
  if (label != "None" & !is.null(clicks) & nrow(df) > 0) {
    all_coords = df
    colnames(all_coords)[colnames(all_coords) == x_name] = "X_plot_coordinate"
    colnames(all_coords)[colnames(all_coords) == y_name] = "Y_plot_coordinate"
    
    curser_x = round(clicks$x, 3)
    curser_y = round(clicks$y, 3)
    
    for (i in 1:nrow(all_coords)) {
      coord_x = round(all_coords$X_plot_coordinate[i], 3)
      coord_y = round(all_coords$Y_plot_coordinate[i], 3)
      
      all_coords$dist[i] = sqrt((curser_x - coord_x)^2 + (curser_y - coord_y)^2)
    }
    
    row = which(grepl(min(all_coords$dist), all_coords$dist))
    selected_row = all_coords[row, ]
    
    # Merge label names if more than one point is present at a given coordinate
    if (nrow(selected_row) > 1) {
      new_label = ""
      for (i in 1:nrow(selected_row)) {
        new_label = paste(new_label, selected_row$cell[i], ", ", sep = "")
      }
      new_label = substr(new_label, 1, nchar(new_label)-2)
      
      selected_row$cell[1] = new_label
      selected_row = selected_row[1,]
    }
    
    x_range = abs(max(all_coords$X_plot_coordinate) - min(all_coords$X_plot_coordinate))
    y_range = abs(max(all_coords$Y_plot_coordinate) - min(all_coords$Y_plot_coordinate))
    
    x_mid = min(all_coords$X_plot_coordinate) + x_range/2
    y_mid = min(all_coords$Y_plot_coordinate) + y_range/2
    h_just =  x_range * ratio
    v_just =  y_range *ratio
    
    quad = 0
    
    if (selected_row$X_plot_coordinate[1] > x_mid) {
      h_just = -h_just
      quad = 1
    }
    if (selected_row$Y_plot_coordinate[1] > y_mid) {
      v_just = -v_just
    }
    
    colnames(selected_row)[colnames(selected_row) == label] = "new_label"
    colnames(all_coords)[colnames(all_coords) == label] = "new_label"
    
    new_type = as.character(selected_row[, "new_label"])
    
    if (grepl(", ", new_type)) {
      new_type = unlist(strsplit(new_type, ","))[1]
    }
    
    # Add text adjustments to plot
    selected_row$h_just = h_just
    selected_row$v_just = v_just
    selected_row$quad = quad
    
    return_plot = return_plot + geom_text(data = selected_row, aes(x=X_plot_coordinate + h_just, y=Y_plot_coordinate + v_just, label = new_label), colour = plot_text_color, size=7, hjust = quad)
    
    if (sm == "s") {
      return_plot = return_plot + geom_point(data = selected_row, aes(x=X_plot_coordinate, y=Y_plot_coordinate), colour = plot_text_color, alpha = 1, size = 7, shape = 1)
    }
    else if (sm == "m") {
      new_data = subset(all_coords, new_label == selected_row$new_label)
      new_data$h_just = h_just
      new_data$v_just = v_just
      return_plot = return_plot + geom_point(data = new_data, aes(x=X_plot_coordinate, y=Y_plot_coordinate), colour = plot_text_color, alpha = 1, size = 7, shape = 1)
    }
  }
  
  return(return_plot)
}
# Label based on data table for fold change
tableLabel = function(p, df, highlight_rows, label, ratio) {
  return_plot = p
  
  if (nrow(df) > 0 & nrow(highlight_rows) > 0) {
    all_coords = df
    
    selected_rows = highlight_rows
    
    x_range = abs(max(all_coords$X1) - min(all_coords$X1))
    y_range = abs(max(all_coords$X2) - min(all_coords$X2))
    
    x_mid = min(all_coords$X1) + x_range/2
    y_mid = min(all_coords$X2) + y_range/2
    h_just =  x_range * ratio
    v_just =  y_range *ratio
    selected_rows$h_just = h_just
    selected_rows$v_just = v_just
    
    colnames(selected_rows)[colnames(selected_rows) == label] = "new_label"
    
    return_plot = return_plot + geom_text(data = selected_rows, aes(x=X1 + h_just, y=X2 + v_just, label = new_label), colour = plot_text_color, size=5, hjust = 0)
    return_plot = return_plot + geom_point(data = selected_rows, aes(x = X1, y = X2), colour = theme_orange, alpha = 0.5, size = 7, shape = 16)
  }
  
  return(return_plot)
}
# findConditionCol takes a dataframe and a string, and changes the name of the column that matches the string to 'current_plot_group'
findConditionCol = function(df, c, genes) {
  if (!is.null(c)) {
    # If any genes are selected
    if (length(genes) > 1) {
      if (c[1] == "Lib_ID") {
        plot_cond_1
        df$plot_cond_1 = df$Lib_ID
      }
      else {
        names(df)[names(df) == c[1]] = "plot_cond_1"
      }
    }
    # If no genes are selected
    else {
      # If one condition is selected
      if (c[1] == "Lib_ID") {
        df$plot_cond_1 = df$Lib_ID
      }
      else {
        names(df)[names(df) == c[1]] = "plot_cond_1"
      }
      # If two conditions are selected
      if (length(c) == 2) {
        if (c[2] == "Lib_ID") {
          df$plot_cond_2 = df$Lib_ID
        }
        else {
          names(df)[names(df) == c[2]] = "plot_cond_2"
        }
      }
    }
  }
  return(df)
}
# Convert gene names to 'Master' name set, if non-mouse data exists in analysis
convertGenes = function(df, o) {
  print(o)
  master = data.frame("Master"=master_name_table$Master, "gene"=master_name_table[colnames(master_name_table) == o])
  colnames(master) = c("Master", "gene")
  data = merge(df, master, by = "gene")
  data$gene = NULL
  colnames(data)[colnames(data) == "Master"] = "gene"
  
  return(data)
}
# Function that creates ggplot based on given coordinates
graphGenes = function(df, gd, x_name, y_name, plot_name, b, b_thresh) {
  plot_data = df
  gene_data = gd
  
  colnames(plot_data)[colnames(plot_data) == x_name] = "X_plot_coordinate"
  colnames(plot_data)[colnames(plot_data) == y_name] = "Y_plot_coordinate"
  
  plot = NULL
  
  # No genes are present
  if (ncol(gene_data) == 1) {
    
    if (any(grepl("plot_cond_1", colnames(plot_data)))) {
      if (any(grepl("plot_cond_2", colnames(plot_data)))) {
        plot = ggplot(plot_data, aes(x = X_plot_coordinate, y = Y_plot_coordinate, colour = factor(plot_cond_1), shape = factor(plot_cond_2))) +
          ggtitle(plot_name) +
          xlab(x_name) +
          ylab(y_name) +
          geom_point(size = dot_size, alpha = dot_alpha) +
          scale_colour_manual(values=cbPalette) +
          scale_shape_manual(values = shapes)
        
      } else {
        plot = ggplot(plot_data, aes(x = X_plot_coordinate, y = Y_plot_coordinate, colour = factor(plot_cond_1))) +
          ggtitle(plot_name) +
          xlab(x_name) +
          ylab(y_name) +
          geom_point(size = dot_size, alpha = dot_alpha) +
          scale_colour_manual(values=cbPalette)
      }
      
    } else {
      plot = ggplot(plot_data, aes(x = X_plot_coordinate, y = Y_plot_coordinate)) +
        ggtitle(plot_name) +
        xlab(x_name) +
        ylab(y_name) +
        geom_point(colour = "grey", size = dot_size, alpha = dot_alpha) +
        scale_colour_manual(values=cbPalette) +
        scale_shape_manual(values = shapes)
    }
  }
  
  # 1 gene is present
  else if (ncol(gene_data) == 2) {
    
    if (b == "On") {
      
      gene_name = colnames(gene_data[2])
      colnames(gene_data) = c("Lib_ID", "Exp")
      plot_data = merge(plot_data, gene_data, by = "Lib_ID", all.x = T)
      
      threshold = b_thresh
      
      current_red = subset(plot_data, Exp >= threshold)
      current_grey = subset(plot_data, Exp < threshold)
      
      colours = c("Exp"="red", "Neg"="grey")
      plot = ggplot() +
        ggtitle(paste(plot_name, gene_name, sep = ": ")) +
        scale_colour_manual(values = colours) +
        scale_shape_manual(values = shapes) +
        xlab(x_name) +
        ylab(y_name)
      
      if (nrow(current_red) > 0) { plot = plot + geom_point(data = current_red, aes(x = X_plot_coordinate, y = Y_plot_coordinate, colour = "Exp", shape = factor(plot_cond_1)), size = dot_size, alpha = dot_alpha) }
      if (nrow(current_grey) > 0) { plot = plot + geom_point(data = current_grey, aes(x = X_plot_coordinate, y = Y_plot_coordinate, colour = "Neg", shape = factor(plot_cond_1)), size = dot_size, alpha = dot_alpha) }
      
    } else {
      
      gene_name = colnames(gene_data[2])
      colnames(gene_data) = c("Lib_ID", "Exp")
      plot_data = merge(plot_data, gene_data, by = "Lib_ID", all.x = T)
      
      current_max = max(plot_data$Exp)
      current_min = min(plot_data$Exp)
      current_mid = mean(plot_data$Exp)
      
      plot = ggplot(plot_data, aes(x = X_plot_coordinate, y = Y_plot_coordinate, shape = factor(plot_cond_1), colour = Exp)) +
        ggtitle(paste(plot_name, gene_name, sep = ": ")) +
        xlab(x_name) +
        ylab(y_name) +
        geom_point(size = dot_size, alpha = dot_alpha) +
        scale_shape_manual(values=shapes) +
        scale_colour_gradient2(low="blue", mid="grey", high=theme_orange, midpoint=current_mid, space="Lab",
                               breaks=c(current_min, current_max),
                               labels=c(round(current_min, 2), round(current_max, 2)),
                               limits=c(current_min, current_max))
    }
    
    
  }
  
  # 2 genes are present
  else if (ncol(gene_data) == 3) {
    
    gene_name_1 = colnames(gene_data[2])
    gene_name_2 = colnames(gene_data[3])
    colnames(gene_data) = c("Lib_ID", "Exp1", "Exp2")
    plot_data = merge(plot_data, gene_data, by = "Lib_ID", all.x = T)
    
    threshold = b_thresh
    
    current_red = subset(plot_data, (Exp1 >= threshold & Exp2 < threshold))
    current_green = subset(plot_data, (Exp2 >= threshold & Exp1 < threshold))
    current_yellow = subset(plot_data, (Exp1 >= threshold & Exp2 >= threshold))
    current_grey = subset(plot_data, (Exp1 < threshold & Exp2 < threshold))
    
    colours = c("Gene1"="red", "Gene2"="green", "Both"="#dfdf02", "Neg"="grey")
    plot = ggplot() +
      ggtitle(paste(plot_name, ": Gene1=", gene_name_1, ", Gene2=", gene_name_2, sep = "")) +
      xlab(x_name) +
      ylab(y_name) +
      scale_colour_manual(values = colours) +
      scale_shape_manual(values = shapes)
    
    if (nrow(current_red) > 0) { plot = plot + geom_point(data = current_red, aes(x = X_plot_coordinate, y = Y_plot_coordinate, colour = "Gene1", shape = factor(plot_cond_1)), size = dot_size, alpha = dot_alpha) }
    if (nrow(current_green) > 0) { plot = plot + geom_point(data = current_green, aes(x = X_plot_coordinate, y = Y_plot_coordinate, colour = "Gene2", shape = factor(plot_cond_1)), size = dot_size, alpha = dot_alpha) }
    if (nrow(current_yellow) > 0) { plot = plot + geom_point(data = current_yellow, aes(x = X_plot_coordinate, y = Y_plot_coordinate, colour = "Both", shape = factor(plot_cond_1)), size = dot_size, alpha = dot_alpha) }
    if (nrow(current_grey) > 0) { plot = plot + geom_point(data = current_grey, aes(x = X_plot_coordinate, y = Y_plot_coordinate, colour = "Neg", shape = factor(plot_cond_1)), size = dot_size, alpha = dot_alpha) }
  }
  
  return(plot)
}

# Load data
datasets = read.table("./data/datasets.txt", sep = "\t", header = T, stringsAsFactors = F)
master_name_table = read.table("./data/Hs.Mm.SymbolMerge.txt", sep = "\t", header = T, stringsAsFactors = F)
file_names = list.files("./sample_info")
annotations = list()

for (i in 1:length(file_names)) {
  path = paste("./sample_info", file_names[i], sep = "/")
  annotations[[i]] = read.table(path, sep = "\t", header = T, stringsAsFactors = F)
}

ui = source("./instance/ui.R",  local = TRUE)$value
server = source("./instance/server.R",  local = TRUE)$value

shinyApp(ui = ui, server = server)