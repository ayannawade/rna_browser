navbarPage("Nordlab: RNA Browser", theme = shinytheme("flatly"), selected = "Plot",
           
           tabPanel("Plot",
                    conditionalPanel(condition = 'output.data_check',
                                     tags$h1("Welcome!"),
                                     tags$h4("This is the Nordlab RNA Browser. It is an interactive tool for viewing rna-seq data. Follow the steps below and happy browsing."),
                                     tags$br(),
                                     tags$hr(),
                                     
                                     tags$h2("Step 1: Choose data"),
                                     tags$hr(),
                                     
                                     tabsetPanel(id = "data_tab",
                                                 tabPanel("Preloaded Data",
                                                          tags$br(),
                                                          tags$h4("Select one or more datasets from the box below and pick a normalization method."),
                                                          tags$h5("Pseudo count will provide more available functions than Relative expression."),
                                                          tags$br(),
                                                          selectizeInput(inputId = "dataset",
                                                                         label = NULL,
                                                                         choices = datasets$Datasets,
                                                                         multiple = T,
                                                                         options = list(placeholder = "Click here to view datasets")),
                                                          
                                                          radioButtons(inputId = "normalize",
                                                                       label = NULL,
                                                                       choices = c("Pseudo count"="r", "Relative expression"="l"),
                                                                       selected = "r",
                                                                       inline = T),
                                                          tags$br()
                                                 ),
                                                 
                                                 tabPanel("Upload Data",
                                                          
                                                          tags$br(),
                                                          tags$h4("Upload your own data in the following format:"),
                                                          tags$br(),
                                                          
                                                          tags$h4("1. Indicate the organism of your given data. Mouse and Human are currently the only compatible organisms with this browser."),
                                                          selectizeInput(inputId = "user_organism",
                                                                         label = NULL,
                                                                         selected = "Mouse",
                                                                         choices = c("Mouse", "Human")),
                                                          tags$br(),
                                                          tags$h4("2. Upload a gene counts matrix in tab delimited text format (.txt), where columns=samples and rows=genes. The first column should be labeled 'gene_id', indicating gene names."),
                                                          fileInput("user_matrix", label = NULL,
                                                                    multiple = F,
                                                                    accept = c("text/csv", "text/comma-separated-values,text/plain")),
                                                          tags$br(),
                                                          tags$h4("3. Metadata table in tab delimited text format (.txt). A column in this table must be named 'Sample_ID', indicating sample names, where samples name match the sample names of the gene counts matrix."),
                                                          fileInput("user_metadata", label = NULL,
                                                                    multiple = F,
                                                                    accept = c("text/csv",
                                                                               "text/comma-separated-values,text/plain")),
                                                          tags$br(),
                                                          tags$h4("3. (Optional) One or more fold change table(s) in tab delimited text format (.txt). See example for formatting."),
                                                          fileInput("user_fold_change", label = NULL,
                                                                    multiple = T,
                                                                    accept = c("text/csv",
                                                                               "text/comma-separated-values,text/plain"))
                                                 )
                                     ),
                                     
                                     tags$br(),
                                     tags$br(),
                                     tags$br(),
                                     
                                     actionButton(inputId = "generate_button",
                                                  label = "Generate")
                    ),
                    
                    conditionalPanel(condition = 'output.data_button',
                                     tags$hr(),
                                     actionButton(inputId = "select_new_data_button",
                                                  label = "Select New Data")
                    ),
                    
                    conditionalPanel(condition = "output.plot_display",
                                     tags$hr(),
                                     tags$h2("Step 2: Check out the plot tabs below"),
                                     tags$h4("Each tab will give you further instructions. All plots are interactive."),
                                     tags$hr(),
                                     
                                     sidebarPanel(
                                       tags$head(
                                         tags$style(type="text/css", "select { width: 350px; }"),
                                         tags$style(type="text/css", ".span4 { width: 350px; }"),
                                         tags$style(type="text/css", ".well { width: 350px; }")),
                                       
                                       tags$h3("- Plot Controls -"),
                                       
                                       tabsetPanel(id = "controller_tab",
                                                   tabPanel("View Options",
                                                            # tags$hr(),
                                                            tags$br(),
                                                            conditionalPanel(condition = 'output.not_fc_check',
                                                                             conditionalPanel(condition = 'output.pc_choose',
                                                                                              tags$h4("PCs (x-axis, y-axis)"),
                                                                                              selectizeInput(inputId = "pcs_chosen",
                                                                                                             label = NULL,
                                                                                                             multiple = T,
                                                                                                             selected = c("PC1", "PC2"),
                                                                                                             choices = c("PC1", "PC2")),
                                                                                              tags$hr()
                                                                             ),
                                                                             
                                                                             conditionalPanel(condition = 'output.genes_yes', tags$h4("Conditions (shape, disabled)")),
                                                                             conditionalPanel(condition = 'output.genes_no', tags$h4("Conditions (color, shape)")),
                                                                             
                                                                             selectizeInput(inputId = "condition",
                                                                                            label = NULL,
                                                                                            multiple = T,
                                                                                            selected = "ExperimentalDesign",
                                                                                            choices = "ExperimentalDesign"),
                                                                             
                                                                             tags$hr(),
                                                                             
                                                                             tags$h4("Gene Expression (color)"),
                                                                             selectizeInput(inputId = "gene",
                                                                                            label = NULL,
                                                                                            multiple = T,
                                                                                            choices = list()),
                                                                             
                                                                             conditionalPanel(condition = 'output.binary_possible', 
                                                                                              radioButtons(inputId = "binary_button",
                                                                                                           label = "Binary Scale",
                                                                                                           choices = c("On", "Off"),
                                                                                                           selected = "On",
                                                                                                           inline = T)
                                                                             ),
                                                                             
                                                                             conditionalPanel(condition = "output.co_expression_check",
                                                                                              tags$h5("Co-expression threshold:"),
                                                                                              sliderInput(inputId = "coex_threshold",
                                                                                                          label = NULL,
                                                                                                          value = 1,
                                                                                                          min = 0,
                                                                                                          max = 100,
                                                                                                          step = 1)
                                                                                              
                                                                                              
                                                                                              
                                                                             )
                                                            ),
                                                            
                                                            conditionalPanel(condition = "output.fc_check",
                                                                             tags$h4("Table Selector"),
                                                                             selectizeInput(inputId = "fc_table_select",
                                                                                            label = NULL,
                                                                                            multiple = T,
                                                                                            choices = list()),
                                                                             
                                                                             conditionalPanel(condition = 'output.fc_sliders',
                                                                                              tags$hr(),
                                                                                              
                                                                                              radioButtons(inputId = "fc_yes",
                                                                                                           label = "logFC Filter",
                                                                                                           choices = c("Yes", "No"),
                                                                                                           selected = "No",
                                                                                                           inline = T),
                                                                                              
                                                                                              conditionalPanel(condition = "input.fc_yes == 'Yes'",
                                                                                                               sliderInput(inputId = "logfc_threshold_1",
                                                                                                                           label = "",
                                                                                                                           value = c(0,1),
                                                                                                                           min = 10,
                                                                                                                           max = 0,
                                                                                                                           step = 1),
                                                                                                               
                                                                                                               conditionalPanel(condition = 'output.multi_fold_data',
                                                                                                                                sliderInput(inputId = "logfc_threshold_2",
                                                                                                                                            label = "",
                                                                                                                                            value = c(0,1),
                                                                                                                                            min = 0,
                                                                                                                                            max = 10,
                                                                                                                                            step = 1)
                                                                                                               )
                                                                                              ),
                                                                                              
                                                                                              radioButtons(inputId = "pval_yes",
                                                                                                           label = "PValue Filter",
                                                                                                           choices = c("Yes", "No"),
                                                                                                           selected = "No",
                                                                                                           inline = T),
                                                                                              
                                                                                              conditionalPanel(condition = "input.pval_yes == 'Yes'",
                                                                                                               sliderInput(inputId = "pval_threshold_1",
                                                                                                                           label = "",
                                                                                                                           value = 1,
                                                                                                                           min = 0,
                                                                                                                           max = 10,
                                                                                                                           step = 1),
                                                                                                               
                                                                                                               conditionalPanel(condition = 'output.multi_fold_data',
                                                                                                                                sliderInput(inputId = "pval_threshold_2",
                                                                                                                                            label = "",
                                                                                                                                            value = 1,
                                                                                                                                            min = 0,
                                                                                                                                            max = 10,
                                                                                                                                            step = 1)
                                                                                                               )
                                                                                              )
                                                                             )
                                                            )
                                                   ),
                                                   
                                                   tabPanel("Click Options",
                                                            # tags$hr(),
                                                            tags$br(),
                                                            tags$h4("Condition (Text shown)"),
                                                            selectizeInput(inputId = "label_type",
                                                                           label = NULL,
                                                                           selected = "none",
                                                                           choices = "none"),
                                                            
                                                            radioButtons(inputId = "multi_label",
                                                                         label = "Highlight",
                                                                         choices = c("only clicked point"="s", "all matching point"="m"),
                                                                         selected = "s",
                                                                         inline = T),
                                                            tags$hr()
                                                   )
                                       ),
                                     
                                       tags$br(),
                                       
                                       downloadButton(outputId="downloadData",
                                                      label = "Download")
                                     
                    ),
                    
                    mainPanel(
                      tabsetPanel(id = "plot_tab",
                                  tabPanel("PCA Scatter Plot",
                                           tags$h3("Step 3: See the Plot Options on the left"),
                                           tags$h4("Each item has a name that identifies the variable, as well as smaller text that identifies what part of the plot
                                                   the variable will control. The options will change depending on which plot is selected. Currently selected is PCA plot."),
                                           tags$hr(),
                                           tags$br(),
                                           plotOutput("pca_plot", height="100%", width = "100%", click = "plot_click_1"),
                                           tags$br(),
                                           tags$hr(),
                                           tags$h4("A PCA (Principal Component Analysis) plot is a form of dimentionality reduction used to
                                                    plot relationships between data points.")
                                  ),
                                  tabPanel("MDS Scatter Plot",
                                           tags$h3("Step 3: See the Plot Options on the left"),
                                           tags$h4("Each item has a name that identifies the variable, as well as smaller text that identifies what part of the plot
                                                   the variable will control. The options will change depending on which plot is selected. Currently selected is MDS plot."),
                                           tags$hr(),
                                           tags$br(),
                                           plotOutput("mds_plot", height="100%", width = "100%", click = "plot_click_2"),
                                           tags$br(),
                                           tags$hr(),
                                           tags$h4("An MDS (multi-dimentional scaling) plot is a form of dimentionality reduction used to
                                                          plot relationships between data points.")
                                  ),
                                  
                                  tabPanel("DGE Box Plot",
                                           tags$h3("Step 3: See the Plot Options on the left"),
                                           tags$h4("In order to activate this plot, you must first select a gene to be plotted. The plot will
                                                    tell you whether the selected gene has sufficient counts to display a boxplot.
                                                    The condition variable will dictate how many boxes are generated, one for each condition."),
                                           tags$hr(),
                                           tags$br(),
                                           plotOutput("box_plot", height="100%", width = "100%", click = "plot_click_3"),
                                           tags$br(),
                                           tags$hr(),
                                           tags$h4("A DGE boxplot measures variation of gene expression across conditions.")
                                  ),
                                  
                                  tabPanel("Fold Change",
                                           tags$br(),
                                           tabsetPanel(id = "fold_change_panel",
                                                       tabPanel("Table",
                                                                tags$h3("Step 3: See the Plot Options on the left"),
                                                                tags$h4("Select one or two tables from the Table Selector. The available tables are based on your selected datasets, generally separated by a condition from the experiment.
                                                                               Your plot controller will give you the option to adjust the thresholds for logFC and P-value for both datasets. Once you've played around with the tables,
                                                                               click on the plot tab."),
                                                                tags$hr(),
                                                                tags$br(),
                                                                conditionalPanel(condition = 'output.fc_tables_yes', tags$h4("No FC tables uploaded.")),
                                                                DT::dataTableOutput("fold_change_table")
                                                       ),
                                                       
                                                       tabPanel("Plot",
                                                                tags$h3("Step 3: See the Plot Options on the left"),
                                                                tags$h4("To activate this plot, two tables must be selected. Adjusting the thresholds for logFC and PValue will automatically add/remove those points from the table."),
                                                                tags$hr(),
                                                                tags$br(),
                                                                plotOutput("fc_plot", height="100%", width = "100%", click = "plot_click_4")
                                                       ),
                                                       
                                                       tabPanel("Heatmap",
                                                                tags$h3("Step 3: See the Plot Options on the left"),
                                                                tags$h4("Press the 'Generate Heatmap' button to generate a heatmap based on the current genes present in the Fold Change, Table display."),
                                                                tags$hr(),
                                                                tags$br(),
                                                                
                                                                actionButton(inputId = "heatmap_button",
                                                                             label = "Generate Heatmap"),
                                                                
                                                                tags$br(),
                                                                tags$br(),
                                                                
                                                                conditionalPanel(condition = 'output.heatmap_yes', plotlyOutput("heatmap", height="800px", width="900px"))
                                                       )
                                                       
                                           ),
                                           tags$br(),
                                           tags$hr(),
                                           tags$h4("A fold change analysis is also a measure of variation across conditions and is quantified with a logFC value and P-value.")
                                           
                                  )
                      )
                    )
           )
           
),

tabPanel("Update Log",
         tags$hr(),
         tags$h2("Update Log"),
         tags$hr(),
         
         tags$h4("V2.01 - Sep 15, 2017"),
         tags$h5("Aestetic update to allow for better table functionality. The table functions are imported from the DT R package. It has search, selection, and sort functionality.
                            Updated Katayama data, which had a few missing counts, and Platt data, which had missing libraries.
                            Human and mouse gene symbol names are now consistent. Genes present in both human and mouse will be registered as the same when analyzing counts.
                            We added a fold change analysis and correlation plot based on log fold change values. All visualizations are currently downloadable.
                            We also are going to implement a GO enrichment analysis soon. Stay tuned for more updates."),
         tags$br(),
         tags$h4("V1.03 - Sep 11, 2017"),
         tags$h5("We added our chd8 data to the data pool (Chd8: Nord). We also added a new condition variable called ExperimentalDesign, 
                            which separates the samples by whichever variables were being tested in the exerpiment. We optimized some of the code so the browser should run faster.
                            Lastly, we updated the browser to allow human and mouse gene symbols to be compatible. It is now possible to view human and mouse data together with counts for both."),
         tags$h5("We're also putting structures in place so eventually users will be able to upload their own data to analyze on their own and compare with the sets
                            we've already provided. A note about that, the browser will have a good deal of flexibility when it comes to the format of user inputted data, but several 
                            columns will be required for the sample table. Sample_ID, Dataset, and ExperimentalDesign will all be required. Organism will be a separate input through the interface, but will also
                            be required. Formatting of gene names in the counts matrix is not required, but certain formats of names will make genes compatible with ours."),
         tags$br(),
         tags$h4("V1.02 - Aug 30, 2017"),
         tags$h5("We added a DGE box plot feature to the browser. Now all of the available data will be able to be viewed in an expression box plot."),
         tags$br(),
         tags$h4("V1.01 - Aug 21, 2017"),
         tags$h5("Welcome to the Nordlab RNA browser. This browser features a flexible interface that allows the user to populate the site with different
                            datasets, then view those sets in a MDS scatter plot. The user can also change which condition variable is viewed on the plot. The color
                            of the dots on the plot will change to reflect the chosen variable. The user can also select a specific gene and view that gene's
                            raw/relative expression values. In this case, the condition variable changes the shape of the dots, and the expression will be shown in color.")
),

tabPanel("Explore Our Other Browsers",
         tags$hr(),
         tags$h2("Explore"),
         tags$hr(),
         
         tags$a(href="https://nordlab.shinyapps.io/base_camp/", "Click here")
),

tags$hr(),
tags$h5("V2.01, kjulim@ucdavis.edu"),

tags$head(tags$style(HTML("
                                    body { 
                                          background-image: url('./DNA.jpeg');
                                          background-repeat: no-repeat;
                                          background-position: left top; 
                                    
                                    }

                                    .progress-bar{background-color:#ffb62f; }
                                    
                                    #generate_button{ background-color:#ffb62f; }
                                    #generate_button:hover{ background-color:#56B4E9; }
                                    
                                    #select_new_data_button{ background-color:#ffb62f; }
                                    #select_new_data_button:hover{ background-color:#56B4E9; }

                                    #user_data_button{ background-color:#ffb62f; }
                                    #user_data_button:hover{ background-color:#56B4E9; }
                                    
                                    #downloadData{ background-color:#f09e21; color:#f0f0f0; }
                                    #downloadData:hover{ background-color:#56B4E9; }

                                    #heatmap_button{ background-color:#f09e21; color:#f0f0f0; }
                                    #heatmap_button:hover{ background-color:#56B4E9; }
                                    
                                    .tabbable > .nav > li[class=active] > a {border-color:#ffb62f; background-color:#ffb62f; color:white }
                                    .tabbable > .nav > li > a:hover { background-color:transparent; color:#56B4E9; border-color:transparent }
                                    .navbar .nav > li > a:hover { color:#ffb62f; }
                                    
                                    hr { border-color:#2b3e50; }
                                    a { color:#f09e21; }
                                    a:hover { color:#56B4E9; }
                                    
                                    .navbar-default .navbar-brand:hover { color:#ffb62f; }")))


)