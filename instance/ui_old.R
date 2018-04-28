navbarPage("Nordlab: RNA Browser", theme = shinytheme("flatly"), selected = "Plot",
           
           tabPanel("Plot",
                    conditionalPanel(condition = 'output.data_check',
                                     tags$h1("Welcome!"),
                                     tags$h4("This is the Nordlab RNA Browser. It is an interactive tool for viewing rna-seq data across multiple datasets. Follow the steps below and happy browsing."),
                                     tags$br(),
                                     tags$hr(),
 
                                     tags$h2("Step 1: Choose data"),
                                     tags$h4("Select one or more datasets from the box below and choose how you would like to view your counts (Pseudo counts will provide more functions), then click Generate."),
                                     tags$hr(),
                                     
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
                                     tags$br(),
                                     
                                     actionButton(inputId = "generate_button",
                                                  label = "Generate")
                    ),
                    
                    conditionalPanel(condition = 'output.data_button',
                                     tags$hr(),
                                     actionButton(inputId = "select_new_data_button",
                                                  label = "Select New Data")
                    ),
                    
                    conditionalPanel(condition = "output.MDS_display",
                                     tags$hr(),
                                     tags$h2("Step 2: Check out the plot tabs below"),
                                     tags$h4("Each tab will give you further instructions. All plots are interactive."),
                                     tags$hr(),
                                     
                                     sidebarPanel(
                                       tags$head(
                                         tags$style(type="text/css", "select { width: 350px; }"),
                                         tags$style(type="text/css", ".span4 { width: 350px; }"),
                                         tags$style(type="text/css", ".well { width: 350px; }")),
                                       
                                       tags$h3("- Plot Controller -"),
                                       tags$hr(),
                                       
                                       conditionalPanel(condition = 'output.not_fc_check',
                                                         tags$h4("Condition selector"),
                                                         selectizeInput(inputId = "condition",
                                                                        label = NULL,
                                                                        selected = "ExperimentalDesign",
                                                                        choices = "ExperimentalDesign"),
                                                         
                                                         tags$hr(),
                                                         
                                                         tags$h4("Gene selector"),
                                                         selectizeInput(inputId = "gene",
                                                                        label = NULL,
                                                                        multiple = T,
                                                                        choices = list()),
                                                         
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
                                                                                                      min = 0,
                                                                                                      max = 10,
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
                                       ),
                                       
                                       # Implemented, but not necessary right now
                                       conditionalPanel(condition = "input.condition == 'fancy'",
                                                        tags$h4("Step 5: Change your label category"),
                                                        selectizeInput(inputId = "label_type",
                                                                       label = NULL,
                                                                       selected = "none",
                                                                       choices = "none"),

                                                        radioButtons(inputId = "multi_label",
                                                                     label = NULL,
                                                                     choices = c("only selected"="s", "all that match"="m"),
                                                                     selected = "s",
                                                                     inline = T),
                                                        tags$hr()
                                       ),
                                       
                                       tags$br(),
                                       
                                       downloadButton(outputId="downloadData",
                                                      label = "Download")
                                       
                                     ),
                                     
                                     mainPanel(
                                       tabsetPanel(id = "plot_tab",
                                         tabPanel("MDS Scatter Plot",
                                                  #tags$br(),
                                                  tags$h3("Step 3: Check out the Plot Controller on the left"),
                                                  tags$h4("The Condition Selector will change the color 
                                                          of the points depending on the chosen condition. The Gene Selector will move the condition denotion
                                                          to shape, and will color the points based on gene expression (we suggest Chd8 as an example gene).
                                                          ExperimentalDesign is the default condition."),
                                                  tags$hr(),
                                                  tags$br(),
                                                  plotOutput("interactive_plot", height="100%", width = "100%", click = "plot_click"),
                                                  tags$br(),
                                                  tags$hr(),
                                                  tags$h4("An MDS (multi-dimentional scaling) plot is a form of dimentionality reduction used to
                                                          plot relationships between data points.")
                                         ),
                                         
                                         tabPanel("DGE Box Plot",
                                                  #tags$br(),
                                                  tags$h3("Step 3: Check out the Plot Controller on the left"),
                                                  tags$h4("In order to activate this plot, you must first select a gene to be plotted. The plot will
                                                          tell you whether the selected gene has sufficient counts to display a boxplot.
                                                          The Condition Selector will dictate how many boxes are generated, one for each condition.
                                                          Again, ExperimentalDesign is the default condition."),
                                                  tags$hr(),
                                                  tags$br(),
                                                  plotOutput("box_plot", height="100%", width = "100%", click = "plot_click_2"),
                                                  tags$br(),
                                                  tags$hr(),
                                                  tags$h4("A DGE boxplot measures variation of gene expression across conditions.")
                                         ),
                                         
                                         tabPanel("Fold Change",
                                                  tags$br(),
                                                  tabsetPanel(id = "fold_change_panel",
                                                              tabPanel("Table",
                                                                       tags$h3("Step 3: Check out the Plot Controller on the left"),
                                                                       tags$h4("Select one or two tables from the Table Selector. The available tables are based on your selected datasets, generally separated by a condition from the experiment.
                                                                               Your plot controller will give you the option to adjust the thresholds for logFC and P-value for both datasets. Once you've played around with the tables,
                                                                               click on the plot tab."),
                                                                       tags$hr(),
                                                                       tags$br(),
                                                                       DT::dataTableOutput("fold_change_table")
                                                              ),
                                                              
                                                              tabPanel("Plot",
                                                                       tags$h3("Step 3: Check out the Plot Controller on the left"),
                                                                       tags$h4("To activate this plot, two tables must be selected. Adjusting the thresholds for logFC and PValue will automatically add/remove those points from the table."),
                                                                       tags$hr(),
                                                                       tags$br(),
                                                                       plotOutput("fc_plot", height="100%", width = "100%", click = "plot_click_3")
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
                            columns will be required for the sample table. Lib_ID, Dataset, and ExperimentalDesign will all be required. Organism will be a separate input through the interface, but will also
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
                                    
                                    #downloadData{ background-color:#f09e21; color:#f0f0f0; }
                                    #downloadData:hover{ background-color:#56B4E9; }
                                    
                                    .tabbable > .nav > li[class=active] > a {border-color:#ffb62f; background-color:#ffb62f; color:white }
                                    .tabbable > .nav > li > a:hover { background-color:transparent; color:#56B4E9; border-color:transparent }
                                    .navbar .nav > li > a:hover { color:#ffb62f; }
                                    
                                    hr { border-color:#2b3e50; }
                                    a { color:#f09e21; }
                                    a:hover { color:#56B4E9; }
                                    
                                    .navbar-default .navbar-brand:hover { color:#ffb62f; }")))
           
           
)