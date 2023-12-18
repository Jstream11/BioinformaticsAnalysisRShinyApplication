library(shiny)
library(DT)
library(ggplot2)
library(DESeq2)
library(knitr)
library(gridExtra)


source("main.R")
# UI
ui <- fluidPage(
  titlePanel("Exploratory Data Analysis"),
  
  # Tabset
  tabsetPanel(
    
    # Sample Information Exploration
    tabPanel("Sample Information Exploration",
              sidebarLayout(
                sidebarPanel(
                  fileInput("sample_info", "Upload Sample Information Matrix (CSV)", accept = "csv", multiple = FALSE),
                  width = 2
                ),
               mainPanel(#creating subtabs
                 tabsetPanel(
                   tabPanel("Summary Table",
                            htmlOutput("sample_summary")),
                   tabPanel("Metadata Table",
                            dataTableOutput("metadata_table")),
                   tabPanel("Graphs",
                            sidebarLayout(
                              sidebarPanel(
                                selectInput("plot_type", "Choose Plot Type", c("Histogram", "Violin", "Density")),
                                radioButtons("plot_column", "Choose Column to Plot", ""),
                                radioButtons("group_by_column", "Choose Column to Group By (Optional)", "")),
                              mainPanel(plotOutput("graph_output")))
                            )
                 )
          )
        )
    
    ),
    
    # Counts Matrix Exploration
    tabPanel("Counts Matrix Exploration",
            sidebarLayout(
              sidebarPanel(
                fileInput("counts_matrix", "Upload Normalized Counts Matrix (CSV)", multiple = FALSE),
                width = 2
              ),
              mainPanel(
                #create subtabs
                tabsetPanel(tabPanel("Counts Data Table", 
                                  dataTableOutput("countsdata_table")),
                        tabPanel("Counts Summary", 
                                 sidebarLayout(
                                   sidebarPanel(
                                     sliderInput("variance_percentile", "Variance Percentile Filter", min = 0, max = 100, value = 50),
                                     sliderInput("non_zero_threshold", "Non-Zero Samples Filter", min = 1, max = 100, value = 10)),
                                   mainPanel(
                                     dataTableOutput("counts_summary_table")))
                                 ),
                        tabPanel("Scatter Plots", plotOutput("scatters")),
                        tabPanel("Heatmap", plotOutput("heatmap")),
                        tabPanel("PCA",
                                 sidebarLayout(
                                   sidebarPanel(
                                     sliderInput("principle_component_x","Principle component(x) to plot: ", min=1, max=100, value =1),
                                     sliderInput("principle_component_y","Principle component(y) to plot: ", min=1, max=100, value =1)),
                                   mainPanel(
                                     plotOutput("pca_plot")))
                                 )
                        )
              )
      )
    ),
    
    # Differential Expression
    tabPanel("Differential Expression",
             sidebarLayout(
               sidebarPanel(
                 fileInput("diff_expr_results", "Upload Differential Expression Results (CSV)", multiple = FALSE, accept = ".csv"),
                 width = 2
                 ),
               mainPanel(
                 tabsetPanel(
                   tabPanel(
                     "Data Table", dataTableOutput("DE_data_table")),
                   tabPanel("Plots",
                     sidebarLayout(
                       sidebarPanel(
                         radioButtons("de_plot_choice", "Choose Plot", choices = c("raw p-value histogram", "log2fc histogram", "volcano plot")),
                         sliderInput("padj_threshold", "Choose padj threshold", min=0, max=300, value= 100)),
                       mainPanel(
                         plotOutput("DE_plots")))
                     )
                 )
               )
             )
            ),
    
    # GSEA
    tabPanel("GSEA",
             sidebarLayout(
               sidebarPanel(
                 fileInput("fgsea_results", "Upload FGSEA Results (CSV)", multiple = FALSE, accept = ".csv"),
                 width = 2
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel(
                     "Top Pathways", 
                     sidebarLayout(
                       sidebarPanel(
                         sliderInput("pathways_padj_threshold", "Choose padj threshold", min=0, max=1, value= 0, step=0.005),
                         width=3),
                       mainPanel(
                         plotOutput("top_pathways_plot")))),
                   tabPanel(
                     "Data Table",
                     sidebarLayout(
                       sidebarPanel(
                         radioButtons("NES_status", "Choose NES Status", choices = c("all", "positive", "negative")),
                         sliderInput("gsea_dt_padj_threshold", "Choose padj threshold", min=0, max=1, value= 1, step = 0.005),
                         width = 3,
                         downloadButton("gsea_dt_download", "Download")),
                       mainPanel(
                         dataTableOutput("GSEA_table")))),
                   tabPanel("Plots",
                            sidebarLayout(
                              sidebarPanel(
                                sliderInput("gsea_scatter_padj_threshold", "Choose padj threshold", min=0, max=1, value= 1, step = 0.005),
                                width = 3),
                              mainPanel(
                                plotOutput("GSEA_scatter"))))
                 )
               )
             )
    )
  )
)

# Server
server <- function(input, output, session) {
  
# Sample Information Exploration
  # Load and preprocess data
  meta_data <- reactive({
    #check if the button 
    req(input$sample_info)
    sample_info <- read.csv(input$sample_info$datapath)
    return(sample_info)
  })
  
  output$sample_summary <- renderTable({
    #check is the dataset is available
    req(meta_data())
    
    #get the column names and data types
    column_names <- names(meta_data())
    column_types <- sapply(meta_data(), class)
    
    #create summary table
    summary_table <- data.frame(
      Column_Name = column_names,
      Type = column_types,
      Means_SD_or_Distinct_Values = sapply(meta_data(), function(col) {
        if(is.numeric(col)){
          paste(round(mean(col), 2), " (", round(sd(col), 2), ")")
        } else{
          paste(length(unique(col)), " distinct values")
        }
      })
    )
    #print summary table
    return(summary_table)
  })
  
  #Reactive function to generate the metadata table
  output$metadata_table <- renderDataTable({
    # check if dataset is available
    req(meta_data())
    metadata_table <- meta_data()
    return(metadata_table)
  })

  #Dynamically generate column options
  observe({
    req(meta_data())
    columns <- names(meta_data())
    
    #update choices for selectInput
    updateRadioButtons(session, "plot_column", choices = columns)
    updateRadioButtons(session, "group_by_column", choices = c("No Grouping", columns))
  })

  # create plot
  output$graph_output <- renderPlot( {
    # check if dataset and plot selections are available
    req(meta_data(), input$plot_column, input$plot_type)
    
    #check if user has selected a column to plot
    if(!is.null(input$plot_column)){
      #set plot data
      plot_data <- meta_data()
      #create the appropriate plot based on the input
      if (input$group_by_column != "No Grouping"){
        #grouped plot
        ggplot(plot_data, aes_string(x = input$group_by_column, y=input$plot_column)) +
          switch(input$plot_type,
                 "Histogram" = geom_histogram(position = "identity", binwidth = 1, fill = "blue", color = "black", alpha = 0.7),
                 "Violin" = geom_violin(),
                 "Density" = geom_density()) +
          labs(title = paste("Distribution of", input$plot_column, 
                             ifelse(input$group_by_column != "No Grouping", paste("Grouped by", input$group_by_column), "")))
      } else{
        #simple plot(no grouping)
        ggplot(plot_data, aes_string(x = input$plot_column)) +
          switch(input$plot_type,
                 "Histogram" = geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7),
                 "Violin" = geom_violin(),
                 "Density" = geom_density()) +
          labs(title = paste("Distribution of", input$plot_column))
      }
    } else {
      ggplot() + theme_minimal()
    }
  })
  
  
  # Counts Matrix Exploration
  # Load and preprocess data
  counts_data <- reactive({
    #check if the button 
    req(input$counts_matrix)
    counts_matrix <- read.csv(input$counts_matrix$datapath)
    return(counts_matrix)
  })
  
  output$countsdata_table <- renderDataTable({
    return(counts_data())
  })
  
  filtering <- reactive({
    var_per_temp = 0.2
    non_zero_temp = 2
    filtered_counts <- counts_filter(counts_data(), var_per_temp, non_zero_temp)
    return(filtered_counts)
  })
  
  output$counts_summary_table <- renderDataTable({
    
    # Calculate summary statistics
    total_samples <- ncol(counts_data())
    total_genes <- nrow(counts_data())
    passed_filter <- nrow(filtering())
    removed_genes <- total_genes - passed_filter
    
    # Calculate percentages
    percent_passed_filter <- (passed_filter / total_genes) * 100
    percent_removed_genes <- (removed_genes / total_genes) * 100
    
    # Create summary table
    summary_table <- data.frame(
      "Total number of samples" = total_samples,
      "Total number of genes" = total_genes,
      "% of genes passed filter" = percent_passed_filter,
      "% of genes removed" = percent_removed_genes,
      check.names = FALSE
    )
    
    #print summary table
    return(summary_table)
  })
  
  output$scatters <- renderPlot({
    plot_var <- scatter_plot_var(filtering(), counts_data())
    plot_zeros <- scatter_plot_zeros(filtering(), counts_data())
    combined_plot <- grid.arrange(plot_var, plot_zeros, ncol = 2)
    
    return(combined_plot)
  })
  
  output$heatmap <- renderPlot({
    clustered_heatmap(filtering())
  })
  
  output$pca_plot <- renderPlot({
    plot_pca(filtering(), input$principle_component_x, input$principle_component_y)
  })
  
  # Differential Expression
  # Load and preprocess data
  DE_data <- reactive({
    #check if file uploaded 
    req(input$diff_expr_results)
    DE_matrix <- read.csv(input$diff_expr_results$datapath)
    return(DE_matrix)
  })
  
  output$DE_data_table <- renderDataTable({
    return(DE_data())
  })
  
  output$DE_plots <- renderPlot(
    #check if user has selected a plot
    if(!is.null(input$de_plot_choice)){
      #set plot data
      plot_data <- DE_data()
      switch(input$de_plot_choice,
             "raw p-value histogram" = plot_pvals(plot_data),
             "log2fc histogram" = plot_log2fc(plot_data, input$padj_threshold),
             "volcano plot" = plot_volcano(plot_data, input$padj_threshold)
      )
    }
  )
  
  #GSEA
  #Load and preprocess data
  fgsea_data <- reactive({
    #check if file uploaded
    req(input$fgsea_results)
    fgsea_matrix <- read.csv(input$fgsea_results$datapath)
    return(fgsea_matrix)
  })
  
  output$top_pathways_plot <- renderPlot(
    return(top_pathways(fgsea_data(), input$pathways_padj_threshold))
  )
  
  output$GSEA_table <- renderDataTable({
    filtered_pathways <- pathway_filter(fgsea_data(), input$gsea_dt_padj_threshold, input$NES_status)
    return(filtered_pathways)
  })
  
  output$gsea_dt_download <- downloadHandler(
    filename = function() {
      # Specify the filename for the downloaded file
      paste("filtered_gsea_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # Write the data to a CSV file
      filtered_pathways <- pathway_filter(fgsea_data(), input$gsea_dt_padj_threshold, input$NES_status)
      write.csv(filtered_pathways, file)
    }
  )
  
  output$GSEA_scatter <- renderPlot({
    plot_data <- fgsea_data()
    return(gsea_scatter(plot_data, input$gsea_scatter_padj_threshold))
  })
  
}

# Run the application
shinyApp(ui, server)
