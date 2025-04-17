library(shiny)
library(biomaRt)
library(Gviz)
library(DT)
library(dplyr)
library(shinycssloaders)

ui <- fillPage(
  titlePanel("Gene Isoform Explorer", windowTitle = "Gene Isoform Visualization"),
  # Input: gene and load button
  fluidRow(
    column(12, align = "center",
           div(style = "max-width: 400px; margin: 0 auto;",
               textInput("gene", "Enter Gene Name:", value = "TP53", width = "100%")
           ),
           actionButton("load", "Load Data", class = "btn-primary",
                        icon = icon("dna"))
    )
  ),
  tabsetPanel(
    # First tab: visualisation of isoforms
    tabPanel("Gene Visualization",  
             div(
               style = "height: 700px; overflow-y: auto; padding: 10px;",
               withSpinner(
                 plotOutput("plot", height = "auto", width = "100%"),
                 type = 6, color = "#4a90e2"
               )
             )
    ),
    # Second tab: transcript/exon data and download button
    tabPanel("Transcript Data",  
             div(
               downloadButton("downloadData", "Download Exon Data",
                              class = "btn-success",
                              icon = icon("download")),
               div(style = "margin-bottom: 12px;"),
               withSpinner(DTOutput("table"), type = 6, color = "#4a90e2")
             )
    )
  )
)



server <- function(input, output, session) {
  # Loading data in
  isoform_data <- eventReactive(input$load, {
    req(input$gene)
    
    showNotification("Fetching data from Ensembl...", type = "message", duration = 5)
    
    
    tryCatch({
      mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      
      transcripts <- getBM(
        attributes = c(
          "ensembl_transcript_id", "external_gene_name",
          "chromosome_name", "strand",
          "exon_chrom_start", "exon_chrom_end",
          "transcript_start", "transcript_end",
          "rank"
        ),
        filters = "external_gene_name",
        values = input$gene,
        mart = mart
      )
      
      # Show error if nothing found
      validate(need(nrow(transcripts) > 0, "No transcripts found for this gene symbol."))
      
      # Creating final version of data table
      transcripts %>%
        mutate(
          chromosome = paste0("chr", chromosome_name),
          strand = ifelse(strand == 1, "+", "-"),
          exon_length = exon_chrom_end - exon_chrom_start
        ) %>%
        select(
          ensembl_transcript_id,
          external_gene_name,
          chromosome,
          strand,
          exon_rank = rank,
          start = exon_chrom_start,
          end = exon_chrom_end,
          exon_length
        ) %>%
        arrange(ensembl_transcript_id, exon_rank)
      
      # If error, notifies and returns an empty df
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      return(data.frame(
        ensembl_transcript_id = character(),
        external_gene_name = character(),
        chromosome = character(),
        strand = character(),
        exon_rank = numeric(),
        start = numeric(),
        end = numeric(),
        exon_length = numeric()
      ))
    })
  })
  
  # Download data as .csv
  output$downloadData <- downloadHandler(
    filename = function() paste(input$gene, "_exons.csv", sep = ""),
    content = function(file) {
      req(nrow(isoform_data()) > 0)
      write.csv(isoform_data(), file, row.names = FALSE)
    }
  )
  
  # Plotting gene isoforms
  output$plot <- renderPlot({
    df <- isoform_data()
    req(nrow(df) > 0)
    
    chr <- unique(df$chromosome)
    if (length(chr) != 1) {
      showNotification("Multiple chromosomes detected - showing first one", type = "warning")
      chr <- chr[1]
    }
    
    gene_start <- min(df$start)
    gene_end <- max(df$end)
    range_width <- gene_end - gene_start
    
    genome_axis <- GenomeAxisTrack(
      fontsize = 14,
      fontcolor = "black",
      labelPos = "above"
    )
    
    gene_track <- GeneRegionTrack(
      df,
      chromosome = chr,
      gene = df$external_gene_name,
      transcript = df$ensembl_transcript_id,
      start = df$start,
      end = df$end,
      strand = df$strand,
      name = "Transcripts",
      showId = TRUE,
      transcriptAnnotation = "transcript",
      fontsize.group = 18,
      col.title = "black"
    )
    
    plotTracks(
      list(genome_axis, gene_track),
      from = gene_start - round(range_width * 0.1),
      to = gene_end + round(range_width * 0.02),
      main = paste("Isoforms of", input$gene),
      cex.main = 1.8
    )
  }, height = function() {
    df <- isoform_data()
    if (nrow(df) == 0) return(400)
    max(400, length(unique(df$ensembl_transcript_id)) * 60)
  })
  
  output$table <- renderDT({
    df <- isoform_data()
    req(nrow(df) > 0)
    
    if ("strand" %in% names(df)) {
      df <- df %>% select(-strand)
    }
    
    datatable(
      df,
      colnames = c("Transcript ID", "Gene Name", "Chromosome", 
                   "Exon Rank", "Start", "End", "Exon Length"),
      options = list(
        scrollX = TRUE,
        pageLength = 10
      ),
      rownames = FALSE
    )
  })
  
}

shinyApp(ui = ui, server = server)
