#
#
#This is a Shiny web application.
#libraries :

library(shiny)
library(plotly)
library(stringr)
library(tidyr)
library(tibble)
library(BiocManager)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(visNetwork)
library(pheatmap)
library(DT)
library(shinycssloaders)
library(bs4Dash)
library(bslib)

#loading expression matrix
expr <- readRDS("exp.mat.rds")


ui <- dashboardPage(
    dark = FALSE,
    header = dashboardHeader(
        title = dashboardBrand("Brain Map : Mapping Brain Genome", color = "lightblue"),
        skin = "light"
    ),
    controlbar = NULL,
    sidebar = dashboardSidebar(disable = TRUE) ,
    body = dashboardBody(
        
        # Custom tab styles
        # CSS to center the input box
        tags$head(tags$style(HTML("
      .center-input {
        display: flex;
        justify-content: center;
        margin-bottom: 20px;
      }
      .shiny-input-container {
        width: 50%;
      }
    "))),
        
        br(),
        
        
        div(class = "center-input",
            textInput("genes", "Enter gene symbols (comma-separated):", value = "CDK1 , CDC20 , CCNB1,BUB1,TOP2A")
        ),
        tabBox(
            width = 12,
            tabPanel("Expression Map",
                     fluidRow(
                     column(12, helpText("This tab shows expression levels of selected genes across brain tissues ")) ,    
                     br(),
                     br(),
                     br(),
                     HTML("<h5 style='text-align:center; color:#333;'>Bar Plot</h5>"),
                     column(12, withSpinner(plotlyOutput("barplot", height = "400px"))),
                     br(),
                     br(),
                     br(),
                     br(),
                     br(),
                     HTML("<h5 style='text-align:center; color:#333;'>Heatmap</h5>"),
                     column(12,withSpinner(plotlyOutput("heatmap", height = "500px"))),
                     br(),
                     HTML("<h5 style='text-align:center; color:#333;'>Missing Genes</h5>"),
                     tableOutput("missingGenes")
            )
            ),
            tabPanel("Gene Co-expression",
                     fluidRow(
                         column(12 , helpText("co expression heatmap of corolated genes")),
                         br(),
                         column (12,withSpinner(plotOutput("cor_heatmap" , height = "500px")))
                         )
                     ),
            
            tabPanel("Enrichment Analysis",
                      fluidRow(
                    column(12, helpText("KEGG pathway enrichment analysis for input genes")),
                    br(),
                    column(12,withSpinner(plotOutput("pathway_plot"))),
                    br(),
                    br(),
                    column(12 ,withSpinner(tableOutput("pathway_table"))))
                     
            
            ),
            tabPanel("Explore Gene Resources",
                     fluidRow(
                         column(12,withSpinner(dataTableOutput("gene_links")))
                     )
                     ),
            tabPanel("About",
                     fluidRow(
                         column(10, offset = 1,
                                h5("Overview"),
                                p("Brain Map is a shiny app for exploring gene expression,
                                  This app focuses on five behaviorally and clinically important brain regions — the hypothalamus,hippocampus, amygdala, frontal cortex, and anterior brain areas — to provide a targeted and interpretable view of gene expression relevant to mood, cognition, and neuroendocrine regulation"),
                                h5("Features"),
                                tags$ul(
                                    tags$li("Expression Map : Explore gene expression patterns in brian tissues, Expression value is illustrated as  bar-plot and heatmap of the log transformed data"),
                                    tags$li("Gene co- expression: Visualise correlation of gene expression "),
                                    tags$li("Pathway enrichment : Perform KEGG enrichment on selected genes "),
                                    tags$li("Gene resources : Access NCBI , Ensembl , RefEX and UniProt gene info "),
                                ),
                                h5("Data Source"),
                                p("The expression matrix used is log-transformed RNA-seq data from brain tissue. Gene identifiers are HGNC symbols. The gene expression data was derived from GTEx (Genotype-Tissue Expression) protal."),
                               
                                h5("Github repository"),
                                HTML('<p>Full source code plus more info available at:<br>
                                     <a href="https://github.com/Yasna81/Brain_map.git"target="_blank">Github_repo</a></p>')
                                
                         )
                     )
                         
       )
       
        ),
    ),
    footer = dashboardFooter(
        left = "Brain Map app 2025"
    ),
    title = "Brain Map : Mapping Brain Genome"
)


                     


    



server <- function(input, output, session) {
    
    selected_genes <- reactive({
        req(input$genes)
        input$genes %>%
            str_split(",") %>% .[[1]] %>%
            str_trim() %>%
            toupper()
    })
    
    output$missingGenes <- renderTable({
        genes <- selected_genes()
        missing <- genes[!genes %in% rownames(expr)]
        if (length(missing) > 0) {
            data.frame(Missing_Genes = missing)
        } else NULL
    })
    
    output$barplot <- renderPlotly({
        genes <- selected_genes()
        found_genes <- genes[genes %in% rownames(expr)]
        req(length(found_genes) > 0)
        
        df <- expr[found_genes, , drop = FALSE] %>%
            rownames_to_column("Gene") %>%
            pivot_longer(-Gene, names_to = "Region", values_to = "Expression")
        
        plot_ly(df, x = ~Region, y = ~Expression, color = ~Gene, type = "bar") %>%
            layout(barmode = "group")
    })
    
    output$heatmap <- renderPlotly({
        genes <- selected_genes()
        found_genes <- genes[genes %in% rownames(expr)]
        req(length(found_genes) > 0)
        
        mat <- expr[found_genes, , drop = FALSE]
        plot_ly(z = as.matrix(mat), x = colnames(mat), y = rownames(mat),
                type = "heatmap", colorscale = "Viridis")
    })
    
    output$cor_heatmap <- renderPlot({
        req(selected_genes())
        genes <- selected_genes()
        found_genes <- genes[genes %in% rownames(expr)]
        req(length(found_genes) > 1)
        
        sub_exp <- expr[rownames(expr) %in% found_genes, , drop = FALSE]
        if (nrow(sub_exp) < 2) return(NULL)
        
        sub_exp <- as.matrix(sub_exp)
        storage.mode(sub_exp) <- "numeric"
        
        cor_mat <- cor(t(sub_exp), use = "pairwise.complete.obs")
        
        tryCatch({
               pheatmap(cor_mat, main = "Gene-Gene Correlation Heatmap")
        }, error = function(e) {
            plot.new()
            title("Error generating heatmap")
        })
    })
    # --- clusterProfiler pathway enrichment ---
    output$pathway_table <- renderTable({
        genes <- selected_genes()
        req(length(genes) > 0)
        
        entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
        if (nrow(entrez) < 2) return(data.frame(Message = "No valid Entrez IDs found."))
        
        kegg <- enrichKEGG(gene = entrez$ENTREZID, organism = 'hsa', pAdjustMethod = "BH", qvalueCutoff = 0.05)
        
        if (is.null(kegg) || nrow(as.data.frame(kegg)) == 0) {
            return(data.frame(Message = "No enriched pathways found."))
        }
        
        as.data.frame(kegg)[, c("Description", "p.adjust", "Count")]
        
    })
    
    output$pathway_plot <- renderPlot({
        genes <- selected_genes()
        req(length(genes) > 0)
        entrez <- bitr(genes, fromType = "SYMBOL",toType = "ENTREZID",OrgDb  = org.Hs.eg.db)
        if (nrow(entrez) < 2) return()
        kegg <-enrichKEGG(gene = entrez$ENTREZID , organism = "hsa",pAdjustMethod = "BH" , qvalueCutoff = 0.05)
        if (is.null(kegg) || nrow(as.data.frame(kegg)) == 0) return()
        dotplot(kegg, showCategory = 10, title = "KEGG Pathway Enrichment")
    })
    
    output$gene_links <- renderDataTable({
        genes <- selected_genes()
        req(length(genes) > 0)
        df_links <- data.frame(
            Gene = genes ,
            NCBI = paste0('<a href="https://www.ncbi.nlm.nih.gov/gene/?term=', genes,'" target="_blank">View in NCBI</a>'),
            Ensembl = paste0('<a href ="https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=',genes,'" target="_blank">View in Ensembl</a>'),
            RefEX = paste0('<a href="https://refex.dbcls.jp/refex_human_gene.php?gene=', 
      genes, '" target="_blank">View in RefEx</a>'),
            UniProt = paste0('<a href="https://www.uniprot.org/uniprot/?query=', genes,'" target="_blank">View in UniProt</a>'),
            stringsAsFactors = FALSE
        )
        datatable(df_links, escape = FALSE, rownames = FALSE , options = list(pageLength = 10)) %>%
            formatStyle(columns = c("Gene","NCBI","Ensembl","RefEX" , "UniProt"),
                        target = 'cell',
                        fontWeight = 'bold',
                        color = 'blue')
    })
    
    observeEvent(input$toggle_dark,{
        toggleDarkMode()
    })
    
}

shinyApp(ui, server)



