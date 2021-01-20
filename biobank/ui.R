library(shiny)
library(shinythemes)

ui <- fluidPage(
  theme = shinytheme("yeti"),
  titlePanel("DFI Commensal Library v4.6"),
  
  mainPanel(uiOutput("url"),
            DT::dataTableOutput("bug_list"),
            tabsetPanel(type="tabs",
                        tabPanel("Phylogenetic Tree",
                                 div(style="display: inline-block;vertical-align:top; width=100px;",
                                     checkboxGroupInput("phyla","Phyla:",
                                                        inline=T,
                                                        c("Firmicutes","Actinobacteria","Bacteroidetes","Proteobacteria"),
                                                        selected= c("Firmicutes","Actinobacteria","Bacteroidetes","Proteobacteria"))),
                                 div(style="display: inline:block;vertical-align:top; width=100px;",
                                     checkboxInput("msk_id","show msk_id",F)),
                                 div(style="display: inline:block;vertical-align:top; width=100px;",
                                     checkboxInput("species","show species",F)),
                                 plotOutput("tree",height="1200px",width="1200px")),
                        tabPanel("BLAST biobank protein",
                                 selectInput("program", "Program:", choices=c("blastn","blastp","blastx"), width="100px",
                                             selected="blastp"),
                                 textAreaInput('queryProt', 'Input sequence:', 
                                               value = "", 
                                               placeholder = "", 
                                               width = "600px", height="200px"),
                                 DT::dataTableOutput("blastProtResults"),
                                 actionButton("blastProt", "Run BLAST"),
                                 downloadButton("downloadProtData", "Download")),
                        tabPanel("BLAST biobank 16s rRNA",
                                 textAreaInput('query16s', 'Input 16s rRNA Sanger sequence:', 
                                               value = "", 
                                               placeholder = "", 
                                               width = "600px", height="200px"),
                                 actionButton("blast16s", "Run BLASTn"),
                                 DT::dataTableOutput("blast16sResults"),
                                 downloadButton("download16sData", "Download")),
                        tabPanel("kraken2 contig contam",
                                 textInput("kraken_mskid","Enter msk_id:",value="MSK.1.10"),
                                 plotOutput("kraken2")),
                        tabPanel("Metabolomics",
                                 checkboxInput("msk_id","show msk_id",F),
                                 plotOutput("metabolomics",height="500px",width="1200px"))
            )
  )
  
)