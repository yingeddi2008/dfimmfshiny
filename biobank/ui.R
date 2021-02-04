library(shiny)
library(shinythemes)

navbarPage("DFI Commensal Library",
           
 # 0 page: welcome ---------------------------------------------------------
           tabPanel("About", 
                    fluidPage(theme = shinytheme("sandstone")),
                    h6("Wiki link:",a("DFI wiki",href="http://128.135.41.103/wiki/index.php")),
                    h3("Mission"),
                    p("The DFI aims to maximize good health and the economic, social, and personal benefits it delivers."),
                    p(""),
                    hr(""),
                    code("v4.6")
           ),
           
  # 1 page: summary of all isolates -----------------------------------------
  navbarMenu("Isolates",
             tabPanel("List",
                      h4("All available isolates in biobank:"),
                      p("You can filter or search for isolate of interet by taxonomy."),
                      DT::dataTableOutput("bug_list")
                      ),
             tabPanel("Phylogenetic Tree",
                      sidebarLayout(
                        sidebarPanel(width = 3,
                            checkboxGroupInput("phyla",
                                             "Select Phyla:",
                                             choices = c("Firmicutes","Actinobacteria",
                                                         "Bacteroidetes","Proteobacteria"),
                                             selected= c("Firmicutes","Actinobacteria",
                                                         "Bacteroidetes","Proteobacteria")),
                          h5(strong("Show:")),
                          checkboxInput("msk_id","msk_id",F),
                          checkboxInput("species","species",F)
                        ),
                        mainPanel(
                          plotOutput("tree",height="1200px",width="800px")
                        )
                      )
                      
                      )
             ),
  
  # 2 page: Taxonomy --------------------------------------------------------
             tabPanel("Taxonomy",
                      sidebarLayout(
                        sidebarPanel(width = 3,
                                     textInput("kraken_mskid",
                                               "Enter msk_id(s):",
                                               value="MSK.1.10")
                                     ),
                        mainPanel(plotOutput("kraken2"))
                       )
                      ),
  
  # 3 page: blast -----------------------------------------------------------
  navbarMenu("Blast",
             tabPanel("Biobank proteins",
                      h4("You are about to compare your own sequences to all Biobank proteins."),
                      selectInput("program", "Program:", choices=c("blastn","blastp","blastx"), width="100px",
                                  selected="blastp"),
                      textAreaInput('queryProt', 'Input sequence:', 
                                    value = "", 
                                    placeholder = "", 
                                    width = "600px", height="200px"),
                      DT::dataTableOutput("blastProtResults"),
                      actionButton("blastProt", "Run BLAST"),
                      downloadButton("downloadProtData", "Download")
                      ),
             tabPanel("Biobank 16s rRNA",
                      h4("You are about to compare your own sequences to all Biobank 16s rRNA sequences"),
                      textAreaInput('query16s', 'Input 16s rRNA Sanger sequence:', 
                                    value = "", 
                                    placeholder = "", 
                                    width = "600px", height="200px"),
                      DT::dataTableOutput("blast16sResults"),
                      actionButton("blast16s", "Run BLASTn"),
                      downloadButton("download16sData", "Download")
                      )
             ),
  
  # 4 page: metabolomics ----------------------------------------------------
             tabPanel("Metabolomics",
                      plotOutput("metabolomics",height="500px",width="1200px"),
                      
                      checkboxInput("metamsk_id","show msk_id",F)
                      
                     ),
  
  # more --------------------------------------------------------------------
             tabPanel("More",
                      h4("Under construction")
                      
                      )
  
             
)



