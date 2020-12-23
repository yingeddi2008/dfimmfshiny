rm(list=ls())
library(tidyverse)
library(reshape2)
library(shiny)
library(DT)
library(shinythemes)
library(broom)
library(data.table)
library(yingtools2)
library(scales)
#setwd("/Volumes/chaubard-lab/shiny_workspace/csvs/")
# functions ---------------------------------------------------------------

twoFold <- function(startCon, 
                    compound = c("Butyrate", "Propionate","Acetate"), 
                    series = 8, 
                    prefix = "CC") {
  require(assertive.types)
  require(assertive.base)
  require(tidyverse)
  # arg match
  #compound <- match.arg(compound)
  # checks
  assert_is_numeric(startCon)
  startCon <- use_first(startCon)
  # body
  conc <- NULL
  conc[1] <- cur <- startCon
  i = 2
  while( i <= series){ 
    conc[i] <- cur/2
    cur <- conc[i]
    i <- i + 1
  }
  # combine to tibble and return
  odf <- tibble(compound, conc, curveLab = paste0(prefix,c(series:1)))
  return(odf)
}

get_all_conc <- function(compounds,starting_cons,series=8){
  
  #for testing:
  # compounds=vec
  # starting_cons=100,10,1
  # series=8
  
  
  #convert to character or numeric vector:
  compounds = unlist(strsplit(compounds, split=","))
  starting_cons = as.numeric(unlist(strsplit(as.character(starting_cons), split=",")))
  
  if(length(starting_cons)==1){
    starting_cons <- rep(starting_cons,length(compounds))
  }else{
    
  }
  
  dflist <- NULL
  for(i in 1:length(compounds)){
    int <- twoFold(starting_cons[i],compound=compounds[i],series=series)
    dflist[[i]] <- int
  }
  cc <- do.call(rbind,dflist) %>%
    dplyr::rename(compound_name=compound,
                  conc_val=conc)
  return(cc)
}

readin_meta_csv <- function(directory,na.value=0,recursive=F){
  print(paste("looking in",directory))
  setwd(as.character(directory))
  files <- dir(pattern="csv$",recursive = recursive)
  if(length(files)==0){
    print("no .csv files found")
  }else{
    cat(paste("found csvs:\n",paste(dir(pattern="csv$",recursive = T),sep="\n")))
    
    dflist <- NULL
    for(i in 1:length(files)){
      print(i)
      #i=1
      int <- read.csv(file=files[i],
                      stringsAsFactors = F) %>%
        select(-Sample,-X) %>%
        filter(!X.1=="Name") %>%
        dplyr::rename(sampleid=X.1,
                      Data.File=X.2,
                      Type=X.3,
                      Level=X.4,
                      Acq.Date.Time=X.5) %>%
        select(-Type,-Level) %>%
        mutate(filename=files[i]) %>%
        reshape2::melt(id.vars=c("sampleid","Data.File",
                                 "Acq.Date.Time","filename")) %>%
        mutate(value=suppressWarnings(as.numeric(value))) %>%
        spread(variable,value,fill=na.value) %>%
        reshape2::melt(id.vars=c("sampleid","Data.File",
                                 "Acq.Date.Time","filename")) %>%
        dplyr::rename(compound=variable,peakarea=value) %>%
        mutate(conc=ifelse(grepl("\\_[Cc]on",sampleid),"concentrated",
                           ifelse(grepl("\\_dil",sampleid),"diluted",
                                  "other")),
               compound=as.character(compound),
               itsd=str_extract(compound,pattern="ITSD"),
               status=str_extract(compound,pattern="\\_.+"),
               status=gsub("\\.Results","",status),
               status=gsub("\\_ITSD","",status),
               status=gsub("^\\_","",status),
               compound_name=gsub("\\_ITSD.+|\\.[Rr]esults","",compound)) %>%
        replace_na(list(status="")) %>%
        select(sampleid,Data.File,Acq.Date.Time,compound,
               compound_name,conc,itsd,peakarea,status,filename) %>%
        mutate(peakarea=as.numeric(peakarea),
               peakarea=ifelse(peakarea < na.value,na.value,peakarea))
      
      dflist[[i]] <- int
    }
    
    meta <- do.call(rbind,dflist)
    rm(dflist)
    return(meta)
  }
}

readin_meta_csv_single_file <- function(filename,na.value=0,recursive=F){
  
  int <- read.csv(file=filename,
                  stringsAsFactors = F) %>%
    select(-Sample,-X) %>%
    filter(!X.1=="Name") %>%
    dplyr::rename(sampleid=X.1,
                  Data.File=X.2,
                  Type=X.3,
                  Level=X.4,
                  Acq.Date.Time=X.5) %>%
    select(-Type,-Level) %>%
    mutate(filename=filename) %>%
    reshape2::melt(id.vars=c("sampleid","Data.File",
                             "Acq.Date.Time","filename")) %>%
    mutate(value=suppressWarnings(as.numeric(value))) %>%
    spread(variable,value,fill=na.value) %>%
    reshape2::melt(id.vars=c("sampleid","Data.File",
                             "Acq.Date.Time","filename")) %>%
    dplyr::rename(compound=variable,peakarea=value) %>%
    mutate(conc=ifelse(grepl("\\_[Cc]on",sampleid),"concentrated",
                       ifelse(grepl("\\_dil",sampleid),"diluted",
                              "other")),
           compound=as.character(compound),
           itsd=str_extract(compound,pattern="ITSD"),
           status=str_extract(compound,pattern="\\_.+"),
           status=gsub("\\.Results","",status),
           status=gsub("\\_ITSD","",status),
           status=gsub("^\\_","",status),
           compound_name=gsub("\\_ITSD.+|\\.[Rr]esults","",compound)) %>%
    replace_na(list(status="")) %>%
    select(sampleid,Data.File,Acq.Date.Time,compound,
           compound_name,conc,itsd,peakarea,status,filename) %>%
    mutate(peakarea=as.numeric(peakarea),
           peakarea=ifelse(peakarea < na.value,na.value,peakarea))
  
  meta <- int
  return(meta)
}

make_norm_conc_tbl <- function(df,
                               dil_compounds="Valerate,Phenol,Valine_D8,Proline_D7, Succinate",
                               conc_compounds="Phenol"){
  
  #make df with average for conc and diluted
  #dil_compounds=c("Valine_D8","Valerate")
  #conc_compounds=c("Succinate","Proline","Phenol")
  #dil_compounds <- "Valine_D8,Valerate"
  
  dil_compounds <- trimws(unlist(strsplit(dil_compounds,split=",")))
  conc_compounds <- trimws(unlist(strsplit(conc_compounds,split=",")))
  
  int_conc <- tibble(compound_name = c(dil_compounds, conc_compounds),
                     conc = c(rep("diluted",length(dil_compounds)),
                              rep("concentrated",length(conc_compounds)))) %>% 
    inner_join(df) %>%
    filter(itsd=="ITSD") %>%
    # dplyr::count(compound_name, conc) %>%
    # dplyr::count(sampleid) %>%
    group_by(sampleid) %>%
    summarize(avg=mean(peakarea))
  
  return(int_conc)
  
}

get_indole_conc <- function(conc,compounds,series=11){
  #compounds <- inputcompounds2
  compounds <- "niacin,tyrosine,phenylalanine,kynurenine,Serotonin,anthranilicacid,tryptophan,5HIAA,tryptamine,kynurenicacid,melatonin"
  series <- 11
  conc <- "909,454.5,227.25,113.625,56.8125,5.68125,0.568125,0.056813,0.014203,0.003551,0.000888"
  
  compounds = unlist(strsplit(compounds, split=","))
  conc = unlist(strsplit(conc,split = ","))
  # # 
  # indole_conc <- read.csv("/Volumes/chaubard-lab/shiny_workspace/other_data/CalibratorConcentrations_ForShinyEric_INDOLE.csv",
  #                       stringsAsFactors = F) %>%
  # select(curveLab=X,
  #        conc_val=Conc.uM) # This produces a df with
  # 
  dflist <- NULL
  for(i in 1:series){
    #print(i)
    curve <- paste0("cc",i)
    dflist[[i]] <- curve
  }
  #class(dflist)
  
  curves <- do.call(rbind,as.list(dflist)) %>%
    as.data.frame() %>%
    mutate(rownum=row_number()) %>%
    arrange(-rownum) %>%
    select(-rownum)
  
  curves <- cbind(curves,"conc_val"=conc)
  colnames(curves)[1] <- "curveLab"  
  
  
  dflist <- NULL
  for(i in 1:length(compounds)){
    #i=1
    int <- cbind("compound_name"=rep(compounds[i],series),curves)
    dflist[[i]] <- int
  }
  indole_conc <- do.call(rbind,as.list(dflist)) %>%
    mutate(conc_val=as.numeric(as.character(conc_val)))
  return(indole_conc)
}



# set up directory and files ----------------------------------------------

#system("open .")
wddir <- "/Volumes/chaubard-lab/shiny_workspace/csvs/"

# ui ----------------------------------------------------------------------

ui <- fluidPage(
  # shinythemes::themeSelector(),
  titlePanel("DFI Metabolomics QC (v1.5)"),
  br(),
  
  # CSV file selector -------------------------------------------------------
  
  fluidRow(column(width = 4,
                  wellPanel(
                    actionButton("refresh_csv", "Refresh CSV files"),
                    selectInput("filename", "Select a CSV file from: ", list.files(wddir, pattern="csv$"))
                  )
  )
  ),
  
  # tabs --------------------------------------------------------------------
  
  tabsetPanel(type="tabs",
              # PFFBr quant QC UI -------------------------------------------------------------
              tabPanel("PFFBr Quant QC", theme = shinytheme("flatly"),
                       sidebarLayout(
                         sidebarPanel(width = 3,
                                      textInput("compounds","Enter ITSD compounds (comma separated):",
                                                value="Acetate,Propionate,Butyrate,Succinate"),
                                      br(),
                                      h4("ITSD dilution calculation"),
                                      numericInput("xfactor","Mult factor:",value = 11),
                                      #numericInput("start","Enter concentration(s):",value = 100),
                                      textInput("start","Enter concentration(s):","100"),
                                      numericInput("series","dilution #",8),
                                      br(),
                                      h4("Filters:"),
                                      textInput("maxcc","Max conc(s) filter:","100,100,100,100"),
                                      textInput("mincc","Min conc(s) filter:","0,0,0,0"),
                         ),
                         mainPanel(
                           plotOutput("quant"),
                           h4("Fitted linear model stats:"),
                           dataTableOutput("model"),
                           h4("Standardized table:"),
                           dataTableOutput("quant_tbl"),
                           downloadButton("downloadData", "Download OFFBr Quant Table")
                         )
                       )
              ),
              
              # PFFBr Normalization UI --------------------------------------------------------
              
              tabPanel("PFFBr Normalization",fluidPage(theme = shinytheme("flatly")),
                       sidebarLayout(
                         sidebarPanel(width = 3,
                                      br(),
                                      h4("Type in ITSD"),
                                      textInput("dil_compounds","diluted standards:", value="Valine_D8,Valerate"),
                                      textInput("conc_compounds","concentrated standards:",value="Proline_D7,Phenol"),
                                      br(),
                                      numericInput("zero_val","minimum value:",value=5000)
                         ),
                         mainPanel(
                           splitLayout(cellWidths = c("25%","75%"),
                                       uiOutput("compound_list"),
                                       plotOutput("raw_boxplots",height="1400px")),
                           h4("Intermediate table:"),
                           dataTableOutput("conc_filter"),
                           h4("Normalized table:"),
                           dataTableOutput("normwide"),
                           checkboxInput("qcfil","RemoveQCs"),
                           downloadButton("downloadData2", "Download PFFBr Normalized Table")
                         )
                       )
              ),

            # indole quant ------------------------------------------------------------
              tabPanel("Indole Quant QC", theme = shinytheme("flatly"),
                       sidebarLayout(
                         sidebarPanel(width = 3,
                                      textInput("compounds2","Enter ITSD compounds (comma separated):",
                                                value="niacin,tyrosine,phenylalanine,kynurenine,Serotonin,anthranilicacid,tryptophan,5HIAA,tryptamine,kynurenicacid,melatonin"),
                                      br(),
                                      h4("ITSD dilution calculation"),
                                      #numericInput("xfactor2","Mult factor:",value = 11),
                                      #numericInput("start","Enter concentration(s):",value = 100),
                                      textInput("start2","Enter concentration(s):","909,454.5,227.25,113.625,56.8125,5.68125,0.568125,0.056813,0.014203,0.003551,0.000888"),
                                      numericInput("series2","dilution #",11),
                                      br(),
                                      h4("Filters:"),
                                      textInput("maxcc2","Max conc(s) filter:","909"),
                                      textInput("mincc2","Min conc(s) filter:","0"),
                         ),
                         mainPanel(
                           plotOutput("quant2"),
                           h4("Fitted linear model stats:"),
                           dataTableOutput("model2"),
                           h4("Standardized table:"),
                           dataTableOutput("quant_tbl2"),
                           downloadButton("downloadData3", "Download Indole Quant Table")
                         )
                       )
              ),
              
              # Indone normalization UI --------------------------------------------------------
              
              tabPanel("Indole Normalization",fluidPage(theme = shinytheme("flatly")),
                       sidebarLayout(
                         sidebarPanel(width = 3,
                                      br(),
                                      h4("Type in ITSD"),
                                      textInput("itsd_compounds","standards:", value="Valine_D8,Valerate"),
                                      #textInput("conc_compounds","concentrated standards:",value="Proline_D7,Phenol"),
                                      br(),
                                      numericInput("zero_val2","minimum value:",value=5000)
                         ),
                         mainPanel(
                           splitLayout(cellWidths = c("25%","75%"),
                                       uiOutput("compound_list2"),
                                       plotOutput("raw_boxplots2",height="1400px")),
                           h4("Intermediate table:"),
                           dataTableOutput("conc_filter2"),
                           h4("Normalized table:"),
                           dataTableOutput("normwide2"),
                           checkboxInput("qcfil2","RemoveQCs"),
                           downloadButton("downloadData4", "Download Normalized Table")
                         )
                       )
              ),

              
              # more --------------------------------------------------------------------
              
              tabPanel("More",fluidPage(theme = shinytheme("flatly")),
                       h3("In construction for some awesome stuff!")
              )
              
  )
  
)

# server ------------------------------------------------------------------

server <- function(input, output, session) {
  
  
  # refresh CSV list when hit button
  observeEvent(input$refresh_csv,ignoreInit = T,ignoreNULL = T, {
    print(list.files(wddir, pattern = "csv$"))
    updateSelectInput(session, 'filename', choices = list.files(wddir, pattern = "csv$"))
  })
  
  
  # PFFBr quant qc tab ------------------------------------------------------------
  
  #make concentration table
  conc_tbl <- reactive({
    get_all_conc(input$start,
                 compounds=input$compounds,series=input$series)
  })
  
  output$conc <- renderDataTable(
    conc_tbl(),
    options = list(pageLength=5)
  )
  
  #make cutoff reactive
  cutoff_df <- reactive({
    compounds = unlist(strsplit(input$compounds, split=","))
    max_cutoffs = unlist(strsplit(input$maxcc, split = ","))
    min_cutoffs = unlist(strsplit(input$mincc, split = ","))
    
    cutoff_df <- cbind(compounds,max_cutoffs,min_cutoffs)
    colnames(cutoff_df) <- c("compound_name","maxcc","mincc")
    
    cutoff_df %>%
      as.data.frame() %>%
      mutate(maxcc=as.character(maxcc),
             mincc=as.character(mincc),
             maxcc=as.numeric(maxcc),
             mincc=as.numeric(mincc))
  })
  
  # read in table as reactive 
  meta <- reactive({ readin_meta_csv_single_file(file.path(wddir,input$filename)) })
  
  #show models and make plots
  modelstart <- reactive({
    
    compounds = unlist(strsplit(input$compounds, split=","))
    compounds = factor(compounds,level=unique(compounds))
    
    model <- meta() %>%
      filter(compound_name %in% compounds,
             conc=="diluted") %>%
      mutate(compound_name=factor(compound_name,levels=compounds)) %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(sampleid+compound_name+conc ~ itsd,value.var="peakarea") %>%
      mutate(#peak = ifelse(peak <= 10000,0,peak),
        norm_peak=peak / ITSD) %>%
      mutate(curveLab=str_extract(sampleid,"\\_CC[1-9][0-9]+\\_|\\_CC[1-9]+\\_"),
             curveLab=gsub("\\_","",curveLab)) %>%
      left_join(conc_tbl()) %>%
      left_join(cutoff_df()) %>%
      filter(conc_val <= maxcc,
             conc_val >= mincc) %>%
      group_by(compound_name) %>%
      summarize(r = cor(norm_peak,conc_val),
                model_list <- broom::tidy(lm(norm_peak ~ conc_val))) %>%
      reshape2::dcast(compound_name+r ~ term,value.var="estimate") %>%
      dplyr::rename(slope_value=conc_val) 
  })
  
  output$model <- renderDataTable(
    modelstart() %>%
      datatable() %>%
      formatRound(c(2:4), 3) %>% 
      formatStyle(columns = c(2:4), 'text-align' = 'center')
  )
  
  
  
  #make quant graph
  quant_plot <- reactive({
    
    compounds = unlist(strsplit(input$compounds, split=","))
    compounds = factor(compounds,level=unique(compounds))
    
    meta() %>%
      filter(compound_name %in% compounds,
             conc=="diluted") %>%
      mutate(compound_name=factor(compound_name,levels=unique(compounds))) %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(sampleid+compound_name+conc ~ itsd,value.var="peakarea") %>%
      mutate(#peak = ifelse(peak <= 10000,0,peak),
        norm_peak=peak / ITSD) %>%
      filter(grepl("\\_CC[0-9]",sampleid)) %>%
      mutate(curveLab=str_extract(sampleid,"\\_CC[1-9][0-9]+\\_|\\_CC[1-9]+\\_"),
             curveLab=gsub("\\_","",curveLab)) %>%
      left_join(conc_tbl()) %>%
      left_join(cutoff_df()) %>%
      filter(conc_val <= maxcc,
             conc_val >= mincc) %>%
      ggplot(aes(x=conc_val,y=norm_peak)) +
      geom_point(size=3) +
      geom_smooth(method="lm") +
      theme(strip.text=element_text(size=15),
            axis.text=element_text(size=13)) +
      facet_wrap("compound_name",scales="free")
  })
  
  output$quant <- renderPlot(
    quant_plot()
  )
  
  #make quant table
  quant_table <- reactive({
    
    compounds = unlist(strsplit(input$compounds, split=","))
    
    
    meta() %>%
      filter(compound_name %in% compounds,
             conc=="diluted") %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(sampleid+compound_name+conc ~ itsd,value.var="peakarea") %>%
      mutate(#peak = ifelse(peak <= 10000,0,peak),
        norm_peak=peak / ITSD) %>%
      filter(!grepl("\\_CC[0-9]",sampleid)) %>%
      left_join(modelstart()) %>%
      mutate(quant_val =  (norm_peak - (`(Intercept)`))/slope_value*as.numeric(input$xfactor)) %>%
      arrange(compound_name) %>%
      mutate(quant_val = ifelse(quant_val < 0,0,quant_val),
             quant_val = round(quant_val,2))
  })
  
  output$quant_tbl <- renderDataTable(
    quant_table() %>%
      datatable() %>%
      formatRound(c(4:ncol(quant_table())), 3) %>% 
      formatStyle(columns = c(4:ncol(quant_table())), 'text-align' = 'center')
  )
  
  quant_table_dl <- reactive({
    
    compounds = unlist(strsplit(input$compounds, split=","))
    #compounds = factor(compounds,level=unique(compounds))
    
    meta() %>%
      filter(compound_name %in% compounds,
             conc=="diluted") %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(sampleid+compound_name+conc ~ itsd,value.var="peakarea") %>%
      mutate(#peak = ifelse(peak <= 10000,0,peak),
        norm_peak=peak / ITSD) %>%
      filter(!grepl("\\_CC[0-9]",sampleid)) %>%
      left_join(modelstart()) %>%
      mutate(quant_val =  (norm_peak - (`(Intercept)`))/slope_value*as.numeric(input$xfactor)) %>%
      arrange(compound_name) %>%
      mutate(quant_val = ifelse(quant_val < 0,0,quant_val),
             quant_val = round(quant_val,2)) %>%
      reshape2::dcast(sampleid ~ compound_name, value.var="quant_val",fill=0) %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="\\_\\_") %>%
      arrange(num) %>%
      select(-num,-date,-batch,-conc)
      
  })
  
  #download handler
  output$downloadData <- downloadHandler(
    
    filename = function(){
      paste0("quant_results_",
             gsub("\\.csv","",input$filename),"_",
             gsub("\\-","",Sys.Date()),".csv")
    },
    
    content = function(file) {
      write.csv(quant_table_dl(),file,
                row.names=F,quote=F)
    }
  )
  
  
  # PFFBr normalization tab -------------------------------------------------------
  
  rawdf <- reactive({ 
    readin_meta_csv_single_file(file.path(wddir,input$filename),na.value = input$zero_val) 
  })
  
  csv <- reactive({
    
    vars <- rawdf() %>%
      mutate(compound_name=as.character(compound_name)) %>%
      arrange(desc(compound_name)) %>%
      filter(is.na(itsd)) %>%
      distinct(compound_name) %>%
      `$`(compound_name)
    
    return(vars)
  })
  
  output$compound_list <- renderUI({
    checkboxGroupInput("check_compounds","Check = diluted",
                       choices=csv(),
                       inline=F)
  })
  
  conc_int <- reactive({
    
    make_norm_conc_tbl(rawdf(),
                       dil_compounds = as.character(input$dil_compounds),
                       conc_compounds = as.character(input$conc_compounds))
  })
  
  # output$conc_int <- renderDataTable(
  #   conc_int()
  # )
  
  #make boxplots:
  raw_boxplots <- reactive({
    
    rawdf() %>%
      filter(is.na(itsd)) %>%
      replace_na(list(itsd="not itsd")) %>%
      ggplot(aes(y=compound_name,x=peakarea)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width=0.1,height=0.1,alpha=0.3) +
      theme_bw() +
      theme(strip.text.x=element_text(size=15),
            axis.text.y=element_text(size=13)) +
      ylab("") +
      xlab("peak area") +
      facet_grid(itsd ~ conc,scales="free_y",space="free") +
      scale_x_continuous(trans=log_epsilon_trans(epsilon=1000))
    
  })
  
  output$raw_boxplots <- renderPlot(
    raw_boxplots()
  )
  
  #subset based on checkboxes.. make a dataframe of compounds in same order as checkboxes with input as a column
  conc_filter <- reactive({
    # print(input$check_compounds)
    # print(class(input$check_compounds))
    # cat(input$check_compounds)
    # print(length(input$check_compounds))
    # print(length(as.list(input$check_compounds)))
    
    if(length(as.list(input$check_compounds)) > 0){
      # print("if")
      
      leftovervars <- csv()[ !(csv() %in% input$check_compounds) ]
      
      tibble(checked = rep("diluted", length(input$check_compounds)),
             compound_name = input$check_compounds) %>%
        bind_rows(tibble(checked = rep("concentrated", length(leftovervars)),
                         compound_name = leftovervars)
        )
    }else{
      
      # print("else")
      tibble(checked = "concentrated",
             compound_name = csv()) 
    }
  })
  
  
  
  subset_conc <- reactive({
    rawdf()
    
    #is 33 rows when should be 66.. the conc averages aren't being included !
    #print(conc_int())
    
    rawdf() %>%
      left_join(conc_filter()) %>%
      # replace_na(list(checked="concentrated")) %>%
      filter(conc==checked) %>%
      select(-checked) %>%
      left_join(conc_int()) %>%
      mutate(norm_peak=peakarea / avg)
  })
  
  output$conc_filter <- renderDataTable(
    subset_conc() %>%
      datatable() %>%
      formatRound(c("norm_peak"), 3) %>% 
      formatStyle(columns = c("norm_peak"), 'text-align' = 'center')
  )
  
  #download table
  #make wide table 
  norm_wide_tbl <- reactive({
    
    if(input$qcfil==F){
    rawdf() %>%
      left_join(conc_filter()) %>%
      # replace_na(list(checked="concentrated")) %>%
      filter(conc==checked, is.na(itsd),
             !grepl("CC[0-9]+", sampleid)) %>%
      #select(-checked) %>%
      left_join(conc_int()) %>%
      mutate(norm_peak=peakarea / avg) %>%
      dplyr::select(sampleid, compound_name, norm_peak) %>%
      spread(compound_name, norm_peak, fill = NA) %>%
      reshape2::melt(id.vars=c("sampleid")) %>%
      dplyr::rename(compound_name=variable,norm_peak=value) %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="\\_\\_") %>%
      filter(!is.na(norm_peak)) %>%
      # filter(sampleid %!in% c("MB_1","MB_2",
      #                        "PooledQC",
      #                        "BHIQC_1",
      #                        "BHIQC_2")) %>%
      dplyr::select(num,sampleid, compound_name, norm_peak) %>%
      mutate(norm_peak = round(norm_peak,5)) %>%
      reshape2::dcast(num+sampleid ~ compound_name,value.var="norm_peak",fun.aggregate=mean) %>%
      arrange(num) %>%
      select(-num)
    }else{
      rawdf() %>%
        left_join(conc_filter()) %>%
        # replace_na(list(checked="concentrated")) %>%
        filter(conc==checked, is.na(itsd),
               !grepl("CC[0-9]+", sampleid)) %>%
        #select(-checked) %>%
        left_join(conc_int()) %>%
        mutate(norm_peak=peakarea / avg) %>%
        dplyr::select(sampleid, compound_name, norm_peak) %>%
        spread(compound_name, norm_peak, fill = NA) %>%
        reshape2::melt(id.vars=c("sampleid")) %>%
        dplyr::rename(compound_name=variable,norm_peak=value) %>%
        separate(sampleid,into=c("num","date","batch","sampleid","conc"),
                 sep="\\_\\_") %>%
        filter(!is.na(norm_peak)) %>%
        filter(sampleid %!in% c("MB_1","MB_2",
                               "PooledQC",
                               "BHIQC_1",
                               "BHIQC_2")) %>%
        dplyr::select(num,sampleid, compound_name, norm_peak) %>%
        mutate(norm_peak = round(norm_peak,5)) %>%
        reshape2::dcast(num+sampleid ~ compound_name,value.var="norm_peak",fun.aggregate=mean) %>%
        arrange(num) %>%
        select(-num)
    }
  })
  
  output$normwide <- renderDataTable(
    norm_wide_tbl() %>%
      datatable() %>%
      formatRound(c(2:ncol(norm_wide_tbl())), 3) %>% 
      formatStyle(columns = c(2:ncol(norm_wide_tbl())), 'text-align' = 'center')
  )
  
  #download handler
  output$downloadData2 <- downloadHandler(
    
    filename = function(){
      paste0("normalized_results_",input$filename,"_",Sys.Date(),".csv")
    },
    
    content = function(file) {
      write.csv(norm_wide_tbl(),file,
                row.names=F,quote=F)
    }
  )
  

# indole quant ------------------------------------------------------------

  
  #make concentration table
  conc_tbl2 <- reactive({
    get_indole_conc(conc=input$start2,
                    compounds=input$compounds2,
                    series=input$series2)
  })
  
  output$conc2 <- renderDataTable(
    conc_tbl2(),
    options = list(pageLength=5)
  )
  
  #make cutoff reactive
  cutoff_df2 <- reactive({
    compounds2 = unlist(strsplit(input$compounds2, split=","))
    max_cutoffs2 = unlist(strsplit(input$maxcc2, split = ","))
    min_cutoffs2 = unlist(strsplit(input$mincc2, split = ","))
    
    cutoff_df2 <- cbind(compounds2,max_cutoffs2,min_cutoffs2)
    colnames(cutoff_df2) <- c("compound_name","maxcc","mincc")
    
    cutoff_df2 %>%
      as.data.frame() %>%
      mutate(maxcc=as.character(maxcc),
             mincc=as.character(mincc),
             maxcc=as.numeric(maxcc),
             mincc=as.numeric(mincc))
  })
  
  # read in table as reactive 
  meta <- reactive({ readin_meta_csv_single_file(file.path(wddir,input$filename)) })
  
  #show models and make plots
  modelstart2 <- reactive({
    
    compounds2 = unlist(strsplit(input$compounds2, split=","))
    compounds2 = factor(compounds2,level=unique(compounds2))
    
      
      meta() %>%
        mutate(compound_name=gsub("\\_[0-9]+$","",compound_name)) %>%
        filter(compound_name %in% compounds2) %>%
        mutate(compound_name=factor(compound_name,levels=unique(compounds2))) %>%
        replace_na(list(itsd="peak")) %>%
        reshape2::dcast(Data.File+compound_name+conc ~ itsd,value.var="peakarea") %>%
        mutate(#peak = ifelse(peak <= 10000,0,peak),
          norm_peak=peak / ITSD) %>%
        filter(grepl("\\_[Cc][Cc][1-9][0-9]+\\_|\\_[Cc][Cc][1-9]+\\_", Data.File)) %>%
        mutate(curveLab=str_extract(Data.File,"\\_[Cc][Cc][1-9][0-9]+\\_|\\_[Cc][Cc][1-9]+\\_"),
               curveLab=gsub("\\_","",curveLab)) %>%
        left_join(conc_tbl2()) %>%
        left_join(cutoff_df2()) %>%
        filter(conc_val <= maxcc,
               conc_val >= mincc) %>%
        group_by(compound_name) %>%
        summarize(r = cor(norm_peak,conc_val),
                  model_list <- broom::tidy(lm(norm_peak ~ conc_val))) %>%
        reshape2::dcast(compound_name+r ~ term,value.var="estimate") %>%
        dplyr::rename(slope_value=conc_val) 
    
  })
  
  output$model2 <- renderDataTable(
    modelstart2() %>%
      datatable() %>%
      formatRound(c(2:4), 3) %>% 
      formatStyle(columns = c(2:4), 'text-align' = 'center')
  )
  
  
  
  #make quant graph
  quant_plot2 <- reactive({
    
    compounds2 = unlist(strsplit(input$compounds2, split=","))
    compounds2 = factor(compounds2,level=unique(compounds2))
    
    meta() %>%
      mutate(compound_name=gsub("\\_[0-9]+$","",compound_name)) %>%
      filter(compound_name %in% compounds2) %>%
      mutate(compound_name=factor(compound_name,levels=unique(compounds2))) %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(Data.File+compound_name+conc ~ itsd,value.var="peakarea") %>%
      mutate(#peak = ifelse(peak <= 10000,0,peak),
        norm_peak=peak / ITSD) %>%
      filter(grepl("\\_[Cc][Cc][1-9][0-9]+\\_|\\_[Cc][Cc][1-9]+\\_", Data.File)) %>%
      mutate(curveLab=str_extract(Data.File,"\\_[Cc][Cc][1-9][0-9]+\\_|\\_[Cc][Cc][1-9]+\\_"),
             curveLab=gsub("\\_","",curveLab)) %>%
      left_join(conc_tbl2()) %>%
      left_join(cutoff_df2()) %>%
      filter(conc_val <= maxcc,
             conc_val >= mincc) %>%
      ggplot(aes(x=conc_val,y=norm_peak)) +
      geom_point(size=3) +
      geom_smooth(method="lm") +
      theme(strip.text=element_text(size=15),
            axis.text=element_text(size=13)) +
      facet_wrap("compound_name",scales="free")
  })
  
  output$quant2 <- renderPlot(
    quant_plot2()
  )
  
  #make quant table
  quant_table2 <- reactive({
    
    compounds2 = unlist(strsplit(input$compounds2, split=","))
    
    meta() %>% 
      mutate(compound_name=gsub("\\_[0-9]+$","",compound_name)) %>%
      filter(compound_name %in% compounds2) %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(Data.File+compound_name ~ itsd,value.var="peakarea") %>%
      mutate(#peak = ifelse(peak <= 10000,0,peak),
        norm_peak=peak / ITSD) %>%
      filter(!grepl("\\_[Cc][Cc][0-9]",Data.File)) %>%
      left_join(modelstart2()) %>%
      mutate(quant_val =  (norm_peak - (`(Intercept)`))/slope_value) %>%
      arrange(compound_name) %>%
      mutate(quant_val = ifelse(quant_val < 0,0,quant_val),
             quant_val = round(quant_val,2))
  

  })
  
  output$quant_tbl2 <- renderDataTable(
    quant_table2() %>%
      datatable() %>%
      formatRound(c(4:ncol(quant_table2())), 3) %>% 
      formatStyle(columns = c(4:ncol(quant_table2())), 'text-align' = 'center')
  )
  
  quant_table_dl2 <- reactive({
    
    compounds2 = unlist(strsplit(input$compounds2, split=","))
    compounds2 = factor(compounds2,level=unique(compounds2))
    
    meta() %>% 
      mutate(compound_name=gsub("\\_[0-9]+$","",compound_name)) %>%
      filter(compound_name %in% compounds2) %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(Data.File+compound_name ~ itsd,value.var="peakarea") %>%
      mutate(#peak = ifelse(peak <= 10000,0,peak),
        norm_peak=peak / ITSD) %>%
      filter(!grepl("\\_[Cc][Cc][0-9]",Data.File)) %>%
      left_join(modelstart2()) %>%
      mutate(quant_val =  (norm_peak - (`(Intercept)`))/slope_value) %>%
      arrange(compound_name) %>%
      mutate(quant_val = ifelse(quant_val < 0,0,quant_val),
             quant_val = round(quant_val,2)) %>%
      reshape2::dcast(Data.File ~ compound_name, value.var="quant_val",fill=0) %>%
      separate(Data.File,into=c("num","date","batch","sampleid","conc"),
               sep="\\_\\_") %>%
      arrange(num) %>%
      select(-num,-date,-batch,-conc)
    
    
  })
  
  #download handler
  output$downloadData3 <- downloadHandler(
    
    filename = function(){
      paste0("indole_quant_results_",
             gsub("\\.csv","",input$filename),"_",
             gsub("\\-","",Sys.Date()),".csv")
    },
    
    content = function(file) {
      write.csv(quant_table_dl2(),file,
                row.names=F,quote=F)
    }
  )
  
  
  # normalization tab -------------------------------------------------------
  
  rawdf <- reactive({ 
    readin_meta_csv_single_file(file.path(wddir,input$filename),na.value = input$zero_val) 
  })
  
  csv <- reactive({
    
    vars <- rawdf() %>%
      mutate(compound_name=as.character(compound_name)) %>%
      arrange(desc(compound_name)) %>%
      filter(is.na(itsd)) %>%
      distinct(compound_name) %>%
      `$`(compound_name)
    
    return(vars)
  })
  
  output$compound_list2 <- renderUI({
    checkboxGroupInput("check_compounds","Check = diluted",
                       choices=csv(),
                       inline=F)
  })
  
  conc_int2 <- reactive({
    
    make_norm_conc_tbl(rawdf(),
                       dil_compounds = as.character(input$itsd_compounds))
                       #conc_compounds = as.character(input$conc_compounds))
  })
  
  # output$conc_int <- renderDataTable(
  #   conc_int()
  # )
  
  #make boxplots:
  raw_boxplots2 <- reactive({
    
    rawdf() %>%
      filter(is.na(itsd)) %>%
      replace_na(list(itsd="not itsd")) %>%
      ggplot(aes(y=compound_name,x=peakarea)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width=0.1,height=0.1,alpha=0.3) +
      theme_bw() +
      theme(strip.text.x=element_text(size=15),
            axis.text.y=element_text(size=13)) +
      ylab("") +
      xlab("peak area") +
      facet_grid(itsd ~ conc,scales="free_y",space="free") +
      scale_x_continuous(trans=log_epsilon_trans(epsilon=1000))
    
  })
  
  output$raw_boxplots2 <- renderPlot(
    raw_boxplots2()
  )
  
  #subset based on checkboxes.. make a dataframe of compounds in same order as checkboxes with input as a column
  conc_filter2 <- reactive({
    # print(input$check_compounds)
    # print(class(input$check_compounds))
    # cat(input$check_compounds)
    # print(length(input$check_compounds))
    # print(length(as.list(input$check_compounds)))
    
    if(length(as.list(input$check_compounds2)) > 0){
      # print("if")
      
      leftovervars <- csv()[ !(csv() %in% input$check_compounds2) ]
      
      tibble(checked = rep("diluted", length(input$check_compounds2)),
             compound_name = input$check_compounds2) %>%
        bind_rows(tibble(checked = rep("concentrated", length(leftovervars)),
                         compound_name = leftovervars)
        )
    }else{
      
      # print("else")
      tibble(checked = "concentrated",
             compound_name = csv()) 
    }
  })
  
  
  
  subset_conc2 <- reactive({
    rawdf()
    
    #is 33 rows when should be 66.. the conc averages aren't being included !
    #print(conc_int())
    
    rawdf() %>%
      left_join(conc_filter2()) %>%
      # replace_na(list(checked="concentrated")) %>%
      filter(conc==checked) %>%
      select(-checked) %>%
      left_join(conc_int2()) %>%
      mutate(norm_peak=peakarea / avg)
  })
  
  output$conc_filter2 <- renderDataTable(
    subset_conc2() %>%
      datatable() %>%
      formatRound(c("norm_peak"), 3) %>% 
      formatStyle(columns = c("norm_peak"), 'text-align' = 'center')
  )
  
  #download table
  #make wide table 
  norm_wide_tbl2 <- reactive({
    
    if(input$qcfil==F){
      rawdf() %>%
        left_join(conc_filter2()) %>%
        # replace_na(list(checked="concentrated")) %>%
        filter(conc==checked, is.na(itsd),
               !grepl("CC[0-9]+", sampleid)) %>%
        #select(-checked) %>%
        left_join(conc_int2()) %>%
        mutate(norm_peak=peakarea / avg) %>%
        dplyr::select(sampleid, compound_name, norm_peak) %>%
        spread(compound_name, norm_peak, fill = NA) %>%
        reshape2::melt(id.vars=c("sampleid")) %>%
        dplyr::rename(compound_name=variable,norm_peak=value) %>%
        separate(sampleid,into=c("num","date","batch","sampleid","conc"),
                 sep="\\_\\_") %>%
        filter(!is.na(norm_peak)) %>%
        # filter(sampleid %!in% c("MB_1","MB_2",
        #                        "PooledQC",
        #                        "BHIQC_1",
        #                        "BHIQC_2")) %>%
        dplyr::select(num,sampleid, compound_name, norm_peak) %>%
        mutate(norm_peak = round(norm_peak,5)) %>%
        reshape2::dcast(num+sampleid ~ compound_name,value.var="norm_peak",fun.aggregate=mean) %>%
        arrange(num) %>%
        select(-num)
    }else{
      rawdf() %>%
        left_join(conc_filter2()) %>%
        # replace_na(list(checked="concentrated")) %>%
        filter(conc==checked, is.na(itsd),
               !grepl("CC[0-9]+", sampleid)) %>%
        #select(-checked) %>%
        left_join(conc_int2()) %>%
        mutate(norm_peak=peakarea / avg) %>%
        dplyr::select(sampleid, compound_name, norm_peak) %>%
        spread(compound_name, norm_peak, fill = NA) %>%
        reshape2::melt(id.vars=c("sampleid")) %>%
        dplyr::rename(compound_name=variable,norm_peak=value) %>%
        separate(sampleid,into=c("num","date","batch","sampleid","conc"),
                 sep="\\_\\_") %>%
        filter(!is.na(norm_peak)) %>%
        filter(sampleid %!in% c("MB_1","MB_2",
                                "PooledQC",
                                "BHIQC_1",
                                "BHIQC_2")) %>%
        dplyr::select(num,sampleid, compound_name, norm_peak) %>%
        mutate(norm_peak = round(norm_peak,5)) %>%
        reshape2::dcast(num+sampleid ~ compound_name,value.var="norm_peak",fun.aggregate=mean) %>%
        arrange(num) %>%
        select(-num)
    }
  })
  
  output$normwide2 <- renderDataTable(
    norm_wide_tbl2() %>%
      datatable() %>%
      formatRound(c(2:ncol(norm_wide_tbl2())), 3) %>% 
      formatStyle(columns = c(2:ncol(norm_wide_tbl())), 'text-align' = 'center')
  )
  
  #download handler
  output$downloadData4 <- downloadHandler(
    
    filename = function(){
      paste0("normalized_results_",input$filename,"_",Sys.Date(),".csv")
    },
    
    content = function(file) {
      write.csv(norm_wide_tbl2(),file,
                row.names=F,quote=F)
    }
  )
}



# run app -----------------------------------------------------------------

#shinyApp(ui, server)
runApp(list(ui=ui,server=server),host="0.0.0.0",port=5000)
