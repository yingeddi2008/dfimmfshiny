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
library(pheatmap)
library(ggsci)
library(zoo)
library(ggrepel) # New for Version 1.8
library(ggsci) # New for Version 1.8
library(gridExtra) # New for Version 1.8
library(ggpmisc) # New for Version 1.8

setwd("/Volumes/chaubard-lab/shiny_workspace/csvs/")

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
anyFold <- function(startCon, 
                    compound = c("Butyrate", "Propionate","Acetate"), 
                    series = 8, 
                    fold = 2,
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
    conc[i] <- cur/fold
    cur <- conc[i]
    i <- i + 1
  }
  # combine to tibble and return
  odf <- tibble(compound, conc, curveLab = paste0(prefix,c(series:1)))
  return(odf)
}

get_all_conc <- function(compounds,starting_cons,series=8,fold=2){
  
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
    int <- anyFold(starting_cons[i],compound=compounds[i],series=series,fold=fold)
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
               peakarea=ifelse(peakarea < na.value,1,peakarea))
      
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
           peakarea=ifelse(peakarea < na.value,1,peakarea),
           compound_name=ifelse(grepl("[0-9]",compound_name),
                                gsub("^X","",compound_name),compound_name))
  
  meta <- int
  return(meta)
}

readin_meta_csv_single_file_tms <- function(filename,na.value=0,recursive=F){
  
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
    mutate(conc=ifelse(grepl("[Ss][Pp][Ll][Ii][Tt]\\d.+",sampleid),"split1to5",
                       ifelse(grepl("\\_[Nn][Oo][Ss][Pp][Ll][Ii][Tt]",sampleid),"nosplit",
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
           peakarea=ifelse(peakarea < na.value,1,peakarea),
           compound_name=ifelse(grepl("[0-9]",compound_name),
                                gsub("^X","",compound_name),compound_name),
           compound_name = gsub("\\.", "_", compound_name))
  
  meta <- int
  return(meta)
}

make_norm_conc_tbl <- function(df,
                               dil_compounds = "Valerate,Phenol,Valine_D8,Proline_D7, Succinate",
                               conc_compounds = "Phenol"){
  
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
    summarize(avg=mean(peakarea),
              med=median(peakarea))
  
  return(int_conc)
  
}

# make_norm_conc_heatmap <- function(df,
#                                    dil_compounds = "Valerate,Phenol,Valine_D8,Proline_D7, Succinate",
#                                    conc_compounds = "Phenol"){
#   
#   #make df with average for conc and diluted
#   #dil_compounds=c("Valine_D8","Valerate")
#   #conc_compounds=c("Succinate","Proline","Phenol")
#   #dil_compounds <- "Valine_D8,Valerate"
#   
#   dil_compounds <- trimws(unlist(strsplit(dil_compounds,split=",")))
#   conc_compounds <- trimws(unlist(strsplit(conc_compounds,split=",")))
#   
#   int_conc <- tibble(compound_name = c(dil_compounds, conc_compounds),
#                      conc = c(rep("diluted",length(dil_compounds)),
#                               rep("concentrated",length(conc_compounds)))) %>% 
#     inner_join(df) %>%
#     filter(itsd=="ITSD") %>%
#     # dplyr::count(compound_name, conc) %>%
#     # dplyr::count(sampleid) %>%
#     group_by(compound_name) %>%
#     summarize(min_peak = ifelse(all(peakarea == 0), 0, min(peakarea[peakarea != 0])))
#   
#   return(int_conc)
#   
# }


get_indole_conc <- function(conc, compounds,series=11){
  #compounds <- inputcompounds2
  compounds <- "niacin,tyrosine,phenylalanine,kynurenine,Serotonin,anthranilicacid,tryptophan,5HIAA,tryptamine,kynurenicacid,melatonin"
  # series <- 11
  # conc <- c(909,454.5,227.25,113.625,56.8125,5.68125,0.568125,0.056813,0.014203,0.003551,0.000888)
  
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

readin_bile_csv_single_file <- function(filename,na.value=0,recursive=F){
  
  #filename="bile_acid_test.csv"
  #na.value=0
  
  int <- read.csv(file=filename) %>%
    select(-CAS.ID) %>%
    reshape2::melt(id.vars=c("Compound.Name","Formula",
                             "Mass","RT")) %>%
    replace_na(list(value=na.value)) %>%
    mutate(letter=str_extract(Compound.Name,"\\_[A-Z]+$"),
           letter=gsub("\\_","",letter),
           itsd=str_extract(Compound.Name,pattern="ITSD"),
           com=gsub("\\_ITSD","",Compound.Name),
           Data.File=variable) %>%
    separate(com,into=c("num","compound_name","letter"),sep="\\_") %>%
    separate(variable,into=c("num2","date_run","batch","sampleid","conc"),sep="\\_\\_") %>%
    mutate(num2 = gsub("[Xx]", "", num2),
           sampleid = paste(num2, sampleid, sep = "_")) %>% 
    select(Data.File,sampleid,date_run,Compound.Name,compound_name,
           batch,letter,itsd,conc,peakarea=value) %>%
    mutate(filename=filename,
           compound_name=gsub("D[0-9]+\\-","",compound_name),
           compound_name=tolower(compound_name),
           conc=ifelse(grepl("dil",conc),"diluted","concentrated"),
           peakarea=as.numeric(peakarea),
           peakarea=ifelse(peakarea < na.value,1,peakarea))
  
  
  return(int)
}

# This function is used to generate heights and widths of PDFs for faceted plots
get_row_col <- function(p) {
  n <- length(unique(ggplot_build(p)$data[[1]]$PANEL))
  par <- ggplot_build(p)$layout$facet$params
  wrap_dims(n, par$nrow, par$ncol)
}

# set up directory and files ----------------------------------------------

#system("open .")
wddir <- "/Volumes/chaubard-lab/shiny_workspace/csvs/"

# ui ----------------------------------------------------------------------

ui <- fluidPage(
  shinythemes::themeSelector(),
  shinytheme("journal"),
  titlePanel("DFI Metabolomics QC (v1.8.11)"),
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
              # PFBBr quant QC UI -------------------------------------------------------------
              tabPanel(HTML("PFBBr<br/>Quant QC"), theme = shinytheme("flatly"),
                       sidebarLayout(
                         sidebarPanel(width = 3,
                                      textInput(inputId = "compounds",
                                                label = "Enter ITSD compounds (comma separated):",
                                                value = "Acetate,Propionate,Butyrate,Succinate"),
                                      br(),
                                      textInput(inputId = "quant_conc",
                                                label = "Quant con/dil:",
                                                value = "dil,dil,dil,conc"),
                                      h4("ITSD dilution calculation"),
                                      numericInput("xfactor","Mult factor:",value = 1),
                                      #numericInput("start","Enter concentration(s):",value = 100),
                                      textInput("start","Enter concentration(s):","100,25,12.5,50"),
                                      numericInput("series","dilution #",10),
                                      br(),
                                      h4("Filters:"),
                                      textInput("maxcc","Max conc(s) filter:","100,25,12.5,50"),
                                      textInput("mincc","Min conc(s) filter:","0,0,0,0"),
                                      br(),
                                      br(),
                                      numericInput("quant_zero_val",
                                                   "Minimum value:",
                                                   value=1000),
                         ),
                         mainPanel(
                           br(),
                           br(),
                           dataTableOutput("conc"),
                           plotOutput("quant"),
                           h4("Calibration Curve Peak Areas"),
                           dataTableOutput("calib_table"),
                           h4("Fitted linear model stats:"),
                           dataTableOutput("model"), # Slope table
                           br(),
                           actionButton("cc_metrics_download", "Store Calibration Curve Metrics"),
                           tags$style(type="text/css", 
                                      "#cc_metrics_download {background-color:green;color: white}"),
                           br(),
                           br(),
                           downloadButton("qc_report_download", "Download PFBBr QC Quant Report", class = "butt"),
                           tags$style(type="text/css", 
                                      "#qc_report_download {background-color:green;color: white}"),
                           br(),
                           br(),
                           plotOutput("plasma_plot",
                                      height = "750px",
                                      width = "100%"),
                           h4("Standardized table:"),
                           dataTableOutput("quant_tbl"),
                           actionButton("downloadData", "Download PFBBr Quant Table"),
                           tags$style(type="text/css", 
                                      "#downloadData {background-color:green;color: white}"),
                           br(),
                           br(),
                           h4("Quantitative Results:"),
                           checkboxInput("qcfil_quant", "Remove QCs"),
                           downloadButton("quant_download", "Download PFBBr Quant Barplot"),
                           plotOutput("quant_barplot", 
                                      height = "800px",
                                      width = "150%")
                         )
                       )
              ),
              
              # PFBBr normalization UI --------------------------------------------------------
              
              tabPanel(HTML("PFBBr<br/>Normalization"),fluidPage(theme = shinytheme("flatly")),
                       sidebarLayout(
                         sidebarPanel(width = 3,
                                      br(),
                                      h4("Type in ITSD"),
                                      textInput("dil_compounds",
                                                "diluted standards:", 
                                                value="Valine_D8,Valerate"),
                                      textInput("conc_compounds",
                                                "concentrated standards:",
                                                value="Proline_D7,Phenol"),
                                      br(),
                                      numericInput("zero_val",
                                                   "minimum value:",
                                                   value=1000)
                         ),
                         mainPanel(
                           br(),
                           downloadButton("norm_qc_report_download", "Download PFBBr QC Norm Report", class = "butt"),
                           tags$style(type="text/css", "#norm_qc_report_download {background-color:green;color: white}"),
                           splitLayout(cellWidths = c("25%","75%"),
                                       uiOutput("compound_list"),
                                       plotOutput("raw_boxplots",height="1515px")),
                           h4("Normalized heatmap"),
                           downloadButton("heatmap_download", "Download PFBBr Normalized Heatmap"),
                           h4("This heatmap shows the log2 fold-change of median-normalized peak areas for each compound. Compounds are clustered on the y-axis and samples are clustered on the x-axis."),
                           plotOutput("heatmap_plot",
                                      height = "1000px",
                                      width = "125%"),
                           br(),
                           br(),
                           h4("Intermediate table:"),
                           dataTableOutput("conc_filter"),
                           h4("Normalized table:"),
                           dataTableOutput("normwide"),
                           checkboxInput("qcfil","Remove QCs"),
                           downloadButton("downloadData2", "Download PFBBr Normalized Table")
                         )
                       )
              ),
              
              # Indole quant QC UI------------------------------------------------------------
              tabPanel(HTML("Indole<br/>Quant QC"), theme = shinytheme("flatly"),
                       sidebarLayout(
                         sidebarPanel(width = 3,
                                      textAreaInput("compounds2","Enter ITSD compounds (comma separated):",
                                                    value="5HIAA,anthranilicacid,kynurenicacid,kynurenine,melatonin,niacin,phenylalanine,Serotonin,tryptamine,tryptophan,tyrosine"),
                                      br(),
                                      h4("ITSD dilution calculation"),
                                      # numericInput("xfactor2","Mult factor:",value = 11),
                                      # numericInput("start","Enter concentration(s):",value = 100),
                                      textAreaInput("start2","Enter concentration(s):","909,454.5,227.25,113.625,56.8125,5.68125,0.568125,0.056813,0.014203,0.003551,0.000888"),
                                      numericInput("series2","dilution #",11),
                                      br(),
                                      h4("Filters:"),
                                      textAreaInput("maxcc2","Max conc(s) filter:","909,909,909,909,909,909,909,909,909,909,909"),
                                      textInput("mincc2","Min conc(s) filter:","0,0,0,0,0,0,0,0,0,0,0"),                                                     br(),
                                      br(),
                                      numericInput("quant_zero_val2",
                                                   "Minimum value:",
                                                   value=100),
                         ),
                         mainPanel(
                           br(),
                           br(),
                           dataTableOutput("conc2"),
                           plotOutput("quant2"),
                           # dataTableOutput("TEST"),
                           # plotOutput("TEST_PLOT",
                           #            height = "800px",
                           #            width = "150%"),
                           checkboxInput("sety","set y?",value = 0),
                           numericInput("yint","Y intercept:",value=0),
                           h4("Calibration Curve Peak Areas"),
                           dataTableOutput("calib_table2"),
                           h4("Fitted linear model stats:"),
                           dataTableOutput("model2"),
                           br(),
                           actionButton("indole_cc_metrics_download", "Store Calibration Curve Metrics"),
                           tags$style(type="text/css", "#indole_cc_metrics_download {background-color:green;color: white}"),
                           br(),
                           br(),
                           downloadButton("indole_qc_report_download", "Download Indole QC Quant Report", class = "butt"),
                           tags$style(type="text/css", "#indole_qc_report_download {background-color:green;color: white}"),
                           br(),
                           br(),
                           plotOutput("indole_plasma_plot",
                                      height = "750px",
                                      width = "100%"),
                           h4("Standardized table:"),
                           dataTableOutput("quant_tbl2"),
                           actionButton("downloadData3", "Download Indole Quant Table"),
                           tags$style(type="text/css", "#downloadData3 {background-color:green;color: white}"),
                           br(),
                           br(),
                           h4("Quantitative Results:"),
                           checkboxInput("indole_qcfil_quant", "Remove QCs"),
                           downloadButton("indole_quant_download", "Download Indole Quant Barplot"),
                           plotOutput("indole_quant_barplot", 
                                      height = "800px",
                                      width = "150%")
                         )
                       )
              ),
              
              # Indole normalization UI --------------------------------------------------------
              
              tabPanel(HTML("Indole<br/>Normalization"),fluidPage(theme = shinytheme("flatly")),
                       sidebarLayout(
                         sidebarPanel(width = 3,
                                      br(),
                                      h4("Type in ITSD"),
                                      textInput("itsd_compounds","standards:", value="Serotonin,melatonin"),
                                      #textInput("conc_compounds","concentrated standards:",value="Proline_D7,Phenol"),
                                      br(),
                                      numericInput("zero_val2","minimum value:",value=100)
                         ),
                         mainPanel(
                           br(),
                           downloadButton("indole_norm_qc_report_download", "Download Indole QC Norm Report", class = "butt"),
                           tags$style(type="text/css", "#indole_norm_qc_report_download {background-color:green;color: white}"),
                           splitLayout(#cellWidths = c("25%","75%"),
                             # uiOutput("compound_list2"),
                             plotOutput("raw_boxplots2",height="750px")),
                           # dataTableOutput("TESTMAP"),
                           downloadButton("heatmap_table2", "Download Indole Heatmap Data"),
                           # plotOutput("TEST_PLOT",
                           #            height = "800px",
                           #            width = "150%"),
                           h4("Normalized heatmap"),
                           downloadButton("heatmap_download2", "Download Indole Heatmap"),
                           h4("This heatmap shows the log2 fold-change of median-normalized peak areas for each compound. Compounds are clustered on the y-axis and samples are clustered on the x-axis."),
                           plotOutput("heatmap_plot2",
                                      height = "1000px",
                                      width = "125%"),
                           h4("Intermediate table:"),
                           dataTableOutput("conc_int2"),
                           dataTableOutput("conc_filter2"),
                           h4("Normalized table:"),
                           dataTableOutput("normwide2"),
                           checkboxInput("qcfil2","RemoveQCs"),
                           downloadButton("downloadData4", "Download Normalized Table")
                         )
                       )
              ),
              
              # Bile acid quant ---------------------------------------------------------
              tabPanel(HTML("Bile Acid<br/>Quant QC"), theme = shinytheme("flatly"),
                       sidebarLayout(
                         sidebarPanel(width = 3,
                                      textInput("compounds5","Enter ITSD compounds (comma separated):",
                                                value="Cholic Acid,Deoxycholic Acid,Lithocholic Acid,Glycocholic Acid,Taurocholic Acid,Isodeoxycholic Acid,Alloisolithocholic Acid,3-Oxolithocholic Acid"),
                                      br(),
                                      textInput("quant_conc5","Quant con/dil:",
                                                value="dil,dil,dil,dil,dil,dil,dil,dil"),
                                      h4("ITSD dilution calculation"),
                                      numericInput("xfactor5","Mult factor:",value = 1),
                                      #numericInput("start","Enter concentration(s):",value = 100),
                                      textInput("start5","Enter concentration(s):","125"),
                                      numericInput("series5","dilution #",10),
                                      br(),
                                      h4("Filters:"),
                                      textInput("maxcc5","Max conc(s) filter:","125,125,125,125,125,125,125,125"),
                                      textInput("mincc5","Min conc(s) filter:","0,0,0,0,0,0,0,0"),
                                      numericInput("quant_zero_val_bile_acid",
                                                   "",
                                                   value=3000),
                         ),
                         mainPanel(
                           br(),
                           br(),
                           dataTableOutput("TEST5"),
                           dataTableOutput("conc5"),
                           #dataTableOutput("conc_tbl5"),
                           #dataTableOutput("cutoff_df5"),
                           plotOutput("quant5"),
                           h4("Fitted linear model stats:"),
                           dataTableOutput("model5"),
                           br(),
                           actionButton("ba_cc_metrics_download", "Store Calibration Curve Metrics"),
                           tags$style(type="text/css", "#ba_cc_metrics_download {background-color:green;color: white}"),
                           br(),
                           br(),
                           downloadButton("ba_qc_report_download", "Download Bile Acid QC Quant Report", class = "butt"),
                           tags$style(type="text/css", "#ba_qc_report_download {background-color:green;color: white}"),
                           br(),
                           br(),
                           plotOutput("ba_plasma_plot",
                                      height = "750px",
                                      width = "100%"),
                           h4("Standardized table:"),
                           dataTableOutput("quant_tbl5"),
                           actionButton("downloadData5", "Download Bile Acid Quant Table"),
                           tags$style(type = "text/css", "#downloadData5 {background-color:green;color: white}"),
                           br(),
                           br(),
                           h4("Quantitative Results:"),
                           checkboxInput("bile_qcfil_quant", "Remove QCs"),
                           downloadButton("bile_quant_download", "Download Barplot"),
                           plotOutput("bile_quant_barplot", 
                                      height = "500px",
                                      width = "175%")
                         )
                       )
              ),
              
              
              # Bile acid normalization UI --------------------------------------------------------
              
              tabPanel(HTML("Bile Acid<br/>Normalization"),fluidPage(theme = shinytheme("flatly")),
                       sidebarLayout(
                         sidebarPanel(width = 3,
                                      br(),
                                      h4("Minimum Value"),
                                      # textInput("dil_compounds_bile_acid",
                                      #           "diluted standards:",
                                      #           value=""),
                                      # textInput("conc_compounds_bile_acid",
                                      #           "concentrated standards:",
                                      #           value=""),
                                      numericInput("zero_val_bile_acid",
                                                   "",
                                                   value=3000)
                         ),
                         mainPanel(
                           br(),
                           # dataTableOutput("BATEST1"),
                           downloadButton("ba_norm_qc_report_download", "Download Bile Acid QC Norm Report", class = "butt"),
                           tags$style(type="text/css", "#ba_norm_qc_report_download {background-color:green;color: white}"),
                           # plotOutput("TEST_PLOT",
                           #            height = "800px",
                           #            width = "150%"),
                           
                           # dataTableOutput("BATEST2"),
                           splitLayout(cellWidths = c("25%","75%"),
                                       uiOutput("compound_list_bile_acid"),
                                       plotOutput("raw_boxplots_bile_acid",height="1072px")),
                           h4("Normalized heatmap"),
                           downloadButton("heatmap_download_bile_acid", "Download Bile Acid Heatmap"),
                           h4("This heatmap shows the log2 fold-change of median-normalized peak areas for each compound. Compounds are clustered on the y-axis and samples are clustered on the x-axis."),
                           plotOutput("heatmap_plot_bile_acid",
                                      height = "1000px",
                                      width = "125%"),
                           # h4("Conc_Int"),
                           # dataTableOutput("conc_int_bile_acid"),
                           h4("Intermediate table:"),
                           dataTableOutput("conc_filter_bile_acid"),
                           h4("Normalized table:"),
                           dataTableOutput("normwide_bile_acid"),
                           checkboxInput("qcfil_bile_acid","Remove QCs"),
                           downloadButton("downloadData_bile_acid", "Download Bile Acid Normalized Table")
                         )
                       )
              ),
              # TMS normalization UI --------------------------------------------------------
              
              tabPanel(HTML("TMS<br/>Normalization"),fluidPage(theme = shinytheme("flatly")),
                       sidebarLayout(
                         sidebarPanel(width = 3,
                                      br(),
                                      h4("Type in ITSD"),
                                      textInput("compounds_tms",
                                                "standards:",
                                                value="D7_proline,13C_palmitate_329"),
                                      br(),
                                      h4("Type in ITSD"),
                                      br(),
                                      numericInput("zero_val_tms",
                                                   "minimum value:",
                                                   value=1000)
                         ),
                         mainPanel(
                           # dataTableOutput("TMSTEST1"),
                           br(),
                           downloadButton("norm_qc_report_download_tms", "Download TMS QC Norm Report", class = "butt"),
                           tags$style(type="text/css", "#norm_qc_report_download_tms {background-color:green;color: white}"),
                           splitLayout(cellWidths = c("25%","75%"),
                                       uiOutput("compound_list_tms"),
                                       plotOutput("raw_boxplots_tms",height="3555px")),
                           h4("Normalized heatmap"),
                           downloadButton("heatmap_download_tms", "Download TMS Heatmap"),
                           h4("This heatmap shows the log2 fold-change of median-normalized peak areas for each compound. Compounds are clustered on the y-axis and samples are clustered on the x-axis."),
                           plotOutput("heatmap_plot_tms",
                                      height = "1300px",
                                      width = "125%"),
                           h4("Intermediate table:"),
                           dataTableOutput("conc_filter_tms"),
                           h4("Normalized table:"),
                           dataTableOutput("normwide_tms"),
                           checkboxInput("qcfil_tms","Remove QCs"),
                           downloadButton("downloadData_tms", "Download TMS Normalized Table")
                         )
                       )
              ),
              
              
              # Instrument QC UI --------------------------------------------------------
              
              tabPanel(HTML("Instrument<br/>QC"),fluidPage(theme = shinytheme("flatly")),
                       sidebarLayout(
                         sidebarPanel(width = 3,
                                      br(),
                                      h4("Select Quant Method for QC Analysis"),
                                      selectInput("method",
                                                "Quant Method:",
                                                choices = c("PFBBr","Indole","BileAcid")),
                                      br()
                         ),
                         mainPanel(
                         # dataTableOutput("TEST"),
                         h4("Calibration Slopes per Batch"),
                         plotOutput("slope_qc_plot", height="360px", width = "110%"),
                         downloadButton("slope_qc_plot_download", "Download Calibration Plot", class = "butt"),
                         tags$style(type="text/css", "#slope_qc_plot_download {background-color:green;color: white}"),
                         br(),
                         br(),
                         h4("Plasma QC Samples per Batch"),
                         plotOutput("plasma_qc_plot", height = "360px", width = "110%"),
                         downloadButton("plasma_qc_plot_download", "Download Plasma QC Plot", class = "butt"),
                         tags$style(type="text/css", "#plasma_qc_plot_download {background-color:green;color: white}")
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
  
  
  # PFBBr quant tab ------------------------------------------------------------
  
  #make concentration table
  conc_tbl <- reactive({
    get_all_conc(input$start,
                 compounds=input$compounds,
                 series=input$series)
  })
  
  output$conc <- renderDataTable(
    conc_tbl(),
    options = list(pageLength=5)
  )
  
  quant_conc_tbl <- reactive({
    
    compounds = unlist(strsplit(input$compounds, split=","))
    concs <- unlist(strsplit(input$quant_conc, split=","))
    
    int <- cbind(compounds,concs) %>% as.data.frame()
    colnames(int) <- c("compound_name","conc")
    
    int %>%
      mutate(conc=ifelse(conc=="dil","diluted","concentrated"))
    
  })
  
  

  
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
  meta <- reactive({ readin_meta_csv_single_file(file.path(wddir,input$filename),na.value = input$quant_zero_val)})

  #show models and make plots
  modelstart <- reactive({
    
    compounds = unlist(strsplit(input$compounds, split=","))
    compounds = factor(compounds,level=unique(compounds))
    
    model <- meta() %>%
      filter(compound_name %in% compounds) %>%
      inner_join(quant_conc_tbl()) %>%
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
      filter(compound_name %in% compounds) %>%
      inner_join(quant_conc_tbl()) %>%
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
  
  #make calibration peak area table
  calib_table <- reactive({
    
    compounds = unlist(strsplit(input$compounds, split=","))
    compounds = factor(compounds,level=unique(compounds))
    
    meta() %>%
      filter(compound_name %in% compounds) %>%
      inner_join(quant_conc_tbl()) %>%
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
      dplyr::rename(slope_value=conc_val) %>% 
      select(sampleid, compound_name, conc, curveLab, ITSD, peak, norm_peak, slope_value, maxcc, mincc)
  })
  
  output$calib_table <- renderDataTable(
    calib_table() %>%
      datatable() %>%
      formatRound(c(5:ncol(calib_table())), 3) #%>% 
      # formatStyle(columns = c(4:ncol(calib_table())), 'text-align' = 'center')
  )
  
  #make quant table
  quant_table <- reactive({
    
    compounds = unlist(strsplit(input$compounds, split=","))
    
    meta() %>%
      filter(compound_name %in% compounds) %>%
      inner_join(quant_conc_tbl()) %>%
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
      filter(compound_name %in% compounds) %>%
      inner_join(quant_conc_tbl()) %>%
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
  
  
  ## Save quant table with and without QCs
  observeEvent( 
    input$downloadData,
    {
    write.csv(quant_table_dl(), 
              paste0("/Volumes/chaubard-lab/shiny_workspace/CLIN_Finals_QCs/quant_results_",
                                gsub("\\.csv","",input$filename),"_",
                                gsub("\\-","",Sys.Date()),".csv"), row.names=F,quote=F)
    write.csv(quant_table_dl() %>% filter(!str_detect(sampleid, "[Mm][Bb]"),
                                          !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
                                          !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
                                          !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
                                          !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
                                          !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
                                          !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
                                          !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
                                          !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
                                          !grepl("CC[0-9]+", sampleid)), 
              paste0("/Volumes/chaubard-lab/shiny_workspace/CLIN_Finals/removed_qcs_quant_results_",
                     gsub("\\.csv","",input$filename), "_",
                     gsub("\\-","",Sys.Date()),".csv"), row.names=F,quote=F)
    } 
  )
  
  #make quant bar graph
  quant_barplot <- function()({
    
    compounds = unlist(strsplit(input$compounds, split=","))
    
    if(input$qcfil_quant==F){
      
      meta() %>%
        filter(compound_name %in% compounds) %>%
        inner_join(quant_conc_tbl()) %>%
        replace_na(list(itsd="peak")) %>%
        reshape2::dcast(sampleid+compound_name+conc ~ itsd,value.var="peakarea") %>%
        mutate(#peak = ifelse(peak <= 10000,0,peak),
          norm_peak=peak / ITSD) %>%
        filter(!grepl("\\_CC[0-9]",sampleid)) %>%
        left_join(modelstart()) %>%
        mutate(quant_val =  (norm_peak - (`(Intercept)`))/slope_value*as.numeric(input$xfactor)) %>%
        arrange(compound_name) %>%
        mutate(quant_val = ifelse(quant_val < 0,0,quant_val),
               quant_val = round(quant_val,2))%>%
        select(sampleid, compound_name, quant_val) %>% 
        separate(sampleid,into=c("num","date","batch","sampleid","conc"),
                 sep="__") %>% 
        mutate(sampleid_conc = paste(sampleid, conc, sep = "_")) %>% 
        ggplot(aes(x = sampleid, y = quant_val, fill = compound_name)) +
        geom_bar(stat="identity") +
        theme_bw() +
        theme(plot.margin=unit(c(5.5, 15, 5.5, 10),"points"),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              legend.text = element_text(size = 14),
              legend.title = element_text(size = 16),
              legend.position = "top",
              strip.text=element_text(size=18),
              axis.text.y =element_text(size=11),
              axis.title = element_text(size = 15),
              axis.text.x =element_text(vjust = 0.5,
                                        hjust = 1,
                                        angle = 90,
                                        size = ifelse(
                                          (((nrow(quant_bar_data())))+(nrow(quant_bar_data()))) / 
                                            (nrow(quant_bar_data())/10 * nrow(quant_bar_data())/200) >= 11, 
                                          11, 
                                          (((nrow(quant_bar_data())))+(nrow(quant_bar_data()))) / 
                                            (nrow(quant_bar_data())/10 * nrow(quant_bar_data())/200)))) +
        facet_grid(~compound_name, space = "free_x") +
        scale_fill_locuszoom() +
        xlab("\nSampleID") +
        ylab("Concentration (mM)\n") +
        guides(fill = guide_legend(title="Compound    "))
      
    } else{
      
      meta() %>%
        filter(compound_name %in% compounds) %>%
        inner_join(quant_conc_tbl()) %>%
        replace_na(list(itsd="peak")) %>%
        reshape2::dcast(sampleid+compound_name+conc ~ itsd,value.var="peakarea") %>%
        mutate(#peak = ifelse(peak <= 10000,0,peak),
          norm_peak=peak / ITSD) %>%
        filter(!grepl("\\_CC[0-9]",sampleid)) %>%
        left_join(modelstart()) %>%
        mutate(quant_val =  (norm_peak - (`(Intercept)`))/slope_value*as.numeric(input$xfactor)) %>%
        arrange(compound_name) %>%
        mutate(quant_val = ifelse(quant_val < 0,0,quant_val),
               quant_val = round(quant_val,2))%>%
        select(sampleid, compound_name, quant_val) %>% 
        separate(sampleid,into=c("num","date","batch","sampleid","conc"),
                 sep="__") %>% 
        mutate(sampleid_conc= paste(sampleid, conc, sep = "_")) %>% 
        filter(!str_detect(sampleid, "[Mm][Bb]"),
               !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
               !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
               !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
               !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
               !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
               !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
               !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
               !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
               !grepl("CC[0-9]+", sampleid)
        ) %>% 
        ggplot(aes(x = sampleid, y = quant_val, fill = compound_name)) +
        geom_bar(stat="identity") +
        theme_bw() +
        theme(plot.margin=unit(c(5.5, 15, 5.5, 10),"points"),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              legend.position = "top",
              legend.text = element_text(size = 14),
              legend.title = element_text(size = 16),
              strip.text=element_text(size=18),
              axis.text.y =element_text(size=11),
              axis.title = element_text(size = 15),
              axis.text.x =element_text(vjust = 0.5,
                                        hjust = 1,
                                        angle = 90,
                                        size = ifelse(
                                          (((nrow(quant_bar_data())))+(nrow(quant_bar_data()))) / 
                                            (nrow(quant_bar_data())/10 * nrow(quant_bar_data())/200) >= 11, 
                                          11, 
                                          (((nrow(quant_bar_data())))+(nrow(quant_bar_data()))) / 
                                            (nrow(quant_bar_data())/10 * nrow(quant_bar_data())/200)))) +
        facet_grid(~compound_name, space = "free_x") +
        scale_fill_locuszoom() +
        xlab("\nSampleID") +
        ylab("Concentration (mM)\n") +
        guides(fill = guide_legend(title="Compound    "))
    }
  })
  
  # Output the bar graph
  output$quant_barplot <- renderPlot(
    quant_barplot()
  )
  
  # Make quant_bar dataframe to assign sizes based on number of samples
  quant_bar_data <- function()({
    compounds = unlist(strsplit(input$compounds, split=","))
    
    if(input$qcfil_quant==F){
      
      meta() %>%
        filter(compound_name %in% compounds) %>%
        #inner_join(quant_conc_tbl()) %>%
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
        select(sampleid, compound_name, quant_val) %>% 
        separate(sampleid,into=c("num","date","batch","sampleid","conc"),
                 sep="__") %>% 
        mutate(sampleid_conc = paste(sampleid, conc, sep = "_")) 
    } else {
      meta() %>%
        filter(compound_name %in% compounds) %>%
        #inner_join(quant_conc_tbl()) %>%
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
        select(sampleid, compound_name, quant_val) %>% 
        separate(sampleid,into=c("num","date","batch","sampleid","conc"),
                 sep="__") %>% 
        mutate(sampleid_conc = paste(sampleid, conc, sep = "_")) %>% 
        filter(!str_detect(sampleid, "[Mm][Bb]"),
               !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
               !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
               !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
               !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
               !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
               !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
               !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
               !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
               !grepl("CC[0-9]+", sampleid))
    }
  })
  
  # Produce file labels based on qcfil_quant
  quant_bar_data_label <- reactive({
    if(input$qcfil_quant==F){
      return("quant_barplot_")
    } else {
      return("removed_qcs_quant_barplot_")
    } 
  })
  
  output$quant_download <- downloadHandler(
    filename = function(){
      paste0(quant_bar_data_label(),gsub(".csv","",input$filename),"_",Sys.Date(),".pdf")
    },
    content = function(file) {
      ggsave(file,plot=quant_barplot(),
             width = nrow(quant_bar_data()) / (nrow(quant_bar_data())*0.025)+1.5,
             height = ncol(quant_bar_data()) / (ncol(quant_bar_data())*0.05))
    }
  )
  
  ##### QC Report ####
  # ITSD Raw Peak Area # 
  # Build joined dataframe to loop over
  meta2_1 <- reactive({
    meta() %>% 
      inner_join(., quant_conc_tbl(), by = c("compound_name", "conc")) %>% 
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>% 
      filter(itsd == "ITSD",
             !str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]")) %>% 
      mutate(peakarea = ifelse(peakarea <= input$zero_val, input$zero_val, peakarea),
             num = as.numeric(num),
             cc_shape = ifelse(grepl("CC[0-9]+", sampleid), "Calibration Curve Sample", "Standard Sample"))
  })
  
  
  # Build another dataframe of summary stats
  meta2_2 <- reactive({
    meta() %>% 
      inner_join(., quant_conc_tbl(), by = c("compound_name", "conc")) %>% 
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>% 
      filter(itsd == "ITSD",
             !str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
             !grepl("[Cc][Cc][0-9]+", sampleid)) %>% 
      mutate(peakarea = ifelse(peakarea <= input$zero_val, input$zero_val, peakarea)) %>% 
      group_by(batch, compound_name) %>% 
      summarise(stdev = sd(peakarea),
                average = mean(peakarea),
                middle = median(peakarea),
                cv = stdev / average,
                cv_med = stdev / median(peakarea)# Don't turn into % here since it will be applied in the y-axis scale
      )
  })
  
  pfbbr_qc_plot1 <-function()({
    meta2_1() %>%
      left_join(meta2_2()) %>% 
      mutate(flag = ifelse(peakarea > average + (1.5 * stdev) | 
                             peakarea < average - (1.5 * stdev), paste(num,batch,sampleid, sep = "_"), NA)) %>% 
      group_by(batch,compound_name) %>%
      ggplot(., aes(x = num, y = peakarea, group = compound_name, label = flag)) +
      geom_point(aes(color = compound_name, shape = cc_shape), fill = "black", alpha = 0.6, size = 2) +
      geom_line(aes(y = average, color = compound_name)) +
      geom_line(aes(y = average + 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_line(aes(y = average - 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_ribbon(aes(ymin = average - stdev, ymax = average + stdev,
                      fill = compound_name), alpha=0.2) +
      ggrepel::geom_label_repel(size = 1.2, max.overlaps = Inf, 
                                min.segment.length = 0.1, label.padding = 0.1) +
      theme_bw() +
      theme(panel.grid.minor= element_blank(),
            panel.grid.major.x = element_blank(),
            legend.text = element_text(color = "black", size = 8),
            legend.title = element_text(color = "black", size = 10),
            legend.position = "top",
            strip.text=element_text(color = "black", size=9),
            axis.text =element_text(color = "black", size=8),
            axis.title = element_text(color = "black", size = 14),
            plot.margin = margin(1,1,0,1.1, unit = 'cm')) +
      ggsci::scale_color_locuszoom() +
      ggsci::scale_fill_locuszoom() +
      guides(color = guide_legend(title = "Internal Standard Compound",
                                  override.aes = list(size = 2.5), nrow = 1,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right"), fill = F,
             shape = guide_legend(title = "",
                                  override.aes = list(size = 2.5), nrow = 2,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right")) +
      scale_shape_manual(values = c(24,16))+
      scale_y_continuous(label = scales::scientific) +
      scale_x_continuous(breaks = seq(0,150,25)) +
      ylab("Raw Peak Area\n") +
      xlab("\nInjection Number")+
      facet_grid(~compound_name, scales="free_x")
  })
  # output$TEST_PLOT <- renderPlot(
  #   pfbbr_qc_plot2()
  # )
  
  pfbbr_qc_plot2 <- function()({
    
    meta2_1() %>%
      left_join(meta2_2()) %>% 
      mutate(flag = ifelse(peakarea > average + (1.5 * stdev) | 
                             peakarea < average - (1.5 * stdev), paste(num,batch,sampleid, sep = "_"), NA)) %>% 
      group_by(batch,compound_name) %>%
      ggplot(., aes(x = num, y = peakarea, group = compound_name, label = flag)) +
      geom_point(aes(color = compound_name, shape = cc_shape), fill = "black", alpha = 0.6, size = 2) +
      geom_line(aes(y = average, color = compound_name)) +
      geom_line(aes(y = average + 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_line(aes(y = average - 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_ribbon(aes(ymin = ifelse(average - stdev <= 1, 1, average - stdev), 
                      ymax = average + stdev,
                      fill = compound_name), alpha=0.2) +
      ggrepel::geom_label_repel(size = 1.2, max.overlaps = Inf, 
                                min.segment.length = 0.1, label.padding = 0.1) +
      theme_bw() +
      theme(panel.grid.minor= element_blank(),
            panel.grid.major.x = element_blank(),
            legend.text = element_text(color = "black", size = 8),
            legend.title = element_text(color = "black", size = 10),
            legend.position = "top",
            strip.text=element_text(color = "black", size=9),
            axis.text =element_text(color = "black", size=8),
            axis.title = element_text(color = "black", size = 14),
            plot.margin = margin(1,1,0,1.1, unit = 'cm')) +
      ggsci::scale_color_locuszoom() +
      ggsci::scale_fill_locuszoom() +
      guides(color = guide_legend(title = "Internal Standard Compound",
                                  override.aes = list(size = 2.5), nrow = 1,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right"), fill = F,
             shape = guide_legend(title = "",
                                  override.aes = list(size = 2.5), nrow = 2,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right")) +
      scale_shape_manual(values = c(24,16))+
      ylab("Raw Peak Area\n(log10 scale)") +
      xlab("\nInjection Number") +
      # coord_trans(y = "log10") +
      scale_y_continuous(trans = "log10", labels = scales::scientific)+
      scale_x_continuous(breaks = seq(0,150,25))+
      facet_grid(~compound_name, scales="free_x")
  })
  
  # Show models and make plots
  model2 <-  reactive({
    
    compounds = unlist(strsplit(input$compounds, split=","))
    compounds = factor(compounds,level=unique(compounds))
    
    meta() %>%
      inner_join(quant_conc_tbl()) %>%
      mutate(compound_name=factor(compound_name,levels=compounds)) %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(sampleid+compound_name+conc ~ itsd,value.var="peakarea",
                      fun.aggregate=mean) %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      mutate(ITSD = ifelse(ITSD <= input$zero_val, input$zero_val, ITSD),
             norm_peak = peak / ITSD,
             curveLab = sampleid) %>%
      filter(grepl(pattern = "CC", sampleid)) %>%
      left_join(conc_tbl()) %>%
      left_join(cutoff_df()) %>%
      filter(conc_val <= maxcc,
             conc_val >= mincc) %>%
      group_by(batch,compound_name,conc) %>%
      summarize(r = cor(norm_peak,conc_val),
                model_list <- broom::tidy(lm(norm_peak ~ conc_val))) %>%
      reshape2::dcast(batch+compound_name+conc+r ~ term,value.var="estimate") %>%
      dplyr::rename(Batch = batch,
                    `Internal Standard` = compound_name,
                    Concentration = conc,
                    Intercept = `(Intercept)`,
                    Slope = conc_val)
  })
  
  # Plot CC models
  pfbbr_qc_plot3 <- function()({
    
    compounds = unlist(strsplit(input$compounds, split=","))
    compounds = factor(compounds,level=unique(compounds))
    
    meta() %>%
      inner_join(quant_conc_tbl()) %>%
      mutate(compound_name=factor(compound_name,levels=compounds)) %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(sampleid+compound_name+conc ~ itsd,value.var="peakarea",
                      fun.aggregate=mean) %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      mutate(ITSD = ifelse(ITSD <= input$zero_val, input$zero_val, ITSD),
             norm_peak = peak / ITSD,
             curveLab = sampleid,
             conc = ifelse(conc == "conc","Concentrated Standard", "Diluted Standard")) %>%
      filter(grepl(pattern = "CC", sampleid)) %>%
      left_join(conc_tbl()) %>%
      left_join(cutoff_df()) %>%
      filter(conc_val <= maxcc,
             conc_val >= mincc) %>%
      group_by(compound_name, conc) %>%
      ggplot(aes(x=conc_val,y=norm_peak)) +
      geom_point(size=3) +
      geom_smooth(method="lm", formula = y~x) +
      ggpmisc::stat_poly_eq(formula = y~x,
                            rr.digits = 5,
                            coef.digits = 5,
                            aes(label = paste(..eq.label..,
                                              ..rr.label..,
                                              sep = "~~~")),
                            parse = TRUE,
                            size = 3.7) +
      theme(strip.text=element_text(size=13, color = "black"),
            axis.text=element_text(size=13, color = "black"),
            plot.title = element_text(hjust = 0.5, size = 16, color = "black", face = "bold"),
            plot.margin = margin(1,0.5,0.5,0.5, unit = 'cm')) +
      facet_wrap(~compound_name+conc,scales="free") +
      xlab("\nConcentration (mM)") +
      ylab("Normalized Peak Area\n") +
      ggtitle("Calibration Curves")
  })
  
  ### ITSD CV Percent ###
  
  # Build CV dataframe
  meta2_3 <- reactive({
    meta() %>% 
      inner_join(., quant_conc_tbl(), by = c("compound_name", "conc")) %>% 
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>% 
      filter(itsd == "ITSD",
             !grepl(pattern = "^CC[0-9]+", sampleid),
             !str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]")) %>% 
      mutate(peakarea = ifelse(peakarea <= input$zero_val, input$zero_val, peakarea)) %>% 
      group_by(batch, compound_name) %>% 
      summarise(stdev = sd(peakarea),
                average = mean(peakarea),
                middle = median(peakarea),
                cv = stdev / average,
                cv_med = stdev / median(peakarea)
      ) %>% 
      dplyr::rename(Batch = batch,
                    `Internal Standard` = compound_name,
                    StDev = stdev,
                    Mean = average,
                    Median = middle,
                    CV = cv,
                    `CV Median` = cv_med) %>% 
      mutate(StDev = round(StDev, digits = 0),
             Mean = round(Mean, digits = 0),
             Median = round(Median, digits = 0),
             `CV (%)` = round(CV * 100, digits = 1),
             `CV Median (%)` = round(`CV Median` * 100, digits = 1)) %>%
      select(Batch,`Internal Standard`,StDev,Mean,Median,`CV (%)`,`CV Median (%)`)
  })
  
  # Build summary table
  tt <- gridExtra::ttheme_default(
    colhead=list(fg_params=list(col="black", fontface=2L)),
    padding = unit(c(0.5,0.75), "cm"))
  
  
  table_list1 <- function()({
    gridExtra::tableGrob(meta2_3(), rows = NULL, theme = tt)
  })
  
  # output$TEST <- renderDataTable(
  #   meta2_3(),
  #   options = list(pageLength=5)
  # )
  
  ### Plasma QC ###
  
  # Build summary specifically for high and low ranges for horizontal lines
  meta2_4 <- reactive({
    compounds = unlist(strsplit(input$compounds, split=","))
    
    meta() %>%
      filter(compound_name %in% compounds) %>%
      inner_join(quant_conc_tbl()) %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(sampleid+compound_name+conc ~ itsd,value.var="peakarea") %>%
      mutate(norm_peak=peak / ITSD) %>%
      filter(!grepl("\\_CC[0-9]",sampleid),
             grepl(pattern = "[Pp][Ll][Aa][Ss][Mm][Aa]", sampleid)) %>%
      left_join(modelstart()) %>%
      mutate(quant_val =  (norm_peak - (`(Intercept)`))/slope_value*as.numeric(input$xfactor)) %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>% 
      arrange(compound_name) %>%
      mutate(quant_val = ifelse(quant_val < 0,0,quant_val),
             quant_val = round(quant_val,2)) %>% 
        group_by(batch, compound_name) %>%
        summarise(stdev = sd(quant_val, na.rm = T),
                  average = mean(quant_val, na.rm = T),
                  cv = (stdev / average)*100,
                  upper_limit = average + 2.5*stdev,
                  lower_limit = average - 2.5*stdev
        ) %>%
      pivot_longer(-c(batch, compound_name), names_to = "variable", values_to = "value") %>%
      mutate(variable = factor(variable, levels = c("upper_limit", "average", "lower_limit","stdev","cv")))
    })
  
  # Need another summary dataframe to plot CVs and averages for geom_label
  meta2_5 <- reactive({
    compounds = unlist(strsplit(input$compounds, split=","))
    
    meta() %>%
      filter(compound_name %in% compounds) %>%
      inner_join(quant_conc_tbl()) %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(sampleid+compound_name+conc ~ itsd,value.var="peakarea") %>%
      mutate(norm_peak=peak / ITSD) %>%
      filter(!grepl("\\_CC[0-9]",sampleid),
             grepl(pattern = "[Pp][Ll][Aa][Ss][Mm][Aa]", sampleid)) %>%
      left_join(modelstart()) %>%
      mutate(quant_val =  (norm_peak - (`(Intercept)`))/slope_value*as.numeric(input$xfactor)) %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>% 
      arrange(compound_name) %>%
      mutate(quant_val = ifelse(quant_val < 0,0,quant_val),
             quant_val = round(quant_val,2)) %>% 
      group_by(batch, compound_name) %>%
      summarise(stdev = sd(quant_val, na.rm = T),
                average = mean(quant_val, na.rm = T),
                cv = (stdev / average)*100,
                upper_limit = average + 2.5*stdev,
                lower_limit = average - 2.5*stdev
      ) 
  })

  # Build another dataframe to plot
  meta2_6 <- reactive({
    compounds = unlist(strsplit(input$compounds, split=","))
    
    meta() %>%
      filter(compound_name %in% compounds) %>%
      inner_join(quant_conc_tbl()) %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(sampleid+compound_name+conc ~ itsd,value.var="peakarea") %>%
      mutate(norm_peak=peak / ITSD) %>%
      filter(!grepl("\\_CC[0-9]",sampleid),
             grepl(pattern = "[Pp][Ll][Aa][Ss][Mm][Aa]", sampleid)) %>%
      left_join(modelstart()) %>%
      mutate(quant_val =  (norm_peak - (`(Intercept)`))/slope_value*as.numeric(input$xfactor)) %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>% 
      separate(sampleid, c("sample","qc","replicate"), sep = "_") %>% 
      arrange(compound_name) %>%
      mutate(quant_val = ifelse(quant_val < 0,0,quant_val),
             quant_val = round(quant_val,2))
  })
  
  # output$TEST <- renderDataTable(
  #   meta2_6(),
  #   options = list(pageLength=5)
  # )

  plasma_plot = reactive({
    ggplot(meta2_6(), aes(x = batch, y = quant_val, shape = replicate, color = compound_name, 
                          label = ifelse(is.na(quant_val),"NA",paste0(round(quant_val, digits = 1), "mM")))) +
      geom_point(size = 3) +
      geom_point(meta2_5(), mapping = aes(x= batch, y = average),
                 color = "black", shape = 3, size = 3, inherit.aes = F) +
      geom_hline(subset(meta2_4(), variable %in%
                          c("average","upper_limit","lower_limit")),
                 mapping = aes(yintercept = value, linetype = variable))+
      ggrepel::geom_label_repel(meta2_5(),
                                mapping = aes(x= batch, y = average,
                                              label = ifelse(is.na(cv),"NA",paste0(round(cv, digits = 1), "%"))),
                                label.padding = 0.1,
                                direction = "y", max.overlaps = 50, min.segment.length = 5,
                                size = 1.2, inherit.aes = F) +
      ggrepel::geom_label_repel(label.padding = 0.1,
                                direction = "x", max.overlaps = 50, min.segment.length = 5,
                                size = 1.2, show.legend = FALSE) +
      theme_bw() +
      theme(panel.grid.minor= element_blank(),
            panel.grid.major.x = element_blank(),
            legend.text = element_text(color = "black", size = 8),
            legend.title = element_text(color = "black", size = 10, hjust = 0.5),
            legend.position = "top",
            plot.title = element_text(color = "black", size = 18, face = "bold", hjust= 0.5),
            strip.text=element_text(color = "black", size=9),
            axis.text.x =element_text(color = "black", size=8),
            axis.title = element_text(color = "black", size = 14),
            plot.margin = margin(2,0.5,2,1.1, unit = 'cm')) +
      ggsci::scale_color_locuszoom() +
      guides(color = guide_legend(title = "Compound",
                                  override.aes = list(size = 2.5), ncol = 1,
                                  title.position="top",
                                  label.position = "right"), fill = F,
             shape = guide_legend(title = "Replicate",
                                  override.aes = list(size = 2.5), ncol = 1,
                                  title.position="top",
                                  label.position = "right"),
             linetype = guide_legend(ncol = 1,
                                     title.position="top",
                                     label.position = "right")) +
      scale_linetype_manual(name = "Metrics", values = c(1, 2, 1),
                            labels = c("+2.5 SD", "Mean", "-2.5 SD"),
                            guide = guide_legend(ncol = 1,
                                                 title.position="top",
                                                 label.position = "right"))+
      xlab("\nBatch Number") +
      ylab("Concentration (mM)\n") +
      ggtitle("Plasma QC Summary\n") +
      facet_wrap(~compound_name, nrow = 1)
  })
  
  output$plasma_plot <- renderPlot(
    plasma_plot()
  )

  
  ## Save report 
  output$qc_report_download <- downloadHandler(
    filename = function(){
      paste0("PFBBr_QC_Report_",unique(meta2_2()$batch),"_",Sys.Date(),".pdf")
    },
    content = function(file) {
      pdf(file, height = 11, width = 8.5, onefile = TRUE)
      print(
        gridExtra::grid.arrange(
          pfbbr_qc_plot1() +
            xlab("") +
            ggtitle(paste("PFBBr Quantitative QC Report\n", unique(meta2_1()$batch))) +
            theme(plot.title = element_text(color = "black",
                                            hjust = 0.5,
                                            size = 20,
                                            face = "bold")),
          pfbbr_qc_plot2() + theme(legend.position = "none"),
          nrow = 2)
      )
      print(plasma_plot())
      print(
        gridExtra::grid.arrange(
          table_list1())
      )
      print(pfbbr_qc_plot3())
      dev.off()
    }
  )
  
  ## Save model metrics in csv
  observeEvent( 
    input$cc_metrics_download,
    {write.csv(model2(), paste0("/Volumes/chaubard-lab/shiny_workspace/calibration_metrics/",
                                gsub(".csv","",input$filename),"_CC_Metrics",".csv"))} 
  )


  # PFBBr normalization tab -------------------------------------------------------

  rawdf <- reactive({
    readin_meta_csv_single_file(file.path(wddir,input$filename),na.value = input$zero_val)
  })
  
  output$dfcheck <- renderDataTable(
    rawdf()
  )

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

  # conc_int_heatmap <- reactive({
  # 
  #   make_norm_conc_heatmap(rawdf(),
  #                          dil_compounds = as.character(input$dil_compounds),
  #                          conc_compounds = as.character(input$conc_compounds))
  # })

  conc_int_sep <- reactive({
    conc_int() %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="__")
  })


  output$conc_int <- renderDataTable(
    conc_int()
  )

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

  #make heatmap:
  heatmap_plot <- function()({
    rawdf() %>%
      left_join(conc_filter()) %>%
      filter(conc==checked, is.na(itsd),
             !grepl("(__CC[0-9]+__)", sampleid)) %>%
      left_join(conc_int()) %>%
      # left_join(conc_int_heatmap()) %>%
      mutate(norm_peak = peakarea / avg) %>%
      group_by(compound_name) %>% 
      mutate(compound_med = median(norm_peak)) %>% 
      ungroup() %>% 
      mutate(heat_val = log((norm_peak / compound_med), base = 2)) %>%
      dplyr::select(sampleid, compound_name, heat_val) %>%
      spread(compound_name, heat_val, fill = NA) %>%
      reshape2::melt(id.vars=c("sampleid")) %>%
      dplyr::rename(compound_name=variable,heat_val=value) %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="\\_\\_") %>%
      filter(!is.na(heat_val)) %>%
      mutate(sampleid = ifelse(str_detect(sampleid, "[Mm][Bb]|[Pp][Oo][Oo][Ll][Ee][Dd]|[Bb][Hh][Ii][Qq][Cc]|[Pp][Ll][Aa][Ss][Mm][Aa]|[Hh][Ee][Xx][Aa][Nn][Ee][Ss]|[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]|50%_[Mm][Ee][Oo][Hh]|[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]|50%[Mm][Ee][Oo][hh]"),
                               paste(num, sampleid, conc, sep = "__"),
                               sampleid)) %>%
      dplyr::select(sampleid, compound_name, heat_val) %>%
      mutate(heat_val = round(heat_val,5)) %>%
      group_by(sampleid, compound_name) %>%
      summarise(heat_val = mean(heat_val)) %>%
      pivot_wider(names_from = compound_name, values_from = heat_val) %>%
      filter(!str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]")) %>%
      drop_na(.) %>%
      column_to_rownames(., var = "sampleid") %>%
      as.matrix(.) %>%
      t(.) %>%
      pheatmap(., scale = "none",
               cellheight = nrow(heatmap_data()) / (nrow(heatmap_data())*0.075),
               cellwidth = ncol(heatmap_data()) / (ncol(heatmap_data())*0.0825)+1.5,
               angle_col = "90",
               color=colorRampPalette(c("navy", "white", "red"))(50),
      )
  })

  output$heatmap_plot <- renderPlot(
    heatmap_plot()
  )

  
  heatmap_data <- function()({
    rawdf() %>%
      left_join(conc_filter()) %>%
      filter(conc==checked, is.na(itsd),
             !grepl("(__CC[0-9]+__)", sampleid)) %>%
      left_join(conc_int()) %>%
      # left_join(conc_int_heatmap()) %>%
      mutate(norm_peak = peakarea / avg) %>%
      group_by(compound_name) %>% 
      mutate(compound_med = median(norm_peak)) %>% 
      ungroup() %>% 
      mutate(heat_val = log((norm_peak / compound_med), base = 2)) %>%
      dplyr::select(sampleid, compound_name, heat_val) %>%
      spread(compound_name, heat_val, fill = NA) %>%
      reshape2::melt(id.vars=c("sampleid")) %>%
      dplyr::rename(compound_name=variable,heat_val=value) %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="\\_\\_") %>%
      filter(!is.na(heat_val)) %>%
      mutate(sampleid = ifelse(str_detect(sampleid, "[Mm][Bb]|[Pp][Oo][Oo][Ll][Ee][Dd]|[Bb][Hh][Ii][Qq][Cc]|[Pp][Ll][Aa][Ss][Mm][Aa]|[Hh][Ee][Xx][Aa][Nn][Ee][Ss]|[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]|50%_[Mm][Ee][Oo][Hh]|[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]|50%[Mm][Ee][Oo][hh]"),
                               paste(num, sampleid, conc, sep = "__"),
                               sampleid)) %>%
      dplyr::select(sampleid, compound_name, heat_val) %>%
      mutate(heat_val = round(heat_val,5)) %>%
      group_by(sampleid, compound_name) %>%
      summarise(heat_val = mean(heat_val)) %>%
      pivot_wider(names_from = compound_name, values_from = heat_val) %>%
      filter(!str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]")) %>%
      drop_na(.) %>%
      column_to_rownames(., var = "sampleid") %>%
      as.matrix(.) %>%
      t(.)
  })

  output$heatmap_download <- downloadHandler(
    filename = function(){
      paste0("normalized_heatmap_",gsub(".csv","",input$filename),"_",Sys.Date(),".pdf")
    },
    content = function(file) {
      pdf(file,
          height = nrow(heatmap_data()) / (nrow(heatmap_data())*0.075),
          width = (ncol(heatmap_data()) / (ncol(heatmap_data())*0.07))+1.5
      )
      heatmap_plot()
      dev.off()
    },
    contentType = 'PDF'
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
      replace_na(list(checked="concentrated")) %>%
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
        replace_na(list(checked="concentrated")) %>%
        filter(conc==checked, is.na(itsd),
               !grepl("(__CC[0-9]+__)", sampleid)) %>%
        select(-checked) %>%
        left_join(conc_int()) %>%
        mutate(norm_peak=peakarea / avg) %>%
        dplyr::select(sampleid, compound_name, norm_peak) %>%
        spread(compound_name, norm_peak, fill = NA) %>%
        reshape2::melt(id.vars=c("sampleid")) %>%
        dplyr::rename(compound_name=variable,norm_peak=value) %>%
        separate(sampleid,into=c("num","date","batch","sampleid","conc"),
                 sep="__") %>%
        filter(!is.na(norm_peak)) %>%
        mutate(sampleid = ifelse(str_detect(sampleid, "[Mm][Bb]|[Pp][Oo][Oo][Ll][Ee][Dd]|[Bb][Hh][Ii][Qq][Cc]|[Pp][Ll][Aa][Ss][Mm][Aa]|[Hh][Ee][Xx][Aa][Nn][Ee][Ss]|[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]|50%_[Mm][Ee][Oo][Hh]|[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]|50%[Mm][Ee][Oo][hh]"),
                                 paste(num, sampleid, conc, sep = "__"),
                                 sampleid)) %>%
        dplyr::select(sampleid, compound_name, norm_peak) %>%
        mutate(norm_peak = round(norm_peak,5)) %>%
        group_by(sampleid, compound_name) %>%
        summarise(norm_peak = mean(norm_peak)) %>%
        #     # reshape2::dcast(sampleid ~ compound_name, value.var = "norm_peak") %>%
        arrange(compound_name) %>%
        pivot_wider(names_from = compound_name, values_from = norm_peak)
      # select(-num)
    }else{
      rawdf() %>%
        left_join(conc_filter()) %>%
        replace_na(list(checked="concentrated")) %>%
        filter(conc==checked, is.na(itsd)) %>%
        select(-checked) %>%
        left_join(conc_int()) %>%
        mutate(norm_peak=peakarea / avg) %>%
        dplyr::select(sampleid, compound_name, norm_peak) %>%
        spread(compound_name, norm_peak, fill = NA) %>%
        reshape2::melt(id.vars=c("sampleid")) %>%
        dplyr::rename(compound_name=variable,norm_peak=value) %>%
        separate(sampleid,into=c("num","date","batch","sampleid","conc"),
                 sep="__") %>%
        mutate(sampleid = ifelse(str_detect(sampleid, "[Mm][Bb]|[Pp][Oo][Oo][Ll][Ee][Dd]|[Bb][Hh][Ii][Qq][Cc]|[Pp][Ll][Aa][Ss][Mm][Aa]|[Hh][Ee][Xx][Aa][Nn][Ee][Ss]|[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]|50%_[Mm][Ee][Oo][Hh]|[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]|50%[Mm][Ee][Oo][hh]"),
                                 paste(num, sampleid, conc, sep = "__"),
                                 sampleid)) %>%
        filter(!is.na(norm_peak)) %>%
        filter(!str_detect(sampleid, "[Mm][Bb]"),
               !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
               !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
               !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
               !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
               !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
               !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
               !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
               !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
               !grepl("CC[0-9]+", sampleid)
        ) %>%
        dplyr::select(num,sampleid, compound_name, norm_peak) %>%
        mutate(norm_peak = round(norm_peak,5)) %>%
        group_by(sampleid, compound_name) %>%
        summarise(norm_peak = mean(norm_peak)) %>%
        # reshape2::dcast(num+sampleid ~ compound_name,value.var="norm_peak",fun.aggregate=mean) %>%
        arrange(compound_name) %>%
        pivot_wider(names_from = compound_name, values_from = norm_peak)
      # select(-num)
    }
  })

  output$normwide <- renderDataTable(
    norm_wide_tbl() %>%
      datatable() %>%
      formatRound(c(2:ncol(norm_wide_tbl())), 3) %>%
      formatStyle(columns = c(2:ncol(norm_wide_tbl())), 'text-align' = 'center')
  )


  norm_wide_tbl_label <- reactive({
    if(input$qcfil==F){
      return("normalized_results_")
    } else {
      return("removed_qcs_normalized_results_")
    }
  })

  #download handler
  output$downloadData2 <- downloadHandler(

    filename = function(){
      paste0(norm_wide_tbl_label(),input$filename,"_",Sys.Date(),".csv")
    },

    content = function(file) {
      write.csv(norm_wide_tbl(),file,
                row.names=F,quote=F)
    }
  )

  ##### QC Report #####
  # ITSD Raw Peak Area #

  # Build a Norm CV dataframe of summary stats
  rawdf2 <- reactive({
    compounds = c(unlist(strsplit(input$dil_compounds, split=",")),unlist(strsplit(input$conc_compounds, split=",")))
    compounds = factor(compounds,level=unique(compounds))

    rawdf() %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      inner_join(., conc_int_sep(), by = c("num","date","batch","sampleid","conc")) %>%
      filter(itsd == "ITSD",
             compound_name %in% compounds,
             !str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
             !grepl("[Cc][Cc][0-9]+", sampleid)) %>%
      mutate(peakarea = ifelse(peakarea <= input$zero_val, input$zero_val, peakarea)) %>%
      group_by(batch, compound_name, conc) %>%
      summarise(stdev = sd(peakarea),
                average = mean(peakarea),
                middle = median(peakarea),
                cv = stdev / average,
                cv_med = stdev / median(peakarea)# Don't turn into % here since it will be applied in the y-axis scale
      )
  })

  pfbbr_norm_plot1 <-function()({

    compounds = c(unlist(strsplit(input$dil_compounds, split=",")),unlist(strsplit(input$conc_compounds, split=",")))
    compounds = factor(compounds,level=unique(compounds))

    rawdf() %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      filter(itsd == "ITSD",
             compound_name %in% compounds,
             !str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]")) %>%
      mutate(peakarea = ifelse(peakarea <= input$zero_val, input$zero_val, peakarea),
             cc_shape = ifelse(grepl("CC[0-9]+", sampleid), "CC Sample", "ITSD")) %>%
      inner_join(., conc_int_sep(), by = c("num","date","batch","sampleid","conc")) %>%
      filter(itsd == "ITSD") %>%
      left_join(rawdf2()) %>%
      mutate(num = as.numeric(num),
             flag = ifelse(peakarea > average + (1.5 * stdev) |
                             peakarea < average - (1.5 * stdev), paste(num,batch,sampleid, sep = "_"), NA)) %>%
      ggplot(., aes(x = num, y = peakarea, group = compound_name, label = flag)) +
      geom_point(aes(color = compound_name, shape = cc_shape), fill = "black", alpha = 0.6, size = 2) +
      geom_line(aes(y = average, color = compound_name)) +
      geom_line(aes(y = average + 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_line(aes(y = average - 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_ribbon(aes(ymin = average - stdev, ymax = average + stdev,
                      fill = compound_name), alpha=0.2) +
      ggrepel::geom_label_repel(size = 1.2, max.overlaps = Inf,
                                min.segment.length = 0.1, label.padding = 0.1) +
      theme_bw() +
      theme(panel.grid.minor= element_blank(),
            panel.grid.major.x = element_blank(),
            legend.text = element_text(color = "black", size = 13),
            legend.title = element_text(color = "black", size = 15),
            legend.position = "top",
            strip.text=element_text(color = "black", size=9),
            axis.text =element_text(color = "black", size=8),
            axis.title = element_text(color = "black", size = 15),
            plot.margin = margin(1,1,0,1.1, unit = 'cm')) +
      ggsci::scale_color_jco() +
      ggsci::scale_fill_jco() +
      guides(color = guide_legend(title = "Internal Standard Compound",
                                  override.aes = list(size = 2.5), nrow = 1,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right"), fill = F,
             shape = guide_legend(title = "",
                                  override.aes = list(size = 2.5), nrow = 2,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right")) +
      scale_shape_manual(values = c(24,16))+
      scale_y_continuous(label = scales::scientific) +
      scale_x_continuous(breaks = seq(0,150,25)) +
      ylab("Raw Peak Area\n") +
      xlab("\nInjection Number")+
      facet_wrap(~compound_name+conc, scales="free_x", nrow = 2)
  })


  pfbbr_norm_plot2 <-function()({

    compounds = c(unlist(strsplit(input$dil_compounds, split=",")),unlist(strsplit(input$conc_compounds, split=",")))
    compounds = factor(compounds,level=unique(compounds))

    rawdf() %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      filter(itsd == "ITSD",
             compound_name %in% compounds,
             !str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]")) %>%
      mutate(peakarea = ifelse(peakarea <= input$zero_val, input$zero_val, peakarea),
             cc_shape = ifelse(grepl("CC[0-9]+", sampleid), "CC Sample", "ITSD")) %>%
      inner_join(., conc_int_sep(), by = c("num","date","batch","sampleid","conc")) %>%
      filter(itsd == "ITSD") %>%
      left_join(rawdf2()) %>%
      mutate(num = as.numeric(num),
             flag = ifelse(peakarea > average + (1.5 * stdev) |
                             peakarea < average - (1.5 * stdev), paste(num,batch,sampleid, sep = "_"), NA)) %>%
      ggplot(., aes(x = num, y = peakarea, group = compound_name, label = flag)) +
      geom_point(aes(color = compound_name, shape = cc_shape), fill = "black", alpha = 0.6, size = 2) +
      geom_line(aes(y = average, color = compound_name)) +
      geom_line(aes(y = average + 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_line(aes(y = average - 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_ribbon(aes(ymin = ifelse(average - stdev <= 1, 1, average - stdev),
                      ymax = average + stdev,
                      fill = compound_name), alpha=0.2) +
      ggrepel::geom_label_repel(size = 1.2, max.overlaps = Inf,
                                min.segment.length = 0.1, label.padding = 0.1) +
      theme_bw() +
      theme(panel.grid.minor= element_blank(),
            panel.grid.major.x = element_blank(),
            legend.text = element_text(color = "black", size = 13),
            legend.title = element_text(color = "black", size = 15),
            legend.position = "top",
            strip.text=element_text(color = "black", size=9),
            axis.text =element_text(color = "black", size=8),
            axis.title = element_text(color = "black", size = 15),
            plot.margin = margin(1,1,0,1.1, unit = 'cm')) +
      ggsci::scale_color_jco() +
      ggsci::scale_fill_jco() +
      guides(color = guide_legend(title = "Internal Standard Compound",
                                  override.aes = list(size = 2.5), nrow = 1,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right"), fill = F,
             shape = guide_legend(title = "",
                                  override.aes = list(size = 2.5), nrow = 2,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right")) +
      scale_shape_manual(values = c(24,16))+
      ylab("Raw Peak Area\n(log10 scale)") +
      xlab("\nInjection Number") +
      # coord_trans(y = "log10") +
      scale_y_continuous(trans = "log10", labels = scales::scientific)+
      scale_x_continuous(breaks = seq(0,150,25))+
      facet_wrap(~compound_name+conc, scales="free_x", nrow = 2)
  })

  # ITSD CV Percent #

  # Build CV dataframe
  rawdf2_1 <- reactive({
    compounds = c(unlist(strsplit(input$dil_compounds, split=",")),unlist(strsplit(input$conc_compounds, split=",")))
    compounds = factor(compounds,level=unique(compounds))

    rawdf() %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      inner_join(., conc_int_sep(), by = c("num","date","batch","sampleid","conc")) %>%
      filter(itsd == "ITSD",
             compound_name %in% compounds,
             !str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
             !grepl("[Cc][Cc][0-9]+", sampleid)) %>%
      mutate(peakarea = ifelse(peakarea <= input$zero_val, input$zero_val, peakarea)) %>% 
      group_by(batch, compound_name, conc) %>%
      summarise(stdev = sd(peakarea),
                average = mean(peakarea),
                middle = median(peakarea),
                cv = stdev / average,
                cv_med = stdev / median(peakarea)# Don't turn into % here since it will be applied in the y-axis scale
      ) %>%
      dplyr::rename(Batch = batch,
                    `Internal Standard` = compound_name,
                    StDev = stdev,
                    Mean = average,
                    Median = middle,
                    CV = cv,
                    `CV Median` = cv_med,
                    Concentration = conc) %>%
      mutate(StDev = round(StDev, digits = 0),
             Mean = round(Mean, digits = 0),
             Median = round(Median, digits = 0),
             `CV (%)` = round(CV * 100, digits = 1),
             `CV Median (%)` = round(`CV Median` * 100, digits = 1)) %>%
      select(Batch,Concentration,`Internal Standard`,StDev,Mean,Median,`CV (%)`,`CV Median (%)`)
  })

  # Build summary table
  norm_table_list1 <- function()({
    gridExtra::tableGrob(rawdf2_1(), rows = NULL, theme = tt)
  })

  ## Save Report ##
  output$norm_qc_report_download <- downloadHandler(
    filename = function(){
      paste0("PFBBr_QC_Norm_Report_",unique(rawdf2()$batch),"_",Sys.Date(),".pdf")
    },
    content = function(file) {
      pdf(file, height = 11, width = 8.5)
      print(pfbbr_norm_plot1() +
            xlab("") +
            ggtitle(paste("PFBBr Qualitative QC Report\n", unique(rawdf2()$batch))) +
            theme(plot.title = element_text(color = "black",
                                            hjust = 0.5,
                                            size = 20,
                                            face = "bold")))
      print(pfbbr_norm_plot2() + theme(legend.position = "none"))
      print(
        gridExtra::grid.arrange(
          norm_table_list1())
      )
      dev.off()
    }
  )
  


  # Indole quant tab ------------------------------------------------------------


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
  output$conc2 <- renderDataTable(
    conc_tbl2(),
    options = list(pageLength=5)
  )


  # read in table as reactive
  indole_meta <- reactive({ readin_meta_csv_single_file(file.path(wddir,input$filename), na.value = input$quant_zero_val2) })

  #show models and make plots
  modelstart2 <- reactive({

    compounds2 = unlist(strsplit(input$compounds2, split=","))
    compounds2 = factor(compounds2,levels=compounds2)

    if(input$sety==0){
      indole_meta() %>%
        mutate(compound_name=gsub("\\_[0-9]+$","",compound_name)) %>%
        filter(compound_name %in% compounds2) %>%
        mutate(compound_name=factor(compound_name,levels=compounds2)) %>%
        replace_na(list(itsd="peak")) %>%
        reshape2::dcast(Data.File+compound_name+conc ~ itsd,value.var="peakarea") %>%
        mutate(#peak = ifelse(peak <= 10000,0,peak),
          norm_peak=peak / ITSD) %>%
        filter(grepl("\\_[Cc][Cc][1-9][0-9]+\\_|\\_[Cc][Cc][1-9]+\\_", Data.File)) %>%
        mutate(curveLab=str_extract(Data.File,"\\_[Cc][Cc][1-9][0-9]+\\_|\\_[Cc][Cc][1-9]+\\_"),
               curveLab=gsub("\\_","",curveLab),
               curveLab = tolower(curveLab)) %>%
        left_join(conc_tbl2()) %>%
        left_join(cutoff_df2()) %>%
        filter(conc_val <= maxcc,
               conc_val >= mincc) %>%
        group_by(compound_name) %>%
        summarize(r = cor(norm_peak,conc_val),
                  model_list <- broom::tidy(lm(norm_peak ~ conc_val))) %>%
        reshape2::dcast(compound_name+r ~ term,value.var="estimate") %>%
        dplyr::rename(slope_value=conc_val)
    }else{
      indole_meta() %>%
        mutate(compound_name=gsub("\\_[0-9]+$","",compound_name)) %>%
        filter(compound_name %in% compounds2) %>%
        mutate(compound_name=factor(compound_name,levels=compounds2)) %>%
        replace_na(list(itsd="peak")) %>%
        reshape2::dcast(Data.File+compound_name+conc ~ itsd,value.var="peakarea") %>%
        mutate(#peak = ifelse(peak <= 10000,0,peak),
          norm_peak=peak / ITSD) %>%
        filter(grepl("\\_[Cc][Cc][1-9][0-9]+\\_|\\_[Cc][Cc][1-9]+\\_", Data.File)) %>%
        mutate(curveLab=str_extract(Data.File,"\\_[Cc][Cc][1-9][0-9]+\\_|\\_[Cc][Cc][1-9]+\\_"),
               curveLab=gsub("\\_","",curveLab),
               curveLab = tolower(curveLab)) %>%
        left_join(conc_tbl2()) %>%
        left_join(cutoff_df2()) %>%
        filter(conc_val <= maxcc,
               conc_val >= mincc) %>%
        group_by(compound_name) %>%
        summarize(r = cor(norm_peak,conc_val),
                  model_list <- broom::tidy(lm(as.formula(paste("norm_peak ~ conc_val +",input$yint))))) %>%
        reshape2::dcast(compound_name+r ~ term,value.var="estimate") %>%
        dplyr::rename(slope_value=conc_val) %>%
        mutate(`(Intercept)`= input$yint)

    }

  })

  output$model2 <- renderDataTable(
    modelstart2() %>%
      datatable() %>%
      formatRound(c(2:4), 3) %>%
      formatStyle(columns = c(2:4), 'text-align' = 'center')
  )

  #make calibration peak area table
  calib_table2 <- reactive({
    
    compounds2 = unlist(strsplit(input$compounds2, split=","))
    compounds2 = factor(compounds2,levels=compounds2)
    
    indole_meta() %>%
      mutate(compound_name=gsub("\\_[0-9]+$","",compound_name)) %>%
      filter(compound_name %in% compounds2) %>%
      mutate(compound_name=factor(compound_name,levels=compounds2)) %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(Data.File+compound_name+conc ~ itsd,value.var="peakarea") %>%
      mutate(#peak = ifelse(peak <= 10000,0,peak),
        norm_peak=peak / ITSD) %>%
      filter(grepl("\\_[Cc][Cc][1-9][0-9]+\\_|\\_[Cc][Cc][1-9]+\\_", Data.File)) %>%
      mutate(curveLab=str_extract(Data.File,"\\_[Cc][Cc][1-9][0-9]+\\_|\\_[Cc][Cc][1-9]+\\_"),
             curveLab=gsub("\\_","",curveLab),
             curveLab = tolower(curveLab)) %>%
      left_join(conc_tbl2()) %>%
      left_join(cutoff_df2()) %>%
      filter(conc_val <= maxcc,
             conc_val >= mincc) %>%
      dplyr::rename(slope_value=conc_val,
                    sampleid = Data.File) %>% 
      select(sampleid, compound_name, conc, curveLab, ITSD, peak, norm_peak, slope_value, maxcc, mincc)
  })
  
  output$calib_table2 <- renderDataTable(
    calib_table2() %>%
      datatable() %>%
      formatRound(c(5:ncol(calib_table2())), 3) #%>% 
    # formatStyle(columns = c(4:ncol(calib_table2())), 'text-align' = 'center')
  )

  #make quant graph
  quant_plot2 <- reactive({

    compounds2 = unlist(strsplit(input$compounds2, split=","))
    compounds2 = factor(compounds2,levels=compounds2)

    indole_meta() %>%
      mutate(compound_name=gsub("\\_[0-9]+$","",compound_name)) %>%
      filter(compound_name %in% compounds2) %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(Data.File+compound_name+conc ~ itsd,value.var="peakarea") %>%
      mutate(#peak = ifelse(peak <= 10000,0,peak),
        norm_peak=peak / ITSD) %>%
      filter(grepl("\\_[Cc][Cc][1-9][0-9]+\\_|\\_[Cc][Cc][1-9]+\\_", Data.File)) %>%
      mutate(curveLab=str_extract(Data.File,"\\_[Cc][Cc][1-9][0-9]+\\_|\\_[Cc][Cc][1-9]+\\_"),
             curveLab=gsub("\\_","",curveLab),
             curveLab = tolower(curveLab)) %>%
      left_join(conc_tbl2()) %>%
      left_join(cutoff_df2()) %>%
      filter(conc_val <= maxcc,
             conc_val >= mincc) %>%
      mutate(compound_name=factor(compound_name,levels=compounds2)) %>%
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

    indole_meta() %>%
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
      dplyr::rename(sampleid = Data.File)
  })

  output$quant_tbl2 <- renderDataTable(
    quant_table2() %>%
      datatable() %>%
      formatRound(c(4:ncol(quant_table2())), 3) %>%
      formatStyle(columns = c(4:ncol(quant_table2())), 'text-align' = 'center')
  )

  quant_table_dl2 <- reactive({

    compounds2 = unlist(strsplit(input$compounds2, split=","))
    compounds2 = factor(compounds2,levels=compounds2)

    indole_meta() %>%
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

## Save quant table with and without QCs
  observeEvent( 
    input$downloadData3,
    {
      write.csv(quant_table_dl2(), 
                paste0("/Volumes/chaubard-lab/shiny_workspace/CLIN_Finals_QCs/quant_results_",
                       gsub("\\.csv","",input$filename),"_",
                       gsub("\\-","",Sys.Date()),".csv"), row.names=F,quote=F)
      write.csv(quant_table_dl2() %>% filter(!str_detect(sampleid, "[Mm][Bb]"),
                                             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
                                             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
                                             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
                                             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
                                             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
                                             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
                                             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
                                             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
                                            !grepl("CC[0-9]+", sampleid)), 
                paste0("/Volumes/chaubard-lab/shiny_workspace/CLIN_Finals/removed_qcs_quant_results_",
                       gsub("\\.csv","",input$filename), "_",
                       gsub("\\-","",Sys.Date()),".csv"), row.names=F,quote=F)
    } 
  )
  
  
  #make quant bar graph
  indole_quant_barplot <- function()({
    
    compounds2 = unlist(strsplit(input$compounds2, split=","))
    
    if(input$indole_qcfil_quant==F){
      
      compounds2 = unlist(strsplit(input$compounds2, split=","))
      compounds2 = factor(compounds2,levels=compounds2)
      
      indole_meta() %>%
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
        separate(Data.File,into=c("num","date","batch","sampleid","conc"),
                 sep="\\_\\_") %>% 
        mutate(sampleid_conc = paste(sampleid, conc, sep = "_")) %>% 
        ggplot(aes(x = sampleid, y = quant_val, fill = compound_name)) +
        geom_bar(stat="identity") +
        theme_bw() +
        theme(plot.margin=unit(c(5.5, 15, 5.5, 10),"points"),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              legend.text = element_text(size = 14),
              legend.title = element_text(size = 16),
              legend.position = "top",
              strip.text=element_text(size=18),
              axis.text.y =element_text(size=11),
              axis.title = element_text(size = 15),
              axis.text.x =element_text(vjust = 0.5,
                                        hjust = 1,
                                        angle = 90,
                                        size = ifelse(
                                          (((nrow(indole_quant_bar_data())))+(nrow(indole_quant_bar_data()))) / 
                                            (nrow(indole_quant_bar_data())/10 * nrow(indole_quant_bar_data())/200) >= 11, 
                                          11, 
                                          (((nrow(indole_quant_bar_data())))+(nrow(indole_quant_bar_data()))) / 
                                            (nrow(indole_quant_bar_data())/10 * nrow(indole_quant_bar_data())/200)))) +
        facet_grid(~compound_name, space = "free_x") +
        scale_fill_igv() +
        xlab("\nSampleID") +
        ylab("Concentration (uM)\n") +
        guides(fill = guide_legend(title="Compound    "))
      
    } else{
      
      indole_meta() %>%
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
        separate(Data.File,into=c("num","date","batch","sampleid","conc"),
                 sep="\\_\\_") %>% 
        mutate(sampleid_conc= paste(sampleid, conc, sep = "_")) %>% 
        filter(!str_detect(sampleid, "[Mm][Bb]"),
               !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
               !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
               !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
               !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
               !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
               !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
               !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
               !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
               !grepl("CC[0-9]+", sampleid)
        ) %>% 
        ggplot(aes(x = sampleid, y = quant_val, fill = compound_name)) +
        geom_bar(stat="identity") +
        theme_bw() +
        theme(plot.margin=unit(c(5.5, 15, 5.5, 10),"points"),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              legend.position = "top",
              legend.text = element_text(size = 14),
              legend.title = element_text(size = 16),
              strip.text=element_text(size=18),
              axis.text.y =element_text(size=11),
              axis.title = element_text(size = 15),
              axis.text.x =element_text(vjust = 0.5,
                                        hjust = 1,
                                        angle = 90,
                                        size = ifelse(
                                          (((nrow(indole_quant_bar_data())))+(nrow(indole_quant_bar_data()))) / 
                                            (nrow(indole_quant_bar_data())/10 * nrow(indole_quant_bar_data())/200) >= 11, 
                                          11, 
                                          (((nrow(indole_quant_bar_data())))+(nrow(indole_quant_bar_data()))) / 
                                            (nrow(indole_quant_bar_data())/10 * nrow(indole_quant_bar_data())/200)))) +
        facet_grid(~compound_name, space = "free_x") +
        scale_fill_igv() +
        xlab("\nSampleID") +
        ylab("Concentration (uM)\n") +
        guides(fill = guide_legend(title="Compound    "))
    }
  })
  
  # Output the bar graph
  output$indole_quant_barplot <- renderPlot(
    indole_quant_barplot()
  )
  
  # Make quant_bar dataframe to assign sizes based on number of samples
  indole_quant_bar_data <- function()({
    compounds2 = unlist(strsplit(input$compounds2, split=","))
    
    if(input$indole_qcfil_quant==F){
      
      indole_meta() %>%
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
        separate(Data.File,into=c("num","date","batch","sampleid","conc"),
                 sep="\\_\\_") %>% 
        mutate(sampleid_conc = paste(sampleid, conc, sep = "_")) 
    } else {
      indole_meta() %>%
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
        separate(Data.File,into=c("num","date","batch","sampleid","conc"),
                 sep="\\_\\_") %>% 
        mutate(sampleid_conc = paste(sampleid, conc, sep = "_")) %>% 
        filter(!str_detect(sampleid, "[Mm][Bb]"),
               !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
               !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
               !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
               !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
               !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
               !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
               !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
               !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
               !grepl("CC[0-9]+", sampleid))
    }
  })
  
  # Produce file labels based on qcfil_quant
  indole_quant_bar_data_label <- reactive({
    if(input$indole_qcfil_quant==F){
      return("quant_barplot_")
    } else {
      return("removed_qcs_quant_barplot_")
    } 
  })
  
  output$indole_quant_download <- downloadHandler(
    filename = function(){
      paste0(indole_quant_bar_data_label(),gsub(".csv","",input$filename),"_",Sys.Date(),".pdf")
    },
    content = function(file) {
      ggsave(file,plot=indole_quant_barplot(),
             width = nrow(indole_quant_bar_data()) / (nrow(indole_quant_bar_data())*0.025)+1.5,
             height = ncol(indole_quant_bar_data()) / (ncol(indole_quant_bar_data())*0.05))
    }
  )

  #### QC Report ####
  indole_meta2_1 <- reactive({

    compounds2 = unlist(strsplit(input$compounds2, split=","))
    compounds2 = factor(compounds2,levels=compounds2)

    indole_meta() %>%
      filter(compound_name %in% compounds2) %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      filter(itsd == "ITSD",
             !str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]")) %>%
      mutate(peakarea = ifelse(peakarea <= input$zero_val2, input$zero_val2, peakarea),
             num = gsub("X","",num),
             num = as.numeric(num),
             cc_shape = ifelse(grepl("cc[0-9]+", sampleid), "Calibration Curve Sample", "Standard Sample"))
  })


  # Build another dataframe of summary stats
  indole_meta2_2 <- reactive({

    compounds2 = unlist(strsplit(input$compounds2, split=","))
    compounds2 = factor(compounds2,levels=compounds2)

    indole_meta() %>%
      filter(compound_name %in% compounds2) %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      filter(itsd == "ITSD",
             !str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
             !grepl("[Cc][Cc][0-9]+", sampleid)) %>%
      mutate(peakarea = ifelse(peakarea <= input$zero_val2, input$zero_val2, peakarea),
             num = gsub("X","",num),
             num = as.numeric(num)) %>%
      group_by(batch, compound_name) %>%
      summarise(stdev = sd(peakarea),
                average = mean(peakarea),
                middle = median(peakarea),
                cv = stdev / average,
                cv_med = stdev / median(peakarea))
  })


  indole_qc_plot1 <-function()({

    indole_meta2_1() %>%
      left_join(as.data.frame(indole_meta2_2())) %>%
      mutate(flag = ifelse(peakarea > average + (1.5 * stdev) |
                             peakarea < average - (1.5 * stdev), paste(num,batch,sampleid, sep = "_"), NA)) %>%
      group_by(compound_name) %>%
      ggplot(., aes(x = num, y = peakarea, group = compound_name, label = flag)) +
      geom_point(aes(color = compound_name, shape = cc_shape), fill = "black", alpha = 0.6, size = 2) +
      geom_line(aes(y = average, color = compound_name)) +
      geom_line(aes(y = average + 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_line(aes(y = average - 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_ribbon(aes(ymin = average - stdev, ymax = average + stdev,
                      fill = compound_name), alpha=0.2) +
      ggrepel::geom_label_repel(size = 1.2, max.overlaps = Inf,
                                min.segment.length = 0.1, label.padding = 0.1) +
      theme_bw() +
      theme(panel.grid.minor= element_blank(),
            panel.grid.major.x = element_blank(),
            legend.text = element_text(color = "black", size = 8),
            legend.title = element_text(color = "black", size = 10),
            legend.position = "top",
            strip.text=element_text(color = "black", size=5),
            axis.text =element_text(color = "black", size=8),
            axis.title = element_text(color = "black", size = 14),
            plot.margin = margin(1,0.5,0,0.6, unit = 'cm')) +
      ggsci::scale_color_igv() +
      ggsci::scale_fill_igv() +
      guides(color = guide_legend(title = "Internal Standard Compound",
                                  override.aes = list(size = 2.5), nrow = 2,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right"), fill = F,
             shape = guide_legend(title = "",
                                  override.aes = list(size = 2.5), nrow = 2,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right")) +
      scale_shape_manual(values = c(24,16))+
      scale_y_continuous(label = scales::scientific) +
      scale_x_continuous(breaks = seq(0,150,25)) +
      ylab("Raw Peak Area\n") +
      xlab("\nInjection Number")+
      facet_wrap(~compound_name, scales="free_x", nrow = 3)
  })


  indole_qc_plot2 <- function()({
    indole_meta2_1() %>%
      left_join(as.data.frame(indole_meta2_2())) %>%
      mutate(flag = ifelse(peakarea > average + (1.5 * stdev) |
                             peakarea < average - (1.5 * stdev), paste(num,batch,sampleid, sep = "_"), NA)) %>%
      group_by(compound_name) %>%
      ggplot(., aes(x = num, y = peakarea, group = compound_name, label = flag)) +
      geom_point(aes(color = compound_name, shape = cc_shape), fill = "black", alpha = 0.6, size = 2) +
      geom_line(aes(y = average, color = compound_name)) +
      geom_line(aes(y = average + 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_line(aes(y = average - 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_ribbon(aes(ymin = ifelse(average - stdev <= 1, 1, average - stdev),
                      ymax = average + stdev,
                      fill = compound_name), alpha=0.2) +
      ggrepel::geom_label_repel(size = 1.2, max.overlaps = Inf,
                                min.segment.length = 0.1, label.padding = 0.1) +
      theme_bw() +
      theme(panel.grid.minor= element_blank(),
            panel.grid.major.x = element_blank(),
            legend.text = element_text(color = "black", size = 8),
            legend.title = element_text(color = "black", size = 10),
            legend.position = "top",
            strip.text=element_text(color = "black", size=5),
            axis.text =element_text(color = "black", size=8),
            axis.title = element_text(color = "black", size = 14),
            plot.margin = margin(3.5,0.5,0,0.6, unit = 'cm')) +
      ggsci::scale_color_igv() +
      ggsci::scale_fill_igv() +
      guides(color = guide_legend(title = "Internal Standard Compound",
                                  override.aes = list(size = 2.5), nrow = 2,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right"), fill = F,
             shape = guide_legend(title = "",
                                  override.aes = list(size = 2.5), nrow = 2,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right")) +
      scale_shape_manual(values = c(24,16))+
      ylab("Raw Peak Area\n(log10 scale)") +
      xlab("\nInjection Number") +
      # coord_trans(y = "log10") +
      scale_y_continuous(trans = "log10", labels = scales::scientific)+
      scale_x_continuous(breaks = seq(0,150,25))+
      facet_wrap(~compound_name, scales="free_x", nrow = 3)
  })


  # Build model dataframe
  indole_model <-  reactive({

    compounds2 = unlist(strsplit(input$compounds2, split=","))
    compounds2 = factor(compounds2,levels=compounds2)

    if(input$sety==0){
      indole_meta() %>%
        mutate(compound_name=gsub("\\_[0-9]+$","",compound_name)) %>%
        filter(compound_name %in% compounds2) %>%
        mutate(compound_name=factor(compound_name,levels=compounds2)) %>%
        separate(sampleid,into=c("num","date","batch","sampleid","conc"), sep="__") %>%
        replace_na(list(itsd="peak")) %>%
        reshape2::dcast(Data.File+batch+compound_name+conc ~ itsd,value.var="peakarea") %>%
        mutate(#peak = ifelse(peak <= 10000,0,peak),
          norm_peak=peak / ITSD) %>%
        filter(grepl("\\_[Cc][Cc][1-9][0-9]+\\_|\\_[Cc][Cc][1-9]+\\_", Data.File)) %>%
        mutate(curveLab=str_extract(Data.File,"\\_[Cc][Cc][1-9][0-9]+\\_|\\_[Cc][Cc][1-9]+\\_"),
               curveLab=gsub("\\_","",curveLab),
               curveLab = tolower(curveLab),
               conc = ifelse(conc == "conc","Concentrated Standard", "Diluted Standard")) %>%
        left_join(conc_tbl2()) %>%
        left_join(cutoff_df2()) %>%
        filter(conc_val <= maxcc,
               conc_val >= mincc) %>%
        group_by(batch,compound_name, conc) %>%
        summarize(r = cor(norm_peak,conc_val),
                  model_list <- broom::tidy(lm(norm_peak ~ conc_val))) %>%
        reshape2::dcast(batch+compound_name+conc+r ~ term,value.var="estimate") %>%
        dplyr::rename(Batch = batch,
                      `Internal Standard` = compound_name,
                      Concentration = conc,
                      Intercept = `(Intercept)`,
                      Slope = conc_val)

    }else{
      indole_meta() %>%
        mutate(compound_name=gsub("\\_[0-9]+$","",compound_name)) %>%
        filter(compound_name %in% compounds2) %>%
        mutate(compound_name=factor(compound_name,levels=compounds2)) %>%
        separate(sampleid,into=c("num","date","batch","sampleid","conc"), sep="__") %>%
        replace_na(list(itsd="peak")) %>%
        reshape2::dcast(Data.File+batch+compound_name+conc ~ itsd,value.var="peakarea") %>%
        mutate(#peak = ifelse(peak <= 10000,0,peak),
          norm_peak=peak / ITSD) %>%
        filter(grepl("\\_[Cc][Cc][1-9][0-9]+\\_|\\_[Cc][Cc][1-9]+\\_", Data.File)) %>%
        mutate(curveLab=str_extract(Data.File,"\\_[Cc][Cc][1-9][0-9]+\\_|\\_[Cc][Cc][1-9]+\\_"),
               curveLab=gsub("\\_","",curveLab),
               curveLab = tolower(curveLab),
               conc = ifelse(conc == "conc","Concentrated Standard", "Diluted Standard")) %>%
        left_join(conc_tbl2()) %>%
        left_join(cutoff_df2()) %>%
        filter(conc_val <= maxcc,
               conc_val >= mincc) %>%
        group_by(batch,compound_name, conc) %>%
        summarize(r = cor(norm_peak,conc_val),
                  model_list <- broom::tidy(lm(norm_peak ~ conc_val))) %>%
        reshape2::dcast(batch+compound_name+conc+r ~ term,value.var="estimate") %>%
        dplyr::rename(Batch = batch,
                      `Internal Standard` = compound_name,
                      Concentration = conc,
                      Intercept = `(Intercept)`,
                      Slope = conc_val)

    }
  })


  # Show models and make plots
  indole_qc_plot3 <-  reactive({

      compounds2 = unlist(strsplit(input$compounds2, split=","))
      compounds2 = factor(compounds2,levels=compounds2)

      if(input$sety==0){
        indole_meta() %>%
          mutate(compound_name=gsub("\\_[0-9]+$","",compound_name)) %>%
          filter(compound_name %in% compounds2) %>%
          mutate(compound_name=factor(compound_name,levels=compounds2)) %>%
          separate(sampleid,into=c("num","date","batch","sampleid","conc"), sep="__") %>%
          replace_na(list(itsd="peak")) %>%
          reshape2::dcast(Data.File+batch+compound_name+conc ~ itsd,value.var="peakarea") %>%
          mutate(#peak = ifelse(peak <= 10000,0,peak),
            norm_peak=peak / ITSD) %>%
          filter(grepl("\\_[Cc][Cc][1-9][0-9]+\\_|\\_[Cc][Cc][1-9]+\\_", Data.File)) %>%
          mutate(curveLab=str_extract(Data.File,"\\_[Cc][Cc][1-9][0-9]+\\_|\\_[Cc][Cc][1-9]+\\_"),
                 curveLab=gsub("\\_","",curveLab),
                 curveLab = tolower(curveLab),
                 conc = ifelse(conc == "conc","Concentrated Standard", "Diluted Standard")) %>%
          left_join(conc_tbl2()) %>%
          left_join(cutoff_df2()) %>%
          filter(conc_val <= maxcc,
                 conc_val >= mincc) %>%
          group_by(batch,compound_name, conc) %>%
          ggplot(aes(x=conc_val,y=norm_peak)) +
          geom_point(size=3) +
          geom_smooth(method="lm", formula = y~x) +
          ggpmisc::stat_poly_eq(formula = y~x,
                                rr.digits = 5,
                                coef.digits = 5,
                                size = 3.2,
                                aes(label = paste(..eq.label..,
                                                  ..rr.label..,
                                                  sep = "~~~~")),
                                parse = TRUE) +
          theme(strip.text=element_text(size=9, color = "black"),
                axis.text=element_text(size=8, color = "black"),
                axis.title = element_text(size = 11, color = "black"),
                plot.title = element_text(hjust = 0.5, size = 16,
                                          color = "black", face = "bold"),
                plot.margin = margin(1,0.5,0.5,0.5, unit = 'cm')) +
          facet_wrap(~compound_name+conc,scales="free", ncol = 2) +
          xlab("\nConcentration (uM)") +
          ylab("Normalized Peak Area\n") +
          ggtitle("Calibration Curves")
      }else{
        indole_meta() %>%
          mutate(compound_name=gsub("\\_[0-9]+$","",compound_name)) %>%
          filter(compound_name %in% compounds2) %>%
          mutate(compound_name=factor(compound_name,levels=compounds2)) %>%
          separate(sampleid,into=c("num","date","batch","sampleid","conc"), sep="__") %>%
          replace_na(list(itsd="peak")) %>%
          reshape2::dcast(Data.File+batch+compound_name+conc ~ itsd,value.var="peakarea") %>%
          mutate(#peak = ifelse(peak <= 10000,0,peak),
            norm_peak=peak / ITSD) %>%
          filter(grepl("\\_[Cc][Cc][1-9][0-9]+\\_|\\_[Cc][Cc][1-9]+\\_", Data.File)) %>%
          mutate(curveLab=str_extract(Data.File,"\\_[Cc][Cc][1-9][0-9]+\\_|\\_[Cc][Cc][1-9]+\\_"),
                 curveLab=gsub("\\_","",curveLab),
                 curveLab = tolower(curveLab),
                 conc = ifelse(conc == "conc","Concentrated Standard", "Diluted Standard")) %>%
          left_join(conc_tbl2()) %>%
          left_join(cutoff_df2()) %>%
          filter(conc_val <= maxcc,
                 conc_val >= mincc) %>%
          group_by(batch,compound_name, conc) %>%
          ggplot(aes(x=conc_val,y=norm_peak)) +
          geom_point(size=3) +
          geom_smooth(method="lm", formula = y~x) +
          ggpmisc::stat_poly_eq(formula = y~x,
                                rr.digits = 5,
                                coef.digits = 5,
                                size = 3.2,
                                aes(label = paste(..eq.label..,
                                                  ..rr.label..,
                                                  sep = "~~~~")),
                                parse = TRUE) +
          theme(strip.text=element_text(size=9, color = "black"),
                axis.text=element_text(size=8, color = "black"),
                axis.title = element_text(size = 11, color = "black"),
                plot.title = element_text(hjust = 0.5, size = 16,
                                          color = "black", face = "bold"),
                plot.margin = margin(1,0.5,0.5,0.5, unit = 'cm')) +
          facet_wrap(~compound_name+conc,scales="free", ncol = 2) +
          xlab("\nConcentration (uM)") +
          ylab("Normalized Peak Area\n") +
          ggtitle("Calibration Curves")

      }

    })

  ## ITSD CV Percent ##

  # Build CV dataframe
  indole_meta2_3 <- reactive({

      compounds2 = unlist(strsplit(input$compounds2, split=","))
      compounds2 = factor(compounds2,levels=compounds2)

      indole_meta() %>%
        filter(compound_name %in% compounds2) %>%
        separate(sampleid,into=c("num","date","batch","sampleid","conc"),
                 sep="__") %>%
        filter(itsd == "ITSD",
               !str_detect(sampleid, "[Mm][Bb]"),
               !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
               !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
               !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
               !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
               !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
               !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
               !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
               !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
               !grepl("[Cc][Cc][0-9]+", sampleid)) %>%
        mutate(peakarea = ifelse(peakarea <= input$zero_val2, input$zero_val2, peakarea),
               num = gsub("X","",num),
               num = as.numeric(num)) %>%
        group_by(batch, compound_name) %>%
        summarise(stdev = sd(peakarea),
                  average = mean(peakarea),
                  middle = median(peakarea),
                  cv = stdev / average,
                  cv_med = stdev / median(peakarea)) %>%
      dplyr::rename(Batch = batch,
                    `Internal Standard` = compound_name,
                    StDev = stdev,
                    Mean = average,
                    Median = middle,
                    CV = cv,
                    `CV Median` = cv_med) %>%
      mutate(StDev = round(StDev, digits = 0),
             Mean = round(Mean, digits = 0),
             Median = round(Median, digits = 0),
             `CV (%)` = round(CV * 100, digits = 1),
             `CV Median (%)` = round(`CV Median` * 100, digits = 1)) %>%
      select(Batch,`Internal Standard`,StDev,Mean,Median,`CV (%)`,`CV Median (%)`)
  })

  # Build summary table
  indole_tt <- gridExtra::ttheme_default(
    colhead=list(fg_params=list(col="black", fontface=2L)),
    padding = unit(c(0.5,0.75), "cm"))

  indole_table_list1 <- function()({
    gridExtra::tableGrob(indole_meta2_3(), rows = NULL, theme = indole_tt)
  })


  ### Plasma QC ###
  
  # Build summary specifically for high and low ranges for horizontal lines
  indole_meta2_4 <- reactive({
    
    compounds2 = unlist(strsplit(input$compounds2, split=","))
    
    indole_meta() %>%
      mutate(compound_name=gsub("\\_[0-9]+$","",compound_name)) %>%
      filter(compound_name %in% compounds2) %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(Data.File+compound_name ~ itsd,value.var="peakarea") %>%
      mutate(#peak = ifelse(peak <= 10000,0,peak),
        norm_peak=peak / ITSD) %>%
      filter(!grepl("\\_[Cc][Cc][0-9]",Data.File),
             grepl(pattern = "[Pp][Ll][Aa][Ss][Mm][Aa]", Data.File)) %>%
      left_join(modelstart2()) %>%
      mutate(quant_val =  (norm_peak - (`(Intercept)`))/slope_value) %>%
      separate(Data.File,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      arrange(compound_name) %>%
      mutate(quant_val = ifelse(quant_val < 0,0,quant_val),
             quant_val = round(quant_val,2)) %>%
      group_by(batch, compound_name) %>%
      summarise(stdev = sd(quant_val, na.rm = T),
                average = mean(quant_val, na.rm = T),
                cv = (stdev / average)*100,
                upper_limit = average + 2.5*stdev,
                lower_limit = average - 2.5*stdev
      ) %>%
      pivot_longer(-c(batch, compound_name), names_to = "variable", values_to = "value") %>%
      mutate(variable = factor(variable, levels = c("upper_limit", "average", "lower_limit","stdev","cv")))
  })
  
  # Need another summary dataframe to plot CVs and averages for geom_label
  indole_meta2_5 <- reactive({
    compounds2 = unlist(strsplit(input$compounds2, split=","))
    
    indole_meta() %>%
      mutate(compound_name=gsub("\\_[0-9]+$","",compound_name)) %>%
      filter(compound_name %in% compounds2) %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(Data.File+compound_name ~ itsd,value.var="peakarea") %>%
      mutate(#peak = ifelse(peak <= 10000,0,peak),
        norm_peak=peak / ITSD) %>%
      filter(!grepl("\\_[Cc][Cc][0-9]",Data.File),
             grepl(pattern = "[Pp]lasma", Data.File)) %>%
      left_join(modelstart2()) %>%
      mutate(quant_val =  (norm_peak - (`(Intercept)`))/slope_value) %>%
      separate(Data.File,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>% 
      arrange(compound_name) %>%
      mutate(quant_val = ifelse(quant_val < 0,0,quant_val),
             quant_val = round(quant_val,2)) %>%
      group_by(batch, compound_name) %>%
      summarise(stdev = sd(quant_val, na.rm = T),
                average = mean(quant_val, na.rm = T),
                cv = (stdev / average)*100,
                upper_limit = average + 2.5*stdev,
                lower_limit = average - 2.5*stdev
      )
  })
  
  # Build another dataframe to plot
  indole_meta2_6 <- reactive({
    compounds2 = unlist(strsplit(input$compounds2, split=","))
    
    indole_meta() %>%
      mutate(compound_name=gsub("\\_[0-9]+$","",compound_name)) %>%
      filter(compound_name %in% compounds2) %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(Data.File+compound_name ~ itsd,value.var="peakarea") %>%
      mutate(#peak = ifelse(peak <= 10000,0,peak),
        norm_peak=peak / ITSD) %>%
      filter(!grepl("\\_[Cc][Cc][0-9]",Data.File),
             grepl(pattern = "[Pp][Ll][Aa][Ss][Mm][Aa]", Data.File)) %>%
      left_join(modelstart2()) %>%
      mutate(quant_val =  (norm_peak - (`(Intercept)`))/slope_value) %>%
      separate(Data.File,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>% 
      mutate(sample = str_extract(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             qc = str_extract(sampleid, "[Qq][Cc]"),
             replicate = as.factor(str_extract(sampleid, "[0-9]+"))) %>% 
      arrange(compound_name) %>%
      mutate(quant_val = ifelse(quant_val < 0,0,quant_val),
             quant_val = round(quant_val,2))
  })
  
  # output$TEST <- renderDataTable(
  #   meta2_6(),
  #   options = list(pageLength=5)
  # )
  
  indole_plasma_plot = reactive({
    ggplot(indole_meta2_6(), aes(x = batch, y = quant_val, shape = replicate, color = compound_name, 
                          label = ifelse(is.na(quant_val),"NA",paste0(round(quant_val, digits = 1), "uM")))) +
      geom_point(size = 3) +
      geom_point(indole_meta2_5(), mapping = aes(x= batch, y = average),
                 color = "black", shape = 3, size = 3, inherit.aes = F) +
      geom_hline(subset(indole_meta2_4(), variable %in%
                          c("average","upper_limit","lower_limit")),
                 mapping = aes(yintercept = value, linetype = variable))+
      ggrepel::geom_label_repel(indole_meta2_5(),
                                mapping = aes(x= batch, y = average,
                                              label = ifelse(is.na(cv),"NA",paste0(round(cv, digits = 1), "%"))),
                                label.padding = 0.1,
                                direction = "y", max.overlaps = 50, min.segment.length = 5,
                                size = 1.2, inherit.aes = F) +
      ggrepel::geom_label_repel(label.padding = 0.1,
                                direction = "both", max.overlaps = 50, min.segment.length = 10,
                                force = 2, nudge_x = 0.3,
                                size = 1.2, show.legend = FALSE) +
      theme_bw() +
      theme(panel.grid.minor= element_blank(),
            panel.grid.major.x = element_blank(),
            legend.text = element_text(color = "black", size = 8),
            legend.title = element_text(color = "black", size = 10, hjust = 0.5),
            legend.position = "top",
            plot.title = element_text(color = "black", size = 18, face = "bold", hjust= 0.5),
            strip.text=element_text(color = "black", size=9),
            axis.text.x =element_text(color = "black", size=8),
            axis.text.y =element_text(color = "black", size=8, angle = 60, hjust = 0.5, vjust = 0.5),
            axis.title = element_text(color = "black", size = 14),
            plot.margin = margin(2,0.5,1.1,0.5, unit = 'cm')) +
      ggsci::scale_color_igv() +
      guides(color = guide_legend(title = "Compound",
                                  override.aes = list(size = 2.5), ncol = 3,
                                  title.position="top", 
                                  label.position = "right"), fill = F,
             shape = guide_legend(title = "Replicate",
                                  override.aes = list(size = 2.5), ncol = 1,
                                  title.position="top",
                                  label.position = "right"),
             linetype = guide_legend(ncol = 1,
                                     title.position="top",
                                     label.position = "right")) +
      scale_linetype_manual(name = "Metrics", values = c(1, 2, 1),
                            labels = c("+2.5 SD", "Mean", "-2.5 SD"),
                            guide = guide_legend(ncol = 1,
                                                 title.position="top",
                                                 label.position = "right"))+
      xlab("\nBatch Number") +
      ylab("Concentration (uM)\n") +
      ggtitle("Plasma QC Summary\n") +
      facet_wrap(~compound_name, nrow = 3, scales = "free_y")
  })
  
  output$indole_plasma_plot <- renderPlot(
    indole_plasma_plot()
  )
  
  
  ## Save report ##
  output$indole_qc_report_download <- downloadHandler(
    filename = function(){
      paste0("Indole_QC_Report_",unique(indole_meta2_2()$batch),"_",Sys.Date(),".pdf")
    },
    content = function(file) {
      pdf(file, height = 11, width = 8.5)
      print(indole_qc_plot1() +
            xlab("") +
            ggtitle(paste("Indole Quantitative QC Report\n", unique(indole_meta2_2()$batch))) +
            theme(plot.title = element_text(color = "black",
                                            hjust = 0.5,
                                            size = 20,
                                            face = "bold"))
            )

      print(indole_qc_plot2() + theme(plot.title = element_blank()))
      print(indole_plasma_plot())
      print(gridExtra::grid.arrange(
          indole_table_list1())
      )
      print(indole_qc_plot3())
      dev.off()
    }
  )

  ## Save model metrics in csv
  observeEvent(
    input$indole_cc_metrics_download,
    {write.csv(indole_model(), paste0("/Volumes/chaubard-lab/shiny_workspace/calibration_metrics/",
                                  gsub(".csv","",input$filename),"_CC_Metrics",".csv"))}
  )




# Indole normalization tab -------------------------------------------------------

  indole_rawdf <- reactive({
  readin_meta_csv_single_file(file.path(wddir,input$filename),na.value = input$zero_val2)
})

  indole_csv <- reactive({

    indole_vars <- indole_rawdf() %>%
    mutate(compound_name=as.character(compound_name)) %>%
    arrange(desc(compound_name)) %>%
    filter(is.na(itsd)) %>%
    distinct(compound_name) %>%
    `$`(compound_name)

  return(indole_vars)
})

# output$compound_list2 <- renderUI({
#   checkboxGroupInput("check_compounds","Check = diluted",
#                      choices=csv(),
#                      inline=F)
# })

indole_conc_int2 <- reactive({

  make_norm_conc_tbl(indole_rawdf(),
                     conc_compounds = as.character(input$itsd_compounds))
  #conc_compounds = as.character(input$conc_compounds))
})

output$conc_int2 <- renderDataTable(
  indole_conc_int2()
)

indole_conc_int_sep2 <- reactive({
  indole_conc_int2() %>%
    separate(sampleid,into=c("num","date","batch","sampleid","conc"),
             sep="__")
})


#make boxplots:
raw_boxplots2 <- reactive({

  indole_rawdf() %>%
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

#make heatmap:
heatmap_plot2 <- function()({
  indole_rawdf() %>%
    mutate(Data.File=as.character(Data.File)) %>%
    filter(#conc==checked,
      is.na(itsd),
      !grepl("[Cc][Cc][0-9]+", Data.File)) %>%
    left_join(indole_conc_int2()) %>%
    mutate(norm_peak = peakarea / avg) %>%
    group_by(compound_name) %>% 
    mutate(compound_med = median(norm_peak)) %>% 
    ungroup() %>% 
    mutate(heat_val = log((norm_peak / compound_med), base = 2),
           compound_name=gsub("\\_[0-9]+","",compound_name)) %>%
    separate(Data.File,into=c("num","date","batch","sampleid","conc"),
             sep="__") %>%
    dplyr::select(num, date, batch, sampleid, conc, compound_name, heat_val) %>%
    ungroup() %>%
    add_count(date,batch,sampleid, compound_name, conc) %>%
    mutate(sampleid = ifelse(n > 1, paste(num, sampleid, sep = "__"), sampleid),
           sampleid = gsub(pattern = "X", replacement = "", sampleid)) %>%
    dplyr::select(sampleid, compound_name, heat_val) %>%
    pivot_wider(names_from = compound_name, values_from = heat_val, values_fill = NA) %>%
    filter(!str_detect(sampleid, "[Mm][Bb]"),
           !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
           !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
           !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
           !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
           !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
           !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
           !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
           !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]")) %>%
    drop_na(.) %>%
    column_to_rownames(., var = "sampleid") %>%
    as.matrix(.) %>%
    t(.) %>%
    pheatmap(., scale = "none",
             cellheight = nrow(heatmap_data2()) / (nrow(heatmap_data2())*0.075),
             cellwidth = ncol(heatmap_data2()) / (ncol(heatmap_data2())*0.0825)+1.5,
             angle_col = "90",
             color=colorRampPalette(c("navy", "white", "red"))(50),
    )
})


output$heatmap_plot2<- renderPlot(
  heatmap_plot2()
)

heatmap_data2 <- function()({
  indole_rawdf() %>%
    mutate(Data.File=as.character(Data.File)) %>%
    filter(#conc==checked,
      is.na(itsd),
      !grepl("[Cc][Cc][0-9]+", Data.File)) %>%
    left_join(indole_conc_int2()) %>%
    mutate(norm_peak = peakarea / avg) %>%
    group_by(compound_name) %>%
    mutate(compound_med = median(norm_peak)) %>%
    ungroup() %>%
    mutate(heat_val = log((norm_peak / compound_med), base = 2),
           compound_name=gsub("\\_[0-9]+","",compound_name)) %>%
    separate(Data.File,into=c("num","date","batch","sampleid","conc"),
             sep="__") %>%
    dplyr::select(num, date, batch, sampleid, conc, compound_name, heat_val) %>%
    ungroup() %>%
    add_count(date,batch,sampleid, compound_name, conc) %>%
    mutate(sampleid = ifelse(n > 1, paste(num, sampleid, sep = "__"), sampleid),
           sampleid = gsub(pattern = "X", replacement = "", sampleid)) %>%
    dplyr::select(sampleid, compound_name, heat_val) %>%
    pivot_wider(names_from = compound_name, values_from = heat_val, values_fill = NA) %>%
    filter(!str_detect(sampleid, "[Mm][Bb]"),
           !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
           !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
           !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
           !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
           !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
           !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
           !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
           !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]")) %>%
    drop_na(.) %>%
    column_to_rownames(., var = "sampleid") %>%
    as.matrix(.) %>%
    t(.)
})

output$heatmap_download2 <- downloadHandler(
  filename = function(){
    paste0("normalized_heatmap_",input$filename,"_",Sys.Date(),".pdf")
  },
  content = function(file) {
    pdf(file,
        height = nrow(heatmap_data2()) / (nrow(heatmap_data2())*0.075),
        width = (ncol(heatmap_data2()) / (ncol(heatmap_data2())*0.07))+1.5
    )
    heatmap_plot2()
    dev.off()
  },
  contentType = 'PDF'
)

#make heatmap:
heatmap_table <- function()({
  indole_rawdf() %>%
    mutate(Data.File=as.character(Data.File)) %>%
    filter(#conc==checked,
      is.na(itsd),
      !grepl("[Cc][Cc][0-9]+", Data.File)) %>%
    left_join(indole_conc_int2()) %>%
    group_by(compound_name) %>% 
    mutate(compound_med = median(peakarea)) %>% 
    ungroup() %>% 
    mutate(heat_val = log((norm_peak / compound_med), base = 2),
           compound_name=gsub("\\_[0-9]+","",compound_name)) %>%
    separate(Data.File,into=c("num","date","batch","sampleid","conc"),
             sep="__") %>%
    dplyr::select(num, date, batch, sampleid, conc, compound_name, norm_peak) %>%
    ungroup() %>%
    add_count(date,batch,sampleid, compound_name, conc) %>%
    mutate(sampleid = ifelse(n > 1, paste(num, sampleid, sep = "__"), sampleid),
           sampleid = gsub(pattern = "X", replacement = "", sampleid)) %>%
    dplyr::select(sampleid, compound_name, norm_peak) %>%
    pivot_wider(names_from = compound_name, values_from = norm_peak, values_fill = NA) %>%
    filter(!str_detect(sampleid, "[Mm][Bb]"),
           !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
           !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
           !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
           !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
           !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
           !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
           !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
           !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]")) %>%
    drop_na(.) %>%
    column_to_rownames(., var = "sampleid") %>%
    as.matrix(.) %>%
    t(.)
})


# output$TESTMAP <- renderDataTable(
#   heatmap_table()
# )

output$heatmap_table2 <- downloadHandler(
  
  filename = function(){
    paste0("heatmap_table",input$filename,"_",Sys.Date(),".csv")
  },
  
  content = function(file) {
    write.csv(heatmap_table(),file,
              row.names=T,quote=F)
  }
)

indole_subset_conc2 <- reactive({
  indole_rawdf()

  #is 33 rows when should be 66.. the conc averages aren't being included !
  #print(conc_int())

  indole_rawdf() %>%
    select(-conc) %>%
    left_join(indole_conc_int2()) %>%
    mutate(norm_peak=peakarea / avg)
})

output$conc_filter2 <- renderDataTable(
  indole_subset_conc2() %>%
    datatable() %>%
    formatRound(c("norm_peak"), 3) %>%
    formatStyle(columns = c("norm_peak"), 'text-align' = 'center')
)

#download table
#make wide table
normwide2 <- reactive({

  if(input$qcfil2==F){

    indole_rawdf() %>%
      mutate(Data.File=as.character(Data.File)) %>%
      filter(#conc==checked,
        is.na(itsd),
        !grepl("[Cc][Cc][0-9]+", Data.File)) %>%
      left_join(indole_conc_int2()) %>%
      mutate(norm_peak=peakarea / avg) %>%
      dplyr::select(Data.File, compound_name, norm_peak) %>%
      spread(compound_name, norm_peak, fill = NA) %>%
      #reshape2::dcast(Data.File ~ compound_name,value.var="norm_peak",fill=NA) %>%
      reshape2::melt(id.vars=c("Data.File")) %>%
      dplyr::rename(compound_name=variable,norm_peak=value) %>%
      separate(Data.File,into=c("num","date","batch","sampleid","conc"),
               sep="\\_\\_") %>%
      mutate(sampleid = ifelse(str_detect(sampleid, "[Mm][Bb]|[Pp][Oo][Oo][Ll][Ee][Dd]|[Bb][Hh][Ii][Qq][Cc]|[Pp][Ll][Aa][Ss][Mm][Aa]|[Hh][Ee][Xx][Aa][Nn][Ee][Ss]|[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]|50%_[Mm][Ee][Oo][Hh]|[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]|50%[Mm][Ee][Oo][hh]"),
                               paste(num, sampleid, conc, sep = "__"),
                               sampleid)) %>%
      mutate(sampleid = gsub(pattern = "conc.d", replacement = "conc", sampleid)) %>%
      dplyr::select(num,sampleid, compound_name, norm_peak) %>%
      mutate(norm_peak = round(norm_peak,5),
             compound_name=gsub("\\_[0-9]+","",compound_name)) %>%
      reshape2::dcast(num+sampleid ~ compound_name,value.var="norm_peak",fun.aggregate=mean) %>%
      arrange(num) %>%
      select(-num)

  }else{
    indole_rawdf() %>%
      mutate(Data.File=as.character(Data.File)) %>%
      filter(#conc==checked,
        is.na(itsd),
        !grepl("[Cc][Cc][0-9]+", Data.File)) %>%
      left_join(indole_conc_int2()) %>%
      mutate(norm_peak=peakarea / avg) %>%
      dplyr::select(Data.File, compound_name, norm_peak) %>%
      spread(compound_name, norm_peak, fill = NA) %>%
      reshape2::melt(id.vars=c("Data.File")) %>%
      dplyr::rename(compound_name=variable,norm_peak=value) %>%
      separate(Data.File,into=c("num","date","batch","sampleid","conc"),
               sep="\\_\\_") %>%
      dplyr::select(num,sampleid, compound_name, norm_peak) %>%
      mutate(norm_peak = round(norm_peak,5),
             compound_name=gsub("\\_[0-9]+","",compound_name)) %>%
      reshape2::dcast(num+sampleid ~ compound_name,value.var="norm_peak",fun.aggregate=mean) %>%
      arrange(num) %>%
      select(-num) %>%
      filter(!str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"))
  }
})

output$normwide2 <- renderDataTable(
  normwide2() %>%
    datatable() %>%
    formatRound(c(2:ncol(normwide2())), 3) %>%
    formatStyle(columns = c(2:ncol(normwide2())), 'text-align' = 'center')
)


normwide2_label <- reactive({
  if(input$qcfil2==F){
    return("normalized_results_")
  } else {
    return("removed_qcs_normalized_results_")
  }
})

#download handler
output$downloadData4 <- downloadHandler(

  filename = function(){
    paste0(normwide2_label(),input$filename,"_",Sys.Date(),".csv")
  },

  content = function(file) {
    write.csv(normwide2(),file,
              row.names=F,quote=F)
  }
)
  ##### QC Report #####
  # ITSD Raw Peak Area #

  # Build a Norm CV dataframe of summary stats
indole_rawdf2 <- reactive({
    compounds = unlist(strsplit(input$itsd_compounds, split=","))
    compounds = factor(compounds,level=unique(compounds))

    indole_rawdf() %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      inner_join(., indole_conc_int_sep2(), by = c("num","date","batch","sampleid","conc")) %>%
      filter(itsd == "ITSD",
             compound_name %in% compounds,
             !str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
             !grepl("[Cc][Cc][0-9]+", sampleid)) %>%
      mutate(peakarea = ifelse(peakarea <= input$zero_val2, input$zero_val2, peakarea)) %>%
      group_by(batch, compound_name) %>%
      summarise(stdev = sd(peakarea),
                average = mean(peakarea),
                middle = median(peakarea),
                cv = stdev / average,
                cv_med = stdev / median(peakarea)# Don't turn into % here since it will be applied in the y-axis scale
      )
  })

  indole_norm_plot1 <-function()({

    compounds = unlist(strsplit(input$itsd_compounds, split=","))
    compounds = factor(compounds,level=unique(compounds))

    indole_rawdf() %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      inner_join(., indole_conc_int_sep2(), by = c("num","date","batch","sampleid","conc")) %>%
      filter(itsd == "ITSD",
             compound_name %in% compounds,
             !str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]")) %>%
      mutate(peakarea = ifelse(peakarea <= input$zero_val2, input$zero_val2, peakarea),
             cc_shape = ifelse(grepl("cc[0-9]+", sampleid), "CC Sample", "ITSD")) %>%
      filter(itsd == "ITSD") %>%
      left_join(indole_rawdf2()) %>%
      mutate(num = as.numeric(num),
             flag = ifelse(peakarea > average + (1.5 * stdev) |
                             peakarea < average - (1.5 * stdev), paste(num,batch,sampleid, sep = "_"), NA)) %>%
      ggplot(., aes(x = num, y = peakarea, group = compound_name, label = flag)) +
      geom_point(aes(color = compound_name, shape = cc_shape), fill = "black", alpha = 0.6, size = 2) +
      geom_line(aes(y = average, color = compound_name)) +
      geom_line(aes(y = average + 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_line(aes(y = average - 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_ribbon(aes(ymin = average - stdev, ymax = average + stdev,
                      fill = compound_name), alpha=0.2) +
      ggrepel::geom_label_repel(size = 1.2, max.overlaps = Inf,
                                min.segment.length = 0.1, label.padding = 0.1) +
      theme_bw() +
      theme(panel.grid.minor= element_blank(),
            panel.grid.major.x = element_blank(),
            legend.text = element_text(color = "black", size = 13),
            legend.title = element_text(color = "black", size = 15),
            legend.position = "top",
            strip.text=element_text(color = "black", size=9),
            axis.text =element_text(color = "black", size=8),
            axis.title = element_text(color = "black", size = 15),
            plot.margin = margin(1,1,0,1.1, unit = 'cm')) +
      ggsci::scale_color_futurama() +
      ggsci::scale_fill_futurama() +
      guides(color = guide_legend(title = "Internal Standard Compound",
                                  override.aes = list(size = 2.5), nrow = 1,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right"), fill = F,
             shape = guide_legend(title = "",
                                  override.aes = list(size = 2.5), nrow = 2,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right")) +
      scale_shape_manual(values = c(24,16))+
      scale_y_continuous(label = scales::scientific) +
      scale_x_continuous(breaks = seq(0,150,25)) +
      ylab("Raw Peak Area\n") +
      xlab("\nInjection Number")+
      facet_grid(~compound_name, scales="free_x")
  })


indole_norm_plot2 <-function()({

    compounds = unlist(strsplit(input$itsd_compounds, split=","))
    compounds = factor(compounds,level=unique(compounds))

    indole_rawdf() %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      inner_join(., indole_conc_int_sep2(), by = c("num","date","batch","sampleid","conc")) %>%
      filter(itsd == "ITSD",
             compound_name %in% compounds,
             !str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]")) %>%
      mutate(peakarea = ifelse(peakarea <= input$zero_val2, input$zero_val2, peakarea),
             cc_shape = ifelse(grepl("cc[0-9]+", sampleid), "CC Sample", "ITSD")) %>%
      filter(itsd == "ITSD") %>%
      left_join(indole_rawdf2()) %>%
      mutate(num = as.numeric(num),
             flag = ifelse(peakarea > average + (1.5 * stdev) |
                             peakarea < average - (1.5 * stdev), paste(num,batch,sampleid, sep = "_"), NA)) %>%
      ggplot(., aes(x = num, y = peakarea, group = compound_name, label = flag)) +
      geom_point(aes(color = compound_name, shape = cc_shape), fill = "black", alpha = 0.6, size = 2) +
      geom_line(aes(y = average, color = compound_name)) +
      geom_line(aes(y = average + 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_line(aes(y = average - 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_ribbon(aes(ymin = ifelse(average - stdev <= 1, 1, average - stdev),
                      ymax = average + stdev,
                      fill = compound_name), alpha=0.2) +
      ggrepel::geom_label_repel(size = 1.2, max.overlaps = Inf,
                                min.segment.length = 0.1, label.padding = 0.1) +
      theme_bw() +
      theme(panel.grid.minor= element_blank(),
            panel.grid.major.x = element_blank(),
            legend.text = element_text(color = "black", size = 13),
            legend.title = element_text(color = "black", size = 15),
            legend.position = "top",
            strip.text=element_text(color = "black", size=9),
            axis.text =element_text(color = "black", size=8),
            axis.title = element_text(color = "black", size = 15),
            plot.margin = margin(1,1,0,1.1, unit = 'cm')) +
      ggsci::scale_color_futurama() +
      ggsci::scale_fill_futurama() +
      guides(color = guide_legend(title = "Internal Standard Compound",
                                  override.aes = list(size = 2.5), nrow = 1,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right"), fill = F,
             shape = guide_legend(title = "",
                                  override.aes = list(size = 2.5), nrow = 2,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right")) +
      scale_shape_manual(values = c(24,16))+
      ylab("Raw Peak Area\n(log10 scale)") +
      xlab("\nInjection Number") +
      # coord_trans(y = "log10") +
      scale_y_continuous(trans = "log10", labels = scales::scientific)+
      scale_x_continuous(breaks = seq(0,150,25))+
      facet_grid(~compound_name, scales="free_x")
  })

  # ITSD CV Percent #

  # Build CV dataframe
indole_rawdf2_1 <- reactive({

    compounds = unlist(strsplit(input$itsd_compounds, split=","))
    compounds = factor(compounds,level=unique(compounds))

    indole_rawdf() %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      inner_join(., indole_conc_int_sep2(), by = c("num","date","batch","sampleid","conc")) %>%
      filter(itsd == "ITSD",
             compound_name %in% compounds,
             !str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
             !grepl("[Cc][Cc][0-9]+", sampleid)) %>%
      mutate(peakarea = ifelse(peakarea <= input$zero_val2, input$zero_val2, peakarea)) %>%
      group_by(batch, compound_name) %>%
      summarise(stdev = sd(peakarea),
                average = mean(peakarea),
                middle = median(peakarea),
                cv = stdev / average,
                cv_med = stdev / median(peakarea)# Don't turn into % here since it will be applied in the y-axis scale
      ) %>%
      dplyr::rename(Batch = batch,
                    `Internal Standard` = compound_name,
                    StDev = stdev,
                    Mean = average,
                    Median = middle,
                    CV = cv,
                    `CV Median` = cv_med) %>%
      mutate(StDev = round(StDev, digits = 0),
             Mean = round(Mean, digits = 0),
             Median = round(Median, digits = 0),
             `CV (%)` = round(CV * 100, digits = 1),
             `CV Median (%)` = round(`CV Median` * 100, digits = 1)) %>%
      select(Batch,`Internal Standard`,StDev,Mean,Median,`CV (%)`,`CV Median (%)`)
  })



  # Build summary table
  indole_norm_table_list1 <- function()({
    gridExtra::tableGrob(indole_rawdf2_1(), rows = NULL, theme = tt)
  })

  ## Save Report ##
  output$indole_norm_qc_report_download <- downloadHandler(
    filename = function(){
      paste0("Indole_QC_Norm_Report_",unique(indole_rawdf2()$batch),"_",Sys.Date(),".pdf")
    },
    content = function(file) {
      pdf(file, height = 11, width = 8.5)
      print(
        gridExtra::grid.arrange(
          indole_norm_plot1() +
            xlab("") +
            ggtitle(paste("Indole Qualitative QC Report\n", unique(indole_rawdf2()$batch))) +
            theme(plot.title = element_text(color = "black",
                                            hjust = 0.5,
                                            size = 20,
                                            face = "bold")),
          indole_norm_plot2() + theme(legend.position = "none"),
          nrow = 2)
      )
      print(
        gridExtra::grid.arrange(
          indole_norm_table_list1())
      )
      dev.off()
    }
  )

# 
  
  # Bile acid quant tab --------------------------------------------------------
  #make concentration table
  conc_tbl5 <- reactive({
    get_all_conc(input$start5,
                 compounds=tolower(input$compounds5),series=input$series5,fold=3)
  })

  output$conc_tbl5 <- renderDataTable(
    conc_tbl5(),
    options=list(pageLength=5)
  )

  quant_conc_tbl5 <- reactive({

    compounds = tolower(unlist(strsplit(input$compounds5, split=",")))
    concs <- unlist(strsplit(input$quant_conc5, split=","))

    int <- cbind(compounds,concs) %>% as.data.frame()
    colnames(int) <- c("compound_name","conc")

    int %>%
      mutate(conc=ifelse(conc=="dil","diluted","concentrated"))

  })


  output$conc5 <- renderDataTable(
    quant_conc_tbl5(),
    options = list(pageLength=5)
  )

  #make cutoff reactive
  cutoff_df5 <- reactive({
    compounds = tolower(unlist(strsplit(input$compounds5, split=",")))
    max_cutoffs = unlist(strsplit(input$maxcc5, split = ","))
    min_cutoffs = unlist(strsplit(input$mincc5, split = ","))

    cutoff_df <- cbind(compounds,max_cutoffs,min_cutoffs)
    colnames(cutoff_df) <- c("compound_name","maxcc","mincc")

    cutoff_df %>%
      as.data.frame() %>%
      mutate(maxcc=as.character(maxcc),
             mincc=as.character(mincc),
             maxcc=as.numeric(maxcc),
             mincc=as.numeric(mincc))
  })

  output$cutoff_df5 <- renderDataTable(
    cutoff_df5(),
    options = list(pageLength=5)
  )

  # read in table as reactive
  meta5 <- reactive({ readin_bile_csv_single_file(file.path(wddir,input$filename), na.value = input$quant_zero_val_bile_acid) })

  #show models and make plots
  modelstart5 <- reactive({

    compounds = tolower(unlist(strsplit(input$compounds5, split=",")))
    compounds = factor(compounds,level=unique(compounds))

    model <- meta5() %>%
      filter(compound_name %in% compounds) %>%
      inner_join(quant_conc_tbl5()) %>%
      mutate(compound_name=factor(compound_name,levels=compounds)) %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(sampleid+compound_name+conc+letter ~ itsd,value.var="peakarea",
                      fun.aggregate=mean) %>%
      group_by(sampleid, letter) %>%
      mutate(ITSD=zoo::na.locf(ITSD),
             norm_peak = peak / ITSD,
             curveLab=str_extract(sampleid,pattern="CC[1-9][0-9]+|CC[1-9]+")) %>%
      inner_join(conc_tbl5()) %>%
      left_join(cutoff_df5()) %>%
      filter(conc_val <= maxcc,
             conc_val >= mincc) %>%
      group_by(compound_name) %>%
      summarize(r = cor(norm_peak,conc_val),
                model_list <- broom::tidy(lm(norm_peak ~ conc_val))) %>%
      reshape2::dcast(compound_name+r ~ term,value.var="estimate") %>%
      dplyr::rename(slope_value=conc_val)
  })

  output$model5 <- renderDataTable(
    modelstart5() %>%
      datatable() #%>%
    # formatRound(c(2:4), 3) %>%
    # formatStyle(columns = c(2:4), 'text-align' = 'center')
  )

  #make quant graph
  quant_plot5 <- reactive({

    compounds5 = tolower(unlist(strsplit(input$compounds5, split=",")))
    compounds5 = factor(compounds5,level=unique(compounds5))

    meta5() %>%
      filter(compound_name %in% compounds5) %>%
      inner_join(quant_conc_tbl5()) %>%
      mutate(compound_name=factor(compound_name,levels=compounds5)) %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(sampleid+compound_name+conc+letter ~ itsd,value.var="peakarea",
                      fun.aggregate=mean) %>%
      group_by(sampleid, letter) %>%
      mutate(ITSD=zoo::na.locf(ITSD),
             norm_peak = peak / ITSD,
             curveLab=str_extract(sampleid,pattern="CC[1-9][0-9]+|CC[1-9]+")) %>%
      inner_join(conc_tbl5()) %>%
      left_join(cutoff_df5()) %>%
      filter(conc_val <= maxcc,
             conc_val >= mincc) %>%
      ggplot(aes(x=conc_val,y=norm_peak)) +
      geom_point(size=3) +
      geom_smooth(method="lm") +
      theme(strip.text=element_text(size=15),
            axis.text=element_text(size=13)) +
      facet_wrap("compound_name",scales="free")
  })

  output$quant5 <- renderPlot(
    quant_plot5()
  )

  #make quant table
  quant_table5 <- function()({

    compounds5 = tolower(unlist(strsplit(input$compounds5, split=",")))

    meta5() %>%
      filter(compound_name %in% compounds5) %>%
      inner_join(quant_conc_tbl5()) %>%
      mutate(compound_name=factor(compound_name,levels=compounds5)) %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(sampleid+compound_name+conc+letter ~ itsd,value.var="peakarea",
                      fun.aggregate=mean) %>%
      group_by(sampleid, letter) %>%
      mutate(ITSD=zoo::na.locf(ITSD),
             norm_peak = peak / ITSD) %>%
      filter(!grepl("^CC[0-9]",sampleid)) %>%
      left_join(modelstart5()) %>%
      mutate(quant_val =  (norm_peak - (`(Intercept)`))/slope_value*as.numeric(input$xfactor5)) %>%
      arrange(compound_name) %>%
      mutate(quant_val = ifelse(quant_val < 0,0,quant_val),
             quant_val = round(quant_val,2))
  })

  output$quant_tbl5 <- renderDataTable(
    quant_table5() %>%
      datatable() #%>%
    # formatRound(c(4:ncol(quant_table5())), 3) %>%
    # formatStyle(columns = c(4:ncol(quant_table5())), 'text-align' = 'center')
  )

  quant_table_dl5 <- function()({

    compounds5 = tolower(unlist(strsplit(input$compounds5, split=",")))

    meta5() %>%
      filter(compound_name %in% compounds5) %>%
      inner_join(quant_conc_tbl5()) %>%
      mutate(compound_name=factor(compound_name,levels=compounds5)) %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(sampleid+compound_name+conc+letter ~ itsd,value.var="peakarea",
                      fun.aggregate=mean) %>%
      group_by(sampleid, letter) %>%
      fill(ITSD) %>%
      mutate(ITSD=zoo::na.locf(ITSD),
             norm_peak = peak / ITSD) %>%
      filter(!grepl("^CC[0-9]",sampleid)) %>%
      left_join(modelstart5()) %>%
      mutate(quant_val =  (norm_peak - (`(Intercept)`))/slope_value*as.numeric(input$xfactor5)) %>%
      arrange(compound_name) %>%
      mutate(quant_val = ifelse(quant_val < 0,0,quant_val),
             quant_val = round(quant_val,2)) %>%
      reshape2::dcast(sampleid ~ compound_name, value.var="quant_val",fill=0) #%>%
    #separate(Data.File,into=c("num","date_run","batch","sampleid","conc"),sep="\\_\\_") %>%
    #arrange(num) %>%
    #select(-num,-date_run,-batch,-conc)

  })
  
  ## Save quant table with and without QCs
  observeEvent( 
    input$downloadData5,
    {
      write.csv(quant_table_dl5(), 
                paste0("/Volumes/chaubard-lab/shiny_workspace/CLIN_Finals_QCs/quant_results_",
                       gsub("\\.csv","",input$filename),"_",
                       gsub("\\-","",Sys.Date()),".csv"), row.names=F,quote=F)
      write.csv(quant_table_dl5() %>% filter(!str_detect(sampleid, "[Mm][Bb]"),
                                             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
                                             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
                                             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
                                             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
                                             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
                                             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
                                             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
                                             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
                                            !grepl("CC[0-9]+", sampleid)), 
                paste0("/Volumes/chaubard-lab/shiny_workspace/CLIN_Finals/removed_qcs_quant_results_",
                       gsub("\\.csv","",input$filename), "_",
                       gsub("\\-","",Sys.Date()),".csv"), row.names=F,quote=F)
    } 
  )

  # Make quant_bar dataframe to assign sizes based on number of samples
  bile_quant_bar_data <- function()({

    compounds5 = tolower(unlist(strsplit(input$compounds5, split=",")))
    compounds5 = factor(compounds5,level=unique(compounds5))

    if(input$bile_qcfil_quant==F){

      meta5() %>%
        filter(compound_name %in% compounds5) %>%
        inner_join(quant_conc_tbl5()) %>%
        mutate(compound_name=factor(compound_name,levels=compounds5)) %>%
        replace_na(list(itsd="peak")) %>%
        reshape2::dcast(sampleid+compound_name+conc+letter ~ itsd,value.var="peakarea",
                        fun.aggregate=mean) %>%
        group_by(sampleid, letter) %>%
        mutate(ITSD=zoo::na.locf(ITSD),
               norm_peak = peak / ITSD) %>%
        filter(!grepl("^CC[0-9]",sampleid)) %>%
        left_join(modelstart5()) %>%
        mutate(quant_val =  (norm_peak - (`(Intercept)`))/slope_value*as.numeric(input$xfactor5)) %>%
        arrange(compound_name) %>%
        mutate(quant_val = ifelse(quant_val < 0,0,quant_val),
               quant_val = round(quant_val,2)) %>%
        select(sampleid, compound_name, quant_val)
    } else {
      meta5() %>%
        filter(compound_name %in% compounds5) %>%
        inner_join(quant_conc_tbl5()) %>%
        mutate(compound_name=factor(compound_name,levels=compounds5)) %>%
        replace_na(list(itsd="peak")) %>%
        reshape2::dcast(sampleid+compound_name+conc+letter ~ itsd,value.var="peakarea",
                        fun.aggregate=mean) %>%
        group_by(sampleid, letter) %>%
        mutate(ITSD=zoo::na.locf(ITSD),
               norm_peak = peak / ITSD) %>%
        filter(!grepl("^CC[0-9]",sampleid)) %>%
        left_join(modelstart5()) %>%
        mutate(quant_val =  (norm_peak - (`(Intercept)`))/slope_value*as.numeric(input$xfactor5)) %>%
        arrange(compound_name) %>%
        mutate(quant_val = ifelse(quant_val < 0,0,quant_val),
               quant_val = round(quant_val,2)) %>%
        select(sampleid, compound_name, quant_val) %>%
        filter(!str_detect(sampleid, "[Mm][Bb]"),
               !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
               !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
               !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
               !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
               !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
               !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
               !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
               !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
               !grepl("CC[0-9]+", sampleid))
    }
  })

  #make quant bar graph
  bile_quant_barplot <- function()({

    compounds5 = tolower(unlist(strsplit(input$compounds5, split=",")))
    compounds5 = factor(compounds5,level=unique(compounds5))

    if(input$bile_qcfil_quant==F){

      meta5() %>%
        filter(compound_name %in% compounds5) %>%
        inner_join(quant_conc_tbl5()) %>%
        mutate(compound_name=factor(compound_name,levels=compounds5)) %>%
        replace_na(list(itsd="peak")) %>%
        reshape2::dcast(sampleid+compound_name+conc+letter ~ itsd,value.var="peakarea",
                        fun.aggregate=mean) %>%
        group_by(sampleid, letter) %>%
        mutate(ITSD=zoo::na.locf(ITSD),
               norm_peak = peak / ITSD) %>%
        filter(!grepl("^CC[0-9]",sampleid)) %>%
        left_join(modelstart5()) %>%
        mutate(quant_val =  (norm_peak - (`(Intercept)`))/slope_value*as.numeric(input$xfactor5)) %>%
        arrange(compound_name) %>%
        mutate(quant_val = ifelse(quant_val < 0,0,quant_val),
               quant_val = round(quant_val,2)) %>%
        select(sampleid, compound_name, quant_val) %>%
        ggplot(aes(x = sampleid, y = quant_val, fill = compound_name)) +
        geom_bar(stat="identity") +
        theme_bw() +
        theme(plot.margin=unit(c(5.5, 15, 5.5, 10),"points"),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              legend.text = element_text(size = 14),
              legend.title = element_text(size = 16),
              legend.position = "top",
              strip.text=element_text(size=16),
              axis.text.y =element_text(size=11),
              axis.title = element_text(size = 15),
              axis.text.x =element_text(vjust = 0.5,
                                        hjust = 1,
                                        angle = 90,
                                        size = ifelse(
                                          (((nrow(bile_quant_bar_data())))+(nrow(bile_quant_bar_data())))*1000 /
                                            (nrow(bile_quant_bar_data()) * nrow(bile_quant_bar_data())) >= 11,
                                          11,
                                          (((nrow(bile_quant_bar_data())))+(nrow(bile_quant_bar_data())))*1000 /
                                            (nrow(bile_quant_bar_data()) * nrow(bile_quant_bar_data())))
              )) +
        facet_grid(~compound_name, space = "free_x") +
        scale_fill_manual(values = c(pal_ucscgb("default", alpha = 0.7)(7), "bisque4"))+
        # scale_fill_locuszoom() +
        xlab("\nSampleID") +
        ylab("Concentration (ug/mL)\n") +
        guides(fill = guide_legend(title="Compound    "))

    } else {

      meta5() %>%
        filter(compound_name %in% compounds5) %>%
        inner_join(quant_conc_tbl5()) %>%
        mutate(compound_name=factor(compound_name,levels=compounds5)) %>%
        replace_na(list(itsd="peak")) %>%
        reshape2::dcast(sampleid+compound_name+conc+letter ~ itsd,value.var="peakarea",
                        fun.aggregate=mean) %>%
        group_by(sampleid, letter) %>%
        mutate(ITSD=zoo::na.locf(ITSD),
               norm_peak = peak / ITSD) %>%
        filter(!grepl("^CC[0-9]",sampleid)) %>%
        left_join(modelstart5()) %>%
        mutate(quant_val =  (norm_peak - (`(Intercept)`))/slope_value*as.numeric(input$xfactor5)) %>%
        arrange(compound_name) %>%
        mutate(quant_val = ifelse(quant_val < 0,0,quant_val),
               quant_val = round(quant_val,2)) %>%
        select(sampleid, compound_name, quant_val) %>%
        filter(!str_detect(sampleid, "[Mm][Bb]"),
               !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
               !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
               !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
               !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
               !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
               !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
               !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
               !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
               !grepl("CC[0-9]+", sampleid)
        ) %>%
        ggplot(aes(x = sampleid, y = quant_val, fill = compound_name)) +
        geom_bar(stat="identity") +
        theme_bw() +
        theme(plot.margin=unit(c(5.5, 15, 5.5, 10),"points"),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              legend.position = "top",
              legend.text = element_text(size = 14),
              legend.title = element_text(size = 16),
              strip.text=element_text(size=18),
              axis.text.y =element_text(size=11),
              axis.title = element_text(size = 15),
              axis.text.x =element_text(vjust = 0.5,
                                        hjust = 1,
                                        angle = 90,
                                        size = ifelse(
                                          (((nrow(bile_quant_bar_data())))+(nrow(bile_quant_bar_data())))*1000 /
                                            (nrow(bile_quant_bar_data()) * nrow(bile_quant_bar_data())) >= 11,
                                          11,
                                          (((nrow(bile_quant_bar_data())))+(nrow(bile_quant_bar_data())))*1000 /
                                            (nrow(bile_quant_bar_data()) * nrow(bile_quant_bar_data())))
              )) +
        facet_grid(~compound_name, space = "free_x") +
        scale_fill_manual(values = c(pal_ucscgb("default", alpha = 0.7)(7), "bisque4"))+
        # scale_fill_locuszoom() +
        xlab("\nSampleID") +
        ylab("Concentration (ug/mL)\n") +
        guides(fill = guide_legend(title="Compound    "))
    }
  })

  # Output the bar graph
  output$bile_quant_barplot <- renderPlot(
    bile_quant_barplot()
  )


  # Produce file labels based on bile_qcfil_quant
  bile_quant_bar_data_label <- reactive({
    if(input$bile_qcfil_quant==F){
      return("quant_barplot_")
    } else {
      return("removed_qcs_quant_barplot_")
    }
  })

  output$bile_quant_download <- downloadHandler(
    filename = function(){
      paste0(bile_quant_bar_data_label(),input$filename,"_",Sys.Date(),".pdf")
    },
    content = function(file) {
      ggsave(file,plot=bile_quant_barplot(),
             width = nrow(bile_quant_bar_data())*33 / (nrow(bile_quant_bar_data()))+1.5,
             height = ncol(bile_quant_bar_data())*7.5 / (ncol(bile_quant_bar_data())))
    }
  )


  #### QC Report ####
  meta5_1 <- reactive({

    compounds = tolower(unlist(strsplit(input$compounds5, split=",")))

    meta5() %>%
      filter(compound_name %in% compounds) %>%
      inner_join(., quant_conc_tbl5(), by = c("compound_name", "conc")) %>%
      separate(Data.File,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      filter(
        compound_name == "isodeoxycholic acid" |
               compound_name == "alloisolithocholic acid" |
               compound_name == "3-oxolithocholic acid" |
               itsd == "ITSD",
        !str_detect(sampleid, "[Mm][Bb]"),
        !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
        !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
        !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
        !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
        !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
        !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
        !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
        !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]")) %>%
      mutate(peakarea = ifelse(peakarea <= input$zero_val_bile_acid, input$zero_val_bile_acid, peakarea),
             num = gsub("X","",num),
             num = as.numeric(num),
             cc_shape = ifelse(grepl("CC[0-9]+", sampleid), "Calibration Curve Sample", "Standard Sample"))
  })

  # Build another dataframe of summary stats
  meta5_2 <- reactive({

    compounds = tolower(unlist(strsplit(input$compounds5, split=",")))

    meta5() %>%
      filter(compound_name %in% compounds) %>%
      inner_join(., quant_conc_tbl5(), by = c("compound_name", "conc")) %>%
      separate(Data.File,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      filter(
        compound_name == "isodeoxycholic acid" |
               compound_name == "alloisolithocholic acid" |
               compound_name == "3-oxolithocholic acid" |
               itsd == "ITSD",
        !str_detect(sampleid, "[Mm][Bb]"),
        !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
        !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
        !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
        !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
        !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
        !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
        !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
        !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
        !grepl("[Cc][Cc][0-9]+", sampleid)) %>%
      mutate(peakarea = ifelse(peakarea <= input$zero_val_bile_acid, input$zero_val_bile_acid, peakarea),
             num = gsub("X","",num),
             num = as.numeric(num)) %>%
      group_by(batch, compound_name) %>%
      summarise(stdev = sd(peakarea),
                average = mean(peakarea),
                middle = median(peakarea),
                cv = stdev / average,
                cv_med = stdev / median(peakarea))
  })

  ba_qc_plot1 <-function()({

    meta5_1() %>%
      left_join(as.data.frame(meta5_2())) %>%
      mutate(flag = ifelse(peakarea > average + (1.5 * stdev) |
                             peakarea < average - (1.5 * stdev), paste(num,batch,sampleid, sep = "_"), NA)) %>%
      group_by(compound_name) %>%
      ggplot(., aes(x = num, y = peakarea, group = compound_name, label = flag)) +
      geom_point(aes(color = compound_name, shape = cc_shape), fill = "black", alpha = 0.6, size = 2) +
      geom_line(aes(y = average, color = compound_name)) +
      geom_line(aes(y = average + 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_line(aes(y = average - 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_ribbon(aes(ymin = average - stdev, ymax = average + stdev,
                      fill = compound_name), alpha=0.2) +
      ggrepel::geom_label_repel(size = 1.2, max.overlaps = Inf,
                                min.segment.length = 0.1, label.padding = 0.1) +
      theme_bw() +
      theme(panel.grid.minor= element_blank(),
            panel.grid.major.x = element_blank(),
            legend.text = element_text(color = "black", size = 8),
            legend.title = element_text(color = "black", size = 10),
            legend.position = "top",
            strip.text=element_text(color = "black", size=5),
            axis.text =element_text(color = "black", size=8),
            axis.title = element_text(color = "black", size = 14),
            plot.margin = margin(1,0.5,0,0.6, unit = 'cm')) +
      ggsci::scale_color_ucscgb() +
      ggsci::scale_fill_ucscgb() +
      guides(color = guide_legend(title = "Internal Standard Compound",
                                  override.aes = list(size = 2.5), nrow = 2,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right"), fill = F,
             shape = guide_legend(title = "",
                                  override.aes = list(size = 2.5), nrow = 2,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right")) +
      scale_shape_manual(values = c(24,16))+
      scale_y_continuous(label = scales::scientific) +
      scale_x_continuous(breaks = seq(0,150,25)) +
      ylab("Raw Peak Area\n") +
      xlab("\nInjection Number")+
      facet_grid(~compound_name, scales="free_x")
  })

  ba_qc_plot2 <- function()({
    meta5_1() %>%
      left_join(as.data.frame(meta5_2())) %>%
      mutate(flag = ifelse(peakarea > average + (1.5 * stdev) |
                             peakarea < average - (1.5 * stdev), paste(num,batch,sampleid, sep = "_"), NA)) %>%
      group_by(compound_name) %>%
      ggplot(., aes(x = num, y = peakarea, group = compound_name, label = flag)) +
      geom_point(aes(color = compound_name, shape = cc_shape), fill = "black", alpha = 0.6, size = 2) +
      geom_line(aes(y = average, color = compound_name)) +
      geom_line(aes(y = average + 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_line(aes(y = average - 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_ribbon(aes(ymin = ifelse(average - stdev <= 1, 1, average - stdev),
                      ymax = average + stdev,
                      fill = compound_name), alpha=0.2) +
      ggrepel::geom_label_repel(size = 1.2, max.overlaps = Inf,
                                min.segment.length = 0.1, label.padding = 0.1) +
      theme_bw() +
      theme(panel.grid.minor= element_blank(),
            panel.grid.major.x = element_blank(),
            legend.text = element_text(color = "black", size = 8),
            legend.title = element_text(color = "black", size = 10),
            legend.position = "top",
            strip.text=element_text(color = "black", size=5),
            axis.text =element_text(color = "black", size=8),
            axis.title = element_text(color = "black", size = 14),
            plot.margin = margin(3.5,0.5,0,0.6, unit = 'cm')) +
      ggsci::scale_color_ucscgb() +
      ggsci::scale_fill_ucscgb() +
      guides(color = guide_legend(title = "Internal Standard Compound",
                                  override.aes = list(size = 2.5), nrow = 2,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right"), fill = F,
             shape = guide_legend(title = "",
                                  override.aes = list(size = 2.5), nrow = 2,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right")) +
      scale_shape_manual(values = c(24,16))+
      ylab("Raw Peak Area\n(log10 scale)") +
      xlab("\nInjection Number") +
      # coord_trans(y = "log10") +
      scale_y_continuous(trans = "log10", labels = scales::scientific)+
      scale_x_continuous(breaks = seq(0,150,25))+
      facet_grid(~compound_name, scales="free_x")
  })

  # Show models and make plots
  ba_model <- function()({

    compounds = tolower(unlist(strsplit(input$compounds5, split=",")))

    meta5() %>%
      filter(compound_name %in% compounds) %>%
      inner_join(quant_conc_tbl5()) %>%
      mutate(compound_name=factor(compound_name,levels=compounds)) %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(batch+sampleid+compound_name+conc+letter ~ itsd,value.var="peakarea",
                      fun.aggregate=mean) %>%
      group_by(sampleid, letter) %>%
      mutate(ITSD=zoo::na.locf(ITSD),
             norm_peak = peak / ITSD,
             curveLab=str_extract(sampleid,pattern="CC[1-9][0-9]+|CC[1-9]+")) %>%
      inner_join(conc_tbl5()) %>%
      left_join(cutoff_df5()) %>%
      filter(conc_val <= maxcc,
             conc_val >= mincc) %>%
      group_by(batch,compound_name,conc) %>%
      summarize(r = cor(norm_peak,conc_val),
                model_list <- broom::tidy(lm(norm_peak ~ conc_val))) %>%
      reshape2::dcast(batch+compound_name+conc+r ~ term,value.var="estimate") %>%
      dplyr::rename(Batch = batch,
                    Concentrated = conc,
                    `Internal Standard` = compound_name,
                    Intercept = `(Intercept)`,
                    Slope = conc_val)
  })

  # Plot CC models
  ba_qc_plot3 <- function()({

    compounds = tolower(unlist(strsplit(input$compounds5, split=",")))

    meta5() %>%
      filter(compound_name %in% compounds,
             !str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]")) %>%
      inner_join(quant_conc_tbl5()) %>%
      mutate(compound_name=factor(compound_name,levels=compounds)) %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(batch+sampleid+compound_name+conc+letter ~ itsd,value.var="peakarea",
                      fun.aggregate=mean) %>%
      group_by(sampleid, letter) %>%
      mutate(ITSD=zoo::na.locf(ITSD),
             norm_peak = peak / ITSD,
             curveLab=str_extract(sampleid,pattern="CC[1-9][0-9]+|CC[1-9]+"),
             conc = ifelse(conc == "conc","Concentrated Standard", "Diluted Standard")) %>%
      inner_join(conc_tbl5()) %>%
      left_join(cutoff_df5()) %>%
      filter(conc_val <= maxcc,
             conc_val >= mincc) %>%
      ungroup() %>%
      group_by(batch,compound_name,conc) %>%
      ggplot(aes(x=conc_val,y=norm_peak)) +
      geom_point(size=3) +
      geom_smooth(method="lm", formula = y~x) +
      ggpmisc::stat_poly_eq(formula = y~x,
                            rr.digits = 5,
                            coef.digits = 5,
                            size = 3.2,
                            aes(label = paste(..eq.label..,
                                              ..rr.label..,
                                              sep = "~~~~")),
                            parse = TRUE) +
      theme(strip.text=element_text(size=12, color = "black"),
            axis.text=element_text(size=9, color = "black"),
            axis.title = element_text(size = 13, color = "black"),
            plot.title = element_text(hjust = 0.5, size = 16,
                                      color = "black", face = "bold"),
            plot.margin = margin(1,0.5,0.5,0.5, unit = 'cm')) +
      facet_wrap(~compound_name+conc,scales="free", ncol = 2) +
      xlab("\nConcentration (ug/mL)") +
      ylab("Normalized Peak Area\n") +
      ggtitle("Calibration Curves")
  })

  ## ITSD CV Percent ##

  # Build CV dataframe
  meta5_3 <- reactive({

    compounds = tolower(unlist(strsplit(input$compounds5, split=",")))

    meta5() %>%
      filter(compound_name %in% compounds) %>%
      inner_join(., quant_conc_tbl5(), by = c("compound_name", "conc")) %>%
      separate(Data.File,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      filter(
        compound_name == "isodeoxycholic acid" |
               compound_name == "alloisolithocholic acid" |
               compound_name == "3-oxolithocholic acid" |
               itsd == "ITSD",
        !str_detect(sampleid, "[Mm][Bb]"),
        !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
        !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
        !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
        !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
        !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
        !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
        !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
        !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
        !grepl("[Cc][Cc][0-9]+", sampleid)) %>%
      mutate(peakarea = ifelse(peakarea <= input$zero_val_bile_acid, input$zero_val_bile_acid, peakarea),
             num = gsub("X","",num),
             num = as.numeric(num),
             cc_shape = ifelse(grepl("CC[0-9]+", sampleid), "CC Sample", "ITSD")) %>%
      group_by(batch, compound_name) %>%
      summarise(stdev = sd(peakarea),
                average = mean(peakarea),
                middle = median(peakarea),
                cv = stdev / average,
                cv_med = stdev / median(peakarea)
      ) %>%
      dplyr::rename(Batch = batch,
                    `Internal Standard` = compound_name,
                    StDev = stdev,
                    Mean = average,
                    Median = middle,
                    CV = cv,
                    `CV Median` = cv_med) %>%
      mutate(StDev = round(StDev, digits = 0),
             Mean = round(Mean, digits = 0),
             Median = round(Median, digits = 0),
             `CV (%)` = round(CV * 100, digits = 1),
             `CV Median (%)` = round(`CV Median` * 100, digits = 1)) %>%
      select(Batch,`Internal Standard`,StDev,Mean,Median,`CV (%)`,`CV Median (%)`)
  })

  # Build summary table
  ba_tt <- gridExtra::ttheme_default(
    colhead=list(fg_params=list(col="black", fontface=2L)),
    padding = unit(c(0.5,0.75), "cm"))

  ba_table_list1 <- function()({
    gridExtra::tableGrob(meta5_3(), rows = NULL, theme = ba_tt)
  })
  
  ### Plasma QC ###
  
  # Build summary specifically for high and low ranges for horizontal lines
  meta5_4 <- reactive({
    
    compounds5 = tolower(unlist(strsplit(input$compounds5, split=",")))
    
    meta5() %>%
      filter(compound_name %in% compounds5) %>%
      inner_join(quant_conc_tbl5()) %>%
      mutate(compound_name=factor(compound_name,levels=compounds5)) %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(sampleid+compound_name+conc+letter ~ itsd,value.var="peakarea",
                      fun.aggregate=mean) %>%
      group_by(sampleid, letter) %>%
      fill(ITSD) %>%
      mutate(ITSD=zoo::na.locf(ITSD),
             norm_peak = peak / ITSD) %>%
      filter(!grepl("^CC[0-9]",sampleid),
             grepl(pattern = "[Pp][Ll][Aa][Ss][Mm][Aa]", sampleid)) %>%
      left_join(modelstart5()) %>%
      mutate(quant_val =  (norm_peak - (`(Intercept)`))/slope_value*as.numeric(input$xfactor5)) %>%
      mutate(sample = str_extract(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             qc = str_extract(sampleid, "[Qq][Cc]"),
             replicate = as.factor(str_extract(sampleid, "[0-9]+")),
             batch = unique(meta5_1()$batch)) %>%
      arrange(compound_name) %>%
      mutate(quant_val = ifelse(quant_val < 0,0,quant_val),
             quant_val = round(quant_val,2)) %>% 
      group_by(batch, compound_name) %>%
      summarise(stdev = sd(quant_val, na.rm = T),
                average = mean(quant_val, na.rm = T),
                cv = (stdev / average)*100,
                upper_limit = average + 2.5*stdev,
                lower_limit = average - 2.5*stdev
      ) %>%
      pivot_longer(-c(batch, compound_name), names_to = "variable", values_to = "value") %>%
      mutate(variable = factor(variable, levels = c("upper_limit", "average", "lower_limit","stdev","cv")))
  })
  
  # output$TEST <- renderDataTable(
  #   meta5_4(),
  #   options = list(pageLength=5)
  # )
  
  # Need another summary dataframe to plot CVs and averages for geom_label
  meta5_5 <- reactive({
    
    compounds5 = tolower(unlist(strsplit(input$compounds5, split=",")))
    
    meta5() %>%
      filter(compound_name %in% compounds5) %>%
      inner_join(quant_conc_tbl5()) %>%
      mutate(compound_name=factor(compound_name,levels=compounds5)) %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(sampleid+compound_name+conc+letter ~ itsd,value.var="peakarea",
                      fun.aggregate=mean) %>%
      group_by(sampleid, letter) %>%
      fill(ITSD) %>%
      mutate(ITSD=zoo::na.locf(ITSD),
             norm_peak = peak / ITSD) %>%
      filter(!grepl("^CC[0-9]",sampleid),
             grepl(pattern = "[Pp][Ll][Aa][Ss][Mm][Aa]", sampleid)) %>%
      left_join(modelstart5()) %>%
      mutate(quant_val =  (norm_peak - (`(Intercept)`))/slope_value*as.numeric(input$xfactor5)) %>%
      mutate(sample = str_extract(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             qc = str_extract(sampleid, "[Qq][Cc]"),
             replicate = as.factor(str_extract(sampleid, "[0-9]+")),
             batch = unique(meta5_1()$batch)) %>%
      arrange(compound_name) %>%
      mutate(quant_val = ifelse(quant_val < 0,0,quant_val),
             quant_val = round(quant_val,2)) %>% 
      group_by(batch, compound_name) %>%
      summarise(stdev = sd(quant_val, na.rm = T),
                average = mean(quant_val, na.rm = T),
                cv = (stdev / average)*100,
                upper_limit = average + 2.5*stdev,
                lower_limit = average - 2.5*stdev
      )
  })
  
  # Build another dataframe to plot
  meta5_6 <- reactive({

    compounds5 = tolower(unlist(strsplit(input$compounds5, split=",")))

    meta5() %>%
      filter(compound_name %in% compounds5) %>%
      inner_join(quant_conc_tbl5()) %>%
      mutate(compound_name=factor(compound_name,levels=compounds5)) %>%
      replace_na(list(itsd="peak")) %>%
      reshape2::dcast(sampleid+compound_name+conc+letter ~ itsd,value.var="peakarea",
                      fun.aggregate=mean) %>%
      group_by(sampleid, letter) %>%
      fill(ITSD) %>%
      mutate(ITSD=zoo::na.locf(ITSD),
             norm_peak = peak / ITSD) %>%
      filter(!grepl("^CC[0-9]",sampleid),
             grepl(pattern = "[Pp][Ll][Aa][Ss][Mm][Aa]", sampleid)) %>%
      left_join(modelstart5()) %>%
      mutate(quant_val =  (norm_peak - (`(Intercept)`))/slope_value*as.numeric(input$xfactor5)) %>%
      mutate(sample = str_extract(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             qc = str_extract(sampleid, "[Qq][Cc]"),
             replicate = as.factor(str_extract(sampleid, "[0-9]+")),
             batch = unique(meta5_1()$batch)) %>%
      arrange(compound_name) %>%
      mutate(quant_val = ifelse(quant_val < 0,0,quant_val),
             quant_val = round(quant_val,2))

  })
  
  ba_plasma_plot = reactive({
    ggplot(meta5_6(), aes(x = batch, y = quant_val, shape = replicate, color = compound_name,
                          label = ifelse(is.na(quant_val),"NA",paste0(round(quant_val, digits = 1), "ug/mL")))) +
      geom_point(size = 3) +
      geom_point(meta5_5(), mapping = aes(x= batch, y = average),
                 color = "black", shape = 3, size = 3, inherit.aes = F) +
      geom_hline(subset(meta5_4(), variable %in%
                          c("average","upper_limit","lower_limit")),
                 mapping = aes(yintercept = value, linetype = variable))+
      ggrepel::geom_label_repel(meta5_5(),
                                mapping = aes(x= batch, y = average,
                                              label = ifelse(is.na(cv),"NA",paste0(round(cv, digits = 1), "%"))),
                                label.padding = 0.1,
                                direction = "y", max.overlaps = 50, min.segment.length = 5,
                                size = 1.2, inherit.aes = F) +
      ggrepel::geom_label_repel(label.padding = 0.1,
                                direction = "x", max.overlaps = 50, min.segment.length = 5,
                                size = 1.2, show.legend = FALSE) +
      theme_bw() +
      theme(panel.grid.minor= element_blank(),
            panel.grid.major.x = element_blank(),
            legend.text = element_text(color = "black", size = 8),
            legend.title = element_text(color = "black", size = 10, hjust = 0.5),
            legend.position = "top",
            strip.text=element_text(color = "black", size=5),
            axis.text =element_text(color = "black", size=8),
            axis.title = element_text(color = "black", size = 14),
            plot.title = element_text(color = "black", size = 15, hjust = 0.5, face = "bold"),
            plot.margin = margin(2,0.5,2,1.1, unit = 'cm')) +
      ggsci::scale_color_ucscgb() +
      guides(color = guide_legend(title = "Compound",
                                  override.aes = list(size = 2.5), ncol = 2,
                                  title.position="top",
                                  label.position = "right"), fill = F,
             shape = guide_legend(title = "Replicate",
                                  override.aes = list(size = 2.5), ncol = 1,
                                  title.position="top",
                                  label.position = "right"),
             linetype = guide_legend(ncol = 1,
                                     title.position="top",
                                     label.position = "right")) +
      scale_linetype_manual(name = "Metrics", values = c(1, 2, 1),
                            labels = c("+2.5 SD", "Mean", "-2.5 SD"),
                            guide = guide_legend(ncol = 1,
                                                 title.position="top",
                                                 label.position = "right"))+
      xlab("\nBatch Number") +
      ylab("Concentration (ug/mL)\n") +
      ggtitle("Plasma QC Summary\n") +
      facet_wrap(~compound_name, nrow = 1)
  })

  output$ba_plasma_plot <- renderPlot(
    ba_plasma_plot()
  )
  

  ## Save report ##
  output$ba_qc_report_download <- downloadHandler(
    filename = function(){
      paste0("BileAcid_QC_Report_",unique(meta5_2()$batch),"_",Sys.Date(),".pdf")
    },
    content = function(file) {
      pdf(file, height = 11, width = 8.5)
      print(
        gridExtra::grid.arrange(
          ba_qc_plot1() +
            xlab("") +
            ggtitle(paste("Bile Acid Quantitative QC Report\n", unique(meta5_1()$batch))) +
            theme(plot.title = element_text(color = "black",
                                            hjust = 0.5,
                                            size = 20,
                                            face = "bold")),
          ba_qc_plot2() + theme(legend.position = "none"),
          nrow = 2)
      )
      print(ba_plasma_plot())
      print(
        gridExtra::grid.arrange(
          ba_table_list1())
      )
      print(ba_qc_plot3())
      dev.off()
    }
  )

  ## Save model metrics in csv
  observeEvent(
    input$ba_cc_metrics_download,
    {write.csv(ba_model(), paste0("/Volumes/chaubard-lab/shiny_workspace/calibration_metrics/",
                                  gsub(".csv","",input$filename),"_CC_Metrics",".csv"))}
  )

  # Bile acid normalization tab -------------------------------------------------------

  rawdf_ba <- reactive({
    readin_bile_csv_single_file(file.path(wddir,input$filename),na.value = input$zero_val_bile_acid)
  })

  csv_ba <- reactive({

    vars_ba <- rawdf_ba() %>%
      mutate(compound_name=as.character(compound_name)) %>%
      arrange(desc(compound_name)) %>%
      filter(is.na(itsd)) %>%
      distinct(compound_name) %>%
      `$`(compound_name)

    return(vars_ba)
  })
  
  output$compound_list_bile_acid <- renderUI({
    checkboxGroupInput("check_compounds_bile_acid","Check = diluted",
                       choices=csv_ba(),
                       inline=F,
                       selected = csv_ba())
  })

  conc_int_bile_acid <- reactive({

    rawdf_ba() %>%
      filter(itsd=="ITSD") %>%
      group_by(Data.File, letter) %>%
      summarize(avg = mean(peakarea),
                med = median(peakarea))
  })
  
  # output$BATEST2 <- renderDataTable(
  #   conc_int_bile_acid()
  # )

  # conc_int_heatmap_bile_acid <- reactive({
  # 
  #   rawdf_ba() %>%
  #     filter(itsd=="ITSD") %>%
  #     group_by(compound_name) %>%
  #     summarize(min_peak = ifelse(all(peakarea == 0), 0, min(peakarea[peakarea != 0])))
  # 
  # })

  # output$conc_int_bile_acid <- renderDataTable(
  #   conc_int_bile_acid()
  # )

  #make boxplots:
  raw_boxplots_bile_acid <- reactive({

    rawdf_ba() %>%
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

  output$raw_boxplots_bile_acid <- renderPlot(
    raw_boxplots_bile_acid()
  )

  #make heatmap:
  heatmap_plot_bile_acid <- function()({
    rawdf_ba() %>%
      left_join(conc_filter_bile_acid()) %>%
      replace_na(list(checked="concentrated")) %>%
      filter(conc==checked, is.na(itsd),
             !grepl("(CC[0-9]+)", sampleid)) %>%
      select(-checked) %>%
      left_join(conc_int_bile_acid()) %>%
      mutate(norm_peak = peakarea / avg) %>%
      group_by(compound_name) %>% 
      mutate(compound_med = median(norm_peak)) %>% 
      ungroup() %>% 
      mutate(heat_val = log((norm_peak / compound_med), base = 2)) %>%
      separate(Data.File,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      mutate(sampleid = ifelse(str_detect(sampleid, "[Mm][Bb]|[Pp][Oo][Oo][Ll][Ee][Dd]|[Bb][Hh][Ii][Qq][Cc]|[Pp][Ll][Aa][Ss][Mm][Aa]|[Hh][Ee][Xx][Aa][Nn][Ee][Ss]|[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]|50%_[Mm][Ee][Oo][Hh]|[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]|50%[Mm][Ee][Oo][hh]"),
                               paste(num, sampleid, conc, sep = "__"),
                               sampleid)) %>%
      mutate(sampleid = gsub(pattern = "X", replacement = "", sampleid)) %>%
      dplyr::select(num, date, batch, sampleid, conc, compound_name, heat_val) %>%
      ungroup() %>%
      add_count(date,batch,sampleid, compound_name, conc) %>%
      mutate(sampleid = ifelse(n > 1, paste(num, sampleid, sep = "__"), sampleid),
             sampleid = gsub(pattern = "X", replacement = "", sampleid)) %>%
      dplyr::select(sampleid, compound_name, heat_val) %>%
      pivot_wider(names_from = compound_name, values_from = heat_val, values_fill = NA) %>%
      filter(!str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]")) %>%
      drop_na(.) %>%
      column_to_rownames(., var = "sampleid") %>%
      as.matrix(.) %>%
      t(.) %>%
      pheatmap(., scale = "none",
               cellheight = nrow(heatmap_data_bile_acid()) / (nrow(heatmap_data_bile_acid())*0.075),
               cellwidth = ncol(heatmap_data_bile_acid()) / (ncol(heatmap_data_bile_acid())*0.0825)+1.5,
               angle_col = "90",
               color=colorRampPalette(c("navy", "white", "red"))(50),
      )
  })

  output$heatmap_plot_bile_acid <- renderPlot(
    heatmap_plot_bile_acid()
  )

  heatmap_data_bile_acid <- function()({
    rawdf_ba() %>%
      left_join(conc_filter_bile_acid()) %>%
      replace_na(list(checked="concentrated")) %>%
      filter(conc==checked, is.na(itsd),
             !grepl("(CC[0-9]+)", sampleid)) %>%
      select(-checked) %>%
      left_join(conc_int_bile_acid()) %>%
      mutate(norm_peak = peakarea / avg) %>%
      group_by(compound_name) %>% 
      mutate(compound_med = median(norm_peak)) %>% 
      ungroup() %>% 
      mutate(heat_val = log((norm_peak / compound_med), base = 2)) %>%
      separate(Data.File,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      mutate(sampleid = ifelse(str_detect(sampleid, "[Mm][Bb]|[Pp][Oo][Oo][Ll][Ee][Dd]|[Bb][Hh][Ii][Qq][Cc]|[Pp][Ll][Aa][Ss][Mm][Aa]|[Hh][Ee][Xx][Aa][Nn][Ee][Ss]|[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]|50%_[Mm][Ee][Oo][Hh]|[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]|50%[Mm][Ee][Oo][hh]"),
                               paste(num, sampleid, conc, sep = "__"),
                               sampleid)) %>%
      mutate(sampleid = gsub(pattern = "X", replacement = "", sampleid)) %>%
      dplyr::select(num, date, batch, sampleid, conc, compound_name, heat_val) %>%
      ungroup() %>%
      add_count(date,batch,sampleid, compound_name, conc) %>%
      mutate(sampleid = ifelse(n > 1, paste(num, sampleid, sep = "__"), sampleid),
             sampleid = gsub(pattern = "X", replacement = "", sampleid)) %>%
      dplyr::select(sampleid, compound_name, heat_val) %>%
      pivot_wider(names_from = compound_name, values_from = heat_val, values_fill = NA) %>%
      filter(!str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]")) %>%
      drop_na(.) %>%
      column_to_rownames(., var = "sampleid") %>%
      as.matrix(.) %>%
      t(.)
  })
  
  output$BATEST1 <- renderDataTable(
    heatmap_data_bile_acid()
  )

  output$heatmap_download_bile_acid <- downloadHandler(
    filename = function(){
      paste0("normalized_heatmap_",input$filename,"_",Sys.Date(),".pdf")
    },
    content = function(file) {
      pdf(file,
          height = nrow(heatmap_data_bile_acid()) / (nrow(heatmap_data_bile_acid())*0.075),
          width = (ncol(heatmap_data_bile_acid()) / (ncol(heatmap_data_bile_acid())*0.07))+1.5
      )
      heatmap_plot_bile_acid()
      dev.off()
    },
    contentType = 'PDF'
  )


  # Intermediate table
  conc_filter_bile_acid <- reactive({

    if(length(as.list(input$check_compounds_bile_acid)) > 0){

      leftovervars_ba <- csv_ba()[ !(csv_ba() %in% input$check_compounds_bile_acid) ]

      tibble(checked = rep("diluted", length(input$check_compounds_bile_acid)),
             compound_name = input$check_compounds_bile_acid) %>%
        bind_rows(tibble(checked = rep("concentrated", length(leftovervars_ba)),
                         compound_name = leftovervars_ba)
        )
    }else{

      tibble(checked = "concentrated",
             compound_name = csv_ba())
    }
  })

  #subset based on checkboxes.. make a dataframe of compounds in same order as checkboxes with input as a column
  subset_conc_bile_acid <- reactive({
    rawdf_ba() %>%
      left_join(conc_filter_bile_acid()) %>%
      replace_na(list(checked="concentrated")) %>%
      filter(conc==checked) %>%
      select(-checked) %>%
      left_join(conc_int_bile_acid()) %>%
      mutate(norm_peak = peakarea / avg)
  })

  output$conc_filter_bile_acid <- renderDataTable(
    subset_conc_bile_acid() %>%
      datatable() %>%
      formatRound(c("norm_peak"), 3) %>%
      formatStyle(columns = c("norm_peak"), 'text-align' = 'center')
  )

  #download table
  #make wide table
  norm_wide_tbl_bile_acid <- reactive({

    if(input$qcfil_bile_acid==F){
      rawdf_ba() %>%
        left_join(conc_filter_bile_acid()) %>%
        replace_na(list(checked="concentrated")) %>%
        filter(conc==checked, is.na(itsd),
               !grepl("(CC[0-9]+)", sampleid)) %>%
        select(-checked) %>%
        left_join(conc_int_bile_acid()) %>%
        mutate(norm_peak = peakarea / avg) %>%
        separate(Data.File,into=c("num","date","batch","sampleid","conc"),
                 sep="__") %>%
        mutate(sampleid = ifelse(str_detect(sampleid, "[Mm][Bb]|[Pp][Oo][Oo][Ll][Ee][Dd]|[Bb][Hh][Ii][Qq][Cc]|[Pp][Ll][Aa][Ss][Mm][Aa]|[Hh][Ee][Xx][Aa][Nn][Ee][Ss]|[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]|50%_[Mm][Ee][Oo][Hh]|[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]|50%[Mm][Ee][Oo][hh]"),
                                 paste(num, sampleid, conc, sep = "__"),
                                 sampleid)) %>%
        mutate(sampleid = gsub(pattern = "X", replacement = "", sampleid)) %>%
        dplyr::select(num, date, batch, sampleid, conc, compound_name, norm_peak) %>%
        ungroup() %>%
        add_count(date,batch,sampleid, compound_name, conc) %>%
        mutate(sampleid = ifelse(n > 1, paste(num, sampleid, sep = "__"), sampleid),
               sampleid = gsub(pattern = "X", replacement = "", sampleid)) %>%
        dplyr::select(sampleid, compound_name, norm_peak) %>%
        pivot_wider(names_from = compound_name, values_from = norm_peak, values_fill = NA)
    }else{
      rawdf_ba() %>%
        left_join(conc_filter_bile_acid()) %>%
        replace_na(list(checked="concentrated")) %>%
        filter(conc==checked, is.na(itsd),
               !grepl("(CC[0-9]+)", sampleid)) %>%
        select(-checked) %>%
        left_join(conc_int_bile_acid()) %>%
        mutate(norm_peak = peakarea / avg) %>%
        separate(Data.File,into=c("num","date","batch","sampleid","conc"),
                 sep="__") %>%
        mutate(sampleid = ifelse(str_detect(sampleid, "[Mm][Bb]|[Pp][Oo][Oo][Ll][Ee][Dd]|[Bb][Hh][Ii][Qq][Cc]|[Pp][Ll][Aa][Ss][Mm][Aa]|[Hh][Ee][Xx][Aa][Nn][Ee][Ss]|[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]|50%_[Mm][Ee][Oo][Hh]|[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]|50%[Mm][Ee][Oo][hh]"),
                                 paste(num, sampleid, conc, sep = "__"),
                                 sampleid)) %>%
        mutate(sampleid = gsub(pattern = "X", replacement = "", sampleid)) %>%
        dplyr::select(num, date, batch, sampleid, conc, compound_name, norm_peak) %>%
        ungroup() %>%
        add_count(date,batch,sampleid, compound_name, conc) %>%
        mutate(sampleid = ifelse(n > 1, paste(num, sampleid, sep = "__"), sampleid),
               sampleid = gsub(pattern = "X", replacement = "", sampleid)) %>%
        dplyr::select(sampleid, compound_name, norm_peak) %>%
        pivot_wider(names_from = compound_name, values_from = norm_peak, values_fill = NA) %>%
        filter(!str_detect(sampleid, "[Mm][Bb]"),
               !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
               !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
               !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
               !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
               !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
               !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
               !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
               !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"))
    }
  })

  output$normwide_bile_acid <- renderDataTable(
    norm_wide_tbl_bile_acid() %>%
      datatable() %>%
      formatRound(c(2:ncol(norm_wide_tbl_bile_acid())), 3) %>%
      formatStyle(columns = c(2:ncol(norm_wide_tbl_bile_acid())), 'text-align' = 'center')
  )


  norm_wide_tbl_label_bile_acid <- reactive({
    if(input$qcfil_bile_acid==F){
      return("normalized_results_")
    } else {
      return("removed_qcs_normalized_results_")
    }
  })

  #download handler
  output$downloadData_bile_acid <- downloadHandler(

    filename = function(){
      paste0(norm_wide_tbl_label_bile_acid(),input$filename,"_",Sys.Date(),".csv")
    },

    content = function(file) {
      write.csv(norm_wide_tbl_bile_acid(),file,
                row.names=F,quote=F)
    }
  )



  ##### QC Report #####
  # ITSD Raw Peak Area #
  
  # Build a Norm CV dataframe of summary stats
  rawdf_ba2 <- reactive({

    rawdf_ba() %>%
      separate(Data.File,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      mutate(num = gsub(pattern = "X", replacement = "", num)) %>%
      filter(itsd == "ITSD",
             !str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
             !grepl("[Cc][Cc][0-9]+", sampleid)) %>%
      mutate(peakarea = ifelse(peakarea <= input$zero_val_bile_acid, input$zero_val_bile_acid, peakarea)) %>%
      group_by(batch, compound_name) %>%
      summarise(stdev = sd(peakarea),
                average = mean(peakarea),
                middle = median(peakarea),
                cv = stdev / average,
                cv_med = stdev / median(peakarea)# Don't turn into % here since it will be applied in the y-axis scale
      )
  })


  ba_norm_plot1 <-function()({

    rawdf_ba() %>%
      separate(Data.File,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      mutate(num = gsub(pattern = "X", replacement = "", num)) %>%
      filter(itsd == "ITSD",
             !str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]")) %>%
      mutate(peakarea = ifelse(peakarea <= input$zero_val_bile_acid, input$zero_val_bile_acid, peakarea),
             cc_shape = ifelse(grepl("CC[0-9]+", sampleid), "CC Sample", "ITSD")) %>%
      left_join(rawdf_ba2()) %>%
      mutate(num = as.numeric(num),
             flag = ifelse(peakarea > average + (1.5 * stdev) |
                             peakarea < average - (1.5 * stdev), paste(num,batch,sampleid, sep = "_"), NA)) %>%
      ggplot(., aes(x = num, y = peakarea, group = compound_name, label = flag)) +
      geom_point(aes(color = compound_name, shape = cc_shape), fill = "black", alpha = 0.6, size = 2) +
      geom_line(aes(y = average, color = compound_name)) +
      geom_line(aes(y = average + 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_line(aes(y = average - 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_ribbon(aes(ymin = average - stdev, ymax = average + stdev,
                      fill = compound_name), alpha=0.2) +
      ggrepel::geom_label_repel(size = 1.2, max.overlaps = Inf,
                                min.segment.length = 0.1, label.padding = 0.1) +
      theme_bw()+
      theme(panel.grid.minor= element_blank(),
            panel.grid.major.x = element_blank(),
            legend.text = element_text(color = "black", size = 8),
            legend.title = element_text(color = "black", size = 10),
            legend.position = "top",
            strip.text=element_text(color = "black", size=5),
            axis.text =element_text(color = "black", size=8),
            axis.title = element_text(color = "black", size = 10),
            plot.margin = margin(1,0.5,0,0.6, unit = 'cm')) +
      ggsci::scale_color_lancet() +
      ggsci::scale_fill_lancet() +
      guides(color = guide_legend(title = "Internal Standard Compound",
                                  override.aes = list(size = 2.5), nrow = 2,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right"), fill = F,
             shape = guide_legend(title = "",
                                  override.aes = list(size = 2.5), nrow = 2,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right")) +
      scale_shape_manual(values = c(24,16))+
      scale_y_continuous(label = scales::scientific) +
      scale_x_continuous(breaks = seq(0,150,25)) +
      ylab("Raw Peak Area\n") +
      xlab("\nInjection Number")+
      facet_wrap(~compound_name, scales="free_x", nrow = 3)
  })

  
  ba_norm_plot2 <-function()({

    rawdf_ba() %>%
      separate(Data.File,into=c("num","date","batch","sampleid","conc"),
  sep="__") %>%
  mutate(num = gsub(pattern = "X", replacement = "", num)) %>%
  filter(itsd == "ITSD",
         !str_detect(sampleid, "[Mm][Bb]"),
         !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
         !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
         !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
         !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
         !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
         !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
         !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
         !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]")) %>%
  mutate(peakarea = ifelse(peakarea <= input$zero_val_bile_acid, input$zero_val_bile_acid, peakarea),
         cc_shape = ifelse(grepl("CC[0-9]+", sampleid), "CC Sample", "ITSD")) %>%
      left_join(rawdf_ba2()) %>%
      mutate(num = as.numeric(num),
             flag = ifelse(peakarea > average + (1.5 * stdev) |
                             peakarea < average - (1.5 * stdev), paste(num,batch,sampleid, sep = "_"), NA)) %>%
      ggplot(., aes(x = num, y = peakarea, group = compound_name, label = flag)) +
      geom_point(aes(color = compound_name, shape = cc_shape), fill = "black", alpha = 0.6, size = 2) +
      geom_line(aes(y = average, color = compound_name)) +
      geom_line(aes(y = average + 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_line(aes(y = average - 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_ribbon(aes(ymin = ifelse(average - stdev <= 1, 1, average - stdev),
                      ymax = average + stdev,
                      fill = compound_name), alpha=0.2) +
      ggrepel::geom_label_repel(size = 1.2, max.overlaps = Inf,
                                min.segment.length = 0.1, label.padding = 0.1) +
      theme_bw() +
      theme(panel.grid.minor= element_blank(),
            panel.grid.major.x = element_blank(),
            legend.text = element_text(color = "black", size = 8),
            legend.title = element_text(color = "black", size = 10),
            legend.position = "top",
            strip.text=element_text(color = "black", size=5),
            axis.text =element_text(color = "black", size=8),
            axis.title = element_text(color = "black", size = 10),
            plot.margin = margin(1,0.5,0,0.6, unit = 'cm')) +
      ggsci::scale_color_lancet() +
      ggsci::scale_fill_lancet() +
      guides(color = guide_legend(title = "Internal Standard Compound",
                                  override.aes = list(size = 2.5), nrow = 2,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right"), fill = F,
             shape = guide_legend(title = "",
                                  override.aes = list(size = 2.5), nrow = 2,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right")) +
      scale_shape_manual(values = c(24,16))+
      ylab("Raw Peak Area\n(log10 scale)") +
      xlab("\nInjection Number") +
      # coord_trans(y = "log10") +
      scale_y_continuous(trans = "log10", labels = scales::scientific)+
      scale_x_continuous(breaks = seq(0,150,25))+
      facet_wrap(~compound_name, scales="free_x", nrow = 3)
  })
  

  # ITSD CV Percent #

  # Build CV dataframe
  rawdf_ba2_1 <- reactive({

    rawdf_ba() %>%
      separate(Data.File,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      mutate(num = gsub(pattern = "X", replacement = "", num)) %>%
      filter(itsd == "ITSD",
             !str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
             !grepl("[Cc][Cc][0-9]+", sampleid)) %>%
      mutate(peakarea = ifelse(peakarea <= input$zero_val_bile_acid, input$zero_val_bile_acid, peakarea)) %>%
      group_by(batch, compound_name) %>%
      summarise(stdev = sd(peakarea),
                average = mean(peakarea),
                middle = median(peakarea),
                cv = stdev / average,
                cv_med = stdev / median(peakarea)# Don't turn into % here since it will be applied in the y-axis scale
      ) %>%
      dplyr::rename(Batch = batch,
                    `Internal Standard` = compound_name,
                    StDev = stdev,
                    Mean = average,
                    Median = middle,
                    CV = cv,
                    `CV Median` = cv_med) %>%
      mutate(StDev = round(StDev, digits = 0),
             Mean = round(Mean, digits = 0),
             Median = round(Median, digits = 0),
             `CV (%)` = round(CV * 100, digits = 1),
             `CV Median (%)` = round(`CV Median` * 100, digits = 1)) %>%
      select(Batch,`Internal Standard`,StDev,Mean,Median,`CV (%)`,`CV Median (%)`)
  })
  
  
  # Build summary table
  ba_norm_table_list1 <- function()({
    gridExtra::tableGrob(rawdf_ba2_1(), rows = NULL, theme = tt)
  })

  ## Save Report ##
  output$ba_norm_qc_report_download <- downloadHandler(
    filename = function(){
      paste0("BileAcid_QC_Report_",unique(rawdf_ba2()$batch),"_",Sys.Date(),".pdf")
    },
    content = function(file) {
      pdf(file, height = 11, width = 8.5)
      print(ba_norm_plot1() +
              xlab("") +
              ggtitle(paste("Bile Acid Qualitative QC Report\n", unique(rawdf_ba2()$batch))) +
              theme(plot.title = element_text(color = "black",
                                              hjust = 0.5,
                                              size = 20,
                                              face = "bold"))
      )
      
      print(ba_norm_plot2() + theme(plot.title = element_blank()))
      
      print(gridExtra::grid.arrange(
        ba_norm_table_list1())
      )
      dev.off()
    }
  )
  
  
  
  # TMS normalization tab -------------------------------------------------------
  
  rawdf_tms <- reactive({
    readin_meta_csv_single_file_tms(file.path(wddir,input$filename),na.value = input$zero_val_tms)
  })
  
  # This generates a list of all compounds except for the internal standards
  csv_tms <- reactive({
    
    vars_tms <- rawdf_tms() %>%
      mutate(compound_name=as.character(compound_name)) %>%
      arrange(desc(compound_name)) %>%
      filter(is.na(itsd)) %>% # This filters FOR samples that are NOT internal standards
      distinct(compound_name) %>%
      `$`(compound_name) 
    
    return(vars_tms)
  })
  
  # This generates the average and median peak areas for all internal standards in a sample
  conc_int_tms <- reactive({
    
    compounds_tms = unlist(strsplit(input$compounds_tms, split=","))
    compounds_tms = factor(compounds_tms,level=unique(compounds_tms))
    
    rawdf_tms() %>%
      filter(itsd=="ITSD",
             compound_name %in% compounds_tms) %>%
      group_by(sampleid) %>%
      summarize(avg=mean(peakarea),
                med=median(peakarea))
    
  })
  
  
  conc_int_tms_sep <- reactive({
    conc_int_tms() %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="__")
  })
  
  
  # This generates the minimum peak areas for all internal standards in a sample
  # conc_int_heatmap_tms <- reactive({
  #   
  #   compounds_tms = unlist(strsplit(input$compounds_tms, split=","))
  #   compounds_tms = factor(compounds_tms,level=unique(compounds_tms))
  #   
  #   rawdf_tms() %>%
  #     filter(itsd=="ITSD",
  #            compound_name %in% compounds_tms) %>%
  #     group_by(compound_name) %>%
  #     summarize(min_peak = ifelse(all(peakarea == 0), 0, min(peakarea[peakarea != 0])))
  # })
  
  
  
  output$compound_list_tms <- renderUI({
    checkboxGroupInput("check_compounds_tms","Check = split1to5",
                       choices=csv_tms(),
                       inline=F)
  })
  
  #make boxplots:
  raw_boxplots_tms <- reactive({
    
    rawdf_tms() %>%
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
  
  output$raw_boxplots_tms <- renderPlot(
    raw_boxplots_tms()
  )
  
  #subset based on checkboxes.. make a dataframe of compounds in same order as checkboxes with input as a column
  conc_filter_tms <- reactive({
    
    if(length(as.list(input$check_compounds_tms)) > 0){
      
      leftovervars <- csv_tms()[ !(csv_tms() %in% input$check_compounds_tms) ]
      
      tibble(checked = rep("split1to5", length(input$check_compounds_tms)),
             compound_name = input$check_compounds_tms) %>%
        bind_rows(tibble(checked = rep("nosplit", length(leftovervars)),
                         compound_name = leftovervars)
        )
    }else{
      
      tibble(checked = "nosplit",
             compound_name = csv_tms())
    }
  })
  
  subset_conc_tms <- reactive({
    
    rawdf_tms() %>%
      left_join(conc_filter_tms()) %>%
      replace_na(list(checked="nosplit")) %>%
      filter(conc==checked) %>%
      select(-checked) %>%
      left_join(conc_int_tms()) %>%
      mutate(norm_peak=peakarea / avg)
  })
  
  output$conc_filter_tms <- renderDataTable(
    subset_conc_tms() %>%
      datatable() %>%
      formatRound(c("norm_peak"), 3) %>%
      formatStyle(columns = c("norm_peak"), 'text-align' = 'center')
  )
  
  
  #make heatmap:
  heatmap_plot_tms <- function()({
    rawdf_tms() %>%
      left_join(conc_filter_tms()) %>%
      replace_na(list(checked="nosplit"))%>%
      filter(conc==checked, is.na(itsd),
             !grepl("(__CC[0-9]+__)", sampleid)) %>%
      select(-checked) %>%
      left_join(conc_int_tms()) %>%
      # left_join(conc_int_heatmap_tms()) %>%
      mutate(norm_peak = peakarea / avg) %>%
      group_by(compound_name) %>% 
      mutate(compound_med = median(norm_peak)) %>%
      ungroup() %>%
      mutate(heat_val = log((norm_peak / compound_med), base = 2)) %>%
      dplyr::select(sampleid, compound_name, heat_val) %>%
      spread(compound_name, heat_val, fill = NA) %>%
      reshape2::melt(id.vars=c("sampleid")) %>%
      dplyr::rename(compound_name=variable,heat_val=value) %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="\\_\\_") %>%
      filter(!is.na(heat_val)) %>%
      mutate(sampleid = ifelse(str_detect(sampleid, "[Mm][Bb]|[Pp][Oo][Oo][Ll][Ee][Dd]|[Bb][Hh][Ii][Qq][Cc]|[Pp][Ll][Aa][Ss][Mm][Aa]|[Hh][Ee][Xx][Aa][Nn][Ee][Ss]|[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]|50%_[Mm][Ee][Oo][Hh]|[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]|50%[Mm][Ee][Oo][hh]"),
                               paste(num, sampleid, conc, sep = "__"),
                               sampleid)) %>%
      filter(!str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
             !grepl("CC[0-9]+", sampleid)) %>%
      dplyr::select(sampleid, compound_name, heat_val) %>%
      mutate(heat_val = round(heat_val,5)) %>%
      group_by(sampleid, compound_name) %>%
      summarise(heat_val = mean(heat_val)) %>%
      pivot_wider(names_from = compound_name, values_from = heat_val) %>%
      # filter_all(all_vars(!is.infinite(.))) %>%
      drop_na(.) %>%
      column_to_rownames(., var = "sampleid") %>%
      as.matrix(.) %>%
      t(.) %>%
      pheatmap(., scale = "none",
               cellheight = nrow(heatmap_data_tms()) / (nrow(heatmap_data_tms())*0.075),
               cellwidth = ncol(heatmap_data_tms()) / (ncol(heatmap_data_tms())*0.0825)+1.5,
               angle_col = "90",
               color=colorRampPalette(c("navy", "white", "red"))(50),
      )
  })
  
  # output$TMSTEST1 <- renderDataTable(
  #   heatmap_plot_tms()
  # )
  
  output$heatmap_plot_tms <- renderPlot(
    heatmap_plot_tms()
  )
  
  heatmap_data_tms <- function()({
    rawdf_tms() %>%
      left_join(conc_filter_tms()) %>%
      replace_na(list(checked="nosplit"))%>%
      filter(conc==checked, is.na(itsd),
             !grepl("(__CC[0-9]+__)", sampleid)) %>%
      select(-checked) %>%
      left_join(conc_int_tms()) %>%
      # left_join(conc_int_heatmap_tms()) %>%
      mutate(norm_peak = peakarea / avg) %>%
      group_by(compound_name) %>% 
      mutate(compound_med = median(norm_peak)) %>%
      ungroup() %>%
      mutate(heat_val = log((norm_peak / compound_med), base = 2)) %>%
      dplyr::select(sampleid, compound_name, heat_val) %>%
      spread(compound_name, heat_val, fill = NA) %>%
      reshape2::melt(id.vars=c("sampleid")) %>%
      dplyr::rename(compound_name=variable,heat_val=value) %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="\\_\\_") %>%
      filter(!is.na(heat_val)) %>%
      mutate(sampleid = ifelse(str_detect(sampleid, "[Mm][Bb]|[Pp][Oo][Oo][Ll][Ee][Dd]|[Bb][Hh][Ii][Qq][Cc]|[Pp][Ll][Aa][Ss][Mm][Aa]|[Hh][Ee][Xx][Aa][Nn][Ee][Ss]|[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]|50%_[Mm][Ee][Oo][Hh]|[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]|50%[Mm][Ee][Oo][hh]"),
                               paste(num, sampleid, conc, sep = "__"),
                               sampleid)) %>%
      filter(!str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
             !grepl("CC[0-9]+", sampleid)) %>%
      dplyr::select(sampleid, compound_name, heat_val) %>%
      mutate(heat_val = round(heat_val,5)) %>%
      group_by(sampleid, compound_name) %>%
      summarise(heat_val = mean(heat_val)) %>%
      pivot_wider(names_from = compound_name, values_from = heat_val) %>%
      # filter_all(all_vars(!is.infinite(.))) %>%
      drop_na(.) %>%
      column_to_rownames(., var = "sampleid") %>%
      as.matrix(.) %>%
      t(.)
  })
  
  # output$TEST4 <- renderDataTable(
  #   heatmap_data_tms()
  # )
  
  output$heatmap_download_tms <- downloadHandler(
    filename = function(){
      paste0("normalized_heatmap_",input$filename,"_",Sys.Date(),".pdf")
    },
    content = function(file) {
      pdf(file,
          height = nrow(heatmap_data_tms()) / (nrow(heatmap_data_tms())*0.035),
          width = (ncol(heatmap_data_tms()) / (ncol(heatmap_data_tms())*0.07))+1.5
      )
      heatmap_plot_tms()
      dev.off()
    },
    contentType = 'PDF'
  )
  
  
  #download table
  #make wide table
  norm_wide_tbl_tms <- reactive({
    
    if(input$qcfil_tms==F){
      rawdf_tms() %>%
        left_join(conc_filter_tms()) %>%
        replace_na(list(checked="nosplit")) %>%
        filter(conc==checked, is.na(itsd),
               !grepl("(__CC[0-9]+__)", sampleid)) %>%
        select(-checked) %>%
        left_join(conc_int_tms()) %>%
        mutate(norm_peak = peakarea / avg) %>%
        dplyr::select(sampleid, compound_name, norm_peak) %>%
        spread(compound_name, norm_peak, fill = NA) %>%
        reshape2::melt(id.vars=c("sampleid")) %>%
        dplyr::rename(compound_name=variable,norm_peak=value) %>%
        separate(sampleid,into=c("num","date","batch","sampleid","conc"),
                 sep="__") %>%
        filter(!is.na(norm_peak)) %>%
        mutate(sampleid = ifelse(str_detect(sampleid, "[Mm][Bb]|[Pp][Oo][Oo][Ll][Ee][Dd]|[Bb][Hh][Ii][Qq][Cc]|[Pp][Ll][Aa][Ss][Mm][Aa]|[Hh][Ee][Xx][Aa][Nn][Ee][Ss]|[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]|50%_[Mm][Ee][Oo][Hh]|[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]|50%[Mm][Ee][Oo][hh]"),
                                 paste(num, sampleid, conc, sep = "__"),
                                 sampleid)) %>%
        dplyr::select(sampleid, compound_name, norm_peak) %>%
        mutate(norm_peak = round(norm_peak,5)) %>%
        group_by(sampleid, compound_name) %>%
        summarise(norm_peak = mean(norm_peak)) %>%
        #     # reshape2::dcast(sampleid ~ compound_name, value.var = "norm_peak") %>%
        arrange(compound_name) %>%
        pivot_wider(names_from = compound_name, values_from = norm_peak)
      # select(-num)
    }else{
      rawdf_tms() %>%
        left_join(conc_filter_tms()) %>%
        replace_na(list(checked="nosplit")) %>%
        filter(conc==checked, is.na(itsd)) %>%
        select(-checked) %>%
        left_join(conc_int_tms()) %>%
        mutate(norm_peak=peakarea / avg) %>%
        dplyr::select(sampleid, compound_name, norm_peak) %>%
        spread(compound_name, norm_peak, fill = NA) %>%
        reshape2::melt(id.vars=c("sampleid")) %>%
        dplyr::rename(compound_name=variable,norm_peak=value) %>%
        separate(sampleid,into=c("num","date","batch","sampleid","conc"),
                 sep="__") %>%
        mutate(sampleid = ifelse(str_detect(sampleid, "[Mm][Bb]|[Pp][Oo][Oo][Ll][Ee][Dd]|[Bb][Hh][Ii][Qq][Cc]|[Pp][Ll][Aa][Ss][Mm][Aa]|[Hh][Ee][Xx][Aa][Nn][Ee][Ss]|[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]|50%_[Mm][Ee][Oo][Hh]|[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]|50%[Mm][Ee][Oo][hh]"),
                                 paste(num, sampleid, conc, sep = "__"),
                                 sampleid)) %>%
        filter(!is.na(norm_peak)) %>%
        filter(!str_detect(sampleid, "[Mm][Bb]"),
               !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
               !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
               !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
               !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
               !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
               !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
               !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
               !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
               !grepl("CC[0-9]+", sampleid)
        ) %>%
        dplyr::select(num,sampleid, compound_name, norm_peak) %>%
        mutate(norm_peak = round(norm_peak,5)) %>%
        group_by(sampleid, compound_name) %>%
        summarise(norm_peak = mean(norm_peak)) %>%
        # reshape2::dcast(num+sampleid ~ compound_name,value.var="norm_peak",fun.aggregate=mean) %>%
        arrange(compound_name) %>%
        pivot_wider(names_from = compound_name, values_from = norm_peak)
      # select(-num)
    }
  })
  
  output$normwide_tms <- renderDataTable(
    norm_wide_tbl_tms() %>%
      datatable() %>%
      formatRound(c(2:ncol(norm_wide_tbl_tms())), 3) %>%
      formatStyle(columns = c(2:ncol(norm_wide_tbl_tms())), 'text-align' = 'center')
  )
  
  
  norm_wide_tbl_label_tms <- reactive({
    if(input$qcfil_tms==F){
      return("normalized_results_")
    } else {
      return("removed_qcs_normalized_results_")
    }
  })
  
  #download handler
  output$downloadData_tms <- downloadHandler(
    
    filename = function(){
      paste0(norm_wide_tbl_label_tms(),input$filename,"_",Sys.Date(),".csv")
    },
    
    content = function(file) {
      write.csv(norm_wide_tbl_tms(),file,
                row.names=F,quote=F)
    }
  )
  
  ##### QC Report #####
  # ITSD Raw Peak Area #
  
  # Build a Norm CV dataframe of summary stats
  rawdf2_tms <- reactive({
    compounds = unlist(strsplit(input$compounds_tms, split=","))
    compounds = factor(compounds,level=unique(compounds))
    
    rawdf_tms() %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      inner_join(., conc_int_tms_sep(), by = c("num","date","batch","sampleid","conc")) %>%
      filter(itsd == "ITSD",
             compound_name %in% compounds,
             !str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
             !grepl("CC[0-9]+", sampleid)) %>%
      mutate(peakarea = ifelse(peakarea <= input$zero_val_tms, input$zero_val_tms, peakarea)) %>%
      group_by(batch, compound_name, conc) %>%
      summarise(stdev = sd(peakarea),
                average = mean(peakarea),
                middle = median(peakarea),
                cv = stdev / average,
                cv_med = stdev / median(peakarea)# Don't turn into % here since it will be applied in the y-axis scale
      )
  })
  
  tms_norm_plot1 <-function()({
    
    compounds = unlist(strsplit(input$compounds_tms, split=","))
    compounds = factor(compounds,level=unique(compounds))
    
    rawdf_tms() %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      inner_join(., conc_int_tms_sep(), by = c("num","date","batch","sampleid","conc")) %>%
      filter(itsd == "ITSD",
             compound_name %in% compounds,
             !str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]")) %>%
      mutate(peakarea = ifelse(peakarea <= input$zero_val_tms, input$zero_val_tms, peakarea),
             cc_shape = ifelse(grepl("[Cc][Cc][0-9]+", sampleid), "CC Sample", "ITSD")) %>%
      filter(itsd == "ITSD") %>%
      left_join(rawdf2_tms()) %>%
      mutate(num = as.numeric(num),
             flag = ifelse(peakarea > average + (1.5 * stdev) |
                             peakarea < average - (1.5 * stdev), paste(num,batch,sampleid, sep = "_"), NA)) %>%
      ggplot(., aes(x = num, y = peakarea, group = compound_name, label = flag)) +
      geom_point(aes(color = compound_name, shape = cc_shape), fill = "black", alpha = 0.6, size = 2) +
      geom_line(aes(y = average, color = compound_name)) +
      geom_line(aes(y = average + 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_line(aes(y = average - 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_ribbon(aes(ymin = average - stdev, ymax = average + stdev,
                      fill = compound_name), alpha=0.2) +
      ggrepel::geom_label_repel(size = 1.2, max.overlaps = Inf,
                                min.segment.length = 0.1, label.padding = 0.1) +
      theme_bw() +
      theme(panel.grid.minor= element_blank(),
            panel.grid.major.x = element_blank(),
            legend.text = element_text(color = "black", size = 13),
            legend.title = element_text(color = "black", size = 15),
            legend.position = "top",
            strip.text=element_text(color = "black", size=9),
            axis.text =element_text(color = "black", size=8),
            axis.title = element_text(color = "black", size = 15),
            plot.margin = margin(1,1,-1,1.1, unit = 'cm')) +
      ggsci::scale_color_uchicago() +
      ggsci::scale_fill_uchicago() +
      guides(color = guide_legend(title = "Internal Standard Compound",
                                  override.aes = list(size = 2.5), nrow = 1,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right"), fill = F,
             shape = guide_legend(title = "",
                                  override.aes = list(size = 2.5), nrow = 2,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right")) +
      scale_shape_manual(values = c(24,16))+
      scale_y_continuous(label = scales::scientific) +
      scale_x_continuous(breaks = seq(0,150,25)) +
      ylab("Raw Peak Area\n") +
      xlab("\nInjection Number")+
      facet_wrap(~compound_name+conc, scales="free_x", nrow = 2)
  })
  
  
  tms_norm_plot2 <-function()({
    
    compounds = unlist(strsplit(input$compounds_tms, split=","))
    compounds = factor(compounds,level=unique(compounds))
    
    rawdf_tms() %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      inner_join(., conc_int_tms_sep(), by = c("num","date","batch","sampleid","conc")) %>%
      filter(itsd == "ITSD",
             compound_name %in% compounds,
             !str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]")) %>%
      mutate(peakarea = ifelse(peakarea <= input$zero_val_tms, input$zero_val_tms, peakarea),
             cc_shape = ifelse(grepl("[Cc][Cc][0-9]+", sampleid), "CC Sample", "ITSD")) %>%
      filter(itsd == "ITSD") %>%
      left_join(rawdf2_tms()) %>%
      mutate(num = as.numeric(num),
             flag = ifelse(peakarea > average + (1.5 * stdev) |
                             peakarea < average - (1.5 * stdev), paste(num,batch,sampleid, sep = "_"), NA)) %>%
      ggplot(., aes(x = num, y = peakarea, group = compound_name, label = flag)) +
      geom_point(aes(color = compound_name, shape = cc_shape), fill = "black", alpha = 0.6, size = 2) +
      geom_line(aes(y = average, color = compound_name)) +
      geom_line(aes(y = average + 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_line(aes(y = average - 1.5*stdev, color = compound_name), linetype = "dashed") +
      geom_ribbon(aes(ymin = ifelse(average - stdev <= 1, 1, average - stdev),
                      ymax = average + stdev,
                      fill = compound_name), alpha=0.2) +
      ggrepel::geom_label_repel(size = 1.2, max.overlaps = Inf,
                                min.segment.length = 0.1, label.padding = 0.1) +
      theme_bw() +
      theme(panel.grid.minor= element_blank(),
            panel.grid.major.x = element_blank(),
            legend.text = element_text(color = "black", size = 13),
            legend.title = element_text(color = "black", size = 15),
            legend.position = "top",
            strip.text=element_text(color = "black", size=9),
            axis.text =element_text(color = "black", size=8),
            axis.title = element_text(color = "black", size = 15),
            plot.margin = margin(1,1,0.1,1.35, unit = 'cm')) +
      ggsci::scale_color_uchicago() +
      ggsci::scale_fill_uchicago() +
      guides(color = guide_legend(title = "Internal Standard Compound",
                                  override.aes = list(size = 2.5), nrow = 1,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right"), fill = F,
             shape = guide_legend(title = "",
                                  override.aes = list(size = 2.5), nrow = 2,
                                  title.position="top", title.hjust = 0.5,
                                  label.position = "right")) +
      scale_shape_manual(values = c(24,16))+
      ylab("Raw Peak Area\n(log10 scale)") +
      xlab("\nInjection Number") +
      # coord_trans(y = "log10") +
      scale_y_continuous(trans = "log10", labels = scales::scientific)+
      scale_x_continuous(breaks = seq(0,150,25))+
      facet_wrap(~compound_name+conc, scales="free_x", nrow = 2)
  })
  
  # ITSD CV Percent #
  
  # Build CV dataframe
  rawdf2_1_tms <- reactive({
    
    compounds = unlist(strsplit(input$compounds_tms, split=","))
    compounds = factor(compounds,level=unique(compounds))
    
    rawdf_tms() %>%
      separate(sampleid,into=c("num","date","batch","sampleid","conc"),
               sep="__") %>%
      inner_join(., conc_int_tms_sep(), by = c("num","date","batch","sampleid","conc")) %>%
      filter(itsd == "ITSD",
             compound_name %in% compounds,
             !str_detect(sampleid, "[Mm][Bb]"),
             !str_detect(sampleid, "[Pp][Oo][Oo][Ll][Ee][Dd]"),
             !str_detect(sampleid, "[Bb][Hh][Ii][Qq][Cc]"),
             !str_detect(sampleid, "[Pp][Ll][Aa][Ss][Mm][Aa]"),
             !str_detect(sampleid, "[Hh][Ee][Xx][Aa][Nn][Ee][Ss]"),
             !str_detect(sampleid, "[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd]"),
             !str_detect(sampleid, "50%_[Mm][Ee][Oo][Hh]"),
             !str_detect(sampleid, "[Ee][Aa]_[Bb][Ll][Aa][Nn][Kk]"),
             !str_detect(sampleid, "50%[Mm][Ee][Oo][Hh]"),
             !grepl("[Cc][Cc][0-9]+", sampleid)) %>%
      mutate(peakarea = ifelse(peakarea <= input$zero_val_tms, input$zero_val_tms, peakarea)) %>%
      group_by(batch, compound_name, conc) %>%
      summarise(stdev = sd(peakarea),
                average = mean(peakarea),
                middle = median(peakarea),
                cv = stdev / average,
                cv_med = stdev / median(peakarea)# Don't turn into % here since it will be applied in the y-axis scale
      ) %>%
      dplyr::rename(Batch = batch,
                    `Internal Standard` = compound_name,
                    StDev = stdev,
                    Mean = average,
                    Median = middle,
                    CV = cv,
                    `CV Median` = cv_med,
                    Concentration = conc) %>%
      mutate(StDev = round(StDev, digits = 0),
             Mean = round(Mean, digits = 0),
             Median = round(Median, digits = 0),
             `CV (%)` = round(CV * 100, digits = 1),
             `CV Median (%)` = round(`CV Median` * 100, digits = 1)) %>%
      select(Batch,Concentration,`Internal Standard`,StDev,Mean,Median,`CV (%)`,`CV Median (%)`)
  })
  
  
  
  # Build summary table
  tms_norm_table_list1 <- function()({
    gridExtra::tableGrob(rawdf2_1_tms(), rows = NULL, theme = tt)
  })
  
  ## Save Report ##
  output$norm_qc_report_download_tms <- downloadHandler(
    filename = function(){
      paste0("TMS_QC_Norm_Report_",unique(rawdf2_tms()$batch),"_",Sys.Date(),".pdf")
    },
    content = function(file) {
      pdf(file, height = 11, width = 8.5)
      print(
        gridExtra::grid.arrange(
          tms_norm_plot1() +
            xlab("") +
            ggtitle(paste("TMS Qualitative QC Report\n", unique(rawdf2_tms()$batch))) +
            theme(plot.title = element_text(color = "black",
                                            hjust = 0.5,
                                            size = 20,
                                            face = "bold")),
          tms_norm_plot2() + theme(legend.position = "none"),
          nrow = 2)
      )
      print(
        gridExtra::grid.arrange(
          tms_norm_table_list1())
      )
      dev.off()
    }
  )
  
# Instrument QC tab -----------------------------------------------------------------
  #### Calibration Curve-Slope Plot ####
  slope_qc_plot <- function()({
    
    if (input$method == "PFBBr") {
      finals_paths <- reactive({
        list.files(path = "/Volumes/chaubard-lab/shiny_workspace/calibration_metrics/",
                   pattern = input$method, full.names = TRUE)
      })
      
      # Read file content
      finals_content <- reactive({
        finals_paths() %>%
          lapply(read.table,
                 header = TRUE,
                 sep = ",",
                 encoding = "UTF-8")
      })
      
      # Read file name
      finals_filenames <- reactive({
        finals_paths() %>%
          basename() %>%
          as.list()
      })
      
      # Combine file content list and file name list
      finals_lists <- reactive({
        mapply(c, finals_content(), finals_filenames(), SIMPLIFY = FALSE)
      })
      
      # Unlist all lists and change column name
      finals_result <- reactive({
        data.table::rbindlist(finals_lists(), fill = T)
      })
      
      # Change column name and separate out to obtain batch info and date
      finals_result2 <- reactive({
        finals_result() %>%
          data.table::setDF() %>% 
          dplyr::rename(filename = V1) %>% 
          separate(filename, c("date","method","batch2","X2", "X3"), sep = "_") %>% 
          dplyr::select(-c(starts_with("X"), batch2)) %>% 
          dplyr::rename(itsd = Internal.Standard) %>% 
          dplyr::rename_all(tolower) %>% 
          dplyr::select(date, method, batch:slope)
      })
      
      p = reactive({
        ggplot(finals_result2(), aes(x = batch, y = slope, color = itsd)) +
          geom_point(size = 3) +
          theme_bw() +
          theme(panel.grid.minor= element_blank(),
                panel.grid.major.x = element_blank(),
                legend.text = element_text(color = "black", size = 12),
                legend.title = element_text(color = "black", size = 14),
                legend.position = "right",
                strip.text=element_text(color = "black", size=12),
                axis.text.y =element_text(color = "black", size=8),
                axis.text.x =element_text(color = "black", size=8, angle = 60, hjust = 1, vjust = 1),
                axis.title = element_text(color = "black", size = 12),
                plot.margin = margin(0.5,0.5,0.5,0.5, unit = 'cm')) +
          ggsci::scale_color_locuszoom() +
          xlab("\nBatch Number") +
          ylab("Calibration Curve Slope\n") +
          facet_wrap(~itsd, nrow = 1, scales = "free_y") +
          guides(color = guide_legend(title = "Compound",
                                      override.aes = list(size = 2.5), ncol = 1,
                                      title.position="top",
                                      label.position = "right"), fill = F)
      })
      
      return(p())
      
    } else if (input$method == "Indole") {      
      finals_paths <- reactive({
        list.files(path = "/Volumes/chaubard-lab/shiny_workspace/calibration_metrics/",
                   pattern = input$method, full.names = TRUE)
      })
      
      # Read file content
      finals_content <- reactive({
        finals_paths() %>%
          lapply(read.table,
                 header = TRUE,
                 sep = ",",
                 encoding = "UTF-8")
      })
      
      # Read file name
      finals_filenames <- reactive({
        finals_paths() %>%
          basename() %>%
          as.list()
      })
      
      # Combine file content list and file name list
      finals_lists <- reactive({
        mapply(c, finals_content(), finals_filenames(), SIMPLIFY = FALSE)
      })
      
      # Unlist all lists and change column name
      finals_result <- reactive({
        data.table::rbindlist(finals_lists(), fill = T)
      })
      
      # Change column name and separate out to obtain batch info and date
      finals_result2 <- reactive({
        finals_result() %>%
          data.table::setDF() %>% 
          dplyr::rename(filename = V1) %>% 
          separate(filename, c("date","method","batch2","X2", "X3"), sep = "_") %>% 
          dplyr::select(-c(starts_with("X"), batch2)) %>% 
          dplyr::rename(itsd = Internal.Standard) %>% 
          dplyr::rename_all(tolower) %>% 
          dplyr::select(date, method, batch:slope)
      })
      
      p = reactive({
        ggplot(finals_result2(), aes(x = batch, y = slope, color = itsd)) +
          geom_point(size = 3) +
          theme_bw() +
          theme(panel.grid.minor= element_blank(),
                panel.grid.major.x = element_blank(),
                legend.text = element_text(color = "black", size = 12),
                legend.title = element_text(color = "black", size = 14),
                legend.position = "right",
                strip.text=element_text(color = "black", size=12),
                axis.text.y =element_text(color = "black", size=8),
                axis.text.x =element_text(color = "black", size=8, angle = 60, hjust = 1, vjust = 1),
                axis.title = element_text(color = "black", size = 12),
                plot.margin = margin(0.5,0.5,0.5,0.5, unit = 'cm')) +
          ggsci::scale_color_igv() +
          xlab("\nBatch Number") +
          ylab("Calibration Curve Slope\n") +
          facet_wrap(~itsd, nrow = 3, scales = "free_y") +
          guides(color = guide_legend(title = "Compound",
                                      override.aes = list(size = 2.5), ncol = 1,
                                      title.position="top",
                                      label.position = "right"), fill = F)
      })
      
      return(p())
      
    } else {
      
      finals_paths <- reactive({
        list.files(path = "/Volumes/chaubard-lab/shiny_workspace/calibration_metrics/",
                   pattern = input$method, full.names = TRUE)
      })
      
      # Read file content
      finals_content <- reactive({
        finals_paths() %>%
          lapply(read.table,
                 header = TRUE,
                 sep = ",",
                 encoding = "UTF-8")
      })
      
      # Read file name
      finals_filenames <- reactive({
        finals_paths() %>%
          basename() %>%
          as.list()
      })
      
      # Combine file content list and file name list
      finals_lists <- reactive({
        mapply(c, finals_content(), finals_filenames(), SIMPLIFY = FALSE)
      })
      
      # Unlist all lists and change column name
      finals_result <- reactive({
        data.table::rbindlist(finals_lists(), fill = T)
      })
      
      # Change column name and separate out to obtain batch info and date
      finals_result2 <- reactive({
        finals_result() %>%
          data.table::setDF() %>% 
          dplyr::rename(filename = V1) %>% 
          separate(filename, c("date","method","batch2","X2", "X3"), sep = "_") %>% 
          dplyr::select(-c(starts_with("X"), batch2)) %>% 
          dplyr::rename(itsd = Internal.Standard) %>% 
          dplyr::rename_all(tolower) %>% 
          dplyr::select(date, method, batch:slope)
      })
      
      p = reactive({
        ggplot(finals_result2(), aes(x = batch, y = slope, color = itsd)) +
          geom_point(size = 3) +
          theme_bw() +
          theme(panel.grid.minor= element_blank(),
                panel.grid.major.x = element_blank(),
                legend.text = element_text(color = "black", size = 12),
                legend.title = element_text(color = "black", size = 14),
                legend.position = "right",
                strip.text=element_text(color = "black", size=12),
                axis.text.y =element_text(color = "black", size=8),
                axis.text.x =element_text(color = "black", size=8, angle = 60, hjust = 1, vjust = 1),
                axis.title = element_text(color = "black", size = 12),
                plot.margin = margin(0.5,0.5,0.5,0.5, unit = 'cm')) +
          ggsci::scale_color_ucscgb() +
          xlab("\nBatch Number") +
          ylab("Calibration Curve Slope\n") +
          facet_wrap(~itsd, nrow = 2, scales = "free_y")+
          guides(color = guide_legend(title = "Compound",
                                      override.aes = list(size = 2.5), ncol = 1,
                                      title.position="top",
                                      label.position = "right"), fill = F)
      })
      
      return(p())
      
    }
  })
  
  output$slope_qc_plot <- renderPlot(
    slope_qc_plot()
  )
  
  output$slope_qc_plot_download <- downloadHandler(
    filename = function(){
      paste0(input$method,"_", "InstrumentQC_CalibrationCurve_", Sys.Date(),".pdf")
    },
    content = function(file) {
      ggsave(file,plot = slope_qc_plot(),
             width =  if (input$method == "PFBBr") {
               get_row_col(plasma_qc_plot())[1] + 
                 length(unique(plasma_qc_plot()$data$batch))*1.25
             } else if (input$method == "Indole") {
               get_row_col(plasma_qc_plot())[1]*3.5 + 
                 length(unique(plasma_qc_plot()$data$batch))*0.5
             } else {
               get_row_col(plasma_qc_plot())[1]*2.5 + 
                 length(unique(plasma_qc_plot()$data$batch))*0.75
             },
             height = if (input$method == "PFBBr") {
               (get_row_col(plasma_qc_plot())[1] + 
                  length(unique(plasma_qc_plot()$data$batch))*1.25)*0.35
             } else if (input$method == "Indole") {
               (get_row_col(plasma_qc_plot())[1]*3.5 + 
                  length(unique(plasma_qc_plot()$data$batch))*0.5)*0.75
             } else {
               (get_row_col(plasma_qc_plot())[1]*2.5 + 
                  length(unique(plasma_qc_plot()$data$batch))*0.75)*0.5
             }
  )
    }
  )
  
  #### Plasma QC Plot ####
  plasma_qc_plot <- function()({
    
  if (input$method == "PFBBr") {
    finals_paths <- reactive({
      list.files(path = "/Volumes/chaubard-lab/shiny_workspace/CLIN_Finals_QCs/",
                 pattern = paste0("quant_results_[0-9]+_",input$method,"_\\w+.csv"), full.names = TRUE)
    })
    
    # Read file content
    finals_content <- reactive({
      finals_paths() %>%
        lapply(read.table,
               header = TRUE,
               sep = ",",
               encoding = "UTF-8")
    })
    
    # Read file name
    finals_filenames <- reactive({
      finals_paths() %>%
        basename() %>%
        as.list()
    })
    
    # Combine file content list and file name list
    finals_lists <- reactive({
      mapply(c, finals_content(), finals_filenames(), SIMPLIFY = FALSE)
    })
    
    # Unlist all lists and change column name
    finals_result <- reactive({
      data.table::rbindlist(finals_lists(), fill = T)
    })
    
    # Change column name and separate out to obtain batch info and date
    finals_result2 <- reactive({
      finals_result() %>%
        data.table::setDF() %>%
        dplyr::rename(filename = V1) %>%
        separate(filename, c("X1","X2","date","method","batch","X3","X4"), sep = "_") %>%
        select(method, date, batch, sampleid:Succinate)
    })
    
    # Build summary dataframe to produce CV values
    finals_summary <- reactive({
      finals_result2() %>% 
      mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
      pivot_longer(!c(method,date,batch,sampleid), names_to = "compound_name", values_to = "concentration") %>%
      filter(grepl("[Pp][Ll][Aa][Ss][Mm][Aa]", sampleid)) %>% 
      separate(sampleid, c("sample","qc","replicate"), sep = "_") %>% 
      group_by(method,date,batch,sample, compound_name) %>% 
      summarise(stdev = sd(concentration, na.rm = T),
                average = mean(concentration, na.rm = T),
                cv = (stdev / average)*100)
    })

    
    # Build summary specifically for high and low ranges for horizontal lines
    finals_summary2 <- reactive({
      finals_result2() %>%
      mutate_if(is.numeric, list(~na_if(., Inf))) %>%
      pivot_longer(!c(method,date,batch,sampleid), names_to = "compound_name", values_to = "concentration") %>%
      filter(grepl("[Pp][Ll][Aa][Ss][Mm][Aa]", sampleid)) %>%
      separate(sampleid, c("sample","qc","replicate"), sep = "_") %>%
      group_by(compound_name) %>%
      summarise(stdev = sd(concentration, na.rm = T),
                average = mean(concentration, na.rm = T),
                cv = (stdev / average)*100,
                upper_limit = average + 2.5*stdev,
                lower_limit = average - 2.5*stdev) %>%
      pivot_longer(-compound_name, names_to = "variable", values_to = "value") %>% 
        mutate(variable = factor(variable, levels = c("upper_limit", "average", "lower_limit","stdev","cv")))
    })
    

    # Build another dataframe to plot
    finals_result3 <- reactive({
      finals_result2() %>%
      mutate_if(is.numeric, list(~na_if(., Inf))) %>%
      pivot_longer(!c(method,date,batch,sampleid), names_to = "compound_name", values_to = "concentration") %>%
      filter(grepl("[Pp][Ll][Aa][Ss][Mm][Aa]", sampleid)) %>%
      separate(sampleid, c("sample","qc","replicate"), sep = "_") %>%
      ungroup()
    })

p = reactive({
    ggplot(finals_result3(), aes(x = batch, y = concentration, shape = replicate, color = compound_name)) +
      geom_point(size = 3) +
      geom_point(finals_summary(), mapping = aes(x= batch, y = average),
                 color = "black", shape = 3, size = 3, inherit.aes = F) +
      geom_hline(subset(finals_summary2(), variable %in%
                          c("average","upper_limit","lower_limit")),
                 mapping = aes(yintercept = value, linetype = variable))+
      ggrepel::geom_label_repel(finals_summary(),
                 mapping = aes(x= batch, y = average,
                               label = ifelse(is.na(cv),"NA",paste0(round(cv, digits = 1), "%"))),
                                label.padding = 0.1,
                                direction = "y", max.overlaps = 50, min.segment.length = 5,
                                size = 3, inherit.aes = F) +
      theme_bw() +
      theme(panel.grid.minor= element_blank(),
            panel.grid.major.x = element_blank(),
            legend.text = element_text(color = "black", size = 12),
            legend.title = element_text(color = "black", size = 14),
            legend.position = "right",
            strip.text=element_text(color = "black", size=12),
            axis.text.y =element_text(color = "black", size=8),
            axis.text.x =element_text(color = "black", size=8, angle = 60, hjust = 1, vjust = 1),
            axis.title = element_text(color = "black", size = 12),
            plot.margin = margin(0.5,0.5,0.5,0.5, unit = 'cm')) +
      ggsci::scale_color_locuszoom() +
      guides(color = guide_legend(title = "Compound",
                                  override.aes = list(size = 2.5), ncol = 1,
                                  title.position="top",
                                  label.position = "right"), fill = F,
             shape = guide_legend(title = "Replicate",
                                  override.aes = list(size = 2.5), ncol = 1,
                                  title.position="top",
                                  label.position = "right"),
             linetype = guide_legend(ncol = 1,
                                     title.position="top",
                                     label.position = "right")) +
      scale_linetype_manual(name = "Metrics", values = c(1, 2, 1),
                            labels = c("+2.5 SD", "Mean", "-2.5 SD"),
                            guide = guide_legend(ncol = 1,
                                                 title.position="top",
                                                 label.position = "right"))+
      xlab("\nBatch Number") +
      ylab("Concentration (mM)\n") +
      facet_wrap(~compound_name, nrow = 1, scales = "free_y")
})

return(p())

  } else if (input$method == "Indole") {
    finals_paths <- reactive({
      list.files(path = "/Volumes/chaubard-lab/shiny_workspace/CLIN_Finals_QCs/",
                 pattern = paste0("quant_results_[0-9]+_",input$method,"_\\w+.csv"), full.names = TRUE)
    })
    
    # Read file content
    finals_content <- reactive({
      finals_paths() %>%
        lapply(read.table,
               header = TRUE,
               sep = ",",
               encoding = "UTF-8")
    })
    
    # Read file name
    finals_filenames <- reactive({
      finals_paths() %>%
        basename() %>%
        as.list()
    })
    
    # Combine file content list and file name list
    finals_lists <- reactive({
      mapply(c, finals_content(), finals_filenames(), SIMPLIFY = FALSE)
    })
    
    # Unlist all lists and change column name
    finals_result <- reactive({
      data.table::rbindlist(finals_lists(), fill = T)
    })
    
    # Change column name and separate out to obtain batch info and date
    finals_result2 <- reactive({
      finals_result() %>%
        data.table::setDF() %>%
        dplyr::rename(filename = V1) %>%
        separate(filename, c("date","method","batch","X3","X4"), sep = "_") %>% 
        select(method, date, batch, sampleid:tyrosine)
    })

    # Build summary dataframe to produce CV values
    finals_summary <- reactive({
      finals_result2() %>% 
        mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
        pivot_longer(!c(method,date,batch,sampleid), names_to = "compound_name", values_to = "concentration") %>%
        filter(grepl("[Pp][Ll][Aa][Ss][Mm][Aa]", sampleid)) %>% 
        separate(sampleid, c("sample","qc","replicate"), sep = "_") %>% 
        group_by(method,date,batch,sample, compound_name) %>% 
        summarise(stdev = sd(concentration, na.rm = T),
                  average = mean(concentration, na.rm = T),
                  cv = (stdev / average)*100)
    })
    
    
    # Build summary specifically for high and low ranges for horizontal lines
    finals_summary2 <- reactive({
      finals_result2() %>%
        mutate_if(is.numeric, list(~na_if(., Inf))) %>%
        pivot_longer(!c(method,date,batch,sampleid), names_to = "compound_name", values_to = "concentration") %>%
        filter(grepl("[Pp][Ll][Aa][Ss][Mm][Aa]", sampleid)) %>%
        separate(sampleid, c("sample","qc","replicate"), sep = "_") %>%
        group_by(compound_name) %>%
        summarise(stdev = sd(concentration, na.rm = T),
                  average = mean(concentration, na.rm = T),
                  cv = (stdev / average)*100,
                  upper_limit = average + 2.5*stdev,
                  lower_limit = average - 2.5*stdev) %>%
        pivot_longer(-compound_name, names_to = "variable", values_to = "value") %>% 
        mutate(variable = factor(variable, levels = c("upper_limit", "average", "lower_limit","stdev","cv")))
    })
    
    
    # Build another dataframe to plot
    finals_result3 <- reactive({
      finals_result2() %>%
        mutate_if(is.numeric, list(~na_if(., Inf))) %>%
        pivot_longer(!c(method,date,batch,sampleid), names_to = "compound_name", values_to = "concentration") %>%
        filter(grepl("[Pp][Ll][Aa][Ss][Mm][Aa]", sampleid)) %>%
        separate(sampleid, c("sample","qc","replicate"), sep = "_") %>%
        ungroup()
    })
    
    p = reactive({
      ggplot(finals_result3(), aes(x = batch, y = concentration, shape = replicate, color = compound_name)) +
        geom_point(size = 3) +
        geom_point(finals_summary(), mapping = aes(x= batch, y = average),
                   color = "black", shape = 3, size = 3, inherit.aes = F) +
        geom_hline(subset(finals_summary2(), variable %in%
                            c("average","upper_limit","lower_limit")),
                   mapping = aes(yintercept = value, linetype = variable))+
        ggrepel::geom_label_repel(finals_summary(),
                                  mapping = aes(x= batch, y = average,
                                                label = ifelse(is.na(cv),"NA",paste0(round(cv, digits = 1), "%"))),
                                  label.padding = 0.1,
                                  direction = "y", max.overlaps = 50, min.segment.length = 5,
                                  size = 3, inherit.aes = F) +
        theme_bw() +
        theme(panel.grid.minor= element_blank(),
              panel.grid.major.x = element_blank(),
              legend.text = element_text(color = "black", size = 12),
              legend.title = element_text(color = "black", size = 14),
              legend.position = "right",
              strip.text=element_text(color = "black", size=12),
              axis.text.y =element_text(color = "black", size=8),
              axis.text.x =element_text(color = "black", size=8, angle = 60, hjust = 1, vjust = 1),
              axis.title = element_text(color = "black", size = 12),
              plot.margin = margin(0.5,0.5,0.5,0.5, unit = 'cm')) +
        ggsci::scale_color_igv() +
        guides(color = guide_legend(title = "Compound",
                                    override.aes = list(size = 2.5), ncol = 1,
                                    title.position="top",
                                    label.position = "right"), fill = F,
               shape = guide_legend(title = "Replicate",
                                    override.aes = list(size = 2.5), ncol = 1,
                                    title.position="top",
                                    label.position = "right"),
               linetype = guide_legend(ncol = 1,
                                       title.position="top",
                                       label.position = "right")) +
        scale_linetype_manual(name = "Metrics", values = c(1, 2, 1),
                              labels = c("+2.5 SD", "Mean", "-2.5 SD"),
                              guide = guide_legend(ncol = 1,
                                                   title.position="top",
                                                   label.position = "right"))+
        xlab("\nBatch Number") +
        ylab("Concentration (uM)\n") +
        facet_wrap(~compound_name, nrow = 3, scales = "free_y")
    })
    
    return(p())
    
    
  } else {
    
    finals_paths <- reactive({
      list.files(path = "/Volumes/chaubard-lab/shiny_workspace/CLIN_Finals_QCs/",
                 pattern = paste0("quant_results_[0-9]+_",input$method,"_\\w+.csv"), full.names = TRUE)
    })
    
    # Read file content
    finals_content <- reactive({
      finals_paths() %>%
        lapply(read.table,
               header = TRUE,
               sep = ",",
               encoding = "UTF-8")
    })
    
    # Read file name
    finals_filenames <- reactive({
      finals_paths() %>%
        basename() %>%
        as.list()
    })
    
    # Combine file content list and file name list
    finals_lists <- reactive({
      mapply(c, finals_content(), finals_filenames(), SIMPLIFY = FALSE)
    })
    
    # Unlist all lists and change column name
    finals_result <- reactive({
      data.table::rbindlist(finals_lists(), fill = T)
    })
    
    # Change column name and separate out to obtain batch info and date
    finals_result2 <- reactive({
      finals_result() %>% 
        data.table::setDF() %>% 
        dplyr::rename(filename = V1) %>% 
        separate(filename, c("date","method","batch","X4"), sep = "_") %>% 
        select(date,method,batch, Sample.Name:X08_3.Oxolithocholic.Acid_E) %>% 
        separate(Sample.Name, c("num","X2","X3","sampleid"), sep = "__") %>% 
        select(-c(X2,X3))
    })
    
    
    # Build summary dataframe to produce CV values
    finals_summary <- reactive({
      finals_result2() %>%
        mutate_if(is.numeric, list(~na_if(., Inf))) %>%
        pivot_longer(!c(date,method,batch,sampleid,num), names_to = "compound_name", values_to = "concentration") %>%
        filter(grepl("[Pp][Ll][Aa][Ss][Mm][Aa]", sampleid)) %>%
        mutate(compound_name = gsub("X[0-9]+_","",compound_name),
               compound_name = gsub("_[A-Z]","",compound_name),
               compound_name = gsub(".", " ", compound_name)) %>%
        separate(sampleid, c("sample","qc","replicate"), sep = "_") %>%
        group_by(method,date,batch,sample, compound_name) %>%
        summarise(stdev = sd(concentration, na.rm = T),
                  average = mean(concentration, na.rm = T),
                  cv = (stdev / average)*100)
    })

    # Build summary specifically for high and low ranges for horizontal lines
    finals_summary2 <- reactive({
      finals_result2() %>%
        mutate_if(is.numeric, list(~na_if(., Inf))) %>%
        pivot_longer(!c(date,method,batch,sampleid,num), names_to = "compound_name", values_to = "concentration") %>%
        filter(grepl("[Pp][Ll][Aa][Ss][Mm][Aa]", sampleid)) %>%
        separate(sampleid, c("sample","qc","replicate"), sep = "_") %>%
        group_by(compound_name) %>%
        summarise(stdev = sd(concentration, na.rm = T),
                  average = mean(concentration, na.rm = T),
                  cv = (stdev / average)*100,
                  upper_limit = average + 2.5*stdev,
                  lower_limit = average - 2.5*stdev) %>%
        pivot_longer(-compound_name, names_to = "variable", values_to = "value") %>%
        mutate(variable = factor(variable, levels = c("upper_limit", "average", "lower_limit","stdev","cv")))
    })


    # Build another dataframe to plot
    finals_result3 <- reactive({
      finals_result2() %>%
        mutate_if(is.numeric, list(~na_if(., Inf))) %>%
        pivot_longer(!c(date,method,batch,sampleid,num), names_to = "compound_name", values_to = "concentration") %>%
        filter(grepl("[Pp][Ll][Aa][Ss][Mm][Aa]", sampleid)) %>%
        separate(sampleid, c("sample","qc","replicate"), sep = "_") %>%
        ungroup()
    })

    p = reactive({
      ggplot(finals_result3(), aes(x = batch, y = concentration, shape = replicate, color = compound_name)) +
        geom_point(size = 3) +
        geom_point(finals_summary(), mapping = aes(x= batch, y = average),
                   color = "black", shape = 3, size = 3, inherit.aes = F) +
        geom_hline(subset(finals_summary2(), variable %in%
                            c("average","upper_limit","lower_limit")),
                   mapping = aes(yintercept = value, linetype = variable))+
        ggrepel::geom_label_repel(finals_summary(),
                                  mapping = aes(x= batch, y = average,
                                                label = ifelse(is.na(cv),"NA",paste0(round(cv, digits = 1), "%"))),
                                  label.padding = 0.1,
                                  direction = "y", max.overlaps = 50, min.segment.length = 5,
                                  size = 3, inherit.aes = F) +
        theme_bw() +
        theme(panel.grid.minor= element_blank(),
              panel.grid.major.x = element_blank(),
              legend.text = element_text(color = "black", size = 12),
              legend.title = element_text(color = "black", size = 14),
              legend.position = "right",
              strip.text=element_text(color = "black", size=12),
              axis.text.y =element_text(color = "black", size=8),
              axis.text.x =element_text(color = "black", size=8, angle = 60, hjust = 1, vjust = 1),
              axis.title = element_text(color = "black", size = 12),
              plot.margin = margin(0.5,0.5,0.5,0.5, unit = 'cm')) +
        ggsci::scale_color_igv() +
        guides(color = guide_legend(title = "Compound",
                                    override.aes = list(size = 2.5), ncol = 1,
                                    title.position="top",
                                    label.position = "right"), fill = F,
               shape = guide_legend(title = "Replicate",
                                    override.aes = list(size = 2.5), ncol = 1,
                                    title.position="top",
                                    label.position = "right"),
               linetype = guide_legend(ncol = 1,
                                       title.position="top",
                                       label.position = "right")) +
        scale_linetype_manual(name = "Metrics", values = c(1, 2, 1),
                              labels = c("+2.5 SD", "Mean", "-2.5 SD"),
                              guide = guide_legend(ncol = 1,
                                                   title.position="top",
                                                   label.position = "right"))+
        xlab("\nBatch Number") +
        ylab("Concentration (ug/mL)\n") +
        facet_wrap(~compound_name, nrow = 2, scales = "free_y")
    })

    return(p())
    
  }
  })
  
  output$plasma_qc_plot <- renderPlot(
    plasma_qc_plot()
  )
  
  output$plasma_qc_plot_download <- downloadHandler(
    filename = function(){
      paste0(input$method,"_", "InstrumentQC_PlasmaQC_", Sys.Date(),".pdf")
    },
    content = function(file) {
      ggsave(file,plot = plasma_qc_plot(),
             width =  if (input$method == "PFBBr") {
               get_row_col(plasma_qc_plot())[1] + 
                 length(unique(plasma_qc_plot()$data$batch))*1.25
             } else if (input$method == "Indole") {
               get_row_col(plasma_qc_plot())[1]*3.5 + 
                 length(unique(plasma_qc_plot()$data$batch))*0.5
             } else {
               get_row_col(plasma_qc_plot())[1]*2.5 + 
                 length(unique(plasma_qc_plot()$data$batch))*0.75
             },
             height = if (input$method == "PFBBr") {
               (get_row_col(plasma_qc_plot())[1] + 
                  length(unique(plasma_qc_plot()$data$batch))*1.25)*0.45
             } else if (input$method == "Indole") {
               (get_row_col(plasma_qc_plot())[1]*3.5 + 
                  length(unique(plasma_qc_plot()$data$batch))*0.5)*0.75
             } else {
               (get_row_col(plasma_qc_plot())[1]*2.5 + 
                  length(unique(plasma_qc_plot()$data$batch))*0.75)*0.5
             }
             )
    }
  )
  
  

}

# run app -----------------------------------------------------------------

# Run locally
# shinyApp(ui, server)

# Run on supermac
runApp(list(ui=ui,server=server),host="0.0.0.0",port=5000)
