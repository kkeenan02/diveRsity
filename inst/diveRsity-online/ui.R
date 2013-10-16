library(shiny)
library(diveRsity)



# Define UI
shinyUI(pageWithSidebar(
  
  #app title
  headerPanel("diveRsity Online"),
  
  #input
  sidebarPanel(
    
    fileInput("file", "Input file", multiple = FALSE, accept = NULL),
            
    numericInput("gp", "Genepop format", value = 3, min = 2, max = 3, step = 1),
    
    checkboxInput("pairwise", "Pairwise stats", FALSE),
    
    checkboxInput("WC_Fst", "Weir & Cockerham Fst", FALSE),
    
    checkboxInput("bs_locus", "Locus Bootstrap", FALSE),
    
    checkboxInput("bs_pairwise", "Pairwise Bootstraps", FALSE),
    
    checkboxInput("divBasic", 
                  "Calculate basic stats (e.g. HWE, allelic richness)", FALSE),
    
    numericInput("bootstraps", "Number of Bootstraps", value = 0, min = 0, 
                 max = 1000, step = 10),
    
    checkboxInput("parallel", "Parallel", FALSE),
    
    helpText(""),
    helpText("If you would like to plot your locus values",
             "using the corPlot function, please indicate below."),
    
    checkboxInput("corplot", "run corPlot?", FALSE),
    
    
    actionButton("goButton", "Calculate"),
    
    helpText(""),
    
    helpText("Written and designed by Kevin Keenan using shiny",
             "from RStudio and Inc. (2012)."),
    helpText("Any suggestion or questions should be directed to,"),
    helpText("<kkeenan02 AT qub.ac.uk>")
  )  
    
    ,
  
  mainPanel(
    tabsetPanel(
    tabPanel("Welcome",
             helpText("NOTE: This web app is designed for small to",
                      " medium data sets only (e.g.< 15 pops x 50 loci).",
                      "Users wishing to analyse larger data sets are",
                      "advised to use the original diveRsity R package.",
                      "This package provides much more flexibility than",
                      "the web app, as well as providing pairwise population",
                      "bootstrapping facilities."),
             helpText(""),
             
             a(href = "http://goo.gl/ahGTo", 
               "Click here to download the help document"),
             
             helpText(""),
             helpText("For more information about the diveRsity package please",
                      "visit"),
             a(href = "http://goo.gl/ikWmh", 
               "http://diversityinlife.weebly.com/software.html"),
             
             helpText(""),
             helpText("For information on the shiny package, please visit the",
                      "RStudio website at"),
             a(href = "http://www.rstudio.com/shiny/", "http://www.rstudio.com/shiny/")                 
    ),
    tabPanel("divBasic", downloadButton("divBdl", "Download as file"),
             tableOutput("divB")),
    tabPanel("Standard", downloadButton("dlstd", "Download as file"),
             tableOutput("std")),
    tabPanel("Estimates", downloadButton("dlest", "Download as file"), 
             tableOutput("est")),
    tabPanel("Pairwise", downloadButton("dlpw", "Download as file"), 
             tableOutput("pw")),
    tabPanel("Locus Bootstrap",  downloadButton("dllcbs", "Download as file"),
             tableOutput("bs_loc")),
    tabPanel("Pairwise Bootstrap", downloadButton("dlpwbs", "Download as file"), 
             tableOutput("pw_bs")),
    tabPanel("corPlot", downloadButton("corplt", "Download as PDF"),
             plotOutput("cor", width = "90%", height = "550px"))
  ))
))