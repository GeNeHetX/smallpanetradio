

runShinyApp=function(){

  require(shiny)
  require(shinyjs)
  require(ggpubr)

      
  
  shiny::runApp( system.file("shinyapp", package = "smalltneradio")) 
}

