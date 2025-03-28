

runShinyApp=function(){

  require(shiny)
  require(shinyjs)
  require(ggpubr)

      
  app_dir <- system.file("shinyapp", package = "smalltneradio")
  shiny::runApp(app_dir) 
}



# runShinyApp()/Users/remy.nicolle/Workspace/PDAC_aux/TNERadio/smalltneradio/inst/shinyapp/app.R
