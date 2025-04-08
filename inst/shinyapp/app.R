


require(shiny)
require(shinyjs)
require(ggpubr)

# runShinyApp("/Users/remy.nicolle/Workspace/PDAC_aux/TNERadio/smallpanetradio/inst/shinyapp")

shinyApp(
  ui <- fluidPage(
    useShinyjs(),
    
    
    titlePanel("Radiomic prediction of small TNE aggressivity"),
    
    
    sidebarLayout(
      
      sidebarPanel(
        h2("Input radiomic features"),
        fileInput("upart", "Upload arterial phase",accept="xlsx"),
        textOutput("artmsg", container = tags$b),
        
        fileInput("uport", "Upload portal phase",accept="xlsx"),
        textOutput("portmsg", container = tags$b),
        
        
        actionButton("go", "predict"),

      downloadButton("downloadData", "Download Results")
        
        
      ),
      
      
      mainPanel(
        
        plotOutput("plotboth"),
        plotOutput("plotart"),
        plotOutput("plotport"),
        textOutput("message", container = tags$b),
        textOutput("error", container = tags$b)
      )
    )
  ),


  server <- function(input, output,session) {
    
    getdata=eventReactive(input$go, {
      
      artfile=input$upart
      portfile=input$uport

      artmsg="No file"
      portmsg="No file"
      

      artdata=NULL
      if(!is.null(artfile$datapath)){
        # artdata=readxl::read_xlsx(artfile$datapath) 


        tryCatch({
              # artdata=radiomicXlHelper(xlfile,"art",throwError=T)
              artdata=radiomicXlHelper(artfile$datapath,"art",throwError=T)
            },  error = function(e) {
              print(paste("Error:", e$message))
              artmsg=paste("Error Arterial phase file:", e$message)
            })   
        if(!is.null(artdata)){
          artmsg="Arterial phase file loaded"
        }
      }

      print(artmsg)
      portdata=NULL

      if(!is.null(portfile$datapath)){
        tryCatch({
              portdata=radiomicXlHelper(portfile$datapath,"port",throwError=T)
            }, error = function(e) {
              print(paste("Error:", e$message))
              artmsg=paste("Error Portal phase file:", e$message)
            })   
        if(!is.null(portdata)){
          portmsg="Portal phase file loaded"
        }
      }

      print(portmsg)
      
      
      refbothpred=NULL
      newbothpred=NULL
      if(!is.null(artdata) & !is.null(portdata)){
        newbotgdata=data.frame(artdata,portdata)
        coms=intersect(rownames(smallpanetradio::refartdata),rownames(smallpanetradio::refportdata))
        bothref=data.frame(smallpanetradio::refartdata[coms,],smallpanetradio::refportdata[coms,])
        
        
        newbothpred=radiopred(newbotgdata,smallpanetradio::MODELS$bothmodel,bothref)
        refbothpred=radiopred(bothref,smallpanetradio::MODELS$bothmodel)
        refbothpred=data.frame(smallpanetradio::refclin[rownames(refbothpred),],refbothpred)

      }

      artpred=NULL;refartpred=NULL
      if(!is.null(artdata)){
        print("Yes there is artdata")
        print(dim(artdata))
        

        if(nrow(artdata)==1){
          print("One row, so treated as a single patient")
          
          artpred=radiopred(artdata[1,],smallpanetradio::MODELS$artmodel,smallpanetradio::refartdata)
        }else{
           artpred=radiopred(artdata,smallpanetradio::MODELS$artmodel,smallpanetradio::refartdata)
        }
        
        refartpred=radiopred(smallpanetradio::refartdata,smallpanetradio::MODELS$artmodel)
        refartpred=data.frame(smallpanetradio::refclin[rownames(refartpred),],refartpred)
           
      }

      portpred=NULL;refportpred=NULL
        if(!is.null(portdata)){
        if(nrow(portdata)==1){
                  print("One row, so treated as a single patient")
            portpred=radiopred(portdata[1,],smallpanetradio::MODELS$portmodel,smallpanetradio::refportdata)
        }else{

          portpred=radiopred(portdata,smallpanetradio::MODELS$portmodel,smallpanetradio::refportdata)
        }
        refportpred=radiopred(smallpanetradio::refportdata,smallpanetradio::MODELS$portmodel)
        refportpred=data.frame(smallpanetradio::refclin[rownames(refportpred),],refportpred)
           
      }

      alldata=list(artdata=artdata,artmsg=artmsg,portdata=portdata,portmsg=portmsg,
                   bothpred=newbothpred,refbothpred=refbothpred,
                   artpred=artpred,refartpred=refartpred,
                   portpred=portpred,refportpred=refportpred
                   )      
      alldata
    })
    
    
    output$plotboth <- renderPlot({ 
       alldata=getdata();
      if(!is.null(alldata$bothpred) & !is.null(alldata$refbothpred)){
        smallpanetradio:::.plotrez(alldata$bothpred,alldata$refbothpred,"both")
      }else{
         par(mar = c(0, 0, 0, 0))    
          plot(x = 0:10, y = 0:10, ann = F,bty = "n",type = "n",
         xaxt = "n", yaxt = "n")
          text(x = 5,y = 5,paste0("Could not retrieve data for Both phases"),cex=2)
      }
    })

    output$plotart <- renderPlot({
      
      alldata=getdata();
      if(!is.null(alldata$artpred)&!is.null(alldata$refartpred)){
       
        smallpanetradio:::.plotrez(alldata$artpred,alldata$refartpred,"art")
      }else{
        par(mar = c(0, 0, 0, 0))    
        plot(x = 0:10, y = 0:10, ann = F,bty = "n",type = "n",
             xaxt = "n", yaxt = "n")
        text(x = 5,y = 5,paste0("Could not retrieve data for Arterial phase"),cex=2)
      }
      
    })
    output$plotport <- renderPlot({  alldata=getdata();
          
      alldata=getdata();
      if(!is.null(alldata$portpred)&!is.null(alldata$refportpred)){
         

        smallpanetradio:::.plotrez(alldata$portpred,alldata$refportpred,"port")
      }else{
        par(mar = c(0, 0, 0, 0))    
        plot(x = 0:10, y = 0:10, ann = F,bty = "n",type = "n",
             xaxt = "n", yaxt = "n")
        text(x = 5,y = 5,paste0("Could not retrive data for Portal phase"),cex=2)
      }
    })
    
    output$error <- renderText({
      alldata <- getdata()
      # return(nrow(alldata$art.num))
      alldata$msg
    })
    
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("data-", Sys.Date(), ".xlsx", sep = "")
      },
      content = function(con) {
        alldata=getdata()

        openxlsx::write.xlsx(list(
          arterial = alldata$artpred,
          portal = alldata$portpred,
          both = alldata$bothpred
        ), con,rowNames=T)
      }
    )
      
    
    output$artmsg <- renderText({getdata()$artmsg})
    output$portmsg <- renderText({getdata()$portmsg })
    
    
    

    
    
    
    
    
    
    
    
  }
)





