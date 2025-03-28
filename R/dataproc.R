

checkRadiomicdataframe=function(df,phasename,verbose=TRUE,throwError=FALSE){

  thisvarnames=unique(c(
    list(art=smallpanetradio::MODELS$artvars,port=smallpanetradio::MODELS$portvars)[[phasename]],
    grep(paste0("^",phasename),smallpanetradio::MODELS$bothvars,value=TRUE)
    ))


  msg=""
  if(is.null(df)){
    msg="no data"
    return(msg)
  }

  missing=setdiff(thisvarnames,colnames(df))
    OK=FALSE

  if(length(missing)==0){
    msg=paste(msg,"\ngot no missing variables")
    if(verbose)print(msg)
  }else{
    msg=paste(msg,"missing vars for ",phasename)
    if(verbose)print(msg)
  }
  
  if(!all(thisvarnames %in% colnames(df))){
     msg=paste(msg,"mssing all variable  in ",phasename)
    if(verbose)print(msg)

    }else if(length(missing)>=(length(thisvarnames)-1)){
      msg=paste(msg,"Missing all variables, file format likely inadequate")
      if(verbose)print(msg)
    }else if(length(missing)>5 ){
      msg=paste(msg,"mssing a llot of variable  in ",phasename)
      if(verbose)print(msg)

    }else if(length(missing)>0){

      msg=paste(msg,"Missing following variables:",
        paste(missing, collapse=" , "))
      if(verbose)print(msg)

    }
    
    dfnum=data.matrix(df[,thisvarnames])
    
    
    if(sum(is.na(dfnum[1,]) )>0 ){
      msg=paste(msg,"Variable names (column names) OK. \
                        Missing radiomic values")
      if(verbose)print(msg)          

    }else{
    
      msg=paste(msg,"Data OK")
      if(verbose)print(msg)          
      OK=TRUE
    }
    if(throwError & !OK){
      stop(msg)
    }

    

    return(dfnum) 

}

radiomicXlHelper=function(file,phasename,xltab=1,verbose=TRUE,throwError=FALSE){
  rawtab=openxlsx::read.xlsx(file,sheet=xltab)
  coli=4:ncol(rawtab)
  df=data.frame(t(rawtab[,coli]))
  # if(is.null(predefineColnames)){
  df=setNames(df,paste0(phasename,rawtab[,2],"_",rawtab[,3]))
  # }else{
  #   df=setNames(df,paste0(phasename,predefineColnames))
  # }
  
  suppressWarnings({
    df[] <- lapply(df, function(x) if(!any(is.na(as.numeric(as.character(x))))) as.numeric(as.character(x)) else x)
  })
  rownames(df) <- colnames(rawtab)[coli]
  
  checkRadiomicdataframe(df,phasename,verbose,throwError)
  
}



# proj2eval=selpcaproj
# refproj=refpcaproj
projAnomalyEval=function(proj2eval,refproj,alpha = c(0.05,10^-(2:5)) ){

  alpha <- c(0.05,10^-(2:5)) 


  if(is.vector(refproj)){

    n=length(refproj) # Sample size
    x_bar <- mean(refproj) # Sample mean
    s <- sd(refproj)       # Sample standard deviation
    df <- n - 1   # Degrees of freedom
    t_crit <- qt(1 - alpha / 2, df = df) # Two-tailed t-critical value
    me_pred <- t_crit * s * sqrt(1 + 1/n)

    centdist=abs((proj2eval - x_bar)) 
    return(do.call(rbind,lapply(centdist, function(x) {

    resultest <- ifelse(any(x <= me_pred), which(x <= me_pred)[1], length(alpha) + 1)

      data.frame(dist = x, isOK=resultest==1 , anomaly = c("Normal >5%",paste0("Anormal <",alpha))[(resultest)])
    })))


  }else{

  

    #anomaly detection
    x_bar <- colMeans(refproj)
    S <- cov(refproj)         
    n=nrow(refproj)
    
    df1 <- ncol(refproj)    
    df2 <- n - df1
    F_crit <- qf(1 - alpha, df1 = df1, df2 = df2) # Critical F-value
    
   return(do.call(rbind,apply(proj2eval, 1, function(x) {
      mahal_dist <- sqrt(t(x - x_bar) %*% solve(S) %*% (x - x_bar))[1,1]
      critical_mahal <- sqrt(df1*(n+1)*(n-1)*F_crit/(n*(n-df1)))
      
      resultest <- ifelse(any(mahal_dist <= critical_mahal), which(mahal_dist <= critical_mahal)[1], length(alpha) + 1)

      data.frame(dist = mahal_dist, isOK=resultest==1 , anomaly = c("Normal >5%",paste0("Anormal <",alpha))[(resultest)])
    })))
  }

}



.pcaProj=function(newdata,pcaglmodel){
  PCAPROJ=NULL
  
  standev  <- pcaglmodel$PCA$call$ecart.type
  centre   <- pcaglmodel$PCA$call$centre
  projwmat <- pcaglmodel$PCA$svd$V
  varnames <- rownames(pcaglmodel$PCA$var$coord)
  
  

  # varnames[which(!varnames %in%colnames(newdata))]
  # colnames(newdata)[which(!colnames(newdata) %in%varnames)]


  if(is.vector(newdata)) {
    if(!all(varnames %in%names(newdata)  )) {
      stop("Missing variables in newdata")
    }
    newdata=newdata[varnames]
    newdata <- newdata - centre
    newdata <- newdata/standev
    PCAPROJ <- newdata %*% projwmat
  }else{
      if(!all(varnames %in%colnames(newdata)  )) {
      stop("Missing variables in newdata")
    }
    newdata <- newdata[, varnames]    
    newdata <- t(t(as.matrix(newdata)) - centre)
    newdata <- t(t(newdata)/standev)
    PCAPROJ <- crossprod(t(newdata), projwmat)
  }
  
  colnames(PCAPROJ) <- paste("Dim", c(1:ncol(PCAPROJ)), sep = ".")
  if(is.vector(newdata)) {
    rownames(PCAPROJ)="givenSample"
  }else{
    rownames(PCAPROJ) <-  rownames(newdata)
  }
  PCAPROJ
}

# newdata=newbotgdata[1,];pcaglmodel=smallpanetradio::MODELS$bothmodel;refdata=bothref

# newdata=newportdata[1,];pcaglmodel=smallpanetradio::MODELS$portmodel;refdata=smallpanetradio::refportdata

# newdata=newportdata;pcaglmodel=smallpanetradio::MODELS$portmodel;refdata=smallpanetradio::refportdata
# radiopred(newdata,pcaglmodel,refdata)

radiopred=function(newdata,pcaglmodel,refdata=NULL){
  require(Matrix)

  pbmsg="Variable names in newdata do not match the PCA model variables.\
     Use checkRadiomicdataframe to check. \
     Alternatively, yse the radiomicXlHelper function to load the data."
  

  if (is.vector(newdata)) {
      names(newdata)=sub("\\.","-",names(newdata))

    if(!all(rownames(pcaglmodel$PCA$var$coord)  %in%names(newdata)  )) {
      stop(pbmsg)
    }
  }else {
    colnames(newdata)=sub("\\.","-",colnames(newdata))
    if (!all(rownames(pcaglmodel$PCA$var$coord)  %in%colnames(newdata) )) {
        stop(pbmsg)
    }
  }
 

  PCAPROJ=smallpanetradio:::.pcaProj(newdata,pcaglmodel)
  
  selcomps=which((pcaglmodel$glmnet$beta)[-1,1]!=0)
  selpcaproj=PCAPROJ[,selcomps]

  if(is.vector(newdata)) {

    selpcaproj <- t(as.matrix(selpcaproj))
    rownames(selpcaproj)="givenSample"
    anomdf=data.frame(mahal_dist=NA,isOK=NA,anomaly="No ref")

  }else{
     anomdf=data.frame(dist=rep(NA,nrow(newdata)),isOK=rep(NA,nrow(newdata)),anomaly=rep("No ref",nrow(newdata)))
  }
  

  if(!is.null(refdata)){
    colnames(refdata)=sub("\\.","-",colnames(refdata))
    
    refproj=smallpanetradio:::.pcaProj(refdata,pcaglmodel)[,selcomps]
    
  
    if(is.vector(selpcaproj) & length(selcomps)>1){
      selpcaproj=t(as.matrix(selpcaproj))
      refproj=refproj[complete.cases(refproj),]
    }else if(length(selcomps)==1){
      refproj=refproj[complete.cases(refproj)]
      selpcaproj=as.matrix(selpcaproj)
    }else{
      refproj=refproj[complete.cases(refproj),]
    }

    anomdf= projAnomalyEval(selpcaproj,refproj)
    #  anomdf= projAnomalyEval(selpcaproj,FactoMineR::predict.PCA(pcaglmodel$PCA,refdata)$coord[,selcomps])
  }

    



  rez=data.frame(response=glmnet::predict.glmnet(pcaglmodel$glmnet,PCAPROJ,type="response")[,1],
  prob=predict(pcaglmodel$glmnet,PCAPROJ,type="response")[,1],
  anomdf)


  rownames(rez)=rownames(newdata)
  rez
}
# predict(pcaglmodel$glmnet,refPCAPROJ,type="response")


# thisvarnames=list(art=smallpanetradio::MODELS$artvars,port=smallpanetradio::MODELS$portvars)[[phasename]]




.oldplotrez=function(adata,amodel,label){
  
  print(names(adata))
  
  
  if(is.null(adata)){
    
    par(mar = c(0, 0, 0, 0))
    
    plot(x = 0:10, y = 0:10, ann = F,bty = "n",type = "n",
         xaxt = "n", yaxt = "n")
    
    text(x = 5,y = 5,paste0("no data for ",label))
  }else{

    
    refdf=data.frame( pred=predict(amodel,refdataset, na.action = na.pass, type = "prob"),
                        ref=refdataset$y)
    
    
    predprob=predict(amodel,adata, na.action = na.pass, type = "prob")
    predclass=predict(amodel,adata, na.action = na.pass, type = "raw")
    
    
    i=4
    
    print(ggviolin(refdf, x = "ref", y = "pred.met",ylab="Prediction of metastatic disease probability",
                   add = "boxplot",
                   title=c(art="Arterial phase",port="Portal phase")[label])+
            geom_hline(yintercept=predprob[i,"met"],
                       color=c("loc"="blue","met"="red")[levels(predclass)[predclass[i]]], size=1.5)+
            geom_hline(yintercept=0.5,color="grey",linetype="dashed",size=1.2)
    )
    
    
  }
  
}
.plotrez=function(newpred,refpred,label){
  
  if(is.null(newpred)){
    par(mar = c(0, 0, 0, 0))    
    plot(x = 0:10, y = 0:10, ann = F,bty = "n",type = "n",
         xaxt = "n", yaxt = "n")
    text(x = 5,y = 5,paste0("no data for ",label))
  }else{

    refpred$ref=c("loc","met")[1+refpred$y]


    print(ggviolin(refpred, x = "ref", y = "prob",ylab="Prediction of metastatic disease probability",
                   add = "boxplot",
                   title=c(art="Arterial phase",port="Portal phase")[label],
                            subtitle = paste("Anomaly evaluation:", newpred$anomaly[1])) +
            geom_hline(yintercept=newpred$prob[1],
                       color=c("loc"="blue","met"="red")[(newpred$prob[1] >0.5)+1], size=1.5)+
            geom_hline(yintercept=0.5,color="grey",linetype="dashed",size=1.2)
    )
    
    
  }
  
}





if(F){
 

  library(smallpanetradio)
  xlfile="/Users/remy.nicolle/Workspace/PDAC_aux/TNERadio/datasets/JulesGregoryTEST2.xlsx"
  # newdata=smallpanetradio:::radiomicXlHelper(xlfile,"port")
  # smallpanetradio:::.readxl(xlfile,"art")

  newartdata=smallpanetradio:::radiomicXlHelper(xlfile,"art")
  newportdata=smallpanetradio:::radiomicXlHelper(xlfile,"port")
  newbotgdata=data.frame(newartdata,newportdata)
  coms=intersect(rownames(smallpanetradio::refartdata),rownames(smallpanetradio::refartdata))
  bothref=data.frame(smallpanetradio::refartdata[coms,],smallpanetradio::refportdata[coms,])
  
  radiopred(newartdata,smallpanetradio::MODELS$artmodel,smallpanetradio::refartdata)
  radiopred(newartdata[1,],smallpanetradio::MODELS$artmodel,smallpanetradio::refartdata)
  
  radiopred(newportdata[1,],smallpanetradio::MODELS$portmodel,smallpanetradio::refportdata)
  radiopred(newportdata,smallpanetradio::MODELS$portmodel,smallpanetradio::refportdata)

  radiopred(newbotgdata[1,],smallpanetradio::MODELS$bothmodel,bothref)
  radiopred(newbotgdata,smallpanetradio::MODELS$bothmodel,bothref)


  radiopred(smallpanetradio::refartdata,smallpanetradio::MODELS$artmodel)

  radiopred(newportdata,smallpanetradio::MODELS$portmodel,smallpanetradio::refportdata)
  
  
  portvar=rownames( smallpanetradio::MODELS$portmodel$PCA$var$contrib)
  artvar=rownames( smallpanetradio::MODELS$artmodel$PCA$var$contrib)
  bothvar=rownames( smallpanetradio::MODELS$bothmodel$PCA$var$contrib)




  mean(bothvar %in% colnames(newbotgdata))
  colnames(newbotgdata)[which(!colnames(newbotgdata) %in% bothvar)]
  bothvar[which(!bothvar %in% colnames(newbotgdata))]
  bothvar[which(!bothvar %in% setdiff(colnames(newartdata),colnames(newportdata)))]

  bothvar[which(!bothvar %in% colnames(newbotgdata))]

  radiopred(newbotgdata[,],smallpanetradio::MODELS$bothmodel)

  
  pcaglmodel =smallpanetradio::MODELS$portmodel 


  radiopred(newdata,pcaglmodel,refdataset)

  # smallpanetradio:::.sheet2data( openxlsx::read.xlsx(xlfile,sheet=1),phase="art")

  # smallpanetradio::MODELS$allvars    smallpanetradio::MODELS$artvars    smallpanetradio::MODELS$bothvars   smallpanetradio::MODELS$portvars
  # smallpanetradio::MODELS$artmodel   smallpanetradio::MODELS$bothmodel  smallpanetradio::MODELS$portmodel  smallpanetradio::MODELS$unamevar


library(matrixStats)
library(openxlsx)

  refportpred=radiopred(smallpanetradio::refportdata,smallpanetradio::MODELS$portmodel)
  refportpred=data.frame(smallpanetradio::refclin[rownames(refportpred),],refportpred)

  refartpred=radiopred(smallpanetradio::refartdata,smallpanetradio::MODELS$artmodel)
        refartpred=data.frame(smallpanetradio::refclin[rownames(refartpred),],refartpred)
       

        coms=intersect(rownames(smallpanetradio::refartdata),rownames(smallpanetradio::refportdata))
        bothref=data.frame(smallpanetradio::refartdata[coms,],smallpanetradio::refportdata[coms,])
        
  
        refbothpred=radiopred(bothref,smallpanetradio::MODELS$bothmodel)
        refbothpred=data.frame(smallpanetradio::refclin[rownames(refbothpred),],refbothpred)

  allrefdf=data.frame(art=refartpred[coms,"prob"],port=refportpred[coms,"prob"],
  both=refbothpred[coms,"prob"],smallpanetradio::refclin[coms,])
  bjref=allrefdf[allrefdf$center=="bjn",]
  bjref$avgprob=rowMeans(bjref[,1:3])
  bjref$sdprob=rowSds(as.matrix(bjref[,1:3]))
  
  write.xlsx(bjref[order(bjref$avgprob),],file="~/Downloads/bjref.xlsx",
             rowNames=T,overwrite=T)


}


#saving code
if(F){
  # setwd("/Users/remy.nicolle/Workspace/PDAC_aux/TNERadio/radiopred")
  # library(openxlsx)
  # source("customglmforest.R")

  # bothmodel <- readRDS(paste0("models/bothmodel_pcaglmnet_0.55_20_1_TRUE.rds"))
  # artmodel <- readRDS(paste0("models/artmodel_pcaglmnet_0.6_20_1_TRUE.rds"))
  # portmodel <- readRDS(paste0("models/portmodel_pcaglmnet_0.6_20_1_TRUE.rds"))

  # clindf <- readRDS("models/cleanclindf.rds")
  # artdf <- readRDS("models/cleanartdf.rds")
  # portdf <- readRDS("models/cleanportdf.rds")



  # lvar=list(both=rownames(bothmodel$PCA$var$coord),art=rownames(artmodel$PCA$var$coord),port=rownames(portmodel$PCA$var$coord))
  # allvar=unique(unlist(lvar))
  # eachvar=unique(sub("^port|^art","",allvar))



  # varnames=readRDS( "models/varnames.rds")
  # refdataset=readRDS( "models/dataset.rds")

  # xlsheetfile="/Users/remy.nicolle/Workspace/PDAC_aux/TNERadio/datasets/Besancon_TNE_RADIOMIC.xlsx"
  # xltab=2#"2art"
  # xldat=read.xlsx(xlsheetfile,sheet=xltab,colNames=F)
  # x=sheet2data(xldat,phase="port")
  # radiopred(x,portmodel)

  # file="/Users/remy.nicolle/Workspace/PDAC_aux/TNERadio/datasets/JulesGregoryTEST2.xlsx";xltab=1; phasename="art"

  # #anomaly detection
      # x_bar <- colMeans(refpcaproj)
    # S <- cov(refpcaproj)         
    # n=nrow(refpcaproj)
    # alpha <- c(0.05,10^-(2:5)) 
    # df1 <- length(selcomps)    
    # df2 <- n - df1
    # F_crit <- qf(1 - alpha, df1 = df1, df2 = df2) # Critical F-value
    
    # anomdf=do.call(rbind,apply(selpcaproj, 1, function(x) {
    #   mahal_dist <- sqrt(t(x - x_bar) %*% solve(S) %*% (x - x_bar))[1,1]
    #   critical_mahal <- sqrt(df1*(n+1)*(n-1)*F_crit/(n*(n-df1)))
      
    #   resultest <- ifelse(any(mahal_dist <= critical_mahal), which(mahal_dist <= critical_mahal)[1], length(alpha) + 1)

    #   data.frame(dist = mahal_dist,  anomaly = c("Normal >5%",paste0("Anormal <",alpha))[(resultest)])
    # }))
  

  .doprediction=function(df){


  predf <- cbind(#df,
    selclin[rownames(df),],
    pred=predict(bothmodel, df, na.action = na.pass, type = "prob"),
    pred_class=predict(bothmodel, df, na.action = na.pass,type = "raw"),

    artpred=predict(artmodel, df, na.action = na.pass, type = "prob"),
    artpred_class=predict(artmodel, df, na.action = na.pass,type = "raw"),

    portpred=predict(portmodel, df, na.action = na.pass, type = "prob"),
    portpred_class=predict(portmodel, df, na.action = na.pass,type = "raw")
    
  )
  predf$dilatation = factor(c("nodilat","dilat")[predf$DilatationMD+1])
  predf$isValidation=predf$DateoffirstCT <= 41254
  predf$error=(predf$pred_class =="met") != predf$lab

  }
  .readxl=function(file,phasename,xltab=1){

    print(paste("Processing ", phasename))

    dataset=list(rawtab=NULL,
                msg="Unable to read data. Expecting an Excel\
                file with first line containing radiomic feature name",
                num=NULL)

    print(paste(file, "in", phasename))

    if(is.null(file)){
      dataset$msg="no file"
      return(dataset)
    }

    # thisvarnames=paste0(phasename,".", varnames)

    thisvarnames=list(art=smallpanetradio::MODELS$artvars,port=smallpanetradio::MODELS$portvars)[[phasename]]
    
    try({
      dataset$rawtab=smallpanetradio:::.sheet2data( openxlsx::read.xlsx(file,sheet=xltab),phase=phasename)

    })
    print(dim(dataset$rawtab))
    
    # dataset$rawtab=data.frame(dataset$rawtab)
    
    if(is.null(dataset$rawtab)){return(dataset)}

    missing=setdiff(thisvarnames,colnames(dataset$rawtab))
    
    if(length(missing)==0){
      print(paste("got no missing variables"))
    }else{
      print(paste("missing vars for ",phasename))
      print(missing)
    }
    
    if(!all(thisvarnames %in% colnames(dataset$rawtab))){
      print(paste("mssing all variable  in ",phasename))
      }else if(length(missing)>=(length(thisvarnames)-1)){
        dataset$msg=  paste("Missing all variables, file format likely inadequate")
        return(dataset)
      }else if(length(missing)>5 ){
        print(paste("mssing a llot of variable  in ",phasename))
        dataset$msg=  paste("Missing many variables")
        return(dataset)
      }else if(length(missing)>0){
        print(paste("mssing a few of variable  in ",phasename,"\n",
                    paste(missing, collapse=" , ")))
        dataset$msg=  paste("Missing following variables:",
                            paste(missing, collapse=" , "))
        return(dataset)
      }
      
    
      dataset$num=data.matrix(dataset$rawtab[,thisvarnames])
      
      print(paste("got a  dataset for ",phasename))
      
      print(paste("got a numeric dataset for ",phasename))
      
      if(sum(is.na(dataset$num[1,]) )>0 ){
        dataset$msg=paste("Variable names (column names) OK. \
                          Missing radiomic values")
        dataset$num=NULL
        return(dataset)
      }
      
      dataset$msg=paste("Data uploaded correctly")
      #colnames(dataset$num)=paste0(phase,".",colnames(dataset$num))
      
      print(names(dataset))
      return(dataset) 
    
    
  }

  .sheet2data=function(rawtab,coli=4:ncol(rawtab),predefineColnames=NULL,phase=""){

        df=data.frame(t(rawtab[,coli]))
        if(is.null(predefineColnames)){
          
          df=setNames(df,paste0(phase,rawtab[,2],"_",rawtab[,3]))
        }else{
          df=setNames(df,paste0(phase,predefineColnames))
        }
        
        df[] <- lapply(df, function(x) if(!any(is.na(as.numeric(as.character(x))))) as.numeric(as.character(x)) else x)
          rownames(df)=colnames(rawtab)[coli]
      df
    
  } 

}