
read.calculations <- function(folders,name="calculations",calculations=list(),chgcarname="",bandsxmlname="",dosxmlname="",update=F){
  class(calculations) <- "calculations"
  if (is.null(calculations[["name"]])) calculations[["name"]]<-name
  calculations[["folders"]]<-unique(c(calculations[["folders"]],folders))
  if (is.null(calculations[["names"]]))  calculations[["names"]]<-character()
  for (folder in folders)
  {    
    foldername <- folder#tail(strsplit(folder,"/")[[1]],1)
    if (!foldername %in% calculations[["names"]] )
    {
      calculations[["names"]]<- c(calculations[["names"]],foldername)
    }
    dirs <- list.dirs(folder)
    print(paste("processing",folder,"( subfolders:",length(dirs)-1,")"))
    if (length(dirs)>1)
    {
      if (is.null(calculations[[foldername]]))
        calculation<- list()
      else
        calculation <- calculations[[foldername]]        
      class(calculation) <- "calculation"
      for (dir in dirs[-1])
      {     
        
        ### Parameters
        splitteddir <- strsplit(dir,"/")[[1]]
        parasstr <- splitteddir[length(splitteddir)]
        parassplitted <-strsplit(parasstr,"+",fixed=T)[[1]]
        if (!parasstr %in% names(calculation) | update)
        {
          print(paste("processing",parasstr))
          if (is.null(calculation[[parasstr]]))
            calcdata<- list()
          else
            calcdata <- calculation[[parasstr]]       
          class(calcdata) <- "calculationdata"
          calcdata[["paralist"]] <- sapply(parassplitted,FUN=function(x){
            para <- strsplit(x,"_")[[1]]
            return(cbind(para[[1]],para[[2]]))
                                                                         })
          for (ps in parassplitted)
          {
            para <- strsplit(ps,"_")[[1]]
            
            if (length(para)==2)
              calcdata[[para[1]]] <- para[2]
          }
          ### Energie
          oszicar <- file.gethighestversion(dir,"OSZICAR")
          if(oszicar$success)
          {
            oszidata <- read.table(oszicar$file,sep="\n",stringsAsFactors=F,strip.white=T)
            calcdata[["energy"]]<- as.numeric(strsplit(oszidata[nrow(oszidata),]," +")[[1]][[5]])
          }
          ### POSCAR
          poscar <- file.gethighestversion(dir,"POSCAR")
          if(poscar$success)
          {
            calcdata[["poscar"]]<- read.poscar(poscar$file)
          }
          ### CONTCAR
          contcar <- file.gethighestversion(dir,"CONTCAR")
          if(contcar$success)
          {
            calcdata[["contcar"]]<- read.poscar(contcar$file)
          }
          ### CHGCAR
          if (chgcarname!="")
          {
            chgcar <- file.gethighestversion(dir,chgcarname)
            if (chgcar$success)
            {
              calcdata[["chgcar"]] <- read.chgcar(chgcar$file)
            }
          }
          ### BANDS
          if (bandsxmlname!="")
          {
            bandsfile <- file.gethighestversion(dir,bandsxmlname)
            if(bandsfile$success)
            {
              calcdata[["banddata"]] <- read.bandsdata(bandsfile$file)
            }
            
          }
          ### DOS
          if (dosxmlname!="")
          {
            dosfile <- file.gethighestversion(dir,dosxmlname)
            if(dosfile$success)
            {
              calcdata[["dosdata"]] <- read.dosdata(dosfile$file)
            }
            
          }
        
          calculation[[parasstr]]<- calcdata
        }
        else
          print(paste("skipping",parasstr))
      }
      calculations[[foldername]]<-calculation
      gc()
    }
    
  }
  return (calculations)
}

print.calculations<-function(calculations,...){

  cat(paste("###",calculations$name,"###\n"))
  cat(paste("Calculations:",length(calculations$names)))
  cat("\n###########\n")
  for (name in calculations$names)
  {
    cat(paste(tail(strsplit(name,"/")[[1]],1),"\n"))
    print(calculations[[name]])
    cat("###########\n")
  }
}

print.calculation<-function(calculation,...){
  for (i in 1:length(calculation))
  {
    cat(paste(" Parameter:",names(calculation)[[i]],"\n"))
    print(calculation[[i]])
    cat(" **********\n")
  }
}

file.gethighestversion <- function(dir,filename)
{
  file <- list()
  file$success <- F
  files <- list.files(dir)
  result <- grep(filename,files,ignore.case=T,value=T)
  if (length(result)>0)
  {
    matches <- regexec("[0-9]+",result)
    number<-as.numeric(regmatches(result,matches))
    number[is.na(number)]<-0
    file$file<-paste(dir,result[which.max(number)],sep="/")
    file$success <- T
  }
  return(file)
}
load.calculations<-function(name,update=F)
{
  filename <- paste0(name,".RData")
  print(paste("loading",filename))
  if (file.exists(filename))
  {
    if (!exists(name)|update)
      load(file=filename,globalenv())
    return (T)
  }
  else
  {
    return (F)
  }
}
ea.fitEOS <- function(fitdata)
{
  getstartEOS<-function(data)
  {
    data$a <- data$a^3
    lfit <- lm(E~a+I(a^2),data)
    a<- lfit$coefficients[[3]]
    b<-lfit$coefficients[[2]]
    c<- lfit$coefficients[[1]]
    result <- list()
    v0 <- -b/2/a
    e0 <- a*(v0)^2+b*v0+c
    b0 <- 2*a*v0
    bP <- 5
    return(c(v0,e0,b0,bP))
  }
  murnaghan<-function(paras,vol)
  {
    v0<-paras[[1]]
    e0 <- paras[[2]]
    b0 <- paras[[3]]
    bP <- paras[[4]]
    E<-e0+b0*vol/bP*(((v0/vol)^bP)/(bP-1)+1)-v0*b0/(bP-1)
    return (E)
  }
  objectivem <- function(paras,data)
  {
    err = sum((data$E-murnaghan(paras,data$a^3))^2)
    return (err)
  }
  starteos <- getstartEOS(fitdata)
  
  obj <- optim(starteos,objectivem,data=fitdata,method="BFGS")
  obj$arange<-range(fitdata$a)
  class(obj)<-"EOS"
  return(obj)
}
predict.EOS<-function(o,x)
{
  murnaghan<-function(paras,vol)
  {
    v0<-paras[[1]]
    e0 <- paras[[2]]
    b0 <- paras[[3]]
    bP <- paras[[4]]
    E<-e0+b0*vol/bP*(((v0/vol)^bP)/(bP-1)+1)-v0*b0/(bP-1)
    return (E)
  }
  data <- data.frame(x,murnaghan(o$par,x^3))
  names(data) <- c("a","E")
  return(data)
}
calculation.getea<-function(calculation)
{
  data <- data.frame(do.call(rbind,lapply(calculation,FUN=function(x)
    if(!is.null(x$poscar$a) & !is.null(x$energy))
    cbind(x$poscar$a,x$energy))))
  names(data)<-c("a","E")
  return(data)
}
plot.calculation.ea<-function(calculation,energyfactor=1,energyshift=0,fit=F,type="p",...)
{
  data <- calculation.getea(calculation)
  data$E <- data$E * energyfactor-energyshift
  plot(data,type=type,...)
  if(fit){
    o <- ea.fitEOS(data)
    plot.EOS.add(o,...)
    return(o)
  }
}

plot.calculation.ea.addpoints<-function(calculation,afactor=1,energyfactor=1,energyshift=0,fit=F,...)
{
  data <- calculation.getea(calculation)
  data$a <- data$a*afactor
  data$E <- data$E * energyfactor-energyshift
  points(data,...)
  if(fit){
    o <- ea.fitEOS(data)
    plot.EOS.add(o,...)
    return(o)
  }
  
}
plot.calculation.ea.addfit<-function(calculation,afactor=1,energyfactor=1,energyshift=0,...)
{
  data <- calculation.getea(calculation)
  data$a <- data$a*afactor
  data$E <- data$E * energyfactor-energyshift
  o <- ea.fitEOS(data)
  plot.EOS.add(o,...)
  return(o)
}
plot.EOS.add<-function(o,...)
{
  x <- seq(o$arange[[1]],o$arange[[2]],length.out=101)
  lines(predict.EOS(o,x),...)
}