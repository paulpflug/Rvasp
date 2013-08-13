#' Reads severaly calculations
#' 
#' \code{read.calculations} reads severaly calculations.
#' Needed folder structure
#' /FolderYouProvide/Parameter_Value/POSCAR
#' For two parameters:
#' /FolderYouProvide/Parameter_Value+Parameter2_Value/POSCAR
#' Will read CONTCAR / POSCAR / (total) Energy (from OSZICAR)
#' To read CHGCAR / Bands / Dos
#' provide CHGCAR-type filename / vasprun.xml / vasprun.xml
#' 
#' @param folders list/vector of folders containing a calculation each
#' @param name name
#' @param calculations provide an existing calculations object for updating or extending purpose
#' @param chgcarname if name is provided, calls \code{\link{read.chgcar}}
#' @param bandsxmlname if name is provided, calls \code{\link{read.bandsdata}}
#' @param dosxmlname if name is provided, calls \code{\link{read.dosdata}}
#' @param update if \code{TRUE} will override existing folders in calculations object provided by \code{calculations}
#' @export
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

#' Custom print for objects of class calculations
#' 
#' \code{print.calculations} custom print for objects of class calculations.
#' @param calculations objects of class calculations
#' @export
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

#' Custom print for objects of class calculation
#' 
#' \code{print.calculation} custom print for objects of class calculation.
#' @param calculation objects of class calculation
#' @export
print.calculation<-function(calculation,...){
  for (i in 1:length(calculation))
  {
    cat(paste(" Parameter:",names(calculation)[[i]],"\n"))
    print(calculation[[i]])
    cat(" **********\n")
  }
}

#' Will give the highest version of provided filename
#' 
#' \code{file.gethighestversion} will give the highest version of provided filename.
#' Searches for example for
#' POSCAR
#' POSCAR1
#' POSCAR2
#' POSCAR3
#' etc.
#' and will return the highest found filename
#' 
#' @param dir folder to search
#' @param filename filename
#' @export
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

#' Loads a calculations object in RData format
#' 
#' \code{load.calculations} loads a calculations object in RData format.
#' Will search current working directory for \code{name}.RData
#' 
#' @param file with calculations object to load
#' @param objectname will check if already loaded
#' @param update if object is already loaded determines, if object will be loaded again
#' @export
load.calculations<-function(file,update=F,objectname=NULL){
  if(length(grep(".RData",file))==0)
  {
    file <- paste0(file,".RData")
  }
  print(paste("loading",file))
  if (file.exists(file))
  {
    if (!is.null(objectname)){
      if(update|!exists(objectname)){
        load(file=file,globalenv())
      }
    }      
    else{
      load(file=file,globalenv())
    }
    return (T)
  }
  print(paste("didn't find",file))
  print("abort")
  return (F)
}


#' Fits equation of state to a given ea object
#' 
#' \code{ea.fitEOS} fits 3D equation of state to a given ea object.
#' 
#' @param fitdata ea object
#' @export
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

#' Predicts values based on a equation of state fit
#' 
#' \code{predict.EOS} predicts values based on a equation of state fit.
#' 
#' @param o object of EOS type
#' @param x values for which a y should be predicted
#' @export
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

#' Will give a vector containing lattice constants and energies
#' 
#' \code{calculation.getea} will give a vector containing lattice constants and energies
#' based on calculation object.
#' Vector will have class ea.
#' 
#' @param calculation object of type calculation
#' @export
calculation.getea<-function(calculation)
{
  data <- data.frame(do.call(rbind,lapply(calculation,FUN=function(x)
    if(!is.null(x$poscar$a) & !is.null(x$energy))
    cbind(x$poscar$a,x$energy))))
  names(data)<-c("a","E")
  class(data)<-"ea"
  return(data)
}

#' Will plot a e over a curve
#' 
#' \code{plot.calculation.ea} will plot a e over a curve.
#' 
#' @param calculation object of type calculation
#' @param energyfactor will be multiplied with energy
#' @param energyshift will be substracted from energy
#' @param fit if \code{TRUE} will fit with \code{\link{ea.fitEOS}} plot and return result
#' @param ... further plotting parameters
#' @export
plot.calculation.ea<-function(calculation,energyfactor=1,energyshift=0,fit=F,type="p",...)
{
  data <- calculation.getea(calculation)
  data$E <- data$E * energyfactor-energyshift
  rdata <- data.frame(a=data$a, E=data$E)
  plot(rdata, type=type, ...)
  if(fit){
    o <- ea.fitEOS(data)
    plot.EOS.add(o,...)
    return(o)
  }
}

#' Adds a curve to an existing e over a curve
#' 
#' \code{plot.calculation.ea.addpoints} adds a curve to an existing e over a curve.
#' 
#' @param calculation object of type calculation
#' @param afactor will be multiplied with lattice constant
#' @param energyfactor will be multiplied with energy
#' @param energyshift will be substracted from energy
#' @param fit if \code{TRUE} will fit with \code{\link{ea.fitEOS}} plot and return result
#' @param ... further plotting parameters
#' @export
plot.calculation.ea.addpoints<-function(calculation,afactor=1,energyfactor=1,energyshift=0,fit=F,...){
  data <- calculation.getea(calculation)
  data$a <- data$a*afactor
  data$E <- data$E * energyfactor-energyshift
  rdata <- data.frame(a=data$a, E=data$E)
  points(rdata,...)
  if(fit){
    o <- ea.fitEOS(data)
    plot.EOS.add(o,...)
    return(o)
  }
  
}

#' Adds a fit to an existing e over a curve
#' 
#' \code{plot.calculation.ea.addfit} adds a fit to an existing e over a curve.
#' 
#' @param calculation object of type calculation
#' @param afactor will be multiplied with lattice constant
#' @param energyfactor will be multiplied with energy
#' @param energyshift will be substracted from energy
#' @param ... further plotting parameters
#' @export
plot.calculation.ea.addfit<-function(calculation,afactor=1,energyfactor=1,energyshift=0,...){
  data <- calculation.getea(calculation)
  data$a <- data$a*afactor
  data$E <- data$E * energyfactor-energyshift
  o <- ea.fitEOS(data)
  plot.EOS.add(o,...)
  return(o)
}

#' Adds a EOS fit to an existing plot
#' 
#' \code{plot.EOS.add} adds a EOS fit to an existing plot.
#' fits with \code{\link{ea.fitEOS}}.
#' 
#' @param o object of EOS class
#' @param ... further plotting parameters
#' @export
plot.EOS.add<-function(o,...){
  x <- seq(o$arange[[1]],o$arange[[2]],length.out=101)
  lines(predict.EOS(o,x),...)
}