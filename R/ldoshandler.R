#' Adds a gaussian smearing to ldos vector
#' 
#' \code{ldosvector.calcsmearing} Adds a gaussian smearing
#' 
#' @param energy vector containing the energy
#' @param dos vector containing the dos
#' @param sigma smearing parameter
#' @export
ldosvector.calcsmearing<-function(energy,dos,sigma=0.1){
  matrix <- sapply(energy,dnorm,energy,sigma) 
  return (rowSums(sweep(matrix,2,colSums(matrix)/dos,"/")))
}

#' Creates a ldos image.
#' Based on base package.
#' 
#' \code{plot.calculation.ldos} Creates a ldos image
#' 
#' @param calculation object of class calculation
#' @param positions object of class positions
#' @param smearing (optional) adds gaussian smearing
#' @param efermi (optional) alignment to the Fermi level
#' @param legend (optional) should an automatic legend should be created?
#' @export
plot.calculation.ldos<- function(calculation,
                                 positions,
                                 col=rainbow(length(positions)),
                                 typ="l",
                                 efermi=0,
                                 smearing=NULL,
                                 xlab=expression("Energy"~epsilon~"(eV)"),
                                 ylab=expression(rho*"(x,"*epsilon*")"~("\u00C5"^-3)),
                                 legend=F,
                                 lty=1,
                                 ...){
  indices <- calculation.getindicesfrompositions(calculation,positions)
  if(!is.null(indices))
  {
  dataset <- calculation.getdatafromindices(calculation,indices,efermi=efermi,smearing=smearing)  
  rng <- do.call(rbind,dataset)
  rng<- apply(rng,2,range)
  plot(rng,xlab=xlab,ylab=ylab,typ="n",...)
  if(typ=="l")
  {
    plot.calculation.addldos(calculation,positions,efermi=efermi,smearing=smearing,dataset=dataset,col=col,lty=lty,...)
  }
  if(legend)
  {
    plot.positions.addlegend(positions,lty=lty,col=col)
  }
  }
}

#' searches for the indices of the value in the CHGCAR for a list of given positions
#' 
#' \code{calculation.getindicesfrompositions} searches for the index of the value in the CHGCAR by a given position
#' 
#' @param calculation object of class calculation
#' @param positions object of class positions
#' @export
calculation.getindicesfrompositions<-function(calculation,positions){
  if(!is.null(calculation[[1]]$chgcar))
  {
  xyz <- calculation[[1]]$chgcar$data[,1:3]
  indices <- sapply(positions,FUN=function(position){  
    centeredxyz <- sweep(as.matrix(xyz),2,position[1:3])
    index <- which.min(rowSums(abs(centeredxyz)))    
    return(index)
  })
  for (i in 1:length(positions))
  print(paste("found position",names(positions)[[i]],"offset:",paste(signif(xyz[indices[[i]],]-positions[[i]],2),collapse=" ")))
  return(indices)
  }
  return(NULL)
}

#' gives the CHGCAR data by a list of given dataindices
#' 
#' \code{calculation.getdatafromindices} gives the CHGCAR data by a given dataindex
#' 
#' @param calculation object of class calculation
#' @param indices dataindices based on positions
#' @param efermi (optional) alignment to the Fermi level
#' @param smearing (optional) adds gaussian smearing
#' @export
calculation.getdatafromindices <- function(calculation,indices,efermi=0,smearing=NULL){
  dataset <- lapply(indices,function(index){
  data <- lapply(calculation,function(para)
  {
    intervall <- as.numeric(strsplit(para$EINT," +")[[1]])
    eint <- mean(intervall)
    np <- ifelse(para$ldosdata$npoints>0,para$ldosdata$npoints,1)
    return (c(eint-efermi,para$chgcar$data[index,4]/abs(intervall[1]-intervall[2])))
  #  return (c(eint-efermi,para$chgcar$data[index,4]/np))
  })
  data <- do.call(rbind,data)
  data <- data[order(data[,1]),]
  if(!is.null(smearing))
  {
    data <- cbind(data[,1],ldosvector.calcsmearing(data[,1],data[,2],smearing))
  }
  return(data)
  })
  return(dataset)
}

#' adds a ldos to an existing plot. either positions or dataset must be given.
#' 
#' \code{plot.calculation.addldos} adds a ldos to an existing plot
#' 
#' @param calculation object of class calculation
#' @param positions object of class positions
#' @param dataset list of ldos data objects
#' @param efermi (optional) alignment to the Fermi level
#' @param smearing (optional) adds gaussian smearing
#' @export
plot.calculation.addldos<-function(calculation
                                   ,positions=NULL
                                   ,efermi=0
                                   ,smearing=NULL
                                   ,dataset=NULL
                                   ,col="black"
                                   ,lty=1
                                   ,...){
  stopifnot(!is.null(dataset)|!is.null(positions))
  if(is.null(dataset))
  {    
    indices <- calculation.getindicesfrompositions(calculation,positions)
    dataset <- calculation.getdatafromindices(calculation,indices,efermi=efermi,smearing=smearing)
  }
  lty <- rep(lty,length.out=length(dataset))
  for (i in 1:length(dataset))
    lines(dataset[[i]],col=col[[i]],lty=lty[[i]],...)
}

#' adds a positions legend to an existing plot.
#' 
#' \code{plot.positions.addlegend} adds a positions legend to an existing plot
#' 
#' @param positions object of class positions
#' @param ... further ploting parameters
#' @export
plot.positions.addlegend<-function(positions,lty=1,col=rainbow(length(positions)),position="topright",...){
  s<-c("")
  p <- do.call(rbind,positions)
  if(length(unique(p[,4]))>1)
   s<-paste("layer",p[,4])
  if(length(unique(p[,5]))>1)
    s<-paste0(s," distance ",p[,5],"Å")
  legend(x=position,legend=s,col=col,lty=lty,...)
}
# load.calculationdata(name="ldosdata1")
#calculation <- ldosdata[[4]]
# positions <- getpositions.poscar(calculation[[1]]$poscar,layers=14,layer=14,zdist=c(1,2))
# plotldos.calculation(calculation,positions)
#addlegend.calculation(calculation,positions)
# position <- getposition.poscar(calculation[[1]]$poscar,layers=14,layer=13)
# plotldos.calculation(calculation,position,ylim=c(0,0.0015))