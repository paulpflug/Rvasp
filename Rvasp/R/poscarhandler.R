#' Creates a new object of class poscar
#' 
#' \code{create.poscar} creates a new object of class poscar.
#' 
#' @export
create.poscar<-function(){
  poscar <- list()
  poscar$firstline <- "new Poscar"
  poscar$a <- 1
  poscar$basis <- cbind(c(1,0,0),c(0,1,0),c(0,0,1))
  poscar$eightline <-"Selective dynamics"
  poscar$ninthline <-"Direct"
  poscar$atoms <- data.frame(0.5,0.5,0.5,"T","T","T","H")
  comps <- c("x","y","z")
  colnames(poscar$atoms)<-c(paste("r",comps,sep=""),paste("move",comps,sep=""),"type")
  class(poscar)<-"poscar"
  return(poscar)
}

#' Reads file in POSCAR format
#' 
#' \code{read.poscar} reads file in POSCAR format.
#' returns object of class poscar
#' 
#' @param file filename of poscar
#' @export
read.poscar<-function(file){
  data <- read.table(file,strip.white=T,sep="\n",stringsAsFactors=F)
  poscar <- parse.poscar(data)
  return(poscar)
}

#' Parses a dataframe into object of class poscar
#' 
#' \code{parse.poscar} parses a dataframe into object of class poscar.
#' The dataframe only consists of one column whereas each row contains a row of the poscar.
#' 
#' @param stringdataframe Dataframe of strings
#' @export
parse.poscar <- function(stringdataframe){
  poscar <- list()
  poscar$firstline <- stringdataframe[1,1]
  poscar$a <- as.numeric(stringdataframe[2,1])
  poscar$basis <- apply(do.call(rbind,strsplit(stringdataframe[3:5,1]," +")),2,as.numeric)
  atoms <- strsplit(stringdataframe[6,1]," +")[[1]]
  atomnumbers <- as.numeric(strsplit(stringdataframe[7,1]," +")[[1]])
  startatompos <- ifelse(length(strsplit(stringdataframe[9,1]," +")[[1]])>2,9,10)
  poscar$eightline <- stringdataframe[8,1]
  if (startatompos==10)
  {
    poscar$ninthline <-stringdataframe[9,1]
    poscar$selectivedynamics <- T
  }
  else
  {
    poscar$selectivedynamics <- F
  }
  block <-do.call(rbind,strsplit(stringdataframe[startatompos:(startatompos+sum(atomnumbers)-1),1]," +"))
  if (nrow(block)>1)
  {
    poscar$atoms <- data.frame(t(apply(block[,1:3],1,as.numeric)),stringsAsFactors=F)
    if(ncol(block)>3)
      poscar$atoms <- cbind(poscar$atoms ,block[,4:6])
    poscar$atoms <- cbind(poscar$atoms ,rep(atoms,atomnumbers))
  }
  else
  {
    poscar$atoms <- data.frame(t(as.numeric(block[1:3])))
    if(ncol(block)>3)
      poscar$atoms <- cbind(poscar$atoms ,t(block[,4:6]))
    poscar$atoms <- cbind(poscar$atoms ,atoms[[1]])
  }
  poscar$length <- (startatompos+sum(atomnumbers)-1)
  comps <- c("x","y","z")
  colnames <- paste("r",comps,sep="")
  if(ncol(poscar$atoms)>4)
    colnames <- c(colnames,paste("move",comps,sep=""))
  colnames <- c(colnames,"type")
  colnames(poscar$atoms)<-colnames
  class(poscar)<-"poscar"
  return(poscar)
}

#' Turnes object of class poscar into POSCAR format
#' 
#' \code{format.poscar} turnes object of class poscar into POSCAR format.
#' 
#' @param poscar object of class poscar
#' @export
format.poscar<-function(poscar,...){
  poscar$atoms <- format(poscar$atoms,scientific=F,nsmall=10)
  poscar$a <- format(poscar$a,scientific=F)
  poscar$basis <- format(poscar$basis,scientific=F,nsmall=10)
  collapsedbasis <- apply(poscar$basis,1,function(x)paste(x,collapse=' '))
  atoms <- c()
  atomnumbers <-c()
  lastatom <- ""
  index <- 0
  selectivdynamics <- poscar$selectivedynamics
  if (!is.null(selectivdynamics) &&selectivdynamics&is.null(poscar$ninthline))
  {
    poscar$ninthline <- poscar$eightline
    poscar$eightline <- "Selective dynamics"
    
    poscar$atoms <- cbind(poscar$atoms[1:3],matrix("T",ncol=3,nrow=nrow(poscar$atoms)),poscar$atoms[4])
  }
  atomcolumn <- 7
  if (is.null(poscar$ninthline))
    atomcolumn <- 4
  poscar$atoms[atomcolumn] <- lapply(poscar$atoms[atomcolumn], as.character)
  for (i in 1:nrow(poscar$atoms))
  {
    atom <-poscar$atoms[i,atomcolumn]
    if(poscar$atoms[i,atomcolumn]!=lastatom)
    {
      index <- index+1
      lastatom<-atom
      atoms <- c(atoms,atom)
      atomnumbers <- c(atomnumbers,0)
    }
    atomnumbers[index]<-atomnumbers[index]+1
  }
  collapsedatoms<- paste(atoms,collapse=' ')
  collapsedatomnumbers<- paste(atomnumbers,collapse=' ')
  collapsedatomstrings <- apply(poscar$atoms[,1:(atomcolumn-1)],1,function(x)paste(x,collapse=' '))
  data <- c(poscar$firstline,poscar$a,collapsedbasis,collapsedatoms,collapsedatomnumbers,poscar$eightline)
  if (!is.null(poscar$ninthline))
    data <- c(data,poscar$ninthline)
  data<-c(data,collapsedatomstrings)
  return(data)
}

#' Writes object of class poscar into file of POSCAR format
#' 
#' \code{write.poscar} writes object of class poscar into file of POSCAR format.
#' 
#' @param poscar object of class poscar
#' @param file name of file
#' @export
write.poscar<-function(poscar,file){  
    cat(paste(format.poscar(poscar),collapse="\n"),file=file)
    return(T)
}

#' Custom print for object of class poscar
#' 
#' \code{print.poscar} custom print for object of class poscar.
#' 
#' @param poscar object of class poscar
#' @export
print.poscar <- function(poscar,...){
  cat("\n###POSCAR###\n")
  for (name in names(poscar))
  {
    if (length(poscar[[name]])>1)
    {
      cat(paste0("   ",name," length:",length(poscar[[name]])," (",class(poscar[[name]]),")","\n"))    
    }
    else
    {
      cat(paste0("   ",name," ",poscar[[name]],"\n"))
    }
  }
  cat("###########\n")
  cat("printraw.poscar(poscar) for copy+paste poscar")
}

#' Raw print for object of class poscar
#' 
#' \code{printraw.poscar} prints object of class poscar in ASCII format.
#' 
#' @param poscar object of class poscar
#' @export
printraw.poscar <- function(poscar){
  cat("###POSCAR###\n")
  cat(paste(format.poscar(poscar),collapse="\n"))
  cat("\n###END of POSCAR###")
}

#' Custom plot for object of class poscar
#' 
#' \code{plot.poscar} custom plot for object of class poscar.
#' Creates only the plotbox for 2d projection plot.
#' 
#' @param poscar object of class poscar
#' @param direction of projection
#' @param basis if activated, draws basis
#' @param unitcell if activated, draws basis
#' @param fullcell if activated, draws full cell, rather than only the volume where atoms are
#' @param ... further plotting parameters
#' @export
plot.poscar <- function(poscar,direction=3,xlab="x",ylab="y",basis=F,unitcell=F,fullcell=F,...){
  dir <- (1:3)[-direction]
  atoms <- poscar.getbasisconvertedatoms(poscar)
  rng<-apply(atoms[,dir],MARGIN=2,FUN=range)
  if (fullcell)
  {
    rng <- (poscar$basis*poscar$a)[dir,dir]
  }
  plot(rng[,1],rng[,2],asp=1,typ="n",axes=F,xlab="",ylab="",xaxs="i",yaxs="i",...)
  if(basis)
    plot.poscar.addbasis(poscar,direction=direction,...)
  if(unitcell)
    plot.poscar.addunitcell(poscar,direction=direction,...)
}

#' Adds basis to existing plot
#' 
#' \code{plot.poscar.addbasis} adds basis to existing plot.
#' 
#' @param poscar object of class poscar
#' @param direction of projection
#' @param xoffin offset in x direction in inch
#' @param yoffin offset in y direction in inch
#' @param basisnames names of basis vectors
#' @param fullcell if activated, draws full cell, rather than only the volume where atoms are
#' @param arrowlength length of basis arrows in user coordinates
#' @param arrowsize size of arrow tip
#' @param ... further plotting parameters
#' @export
plot.poscar.addbasis<-function(poscar
                              ,xoffin=0.1
                              ,yoffin=0.1
                              ,direction=3
                              ,basisnames=c("a","b","c")
                              ,arrowlength=0
                              ,arrowsize=0.1
                              ,fullcell=F
                              ,...){
  x1 <- grconvertX(xoffin,from="inches")
  y1 <- grconvertY(yoffin,from="inches")
  dir <- (1:3)[-direction]
  basis <- poscar$basis*poscar$a  
  if (fullcell)
  {
    rng <- basis[dir,dir]
  }
  fac1 <- 1
  fac2 <- 1
  if (arrowlength>0)
  {
  l1 <- sqrt(sum(basis[dir[1],dir]^2))
  l2 <- sqrt(sum(basis[dir[2],dir]^2))
  fac1 <- arrowlength/l1
  fac2 <- arrowlength/l2
  }
  target1 <- c(x1,y1)+basis[dir[1],dir]*fac1
  target2 <- c(x1,y1)+basis[dir[2],dir]*fac2
  arrows(x1,y1,target1[[1]],target1[[2]],length=arrowsize,xpd=T)
  arrows(x1,y1,target2[[1]],target2[[2]],length=arrowsize,xpd=T)
  text(target1[1],target1[2],basisnames[dir[1]],xpd=T,pos=4,...)
  text(target2[1],target2[2],basisnames[dir[2]],xpd=T,pos=3,...)
}

#' Adds unitcell to existing plot
#' 
#' \code{plot.poscar.addunitcell} adds unitcell to existing plot.
#' 
#' @param poscar object of class poscar
#' @param direction of projection
#' @param ... further plotting parameters
#' @export
plot.poscar.addunitcell<-function(poscar,direction=3,lty=2,...){
  dir <- (1:3)[-direction]
  basis <- poscar$basis[dir,dir]*poscar$a
  x<- c(0,basis[1,1],basis[1,1]+basis[2,1],basis[2,1])
  y <- c(0,basis[1,2],basis[1,2]+basis[2,2],basis[2,2])
  polygon(x=x,y=y,lty=lty,xpd=T,...)
}

#' Adds atoms to existing plot
#' 
#' \code{plot.atoms.add} adds atoms to existing plot.
#' 
#' @param atoms dataframe of atoms
#' @param direction of projection
#' @param basis basis if atoms are in direct coordinates
#' @param ... further plotting parameters
#' @export
plot.atoms.add<-function(atoms,basis=NULL,direction=3,col="white",cex=3,lwd=1,lty=1,...){
  dir <- (1:3)[-direction]
  if (!is.null(basis))
    atoms <- atoms.convertbasis(atoms,basis,forth=F)
  atoms <- atoms[order(atoms[,direction]),]  
  points(atoms[,dir],col="black",bg=col,pch=21,cex=cex,lwd=lwd,xpd=T,lty=1,...)
}

#' Adds layers of atoms to existing plot
#' 
#' \code{plot.poscar.addlayers} adds layers of atoms to existing plot.
#' 
#' @param poscar object of class poscar
#' @param direction of projection
#' @param layer vector of indices of layers to plot
#' @param layers total layer count in poscar
#' @param ... further plotting parameters
#' @export
plot.poscar.addlayers<-function(poscar
                                ,layer
                                ,layers
                                ,direction=3
                                ,color=rainbow(length(layer))
                                ,size=rep(1,length(layer))
                                ,lwd=rep(1,length(layer))
                                ,...){
  dir <- (1:3)[-direction]  
  layerindices <- poscar.getatomlayerindices(poscar,layers)
  atoms <- poscar.getbasisconvertedatoms(poscar)
  color <- rep(color,length.out=length(layer))
  size <- rep(size,length.out=length(layer))
  lwd <- rep(lwd,length.out=length(layer))
  
  data <- cbind(atoms[layerindices%in%layer,1:3],layerindices[layerindices%in%layer])
  data[,4] <- sapply(data[,4],FUN=function(x)which(x==layer))
  data <- cbind(data,color[data[,4]],size[data[,4]],lwd[data[,4]])
  data <- data[order(data[,direction]),]
  points(data[,dir],col="black",bg=as.character(data[,5]),pch=21,cex=data[,6],lwd=data[,7],xpd=T,...)
}

#' Adds distance between layers to existing plot
#' 
#' \code{plot.poscar.addlayers} adds distance between layers of atoms to existing plot.
#' 
#' @param poscar object of class poscar
#' @param direction of projection
#' @param layer vector of indices of layers to plot
#' @param layers total layer count in poscar
#' @param color color of lines
#' @param ... further plotting parameters
#' @export
plot.poscar.addlayerdistance<-function(poscar,layer,layers,color="black",direction=1,...){
  stopifnot(length(layer)>1)
  dir <- (1:3)[-direction]
  atoms <- poscar.getbasisconvertedatoms(poscar) 
  rawrng<-apply(atoms[,dir],MARGIN=2,FUN=range)
  layerindices <- poscar.getatomlayerindices(poscar,layers)
  plotrng <- matrix(par("usr"),nrow=2)
  off <- (plotrng[2,]-plotrng[1,])*0.05  
  offset <- rbind(-off,+off)
  rng <- rawrng+offset
  z <- sapply(layer,FUN=function(x)mean(atoms[layerindices==x,3]))
  color <- rep(color,length.out=length(layer))
  for (i in 1:length(layer))
{
  lines(rawrng[2,1]+rbind(-off[[1]]*2,off[[1]]*0.5),c(z[[i]],z[[i]]),lty=2,col=color[[i]],xpd=T,...)
  }
  for (i in 1:(length(z)-1))
  {
    dist<- round(abs(z[i]-z[i+1]),2)
    legend((rawrng+offset*0.15)[2,1],min(c(z[i+1],z[i]))+dist/2,legend=paste(dist,"\u00C5"),xjust=0,yjust=0.6,x.intersp=-0.4,y.intersp=0.3,bty="n",bg=NULL,xpd=T,...)
  }
}

#' Adds distance between atoms to existing plot
#' 
#' \code{plot.atoms.adddistance} adds distance between atoms to existing plot.
#' 
#' @param atoms dataframe of atoms, will only use first two atoms
#' @param direction of projection
#' @param basis basis if atoms are in direct coordinates
#' @param length arrow tip length
#' @param lwd arrow line width
#' @param col arrow col
#' @param selectedalpha transparency level of highlighted point over selected atom
#' @param selectedsize if larger than 0 will draw white transparent point over selected atom for highlighting purpose
#' @param ... further plotting parameters
#' @export
plot.atoms.adddistance<-function(atoms
                                 ,basis=NULL
                                 ,direction=3
                                 ,length=0.1
                                 ,lwd=1
                                 ,col="black"
                                 ,selectedalpha=150
                                 ,selectedsize=1.5
                                 ,...){
  dir <- (1:3)[-direction]
  if (!is.null(basis))
    atoms <- atoms.convertbasis(atoms,basis,forth=F)
  atoms <- atoms[c(1,2),dir]
  if(selectedsize>0)
    points(atoms,pch=21,bg=makeTransparent("white",selectedalpha),col=makeTransparent("white",selectedalpha),cex=selectedsize)
  arrows(atoms[1,1],atoms[1,2],atoms[2,1],atoms[2,2],length=length,code=3,lwd=lwd,col=col)
  rdist<-c((atoms[2,1]-atoms[1,1]),(atoms[2,2]-atoms[1,2]))
  dist<- round(sqrt(sum(rdist^2)),2)
  r <-atoms[1,]+rdist/2 
  legend(r[1],r[2],legend=paste(dist,"Å"),bg="white",xjust=0.5,yjust=0.5,x.intersp=-0.4,y.intersp=0.3 ,...)
}

#' Adds numbers on atoms to existing plot
#' 
#' \code{plot.poscar.addnumbers} adds numbers on atoms, based on layers, to existing plot.
#' 
#' @param poscar object of class poscar
#' @param direction of projection
#' @param layer vector of indices of layers to plot
#' @param layers total layer count in poscar
#' @param absolutenumber determines if index of atom in poscar(absolute) or index of atom in layer is printed
#' @param ... further plotting parameters
#' @export
plot.poscar.addnumbers<-function(poscar,layers=1,layer=1,direction=3,absolutenumber=T,...){
  atomselector=NULL
  indices <- poscar.getatomlayerindices(poscar,layers=layers)
  if(absolutenumber)
  {
    atoms <- poscar$atoms
    atomselector <- indices%in%layer
  }
  else
    atoms<-poscar$atoms[indices%in%layer,]
  plot.atoms.addnumbers(atoms,basis=poscar$a*poscar$basis,direction=direction,atomselector=atomselector,...)
}

#' Adds numbers on atoms to existing plot
#' 
#' \code{plot.atoms.addnumbers} adds numbers on atoms to existing plot.
#' 
#' @param atoms dataframe of atoms
#' @param direction of projection
#' @param basis basis if atoms are in direct coordinates
#' @param atomselector vector of atom indices which should be used
#' @param ... further plotting parameters
#' @export
plot.atoms.addnumbers<-function(atoms,basis=NULL,direction=3,atomselector=NULL,...){
  dir <- (1:3)[-direction]
  if (!is.null(basis))
    atoms <- atoms.convertbasis(atoms,basis,forth=F)
  atoms <- atoms[,dir]  
  if (is.null(atomselector))
    atomselector<- 1:nrow(atoms)
  text(atoms[atomselector,1],atoms[atomselector,2],labels=(1:nrow(atoms))[atomselector],xpd=T,...)
}

#' Adds arrows from specific atoms to other atoms to existing plot
#' 
#' \code{plot.atoms.addarrows} adds arrows from specific atoms to other atoms to existing plot.
#' Arrows will be drawn pairwise.
#' 
#' @param atomsold dataframe of atoms, used for arrow start
#' @param atomsnew dataframe of atoms, used for arrow end
#' @param direction of projection
#' @param basisold basis if atomsold are in direct coordinates
#' @param basisnew basis if atomsold are in direct coordinates
#' @param length size of arrow tips
#' @param ... further plotting parameters
#' @export
plot.atoms.addarrows<-function(atomsold,atomsnew,basisold=NULL,basisnew=NULL,direction=3,length=0.1,...){
  dir <- (1:3)[-direction]
  if (!is.null(basisold))
    atomsold <- atoms.convertbasis(atomsold,basisold,forth=F)
  if (!is.null(basisnew))
    atomsnew <- atoms.convertbasis(atomsnew,basisnew,forth=F)
  atoms <- atoms[c(1,2),dir]  
  arrows(atomsold[,dir[1]],atomsold[,dir[2]],atomsnew[,dir[1]],atomsnew[,dir[2]],length=length,...)
}

poscar.getbasisconvertedatoms<-function(poscar,forth=F){
  return(atoms.convertbasis(poscar$atoms, poscar$basis*poscar$a,forth))
}

#' Calculates basis converted atoms
#' 
#' \code{atoms.convertbasis} calculates basis converted atoms.
#' 
#' @param atoms dataframe of atoms
#' @param basis to use
#' @param forth from cartesian to direct
#' @export
atoms.convertbasis<-function(atoms,basis,forth=T){
  if (forth)
    basis <- solve(basis)
  if (ncol(atoms)>4)
    return (cbind(as.matrix(atoms[,1:3])%*%basis,atoms[,4:7]))
  else
    return (cbind(as.matrix(atoms[,1:3])%*%basis,atoms[,4]))
}

#' Gives the vacuum
#' 
#' \code{poscar.getvacuum} calculates the vacuum.
#' Vacuum is defined by the distance between atoms across unitcell borders.
#' 
#' @param poscar object of class poscar
#' @export
poscar.getvacuum <- function(poscar){
  # get data
  basis <- data.matrix(poscar$basis*poscar$a)
  atoms <- poscar$atoms
  # length of basis vectors
  size <- apply(basis,2,FUN=function(x)norm(as.matrix(x),"F"))
  # range of atoms for old vacuum
  rng <- apply(poscar$atoms[1:3],2,range)
  vacuumold<- c(0,0,0)
  for(i in 1:3)
  {
    span <- abs(rng[2,i]-rng[1,i])
    vacuumold[i] <- (1-span)*size[[i]]
  }
  return(vacuumold)
}

#' Sets the vacuum
#' 
#' \code{poscar.setvacuum} changes the vacuum.
#' Vacuum is defined by the distance between atoms across unitcell borders.
#' 
#' @param poscar object of class poscar
#' @param vacuum vector of new vacuum in all three directions. \code{c(0,0,10)} will change the vacuum in third direction to 10 (and leave the others).
#' @param center \code{TRUE} will center the cell after change
#' @export
poscar.setvacuum<- function(poscar,vacuum=c(0,0,0),center=T){
  if(center)
    poscar <- poscar.centeratoms(poscar)
  # get data
  vacuumold <- poscar.getvacuum(poscar)
  basis <- data.matrix(poscar$basis*poscar$a)
  atoms <- poscar$atoms
  # length of basis vectors
  size <- apply(basis,2,FUN=function(x)norm(as.matrix(x),"F"))
  # range of atoms for old vacuum
  rng <- apply(poscar$atoms[1:3],2,range)
  factor <- c(1,1,1)
  for(i in 1:3)
  {
    a <- atoms[,i]
    # calculate old vacuum
    span <- abs(rng[2,i]-rng[1,i])
    # if no new vacuum is given, use old one
    if (vacuum[[i]]==0)
      vacuum[[i]]<- vacuumold[[i]]
    # center all atoms in respect to the basis
    a <- a-rng[1,i]-span/2
    # calculate factor for each component
    factor[[i]] <- (vacuum[[i]]+span*size[[i]])/(vacuumold[[i]]+span*size[[i]])
    # apply new vacuum on atoms
    atoms[,i] <- a/factor[[i]]+rng[1,i]+span/2
  }
  # apply new vacuum on basis
  basis <- t(t(basis)*factor)
  # give data
  poscar$basis<-basis /poscar$a
  poscar$atoms <- atoms
  return(poscar)
}

#' Centers atoms
#' 
#' \code{atoms.centeratoms} centers the atoms relativ to a specific position.
#' 
#' @param atomsdirect atoms in direct coordinates
#' @param directions subset of 1,2,3 determines which directions should be centered
#' @param position relativ to which the atoms should be aranged (3d vector)
#' @export
atoms.centeratoms <- function(atomsdirect,directions=1:3,position=rep(0,3)){
  directions <- round(directions)
  if(any(directions<1|directions>3))
    warning("wrong direction encountered in atoms.centeratoms")
  # range of atoms for old vacuum
  rng <- apply(atomsdirect[1:3],2,range)
  position <- rep(position,length.out=length(directions))
  for(i in directions[directions>0&directions<4])
  {
    span <- abs(rng[2,i]-rng[1,i])
    atomsdirect[,i]<- atomsdirect[,i]-rng[1,i]-span/2+position[[which(directions==i)]]
  }
  return(atomsdirect)
}

#' Centers atoms in object of class poscar
#' 
#' \code{poscar.centeratoms} centers the atoms relativ to a specific position in object of class poscar.
#' 
#' @param poscar object of class poscar
#' @param directions subset of 1,2,3 determines which directions should be centered
#' @param position relativ to which the atoms should be aranged (3d vector)
#' @export
poscar.centeratoms <- function(poscar,directions=1:3,position=rep(0.5,3)){
  poscar$atoms <- atoms.centeratoms(atomsdirect=poscar$atoms,directions=directions,position=position)
  return(poscar)
}

#' Sorts atoms
#' 
#' \code{atoms.sort} sorts atoms.
#' 
#' @param atoms dataframe of atoms 
#' @param sortindices vector of indices after which to sort. \code{c(1,3)} will sort by rx then by rz
#' @export
atoms.sort<-function(atoms,sortindices=c(7,3,1,2)){
  return(atoms[do.call(order,lapply(sortindices,FUN=function(x)atoms[,x])),])  
}

#' Sorts atoms in object of class poscar
#' 
#' \code{poscar.sortatoms} sorts atoms in object of class poscar
#' 
#' @param poscar object of class poscar
#' @param sortindices vector of indices after which to sort. \code{c(1,3)} will sort by rx then by rz
#' @export
poscar.sortatoms <- function(poscar,sortindices=c(7,3,1,2))
{
  poscar$atoms <- atoms.sort(atoms=poscar$atoms,sortindices=sortindices)
  return(poscar)
}

poscar.extractatoms <- function(poscar,atomindices,vacuum=c(0,0,0),center=T)
{
  stopifnot(all(atomindices <= nrow(poscar$atoms)))
  poscar$atoms <- poscar$atoms[atomindices,]
  if (any(vacuum != c(0,0,0)))
    poscar <- poscar.setvacuum(poscar,vacuum,center=center)
  return(poscar)
}

poscar.extractlayers <- function(poscar,layer,layers,vacuum=c(0,0,0),center=T){
  layerindices <- poscar.getatomlayerindices(poscar,layers)
  return(poscar.extractatoms(poscar,which(layerindices %in% layer),vacuum,center=center))
}

poscar.removelayers <- function(poscar,layer,layers,vacuum=c(0,0,0),center=T){
  oldvacuum <- poscar.getvacuum(poscar)
  vacuum <- ifelse(vacuum==0,oldvacuum,vacuum)
  layerindices <- poscar.getatomlayerindices(poscar,layers)
  atoms <- poscar$atoms
  layer <- sort(layer,decreasing=T)
  for(l in layer){  
    if(!(l==1 | layer==layers)){
      zlayer <- mean(atoms[layerindices==l,3])
      zlayerabove <- mean(atoms[layerindices==(l+1),3])
      zlayerbelow <- mean(atoms[layerindices==(l-1),3])
      zoffset <- mean(abs(zlayer-c(zlayerabove,zlayerbelow)))      
      atoms[layerindices>l,3] <- atoms[layerindices>l,3]-zoffset      
    }
    atoms <- atoms[!(layerindices==l),]
    layerindices <- layerindices[!(layerindices==l)]
    layers <- layers-1 
    layerindices[layerindices>l] <- layerindices[layerindices>l]-1
  }
  poscar$atoms <- atoms
  poscar <- poscar.setvacuum(poscar=poscar,vacuum=vacuum,center=center)
  return(poscar)
}

poscar.getlayertranslation<-function(poscar,fromlayer,tolayer,layers){
  layerindices <- poscar.getatomlayerindices(poscar,layers)
  fatom <- poscar$atoms[(layerindices==fromlayer),1:2]  
  tatom <- poscar$atoms[(layerindices==tolayer),1:2]
  xmean <- mean(c(fatom[[1]],tatom[[1]]))
  ymean <- mean(c(fatom[[2]],tatom[[2]]))
  fatom <- fatom[which.min((xmean-fatom[[1]])^2+(ymean-fatom[[2]])^2),]
  tatom <- tatom[which.min((fatom[[1]]-tatom[[1]])^2+(fatom[[2]]-tatom[[2]])^2),]
  return(as.numeric(tatom-fatom))
}

poscar.translatelayers <- function(poscar,layer,layers,directtranslation){
  layerindices <- poscar.getatomlayerindices(poscar,layers)
  poscar$atoms[layerindices%in%layer,1:2] <- sweep(poscar$atoms[layerindices%in%layer,1:2],2,directtranslation,"+")%%1
  return(poscar)
}

poscar.getatomlayerindices<-function(poscar,layers){
  return (atoms.getlayerindices(poscar$atoms,poscar$basis*poscar$a,layers))
}

poscar.rotatelayer.deg<-function(poscar,layer,layers,angle,offset=NULL){
  return(poscar.rotatelayer.rad(poscar=poscar,layer=layer,layers=layers,angle=angle/180*pi,offset=offset))
}

poscar.rotatelayer.rad<-function(poscar,layer,layers,angle,offset=NULL){
  layerindices <- poscar.getatomlayerindices(poscar,layers)
  basis <- poscar$basis*poscar$a
  atoms <- poscar$atoms[layerindices%in%layer,]
  oldatoms <- atoms
  if(is.null(offset))
  {
    offset <- as.numeric(atoms[which.min(atoms[,1]^2+atoms[,2]^2),1:2])
  }
  else{
    offset <- as.numeric(offset)
  }
  stopifnot(length(offset)==2)
  acount<- nrow(atoms)
  atoms[,1:2] <- sweep(atoms[,1:2],2,offset,"-")
  atoms <- atoms.createsupercell(atomsdirect=atoms,super=c(5,5,1),center=F)  
  atoms[,1:2] <- atoms[,1:2]-2
  a <- atoms.convertbasis(atoms=atoms,basis=basis,forth=F)  
  matr <- rbind(c(cos(angle),-sin(angle)),c(sin(angle),cos(angle)))
  a[,1:2] <-t(matr%*%t(data.matrix(a[,1:2])))  
  atoms <- atoms.convertbasis(atoms=a,basis=basis,forth=T)
  atoms[,1:2] <- sweep(atoms[,1:2],2,offset,"+")  
  selector <- atoms[,1]>=0&atoms[,1]< 1&atoms[,2]>=0&atoms[,2]< +1
  atoms <- atoms[selector,]
  stopifnot(nrow(atoms)==acount)
  poscar$atoms[layerindices%in%layer,] <- atoms
  return(poscar)
}


atoms.getlayerindices<-function(atoms,basis=NULL,layers){
  if(layers==1)
    return(rep(1,nrow(atoms)))
  if(layers==nrow(atoms))
    return(order(atoms[,3]))
  stopifnot(nrow(atoms)>layers)
  stopifnot(nrow(atoms)>0)
  stopifnot(layers>0)
  if(!is.null(basis))
    atoms <- atoms.convertbasis(atoms,basis,forth=F)
  centers <- sort(unique(signif(atoms[,3],4)))
  while(length(centers)>layers)
  {
    o <- abs(outer(centers,centers,"-"))
    o[o==0]<-1000
    indices <- which(o==min(o),arr.ind=T)[1,]
    centers[indices[[1]]]<- mean(centers[indices])
    centers<- centers[-indices[[2]]]
  }
  k<-kmeans(atoms[,3],centers) 
  return (k$cluster)  
}

basis.rotate2d <- function(basis,angle=NULL){
  ### rotate
  if (is.null(angle))
    angle <- atan2(basis[1,2],basis[1,1])
  if (abs(angle) >1e-6)
  {
    print("rotating")
    matr <- rbind(c(cos(angle),-sin(angle)),c(sin(angle),cos(angle)))
    basis[,1:2]<-basis[,1:2]%*%matr
  }
  return(basis)
}

poscar.rotate2d <- function(poscar,angle=NULL){
  poscar$basis <- basis.rotate2d(poscar$basis,angle)
  return(poscar)
}

atoms.createsupercell <- function(atomsdirect,super=c(1,1,1),center=T,center.directions=1:3,center.position=rep(0,3)){
  super <- round(super-1)
  super <- rep(super,length.out=3)
  combinations <- expand.grid(lapply(super,seq,0))
  atoms <- apply(combinations,1,FUN=function(x)
  { 
    a <- atomsdirect[1:3]
    a <- sweep(a,2,x,"+")
    return(cbind(a,atomsdirect[4:ncol(atomsdirect)]))
  })
  atoms <- do.call(rbind,atoms)
  if(center)
    atoms <- atoms.centeratoms(atomsdirect=atoms,directions=center.directions,position=center.position)
  return(atoms)
}

poscar.createsupercell<-function(poscar,super=c(1,1,1),center=T,center.directions=1:3,center.position=rep(0.5,3)){
  poscar$atoms <- atoms.createsupercell(atomsdirect=poscar$atoms,super=super,center=F)
  atomsreal <- poscar.getbasisconvertedatoms(poscar)
  poscar$basis <- diag(super)%*%poscar$basis
  poscar$atoms <- atoms.convertbasis(atoms=atomsreal,basis=poscar$basis*poscar$a,forth=T)
  if(center)
    poscar <- poscar.centeratoms(poscar,directions=center.directions,position=center.position)
  return(poscar)
}

poscar.getpositionbyatom<-function(poscar,atomselector,zdist=c(1),layername=-1){
  positions <- list()
  atoms <- poscar.getbasisconvertedatoms(poscar)
  for(z in zdist)
  {
    atompos <- atoms[atomselector,1:3][1,]
    atompos[,3]<-atompos[,3]-z  
    d <- c(as.numeric(atompos),layername,z)
    names(d)<- c("x","y","z","layer","zdist")
    positions[[paste0("layer=",layername,"&zdist=",z)]]<-d
  }
  class(positions) <-"positions"
  return(positions)
}

poscar.getpositionbylayer<-function(poscar,layers,zdist=c(1),layer=layers){
  positions <- list()
  indices<-poscar.getatomlayerindices(poscar,layers=layers)  
  for(l in layer)
  {
    positions <- c(positions,poscar.getpositionbyatom(poscar,l==indices,zdist))  
  }
  class(positions) <-"positions"
  return (positions)
}

poscar.getreciprocalbasis<-function(poscar){
  return(basis.getreciprocal(poscar$basis,poscar$a))
}

basis.getreciprocal<-function(basis,a){
  basis <- basis*a
  volume <- c(crossprod.vec(basis[1,],basis[2,])%*%basis[3,])
  rbasis <- rbind(crossprod.vec(basis[2,],basis[3,]),
                 crossprod.vec(basis[3,],basis[1,]),
                 crossprod.vec(basis[1,],basis[2,]))
  rbasis <- rbasis/volume*2*pi
  
  return(rbasis)
}

#' gives you a poscar with all atoms in mirrorlayer
#' mirrored by the diagonal going trough a atom of the baselayer
#' 
#' \code{poscar.mirrorlayers} test
#' 
#' @param poscar the input poscar.
#' @param layers layercount of poscar
#' @param baselayer layer in which the mirror-diagonal will be layed
#' @param mirrorlayers layers which will be mirrored by the diagonal
#' @export
poscar.mirrorlayers<-function(poscar,layers,baselayer,mirrorlayers){
  indices <- poscar.getatomlayerindices(poscar,layers)
  ### find offset
  basisindex <- which(indices%in%baselayer)
  positionorder <- order(poscar$atoms[,1]+poscar$atoms[,2])
  index<-match(basisindex,positionorder)[[1]]
  offset <- as.numeric(poscar$atoms[positionorder[[index]],1:2])
  ### sweep offset
  atoms <- poscar$atoms[indices%in%mirrorlayers,1:2]
  atoms <- sweep(atoms,2,offset)
  atoms <- atoms[,2:1]
  atoms <- sweep(atoms,2,-offset)
  atoms <- atoms %%1
  poscar$atoms[indices%in%mirrorlayers,1:2] <- atoms
  return(poscar)
}

#' gives you a poscar with atoms in area of given size and shape
#' 
#' \code{poscar.getshapedposcar} test
#' 
#' @param poscar the input poscar.
#' @param x max x (range will be c(0,x))
#' @export
poscar.getshapedposcar <- function(poscar  # Input poscar
                                  ,x  # Max x (range will be c(0,x))
                                  ,y  # Max y (range will be c(0,y))
                                  ,shape=c("rectangular")  # Shape of area atoms are allowed
                                  ){  
  shape <- match.arg(shape)
  poscar <- poscar.rotate2d(poscar)
  basis <- poscar$basis*poscar$a
  atoms <- poscar.getbasisconvertedatoms(poscar)  
  super <- ceiling(max(x/basis[1,1],y/basis[2,2]))+1
  pairs <- t(do.call(cbind,lapply(-super:super,FUN=function(x)sapply(-super:super,FUN=function(y)c(x,y)))))
  atomcoords <- apply(pairs,1,FUN=function(x)
  {
    addvec <- x%*%basis[1:2,1:2]
    xatom <- atoms[,1]+addvec[1]
    yatom <- atoms[,2]+addvec[2]  
    return(list(cbind(xatom,yatom,atoms[,3:ncol(atoms)])))
  })
  atomcoords<-data.frame(do.call(rbind,lapply(atomcoords,FUN=function(x)x[[1]])))  
  names(atomcoords) <- c("x","y","z","movex","movey","movez","type")
  selector <- rep(T,nrow(atomcoords))
  selector <- selector & atomcoords$x<=x & atomcoords$x>0& atomcoords$y<=y & atomcoords$y>0
  print(sum(selector))
  
  poscar$atoms <- atoms.convertbasis(atoms=atomcoords[selector,],basis=basis,forth=T)
  return(poscar)
}
