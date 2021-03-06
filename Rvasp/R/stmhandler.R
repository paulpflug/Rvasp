require(akima)
require(snowfall)
require(lattice)

#' Plots a constant-current stm image
#' 
#' \code{plot.stm} Creates a constant-current stm image.
#' Based on lattice package.
#' 
#' @param stm object of class stm
#' @param super (optional) plots a supercell (symmetric in x and y)
#' @param ... further plotting parameters
#' @export
plot.stm<- function(stm,super=1,xlim=NULL,ylim=NULL,more=F,useRaster=F                       
                    ,...){ 
  base <- stm$poscar$basis*stm$poscar$a
  ##supercell
  if (super>1)
  {
    print("calculating supercell")    
    tmpx <- lapply(1:(super-1),FUN=function(x){stm$x+x*base[1,1]})
    stm$x <- do.call(c,c(stm$x,tmpx))
    tmpy <- lapply(1:(super-1),FUN=function(x){stm$y+x*base[2,2]})
    stm$y <- do.call(c,c(stm$y,tmpy))
    stm$z <- do.call(rbind,rep(list(stm$z),super))  
    stm$z <- do.call(cbind,lapply(0:(super-1),FUN=function(x)
    {
      index <- nrow(stm$z)-which.min(abs(stm$x-x*base[2,1]))
      i <- 1:nrow(stm$z)
      if (x>0)
        i<- i[c((index):length(i),1:(index-1))]
      stm$z[i,]
    }))
  }
  ############
  if(is.null(xlim)) xlim <- range(stm$x)/10
  if(is.null(ylim)) ylim <- range(stm$y)/10
  xlim[2] <- stm$x[which.min(abs(stm$x-xlim[2]))]
  ylim[2] <- stm$y[which.min(abs(stm$y-ylim[2]))]
  
  grid <- expand.grid(x=stm$x/10, y=stm$y/10)
  grid$z <- as.vector(stm$z)
  zrange <- range(stm$z)
  zat <-pretty(zrange,n=5)
  zlabels <- paste(format(zat,digits=2),"\u00C5")
  fig <- levelplot(z~x+y,grid,pretty=F,panel= panel.levelplot.raster,cuts=100,
                   colorkey=list(T,raster=T,interpolate=T,labels=list(at=zat,labels=zlabels),width=1),
                   aspect="iso",contour=F,col.regions=colorRampPalette(c("orangered3","yellow"))(400),useRaster=useRaster,
                   xlim=xlim,ylim=ylim,...)
  print(fig,more=more,...)
}

#' Adds atoms
#' 
#' \code{plot.stm.addatoms} Adds atoms to a stm image
#' 
#' @param stm object of class stm
#' @param atomselector selector to filter atoms which should be added
#' @param super (optional) makes a supercell of filtered atoms (symmetric in x and y)
#' @param ... further plotting parameters
#' @export
plot.stm.addatoms<-function(stm,atomselector,col="black",atomsize=2*super,
                       super=1,xlim=NULL,ylim=NULL,...){
  ##atoms
  base <- stm$poscar$basis*stm$poscar$a
  atoms <- poscar.getbasisconvertedatoms(stm$poscar)   
  atoms <- data.matrix(atoms[atomselector,])
  if (ncol(atoms)==1)
    atoms <- t(atoms)
  pairs <- t(do.call(cbind,lapply(-super:super,FUN=function(x)sapply(-super:super,FUN=function(y)c(x,y)))))
  atomcoords <- apply(pairs,1,FUN=function(x)
  {
    addvec <- x%*%base[1:2,1:2]
    xatom <- atoms[,1]+addvec[1]
    yatom <- atoms[,2]+addvec[2]  
    return(list(cbind(xatom,yatom)))
  })
  atomcoords<-data.frame(do.call(rbind,lapply(atomcoords,FUN=function(x)x[[1]])))  
  names(atomcoords) <- c("x","y")
  selector <- rep(T,nrow(atomcoords))
  if (!is.null(xlim))
    selector <- selector & atomcoords$x<=xlim[2]*10 & atomcoords$x>xlim[1]*10
  if (!is.null(ylim))
    selector <- selector & atomcoords$y<=ylim[2]*10 & atomcoords$y>ylim[1]*10
  if(any(selector))
  {
    if (sum(selector)>1)
    coords <- atomcoords[selector,]/10
    else
      coords<- atomcoords/10
  trellis.focus("panel", 1, 1, highlight=F) 
  lpoints(coords,fill=col,col=grey(0.4),pch=21,lwd=1,cex=atomsize/super,...)
  trellis.unfocus()
  }
}

#' Adds unitcell
#' 
#' \code{plot.stm.addunitcell} Adds the unitcell to a stm image indicated by a dashed line
#' 
#' @param stm object of class stm
#' @param atomnumber (optional) number of the atom on which the bottom left corner of the unitcell will be positioned
#' @param ... further plotting parameters
#' @export
plot.stm.addunitcell<-function(stm,atomnumber=NULL,col="black",...){  
  offset <- c(0,0)
  if(!is.null(atomnumber))
  {
    atoms <- poscar.getbasisconvertedatoms(stm$poscar)  
    offset <-atoms[atomnumber,1:2]
  }
  base <- (stm$poscar$basis*stm$poscar$a)[1:2,1:2]  
  trellis.focus("panel", 1, 1, highlight=F) 
  points <- t(cbind(c(0,0),base[1,],base[1,]+base[2,],base[2,],c(0,0))+offset)/10
  llines(points,lty=2,col=col,...)
  trellis.unfocus()

}

#' Calculates a constant-current stm
#' 
#' 
#' \code{stm} Calculates a constant-current stm by a given chgcar.
#' preferred orientation: z-direction 
#'  
#' @param chgcar object of class chgcar
#' @param emax cut-off electron density
#' @param direction direction of stm creation (negativ will go from small values to big, positiv vice versa)
#' @param cpus uses snowfall package to parallelize calculation
#' @param interpolation only linear implemented
#' @export
stm<- function(chgcar,emax,direction=3,cpus=4,interpolation=c("linear")){  
  stm <- list()
  stm$direction <- direction
  topdown <- sign(direction)==1
  direction <- abs(round(direction))  
  dir <- (1:3)[-direction]
  if(chgcar$poscar$basis[dir[[1]],]%*%chgcar$poscar$basis[direction,]>=1e-6|chgcar$poscar$basis[dir[[2]],]%*%chgcar$poscar$basis[direction,]>=1e-6){
    stop("chosen direction is not perpendicular on other two directions")
  }
  print ("calculating stm")
  interpolation<-match.arg(interpolation)
  width <- switch(interpolation,linear=1)
  formula <- switch(interpolation,linear=formula(e~z))
  step <- c(-width:width)
  if(direction>1){
    step <- step*chgcar$n[1]
    if(direction>2){
      step <- step*chgcar$n[2]
    }
  }
  xy <- unique(chgcar$data[,dir])
  fac <- ceiling(log10(max(xy[,2])))
  nd <- chgcar$data[chgcar$data[,4]>=emax,]  
  nd <- nd[order(nd[,direction],decreasing=!topdown),]
 
  print("searching xy for given z")
  i <-match(xy[,1]+10^fac*xy[,2],nd[,dir[[1]]]+10^fac*nd[,dir[[2]]])
  nd<-nd[i,]
  fitfunc <- function(x)
  {
    x <- as.numeric(x)
    s <- step+x
    xy <- chgcar$data[x,dir]
    z <-  chgcar$data[x,direction]
    data <- chgcar$data[s,c(direction,4)]
    names(data)<-c("z","e")
    data[,2] <- data[,2]-emax    
    o <- lm(formula,data)
    zeros <- polyroot(o$coefficients)
    zeros <- zeros[abs(Im(zeros))<1e-4]
    index <- which.min(abs(Re(zeros)-z))        
    zero <- Re(zeros[[index]])
    return(data.frame(xy,zero))
  }
  print("interpolating z")
  if (cpus>1)
  {
    sfInit(parallel=TRUE, cpus=cpus, type="SOCK")
    sfExport("chgcar")
    data <- sfLapply(rownames(nd),fun=fitfunc)
    sfRemoveAll()
    sfStop()   
  }
  else
  {
    print("on single cpu")
    data <- lapply(rownames(nd),FUN=fitfunc)
  }
  data <- do.call(rbind,data)
  print("interpolating z finished")
  # ranges of z to find corresponding atoms for plotting purpose
  atomrange<- range(data[,direction])
  # lowest to zero
  data[,direction]<-data[,direction]-atomrange[[1]]
  # reverse in z direction
  if(topdown){
    data[,direction]<-abs(data[,direction]-max(data[,direction]))
  }
  
  print("interpolating quadratic area")
  base <- chgcar$poscar$basis*chgcar$poscar$a  
  ### rotate
  angle <- -acos(base[1,1]/sqrt(base[1,1]^2+base[1,2]^2))
  if (abs(angle) >1e-6)
  {
    print("rotating")
    matr <- rbind(c(cos(angle),-sin(angle)),c(sin(angle),cos(angle)))
    base[dir,dir]<-base[dir,dir]%*%matr
    data[,dir]<-data.matrix(data[,dir])%*%matr    
  }
  ## preparing quadratic area
  pairs <- rbind(c(1,0),c(0,1),c(0,0),c(-1,0),c(0,-1),c(-1,-1),c(1,1),c(-1,1),c(1,-1))
  addfunc <- function(x)
  {
    addvec <- x%*%base[dir,dir]
    d2 <- data
    d2[,dir[[1]]] <- data[,dir[[1]]]+addvec[1]
    d2[,dir[[2]]] <- data[,dir[[2]]]+addvec[2]
    return (d2)
  }
  dtmp <- apply(pairs,1,FUN=addfunc)  
  dtmp <- do.call(rbind,dtmp)
  xmin <- 0#min(data[,1])
  ymin <- 0#min(data[,2])
  xmax <- xmin+abs(base[1,dir[[1]]])
  ymax <- ymin+abs(base[2,dir[[2]]])
  xo=seq(xmin,xmax,length=201)[1:200]
  yo=seq(ymin,ymax,length=201)[1:200] 
  ## calculation of quadratic area
  o<-interp(dtmp[,dir[[1]]],dtmp[,dir[[2]]],dtmp[,direction],xo=xo,yo=yo)
  print("calculating stm successfull")
  
  stm$poscar <- chgcar$poscar
  stm$poscar$basis <- base/chgcar$poscar$a
  stm$x <- o$x
  stm$y <- o$y
  stm$z <- o$z
  stm$emax <- emax
  stm$atomrange<-atomrange
  class(stm)<-"stm"
  gc()
  return (stm)
}
