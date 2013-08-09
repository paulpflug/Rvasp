#' Plots brillouinzone
#' 
#' \code{plot.brillouinzone} plots brillouinzone.
#' 
#' @param brillouinzone object of class brillouinzone
#' @param ... further plotting parameters
#' @export
plot.brillouinzone<- function(brillouinzone,
                              xlab=NA,
                              ylab=NA,
                              xaxt="n",
                              yaxt="n",
                              col="blue",
                              typ="l",
                              ...){
  plot(apply(do.call(rbind,brillouinzone),2,range),asp=1,typ="n",xlab=xlab,ylab=ylab,xaxt=xaxt,yaxt=yaxt,...)
  if(typ!="n")
    plot.brillouinzone.add(brillouinzone,col=col,...)
}

#' Gets hexagonal brilloinzone vector
#' 
#' \code{brillouinzone.gethexagonalvector} gets hexagonal brilloinzone vector.
#' 
#' @param poscar object of class poscar
#' @param rotate rotates brillouinzone (in degrees)
#' @param extend creates supercell of brillouinzone 
#' @param strain applied to brillouinzone
#' @export
poscar.getbrillouinzone.hexagonal<-function(poscar,rotate=0,extend=1,strain=0){  
  basicpattern <- list(c(1,1,0),c(0,1,1))
  poscar$a <- poscar$a*(1-strain)
  rbase <- reciprocal.poscar(poscar)[1:2,1:2]%*%get2droationmatrix.deg(rotate)
  vec2d <- rbind(rbase,(rbase[1,]+rbase[2,]),-rbase,-(rbase[1,]+rbase[2,]))
  vec2d <- vec2d[c(5,6,4,2,3,1),]
  if (extend==1)
  vec2d<- list(vec2d)
  else
  {
    patterns <- unlist(sapply(1:(extend),FUN=function(e){
      e1 <- 1:e
      e2 <- (e:1)-1
      l <- lapply(1:length(e1),FUN=function(i){
        e1[i]*basicpattern[[1]]+e2[i]*basicpattern[[2]]
        
      })
      return(l)
      }),F)
    vec2d<- c(vec2d,sapply(patterns[1:(sum(1:(extend-1)))],FUN=function(pattern){
      lapply(1:6,FUN=function(y)
      {    
        z <- rep(0,6)
        z[((y:(y+length(pattern)-1))%%6)+1]<-pattern
        add.vectomatrix(vec2d,apply(vec2d*z,2,sum)) 
      })})
                           
                           )
  }
  class(vec2d)<-"brillouinzone"
  return(vec2d)
}

#' Adds brilloinzonevector to existing plot
#' 
#' \code{plot.brillouinzone.addvector} adds brilloinzonevector to existing plot based on poscar.
#' 
#' @param brillouinzone object of class brillouinzone
#' @param ... further plotting parameters
#' @export
plot.brillouinzone.add<-function(brillouinzone,col="blue",...){  
  lapply(brillouinzone,FUN=function(x)polygon(x,border=col,...))
}

#' Adds high symmetry points to existing plot
#' 
#' \code{plot.brillouinzone.addsympoints} adds high symmetry points to existing plot based on poscar.
#' 
#' @param brillouinzone object of class brillouinzone
#' @param directcoordinates positions of highsymmetry points in direct coordinates
#' @param labels labels of highsymmetry points in order given by \code{directcoordinates}
#' @param vectors indices of two vectors from brilloinzone, which will be used as basis for calculation
#' @param textpos position of labels, see \code{\link{text}} for futher information
#' @param ... further plotting parameters
#' @export
plot.brillouinzone.addsympoints<-function(brillouinzone
                                                 ,directcoordinates=list(c(0,0))
                                                 ,labels=1:length(directcoordinates)
                                                 ,vectors=1:2
                                                 ,textpos=2
                                                 ,col="blue"
                                                 ,pch=16
                                                 ,typ="p"
                                                 ,xoffset=0
                                                 ,yoffset=0
                                                 ,...){  
  if (length(xoffset)<length(directcoordinates))
    xoffset <- rep(xoffset,ceiling(length(directcoordinates)/length(xoffset)))
  if (length(yoffset)<length(directcoordinates))
    yoffset <- rep(yoffset,ceiling(length(directcoordinates)/length(yoffset)))
  library(grid)
  brillouinzone<-do.call(rbind,brillouinzone)
  lapply(1:length(directcoordinates),FUN=function(i){
    vec <- directcoordinates[[i]]%*%brillouinzone[vectors,]  
    #legend(vec[[1]],vec[[2]],legend=labels[i],seg.len=0.1,text.col=col,xjust=0.5,yjust=0.5,x.intersp=-0.4,y.intersp=0.3)
    text(vec[[1]]+xoffset[i],vec[[2]]+yoffset[i],labels=labels[i],pos=textpos,col=col,...)
    #grid.text(labels[i],vec[[1]],vec[[2]],default.units="native",...)
    if(typ=="p")
      points(vec[[1]],vec[[2]],col=col,pch=pch)
  })
  
}

