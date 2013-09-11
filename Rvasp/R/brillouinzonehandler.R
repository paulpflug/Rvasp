#' Plots brillouinzones
#' 
#' \code{plot.brillouinzones} plots brillouinzones.
#' 
#' @param brillouinzones object of class brillouinzones
#' @param ... further plotting parameters
#' @export
plot.brillouinzones<- function(brillouinzone,
                              xlab=NA,
                              ylab=NA,
                              xaxt="n",
                              yaxt="n",
                              col="blue",
                              typ="l",
                              ...){
  plot(apply(do.call(rbind,brillouinzone),2,range),asp=1,typ="n",xlab=xlab,ylab=ylab,xaxt=xaxt,yaxt=yaxt,...)
  if(typ!="n")
    plot.brillouinzones.add(brillouinzone,col=col,...)
}

#' Adds brilloinzonevector to existing plot
#' 
#' \code{plot.brillouinzonse.addvector} adds brilloinzonevector to existing plot based on poscar.
#' 
#' @param brillouinzones object of class brillouinzones
#' @param ... further plotting parameters
#' @export
plot.brillouinzones.add<-function(brillouinzones,col="blue",...){  
  invisible(lapply(brillouinzone,FUN=function(x)polygon(x,border=col,...)))
}

#' Adds high symmetry points to existing plot
#' 
#' \code{plot.brillouinzones.addsympoints} adds high symmetry points to existing plot based on poscar.
#' 
#' @param brillouinzones object of class brillouinzones
#' @param directcoordinates positions of highsymmetry points in direct coordinates
#' @param labels labels of highsymmetry points in order given by \code{directcoordinates}
#' @param vectors indices of two vectors from brilloinzone, which will be used as basis for calculation
#' @param textpos position of labels, see \code{\link{text}} for futher information
#' @param ... further plotting parameters
#' @export
plot.brillouinzones.addsympoints<-function(brillouinzones
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
  invisible(lapply(1:length(directcoordinates),FUN=function(i){
    vec <- directcoordinates[[i]]%*%brillouinzone[vectors,]  
    #legend(vec[[1]],vec[[2]],legend=labels[i],seg.len=0.1,text.col=col,xjust=0.5,yjust=0.5,x.intersp=-0.4,y.intersp=0.3)
    text(vec[[1]]+xoffset[i],vec[[2]]+yoffset[i],labels=labels[i],pos=textpos,col=col,...)
    #grid.text(labels[i],vec[[1]],vec[[2]],default.units="native",...)
    if(typ=="p")
      points(vec[[1]],vec[[2]],col=col,pch=pch)
  }))  
}

#' Gets brillouinzone vector
#' 
#' \code{basis.getbrillouinzone} gets brillouinzone vector.
#' 
#' @param basis basis in real space
#' @export
basis.getbrillouinzone<-function(basis){
  recbasis <- basis.getreciprocal(basis)
  return(reciprocalbasis.getbrillouinzone(recbasis))
}

#' Gets brillouinzone vector
#' 
#' \code{reciprocalbasis.getbrillouinzone} gets brillouinzone vector.
#' 
#' @param recbasis basis in reciprocal space
#' @export
reciprocalbasis.getbrillouinzone<-function(recbasis){
  baseanglexy <- round(vectors.calcangle.degree(recbasis[1,],recbasis[2,],period=180),2)
  type <- reciprocalbasis.getbrillouinzonetype(recbasis)
  if(!is.na(type)){
    rbasis <- recbasis[1:2,1:2]
    if(type=="rectangular"){
      vec2d <- rbind(c(0,0),rbasis,rbasis[1,]+rbasis[2,])
      vec2d <- scale(vec2d,scale=F)
      vec2d <- vec2d[c(1,3,4,2),]
      class(vec2d)<-"brillouinzone"
      return(vec2d)
    }
    if(type=="hexagonal"){
      if(baseanglexy==120){
        rbasis <- rbind(rbasis[1,],rbasis[1,]+rbasis[2,],rbasis[2,])
      }
      else{
        rbasis <- rbind(rbasis[1,],rbasis[2,],rbasis[1,]-rbasis[2,])
      }
      rbasis <- (rbasis/2/cos(30/180*pi))%*%matrix.rotation2d.degree(30)
      vec2d <- rbind(rbasis,-rbasis)
      class(vec2d)<-"brillouinzone"
      return(vec2d)
    }
  }
  print("calculation of your brillouinzone is not implemented")
  return(numeric(0))
}

#' Gets type of brillouinzone
#' 
#' \code{reciprocalbasis.getbrillouinzone} gets type of brillouinzone.
#' 
#' @param recbasis basis in reciprocal space
#' @export
reciprocalbasis.getbrillouinzonetype<-function(recbasis){
  baseanglexy <- round(vectors.calcangle.degree(recbasis[1,],recbasis[2,],period=180),2)
  baseanglexz <- round(vectors.calcangle.degree(recbasis[1,],recbasis[3,],period=180),2)
  baseangleyz <- round(vectors.calcangle.degree(recbasis[2,],recbasis[3,],period=180),2)
  print(paste0("found following angles: xy=",baseanglexy,"° xz=",baseanglexz,"° yz=",baseangleyz,"°"))
  if(baseanglexz==90 & baseangleyz==90){
    rbasis <- recbasis[1:2,1:2]
    if(baseanglexy==90){
      print("found rectangular brillouinzone") 
      return("rectangular")
    }
    if(baseanglexy==60 | baseanglexy==120){
      print("found hexagonal brillouinzone")
      return("hexagonal")
    }
  }
  print("calculation of your brillouinzone is not implemented")
  return(NA)
}

#' Gets brillouinzone vector
#' 
#' \code{poscar.getbrillouinzones} gets vector containing the vertices of the brillouinzones.
#' 
#' @param poscar object of class poscar
#' @param rotate rotates brillouinzone (in degrees)
#' @param extend creates supercell of brillouinzone 
#' @param strain applied to brillouinzone
#' @export
poscar.getbrillouinzones<-function(poscar,rotate=0,extend=1,strain=0){ 
  poscar$a <- poscar$a*(1-strain)
  poscar$basis[1:2,1:2] <- poscar$basis[1:2,1:2]%*%matrix.rotation2d.degree(-rotate)
  rbase <- poscar.getreciprocalbasis(poscar)
  vec2d <- list(reciprocalbasis.getbrillouinzone(rbase))
  if (extend>1){
    type <- reciprocalbasis.getbrillouinzonetype(rbase)
    if(type=="hexagonal"){
      basicpattern <- list(c(1,1,0),c(0,1,1))    
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
          return(sweep(vec2d[[1]],2,apply(vec2d[[1]]*z,2,sum),"+")) 
        })}))
    }
    if(type=="rectangular"){
      ex <- seq(0,(2*extend-1),by=2)
      p <- -(ex[extend]):(ex[extend])
      patterns <- unlist(apply(expand.grid(p,p),1,function(x){
        if(sum(abs(x))%in%ex){
          return(list(as.numeric(x)))
        }else{
          return(NULL)      
        }
      }),recursive=F)
      vec2d<- lapply(patterns,FUN=function(pattern){
        m <-vec2d[[1]]
        p <- as.numeric(pattern%*%m[c(1,2),])
        m <- sweep(m,2,p)
        return(m)
      })
    }
  }
  class(vec2d)<-"brillouinzones"
  return(vec2d)
}

#' Gives brillouinzone projected kpoints.
#' 
#' \code{brillouinzone.projectkpoints} projects arbitrary kpoints back in the first brillouinzone.
#' 
#' @param brillouinzone object of class brillouinzone
#' @param kpoints object with arbitrary kpoints
#' @export
brillouinzone.projectkpoints <- function(brillouinzone,kpoints){
  newpoints <- apply(kpoints,1,FUN=function(p){
    distances <- order(sqrt((p[[1]] - brillouinzone[,1])^2 + (p[[2]] - brillouinzone[,2])^2)) # vertex distances from point
    v1 <- brillouinzone[distances[[1]],] # vertex 1
    v2 <- brillouinzone[distances[[2]],] # vertex 2
    dp <- vector.length(p) # distance 0P
    v3 <- v1+sum((p-v1)*(v2-v1))/vector.length(v2-v1)*(v2-v1)/vector.length(v2-v1) # vertex of the height of the triangle of v1,v2 and p 
    d3 <- vector.length(v3)
    if(dp>d3){                
      pnew <- 2*(v3-p)+p  # Two times the height, to get in the BZ
      return(pnew)
    }else{
      return(p)
    }    
  })
  return(t(newpoints))
}

#' Select kpoints in first brillouinzone.
#' 
#' \code{brillouinzone.selectkpoints} selects kpoints in first brillouinzone.
#' 
#' @param brillouinzone object of class brillouinzone
#' @param kpoints object with arbitrary kpoints
#' @param maxdistance (optional) allowed distance to first brillouinzone (shape-conserving)
#' @export
brillouinzone.selectkpoints <- function(brillouinzone,kpoints,maxdistance=0){
  newpoints <- apply(kpoints,1,FUN=function(p){
    distances <- order(sqrt((p[[1]] - brillouinzone[,1])^2 + (p[[2]] - brillouinzone[,2])^2)) # vertex distances from point
    v1 <- brillouinzone[distances[[1]],] # vertex 1
    v2 <- brillouinzone[distances[[2]],] # vertex 2
    dp <- vector.length(p) # distance 0P
    v3 <- v1+sum((p-v1)*(v2-v1))/vector.length(v2-v1)*(v2-v1)/vector.length(v2-v1) # vertex of the height of the triangle of v1,v2 and p 
    d3 <- vector.length(v3)
    if(dp>d3){   
      if(vector.length(v3-p)<maxdistance)
        return(p)
      else
        return(NULL)
    }else{
      return(p)
    }    
  })
  return(t(newpoints))
}

