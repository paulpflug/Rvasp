#' Plots brillouinzones
#' 
#' \code{plot.brillouinzones} plots brillouinzones.
#' 
#' @param brillouinzones object of class brillouinzones
#' @param ... further plotting parameters
#' @export
plot.brillouinzones<- function(brillouinzones,
                              xlab=NA,
                              ylab=NA,
                              xaxt="n",
                              yaxt="n",
                              col="blue",
                              typ="l",
                              ...){
  plot(apply(do.call(rbind,brillouinzones),2,range),asp=1,typ="n",xlab=xlab,ylab=ylab,xaxt=xaxt,yaxt=yaxt,...)
  if(typ!="n")
    plot.brillouinzones.add(brillouinzones,col=col,...)
}

#' Adds brilloinzonevector to existing plot
#' 
#' \code{plot.brillouinzonse.addvector} adds brilloinzonevector to existing plot based on poscar.
#' 
#' @param brillouinzones object of class brillouinzones
#' @param ... further plotting parameters
#' @export
plot.brillouinzones.add<-function(brillouinzones,col="blue",...){  
  invisible(lapply(brillouinzones,FUN=function(x)polygon(x,border=col,...)))
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
  brillouinzones<-do.call(rbind,brillouinzones)
  invisible(lapply(1:length(directcoordinates),FUN=function(i){
    vec <- directcoordinates[[i]]%*%brillouinzones[vectors,]  
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
#' @param zone specific the wished zone. \code{bz} will calculate the first brillouinzone, \code{basis} will give the simple zone conducted by the reciprocal basis vectors
#' @export
reciprocalbasis.getbrillouinzone<-function(recbasis,zone=c("bz","basis")){
  zone <- match.arg(zone)
  if(zone=="basis"){
    vec2d <- rbind(c(0,0),recbasis[1,1:2],recbasis[1,1:2]+recbasis[2,1:2],recbasis[2,1:2])
    vec2d <- scale(vec2d,center=colMeans(vec2d),scale=F)
    class(vec2d)<-"basiszone"
    return(vec2d)
  }
  if(zone=="bz"){ 
    baseanglexy <- round(vectors.calcangle.degree(recbasis[1,],recbasis[2,],period=180),1)
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
        rbasis <- rbasis/2
        lineintercept2d <- function(p1,d1,p2,d2){
          if(d1[[2]]==0){
            p3 <- p2
            p2 <- p1
            p1 <- p3
            d3 <- d2
            d2 <- d1
            d1 <- d3          
          }
          h <- d1[[1]]/d1[[2]]
          s <- (h*(p1[[2]]-p2[[2]])+p2[[1]]-p1[[1]])/(h*d2[[2]]-d2[[1]])
          return(p2+s*d2)
        }
        i1 <- lineintercept2d(rbasis[1,],rbasis[1,2:1]*c(1,-1),rbasis[2,],rbasis[2,2:1]*c(1,-1))
        i2 <- lineintercept2d(rbasis[3,],rbasis[3,2:1]*c(-1,1),rbasis[2,],rbasis[2,2:1]*c(-1,1))
        i3 <- lineintercept2d(-rbasis[1,],-rbasis[1,2:1]*c(1,-1),rbasis[3,],rbasis[3,2:1]*c(1,-1))
        is <- rbind(i1,i2,i3)
        is <- unique(round(is,7))
        vec2d <- rbind(is,(is%*%matrix.rotation2d.degree(180))[,])
        class(vec2d)<-"brillouinzone"
        return(vec2d)
      }
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
  baseanglexy <- round(vectors.calcangle.degree(recbasis[1,],recbasis[2,],period=180),1)
  baseanglexz <- round(vectors.calcangle.degree(recbasis[1,],recbasis[3,],period=180),1)
  baseangleyz <- round(vectors.calcangle.degree(recbasis[2,],recbasis[3,],period=180),1)
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
    if(nrow(vec2d[[1]])==6){      
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
    if(nrow(vec2d[[1]])==4){
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
#' @export
brillouinzone.selectkpoints <- function(brillouinzone,kpoints){
  lengths <- apply(brillouinzone,1,vector.length)
  newpoints <- apply(kpoints,1,FUN=function(point){
    p <- as.numeric(point[1:2])
    anglesum <-0
    for(i in 1:nrow(brillouinzone)){
      i2 <- (i+1)
      if(i==nrow(brillouinzone)) i2 <- 1
      v1 <- brillouinzone[i,] # vertex 1
      v2 <- brillouinzone[i2,] # vertex 2
      if(all(v1==p)|all(v2==p)){
        return(point)
      }
      angle <- suppressWarnings(vectors.calcangle.degree(p-v1,p-v2,period=180))
      if(is.nan(angle)) return(point)
      if(abs(angle-180)<1e-6){ ## v1 v2 and p on one line
        a <- ((p-v1)/(v2-v1))
        a <- a[a!=0&is.finite(a)]
        if(!is.null(a)){
          a <- a[[1]]
          if(0<=a&a<=1)
            return(point) 
        }                      
      } 
      anglesum <- anglesum+angle
    }
    if(abs(anglesum-360)<1e-6){
      return(point)
    }else{
      return(NULL)
    }
  })
  newpoints <- do.call(rbind,newpoints)
  return(newpoints)
}

#' Extends kpoints to all second brillouinzones.
#' 
#' \code{brillouinzone.extendkpoints} translates kpoints to all second brillouinzones.
#' 
#' @param brillouinzone object of class brillouinzone
#' @param kpoints object with arbitrary kpoints
#' @export
brillouinzone.extendkpoints <- function(brillouinzone,kpoints){
  names <- colnames(kpoints)
  len <- nrow(brillouinzone)
  permutations <- outer(1:2,(1:len)-1,FUN="+")%%len+1
  newkpoints <- lapply(1:len,FUN=function(i){
    vecindices <- permutations[,i]
    vec <- colSums(brillouinzone[vecindices,])
    np <- sweep(kpoints[,1:2],2,vec,FUN="+")
    if(ncol(kpoints)>2)
      np <- cbind(np,kpoints[,3:ncol(kpoints)])
    return(np)
  })
  newkpoints <- do.call(rbind,newkpoints)
  colnames(newkpoints)<-names
  return(rbind(kpoints,newkpoints))
}
