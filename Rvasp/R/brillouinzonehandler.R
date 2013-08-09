
plot.brillouinzone.poscar<- function(poscar,rotate=0,extend=1,strain=0,xlab=NA,ylab=NA,xaxt="n",yaxt="n",col="blue",typ="l",...){
  vec2d<-poscar.brillouinzone.gethexagonalvector(poscar,rotate,extend,strain)  
  plot(apply(do.call(rbind,vec2d),2,range),asp=1,typ="n",xlab=xlab,ylab=ylab,xaxt=xaxt,yaxt=yaxt,...)
  if(typ!="n")
    plot.brillouinzone.hexagonalvector.add(vec2d,col=col,...)
}

poscar.brillouinzone.gethexagonalvector<-function(poscar,rotate=0,extend=1,strain=0){  
  basicpattern <- list(c(1,1,0),c(0,1,1))
  
  #poscar <- rotate2d.poscar(poscar)
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
  return(vec2d)
}

plot.brillouinzone.poscar.add<-function(poscar,rotate=0,extend=1,strain=0,col="blue",...){  
  plot.brillouinzone.hexagonalvector.add(poscar.brillouinzone.gethexagonalvector(poscar,rotate,extend,strain),col=col,...)
}

plot.brillouinzone.hexagonalvector.add<-function(hexvector,col="blue",...){  
  lapply(hexvector,FUN=function(x)polygon(x,border=col,...))
}

plot.brillouinzone.poscar.addsympoints<-function(poscar
                                                 ,rotate=0
                                                 ,extend=1
                                                 ,strain=0
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
  vec2d<-do.call(rbind,getbrillouinzonehexvector.poscar(poscar,rotate,extend,strain))
  lapply(1:length(directcoordinates),FUN=function(i){
    vec <- directcoordinates[[i]]%*%vec2d[vectors,]  
    #legend(vec[[1]],vec[[2]],legend=labels[i],seg.len=0.1,text.col=col,xjust=0.5,yjust=0.5,x.intersp=-0.4,y.intersp=0.3)
    text(vec[[1]]+xoffset[i],vec[[2]]+yoffset[i],labels=labels[i],pos=textpos,col=col,...)
    #grid.text(labels[i],vec[[1]],vec[[2]],default.units="native",...)
    if(typ=="p")
      points(vec[[1]],vec[[2]],col=col,pch=pch)
  })
  
}

