#' prints a two python vectors for x and y direction with equidistant points within and on the left and bottom border of a parallelogram
#' 
#' \code{getrasteredparallelogram} calculates equidistant points within a parallelogram
#' 
#' @param direction1 first 2d direction
#' @param direction2 second 2d direction
#' @param discretisation1 number of steps in first direction
#' @param discretisation2 number of steps in second direction
#' @param significant prettier printing by use of signif
#' @export
getrasteredparallelogram<-function(direction1,direction2,discretisation1=5,discretisation2=5,significant=5){
  disc1 <- seq(0,1,len=discretisation1+1)
  disc2 <- seq(0,1,len=discretisation2+1)
  disc1 <- disc1[-length(disc1)]
  disc2 <- disc2[-length(disc2)]
  i <- rep(disc1,discretisation2)
  j <- as.numeric(sapply(disc2,rep,discretisation1))
  x<-sapply(i,FUN=function(x)
    as.numeric(direction1)*x
    )
  y<-sapply(j,FUN=function(y)
    as.numeric(direction2)*y
  )
  print(apply(signif(x+y,significant),MARGIN=1,FUN=function(z)paste("[",paste(z,collapse=","),"]")))
}

#' adds a label to a plot with the help of the package grid
#'  
#' \code{plot.addlabel} adds a label to a plot
#' 
#' @param label label
#' @param top distance to top in mm
#' @param left distance to left in mm
#' @export
plot.addlabel<-function(label,top=0.5,left=0){
  library(grid) 
  grid.text(label,
           x=unit(left, "mm"), y=unit(1, "npc")-unit(top, "mm"), 
            just=c("left", "top")) 
}

#' Calculated the crossproduct of two given 3d vectors
#'  
#' \code{crossprod.vec} calculated the crossproduct
#' 
#' @param vector1 first vector
#' @param vector2 second vector
#' @export
vectors.crossproduct<-function(vector1,vector2){
  vec3 <- c(vector1[2]*vector2[3]-vector1[3]*vector2[2],
            vector1[3]*vector2[1]-vector1[1]*vector2[3],
            vector1[1]*vector2[2]-vector1[2]*vector2[1]
  )
  return(vec3)
}

#' Gets the euclidean length of a vector
#' 
#' \code{vector.length} gets the euclidean length of a vector.
#' @param vector the vector
#' @export
vector.length <- function(vector){   
  return(sqrt(sum(vector*vector)))
}

#' Gets the angle between two vectors
#' 
#' \code{vectors.calcangle} calcultes the angle between two vectors.
#' Gives angle in radians.
#' 
#' @param vector1 first vector
#' @param vector2 second vector
#' @param period (optional) if provided folds the angle back to its period. Period in radians.
#' @export
vectors.calcangle <- function(vector1,vector2,period=NA){
  a <- acos( sum(vector1*vector2) / ( vector.length(vector1) * vector.length(vector2)) ) 
  if(!is.na(period)&&!is.na(a)){
    while(a<0||a>period){
      if(a<0)
        a <- a+period
      if(a>period)
        a <- a-period
    }   
  }
  return(a)
}

#' Gets the angle between two vectors
#' 
#' \code{vectors.calcangle.degree} calcultes the angle between two vectors.
#' Gives angle in degrees.
#' 
#' @param vector1 first vector
#' @param vector2 second vector
#' @param period (optional) if provided folds the angle back to its period. Period in degrees.
#' @export
vectors.calcangle.degree <- function(vector1,vector2,period=NA){
  return(vectors.calcangle(vector1,vector2,period/180*pi)/pi*180 )
}

#' Checks two vectors for collinearity
#' 
#' \code{vectors.arecollinear} checks two vectors for collinearity.
#' 
#' @param vector1 first vector
#' @param vector2 second vector
#' @param prec (optional) precission
#' @export
vectors.arecollinear <- function(vector1,vector2,prec=1e-6){
  facs <- vector1/vector2
  return(all(abs(1-facs/facs[[1]])<prec))
}

#' Gets 2d rotation matrix
#' 
#' \code{matrix.rotation2d} gets 2d rotation matrix for given angle.
#' 
#' @param angle angle in radians
#' @export
matrix.rotation2d<-function(angle)
  return(rbind(c(cos(angle),-sin(angle)),c(sin(angle),cos(angle))))

#' Gets 2d rotation matrix
#' 
#' \code{matrix.rotation2d.degree} gets 2d rotation matrix for given angle.
#' 
#' @param angle angle in degrees
#' @export
matrix.rotation2d.degree<-function(angle)
  return(matrix.rotation2d(angle/180*pi))

#' Gets matrix for reflection across a line
#' 
#' \code{matrix.reflection2d} gets matrix for reflection across a line which goes through zero.
#' 
#' @param slope given in degrees
#' @export
matrix.reflection2d.degree<-function(slope){
  return(matrix.reflection2d(slope/180*pi))
}

#' Gets matrix for reflection across a line
#' 
#' \code{matrix.reflection2d} gets matrix for reflection across a line which goes through zero.
#' 
#' @param slope given in radians
#' @export
matrix.reflection2d<-function(slope){
  return(cbind(c(cos(2*slope),sin(2*slope)),c(sin(2*slope),-cos(2*slope))))
}


#' makes a vector of colors transparent
#' snippet from stackoverflow.com
#'  
#' \code{makeTransparent} makes a vector of colors transparent
#' 
#' @param someColor colors
#' @param alpha transparency
#' @export
makeTransparent<-function(someColor, alpha=100){
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

#' gets the fig variables of a plotinplot
#'  
#' \code{getplotinplotfig} gets the fig variables of a plotinplot
#' 
#' @param x x position
#' @param y y position
#' @param coords coordinate type of positions
#' @export
getplotinplotfig <- function(x,y,coords="ndc"){
  x <- range(x)
  y <- range(y)
  if(exists("paraold"))
  {
    usr <- paraold$usr
  }
  else
  {
    usr <- par("usr")
  }
  xstep <- (usr[2]-usr[1])/100
  ystep <- (usr[4]-usr[3])/100
  x <- c(usr[1]+x[[1]]*xstep,usr[1]+x[[2]]*xstep)
  y <- c(usr[3]+y[[1]]*ystep,usr[3]+y[[2]]*ystep)
  fig <- c(grconvertX(x, from="user", to=coords),
           grconvertY(y, from="user", to=coords))
  return(fig)
}

#' creates a new plot within a plot
#' uses base package
#'  
#' \code{plotinplot} creates a new plot within a plot
#' 
#' @param x x position
#' @param y y position
#' @param mar margin for new plot
#' @export
plotinplot <- function(x,y,mar=c(0,0,0,0)){
  if(!exists("paraold"))
  {
    print("saving parameters")
    paraold <<- par(no.readonly=T)
  }
  par(paraold)
  fig <- getplotinplotfig(x,y)
  par(fig = fig,
      mar = mar,
      new = TRUE)
}

#' sets the focus back to old plot
#'  
#' \code{endplotinplot} sets the focus back to old plot
#' @export
endplotinplot <- function(){
  par(paraold)
  rm(paraold,pos=1)
}

#' creates a new zoomplot within a plot
#' uses base package
#'  
#' \code{zoomplot} creates a new plot within a plot
#' 
#' @param xpos x position of the plot
#' @param ypos y position of the plot
#' @param xmean new plotting region will be centered around xmean
#' @param xwidth width of the new plotting region in x direction
#' @param ymean new plotting region will be centered around ymean
#' @param ywidth width of the new plotting region in y direction (optional if asp is given)
#' @param lines will draw lines from old plot to new plotting region
#' @param box will draw box around old plot
#' @param asp (optional) is used to calculate ywidth if not given. if both are not given the asp of original plot is used.
#' @param mar margin for new plot
#' @param ... further plotting parameters
#' @export
zoomplot<-function(xpos,
                   ypos,
                   xmean,
                   xwidth,
                   ymean,
                   ywidth=NULL,
                   lines=T,
                   box=T,
                   mar=c(0,0,0,0),
                   xlab="",
                   ylab="",
                   bg="",
                   asp=NULL,
                   xticks=2,
                   yticks=2,
                   ...){
  xlim <- c(xmean-xwidth/2,xmean+xwidth/2)
  usr <- par("usr")
  pin <- par("pin")
  w <- pin[1]/diff(usr[1:2])
  h <- pin[2]/diff(usr[3:4])
  
  if(is.null(ywidth))
  {
    if(is.null(asp))
      asp <- w/h
    ylim <- ymean + (xlim-xmean)/asp
  }
  else 
    ylim <- c(ymean-ywidth/2,ymean+ywidth/2) 
  fig <- getplotinplotfig(xpos,ypos,"user")
  marinfig <- (mar*c(1,1,-1,-1))[c(2,4,1,3)]*0.2
  boxsize <- fig+marinfig/c(w,w,h,h)
  if(box)
    rect(xlim[[1]],ylim[[1]],xlim[[2]],ylim[[2]],xpd=T)
  if(lines)
    for (i in 1:2)
      for (j in 3:4)
        lines(rbind(c(xlim,ylim)[c(i,j)],boxsize[c(i,j)]))
  if(bg!="")
  {
    rect(boxsize[1],boxsize[3],boxsize[2],boxsize[4],col=bg)
  }
  plotinplot(xpos,ypos,mar)
  plot(1, type="n",xlab="",ylab="",xlim=xlim,ylim=ylim,xaxs="i",yaxs="i",yaxp=signif(c(ylim[1],ylim[2],2),yticks),
       xaxp=signif(c(xlim[1],xlim[2],2),xticks),...)
}

#' sets the focus back to old plot
#'  
#' \code{endzoomplot} sets the focus back to old plot
#' @export
endzoomplot<-function(){
  endplotinplot()
}

#' calculates nested range in list
#'  
#' \code{calcrange} calculates nested range in list
#' @param listofobjects a list of object all containing the nameofnestedelement
#' @param nameofnestedelement name of the nested element which is in each item of listofobjects
#' @export
calcrange<-function(listofobjects,nameofnestedelement="range")
{
  ranges <-do.call("rbind",lapply(listofobjects,FUN=function(x)return(x[[nameofnestedelement]])))
  if (length(ranges)>0)
  {
    ranges <- apply(FUN=range,ranges,MARGIN=2)
  }
  else
  {
    print("error: no ranges")    
  }
  return(ranges)
}
