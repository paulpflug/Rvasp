require(XML)

#' Read a dosdata object
#' 
#' \code{read.dosdata} read a dosdata object form a vasprun.xml file.
#' 
#' @param xmlfile vasprun.xml file
#' @export
read.dosdata <- function(xmlfile){  
  result <- list()
  class(result)<-"dosdata"
  print(paste("processing:",xmlfile))
  ### parsing xml tree 
  xmlobj <- xmlTreeParse(xmlfile, getDTD=FALSE, useInternalNodes=TRUE)
  # find atoms
  atoms <- xpathSApply(xmlobj, "/modeling/atominfo/array[@name='atoms']/set/rc/*",
                       function(x) xmlValue(x))
  result$atoms <- atoms[seq(1,length(atoms),by=2)]
  result$natoms<-length(result$atoms)
  # find fermi energy
  result$efermi <-  xpathSApply(xmlobj, "/modeling/calculation/dos/i",
                                function(x) as.numeric(xmlValue(x)))
  result$ngrid <-  xpathSApply(xmlobj, "/modeling/parameters/separator[@name='dos']/i[@name='NEDOS']",
                               function(x) as.numeric(xmlValue(x)))
  free(xmlobj)
  ### get energie and banddata matrices
  total <-matrix(0,result$ngrid,3)
  partial <- matrix(0,result$natoms*result$ngrid,10)
  counter <-1
  counter2<-1

  state <- c(F,F,F,F,F,F)  
  xmlEventParse(xmlfile,useTagName=F,addContext=F,handlers=list(
    startElement=function(name,content,...){
      if(name=="set" & length(content)>0)
      {
        if(content=="spin1" | content=="spin 1")
        {
          state[4]<<-T
        }
      }
      if(state[1])
      {
        if(state[2] | state[3])
        {
          if(state[2] & state[4] & name=="r") state[5]<<- T  else {}
          if(state[3] & state[4] & name=="r") state[6]<<- T  else {}
        }
        else
        {
          if(name=="total") state[2]<<-T   else{}
          if(name=="partial") state[3]<<-T   else{}
        }
      }
      else
        if(name=="dos") state[1]<<-T   else{}   
    },
    text=function(name,...){
      if(state[5])
      {
        total[counter2,]<<-as.numeric(strsplit(name," +")[[1]])
        counter2<<-counter2+1
        state[5]<<-F
      }
      if(state[6])
      {
        partial[counter,]<<-as.numeric(strsplit(name," +")[[1]])
        counter<<-counter+1
        state[6]<<-F
      }
    },
    endElement=function(name,...){
      if(name=="set")
        state[4]<<-F
      if(state[1])
      {
        if(state[2] & name=="total") state[2]<<-F else{}
        if(state[3] & name=="partial")  state[3]<<-F else{}
      }
      else
        if(name=="dos")  state[1]<<-F else{}
    }))
  # preparation for sorting data
  efermi <- result$efermi
  
  
  natoms <- result$natoms
  nenergies <- result$ngrid
  atoms <- result$atoms  
  atomindex <- as.vector(sapply(1:natoms,rep,nenergies))
  names <- c("energy","s","py","pz","px","dxy","dyz","dz2", "dxz", "dx2","atomnr","atomtype")
  
  partial[,1]<-partial[,1]-efermi
  result$range <- apply(FUN=range,partial,MARGIN=2)
  result$partial <- data.frame(partial,atomindex,atoms[atomindex])
  colnames(result$partial)<-names
  result$total <- data.frame(total[,1:2])
  colnames(result$total) <- c("energy","total")
  result$total[,1] <- result$total[,1]-efermi
  return (result)
}

#' Print a dosdata object
#' 
#' \code{print.dosdata} print a dosdata object.
#' 
#' @param dosdata object of type dosdata
#' @export
print.dosdata <- function(dosdata,...){
  for (name in names(dosdata))
  {
    if (length(dosdata[[name]])>1)
    {
      cat(paste0("   ",name," length:",length(dosdata[[name]])," (",class(dosdata[[name]]),")","\n"))
      
    }
    else
    {
      cat(paste0("   ",name," ",dosdata[[name]],"\n"))
    }
  }
}

#' Add smearing to a dosdata object
#' 
#' \code{dosdata.addsmearing} Add gaussian smearing to a dosdata object.
#' 
#' @param dosdata object of type dosdata
#' @param sigma standard deviation
#' @export
dosdata.addsmearing<-function(dosdata,sigma=0.2)
{  
  dosdata$total[[2]] <-  rowSums(
    apply(dosdata$total,FUN=function(x)
    {
      d <- dnorm(dosdata$total[[1]],mean=x[[1]],sd=sigma)
      d <- d / sum(d)*x[[2]]
    }
          ,MARGIN=1
    )
  )
  
  for (ion in 1:dosdata$natoms)
  {        
    for (orbital in (1:9)+1)
    {
      d <- dosdata$partial[dosdata$partial$atomnr==ion,c(1,orbital)]
      dosdata$partial[dosdata$partial$atomnr==ion,orbital]<- rowSums(
        apply(d,FUN=function(x)
        {
          d <- dnorm(d[[1]],mean=x[[1]],sd=sigma)
          d <- d / sum(d)*x[[2]]
        }
              ,MARGIN=1
        )
      )
    }
  }
  dosdata$range <- apply(FUN=range,dosdata$partial,MARGIN=2)
  return(dosdata)
}

#' Add smearing to a dosvector
#' 
#' \code{dosvector.calcsmearing} add gaussian smearing to a dosvector.
#' 
#' @param energy energy vector
#' @param dos dos vector
#' @param sigma standard deviation
#' @export
dosvector.calcsmearing<-function(energy,dos,sigma=0.1)
{
  matrix <- sapply(energy,dnorm,energy,sigma) 
  return (rowSums(sweep(matrix,2,colSums(matrix)/dos,"/")))
}

#' Creates a plot for dosdata
#' 
#' \code{plot.dosdata} creates a plot for dosdata.
#' Will return object of type dosdata where plotting parameters are saved.
#' For adding curves see \code{\link{plot.dosdata.add}}
#' 
#' @param dosdata object of type dosdata
#' @param smearing if greater zero will add a gaussian smearing with \code{smearing} as standard deviation
#' @param flip will exchange x and y
#' @param norm will divide density by atom count
#' @param fermi will draw fermi level
#' @param col.fermi color of fermi level
#' @param ... further plotting parameters
#' @export
plot.dosdata <- function(dosdata,smearing=0,flip=F,norm=F,fermi=F,col.fermi="blue",xaxs="i",yaxs="i",xlab="Energy [eV]",ylab="Electron density",...)
{
  data <- dosdata
  if (smearing>0)
    data$total[,2] <- dosvector.calcsmearing(data$total[,1],data$total[,2],sigma=smearing)
  dosdata$flip<-flip
  dosdata$norm<-norm
  energy <- data$total[,1]
  density <- data$total[,2]
  if (norm)
  {
    density <- density/data$natoms
    if(ylab == "Electron density")
      ylab <- paste(ylab,"per atom")
  }
  
  if (flip)
  {
    oldmar <- par("mar")
    par(mar=oldmar[c(1,4,3,2)])
    plot(density,energy,type="n",xaxs=xaxs,yaxs=yaxs,yaxt="n",xlab=ylab,...)
    axis(4)
    mtext(xlab,4,3)
  }
  else
    plot(energy,density,type="n",xaxs=xaxs,yaxs=yaxs,xlab=xlab,ylab=ylab,...)
  if (fermi)
  {
    plot.dosdata.addfermi(dosdata,col=col.fermi,...)
  }
  dosdata$lastpolygon <- rep(0,length(density))
  return(dosdata)
}

#' Will add dosdata to existing plot.
#' 
#' \code{plot.dosdata.add} will add dosdata to existing plot.
#' Needs partial dos. Offers two modes.
#' Will return object of type dosdata where plotting parameters are saved.
#' 
#' @param dosdata object of type dosdata
#' @param factor is multiplyed with density
#' @param smearing if greater zero will add a gaussian smearing with \code{smearing} as standard deviation
#' @param orbitals vector of indices of orbitals which are seperatly plotted. Use \code{all} to sum all orbitals
#' @param atomindices vector of indices of atoms which are summed
#' @param type \code{line} or \code{polygon} mode
#' @param ... further plotting parameters
#' @export
plot.dosdata.add <-function(dosdata,factor=1,smearing=0,orbitals=NULL,atomindices=NULL,border=c(NA),col=c("grey"),type=c("line","polygon"),...)
{
  type <- match.arg(type)
  flip <- dosdata$flip
  if (is.null(flip))
    flip<-F
  norm <- dosdata$norm
  if (is.null(norm))
    norm<-F
  if(norm & type=="polygon")
  {
    warning("polygon with normalized dos makes no sense: switching to line mode")
    type <- "line"
  }
  data <- dosdata
  energy <- data$total[,1]
  density <- data.frame(data$total[,2])
  if (norm)
    density <- density/data$natoms
  total <- T
  allorbitals <- F
  if(length(border)<length(orbitals))
    border <- rep(border,ceiling(length(orbitals)/length(border)))
  if(length(col)<length(orbitals))
    col <- rep(col,ceiling(length(orbitals)/length(col)))
  if(!is.null(orbitals)|!is.null(atomindices))
  {
    tmpdensity <- data$partial[,1:11]
    total<-F
    if (!is.null(atomindices))
    {
      #print(data)
      tmpdata <-sapply(energy,FUN=function(x)
      {
        
        selector <- tmpdensity[,1]==x & tmpdensity[,11]%in%atomindices
        tmp <- tmpdensity[selector,2:10]
        return(colSums(tmp))
      })
      tmpdensity <- t(tmpdata)
      if(norm)
        tmpdensity <-tmpdensity/length(atomindices)
    }
    else { 
      tmpdata <-sapply(energy,FUN=function(x)
      {
        
        selector <- tmpdensity[,1]==x
        tmp <- tmpdensity[selector,2:10]
        return(colSums(tmp))
      })
      tmpdensity <- t(tmpdata)
      if(norm)
        tmpdensity <-tmpdensity/data$natoms
      }
    if (!is.null(orbitals))
    {
      if ("all" %in% orbitals)
      {  
        allorbitals <- T
        orbitals <- as.numeric(orbitals[-which(orbitals=="all")])
        density <- data.frame(cbind(tmpdensity[,orbitals],rowSums(tmpdensity)))
      }
      else
      {
        density <- data.frame(tmpdensity[,orbitals])
      }
    }
    else
    {
      density <- data.frame(rowSums(tmpdensity))
    }
  }
  density <- density*factor
  for (i in 1:ncol(density))
  {
    if (smearing>0)
      density[,i] <- dosvector.calcsmearing(energy,density[,i],sigma=smearing)
    if (type=="line")
      {if (flip)
        lines(density[,i],energy,col=col[[i]],...)
    else
      lines(energy,density[,i],col=col[[i]],...)
    }
    if(type=="polygon")
    {
      lastdensity<- dosdata$lastpolygon
      if(is.null(lastdensity))
        lastdensity <- rep(0,length(density[,i]))
      data <- density[,i]
      if(!total)
        {data <- data+lastdensity}
      x <- c(energy,rev(energy))
      y <- c(lastdensity,rev(data))
      if (flip)
      {
        h <- x
        x <- y
        y <- h
      }
      polygon(x,y ,border=border[[i]],col=col[[i]],...)
      dosdata$lastpolygon <- data
    }
  }
  return(dosdata)
}


#' Will add fermi level to existing plot.
#' 
#' \code{plot.dosdata.addfermi } will add fermi level to existing plot.
#' 
#' @param dosdata object of type dosdata
#' @param ... further plotting parameters
#' @export
plot.dosdata.addfermi <- function(dosdata,col="blue",lty=3,...)
{
  flip <- dosdata$flip
  if (is.null(flip))
    flip<-F
  rng<- par("usr")[3:4]
  if (flip)
  {
  abline(h=0,col=col,lty=lty,...)
  text(rng[[2]]-0.1*abs(rng[[2]]-rng[[1]]),0,labels="fermi",col=fermicol,pos=2,...)
  }
  else{
    abline(v=0,col=col,lty=lty,...)
    text(0,rng[[2]]-0.1*abs(rng[[2]]-rng[[1]]),labels="fermi",col=fermicol,pos=2,...)
  }
  
}
