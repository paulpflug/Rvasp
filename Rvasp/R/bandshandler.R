require(XML)
require(lattice)
require(grid)
require(snowfall)
require(akima)
require(rgl)
#' Reads bandsdata
#' 
#' \code{read.bandsdata} Reads bandsdata from vasprun.xml.
#' 
#' @param xmlfile vasprun.xml file
#' @export
read.bandsdata <- function(xmlfile){  
  result <- list()
  class(result)<-"bandsdata"
  print(paste("processing:",xmlfile))
  ### parsing xml tree 
  xmlobj <- xmlTreeParse(xmlfile, getDTD=FALSE, useInternalNodes=TRUE)
  # find k basis in xml
  result$kbasis <- t(xpathSApply(xmlobj, "/modeling/structure[@name='initialpos']/crystal/varray[@name='rec_basis']/*",
                               function(x) as.numeric(strsplit(xmlValue(x)," +")[[1]][-1])))
  # find k points in xml
  result$kpoints <- t(xpathSApply(xmlobj, "/modeling/kpoints/varray[@name='kpointlist']/*",
                                  function(x) as.numeric(strsplit(xmlValue(x)," +")[[1]][-1])))
  result$kweights <- t(xpathSApply(xmlobj, "/modeling/kpoints/varray[@name='weights']/*",
                                  function(x) as.numeric(strsplit(xmlValue(x)," +")[[1]][-1])))
  result$nkpoints<-nrow(result$kpoints)
  # calculate high symetric points and project k points on 1-D distance function
  result$sympoints<-c(0)
  result$sympointindices<-c(1)
  result$sympointindicesrev<-c(1)
  result$kpointdistances <- c(0)
  if (result$nkpoints>1)
  {
    for (i in 1:(result$nkpoints-1))
    {
      point1 <- result$kpoints[i,]%*%t(result$kbasis)
      point2 <- result$kpoints[i+1,]%*%t(result$kbasis)
      d <- dist(rbind(point1,point2))[[1]]*2*pi
      if(d == 0 | i == (result$nkpoints-1))
      {
        result$sympoints<-c(result$sympoints,result$kpointdistances[[i]]+d)
        if(d == 0)
        {
        result$sympointindices<-c(result$sympointindices,i)
        result$sympointindicesrev<-c(result$sympointindicesrev,i+1)
        }
      }
      result$kpointdistances <- c(result$kpointdistances,result$kpointdistances[[i]]+d)
    }
  }
  result$sympointindices<-c(result$sympointindices,result$nkpoints)
  result$sympointindicesrev<-c(result$sympointindicesrev,result$nkpoints)
  result$kpointsflat <- result$kpointdistances/max(result$kpointdistances)
  result$sympoints <- result$sympoints/max(result$kpointdistances)
  # find atoms
  atoms <- xpathSApply(xmlobj, "/modeling/atominfo/array[@name='atoms']/set/rc/*",
                       function(x) xmlValue(x))
  result$atoms <- atoms[seq(1,length(atoms),by=2)]
  result$natoms<-length(result$atoms)
  # find fermi energy
  result$efermi <-  xpathSApply(xmlobj, "/modeling/calculation/dos/i",
                                function(x) as.numeric(xmlValue(x)))
  # find band count
  result$nbands <- xpathSApply(xmlobj, "/modeling/parameters/separator[@name='electronic']/i[@name='NBANDS']",
                               function(x) as.numeric(xmlValue(x)))
  free(xmlobj)
  gc()
  ### get energie and banddata matrices
  maxproj <- result$natoms*result$nbands*result$nkpoints
  maxeigen <- result$nbands*result$nkpoints
  maxprojpotenz <- 10^floor(log10(maxproj))
  maxeigenpotenz <- 10^floor(log10(maxeigen))
  energies2 <-matrix(0,maxeigen,1)
  bands2 <- matrix(0,maxproj,9)
  counter <-1
  counter2<-1
  setcount <- 0
  state <- c(F,F,F,F,F)  
  xmlEventParse(xmlfile,useTagName=F,addContext=F,handlers=list(
    startElement=function(name,content,...){                
                if(name=="projected")
                {
                  print("reading projected values")
                  print(paste("total: ",maxproj,sep=""))
                  state[1] <<- T
                }
                else
                {
                  if(name=="eigenvalues")
                  {
                    if(!state[1])
                      {print("reading eigenvalues")
                    print(paste("total: ",maxeigen,sep=""))
                    }
                    state[2]<<-T
                  }
                  else
                  {
                    if(name=="set" & length(content)>0)
                    {
                      if(content=="spin1" | content=="spin 1")
                      {
                        state[3]<<-T
                      }
                      if(state[3])
                        setcount<<-setcount+1
                    }
                    else
                    {
                      if(name=="r"&state[3])
                      {
                        if (state[1]&!state[2])
                          state[4]<<-T
                        if(!state[1]&state[2])
                          state[5]<<-T
                      }
                    }
                  }
                }
              },
    text=function(name,...){
                if(state[5])
                {
                  
                  if(counter2%%maxeigenpotenz==0)
                  {
                    print(paste(counter2,maxeigen,sep="/"))
                  }
                  energies2[counter2,]<<-as.numeric(strsplit(name," +")[[1]][[1]])
                  counter2<<-counter2+1
                }
                if(state[4])
                {
                  if(counter%%maxprojpotenz==0)
                  {
                    print(paste(counter,maxproj,sep="/"))
                  }
                  bands2[counter,]<<-as.numeric(strsplit(name," +")[[1]])                  
                  counter<<-counter+1
                }
              },
    endElement=function(name,...){
      if(name=="projected")
      {
        state[1] <<- F
      }
      else
      {
        if(name=="eigenvalues")
        {
          state[2]<<-F
        }
        else
        {
          if(name=="set")
          {
            if(state[3])
              setcount<<-setcount-1
            if(setcount==0)
            {
              state[3]<<-F
            }
          }
          else
          {
            if(name=="r")
            {
                state[4]<<-F
                state[5]<<-F
            }
          }
        }
      }
              }))
  bands<-bands2
  energies<-energies2
  # preparation for sorting data
  natoms <- result$natoms
  nbands <- result$nbands
  nkpoints <- result$nkpoints
  nenergies <- length(energies)/nbands
  atoms <- result$atoms
  efermi <- result$efermi
  deltakpunkt<- natoms*nbands
  deltaband <- natoms
  names <- c("kpoint","energy","s","py","pz","px","dxy","dyz","dz2", "dxz", "dx2","atomnr","atomtype")
  ### sorted atomnumbers
  index <- as.vector(sapply(1:natoms,rep,nenergies))
  ##########################################
  # K1B1A1 -> 1
  # K1B1A2 -> NK+1
  # K1B1A3 -> 2*NK+2
  # ...
  # K1B2A1 -> 1+NA
  # K1B2A2 -> NK+1+NA
  # ...
  # K2B1A1 -> 2
  # K2B1A2 -> NK+2
  ### spacing between
  dk <- (1:(nkpoints)-1)*deltakpunkt
  sortindex <- rep(dk,natoms)+as.vector(sapply(1:natoms-1,rep,nkpoints))
  #print(sortindex)

 # print(ksortindex)
  #dataindex <- 
  #print(dataindex)
  
  energieindex <- seq(1,length(energies),by=nbands)
  # sort data
  #sfInit(parallel=TRUE,cpus=4,type="SOCK")
  result$bands <- lapply( 0:(nbands-1),FUN=function(band)
  {
    banddata <- list()
    class(banddata) <- "banddata"
    banddata$bandnr <- band+1
    data <- data.frame(result$kpointsflat,
                       energies[band+energieindex]-efermi
    )
    if(counter>1)
    {
    projecteddata <- data.frame(
      rep(result$kpointsflat,natoms), # kpoints
      rep(energies[band+energieindex]-efermi,natoms), # Energy
      bands[sortindex+band*natoms+1,], #s py .. dx2
      index,        #atomnr
      atoms[index]          #atomtype
    )
    
    names(projecteddata)<-names
    banddata$data <-   projecteddata    
    banddata$range <- apply(FUN=range,banddata$data[-13],MARGIN=2)
    }
    names(data)<-names[1:2]
    
    banddata$simpledata <- data
    
    return(banddata)  
  })  
  names(result$bands)<- paste0("band",1:nbands)
  if(counter>1)
  {
  result$range <- calcrange(result$bands)  
  }
  return (result)
}

#' Custom print for object of class bandsdata
#' 
#' \code{print.bandsdata} custom print for object of class bandsdata.
#' 
#' @param bandsdata object of class bandsdata
#' @export
print.bandsdata <- function(bandsdata,...){
  for (name in names(bandsdata))
  {
    if (length(bandsdata[[name]])>1)
    {
      if (class(bandsdata[[name]])=="banddata")
        print(bandsdata[[name]])
      else
      cat(paste0("   ",name," length:",length(bandsdata[[name]])," (",class(bandsdata[[name]]),")","\n"))
      
    }
    else
    {
      cat(paste0("   ",name," ",bandsdata[[name]],"\n"))
    }
  }
}

#' Custom print for object of class banddata
#' 
#' \code{print.banddata} custom print for object of class banddata.
#' 
#' @param banddata object of class banddata
#' @export
print.banddata <- function(banddata,...){
  for (name in names(banddata))
  {
    if (length(banddata[[name]])>1)
    {
      cat(paste0("   ",name," length:",length(banddata[[name]])," (",class(banddata[[name]]),")","\n"))
      
    }
    else
    {
      cat(paste0("   ",name," ",banddata[[name]],"\n"))
    }
  }
}

#' Adds a path of high symmetry points to bandsdata
#' 
#' \code{bandsdata.calcsympointpath} adds a custom path of high symmetry points to bandsdata.
#'  ex. sympointpath=list(c(1,2),c(3,4)) 
#'  goes from 1 to 2 and 3 to 4
#' returns class of type bandsdata
#' 
#' @param bandsdata object of class bandsdata
#' @param sympointpath list of vectors containing indices of high symmetry points
#' @export
bandsdata.calcsympointpath<-function(bandsdata,sympointpath){
  print("calculating sympoints")
  xmax <- 0
  sympath <- list()
  sympath$sympoints<- rep(-1,length(bandsdata$sympoints))
  sympath$data <- do.call(rbind, lapply(sympointpath,FUN=function(x)
    {
    i1 <- bandsdata$sympointindicesrev[x[which.min(x)]]
    i2 <- bandsdata$sympointindices[x[which.max(x)]]    
    i <- c(i1,i2) 
    xrng <- range(bandsdata$kpointsflat[i])
    xrng <- xrng-min(xrng)+xmax
    xmax <<- xrng[2]
    i<- i[order(x)]
    sympath$sympoints[x] <<-xrng
    indices <- i[1]:i[2]
    return(cbind(indices,seq(xrng[1],xrng[2],along.with=indices),bandsdata$kpointsflat[indices]))
  }))
  class(sympath)<-"sympath"
  bandsdata$sympath <- sympath
  return(bandsdata)
}

#' Adds custom high symmetry points to bandsdata
#' 
#' \code{bandsdata.addsympoint} adds custom high symmetry points to bandsdata.
#' returns class of type bandsdata
#' 
#' @param bandsdata object of class bandsdata
#' @param indices of high symmetry points in \code{bandsdata$kpoints}
#' @export
bandsdata.addsympoint <- function(bandsdata,indices){
  bandsdata$sympoints <-sort(c(bandsdata$sympoints,bandsdata$kpointsflat[indices]))
  bandsdata$sympointindices <- sort(c(bandsdata$sympointindices,indices))
  bandsdata$sympointindicesrev <- sort(c(bandsdata$sympointindicesrev,indices))
  print("added sympoints to bandsdata")
  print(paste("new sympointindices:",paste(bandsdata$sympointindices,collapse=" ")))
  return(bandsdata)
}

#' Plots bandsdata
#' 
#' \code{plot.bandsdata} plots bandsdata.
#' returns class of type bandsdata
#' 
#' @param bandsdata object of class bandsdata
#' @param bands limits plotting to specified bands
#' @param col.bands color of bands
#' @param sympointpath calls \code{\link{bandsdata.calcsympointpath}}
#' @param fermi adds line at Fermi level
#' @param fermi.col color of line at Fermi level
#' @param sym.labels adds labels at high symmetry points
#' @param sym.col.labels color of high symmetry point labes
#' @param sym.lty line typ of high symmetry point lines
#' @param sym.col.lines color of high symmetry point lines
#' @param energyoffset will be added to energy of all bands
#' @param ... further plotting parameters
#' @export
plot.bandsdata <- function(bandsdata
                           ,bands=1:length(bandsdata$bands)
                           ,sympointpath=NULL
                           ,col.bands="black"
                           ,type="l"
                           ,fermi=F
                           ,fermi.col="blue"
                           ,sym.labels=NULL
                           ,sym.col.labels="black"
                           ,sym.lty=3
                           ,sym.col.lines="red"                           
                           ,xlim=NULL
                           ,xaxs="i"
                           ,yaxs="i"
                           ,xlab="Wave vector"
                           ,ylab="Energy (eV)"
                           ,energyoffset=0
                           ,...){
  print("plotting bandsdata")
  rng <- bandsdata$range
  if(!is.null(sympointpath))
  {
    bandsdata <- bandsdata.calcsympointpath(bandsdata,sympointpath)    
  }
  bandsdata$energyoffset <- energyoffset  
  rng[1:2,2]<- rng[1:2,2]+energyoffset
  if(!is.null(bandsdata$sympath))
  {
    rng[1:2,1]<-range(bandsdata$sympath$data[,2])
  }
  if (is.null(xlim))
    xlim <- rng[1:2,1]
  plot(rng[,1:2],type="n",xaxt="n",xlim=xlim,xaxs=xaxs,yaxs=yaxs, xlab=xlab, ylab=ylab,...)
  if (!type=="n") plot.bandsdata.addbands(bandsdata,bands,col=col.bands,...)
  if (fermi)
  {
    plot.bandsdata.addfermi(bandsdata,xlim=xlim,col=fermi.col,...)
  }
  if(!is.null(sym.labels))
  {
    plot.bandsdata.addsymnnames(bandsdata,labels=sym.labels,col.labels=sym.col.labels,
                                lty=sym.lty,col.lines=sym.col.lines)
  }
  return(bandsdata)
}

#' Adds high symmetry point labels to existing plot
#' 
#' \code{plot.bandsdata.addsymnnames} adds high symmetry point labels to existing plot
#' 
#' @param bandsdata object of class bandsdata
#' @param labels adds labels at high symmetry points
#' @param col.labels color labels
#' @param lty line typ of high symmetry point lines
#' @param col.lines color of high symmetry point lines
#' @param ... further plotting parameters
#' @export
plot.bandsdata.addsymnnames <- function(bandsdata,labels,col.labels="black",lty=3,col.lines="red",...){
  print("adding symnames")
  sp <- bandsdata$sympoints
  if(!is.null(bandsdata$sympath))
  {
    sp<-bandsdata$sympath$sympoints
  }
  select <- !duplicated(sp)
  abline(v=sp[select],col=col.lines,lty=lty,...)
  axis(1,at=sp[select],labels=labels[select],col=col.labels,...)
} 

#' Adds line at Fermi level to existing plot
#' 
#' \code{plot.bandsdata.addsymnnames} adds line at Fermi level to existing plot
#' 
#' @param bandsdata object of class bandsdata
#' @param col color of line at Fermi level
#' @param lty line typ of line at Fermi level
#' @param ... further plotting parameters
#' @export
plot.bandsdata.addfermi<-function(bandsdata,col="blue",lty=3,...){
  abline(h=0,col=col,lty=lty,...)
}

#' Adds bands to existing plot
#' 
#' \code{plot.bandsdata.addsymnnames} adds bands to existing plot
#' 
#' @param bandsdata object of class bandsdata
#' @param bands limits plotting to specified bands
#' @param energyoffset will be added to energy of all bands
#' @param ... further plotting parameters
#' @export
plot.bandsdata.addbands <- function(bandsdata,bands=1:length(bandsdata$bands),energyoffset=NULL,...){
  if(is.null(energyoffset))
  {
    if(!is.null(bandsdata$energyoffset))
    {
      energyoffset <- bandsdata$energyoffset
    }
    else
    {
      energyoffset<-0
    }
  }
  print("plotting bands")
  k <- bandsdata$kpointsflat
  i <- 1:bandsdata$nkpoints
  if(!is.null(bandsdata$sympath))
  {
    k <- bandsdata$sympath$data[,2]
    i <- bandsdata$sympath$data[,1]
  }
  
  lapply(bandsdata$bands,FUN=function(band)
  {
    if (band$bandnr %in% bands)
    {
    lines(cbind(k,band$simpledata[i,2]+energyoffset),...)
    }
  })
}

#' Adds number of all bands to existing plot
#' 
#' \code{plot.bandsdata.addsymnnames} adds number of all bands to existing plot.
#' 
#' @param bandsdata object of class bandsdata
#' @param kpoints at which numbers will be plotted, is recycled
#' @param energyoffset will be added to energy of all bands
#' @param ... further plotting parameters
#' @export
plot.bandsdata.addnumbers <- function(bandsdata
                                      ,kpoints=NULL
                                      ,energyoffset=NULL
                                      ,...){
  if(is.null(energyoffset))
  {
    if(!is.null(bandsdata$energyoffset))
    {
      energyoffset <- bandsdata$energyoffset
    }
    else
    {
      energyoffset<-0
    }
  }
  print("plotting numbers")
  k <- 1:bandsdata$nkpoints
  kpointsflat <- bandsdata$kpointsflat
  if(!is.null(bandsdata$sympath))
  {
    k <- bandsdata$sympath$data[,1]
    kpointsflat <- bandsdata$sympath$data[,2]
  }
  if(is.null(kpoints))
  {
  kpoints <-if(length(k)>25)
    seq(5,length(k)-5,by=10)+min(k)-1 else k
  }
  if(length(kpoints)<length(bandsdata$bands))
    kpoints<-rep(kpoints,length.out=length(bandsdata$bands))
  lapply(1:length(bandsdata$bands),FUN=function(i){
    text(kpointsflat[k==kpoints[[i]]],bandsdata$bands[[i]]$simpledata[kpoints[[i]],2]+energyoffset,labels=i,...)
  })
  return(0)
}

#' Will give a object of class bulkbands
#' 
#' \code{calculation.getbulkbands} will give a object of class bulkbands
#' based on calculation object.
#' 
#' @param calculation object of class calculation
#' @export
calculation.getbulkbands<-function(calculation,...){
  print("calculating bulk polygons")
  kpointsbulk <-calculation[[1]]$banddata$kpointsflat
  bandcountbulk <- calculation[[1]]$banddata$nbands
  poly<-list()
  for(i in 1:bandcountbulk)
  {
    energies <- numeric()
    for(j in 1:length(calculation))
  {
    d <- calculation[[j]]$banddata$bands[[i]]$simpledata
    energies <- cbind(energies,d[,2])      
  }
  energies <- apply(FUN=range,energies,MARGIN=1)
  poly[[i]]<-list()
  poly[[i]]$x <- c(kpointsbulk,rev(kpointsbulk))
  poly[[i]]$y <- c(energies[1,],rev(as.vector(energies[2,])))
  }
  bb <- list()
  class(bb)<-"bulkbands"
  bb$poly <- poly
  return(bb)
}

#' Adds bulkband polygons to existing plot
#' 
#' \code{plot.bulkbands.add} adds bulkband polygons to existing plot.
#' 
#' @param bulkbands object of class bulkbands
#' @param ... further plotting parameters
#' @export
plot.bulkbands.add<-function(bulkbands,col="grey",...){
  print("plotting bulk polygons")
  for(p in bulkbands$poly)
  {
    polygon(p$x,p$y,col=col,border=NA)
  }
}

#' Will give a object of class projectedbands
#' 
#' \code{bandsdata.getprojecteddata} will give a object of class projectedbands
#' 
#' 
#' @param bandsdata object of class bandsdata
#' @param atomindices indices of atoms over which will be summed
#' @param bands which will be included
#' @param energyintervall in which bands will be included (overwrites \code{bands})
#' @export
bandsdata.getprojecteddata <- function(bandsdata,
                                       atomindices=1:bandsdata$natoms,
                                       bands=1:bandsdata$nbands,
                                       energyintervall=NULL,
                                       cpus=1){
  print("calculating projecteddata")
  print("sum over atoms")
  nkpoints <- length(bandsdata$kpointsflat)
  if(!is.null(energyintervall))
  {
    energyintervall <- range(energyintervall)
    bands <- sapply(bandsdata$bands,FUN=function(band){
      return(any(energyintervall[[1]]<band$simpledata$energy& band$simpledata$energy<energyintervall[[2]]))
    })
  }
  d <- list()
  sumoveratoms <- function(band)
  {
    do.call(rbind,lapply(1:nkpoints,FUN=function(i){
      return(c(band$simpledata[i,],colSums(band$data[i+(atomindices-1)*nkpoints,3:11])))
    }))
  }
  if (cpus>1 & length(bands) >cpus*10)
  {
    sfInit(parallel=TRUE, cpus=cpus, type="SOCK")
    sfExport("atomindices","nkpoints")
    d$bands <- sfLapply(bandsdata$bands[bands],fun=sumoveratoms)
    sfRemoveAll()
    sfStop()   
  }
  else
  {
    print("on single cpu")
    d$bands <- lapply(bandsdata$bands[bands],FUN=sumoveratoms)
  }
  if(!is.null(bandsdata$sympath))
    d$sympath <- bandsdata$sympath
  if(!is.null(bandsdata$energyoffset))
    d$energyoffset <- bandsdata$energyoffset
  class(d)<-"projectedbands"
  return (d)
}

#' Adds projected orbitals to existing plot
#' 
#' \code{plot.projectedbands.add} adds projected orbitals to existing plot.
#' 
#' @param projectedbands object of class projectedbands
#' @param bands limits plotting to specified bands
#' @param orbitals list of orbitals to plot. To sum over orbitals 2 and 3 use \code{list(1,c(2,3),4)}
#' @param col colors for orbitals, will be recycled by the length of orbitals, if not provided, will use col.platte for generation.
#' @param col.palette color palette for orbitals
#' @param legend position of legend, \code{NULL} will supress plotting
#' @param legendcex size of legend
#' @param energyoffset will be added to energy of all bands
#' @param ... further plotting parameters
#' @export
plot.projectedbands.add <- function(projectedbands
                                    ,bands=1:length(projectedbands$bands)
                                    ,orbitals=list(1,2,3,4)
                                    ,col=NULL
                                    ,col.palette=colorRampPalette(c("red","blue","green"))
                                    ,pch=15:(14+length(orbitals))
                                    ,cex=0.8
                                    ,legend="topright"
                                    ,legendcex=0.8
                                    ,energyoffset=NULL
                                    ,...){
  if(is.null(energyoffset))
  {
    if(!is.null(projectedbands$energyoffset))
    {
      energyoffset <- projectedbands$energyoffset
    }
    else
    {
      energyoffset<-0
    }
  }
  print("plotting projecteddata")
  k <- projectedbands$bands[[1]][,1]
  i <- 1:nrow(projectedbands$bands[[1]])
  if(!is.null(projectedbands$sympath))
  {
    k <- projectedbands$sympath$data[,2]
    i <- projectedbands$sympath$data[,1]
  }
  projbands <- lapply(projectedbands$bands,FUN=function(band){
    band <-apply(band[,-(1:2)],MARGIN=2,FUN=as.numeric)
    do.call(cbind,lapply(orbitals,FUN=function(orb){      
      if(length(orb)>1)
        b <- rowSums(band[,(orb)])
      else
        b <- as.numeric(band[,(orb)])
      return (b)
    }))
  })
  maxis <- apply(do.call(rbind,projbands),2,max)
  if(is.null(col)){
  col <- col.palette(length(orbitals))
  }else{
    col <- rep(col,length.out=length(orbitals))
  }
  print("plotting")
  lapply(bands,FUN=function(band)
  {
    lapply(1:length(orbitals),FUN=function(orb){
      da <- projbands[[band]][,orb]
      if(maxis[orb]>0)
      {
        da <- da/maxis[orb]
      }
      #da <- da^(3/4)
      trns <- round(da*255)
     points(cbind(k,simplify2array(projectedbands$bands[[band]][i,2])+energyoffset),col=makeTransparent(col[orb],trns)[i],cex=(da*cex)[i],pch=pch[orb],...)
    })
  })
  names<-colnames(projectedbands$bands[[1]][,-c(1:2)])
  n1 <- substring(names,1,1)
  n2<- substring(names,2)
  rawnames <- paste0(n1,"[",n2,"]")
  names <- sapply(orbitals,FUN=function(orb)
    {
    if (length(unique(n1[orb]))==1&length(orb)!=1)
    {
      return(paste0(n1[orb[[1]]],"[",paste0(n2[orb],collapse=""),"]"))
    }
    else
    {
     return(paste(rawnames[orb],collapse="+"))
    }
  })
  if(!is.null(legend))
    legend(legend,legend=parse(text=names),bg="white",col=col[1:length(orbitals)],pch=pch,cex=legendcex,seg.len=1,x.intersp=0.5,y.intersp=0.8,...)

}

#' Gives a intervall around a high symmetry point
#' 
#' \code{bandsdata.getintervallaroundsympoint} gives a intervall around a high symmetry point.
#' uses flat (2d) kpoints.
#' 
#' @param bandsdata object of class bandsdata
#' @param sympointnumber number of high symmetry point to use
#' @param intervall in per cent
#' @export
bandsdata.getintervallaroundsympoint<-function(bandsdata,sympointnumber,intervall=0.2){
  rng <- range(bandsdata$sympoints)
  #print(rng)
  sympoint <- bandsdata$sympoints[[sympointnumber]]
  #print(sympoint)
  if (sympoint==rng[[1]]){
    min<-sympoint
  }else{
    min <- sympoint-(sympoint-bandsdata$sympoints[[sympointnumber-1]])*intervall
  }
  if (sympoint==rng[[2]]){
    max<-sympoint
  }else{
    
    max <- sympoint-(sympoint-bandsdata$sympoints[[sympointnumber+1]])*intervall
  }
  return(c(min,max))
}

#' Gives the distance between to bands at one kpoint
#' 
#' \code{bandsdata.getbanddistance} gives the distance between to bands at one kpoint.
#' 
#' @param bandsdata object of class bandsdata
#' @param kpoint index of kpoint to use
#' @param bands vector containing two indices of the bands to use
#' @export
bandsdata.getbanddistance<-function(bandsdata,kpoint,bands){
  e1 <- bandsdata$bands[[bands[[1]]]]$simpledata[kpoint,2]
  e2 <- bandsdata$bands[[bands[[2]]]]$simpledata[kpoint,2]
  return(abs(e1-e2))
}

#' Gives the energy of a band at one kpoint
#' 
#' \code{bandsdata.getenergy} gives the energy of a band at one kpoint.
#' 
#' @param bandsdata object of class bandsdata
#' @param kpoint index of kpoint to use
#' @param band index of band to use
#' @export
bandsdata.getenergy<-function(bandsdata,kpoint,band){
  return(bandsdata$bands[[band]]$simpledata[kpoint,2])
}

#' Fit a function to a band
#' 
#' \code{bandsdata.fit} fit a function to a band.
#' will return a object of class bandsfit.
#' 
#' @param bandsdata object of class bandsdata
#' @param bandnr index of band to use
#' @param kpoints indices of kpoints to use
#' @param fitname used fitfunction
#' @param startingparamters for use in fitfunction
#' @param constants for use in fitfunction
#' @export
#' @seealso \code{\link{bandsdata.fit.dirac.function}} and \code{\link{bandsdata.fit.quadratic.function}}
bandsdata.fit <- function(bandsdata,
                          bandnr,
                          kpoints,
                          fitname=c("dirac","quadratic"),
                          startingparameters,
                          constants=NULL){
  fitname <- match.arg(fitname)
  if (!is.null(constants)&&!is.null(constants$k0)&&constants$k0>1)
  {
    constants$k0<- bandsdata$kpointdistances[constants$k0]
  }
  objectivem <- function(parameters,fitfunction,data,constants=list())
  {
    err <- sum((data$energy-fitfunction(parameters,data$kpoint,constants))^2)
    return (err)
  }
  fitfunction<-get(paste("bandsdata.fit",fitname,"function",sep="."))
  data <- bandsdata$bands[[bandnr]]$simpledata[kpoints,]
  data$kpoint <- bandsdata$kpointdistances[kpoints]
  if(is.null(constants))
    obj <- optim(startingparameters,objectivem,fitfunction=fitfunction,data=data,method="BFGS")
  else
    obj <- optim(startingparameters,objectivem,fitfunction=fitfunction,data=data,method="BFGS",constants=constants)
  class(obj)<-"bandsfit"
  obj$ap <-constants
  obj$fitname <- fitname
  obj$fitfunction <- fitfunction
  obj$kdist <- bandsdata$kpointdistances
  obj$korig <- bandsdata$bands[[bandnr]]$simpledata$kpoint
  obj$k <-data$kpoint
  obj$knumbers <- kpoints
  if(!is.null(bandsdata$sympath))
    obj$sympath <- bandsdata$sympath
  if(!is.null(bandsdata$energyoffset))
    obj$energyoffset <- bandsdata$energyoffset
  return(obj)
}

#' Creates startingparamters for dirac fit
#' 
#' \code{bandsdata.fit.dirac.makeparameters} creates startingparamters for dirac fit.
#' 
#' @param vF starting Fermi velocity
#' @export
bandsdata.fit.dirac.makeparameters<-function(vF){
  p <- c(vF)
  names(p) <- c("vF")
  return (p)
}

#' Creates constants for dirac fit
#' 
#' \code{bandsdata.fit.dirac.makeconstants} creates constants for dirac fit.
#' 
#' @param bandsdata object of class bandsdata
#' @param bands indices of bands forming upper and lower dirac cone
#' @param k0 index of high symmetry point where the dirac cone appears
#' @param factor \code{-1}/\code{1} for lower/upper dirac cone 
#' @export
bandsdata.fit.dirac.makeconstants<-function(bandsdata,bands,k0=0,factor=1){
  d <- bandsdata.getbanddistance(bandsdata,k0,bands)
  evalence <- bandsdata.getenergy(bandsdata,k0,bands[[1]])
  ap <- list()
  ap$EF<-d/2+evalence
  ap$k0<-k0
  ap$Eg<-d
  ap$factor<-factor
  return (ap)
}

#' Fit function for dirac fit
#' 
#' \code{bandsdata.fit.dirac.function} contains fit function for dirac fit.
#' 
#' @param parameters list containing all paramters
#' @param kpoints indices of kpoints to use
#' @param constants for use in fitfunction
#' @export
bandsdata.fit.dirac.function<-function(parameters,kpoints,constants){
  p<-parameters ## p[[1]] vF
  ap <- constants  
  return(sign(ap$factor)*sqrt((ap$Eg/2)^2+(abs(kpoints-ap$k0)*p[[1]])^2)+ap$EF)
}

#' Creates startingparamters for quadratic fit
#' 
#' \code{bandsdata.fit.quadratic.makeparameters} creates startingparamters for quadratic fit.
#' used in \code{\link{bandsdata.fit}}.
#' 
#' @param m starting mass
#' @export
bandsdata.fit.quadratic.makeparameters<-function(m){
  p <- c(m)
  names(p) <- c("m*")
  return (p)
}

#' Creates constants for quadratic fit
#' 
#' \code{bandsdata.fit.quadratic.makeconstants} creates constants for quadratic fit.
#' 
#' @param bandsdata object of class bandsdata
#' @param band index of band to use
#' @param k0 index of high symmetry point to use
#' @export
bandsdata.fit.quadratic.makeconstants<-function(bandsdata,band,k0=0){
  ap <- list()
  ap$E<-bandsdata.getenergy(bandsdata,k0,band)
  ap$k0<-k0
  return (ap)
}

#' Fit function for quadratic fit
#' 
#' \code{bandsdata.fit.quadratic.function} contains fit function for quadratic fit.
#' used in \code{\link{bandsdata.fit}}.
#' 
#' @param parameters list containing all paramters
#' @param kpoints indices of kpoints to use
#' @param constants for use in fitfunction
#' @export
bandsdata.fit.quadratic.function<- function(parameters,kpoints,constants){
  p<-parameters ## p[[1]] m*
  ap <- constants  
  return(p[[1]]*(kpoints-ap$k0)^2+ap$E)
}

#' Custom predict for object of class bandsfit
#' 
#' \code{predict.bandsfit} custom predict for object of class bandsfit.
#' 
#' @param o object
#' @param x values to use for prediction
#' @export
predict.bandsfit<-function(o,x){
  data <- data.frame(x,o$fitfunction(o$par,x,o$ap))
  names(data) <- c("kpoint","energy")
  return(data)
}

#' Calculation of Fermi velocity
#' 
#' \code{bandsfit.dirac.getv} calculation of Fermi velocity.
#' Will be in 10e6 m/s
#' 
#' @param bandsfit object of class bandsfit, generated by dirac fit
#' @export
bandsfit.dirac.getv<-function(bandsfit){
  p <- bandsfit$par[[1]]
  hbar <- 6.582119*10^-16 ## eVs
  kfac <- 10^-10
  return(p*kfac/hbar)
}

#' Calculation of effectiv mass
#' 
#' \code{bandsfit.quadratic.getm} calculation of effectiv mass.
#' Will be in units of electron mass
#' 
#' @param bandsfit object of class bandsfit, generated by quadratic fit
#' @export
bandsfit.quadratic.getm<-function(bandsfit){
  p <- bandsfit$par[[1]]
  hbar <- 6.582119*10^-16 ## eVs
  kfac <- 10^-10
  me <- 0.510998928*10^3
  c <- 3*10^8
  return(hbar^2*c^4/p/2/me)
} 

#' Custom print for object of class bandsfit
#' 
#' \code{bandsfit.quadratic.getm} custom print for object of class bandsfit.
#' 
#' @param bandsfit object of class bandsfit
#' @export
print.bandsfit<-function(bandsfit,...){
  switch(bandsfit$fitname,
         "dirac"={
            v <- bandsfit.dirac.getv(bandsfit)
            print(paste("v=",signif(v,3),"m/s"))
          },
         "quadratic"={
           m <- bandsfit.quadratic.getm(bandsfit)
           print(paste("m=",signif(m,3),"m_e"))
         },
          {
              print.default(bandsfit)
          })
}

#' Adds a bandsfit to an existing plot
#' 
#' \code{plot.bandsfit.add} adds a bandsfit to an existing plot
#' 
#' @param bandsfit object of class bandsfit
#' @param kpoints indices of kpoints, where the fit is plotted. All in between points will also be used.
#' @param n count of points used for drawing
#' @param energyoffset will be added to energy the fit
#' @export
plot.bandsfit.add<-function(bandsfit,kpoints=bandsfit$knumbers,n=201,energyoffset=NULL,...){
  if(is.null(energyoffset))
  {
    if(!is.null(bandsfit$energyoffset))
    {
      energyoffset <- bandsfit$energyoffset
    }
    else
    {
      energyoffset<-0
    }
  }
  print("plotting bandsfit")
  x <- bandsfit$kdist[range(kpoints)]
  kdist <- seq(x[[1]],x[[2]],length.out=n)
  data <- predict(bandsfit,kdist)
  if(!is.null(bandsfit$sympath))
  {
    x <- range(bandsfit$sympath$data[bandsfit$sympath$data[,1]%in%kpoints,2])
  }
  else
  {
    x <- bandsfit$korig[range(kpoints)]    
  }  
  k <- seq(x[[1]],x[[2]],length.out=n)
  lines(k,data$energy+energyoffset,col="red",...)
}

#' Plots the 2d grid used in bandsdata
#' 
#' \code{plot.bandsdata.grid} plots the 2d grid used in bandsdata.
#' Main purpose is to test sym operations, to come to a statisfying grid.
#' 
#' @param bandsdata object of class bandsdata
#' @param sym See \code{\link{dataframe.applysymoperations}} for usage.
#' @param foldback folding of points back in the BZ can be disabled
#' @param ... further plotting parameters
#' @export
plot.bandsdata.grid <- function(bandsdata,sym=NA,foldback=T,xlim=NA,ylim=NA,...){
  rbase <- bandsdata$kbasis[1:2,1:2]
  vec2d <- reciprocalbasis.getbrillouinzone(bandsdata$kbasis)
  kpoints <- (bandsdata$kpoints[,1:2]%*%rbase)
  colnames(kpoints)<-c("kx","ky")
  if(foldback){    
    kpoints<-brillouinzone.extendkpoints(vec2d,kpoints)
    kpoints <- brillouinzone.selectkpoints(vec2d,kpoints)    
  }
  if(length(sym)>1 || !is.na(sym)){
    kpoints <- dataframe.applysymoperations(kpoints,sym)
  }
  kpoints <- unique(round(kpoints,6))
  if(any(is.na(xlim))) xlim <- range(vec2d[,1])
  if(any(is.na(ylim))) ylim <- range(vec2d[,2])
  plot(kpoints,asp=1,xlim=xlim,ylim=ylim,...)
  polygon(vec2d)
}

#' Plots the contour bandsdata
#' 
#' \code{plot.bandsdata.contour} plots the contour bandsdata of a single band.
#' Projected states can be added as additional datalayer.
#' First make sure, you have a statisfying grid. See \code{\link{plot.bandsdata.grid}}.
#' 
#' @param bandsdata object of class bandsdata
#' @param band index of band to plot
#' @param n resolution of datalayers
#' @param sym See \code{\link{dataframe.applysymoperations}} for usage and \code{\link{plot.bandsdata.grid}} for testing.
#' @param breaks number of steps in color coding
#' @param colorpalette colors which are used if \code{projected} is \code{False}
#' @param zone see \code{\link{reciprocalbasis.getbrillouinzone}} for usage.
#' @param linear see \code{\link{akima::interp}} for usage.
#' @param projected activate for additional datalayer with projected states.
#' @param projected.atoms indices of atoms over which will be summed
#' @param projected.energyintervall in which bands will be included (overwrites \code{projected.bands})
#' @param projected.bands used for normating of color. searches for maximum projected value in these bands
#' @param projected.orbitals list of orbitals to plot. To sum over orbitals 2 and 3 use \code{list(1,c(2,3),4)}
#' @param projected.colorpalette color palette for orbitals
#' @param ... further plotting parameters
#' @export
plot.bandsdata.contour<-function(bandsdata,band,sym=NA,n=201,breaks=12,colorpalette=colorRampPalette(c("red","blue","green"))
                                 ,zone=c("bz","basis"),linear=F
                                 ,projected=F,projected.atoms=1:bandsdata$natoms,projected.energyintervall=NULL,projected.bands=1:bandsdata$nbands
                           ,projected.orbitals=list(1,2,3,4)
                           ,projected.colors=colorRampPalette(c("red","blue","green"))(length(projected.orbitals))
                           ,xlab=expression(k[x]),ylab=expression(k[y]),...
                           ){ 
  kpoints <- (bandsdata$kpoints[,1:2])%*%(bandsdata$kbasis[1:2,1:2])
  data <- cbind(kpoints,1:nrow(kpoints))
  vec2d <- reciprocalbasis.getbrillouinzone(bandsdata$kbasis)  
  data<-brillouinzone.extendkpoints(vec2d,data)   
  if(length(sym)>1 || !is.na(sym)){
    data <- dataframe.applysymoperations(data,sym)
  }
  if(zone!="bz"){
    vec2d <- reciprocalbasis.getbrillouinzone(bandsdata$kbasis,zone=zone)
    center <- attr(vec2d,"scaled:center")
    if(!is.null(center)){
      data[,1:2] <- scale(data[,1:2],center=center,scale=F)
    }
    data<-brillouinzone.extendkpoints(vec2d,data) 
  }  
  data <- brillouinzone.selectkpoints(vec2d,data)
  data <- unique(round(data,6))
  data <- data[order(data[,1],data[,2],data[,3]),]
  data<- as.data.frame(data)
  colnames(data) <-c("kx","ky","index")
  xo <- seq(min(data[,1]),max(data[,1]),length=n)
  yo <- seq(min(data[,2]),max(data[,2]),length=n)
  e<-interp(data$kx,data$ky,bandsdata$bands[[band]]$simpledata[data$index,2],linear=linear,xo=xo,yo=yo)
  if(all(is.na(e$z))){
    warning("collinearity found, setting linear interpolation to TRUE")
    linear<-T  
    e <-interp(data$kx,data$ky,bandsdata$bands[[band]]$simpledata[data$index,2],linear=linear,xo=xo,yo=yo)
  }  
  plot(1,xlim=range(vec2d[,1]),ylim=range(vec2d[,2]),type="n",asp=1,xlab=xlab,ylab=ylab)
  levels <- pretty(range(e$z,finite=T),breaks)
  if(projected){
    colors <- lapply(projected.colors,FUN=function(x)colorRampPalette(c("white",x))(breaks))
    proj<- bandsdata.getprojecteddata(bandsdata,bands=projected.bands,atomindices=projected.atoms,energyintervall=projected.energyintervall)
    projbands <- lapply(proj$bands,FUN=function(band){
      band <-apply(band[,-(1:2)],MARGIN=2,FUN=as.numeric)
      do.call(cbind,lapply(projected.orbitals,FUN=function(orb){      
        if(length(orb)>1)
          b <- rowSums(band[,(orb)])
        else
          b <- as.numeric(band[,(orb)])
        return (b)
      }))
    })
    maxis <- max(do.call(rbind,projbands))
    pbreaks <- pretty(c(0,maxis),breaks)
    cols <- lapply(1:length(projected.orbitals),FUN=function(orb){
      da <- projbands[[paste0("band",band)]][,orb]
      p<-interp(data$kx,data$ky,da[data$index],linear=linear,xo=xo,yo=yo)        
      zcol <- cut(p$z,pbreaks)
      rgbcol <- col2rgb(colors[[orb]][zcol])
      return(rbind(rgbcol,as.numeric(zcol)/length(pbreaks)))
    })
    rgb <- lapply(1:3,FUN=function(rgb){
      rgbcol <- lapply(1:length(cols),FUN=function(x)cols[[x]][c(rgb,4),])        
      rgbcol <- do.call(rbind,rgbcol)
      rgbcol <- apply(rgbcol,2,FUN=function(column){
        value <-sum(column[c(1,3,5)]*column[c(2,4,6)], na.rm = T)/sum(column[c(2,4,6)], na.rm = T)
        if(is.na(value)) return(255)
        else return(value)
      })
      return(rgbcol)
    })

    col <- rgb(rgb[[1]],rgb[[2]],rgb[[3]],maxColorValue=255)
    ucol <- unique(col)
    m <- matrix(match(col,ucol),nrow=length(e$y),ncol=length(e$x))
    image(e$x,e$y,m,col=ucol,breaks=(0:(length(ucol))),useRaster=T,add=T)
    
  }else{
    image(e$x,e$y,e$z,add=T,breaks=levels,col=colorpalette(length(levels)-1),useRaster=T)
  }
  contour(e$x,e$y,e$z,add=T,levels=levels)
  polygon(vec2d)
}

#' Plots 3d bandsdata
#' 
#' \code{plot.bandsdata.3d} plots 3d bandsdata.
#' Projected states can be added as additional datalayer.
#' First make sure, you have a statisfying grid. See \code{\link{plot.bandsdata.grid}}.
#' 
#' @param bandsdata object of class bandsdata
#' @param bands indices of bands to plot
#' @param n resolution of datalayers
#' @param sym See \code{\link{dataframe.applysymoperations}} for usage and \code{\link{plot.bandsdata.grid}} for testing.
#' @param colorpalette colors which are used if \code{projected} is \code{False}
#' @param asp aspect ration between energy and k-space
#' @param breaks number of steps in color coding
#' @param zone see \code{\link{reciprocalbasis.getbrillouinzone}} for usage.
#' @param linear see \code{\link{akima::interp}} for usage.
#' @param projected activate for additional datalayer with projected states.
#' @param projected.atoms indices of atoms over which will be summed
#' @param projected.bands (optional) used for normating of color. searches for maximum projected value in these bands. Make sure your desired bands are included.
#' @param projected.energyintervall (optional) in which bands will be included for norming projected. Make sure your desired bands are in this intervall.
#' @param projected.orbitals list of orbitals to plot. To sum over orbitals 2 and 3 use \code{list(1,c(2,3),4)}
#' @param projected.colorpalette color palette for orbitals
#' @param ... further plotting parameters
#' @seealso \code{\link{rgl}}
#' @export
plot.bandsdata.3d<-function(bandsdata,bands,sym=NA,n=201,colorpalette=colorRampPalette(c("red","white","blue")),asp=0.1
                            ,breaks=256,zone=c("bz","basis")
                            ,linear=T
                            ,projected=F,projected.atoms=1:bandsdata$natoms,projected.orbitals=list(1,2,3,4)
                            ,projected.bands=1:bandsdata$nbands,projected.energyintervall=NULL
                            ,projected.colors= colorRampPalette(c("red","blue","green"))(length(projected.orbitals))
                                 ,xlab=expression(k[x]),ylab=expression(k[y]),zlab="energy (eV)",back="filled",front="filled",box=F,...
){  
  kpoints <- (bandsdata$kpoints[,1:2])%*%(bandsdata$kbasis[1:2,1:2])
  data <- cbind(kpoints,1:nrow(kpoints))
  vec2d <- reciprocalbasis.getbrillouinzone(bandsdata$kbasis)  
  data<-brillouinzone.extendkpoints(vec2d,data)   
  if(length(sym)>1 || !is.na(sym)){
    data <- dataframe.applysymoperations(data,sym)
  }
  if(zone!="bz"){
    vec2d <- reciprocalbasis.getbrillouinzone(bandsdata$kbasis,zone=zone)
    center <- attr(vec2d,"scaled:center")
    if(!is.null(center)){
      data[,1:2] <- scale(data[,1:2],center=center,scale=F)
    }
    data<-brillouinzone.extendkpoints(vec2d,data) 
  }  
  data <- brillouinzone.selectkpoints(vec2d,data)
  data <- unique(round(data,6))
  data <- data[order(data[,1],data[,2],data[,3]),]
  data<- as.data.frame(data)
  colnames(data) <-c("kx","ky","index")
  zrange <- range(sapply(bandsdata$bands[bands],FUN=function(band){
    return (band$simpledata[,2])
  }))
  pbreaks <- pretty(zrange,breaks)
  xo <- seq(min(data[,1]),max(data[,1]),length=n)
  yo <- seq(min(data[,2]),max(data[,2]),length=n)
  m <- cbind(vec2d[,1],vec2d[,2],rep(0,nrow(vec2d)))
  m <- rbind(m,m[1,])
  energiefac <- sum(abs(range(sapply(bands,FUN=function(band)bandsdata$bands[[band]]$simpledata[,2]))))
  xfac <- sum(abs(range(data[,1])))
  yfac <- sum(abs(range(data[,2])))
  plot3d(m,lwd=3,xlab=xlab,ylab=ylab,zlab=zlab,type="l",aspect=c(1,yfac/xfac,asp*energiefac/xfac),box=box)  
  if(projected){
    colors <- lapply(projected.colors,FUN=function(x)colorRampPalette(c("white",x))(breaks))
    proj<- bandsdata.getprojecteddata(bandsdata,bands=projected.bands,atomindices=projected.atoms,energyintervall=projected.energyintervall)
    projbands <- lapply(proj$bands,FUN=function(band){
      band <-apply(band[,-(1:2)],MARGIN=2,FUN=as.numeric)
      do.call(cbind,lapply(projected.orbitals,FUN=function(orb){      
        if(length(orb)>1)
          b <- rowSums(band[,(orb)])
        else
          b <- as.numeric(band[,(orb)])
        return (b)
      }))
    })
    maxis <- max(do.call(rbind,projbands))
    pbreaks <- pretty(c(0,maxis),breaks)
  } 
  for(iband in 1:length(bands)){
    e<-interp(data$kx,data$ky,bandsdata$bands[[bands[[iband]]]]$simpledata[data$index,2],linear=linear,xo=xo,yo=yo)
    if(all(is.na(e$z))){
      warning("collinearity found, setting linear interpolation to TRUE")
      linear<-T 
      e<-interp(data$kx,data$ky,bandsdata$bands[[bands[[iband]]]]$simpledata[data$index,2],linear=linear,xo=xo,yo=yo)
    }
    zcol <- cut(e$z,pbreaks)
    col <- colorpalette(length(pbreaks))[zcol]    
    if(projected){      
      cols <- lapply(1:length(projected.orbitals),FUN=function(orb){
        da <- projbands[[paste0("band",bands[[iband]])]][,orb]
        p<-interp(data$kx,data$ky,da[data$index],linear=linear,xo=xo,yo=yo)       
        zcol <- cut(p$z,pbreaks)
        rgbcol <- col2rgb(colors[[orb]][zcol])
        return(rbind(rgbcol,as.numeric(zcol)/length(pbreaks)))
      })
      rgb <- lapply(1:3,FUN=function(rgb){
        rgbcol <- lapply(1:length(cols),FUN=function(x)cols[[x]][c(rgb,4),])        
        rgbcol <- do.call(rbind,rgbcol)
        rgbcol <- apply(rgbcol,2,FUN=function(column){
          value <-sum(column[c(1,3,5)]*column[c(2,4,6)], na.rm = T)/sum(column[c(2,4,6)], na.rm = T)
          if(is.na(value)) return(255)
          else return(value)
          })
        return(rgbcol)
      })
      col <- rgb(rgb[[1]],rgb[[2]],rgb[[3]],maxColorValue=255)
    }  
    persp3d(e$x,e$y,e$z,lit=T,col=col,back=back,add=T,front=front,...)
  }  
}

#' Applies chain of symmetric operation on a dataframe
#' 
#' \code{dataframe.applysymoperations} applied a chain of 2d symmetric operations on the first two columns of the dataframe
#' 
#' @param dataframe dataframe with at least two columns
#' @param symoperations list of symmetric operations. \code{list(c("reflection",-30),c("rotation",120,60))} 
#' for example will first reflect across a line with the slope of -30°, afterwards the new dataset (old+reflected) will be rotated by 120° and 60°
#' the resulting dataset will consist of 6 combined datasets. The old, the reflected and the old+reflected, rotated by 120° and by 60°.
#' @param center (optional) 2d point, which will be used as center for all symoperations
#' @export
dataframe.applysymoperations<-function(dataframe,symoperations,center=NULL){
  stopifnot(ncol(dataframe)>=2)
  names <- colnames(dataframe)  
  for(sym in symoperations){
    if(length(sym)>1){
      for(i in 2:length(sym)){
        rmatrix <- switch(sym[[1]],
                          "reflection" = matrix.reflection2d.degree(as.numeric(sym[[i]])),
                          "rotation" = matrix.rotation2d.degree(as.numeric(sym[[i]])),
                          default = NA)
        if(!is.null(center)){
          d2 <- scale(scale(dataframe[,1:2],center=center,scale=F)%*%t(rmatrix),center=-center,scale=F)
        }
        else{
          d2 <-data.matrix(dataframe[,1:2])%*%t(rmatrix)
        }
        if(ncol(dataframe)>2){
          d2 <- cbind(d2,dataframe[,3:ncol(dataframe)])
        }
        if(i==2)
          d <- d2
        else
          d <- rbind(d,d2)        
      }
      if(ncol(d)==ncol(dataframe)){

        colnames(d)<-names
        dataframe <- rbind(dataframe,d)
      }
    }
  }
  dataframe <- dataframe[order(dataframe[,1],dataframe[,2]),]
  dataframe<-unique(dataframe)  
  return(dataframe)
}

