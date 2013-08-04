require(XML)
require(lattice)
require(grid)
require(snowfall)

read.bandsdata <- function(xmlfile){  
  result <- list()
  class(result)<-"bandsdata"
  print(paste("processing:",xmlfile))
  ### parsing xml tree 
  xmlobj <- xmlTreeParse(xmlfile, getDTD=FALSE, useInternalNodes=TRUE)
  # find k basis in xml
  result$kbasis <- xpathSApply(xmlobj, "/modeling/structure[@name='initialpos']/crystal/varray[@name='rec_basis']/*",
                               function(x) as.numeric(strsplit(xmlValue(x)," +")[[1]][-1]))
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

#custom print for bandsdata
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

#custom print for banddata
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


bandsdata.calcsympointpath<-function(bandsdata,sympointpath){
  # adds a path from sympoints to bandsdata
  # ex. sympointpath=list(c(1,2),c(3,4)) 
  # goes from 1 to 2 and 3 to 4
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

bandsdata.addsympoint <- function(bandsdata,index){
  bandsdata$sympoints <-sort(c(bandsdata$sympoints[index],bandsdata$kpointsflat))
  bandsdata$sympointindices <- sort(c(bandsdata$sympointindices,index))
  bandsdata$sympointindicesrev <- sort(c(bandsdata$sympointindicesrev,index))
  print("added sympoint to bandsdata")
  print(paste("new sympointindices:",paste(bandsdata$sympointindices,collapse=" ")))
  return(bandsdata)
}

plot.bandsdata <- function(bandsdata
                           ,bands=1:length(bandsdata$bands)
                           ,sympointpath=NULL
                           ,col="black"
                           ,type="l"
                           ,fermi=F
                           ,symnames=NULL
                           ,symcolor="red"
                           ,symlty=3
                           ,xlim=NULL
                           ,xaxs="i"
                           ,yaxs="i"
                           ,energyoffset=0
                           ,...){
  print("plotting bandsdata")
  rng <- bandsdata$range
  if(!is.null(sympointpath))
  {
    bandsdata <- bandsdata.calcsympointpath(bandsdata,sympointpath)    
  }
  bandsdata$energyoffset <- energyoffset  
  if(!is.null(bandsdata$sympath))
  {
    rng[1:2,1]<-range(bandsdata$sympath$data[,2])
  }
  if (is.null(xlim))
    xlim <- rng[1:2,1]
  plot(rng[,1:2],type="n",xaxt="n",xlim=xlim,xaxs=xaxs,yaxs=yaxs ,...)
  if (!type=="n") plot.bandsdata.addbands(bandsdata,bands,col=col,...)
  if (fermi)
  {
    plot.bandsdata.addfermi(bandsdata,xlim=xlim,...)
  }
  if(!is.null(symnames))
  {
    plot.bandsdata.addsymnnames(bandsdata,symnames,symcolor=symcolor,symlty=symlty)
  }
  return(bandsdata)
}


plot.bandsdata.addsymnnames <- function(bandsdata,symnames,symcolor="red",symlty=3,...){
  print("adding symnames")
  sp <- bandsdata$sympoints
  if(!is.null(bandsdata$sympath))
  {
    sp<-bandsdata$sympath$sympoints
  }
  select <- !duplicated(sp)
  abline(v=sp[select],col=symcolor,lty=symlty,...)
  axis(1,at=sp[select],labels=symnames[select],...)
}

plot.bandsdata.addfermi<-function(bandsdata,fermicolor="blue",lty=3,xlim=NULL,...){
  e <- 0
  if(is.null(bandsdata$energyoffset))
    e <- bandsdata$energyoffset
  rng <- bandsdata$range
  if (is.null(xlim))
    xlim <- rng[1:2,1]
  abline(h=e,col=fermicolor,lty=lty,...)
  #text(xlim[1],0,"fermi",col=fermicolor,adj=c(-0.2,-0.5),...)
}

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

plot.bulkbands.add<-function(bulkbands,col="grey",...){
  print("plotting bulk polygons")
  for(p in bulkbands$poly)
  {
    polygon(p$x,p$y,col=col,border=NA)
  }
}

bandsdata.getprojecteddata <- function(bandsdata,atomindices=1:bandsdata$natoms,bands=1:bandsdata$nbands,cpus=1){
  print("calculating projecteddata")
  print("sum over atoms")
  nkpoints <- length(bandsdata$kpointsflat)
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

plot.projectedbands.add <- function(projectedbands
                                    ,bands=1:length(projectedbands$bands)
                                    ,orbitals=list(1,2,3,4)
                                    ,col.palette=colorRampPalette(c("red","blue","green"))
                                    ,pch=15:(14+length(orbitals))
                                    ,cex=0.8
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
  k <- projectedbands$bands$band1[,1]
  i <- 1:nrow(projectedbands$bands$band1)
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
  col <- col.palette(length(orbitals))
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
  legend("topright",legend=parse(text=names),bg="white",col=col[1:length(orbitals)],pch=pch,cex=legendcex,seg.len=1,x.intersp=0.5,y.intersp=0.8,...)

}

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

bandsdata.getbanddistance<-function(bandsdata,kpoint,bands){
  e1 <- bandsdata$bands[[bands[[1]]]]$simpledata[kpoint,2]
  e2 <- bandsdata$bands[[bands[[2]]]]$simpledata[kpoint,2]
  return(abs(e1-e2))
}

bandsdata.getenergy<-function(bandsdata,kpoint,band){
  return(bandsdata$bands[[band]]$simpledata[kpoint,2])
}

bandsdata.fit <- function(bandsdata,bandnr,kpoints,fitname="dirac",startingparameters,additionalparameters=NULL)
{
  if (!is.null(additionalparameters)&&!is.null(additionalparameters$k0)&&additionalparameters$k0>1)
  {
    additionalparameters$k0<- bandsdata$kpointdistances[additionalparameters$k0]
  }
  objectivem <- function(parameters,fitfunction,data,additionalparameters=list())
  {
    err <- sum((data$energy-fitfunction(parameters,data$kpoint,additionalparameters))^2)
    return (err)
  }
  fitfunction<-get(paste("fitfunction",fitname,sep="."))
  data <- bandsdata$bands[[bandnr]]$simpledata[kpoints,]
  data$kpoint <- bandsdata$kpointdistances[kpoints]
  if(is.null(additionalparameters))
    obj <- optim(startingparameters,objectivem,fitfunction=fitfunction,data=data,method="BFGS")
  else
    obj <- optim(startingparameters,objectivem,fitfunction=fitfunction,data=data,method="BFGS",additionalparameters=additionalparameters)
  class(obj)<-"bandsfit"
  obj$ap <-additionalparameters
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
bandsdata.fit.dirac.makeparameters<-function(vF)
{
  p <- c(vF)
  names(p) <- c("vF")
  return (p)
}

bandsdata.fit.dirac.makeconstants<-function(bandsdata,bands,k0=0,factor=1)
{
  d <- bandsdata.getbanddistance(bandsdata,k0,bands)
  evalence <- bandsdata.getenergy(bandsdata,k0,bands[[1]])
  ap <- list()
  ap$EF<-d/2+evalence
  ap$k0<-k0
  ap$Eg<-d
  ap$factor<-factor
  return (ap)
}
bandsdata.fit.dirac.function<-function(parameters,kpoints,additionalparameters)
{
  p<-parameters ## p[[1]] vF
  ap <- additionalparameters  
  return(sign(ap$factor)*sqrt((ap$Eg/2)^2+(abs(kpoints-ap$k0)*p[[1]])^2)+ap$EF)
}
bandsdata.fit.quadratic.makeparameters<-function(m)
{
  p <- c(m)
  names(p) <- c("m*")
  return (p)
}
bandsdata.fit.quadratic.makeconstants<-function(bandsdata,band,k0=0)
{
  ap <- list()
  ap$E<-bandsdata.getenergy(bandsdata,k0,band)
  ap$k0<-k0
  return (ap)
}
bandsdata.fit.quadratic.function<- function(parameters,kpoints,additionalparameters)
{
  p<-parameters ## p[[1]] m*
  ap <- additionalparameters  
  return(p[[1]]*(kpoints-ap$k0)^2+ap$E)
}
predict.bandsfit<-function(o,x)
{
  data <- data.frame(x,o$fitfunction(o$par,x,o$ap))
  names(data) <- c("kpoint","energy")
  return(data)
}
bandsfit.dirac.getv<-function(bandsfit)
{
  p <- bandsfit$par[[1]]
  hbar <- 6.582119*10^-16 ## eVs
  kfac <- 10^-10
  return(p*kfac/hbar)
}
bandsfit.quadratic.getm<-function(bandsfit)
{
  p <- bandsfit$par[[1]]
  hbar <- 6.582119*10^-16 ## eVs
  kfac <- 10^-10
  me <- 0.510998928*10^3
  c <- 3*10^8
  return(hbar^2*c^4/p/2/me)
} 
print.bandsfit<-function(bandsfit)
{
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
# add.bandsfit<-function(bandsfit,kpoints=bandsfit$knumbers,n=201,energyoffset=NULL,...){
#   if(is.null(energyoffset))
#   {
#     if(!is.null(bandsfit$energyoffset))
#     {
#       energyoffset <- bandsfit$energyoffset
#     }
#     else
#     {
#       energyoffset<-0
#     }
#   }
#   print("plotting bandsfit")
#   x <- bandsfit$kdist[range(kpoints)]
#   kdist <- seq(x[[1]],x[[2]],length.out=n)
#   data <- predict(bandsfit,kdist)
#   if(!is.null(bandsfit$sympath))
#   {
#     x <- range(bandsfit$sympath$data[bandsfit$sympath$data[,1]%in%kpoints,2])
#   }
#   else
#   {
#     x <- bandsfit$korig[range(kpoints)]    
#   }  
#   k <- seq(x[[1]],x[[2]],length.out=n)
#   lines(k,data$energy+energyoffset,col="red",...)
# }
# bandsdata <- getbandsdata("vasprun.xml")
# bandsdata <- plot(bandsdata,fermi=T,symnames=c(1,"L",expression(Gamma),"X"),ylim=c(-10,5),sympointpath=list(c(2,3),c(3,4)),energyoffset=-getenergy.bandsdata(bandsdata,kpoint=113,band=8))
# sp <- fitmakeparameters.quadratic(m=5)
# ap <- fitmakeaddparameters.quadratic.bandsdata(bandsdata,k0=113,band=8)
# addnumbers.bandsdata(bandsdata)
# o <- fit.bandsdata(bandsdata,8,105:113,"quadratic",sp,ap)
# add.bandsfit(o)
# ap <- fitmakeaddparameters.quadratic.bandsdata(bandsdata,k0=113,band=9)
# o <- fit.bandsdata(bandsdata,9,113:118,"quadratic",sp,ap)
# add.bandsfit(o)
####debug
##fit
#  bandsdata <- getbandsdata("./d_silicene_simple/d05_bands/a_2.240/vasprun.xml1")
#  plot(bandsdata)
#  addnumbers.bandsdata(bandsdata)
#  sp <- fitmakeparameters.dirac(1000)
#  ap <- fitmakeaddparameters.dirac(k0=80,factor=-1)
#  o <- fit.bandsdata(bandsdata,4,kpoints=75:80,"dirac",sp,ap)
#  add.bandsfit(o)
# print(o)
#   o <- fit.bandsdata(bandsdata,4,kpoints=80:86,"dirac",sp,ap)
#   add.bandsfit(o)
#   print(o)
#   sp <- fitmakeparameters.dirac(1000)
#   ap <- fitmakeaddparameters.dirac(k0=80)
#   o <- fit.bandsdata(bandsdata,5,kpoints=75:80,"dirac",sp,ap)
# add.bandsfit(o)
# print(o)
# o <- fit.bandsdata(bandsdata,5,kpoints=80:86,"dirac",sp,ap)
#  add.bandsfit(o)
#  print(o)
# ##other stuff
#banddata <- getbandsdata("./d_argentum_slab/d20_7_bands/kline_0.0000+kline_0.0000/vasprun.xml1")
#banddata <- getbandsdata("./d_argentum_bulk/d11_bands_111/a_4.030/vasprun.xml1")
#load.calculationdata("banddata1")
#plot(bandsdata$"./d_argentum_bulk/d11_bands_111/"$a_4.030$banddata,sympoints=list(c(3,4),c(1,2)),symnames= c("Γ","M","K","Γ"))
#load.calculationdata("banddata2")
#plot(bandsdata$"./d_dual_ar_si/d16_8_on_7_bands/"$a_4.030$banddata,col=gray(0.4))
#asel <- which(bandsdata$"./d_dual_ar_si/d16_8_on_7_bands/"$a_4.030$poscar$atoms[,7]=="Si")
#addprojected.bandsdata(bandsdata$"./d_dual_ar_si/d16_8_on_7_bands/"$a_4.030$banddata,asel,cpus=4)