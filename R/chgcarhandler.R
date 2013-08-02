#source("./R/poscarhandler.R")

read.chgcar<-function(filename){
  chgcar <- list()
  class(chgcar)<-"chgcar"
  chgcar[["name"]]<-foldername <- tail(strsplit(filename,"/")[[1]],1)
  data <- read.table(filename,strip.white=T,sep="\n",stringsAsFactors=F)
  chgcar$poscar <- parse.poscar(data)
  pl <- chgcar$poscar$length
  coords <- as.numeric(strsplit(data[(pl+1),1]," +")[[1]])
  coordcount <- prod(coords)
  block <-do.call(c,strsplit(data[(pl+2):nrow(data),1]," +"))
  rm(data)
  block <- as.numeric(block[1:coordcount])
  x <- seq(0,1,len=coords[1]+1)[1:coords[1]]
  y <- seq(0,1,len=coords[2]+1)[1:coords[2]]
  z <- seq(0,1,len=coords[3]+1)[1:coords[3]]
  pos <- cbind(rep(x,coords[2]*coords[3]),
  rep(do.call(c,lapply(y,rep,coords[1])),coords[3]),
  do.call(c,lapply(z,rep,coords[1]*coords[2])))
  poscart <- pos %*%(chgcar$poscar$basis*chgcar$poscar$a)
  chgcar$cellvolume <- det(chgcar$poscar$basis)*chgcar$poscar$a
  chgcar$data <- data.frame(poscart, block / chgcar$cellvolume)
  chgcar$xy <- unique(chgcar$data[,c(1,2)])
  chgcar$z <- unique(chgcar$data[,3])
  chgcar$n <- coords  
  chgcar$pointvolume <- det(diag(c(x[2],y[2],z[2]))%*%(chgcar$poscar$basis*chgcar$poscar$a))
  return(chgcar)
}

print.chgcar <- function(chgcar,...){
    cat("chgcar\n")
    for (name in names(chgcar))
    {
      if (length(chgcar[[name]])>1)
      {
        cat(paste0(" ",name," length:",length(chgcar[[name]])," (",class(chgcar[[name]]),")","\n"))
      }
      else
      {
        cat(paste0(" ",name," ",chgcar[[name]],"\n"))
      }
    }
  
}

chgcar.sumoverlayer <- function(chgcar,layer,layers){
  layer <- sort(layer)
  diff <- layer[-1]-layer[-length(layer)]
  cuts <- which(diff>1)
  stopifnot(length(cuts)==0)
  pos <- chgcar$poscar
  xyz <- chgcar$data[,1:3]
  vol <- chgcar$cellvolume/prod(chgcar$n)
  layerindices <- poscar.getatomlayerindices(pos,layers)
  atoms <- poscar.getbasisconvertedatoms(pos)
  minlayer <- min(layer)
  if(minlayer==1){
    minz <- 0
  }
  else{
    minz <- mean(c(mean(atoms[layerindices==(minlayer-1),3]),mean(atoms[layerindices==(minlayer),3])))
  }
  maxlayer<-max(layer)
  if(maxlayer==layers){
    maxz <- max(xyz[,3])
  }
  else{
    maxz <- mean(c(mean(atoms[layerindices==(maxlayer+1),3]),mean(atoms[layerindices==(maxlayer),3])))
  }
  selector <- xyz[,3]>=minz & xyz[,3]<=maxz
  return(sum(chgcar$data[selector,4]*vol))
}

chgcar.sum<- function(chgcar){
  return(chgcar.sumoverlayer(chgcar,1,1))
}