Rvasp
=====
Tools for loading, manipulating and plotting VASP files within R

### Workflow
To be somewhat efficent, even in research, I developed a workflow for scientific computation, which works quite well for me.
![scientific computation - workflow](../../raw/master/examples/workflow-01.png "scientific computation - workflow")
VASPmanager can be found [here](https://github.com/paulpflug/VASPmanager) 

### Important
Many features originated by a need to solve a specific problem and are propably not generalized enough to declare the general problem fully solved.
All functions are implemented in pure R and are mostly easy to understand. Feel free to read the [source](../../tree/master/Rvasp/R), improve and submit changes.

## Index

- [Install](#install)
- [What it does](#what-it-does)
- [Simple examples](#simple-examples)
  - [Raw manipulation of a POSCAR](#raw-manipulation-of-a-poscar)
  - [Plotting a single 2d surface from top an side](#plotting-a-single-2d-surface-from-top-an-side)
  - [Plotting of a STM image](#plotting-of-a-stm-image)
  - [Plotting of a bandstructure with corresponding DOS](#plotting-of-a-bandstructure-with-corresponding-dos)
  - [Plotting of an elevated bandstructure](#plotting-of-an-elevated-bandstructure)
  - [Plotting of a E over a curve](#plotting-of-a-e-over-a-curve)
  - [Plotting of bulk bands](#plotting-of-bulk-bands)
  - [Plotting a brillouin zone](#plotting-a-brillouin-zone)
  - [Plotting a periodic table](#plotting-a-periodic-table)
- [Complex examples](#complex-examples)
  - [Plotting of 2d surfaces from top an side](#plotting-of-2d-surfaces-from-top-an-side)
  - [Plotting of bandstructures, incl. fitting of direction-depended dirac-cones](#plotting-of-bandstructures-incl-fitting-of-direction-depended-dirac-cones)
  - [Plotting of STM images with top layers atom positions](#plotting-of-stm-images-with-top-layers-atom-positions)
  - [Call for other scenarios](#call-for-other-scenarios)
- [Overview over all functions](#overview-over-all-functions)
  - [POSCAR based](#poscar-based)
  - [CHGCAR based](#chgcar-based)
  - [Vasprun.xml based](#vasprunxml-based)
  - [Calculation based](#calculation-based)
  - [Miscellaneous](#miscellaneous-1)
  
## Install
##### Requirements
* (Linux) `libxml2-dev` as a requirement of the R XML package and `xorg-dev`, `libglu1-mesa-dev` for the rgl package
* the R packages snowfall, akima, XML and rgl

```S
install.packages(pkgs=c("snowfall","akima","XML","rgl"),repos="http://cran.rstudio.com/")
```

##### Install as a package
* download [Rvasp package](https://www.dropbox.com/s/gtuz0o79zurro6a/Rvasp_0.2.tar.gz)
* install in R

```S
install.packages(pkgs="Rvasp_0.2.tar.gz",repos=NULL, type="source")
```

##### Use the sourcecode
* clone or download & extract zip from the right
* open `Rvasp.Proj` in [RStudio](http://www.rstudio.com/)
* build

##### Load in R

```S
library(Rvasp)
```

## What it does
Contains a huge set of functions to work with the VASP output, namely the POSCAR, CHGCAR and vasprun.xml. This is extended by adding a wrapper to manage several outputs at the same time. Beside of the VASP thing, there are few miscellaneous functions

##### POSCAR
* read
* write
* manipulate
* plot

##### CHGCAR
* read
* calculate and plot stm

##### Vasprun.xml 
* read and plot (projected) dosdata
* read and plot (projected) bandsdata and fit a dirac cone or a quadratic function to a band

##### Several outputs at the same time
The output has to be organized in the following scheme:   

```
/manyCalculations/Calculation1/Parameter1/VASPfile1   
/manyCalculations/Calculation1/Parameter1/VASPfile2   
..   
/manyCalculations/Calculation1/Parameter2/VASPfile1   
/manyCalculations/Calculation1/Parameter2/VASPfile2   
..   
/manyCalculations/Calculation2/Parameter1/VASPfile1   
/manyCalculations/Calculation2/Parameter1/VASPfile2   
etc.
```

with these three levels:
* Calculations
* Calculation
* Parameter

On calculation level there are the following functions:
* get / plot / fit an e over a curve
* get / plot a local dos curve
* get / plot bulk band data

##### miscellaneous
* plotting a periodic table

## Simple examples

##### Raw manipulation of a POSCAR
First read a poscar `poscar <- read.poscar(file="POSCAR")`.
The most important objects can be accessed like this:

```S
poscar$a
  [1] 4.03

poscar$basis
       [,1] [,2] [,3]
  [1,]  0.0  0.5  0.5
  [2,]  0.5  0.0  0.5
  [3,]  0.5  0.5  0.0

poscar$atoms
    rx ry rz movex movey movez type
  1  0  0  0     T     T     T   Ag
```

R built-in tools are powerfull:

```S
poscar$atoms$rz[poscar$atoms$type=="Ag"]
```

or

```S
poscar$atoms[poscar$atoms$type=="Ag",3]
```

for example will give you all z components of all silver atoms.

For manipulation, changes have to be saved,

```S
poscar$atoms[poscar$atoms$type=="Ag",3]<-poscar$atoms[poscar$atoms$type=="Ag",3]+0.2
```

will increase all third coordinates of all silver atoms by 0.2.
After desired manipulation write the poscar to a file:

```S
write.poscar(poscar=poscar,file="POSCAR")
```

##### Plotting a single 2d surface from top an side

```S
data(silverslab)
temppos <- poscar.getshapedposcar(silverslab,x=10,y=10)
# left plot - top view
plot(temppos,basis=F,unitcell=T)
plot.poscar.addlayers(poscar=temppos,layer=1:2,layers=12,color=c("blue","red"),size=4)
plot.poscar.addnumbers(poscar=temppos,layer=1:2,layers=12)
plot.atoms.adddistance(atoms=temppos$atoms[c(17,161),],basis=temppos$basis*temppos$a,selectedsize=4)
# right plot - side view
plot(temppos,basis=F,unitcell=F,direction=1)
plot.poscar.addlayers(poscar=temppos,layer=1:2,layers=12,color=c("blue","red"),size=4,direction=1)
plot.poscar.addlayerdistance(poscar=temppos,layer=1:2,layers=12,direction=1)
plot.poscar.addlayers(poscar=temppos,layer=3:12,layers=12,color="grey",size=4,direction=1)
```
![POSCAR picture - top](../../raw/master/examples/poscar_top.png "Structure of a silverslab plotted in R - topview")
![POSCAR picture - side](../../raw/master/examples/poscar_side.png "Structure of a silverslab plotted in R - sideview")

##### Plotting of a STM image

```S
data(silverstm)
plot.stm(silverstm,super=5,xlab="x (nm)",ylab="y (nm)")
plot.stm.addunitcell(silverstm,atomnumber=1)
plot.stm.addatoms(silverstm,super=5,xlim=c(0,0.75),ylim=c(0,0.75),atomselector=1,atomsize=10,col="red")
plot.stm.addatoms(silverstm,super=5,xlim=c(0,0.75),ylim=c(0,0.75),atomselector=2,atomsize=10,col="blue")
```
![STM picture](../../raw/master/examples/stm.png "STM of a silver surface created in R")

##### Plotting of a bandstructure with corresponding DOS

```S
data(silverbands)
data(silverdos)
newfermi <- silverbands$efermi-silverdos$efermi
bands<-plot.bandsdata(silverbands,sym.labels=c(expression(Gamma),"L","W","X",expression(Gamma),"K","X"),fermi=T,ylim=c(-10,20),energyoffset=newfermi)
proj <- bandsdata.getprojecteddata(bands)
plot.projectedbands.add(proj,orbitals=c(1,2,3,4),cex=1,legendcex=1.2)
dosdata <- plot.dosdata(silverdos,flip=T,fermi=T,col.fermi="Blue",xlim=c(0,0.5),ylim=c(-10,20))
plot.dosdata.add(dosdata,orbitals=c(1:4,"all"),type="polygon",col=c(colorRampPalette(c("red","blue","green"))(4),"grey"))
```
![BANDS picture](../../raw/master/examples/bands.png "Bands of a silver created in R")
![DOS picture](../../raw/master/examples/dos.png "DoS of a silver created in R")
##### Plotting of an elevated bandstructure

```S
data(silverbands)
data(silverdos)
newfermi <- silverbands$efermi-silverdos$efermi
bands<-plot.bandsdata(silverbands,sympointpath=list(c(1,2),c(2,3),c(3,4),c(7,6)),sym.labels=c(expression(Gamma),"L","W","X",expression(Gamma),"K","X"),fermi=T,ylim=c(-10,20),energyoffset=newfermi)
proj <- bandsdata.getprojecteddata(bands)
plot.projectedbands.add(proj,orbitals=c(1,2,3,4),cex=1,legendcex=1.2)
axis(side=3)
zoomplot(xpos=c(38,55),ypos=c(65,82),xmean=0.35,xwidth=0.025,ymean=6.1,ywidth=1,xaxt="n",lines=F)
plot.bandsdata.addbands(bands)
plot.bandsdata.addsymnnames(bands,labels=c(expression(Gamma),"L","W","X",expression(Gamma),"K","X"))
plot.projectedbands.add(proj,orbitals=c(1,2,3,4),cex=1,legend=NULL)
endzoomplot()
```
![Elevated BANDS picture](../../raw/master/examples/ebands.png "Bands of a silver created in R")
##### Plotting of a E over a curve

```S
data(silverea)
plot.calculation.ea(silverea$LDA,energyshift=silverea$LDAsingleAtom$a_20.000$energy,ylim=c(-3.8,-2.1))
plot.calculation.ea.addpoints(silverea$GGA,energyshift=silverea$GGAsingleAtom$a_20.000$energy,pch=3)
```
![E over a plot](../../raw/master/examples/eaplot.png "E over a plot created in R")

##### Plotting of bulk bands
Download silverbulkbands.RData from examples

```S
load("silverbulkbands.RData")
bandsdata <- silverbulkbands$slab$a_4.030$banddata
bulkbands <- calculation.getbulkbands(silverbulkbands$bulk)
plot.bandsdata(bandsdata,type="n",sym.lty=NA,fermi=F,ylim=c(-10,10))
plot.bulkbands.add(bulkbands)
plot.bandsdata.addbands(bandsdata=bandsdata,col="black")
plot.bandsdata.addfermi(bandsdata)
plot.bandsdata.addsymnnames(bandsdata,labels=c(expression(Gamma),"M","K",expression(Gamma)))
```
![Bulk bands](../../raw/master/examples/bulkbands.png "Silver slab bands underlayed by silver bulk bands, created in R")

##### Plotting a brillouin zone

```S
data(silverslab)
bzs<-poscar.getbrillouinzones(silverslab)
plot(bzs)
superslab <- poscar.createsupercell(silverslab,A=diag(c(2,2,1))) # 2x2x1 supercell
bzs2<-poscar.getbrillouinzones(superslab,extend=2)
plot.brillouinzones.add(bzs2,col="red")
plot.brillouinzones.addsympoints(brillouinzone=bzs,directcoordinates=list(c(1,0),c(0,0),c(1/2,1/2)),labels=c("K",expression(Gamma),"M"))
plot.brillouinzones.addsympoints(brillouinzone=bzs,directcoordinates=list(c(1,0),c(0,0),c(1/2,1/2)),labels=c("K",expression(Gamma),expression(Gamma)),textpos=4,col="red")
```
![Brillouinzone](../../raw/master/examples/brillouinzone.png "Brillouinzone of 1x1 and sqrt(3)xsqrt(3) silver, generated by R")

##### Plotting a periodic table

```S
plot.periodictable(highlights=1:2,highlighttexts=list(1.008,4.0026),underlay=3:10)
```
![Periodic table](../../raw/master/examples/periodic_table.png "Periodic table generated by R")


## Complex examples

##### Plotting of 2d surfaces from top and side

```S
require(Rvasp)
folders <- c(
  "./Folder1/",
  "./Folder2/",
  "./etc/"
  )
layers <- c(10,10,11,...) # count of layers in the respective structures
layerplot <- 3 # number of layers you want to plot
if(!exists(calc))
  calc <- read.calculations(folders=folders)
for(i in (1:length(folders))){
  pos <- calc[[folders[[i]]]][[1]]$contcar
  pos <- poscar.rotate2d(pos)
  l <- layers[[i]]
  pos <- poscar.alignbylayer(pos,l,l,alignto=c(0.01,0.01))
  superpos <- poscar.getshapedposcar(pos,20,20) 
  plot.poscar(superpos,basis=F,unitcell=T,ylim=c(0,20),xlim=c(0,29))
  plot.poscar.addlayers(superpos,layers=l,layer=(l-layerplot+1):l,size=2)
  superpos <- poscar.extractlayers(superpos,layer=(l-layerplot+1):l,layers=l,vacuum=c(0,0,2))
  plot(superpos,basis=F,unitcell=F,direction=1)
  plot.poscar.addlayers(superpos,layers=layerplot,layer=1:layerplot,size=2 ,direction=1)
  plot.poscar.addlayerdistance(superpos,layers=layerplot,layer=1:layerplot,color=c("black"),cex=0.7)
}
```

##### Plotting of bandstructures, incl. fitting of direction-depended dirac-cones

```S
require(Rvasp)
folders <- c(
  "./Folder1/",
  "./Folder2/",
  "./etc/"
  )
hexsym <- c(expression(Gamma),"M","K",expression(Gamma))
efermi <- c(0.5241,0.62341,...)
### fitting parameters
sp <- bandsdata.fit.dirac.makeparameters(vF=20)
bands <- list(c(28,29),c(28,29),...) # Bands which to fit
k0 <- list(c(80,81),c(80,81),...)# Kpoints where to fit
kwide<- 5 # how many kpoints to include
factor <- c(-1,1)
#######
filename<-"banddata.RData"
if((!load.calculations(filename))){
  data <- read.calculations(folders=folders,bandsxmlname="vasprun.xml")
  save(data,file=filename)
}
for (i in (1:length(folders))){   
  para <- data[[folders[[i]]]][[1]]
  if(!is.null(para)){
    eoff<-0
    if(!is.na(efermi[i])){
      eoff <- efermi[i]-para$banddata$efermi
    }
    projectedbands <- bandsdata.getprojecteddata(para$banddata,energyintervall=ylim+eoff)     
    plot(para$banddata,ylim=ylim,sym.labels=hexsym,fermi=T,energyoffset=eoff)
    plot.projectedbands.add(projectedbands,orbitals=list(1,c(4,2),3),usetransparent=F,cex.cutoff=0.35,energyoffset=eoff)
    plot.bandsdata.addnumbers(para$banddata)  
    ## DIRACFIT
    for (j in 1:2){
      for(k in 1:2){
        know <- k0[[i]][[k]]
        ap <- bandsdata.fit.dirac.makeconstants(para$banddata,bands=bands[[i]],k0=know,factor=factor[[j]])
        kpoints <- sort(know:((factor[[k]]*kwide)+know))
        o <- bandsdata.fit(para$banddata,bands[[i]][[j]],kpoints=kpoints,"dirac",sp,ap)
        plot.bandsfit.add(o)
      }
    }
  }
}
```


##### Plotting of STM images with top layers atom positions

```S
require(Rvasp)
folders <- c(
  "./Folder1/",
  "./Folder2/",
  "./etc/"
  )
filename <- "stmdata.RData" # savefile for calculated STMs
cutoff <- 0.01 # cutoff for stm calculation
layers <- c(10,10,11,...) # count of layers in the respective structures
for (i in 1:length(folders)){  
  if(!load.calculations(filename) ){
    stmdata <- list()
    rawstmdata <- read.calculations(folders[[i]],chgcarname="PARCHG",update=update)
    for (name in folders){ 
      for (para in names(rawstmdata[[name]])){
        if (!is.null(rawstmdata[[name]][[para]][["chgcar"]])){
          stmdata[[name]]<-list()
          stmdata[[name]][[para]]<- stm(rawstmdata[[name]][[para]][["chgcar"]],cutoff,cpus=1)
        }
      }      
    }
    save(stmdata,file=filename)
  }
  for (name in folders)
  {      
    for (para in names(stmdata[[name]]))
    {
      para2 <- strsplit(para,"_")[[1]][[2]]
      parsedpara <- round(as.numeric(strsplit(para2," +")[[1]]))
      parsedpara <- paste(parsedpara,collapse="_")
      stm <- stmdata[[name]][[para]] 
      print(paste("creating stm from",folders[[i]],"with",para))
      xlim <- ylim <- c(0,4)
      super <- 7
      plot.stm(stm,super=super,xlim=xlim,ylim=ylim)
      layerindices <- poscar.getatomlayerindices(stm$poscar,layers[[i]])
      plot.stm.addatoms(stm,layerindices==1,super=super,xlim=xlim/2,ylim=ylim/2)
      plot.stm.addatoms(stm,layerindices==2,super=super,xlim=xlim/2,ylim=ylim/2)
      plot.stm.addunitcell(stm,col="black",lwd=2)
    }    
  }  
}
```


##### Call for other scenarios

if you encounter scenarios you think sharing will help others, feel free to submit these.

## Overview over all functions
To get information about a function use `help(functionname)` or `?functionname` in R.
If you want to look under the hood of a function just use `functionname` (without brackets) in R or look directly in the sourcecode here on Github.
All functions are given with a typical call
### POSCAR based

```S
## Basics
poscar <- read.poscar(file="POSCAR")
print(poscar) # real function name is print.poscar
printraw.poscar(poscar)
write.poscar(poscar=poscar,file="POSCAR")


## Manipulation
atoms <- poscar.getbasisconvertedatoms(poscar=poscar)
poscar$atoms <- atoms.convertbasis(atoms=atoms,basis=poscar$basis*poscar$a)
poscar <- poscar.rotate2d(poscar) # poscar$basis[1,2] will equal 0 afterwards
recbasis <- poscar.getreciprocalbasis(poscar)
vacuum <- poscar.getvacuum(poscar)
poscar <- poscar.setvacuum(poscar,vacuum=c(0,0,10),center=T)
poscar <- poscar.centeratoms(poscar,direction=3,position=0.7)
poscar <- poscar.sortatoms(poscar,sortindices=c(7,3)) # first by atom sort, then by z-direction
newposcar <- poscar.extractatoms(poscar=poscar,atomindices=1:4) # will isolate the first 4 atoms
newposcar <- poscar.createsupercell(poscar=poscar,super=diag(c(2,2,1)) # will create a 2x2x1 supercell


## Slab manipulation
indices <- poscar.getatomlayerindices(poscar=poscar,layers=5)
poscar <- poscar.extractlayers(poscar,layer=1:2,layers=5) # will isolate first 2 layers
poscar <- poscar.removelayers(poscar,layer=1:2,layers=5) # will remove first 2 layers
translationvector <- poscar.getlayertranslation(poscar,fromlayer=1,tolayer=2,layers=5) # will give the 2d offset (x,y) of two layers
poscar <- poscar.translatelayers(poscar,layer=1,layers=5,directtranslation=translationvector)
poscar <- poscar.rotatelayer.deg(poscar,layer=1,layers=5,angle=90)
poscar <- poscar.rotatelayer.rad(poscar,layer=1,layers=5,angle=pi)
poscar <- poscar.mirrorlayers(poscar,layers=1,baselayer=2,mirrorlayers=1) 


## Plotting
poscar <- poscar.getshapedposcar(poscar,20,20)
plot(poscar)
plot.poscar.addbasis(poscar)
plot.poscar.addunitcell(poscar)
plot.atoms.addpositions(atoms)
plot.atoms.adddistance(atoms)
plot.atoms.addnumbers(atoms)
plot.atoms.addarrows<(atomsold,atomsnew)


## Slab plotting
plot.poscar.addlayers(poscar,layer=1:2,layers=5)
plot.poscar.addlayerdistance(poscar,layer=1:2,layers=5) # for side view
plot.poscar.addnumbers(poscar,layers=1:2,layer=5) # atom numbers


## Brillouinzone
bz <- basis.getbrillouinzone(poscar$basis)
bz <- reciprocalbasis.getbrillouinzone(recbasis)
bztype <- reciprocalbasis.getbrillouinzonetype(recbasis)
bzs <- list(bz1,bz2,bz3)
bzs <- poscar.getbrillouinzones(poscar)

plot.brillouinzones(bzs)
plot.brillouinzones.add<-function(bzs)
plot.brillouinzones.addsympoints(brillouinzones,labels=expression(Gamma)

newkpoints <- brillouinzone.projectkpoints(bz,oldkpoints)
newkpoints <- brillouinzone.selectkpoints(bz,oldkpoints)
newkpoints <- brillouinzone.extendkpoints(bz,oldkpoints)
```

### CHGCAR based

```S
## Basics
chgcar <- read.chgcar("CHGCAR")
print(chgcar)
chgcar.sumoverlayer(chgcar,layer=1:2,layers=5) 
chgcar.sum(chgcar) # should give electron count


## STM
stm <- stm(chgcar,emax=0.01,cpus=1,interpolation="linear")
plot.stm(stm,super=4)
plot.stm.addatoms(stm,atomselector=1:10,super=4)
plot.stm.addunitcell(stm)
```
### Vasprun.xml based

```S
## Basics for bandstructures
bandsdata <- read.bandsdata("vasprun.xml")
print(bandsdata)
print(bandsdata$band1)
bandsdata <- bandsdata.addsympoint(bandsdata,c(40,80))
bandsdata <- bandsdata.calcsympointpath(bandsdata,sympointpath=list(c(1,2),c(3,4)))
xlim <- bandsdata.getintervallaroundsympoint(bandsdata,sympointnumber=2)
energydistance <- bandsdata.getbanddistance(bandsdata,kpoint=40,bands=c(4,6))
energy <- bandsdata.getenergy(bandsdata,kpoint=40,band=4)


## Plotting bandstructures
bandsdata <- plot.bandsdata(bandsdata,sympointpath=list(c(1,2),c(3,4)),energyoffset=-1)
plot.bandsdata.addsymnnames(bandsdata,labels=c(expression(Gamma),"M","K",expression(Gamma)))
plot.bandsdata.addfermi(bandsdata)
plot.bandsdata.addbands(bandsdata,bands=1:10)
plot.bandsdata.addnumbers(bandsdata)


## Projected bands data
projecteddata <- bandsdata.getprojecteddata(bandsdata)
plot.projectedbands.add(projecteddata)


## Fitting bandstructure
sp <- bandsdata.fit.dirac.makeparameters(vF=10)
sp <- bandsdata.fit.quadratic.makeparameters(m=10)
ap <- bandsdata.fit.dirac.makeconstants(bandsdata,bands=4:5,k0=40,factor=-1)
ap <- bandsdata.fit.quadratic.makeconstants(bandsdata,band=4,k0=40)
bfit <- bandsdata.fit(bandsdata,bandnr=4,kpoints=36:40,fitname="dirac",startingparameters=sp,constants=ap)
fitfunc <- bandsdata.fit.dirac.function(parameters=sp,kpoints=36:40,constants=ap)
fitfunc <- bandsdata.fit.quadratic.function(parameters=sp,kpoints=36:40,constants=ap)
predicteddata <- predict(bfit,seq(bfit$kdist[[1]],bfit$kdist[[length(bfit$kdist)]],length.out=101))
bandsfit.dirac.getv(bfit)
bandsfit.quadratic.getm(bfit)
print(bfit)
plot.bandsfit.add(bfit,energyoffset=-1)


## Plotting 3d bandstructure
plot.bandsdata.grid(bandsdata)
plot.bandsdata.contour(bandsdata,band=4)
plot.bandsdata.3d(bandsdata,bands=4:6,projected=T)


## Basics DOSdata
dosdata <- read.dosdata("vasprun.xml") 
print(dosdata)
dosdata <- dosdata.addsmearing(dosdata)
dosvector <- dosvector.calcsmearing(energy=energyvector,dos=dosvector)


## Plotting DOS
dosdata <- plot.dosdata(dosdata,smearing=0.1,flip=T,fermi=T)
dosdata <- plot.dosdata.add(dosdata,type="polygon",smearing=0.1,orbitals=c(1,2,3,4,"all"),atomindices=1:2)
plot.dosdata.addfermi(dosdata)
```

### Calculation based

```S
## Basics
calculation <- read.calculations(folders)
print(calculations)
print(calculation)
load.calculations(file="filename")


## E over a
ea <- calculation.getea(calculation)
fiteos <- ea.fitEOS(ea)
data <- predict(fiteos,newa)
plot.calculation.ea(calculation)
plot.calculation.ea.addpoints(calculation)
plot.calculation.ea.addfit(calculation)
plot.EOS.add(fiteos)


## Bulkbands
bulkbands <- calculation.getbulkbands(calculation)
plot.bulkbands.add(bulkbands)


## LDOS
position <- poscar.getpositionbyatom(poscar,1)
position <- poscar.getpositionbylayer(poscar,layers=5,layer=1)
plot.calculation.ldos(calculation,positions)
plot.calculation.addldos(calculation,positions)
plot.positions.addlegend(position)

indices <- calculation.getindicesfrompositions(calculation,positions)
data <- calculation.getdatafromindices(calculation,indices)
ldosvector <- ldosvector.calcsmearing(energyvector,ldosvector)
```

### miscellaneous
```S
## Periodic table
plot.periodictable(highlights=c("H","He"),highlighttexts=list(1,2),underlay=c("Si"))
plot.periodictable.addelements(c("H","He"),texts=list(1,2)) # usefull if first plot with typ="n"
plot.periodictable.addunderlay(c("Si"))

elementnames <- periodictable.getelementnames(element=1:2)
elementselector <- periodictable.getelementselector(element=c("H","He"))
remainingelements <- periodictable.getremainingelements(element=3:109)
elementpositions <- periodictable.getelementpositions(element=c("H","He"))


## other
plot.addlabel("(a)")
df <- dataframe.applysymoperations(df,symoperations=list(c("rotation",180)))
newfilename <- file.gethighestversion(dir,filename)
vector <- vectors.crossproduct(vector1,vector2)
length <- vector.length(vector)
angle <- vectors.calcangle(vector1,vector2)
angle <- vectors.calcangle.degree(vector1,vector2)
matrix <- matrix.rotation2d(angle)
matrix <- matrix.rotation2d.degree(angle)
matrix <- matrix.reflection2d.degree(slope)
matrix <- matrix.reflection2d(slope)
newcolor <- makeTransparent(color, alpha=100)

fig <- getplotinplotfig(x,y,coords="ndc")
plotinplot(x,y)
endplotinplot()
zoomplot(x,y,xmean,xwidth,ymean,ywidth)
endzoomplot()
range <- calcrange(listofobjects)
```