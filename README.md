Rvasp
=====
Tools for loading, manipulating and plotting VASP files within R

### Workflow
To be somewhat efficent, even in research, I developed a workflow for scientific computation, which works quite well for me.
![scientific computation - workflow](/raw/master/examples/workflow-01.png "scientific computation - workflow")
VASPmanager can be found [here](https://github.com/paulpflug/VASPmanager) 

### Important
Many features originated by a need to solve a specific problem and are propably not generalized enough to declare the general problem fully solved.
All functions are implemented in pure R and are mostly easy to understand. Feel free to read the [source](../../tree/master/Rvasp/R), improve and submit changes.

## Index

- [Install](#install)
- [What it does](#what-it-does)
- [POSCAR](#poscar-1)
  - [Basics](#basics)
  - [Manipulation](#manipulation)
    - [Raw manipulation](#raw-manipulation)
    - [General manipulation](#general-manipulation)
    - [Slab manipulation](#slab-manipulation)
  - [Plotting](#plotting)
    - [General plotting](#general-plotting)
    - [Slab plotting](#slab-plotting)
- [CHGCAR](#chgcar-1)
  - [STM](#stm)
- [Vasprun.xml](#vasprunxml-1)
  - [BANDS and DOS](#bands-and-dos)
- [Calculation](#calculation)
  - [E over a](#e-over-a)
  - [Local dos](#local-dos)
  - [Bulk bands](#bulk-bands)
- [Miscellaneous](#miscellaneous)
  - [Brillouin zone](#brillouin-zone)
  - [Periodic table](#periodic-table)
  
## Install
##### Install as a package
* download [Rvasp package](https://www.dropbox.com/s/gtuz0o79zurro6a/Rvasp_0.2.tar.gz)
* install in R

```R
install.packages(pkgs="Rvasp_0.2.tar.gz",repos=NULL, type="source")
```

* load in R

```R
library(Rvasp)
```

##### Use the sourcecode
* clone or download & extract zip from the right
* open `Rvasp.Proj` in [RStudio](http://www.rstudio.com/)
* build

## What it does
Contains functions to work with POSCAR, CHGCAR, vasprun.xml and more

##### POSCAR
* [read](#read) 
* [write](#write) 
* [manipulate](#manipulation)
* [plot](#plotting)

##### CHGCAR
* [read](#read-1)
* calculate and plot [stm](#stm)

##### vasprun.xml 
* read and [plot](../../wiki/example-Plots#dos) (projected) [dosdata](../../wiki/DOS)
* read and [plot](../../wiki/example-Plots#bands) (projected) [bandsdata](../../wiki/BANDS) and fit a dirac cone or a quadratic function to a band

#### more
Contains a wrapper for [calculations](#calculation) organized in the following scheme:   

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
* get / [plot](#e-over-a) / fit an e over a curve
* get / [plot](#local-dos) a local dos curve
* get / [plot](#bulk-bands) bulk band data

Furthermore there are some [miscellaneous](#miscellaneous) functions like
plotting a [periodic table](#periodic-table)


## POSCAR

### Basics

##### Read
POSCAR like files can be read by using `poscar <- read.poscar(file="POSCAR")`

##### Printing
There are two versions of print. `print(poscar)` will produce output similar to this:

```R
###POSCAR###
   firstline Argentum
   a 4.03
   basis length:9 (matrix)
   eightline Selective dynamics
   ninthline Direct
   selectivedynamics TRUE
   atoms length:7 (data.frame)
   length 10
###########
```

by using `printraw.poscar(poscar)` you'll get the ascii version of the file

```R
Argentum
4.03
0.0000000000 0.5000000000 0.5000000000
0.5000000000 0.0000000000 0.5000000000
0.5000000000 0.5000000000 0.0000000000
Ag
1
Selective dynamics
Direct
0.0000000000 0.0000000000 0.0000000000 T T T
```

##### Write
To write the poscar to file use

```R
write.poscar(poscar=poscar,file="POSCAR")
```

### Manipulation
#### Raw manipulation
If your need is very elevated you have to fall back to raw manipulation.
Most important are the following objects

```R
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

With R built-in tools this can be quite powerfull

```R
poscar$atoms$rz[poscar$atoms$type=="Ag"]
```

or

```R
poscar$atoms[poscar$atoms$type=="Ag",3]
```

for example will give you all z components of all silver atoms.

For manipulation, changes have to be saved,

```R
poscar$atoms[poscar$atoms$type=="Ag",3]<-poscar$atoms[poscar$atoms$type=="Ag",3]+0.2
```

will increase all third coordinates of all silver atoms by 0.2 .

#### General manipulation

##### Basis conversion

To get the atoms in cartesian coordinates use:

```R
atoms <- poscar.getbasisconvertedatoms(poscar=poscar)
```

The inversion works like this:

```R
poscar <- atoms.convertbasis(atoms=atoms,basis=poscar$basis*poscar$a)
```

##### Basis rotating

Mainly used to align the first vector with the x direction.

```R
poscar <- poscar.rotate2d(poscar)
```

will rotate `poscar$basis[1:2,1:2]` so that `poscar$basis[1,2]` equals 0

##### Get reciprocal basis

Rarely needed but quite usefull:

```R
recbasis <- poscar.getreciprocalbasis(poscar)
```

will calculate the reciprocal basis.

##### Setting vacuum

The current vacuum (distance of atoms in basis direction over boundary of the unit cell) can be read by

```R
vacuum <- poscar.getvacuum(poscar)
```

and set by

```R
poscar <- poscar.setvacuum(poscar,vacuum=c(0,0,10),center=T)
```

This command will set the vacuum in z direction to 10 Å and center the atoms in the unitcell in all three directions.

##### Center atoms

To center atoms you can use:

```R
poscar <- poscar.centeratoms(poscar,direction=3,position=0.7)
```

This will move the center of all atoms in z-direction (default: all three directions) to the position of 0.7 (default 0.5)

##### Sorting atoms

Sorting is only for readability. It can be accomplished by

```R
poscar <- poscar.sortatoms(poscar,sortindices=c(1,2))
```

which will sort by x-direction followed by a sort in y-direction. (default c(7,3,1,2) equals sort by atom type, z, x and y)

##### Extracting atoms

For simple extraction of certain atoms followed by a optional vacuum change and optional centering you can use

```R
newposcar <- poscar.extractatoms(poscar=poscar,atomindices=1:4,vacuum=c(0,0,15),center=T)
```

which will produce a new poscar with the first four atoms of poscar and a vacuum in z-direction of 15 Å

##### Creating supercell

To create a supercell simply use:

```R
newposcar <- poscar.createsupercell(poscar=poscar,super=c(2,2,1),center=T,center.directions=3,center.position=0.5)
```

this will also center the newly created cell in z-direction. The default chained centering can be disabled by `center=F`

#### Slab manipulation

```R
poscar.getatomlayerindices(poscar=poscar,layers=12)

poscar.extractlayers <- function(poscar,layer,layers,vacuum=c(0,0,0),center=T)
poscar.removelayers <- function(poscar,layer,layers,vacuum=c(0,0,0),center=T)
poscar.getlayertranslation<-function(poscar,fromlayer,tolayer,layers)
poscar.translatelayers <- function(poscar,layer,layers,directtranslation)
poscar.rotatelayer.deg<-function(poscar,layer,layers,angle,offset=NULL)
poscar.rotatelayer.rad
poscar.mirrorlayers<-function(poscar,layers,baselayer,mirrorlayers)
```

### Plotting

```R
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

#### General plotting

```R
plot.poscar <- function(poscar,direction=3,xlab="x",ylab="y",basis=F,unitcell=F,fullcell=F,...)

poscar.getshapedposcar <- function(poscar  # Input poscar
                                  ,x  # Max x (range will be c(0,x))
                                  ,y  # Max y (range will be c(0,y))
                                  ,shape=c("rectangular")  # Shape of area atoms are allowed
                                  )

plot.poscar.addbasis<-function(poscar
                              ,xoffin=0.1
                              ,yoffin=0.1
                              ,direction=3
                              ,basisnames=c("a","b","c")
                              ,arrowlength=0
                              ,arrowsize=0.1
                              ,fullcell=F
                              ,...)

plot.poscar.addunitcell<-function(poscar,direction=3,lty=2,...)

plot.atoms.addpositions<-function(atoms,basis=NULL,direction=3,col="white",cex=3,lwd=1,lty=1,...)

plot.atoms.adddistance<-function(atoms
                                 ,basis=NULL
                                 ,direction=3
                                 ,length=0.1
                                 ,lwd=1
                                 ,col="black"
                                 ,selectedalpha=150
                                 ,selectedsize=1.5
                                 ,...)

plot.atoms.addnumbers<-function(atoms,basis=NULL,direction=3,atomselector=NULL,...)


plot.atoms.addarrows<-function(atomsold,atomsnew,basisold=NULL,basisnew=NULL,direction=3,length=0.1,...)

```

#### Slab plotting

```R
plot.poscar.addlayers<-function(poscar
                                ,layer
                                ,layers
                                ,direction=3
                                ,color=rainbow(length(layer))
                                ,size=rep(1,length(layer))
                                ,lwd=rep(1,length(layer)))

plot.poscar.addlayerdistance<-function(poscar,layer,layers,color=rainbow(length(layer)),length=0.05,direction=1,...)

plot.poscar.addnumbers<-function(poscar,layers=1,layer=1,direction=3,absolutenumber=T,...)
```
## CHGCAR

### read

### STM

```R
data(silverstm)
plot.stm(silverstm,super=5,xlab="x (nm)",ylab="y (nm)")
plot.stm.addunitcell(silverstm,atomnumber=1)
plot.stm.addatoms(silverstm,super=5,xlim=c(0,0.75),ylim=c(0,0.75),atomselector=1,atomsize=10,col="red")
plot.stm.addatoms(silverstm,super=5,xlim=c(0,0.75),ylim=c(0,0.75),atomselector=2,atomsize=10,col="blue")
```
![STM picture](../../raw/master/examples/stm.png "STM of a silver surface created in R")

## Vasprun.xml

### BANDS and DOS

```R
data(silverbands)
proj <- bandsdata.getprojecteddata(silverbands)
plot.bandsdata(silverbands,symnames=c(expression(Gamma),"L","W","X",expression(Gamma),"K","X"),fermi=T,ylim=c(-10,20))
plot.projectedbands.add(proj,orbitals=c(1,2,3,4),cex=1,legendcex=1.2)
data(silverdos)
dosdata <- plot.dosdata(silverdos,flip=T,fermi=T,col.fermi="Blue",xlim=c(0,0.5),ylim=c(-10,20))
plot.dosdata.add(silverdos,orbitals=c(1:4,"all"),type="polygon",col=c(colorRampPalette(c("red","blue","green"))(4),"grey"))
```
![BANDS picture](../../raw/master/examples/bands.png "Bands of a silver created in R")
![DOS picture](../../raw/master/examples/dos.png "DoS of a silver created in R")

## Calculation

### E over a

```R
data(silverea)
plot.calculation.ea(silverea$LDA,energyshift=silverea$LDAsingleAtom$a_20.000$energy,ylim=c(-3.8,-2.1))
plot.calculation.ea.addpoints(silverea$GGA,energyshift=silverea$GGAsingleAtom$a_20.000$energy,pch=3)
```
![E over a plot](../../raw/master/examples/eaplot.png "E over a plot created in R")
### Local dos
### Bulk bands
Download silverbulkbands.RData from examples

```R
load(silverbulkbands.RData)
bandsdata <- silverbulkbands$slab$a_4.030$banddata
bulkbands <- calculation.getbulkbands(silverbulkbands$bulk)
plot.bandsdata(bandsdata,type="n",symlty=NA,fermi=F,ylim=c(-10,10))
plot.bulkbands.add(bulkbands)
plot.bandsdata.addbands(bandsdata=bandsdata,col="black")
plot.bandsdata.addfermi(bandsdata)
plot.bandsdata.addsymnnames(bandsdata,symnames=c(expression(Gamma),"M","K",expression(Gamma)))
```
![Bulk bands](../../raw/master/examples/bulkbands.png "Silver slab bands underlayed by silver bulk bands, created in R")

## Miscellaneous
### Brillouin zone

```R
data(silverslab)
brillouinzone<-poscar.getbrillouinzone.hexagonal(silverslab)
plot(brilloinzone)
silverslab$a <- silverslab$a*sqrt(3)
brilloinzonesqrt3<-poscar.getbrillouinzone.hexagonal(silverslab,extend=2,rotate=30)
plot.brillouinzone.add(brilloinzonesqrt3,col="red")
plot.brillouinzone.addsympoints(brillouinzone=brilloinzone,directcoordinates=list(c(1,0),c(0,0),c(1/2,1/2)),labels=c("K",expression(Gamma),"M"))
plot.brillouinzone.addsympoints(brillouinzone=brilloinzone,directcoordinates=list(c(1,0),c(0,0),c(1/2,1/2)),labels=c(expression(Gamma),expression(Gamma),"M"),textpos=4,col="red")
```
![Brillouinzone](../../raw/master/examples/brillouinzone.png "Brillouinzone of 1x1 and sqrt(3)xsqrt(3) silver, generated by R")
### Periodic table

```R
plot.periodictable(highlights=1:2,highlighttexts=list(1.008,4.0026),underlay=3:10)
```
![Periodic table](../../raw/master/examples/periodic_table.png "Periodic table generated by R")
