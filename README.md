Rvasp
=====
Tools for loading, manipulating and plotting VASP files within R

### Install
* download [Rvasp package](../../raw/master/Rvasp_0.1.tar.gz)
* install in R

```R
install.packages(pkgs="Rvasp_0.1.tar.gz")
```

* load in R

```R
library(Rvasp)
```

Contains functions to work with

### [POSCAR](../../wiki/POSCAR)
* [read](../../wiki/Basics-(POSCAR\)) 
* [write](../../wiki/Basics-(POSCAR\)) 
* [manipulate](../../wiki/Manipulation-(POSCAR\))
* plot

### CHGCAR 
* read 
* calculate and plot stm

### vasprun.xml 
* read and plot (projected) bandsdata and fit a dirac cone or a quadratic function to a band
* read and plot (projected) dosdata 

Contains a wrapper for calculations organized in the following scheme:   
   
/manyCalculations/Calculation1/Parameter1/VASPfile1   
/manyCalculations/Calculation1/Parameter1/VASPfile2   
..   
/manyCalculations/Calculation1/Parameter2/VASPfile1   
/manyCalculations/Calculation1/Parameter2/VASPfile2   
..   
/manyCalculations/Calculation2/Parameter1/VASPfile1   
/manyCalculations/Calculation2/Parameter1/VASPfile2   
etc.

so there are three levels:
* Calculations
* Calculation
* Parameter

On calculation level there are the following functions:
* get / plot / fit an e over a curve
* get / plot a local dos curve
* get / plot bulk band data