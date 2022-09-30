# Windows_Rjml 
This is Rjml windows version   
Rjml is a program that implements the JML software function (detecting hybridization sequences) in R language.   
RJml can use multiple threads to increase its speed compared to JML.  
## 1. Package summary
JML program uses posterior predictive checking to test whether the minimum distance between sequences of two species is smaller than expected under a scenario that does not account for hybridization.
If the observed distance is smaller than, the simulated valued, then we can conclude that lineage sorting cannot explain the data and a hypothesis of hybridization can be accepted (Joly 2012).
In order to implement the above procedure, the program needs to perform the following steps:  
(1) Posterior distributions of species trees, population sizes, and divergence times generated from other software such as BEAST and BPP.  
(2) JML samples each species tree to simulate gene trees through MCMCCOAL (Rannala & Yang, 2003).  
(3) Sequences are then simulated on the gene tree. This was implemented by adapting the code of the software SEQ-GEN (Rambaut & Grassly 1997).   
(4) JML estimates the minimum distance between sequences for all pairs of species.  
This will generate posterior predictive distributions that can then be used to estimate the P-value of empirical minimum sequences distances between species.
To implement all the functionalities of JML in Rjml, we first modified the source code of JML so that it can only simulate DNA sequences based on the posterior distribution of the input species tree, population size, and divergence time.
The modified software, named jml-sim, is used in exactly the same way as JML, but can only read species trees and output simulated DNA sequences. 
After that, the minimum distance between sequences of all species pairs is estimated and the P-value posterior predictive distributions are calculated using R-based code, Rjml.R, working in parallel, producing in a result file in the same format as JML. 
Rjml can afford users four-fold reductions in time over a traditional JML analysis in detecting hybridization. 
For convenience, we assemble all the necessary documents (jml-sim code, R functions and associated packages) into a Windows pipeline, winRjml.bat, for Windows users. 
The user only needs to specify the input sequences, the control file, posterior species trees (include population size and divergence time), the significance level of the P-value, and the number of threads to obtain the result of hybridization detection.
## 2. Installation of Rjml on Windos system
The prerequisite of installing and using Rjml without troubles, regardless of operation system used, is to have R software installed in users’ computers. 
Then, after downloading and unzipping the compressed files: Rjml-windows.rar, users can use a batch scrip file inside the uncompressed folder to finish the installation of Rjml in Windows system. 
The uncompressed folder contains all the necessary files for Rjml installation and computation, including:  
(1) batch script file :winRjml.bat  
(2) Pre-compiled Windows executable file (jml-sim.exe, libgcc_s_seh-1.dll and libstdc++-6.dll)  
(3) R packages argparser, foreach, doparallel, iterators, Rcpp, RcppArmadillo and Cjml1.  
(4) A example dataset include three files：jml.tpi.ctl rosa.species.trees tpi.phy.  
After decompressing Rjml-windows.rar file, no specific installation was required.  
Users can double click the batch file winRjml.bat or open CMD terminal and type the calculation command as below to finish the installation and analyses automatically (users should make sure all the uncompressed files in the same directory by decompressing Rjml-windows.rar files).
##3.Windows usage of Rjml  
We demonstrate Rjml execution on a Windows system using the example data provided by JML. 
This dataset includes three files: rosa. species. trees, tpi. phy and jml. tpi. ctl.   
The rosa.species.trees file contains a posterior distribution of 1000 species trees (with divergence times and population sizes) of 8 species of Rosa genus  
The tpi. phy file contains 92 sequences of the 8 species of Rosa genus (810 bp in length).  
The jml.tpi.ctl file is a control file that provides necessary information for DNA sequence simulation.  
More information about the creation and formatting of the control file, please refer to the usage document of JML (https://github.com/simjoly/jml).  
In fact, Rjml runs in two steps. 
The first step requires generating the simulation sequence using jml-sim software.  
After generating the simulation sequence, users are expected to use Rjml.R calculate the posterior predictive distributions for the uncorrected minimum distance between sequences for all pair of species.   
However, in our Windows implementation of Rjml, we provide batch scripts for helping users to implement all the simulations and calculations simultaneously with a single command line.  
Users should open the CMD or PowerShell program, enter the folder where the Rjml-windows.rar are unzipped.  
Then, for using Rosa example as a demonstration, one can type the following command:  
~$winRjml jml.tpi.ctl rosa.species.trees tpi.phy 0.1 8~
