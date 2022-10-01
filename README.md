# Windows_Rjml 
This is **Rjml** windows version.   
**Rjml** is a program that implements the **JML** software function (detecting hybridization sequences) in R language.   
**Rjml** uses the strategy of R/C++ integration, array programming and multiple threads to increase its speed compared to **JML**. The computational speed can be two to five folds faster than **JML**, thus being capable to handle large data sets.
## 1. Package summary
**JML** program uses posterior predictive checking to test whether the minimum distance between sequences of two species is smaller than expected under a scenario that does not account for hybridization.
If the observed distance is smaller than, the simulated valued, then we can conclude that lineage sorting cannot explain the data and a hypothesis of hybridization can be accepted (Joly 2012).  
In order to implement the above procedure, the program needs to perform the following steps:  
(1) Posterior distributions of species trees, population sizes, and divergence times generated from other software such as BEAST and BPP (Drummond & Rambaut et al., 2007; Flouri et al., 2018).  
(2) **JML** samples each species tree to simulate gene trees through MCMCCOAL (Rannala & Yang, 2003).  
(3) Sequences are then simulated on the gene tree. This was implemented by adapting the code of the software SEQ-GEN (Rambaut & Grassly 1997).    
(4) **JML** estimates the minimum distance between sequences for all pairs of species. This will generate posterior predictive distributions that can then be used to estimate the P-value of empirical minimum sequences distances between species.  
To implement all the functionalities of **JML** in **Rjml**, we first modified the source code of **JML** so that it can only simulate DNA sequences based on the posterior distribution of the input species tree, population size, and divergence time.
The modified software, named **jml-sim**, is used in exactly the same way as **JML**, but can only read species trees and output simulated DNA sequences. 
After that, the minimum distance between sequences of all species pairs is estimated and the P-value posterior predictive distributions are calculated using R-based code, **Rjml.R**, producing in a result file in the same format as **JML**.  
The computational speed of **Rjml** is much faster in comparison to **JML**. To achieve this, we took advantage of C++ seamless integration through Rcpp, parallel computation and array programming in R. Rcpp functions were developed to read sequences and calculate genetic distances between individuals of species and then calculate minimal genetic distances between species. In addition, R packages DOPARALLEL (Microsoft Corporation & Weston, 2017) and FOREACH (Calaway et al., 2015) were used for reading simulated DNA sequence files and compute genetic distances between individuals of species in parallel. These practices can greatly speed up the computation and detection of introgressed sequences between species from many posterior species trees. Our tests show that **Rjml** can afford users **four-fold reductions in time** over a traditional **JML** analysis in detecting hybridization. 
## 2. Installation of **Rjml** on Windos system
The prerequisite of installing and using **Rjml** without troubles, regardless of operation system used, is to have R software installed in users’ computers. 
Then, after downloading and unzipping the compressed files: **Rjml-windows.rar**, users can use a batch scrip file inside the uncompressed folder to finish the installation of **Rjml** in Windows system. 
The uncompressed folder contains all the necessary files for **Rjml** installation and computation, including:  
(1) batch script file: **winRjml.bat**  
(2) Pre-compiled Windows executable file :**jml-sim.exe**, **libgcc_s_seh-1.dll** and **libstdc++-6.dll**  
(3) R packages argparser, foreach, doparallel, iterators, Rcpp, RcppArmadillo and Cjml1.  
(4) A example dataset include three files：**jml.tpi.ctl** **rosa.species.trees** **tpi.phy**.  
(5) **Rjml.R**: R code that performs the fourth step of **JML** analysis.  
After decompressing **Rjml-windows.rar** file, no specific installation was required.  
Users can double click the batch file **winRjml.bat** or open CMD terminal and type the calculation command as below to finish the installation and analyses automatically (users should make sure all the uncompressed files in the same directory by decompressing **Rjml-windows.rar** files).
## 3.Usage of **Rjml** on Windows system
We demonstrate **Rjml** execution on a linux system using the example data provided by **JML**. 
This dataset includes three files: **rosa. species. trees**, **tpi. phy** and **jml. tpi. ctl**.   
**rosa.species.trees** contains a posterior distribution of 1000 species trees (with divergence times and population sizes) of 8 species of Rosa genus.   
The **tpi. phy** file contains 92 sequences of the 8 species of Rosa genus (810 bp in length).  
The **jml.tpi.ctl** file is a control file that provides necessary information for DNA sequence simulation.  
More information about the creation and formatting of the control file, please refer to the usage document of **JML** (https://github.com/simjoly/jml).  
In fact, **Rjml** runs in two steps. 
The first step requires generating the simulation sequence using **jml-sim** software.  
After generating the simulation sequence, users are expected to use **Rjml.R** calculate the posterior predictive distributions for the uncorrected minimum distance between sequences for all pair of species.   
However, in our Windows implementation of **Rjml**, we provide batch scripts for helping users to implement all the simulations and calculations simultaneously with a single command line.  
First, users should open the CMD or PowerShell program, enter the folder where the **Rjml-windows.rar** are unzipped.  
Then, for using Rosa example as a demonstration, one can type the following command:  
`winRjml jml.tpi.ctl rosa.species.trees tpi.phy 0.1 8`    
Alternatively, users can simply double click **winRjml.bat** file to type the above command for implementing Rosa example.
The three files used in here was the: jml-sim control file, the posterior species tree file, the original sequence file.  
The last two numbers represent the P-value and thread number for parallel computation.  
It is noted here, Rscript.exe is recommended to be listed in environment variable when it was used. To do so, users should find the path to **R.exe** or **Rscript.exe** on their own computers and add the path into the environmental variables.  
## 4. Output files
During the calculation, users may observe that some intermediate output files have been generated, that is, the generated simulation sequences of **jml-sim** with the file name prefixed with “RepSeqs” followed by a number. It also produces procedural files called out.trees and rep.dat, which are repeatedly written and read during the simulation. User should not remove them while the program is running, and all of "RepSeqs” files will be removed automatically when all the required simulation and calculation are done.  
**Rjml** produces three result files:  
**Distributions_Rjml.csv**, **Probabilities_Rjml.csv** and **Results_Rjml.csv**.  
They correspond to the three files in the original **JML** software: **Distributions.txt**, **Probabilities.txt** and **Results.txt**.  
**Distributions_Rjml.csv** file contains the posterior predictive distributions, which can be used to calculate the P-value of hybridization hypotheses.  
**The Probabilities_Rjml.csv** file contains the probability of observing the minimum distance between all pair of species according to a scenario without hybridization.  
**The Results_Rjml.csv** contains all the pairwise sequence distances that have a P-value smaller than the significance level specified in the control file. Note that if the file is empty or unavailable, this means that no hybridization events were found with a P-value below the threshold.
## 5. References
Calaway, R., Weston, S., & Calaway, M. R. (2015). Package ‘foreach’. Retrieved from https://CRAN.R-project.org/package=foreach.  
Corporation, M., & Weston, S. (2017). doParallel: Foreach parallel adaptor for the “parallel” package. Retrieved from https://CRAN.R-project.org/package=doParallel.  
Drummond, A. J., & Rambaut, A. (2007). BEAST: Bayesian evolutionary analysis by sampling trees. BMC evolutionary biology, 7(1), 1-8.  
Flouri, T., Jiao, X., Rannala, B., & Yang, Z. (2018). Species tree inference with BPP using genomic sequences and the multispecies coalescent. Molecular biology and evolution, 35(10), 2585-2593.  
Joly, S. (2012). JML: testing hybridization from species trees. Molecular Ecology Resources, 12(1), 179-184.  
Rannala, B., & Yang, Z. (2003). Bayes estimation of species divergence times and ancestral population sizes using DNA sequences from multiple loci. Genetics, 164(4), 1645-1656.  
Rambaut, A., & Grass, N. C. (1997). Seq-Gen: an application for the Monte Carlo simulation of DNA sequence evolution along phylogenetic trees. Bioinformatics, 13(3), 235-238.  
****
If you have any questions, please contact:  
Qi Xiao, email: m18033705433@163.com; Address: CAS Key Laboratory of Mountain Ecological Restoration and Bioresource Utilization & Ecological Restoration and Biodiversity Conservation Key Laboratory of Sichuan Province, Chengdu Institute of Biology, Chinese Academy of Sciences, Chengdu, 610041.  
Prof. Youhua Chen, email: haydi@126.com; Address: CAS Key Laboratory of Mountain Ecological Restoration and Bioresource Utilization & Ecological Restoration and Biodiversity Conservation Key Laboratory of Sichuan Province, Chengdu Institute of Biology, Chinese Academy of Sciences, Chengdu, 610041.
