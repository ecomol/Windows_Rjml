::Rjml. R and the input file must be in the same path, otherwise the input file must have a path added
::Enter separately in turn control filel,species tree,sequence,significance level and thread 

@echo off
setlocal enabledelayedexpansion
echo 		******************************************************************************************
echo 		*                                                                                        *
echo 		*                    This is Windows version of Rjml                                     *
echo 		*       a fast program for detecting hybridization from species trees                    *
echo 		*   A Windows command example: winRjml jml.tpi.ctl rosa.species.trees tpi.phy 0.1 8      *
echo 		*                                                                                        *
echo 		*                         Coded by Chen's Lab at CIB                                     *
echo 		*                                                                                        *
echo 		******************************************************************************************
echo. 
set /p var=Please input the command for analysis:
echo.
set /A n=0
for %%a in (%var%) do (
set c[!n!]=%%a
set /A n+=1
Rem echo %%a
)
set R=Rjml.R
.\jml-sim.exe -c %c[1]% -t %c[2]%
Rscript.exe %R% %c[3]% %c[4]% %c[5]%