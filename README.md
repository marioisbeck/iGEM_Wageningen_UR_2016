#iGEM_Wageningen_UR_2016
As part of the iGEM (www.igem.org) competition, I created a dynamic model of the optogenetic tools pDusk and pDawn (described by Ohlendorf et al. 2012) and extended it's functionality towards an optogenetic kill switch by implementing the mazEF toxin-antitoxin system. All analysis was conducted with a Mac on OS El Capitan 10.11.6 in Matlab 2015b (student version). Further documentation can be found in the files themselves.

## List of Files
* Data from Ohlendorf et al. (2012) describing their system with lab data can be found in these files: DataPoints_pDawn.csv, DataPoints_pDusk.csv, Standard_Deviation_Ohlendorf_pDusk.csv, Standard_Deviation_Ohlendorf_pDawn.csv
* The equations for the analysis of pDusk and pDawn are described in the pDusk_function.m, pDawn_function.m, pDusk_function_const_mazF.m, and pDawn_function_const_mazE.m files.
* Output files from the pDusk_pDawn_parameter_estimation_weighted_sum.m parameter estimation start with a 10000_0.010_...
* Output files from the pDusk_pDawn_mazEF_parameter_estimation.m parameter estimation start with a 1000_0.010_...

## pDusk/ pDawn Simulation
To reproduce the data shown for pDusk and pDawn on: [BeeT (iGEM 2016) - Optogenetic Kill Switch Model](http://2016.igem.org/Team:Wageningen_UR/Model#light "BeeT (iGEM 2016) - Optogenetic Kill Switch Model") the pDusk_pDawn_parameter_estimation_weighted_sum.m file is the file of your choice. The analysis was already done and to run and see the evaluation of the 10,000 parameter set run, please load and run the pDusk_pDawn_evaluation.m file.

## pDusk + constitutive mazF/ pDawn + constitutive mazE Simulation
Further analysis to implement the mazEF toxin-antitoxin system was done by using the pDusk_pDawn_mazEF_parameter_estimation.m file for the parameter estimation and the pDusk_pDawn_mazEF_evaluation.m can be used to evaluate the 1000 parameter set results.
