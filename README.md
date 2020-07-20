# CLM-R
These files can be used to create and simulate a compartmental logical model (CLM) of mitotic exit control, using the BoolNet package. 
File description:
- CreateCLM.R is a function that can be used to create a CLM from activity and localization network files (contained in networkData file).
- CLMshell.R is a shell file that can be used to create a CLM and run some basic simulations. 
- CLMfunctions.R provides a number of functions used in other scripts. 
- extractSSlocdata.R is a function that can be used to extract and visually represent steady state localization data from a model. 
- simulateMutantSync.R provides functions to simulate and find steady states of mutants under a synchronous update scheme.
- simulateMutantAsync.R provides functions to simulate mutants under an asynchronous update scheme.
- model5.txt, model5OE.txt and model5SPB.txt are BoolNet format representations of the CLM of mitotic exit. The OE and SPB variants can be used to simulate overexpression and forced localization at the SPB respectively.

These files can be used to train a logical model of the FEAR network, using CellNOptR:
- CreateFEARPKN.R is a script to create a prior knowledge network in CellNOptR format.
- ReadFEARphenotypes.R is a script to read a list of FEAR phenotypes and create a phenotype list (CNO file) in CellNOptR format. 
- TrainFEARnetwork.R is a script that train the FEAR PKN against the phenotype list using a genetic algorithm (CellNOptR).
- TrainingFuncs.R provides various functions to assist in training the model.

Additional misc files:
- LoadModels.R
- Mutants.R
- SimulateCaydasiModel.R
