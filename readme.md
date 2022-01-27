README
===============
**As part of publication: Integrated model of yeast metabolism and replicative ageing reveals importance of metabolic phases during ageing, Schnitzer and Österberg et al., 2022**

**Abstract**




## Structure of the repository

- **ExperimentalData** contains the experimental data that was used for the calibration and validation of the model.

- **Functions** contains all functions needed to run the model, as well as a schematic picture of the code structure. To run a lifespan simluation, only the file "IntegratedModel.jl" is needed to be imported.

- **ModelFiles** contains model files for the three modules of the integrated model: Boolean model for the signalling, FBA model for the central carbon metabolism and ODE model for the dynamic damage accumulation.

- **Simulations** contains all simulations relevant to the publication, sorted by the figure numbers. The file type defines the purpose of the file (executable or parameter definitions: .jl, result: .txt, plotting function .R, plots .svg). For a simple individual lifespan simulation see Fig2B. 

Each folder contains another readme file for more details.

## Required software

The simulations were tested with the **Julia** programming language, version 1.6.1, including important packages MAT (reading in FBA model), JuMP (optimisation toolbox) and Gurobi (solver for linear programs).
