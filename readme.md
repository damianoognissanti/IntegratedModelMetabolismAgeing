README
===============
**The repository is part of publication: Multi-scale model suggests the trade-off between protein and ATP demand as a driver of metabolic changes during yeast replicative ageing, Schnitzer and Österberg et al., 2022** (https://www.biorxiv.org/content/10.1101/2022.03.07.483339v1)


**Abstract**

The accumulation of protein damage is one of the major drivers of replicative ageing, describing a cell’s reduced ability to reproduce over time even under optimal conditions. Reactive oxygen and nitrogen species are precursors of protein damage and therefore tightly linked to ageing. At the same time, they are an inevitable by-product of the cell’s metabolism. Cells are able to sense high levels of reactive oxygen and nitrogen species and can subsequently adapt their metabolism through gene regulation to slow down damage accumulation. However, the older or damaged a cell is the less flexibility it has to allocate enzymes across the metabolic network, forcing further adaptions in the metabolism.
To investigate changes in the metabolism during replicative ageing, we developed an multi-scale mathematical model using budding yeast as a model organism. The model consists of three interconnected modules: a Boolean model of the signalling network, an enzyme-constrained flux balance model of the central carbon metabolism and a dynamic model of growth and protein damage accumulation with discrete cell divisions.
The model can explain known features of replicative ageing, like average lifespan and increase in generation time during successive division, in yeast wildtype cells by a decreasing pool of functional enzymes and an increasing energy demand for maintenance. We further used the model to identify three consecutive metabolic phases, that a cell can undergo during its life, and their influence on the replicative potential, and proposed an intervention span for lifespan control.


## Structure of the repository

- **ExperimentalData** contains the experimental data that was used for the calibration and validation of the model.

- **Functions** contains all functions needed to run the model, as well as a schematic picture of the code structure. To run a lifespan simluation, only the file "IntegratedModel.jl" is needed to be imported.

- **ModelFiles** contains model files for the three modules of the integrated model: Boolean model for the signalling, FBA model for the central carbon metabolism and ODE model for the dynamic damage accumulation.

- **Simulations** contains all simulations relevant to the publication, sorted by the figure numbers. The file type defines the purpose of the file (executable or parameter definitions: .jl, result: .txt, plotting function .R, plots .svg). For a simple individual lifespan simulation see Fig2B. 

Each folder contains another readme file for more details.

## Required software

The simulations were tested with the **Julia** programming language, version 1.6.1, including important packages MAT (reading in FBA model), JuMP (optimisation toolbox) and Gurobi (solver for linear programs).
