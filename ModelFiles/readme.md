README
===============
**as part of publication: Schnitzer, Österberg et al., 2022**


## This folder contains all individual model files that are used in the integrated model simulation:

**Boolean** includes the original Boolean model of nutrient signalling (Österberg et al., PLOS Computational Biology 2021) in a more readable format that can be processed directly from our program, as well as the modified version that additionally includes oxidative stress signalling (according to Fig 1B, Supplementart Text 1). Transcription factor targets are obtained by the Yeastract database, as well as Lee et al. 1999, and are matched to the enzymes included in the ecFBA model.
OBS: If the model files are opened with Excel it will replace all true/false with TRUE/FALSE which cannot be interpreted from Julia. Make sure to not save those versions or import the files as 'text' not as 'general' to prevent the conversion.

**ecFBA** includes the original reduced enzyme-constraint flux balance model of the central carbon metabolism (Österberg et al., PLOS Computational Biology 2021), as well as a modified version that additionally includes reactive oxygen and nitrogen producing reactions (according to Fig 1A, Supplementart Text 1 and damageReactions.xlsx).

**ODE** includes a simple ODE model for cellular protein damage accumulation (equations (2)-(4)). Since the ODEs are very simple, the solutions are hardcoded in the simulations.
