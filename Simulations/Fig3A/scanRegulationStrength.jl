################################################################################
################################################################################
# simulate the lifespan of a cell with the integrated model 
# (ecFBA for the metabolism, Boolean for the signalling, ODE
# for damage accumulation) for a range of regulation factors for
# parameter combinations that lead to wildtypes (23 divisions)
################################################################################
################################################################################

using Printf
using DelimitedFiles
using Statistics
include("../../Functions/IntegratedModel.jl")


################################################################################
@printf("Define parameters ... \n")
################################################################################
include("parameters.jl")


################################################################################
@printf("Read in parameter combinations ... \n")
################################################################################
parameters = readdlm(parameterFile, '\t', Any, '\n', skipstart = 5)

# take only the combinations with 23 divisions
parameters = parameters[parameters[:, 4] .== 23, :]
nCombinations = size(parameters, 1)


################################################################################
@printf("Initialise the output ... \n")
################################################################################
outputFile = "regulationFactor.txt"
nRegulations = length(regulationFactor)
output = Array{Any, 2}(undef, 8, nRegulations * nCombinations * 3)


################################################################################
@printf("Simulate cell for different regulation factors ... \n")
################################################################################
for r = 1 : nRegulations
    for n = 1 : nCombinations

        @printf("\tregulation factor %.3f, parameter set %i/%i \n", 
                regulationFactor[r], n, nCombinations)

        # extract parameters
        formationRate0 = parameters[n, 1]
        repairRate0 = parameters[n, 2]

        # initialise cell 
        cell = createCell(fbaPath, booleanSpeciesPath, booleanRulesPath, TFpath,
                          nDelay, maxGrowthRate, formationRate0, repairRate0,
                          P, D, mass, [glucoseThreshold, damageThreshold, 
                          trxThreshold])

        # iterate nDelay times to fill up the recent Boolean input states
        simulateLife!(cell, mass, sizeProportion, retention, timestep, nDelay,  
                      regulationFactor[r], growthFlexibility)

        # simulate life
        results = simulateLife!(cell, mass, sizeProportion, retention, timestep, 
                                maxTimesteps, regulationFactor[r], 
                                growthFlexibility)
        status, rls, averageGenerationTime, phaseSplit = results

        # save in output
        idx = (r - 1) * 3 * nCombinations + (n - 1) * 3
        output[:, idx + 1] = [regulationFactor[r]; "phase 1"; phaseSplit[1]; 
                              formationRate0; repairRate0; status]
        output[:, idx + 2] = [regulationFactor[r]; "phase 2"; phaseSplit[2];
                              formationRate0; repairRate0; status]
        output[:, idx + 3] = [regulationFactor[r]; "phase 3"; phaseSplit[3];
                              formationRate0; repairRate0; status]

    end
end


################################################################################
@printf("Save in output file ...\n")
################################################################################
pretext = 
"# REGULATION FACTOR GRID
# FBA model $fbaPath with objective max. growth (parsimonious)
# Boolean model $booleanRulesPath, $booleanSpeciesPath
# TF layer $TFpath
regulation factor\tphase\ttime\tdivisions\tdamage at end\tdamage formation\t" *
"damage repair\tstatus
"
          
open(outputFile, "w") do file
    write(file, pretext)
    writedlm(file, permutedims(output), "\t")
end

