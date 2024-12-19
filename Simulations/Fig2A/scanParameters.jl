################################################################################
################################################################################
# simulate the lifespan of a cell with the integrated model 
# (ecFBA for the metabolism, Boolean for the signalling, ODE
# for damage accumulation) for a grid of non-metabolic damage formation and repair
# rates for different regulation factors and retention values
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
parameterCombinationIdx = [[i, j, k] for i = 1:length(formationRate0), 
                                         j = 1:length(repairRate0),
                                         k = 1:length(regulationFactor)]
parameterCombinationIdx = parameterCombinationIdx[:]
nCombinations = size(parameterCombinationIdx, 1)


################################################################################
@printf("Initialise the output ... \n")
################################################################################
outputFile = "f0r0Grid.txt"
output = Array{Any, 2}(undef, 15, nCombinations)


################################################################################
@printf("Go through values for damage formation and repair and simulate cell ... \n")
################################################################################
for n = 1 : nCombinations

    tmpFormationRate = formationRate0[parameterCombinationIdx[n][1]]
    tmpRepairRate = repairRate0[parameterCombinationIdx[n][2]]
    tmpRegulation = regulationFactor[parameterCombinationIdx[n][3]]

    @printf("\tcombination %i/%i\n", n, nCombinations)

    # initialise cell 
    cell = createCell(fbaPath, booleanSpeciesPath, booleanRulesPath, TFpath,
                      nDelay, maxGrowthRate, tmpFormationRate, tmpRepairRate,
                      P, D, mass, [glucoseThreshold, damageThreshold, 
                      trxThreshold])

    # iterate nDelay times to fill up the recent Boolean input states
    simulateLife!(cell, mass, sizeProportion, retention, timestep, nDelay, 
                  tmpRegulation, growthFlexibility)

    # simulate life
    results = simulateLife!(cell, mass, sizeProportion, retention, timestep, 
                            maxTimesteps, tmpRegulation, growthFlexibility)
    status, rls, averageGenerationTime, phaseSplit = results

    # save in output
    output[:, n] = [tmpFormationRate; tmpRepairRate; tmpRegulation;
                    rls; averageGenerationTime; phaseSplit[1]; 
                    phaseSplit[2]; phaseSplit[3]; status]
    println(output[:,n])
end


################################################################################
@printf("Save in output file ...\n")
################################################################################
pretext = 
"# F0-R0 PARAMETER SCAN
# FBA model $fbaPath with objective max. growth (parsimonious)
# Boolean model $booleanRulesPath, $booleanSpeciesPath
# TF layer $TFpath
damage formation\tdamage repair\tregulation factor\trls\t" *
"average generation time\ttime phase 1\trls phase 1\tdamage phase 1\t" *
"time phase 2\trls phase 2\tdamage phase 2\ttime phase 3\trls phase 3\t"*
"damage phase 3\tstatus
"

open(outputFile, "w") do file
    write(file, pretext)
    writedlm(file, permutedims(output), "\t")
end

