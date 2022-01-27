################################################################################
################################################################################
# simulate the lifespan of a cell with the integrated model 
# (ecFBA for the metabolism, Boolean for the signalling, ODE
# for damage accumulation) for knockouts of signalling proteins for
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
nKnockouts = length(knockouts)
outputFile = "proteinKnockouts.txt"
output = Array{Any, 2}(undef, 9, (nKnockouts + 1) * nCombinations)

# add wildtype to data
relevantIdx = [4, 5, 7, 10, 13, 1, 2, 15]
wildtypeOutput = [["wildtype" for i = 1:nCombinations] parameters[:, relevantIdx]]
output[:, 1:nCombinations] = permutedims(wildtypeOutput)


################################################################################
@printf("Perturb each enzyme ... \n")
################################################################################
for i = 1 : nKnockouts
    for j = 1 : nCombinations

        @printf("%s %i/%i, parameter set %i/%i\n", knockouts[i], i, nKnockouts, 
                j, nCombinations)

        # extract parameters
        formationRate0 = parameters[j, 1]
        repairRate0 = parameters[j, 2]

        # initialise cell
        cell = createCell(fbaPath, booleanSpeciesPath, booleanRulesPath, TFpath,
                          nDelay, maxGrowthRate, formationRate0, repairRate0,
                          P, D, mass, [glucoseThreshold, damageThreshold, 
                          trxThreshold])

        # knockout the respective signalling proteins in the Boolean model
        knockout!(knockouts[i], cell.boolean.components)

        # iterate nDelay times to fill up the recent Boolean input states
        simulateLife!(cell, mass, sizeProportion, retention, timestep, nDelay, 
                      regulationFactor, growthFlexibility)

        # simulate life and perturb enzyme i in its usage
        status, rls, genTime, phaseSplit = simulateLife!(cell, mass, 
                                                         sizeProportion, 
                                                         retention, timestep, 
                                                         maxTimesteps, regulationFactor, 
                                                         growthFlexibility)

        # save in output
        idx = i * nCombinations + j
        output[:, idx] = ["∆" * replace(knockouts[i], "+" => "∆"), 
                          rls, genTime, phaseSplit[1][2], 
                          phaseSplit[2][2], phaseSplit[3][2], 
                          formationRate0, repairRate0, status]

    end
end


################################################################################
@printf("Save in output file ...\n")
################################################################################
pretext =
"# SIGNALLING PROTEIN DELETIONS
# FBA model $fbaPath with objective max. growth (parsimonious)
# Boolean model $booleanRulesPath, $booleanSpeciesPath
# TF layer $TFpath
knockout\trls\taverage generation time\trls phase 1\trls phase 2\trls phase 3" *
"\tdamage formation\tdamage repair\tstatus
"
open(outputFile, "w") do file
    write(file, pretext)
    writedlm(file, permutedims(output), "\t")
end 