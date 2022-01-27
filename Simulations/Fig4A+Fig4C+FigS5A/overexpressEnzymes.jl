################################################################################
################################################################################
# simulate the lifespan of a cell with the integrated model 
# (ecFBA for the metabolism, Boolean for the signalling, ODE
# for damage accumulation) for positive perturbations (overexpressions) 
# of all enzymes (individually)
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
@printf("Get wildtype rls... \n")
################################################################################
# initialise cell
wildtypeCell = createCell(fbaPath, booleanSpeciesPath, booleanRulesPath, TFpath,
                          nDelay, maxGrowthRate, formationRate0, repairRate0,
                          P, D, mass, [glucoseThreshold, damageThreshold, 
                          trxThreshold])
enzymeNames = wildtypeCell.fba.enzymeNames

# iterate nDelay times to fill up the recent Boolean input states
simulateLife!(wildtypeCell, mass, sizeProportion, retention, timestep, nDelay,
              regulationFactor, growthFlexibility)

# simulate life
status, wtRls, wtGenTime, wtPhases, wtUsage = simulateLife!(wildtypeCell, mass, 
                                                sizeProportion, retention,
                                                timestep, maxTimesteps, 
                                                regulationFactor, growthFlexibility,
                                                wildtypeCell.fba.model[:enzymes])


################################################################################
@printf("Create list of all enzymes plus all isoenzyme combinations ...\n")
################################################################################
nEnzymes = size(enzymeNames, 1)
enzymeCombinations = [[i] for i = 1 : size(enzymeNames, 1)]
isoenzymes = getIsoenzymes(wildtypeCell.fba)
enzymeCombinations = [enzymeCombinations; isoenzymes]
nOverexpressions = length(enzymeCombinations)


################################################################################
@printf("Map all enzymes to the pathways where they are included ...\n")
################################################################################
enzymesToComponents, componentsToFluxes = mapEnzymeIdx(wildtypeCell.fba)
idxArray = collect(1 : nEnzymes)
pathways = ["" for i = 1 : nEnzymes]
for n = 1 : nEnzymes
    comp = enzymesToComponents[idxArray[n]]
    flux = componentsToFluxes[comp]
    pathway = wildtypeCell.fba.reactionPathways[flux]
    pathways[n] = join(unique(pathway), ", ")
end

# add also the isoenzyme combinations
combinationPathways = ["" for i = 1 : length(isoenzymes)]
for n = 1 : length(isoenzymes)
    combinationPathways[n] = join(unique(pathways[isoenzymes[n]]), ", ")
end

# combine both parts
pathways = [pathways; combinationPathways]


################################################################################
@printf("Initialise the output ... \n")
################################################################################
outputFile = "enzymeOverexpressions.txt"
output = Array{Any, 2}(undef, 18, 4 * nOverexpressions + 1)

output[:, 1] = ["wildtype"; "wildtype"; ""; ""; wtRls; 0.0; 0.0; 0.0; 
                wtPhases[1]; wtPhases[2]; wtPhases[3]; status]


################################################################################
@printf("Perturb each enzyme ... \n")
################################################################################
for i = 1 : nOverexpressions

    global output
    local status

    idx = enzymeCombinations[i]
    @printf("\t%s %i/%i\n", join(enzymeNames[idx, 2], ","), i, nOverexpressions)

    overexpressionSpan = ["complete", "phase 1", "phase 2", "phase 3"]

    for j = 1 : length(overexpressionSpan)

        # initialise cell
        cell = createCell(fbaPath, booleanSpeciesPath, booleanRulesPath, TFpath,
                          nDelay, maxGrowthRate, formationRate0, repairRate0,
                          P, D, mass, [glucoseThreshold, damageThreshold, 
                          trxThreshold])

        # iterate nDelay times to fill up the recent Boolean input states
        simulateLife!(cell, mass, sizeProportion, retention, timestep, nDelay, 
                      regulationFactor, growthFlexibility)

        # simulate life and perturb enzyme i in its usage
        status, rls, genTime, phases, usage = simulateLife!(cell, mass, 
                                                    sizeProportion, retention,
                                                    timestep, maxTimesteps, 
                                                    regulationFactor, 
                                                    growthFlexibility, 
                                                    cell.fba.model[:enzymes][idx], 
                                                    "overexpression", 
                                                    overexpressionSpan[j])

        # calculate the LCC
        rlsChange = rls / wtRls
        enzymeChange = sum(usage) / sum(wtUsage[idx])
        lcc = (rls - wtRls) / wtRls * 
              abs(sum(wtUsage[idx]) / (sum(usage) - sum(wtUsage[idx])))

        # add to output
        enzyme = join(enzymeNames[idx, 2], ",")
        enzymeStandard = join(enzymeNames[idx, 3], ",")
        output[:, 1 + (i-1)*4 + j] = [enzyme; enzymeStandard; pathways[i];
                                      overexpressionSpan[j]; rls; rlsChange; enzymeChange;
                                      lcc; phases[1]; phases[2]; phases[3]; status]
    end
end

# sort by rls
output = sortslices(output, dims = 2, by = x -> x[4])


################################################################################
@printf("Save in output file ...\n")
################################################################################
pretext = 
"# OVEREXPRESSION SIMULATIONS
# FBA model $fbaPath with objective max. growth (parsimonious)
# Boolean model $booleanRulesPath, $booleanSpeciesPath
# TF layer $TFpath
# non-metabolic damage formation f0 = $formationRate0
# damage repair r0 = $repairRate0
enzyme(s)\tstandard name\tpathway\toverexpression span\trls\trls change\tenzyme change\t" *
"lcc\ttime phase 1\trls phase 1\tdamage phase 1\ttime phase 2\trls phase 2\t" *
"damage phase 2\ttime phase 3\trls phase 3\tdamage phase 3\tstatus
"

open(outputFile, "w") do file
    write(file, pretext)
    writedlm(file, permutedims(output), "\t")
end 
