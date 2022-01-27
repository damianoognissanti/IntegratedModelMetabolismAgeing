################################################################################
################################################################################
# simulate the life of a cell with the integrated model 
# (ecFBA for the metabolism, Boolean for the signalling, ODE
# for damage accumulation)
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
@printf("Initialise the output ... \n")
################################################################################
outputFile = "life.txt"
output = zeros(Float64, 20, maxTimesteps)


################################################################################
@printf("Initialise cell ... \n")
################################################################################
cell = createCell(fbaPath, booleanSpeciesPath, booleanRulesPath, TFpath, nDelay,
                  maxGrowthRate, formationRate0, repairRate0, P, D, mass,
                  [glucoseThreshold, damageThreshold, trxThreshold])


################################################################################
@printf("Get exchange reactions of interest... \n")
################################################################################
# extract few more reactions that are outputed over time
exchangeReactions = ["Uptake of O2", "Production of CO2", "Production of ethanol",
                     "Production of ethanol (reversible)", "Production of acetate", 
                     "Production of glycerol", "Production of pyruvate"]
exchangeFluxes = cell.fba.model[:fluxes][map(x -> findfirst(y -> y == x, 
                                         cell.fba.reactionNames), 
                                         exchangeReactions)]
damageReactions = ["superoxide production (mitochondria)", 
                   "superoxide oxidoreductase (arm, mitochondria)"]
damageFluxes = cell.fba.model[:fluxes][map(x -> findfirst(y -> y == x, 
                                       cell.fba.reactionNames), 
                                       damageReactions)]
otherReactions = ["ATP hydrolysis"]
otherFluxes = cell.fba.model[:fluxes][map(x -> findfirst(y -> y == x, 
                                      cell.fba.reactionNames), 
                                      otherReactions)]
enzymePool = cell.fba.model[:pool]


################################################################################
@printf("Iterate over time ... \n")
################################################################################
# iterate nDelay times to fill up the recent Boolean input states
simulateLife!(cell, mass, sizeProportion, retention, timestep, nDelay, 
               regulationFactor, growthFlexibility)

# initialise time and division times
time = 0.0
divisionTimes = [0.0]
rls = 0
it = 1

# go through time steps individually and save output fluxes over time
# NOTE : since we go only one time step at a time the simulateLife!() 
# function always start from 0 in each time step, while we count time and rls
# here instead
for n = 1 : maxTimesteps

    global cell, time, rls, divisionTimes, it

    # save the currently used inputs for the Boolean model
    tmpBooleanInputs = cell.latestBooleanInputs[1]

    # run one time step (maximalTimesteps = 1)
    status, tmpRls = simulateLife!(cell, mass, sizeProportion, retention,
                                   timestep, 1, regulationFactor, 
                                   growthFlexibility)

    # break if there is no solution anymore
    !(status == "alive") ? break : nothing

    # save variables of interest in output
    output[:, it] = [time; value(cell.relevantRefs["Growth"]); 
                    value(cell.relevantRefs["Uptake of glucose"]);
                    value.(exchangeFluxes[1:2]); 
                    value(exchangeFluxes[3]) - value(exchangeFluxes[4]);
                    value.(exchangeFluxes[5:end]);
                    value.(otherFluxes)
                    cell.ode.fixedParameters[1] +
                    cell.ode.variableParameters[1];
                    cell.ode.state; 
                    value(enzymePool); 
                    tmpBooleanInputs;
                    value.(damageFluxes)]

    # update time, rls and division times
    if tmpRls == 1
        rls += 1
        push!(divisionTimes, time)
    end
    time += timestep
    it += 1

end 

# get average generation time
if rls > 0
    averageGenerationTime = mean(divisionTimes[2:end-1] .- divisionTimes[1:end-2])
else
    averageGenerationTime = 0.0
end


################################################################################
@printf("Save in output file ...\n")
################################################################################
pretext = 
"# LIFESPAN SIMULATION WITH INTEGRATED MODEL
# FBA model $fbaPath with objective max. growth (parsimonious)
# Boolean model $booleanRulesPath, $booleanSpeciesPath
# TF layer $TFpath
# non-metabolic damage formation f0 = $formationRate0
# damage repair r0 = $repairRate0
# $rls division(s) with average generation time $averageGenerationTime
time\tgrowth\tglucose\tO2\tCO2\tethanol\tacetate\tglycerol\t" *
"pyruvate\tATP\tDformation\tP\tD\tM\tenzymeUsage\tglucose(bool)\t" * 
"H2O2(bool)\tTrx1/2(bool)\tsuperoxide\tH2O2 production
"

open(outputFile, "w") do file
    write(file, pretext)
    writedlm(file, transpose(output[:, 1:it-1]), "\t")
end 