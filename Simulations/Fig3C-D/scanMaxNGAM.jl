################################################################################
################################################################################
# simulate the lifespan of a cell with the integrated model 
# (ecFBA for the metabolism, Boolean for the signalling, ODE
# for damage accumulation) for different maximal NGAM values (repair costs)
# for parameter combinations that lead to wildtypes (23 divisions)
################################################################################
################################################################################

using Printf
using DelimitedFiles
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
outputFile = "ngamEffects.txt"
nValues = length(ngamValues)
output = Array{Any, 2}(undef, 13, nValues * 3 * nCombinations)


################################################################################
@printf("Simulate cells for different maximal NGAM values ... \n")
################################################################################
for i = 1 : nValues
    for n = 1 : nCombinations

        @printf("max NGAM : %.2f, parameter set %i/%i\n", ngamValues[i], n, 
                nCombinations)

        # extract parameters
        formationRate0 = parameters[n, 1]
        repairRate0 = parameters[n, 2]

        # initialise cell 
        cell = createCell(fbaPath, booleanSpeciesPath, booleanRulesPath, TFpath,
                          nDelay, maxGrowthRate, formationRate0, repairRate0,
                          P, D, mass, [glucoseThreshold, damageThreshold, 
                          trxThreshold])

        # set max. NGAM value
        set_upper_bound(cell.relevantRefs["NGAM"], ngamValues[i])

        # iterate nDelay times to fill up the recent Boolean input states
        status = simulateLife!(cell, mass, sizeProportion, retention, timestep,  
                               nDelay, regulationFactor, growthFlexibility)

        # simulate life
        # initialise time and division times
        time = 0.0
        divisionTimes = [0.0]
        rls = 0
        it = 1

        # prepare finding the phase shifts
        initialGrowth = maxGrowthRate
        phaseShift = [false, false]
        shiftIdx = [0, 0]
        phaseSplit = [[0.0, 0.0, NaN] for i = 1 : 3]

        # get relevant fluxes
        modelFluxes = cell.fba.model[:fluxes]
        growthFlux = cell.relevantRefs["Growth"]
        glucoseUptakeFlux = cell.relevantRefs["Uptake of glucose"]
        ethanolUptakeFlux = cell.relevantRefs["Production of ethanol (reversible)"]
        ethanolProductionFlux = cell.relevantRefs["Production of ethanol"]
        O2UptakeFlux = modelFluxes[findfirst(x -> x == "Uptake of O2", 
                                             cell.fba.reactionNames)]
        CO2ProductionFlux = modelFluxes[findfirst(x -> x == "Production of CO2", 
                                                  cell.fba.reactionNames)]
        acetateProductionFlux = modelFluxes[findfirst(x -> x == "Production of acetate", 
                                                      cell.fba.reactionNames)]
        exchangeFluxes = [glucoseUptakeFlux, ethanolProductionFlux, ethanolUptakeFlux, 
                          O2UptakeFlux, CO2ProductionFlux, acetateProductionFlux]
        outputFluxes = zeros(Float64, 6, maxTimesteps)

        # go through time steps individually and save output fluxes over time
        # NOTE : since we go only one time step at a time the simulateLife!() 
        # function always start from 0 in each time step, while we count time and rls
        # and the phases here instead
        for n = 1 : maxTimesteps

            # save the currently used inputs for the Boolean model
            tmpBooleanInputs = cell.latestBooleanInputs[1]

            # run one time step (maximalTimesteps = 1)
            status, tmpRls = simulateLife!(cell, mass, sizeProportion, retention, 
                                           timestep, 1, regulationFactor, 
                                           growthFlexibility)

            # break if there is no solution anymore
            !(status == "alive") ? break : nothing

            # save fluxes (normalised to glucose intake)
            outputFluxes[:, it] = value.(exchangeFluxes)
            outputFluxes[2:end, it] ./= outputFluxes[1, it] 

            # update rls and division times
            if tmpRls == 1
                rls += 1
                push!(divisionTimes, time)
            end

            # distinguish between phases (respect precision of solver)
            it == 1 ? initialGrowth = value(growthFlux) : nothing
            highestGrowth = value(growthFlux) - initialGrowth * 0.95 >= -precision
            ethanolUptake = (value(ethanolUptakeFlux) - value(ethanolProductionFlux)) > 
                            precision
            if highestGrowth && !ethanolUptake && phaseShift[1] == false
                phaseSplit[1][1] += timestep
                phaseSplit[1][2] = rls
                phaseSplit[1][3] = cell.ode.state[2]
            elseif !highestGrowth && !ethanolUptake && phaseShift[2] == false
                if  phaseShift[1] == false
                    phaseShift[1] = true
                    shiftIdx[1] = it
                end
                phaseSplit[2][1] += timestep
                phaseSplit[2][2] = rls - phaseSplit[1][2] 
                phaseSplit[2][3] = cell.ode.state[2]
            elseif ethanolUptake
                if phaseShift[2] == false
                    phaseShift[2] = true
                    shiftIdx[2] = it
                end
                phaseSplit[3][1] += timestep
                phaseSplit[3][2] = rls - phaseSplit[1][2] - phaseSplit[2][2] 
                phaseSplit[3][3] = cell.ode.state[2]
            end

            # update time
            time += timestep
            it += 1

        end 
        
        # get ethanol exchange as difference between production and uptake
        outputFluxes[2, :] .-= outputFluxes[3, :]
        outputFluxes = outputFluxes[[1, 2, 4, 5, 6], :]

        # get exchange fluxes at end of phases
        # (in some cases there is no phase 3)
        phase1Exchange = outputFluxes[:, shiftIdx[1] - 1]
        if phaseShift[2] == true
            phase2Exchange = outputFluxes[:, shiftIdx[2] - 1]
            phase3Exchange = outputFluxes[:, it - 1]
        else
            phase2Exchange = outputFluxes[:, it - 1]
            phase3Exchange = [NaN for k = 1:5]
        end

        # save in output
        idx = (i - 1) * 3 * nCombinations + (n - 1) * 3
        output[:, idx + 1] = [ngamValues[i]; "phase 1"; phaseSplit[1]; 
                              phase1Exchange; formationRate0; repairRate0; status]
        output[:, idx + 2] = [ngamValues[i]; "phase 2"; phaseSplit[2]; 
                              phase2Exchange; formationRate0; repairRate0; status]
        output[:, idx + 3] = [ngamValues[i]; "phase 3"; phaseSplit[3]; 
                              phase3Exchange; formationRate0; repairRate0; status]
                              
    end
end


################################################################################
@printf("Save in output file ...\n")
################################################################################
pretext =
"# NGAM GRID
# FBA model $fbaPath with objective max. growth (parsimonious)
# Boolean model $booleanRulesPath, $booleanSpeciesPath
# TF layer $TFpath
NGAM\tphase\ttime\tdivisions\tdamage\tglucose uptake (end)\t" *
"norm. ethanol production (end)\tnorm. O2 uptake (end)\tnorm. CO2 production (end)\t" *
"norm. acetate production (end)\tdamage formation\tdamage repair\tstatus
"

open(outputFile, "w") do file
    write(file, pretext)
    writedlm(file, permutedims(output), "\t")
end 