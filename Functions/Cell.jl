################################################################################
################################################################################
# Cell.jl includes functions for simulating a cell's lifespan with the 
# integrated model (ecFBA for the metabolism, Boolean for the signalling, ODE
# for damage accumulation):
#   - cell = createCell(fbaPath, booleanSpeciesPath, booleanRulesPath, 
#                       TFpath, maxGrowthRate, damageFormation, damageRepair, 
#                       P, D, mass, signallingThresholds)
#   - relevantRefs = getRelevantReferences(fba, boolean)
#   - outputs = simulateLife!(cell, initialMass, timestep, maximalTimesteps, 
#                             regulationFactor, growthFlexibility, 
#                             perturbedEnzymes, perturbationDirection, 
#                             perturbationPhase)
################################################################################
################################################################################


################################################################################
# cell struct
################################################################################
mutable struct Cell 
    fba::FBAModel
    boolean::BooleanModel
    tfTargets::Array{typejoin(Expr, Array{Int, 1}), 2}
    latestBooleanInputs::Array{Array{Bool, 1}, 1}
    ode::ODEModel
    cellcycle::CellCycleModel
    signallingThresholds::Array{Float64, 1}
    relevantRefs::Dict{String, typejoin(VariableRef, ConstraintRef, BooleanComponent)}
end


################################################################################
# initialise a Cell struct including the three models as well as some parameters
# and references to relevant model parts
#
# input parameters:
#   - fbaPath : path to the .mat file of the ecFBA model
#   - booleanSpeciesPath : path to the .txt file specifying the components in the 
#                          Boolean model
#   - booleanRulesPath : path to the .txt file specifying the rules in the 
#                        Boolean model
#   - TFpath : path to the .txt file specifying the gene regulation of enzymes 
#              in the ecFBA model triggered by transcription factors in the 
#              Boolean model
#   - nDelay : indicates how many time steps it takes until the a signal actually
#              affects the metabolic network
#   - maxGrowthRate : maximal growth rate allowed in the ecFBA model
#   - damageFormation : rate at which damage forms, excluding the damage caused 
#                       by the metabolic network
#   - damageRepair : rate at which damage is repaired
#   - P : intact protein fraction
#   - D : damaged protein fraction
#   - mass : cell mass
#   - signallingThresholds : threshold values for the switching on of glucose, 
#                            hydrogen peroxide and Trx in the Boolean model
#
# output parameters:
#   - cell : Cell struct with all relevant information
#
################################################################################
function createCell(fbaPath::String, booleanSpeciesPath::String, 
                    booleanRulesPath::String, TFpath::String, nDelay::Int,
                    maxGrowthRate::Float64, damageFormation::Float64, 
                    damageRepair::Float64, P::Float64, D::Float64, mass::Float64,
                    signallingThresholds::Array{Float64, 1})::Cell

    # initialise ODE model
    ode = initialiseODEModel(maxGrowthRate, damageFormation, damageRepair,
                             P, D, mass)

    # initialise ecFBA model
    fba = initialiseJuMPModel(fbaPath)
    prepareLifespanExperiment!(fba, maxGrowthRate)

    # initialise Boolean model
    boolean = initialiseBooleanModel(booleanSpeciesPath, booleanRulesPath)
    targets = readTranscriptionalTargets(TFpath, fba.enzymeNames)

    # initialize cell cycle model
    cellcycle = initialiseCellCycleModel()

    # extract and save relevant reactions and components that are used very often
    relevantRefs = getRelevantReferences(fba, boolean)

    # use the current Boolean model for the regulation in first few timesteps
    booleanInputs = [relevantRefs["exGlc"].present, 
                     relevantRefs["H2O2"].present,
                     relevantRefs["Trx1_2"].present]
    latestBooleanInputs = [booleanInputs for i = 1 : nDelay + 1]

    # update the GAM component in the FBA according to the maximal growth rate
    updateGAMValue!(relevantRefs["Growth"], 
                    relevantRefs["Maintainance for growth [cytosol]"], 
                    maxGrowthRate)

    # initialise and return cell
    return Cell(fba, boolean, targets, latestBooleanInputs, ode, 
                cellcycle, signallingThresholds, relevantRefs)
end


################################################################################
# map relevant reaction names in an ecFBA and Boolean model to the variables and 
# contraints in the model, in order to facilitate their usage
#
# input parameters:
#   - fba : FBA model struct with all important information
#   - boolean : Boolean model struct with all relevant information
#
# output parameters:
#   - relevantRefs : dictionary that translates each reaction or component name
#                    to the actual variable in the models
# 
################################################################################
function getRelevantReferences(fba::FBAModel, boolean::BooleanModel)::Dict{String, 
                               typejoin(VariableRef, ConstraintRef, BooleanComponent)}

    # extract relevant FBA reactions
    reactionsFBA = ["Growth", "NGAM", "Uptake of glucose", 
                    "Production of ethanol",
                    "Production of ethanol (reversible)",
                    "protein damage exchange (cytosol)",
                    "protein damage exchange (mitochondria)"]
    refsReactions = fba.model[:fluxes][map(x -> findfirst(y -> y == x, 
                              fba.reactionNames), reactionsFBA)]

    # extract relevant FBA enzymes
    enzymesFBA = ["draw_prot_P22217", "draw_prot_P22803"]
    refsEnzymes = fba.model[:enzymes][map(x -> findfirst(y -> y == x, fba.reactionNames), 
                            enzymesFBA) .- length(fba.model[:fluxes])]

    # extract relevant FBA components
    componentsFBA = ["Maintainance for growth [cytosol]"]
    refsComponents = fba.model[:stochiometry][map(x -> findfirst(y -> y == x, 
                               fba.componentNames), componentsFBA)]

    # extract relevant Boolean components for signalling
    componentsBool = ["exGlc", "H2O2", "Trx1_2"]
    refsBool = boolean.components[map(x -> findfirst(y -> y.name == x, 
                                  boolean.components), componentsBool)]

    # create dictionary with the important variable and constraint references 
    allRefs = [refsReactions; refsEnzymes; refsComponents; refsBool]
    allNames = [reactionsFBA; enzymesFBA; componentsFBA; componentsBool]
    relevantRefs = Dict(allNames .=> allRefs)

    return relevantRefs
end 


################################################################################
# simulate the life of a cell with the integrated model (ecFBA for the 
# metabolism, Boolean for the signalling, ODE for damage accumulation) until 
# the maximal number of timesteps are reached or the cell died
# 
# a cells life is divided into three phases:
#   1. maximal growth phase (via fermentation)
#   2. switch to respiration due to increased ATP demand for repair
#   3. ethanol respiration 
#
# optional is the perturbation of a specific enzyme to perform LCC calculations, 
# which can be a deletion or an overexpression of an enzyme during either 
# the whole lifespan or only a specific phase
#
# input parameters:
#   - cell : Cell struct with all relevant information
#   - minimalMass : minimal mass of the mother cell (important for division)
#   - sizeProportion : fraction of biomass that remains in the mother cell at 
#                      cell divison
#   - retention : fraction of damage in the mother that is caused by 
#                 damage retention
#   - timestep : time of the iteration
#   - maximalTimesteps : maximal number of timesteps before the simulation stops
#   - regulationFactor : factor determining how much the Boolean signalling layer
#                        affects the enzyme usages in the ecFBA model
#   - growthFlexibility : defines how much the growth rate can deviate from the 
#                         optimum
#   - perturbedEnzymes : enzyme that is perturbed in its usage 
#   - perturbationDirection : "deletion" or "overexpression"
#   - perturbationPhase : "complete", "phase 1", "phase 2" or "phase 3"
#
# output parameters:
#   - status : indicates if the cell is still alive, or the reason of death
#   - rls : replicative lifespan of the cell
#   - averageGenerationTime : average time between divisions
#   - phaseSplit : times in, number of divisions during and damage levels at the 
#                  end of the three phases described above
#   - perturbedEnzymeUsage : total amount of the perturbed enzyme(s) that was used,
#                            needed for LCC calculations, in particular in the 
#                            case of overexpression
# 
################################################################################
function simulateLife!(cell::Cell, minimalMass::Float64, sizeProportion::Float64, 
                       retention::Float64, timestep::Float64, 
                       maximalTimesteps::Int, regulationFactor::Float64, 
                       growthFlexibility::Float64, perturbedEnzymes::typejoin(VariableRef, 
                       Array{VariableRef, 1}) = [], 
                       perturbationDirection::String = "",
                       perturbationPhase::String = ""
                       )::Tuple{String, Int64, Float64, 
                       Array{Array{Float64, 1}, 1}, Array{Float64, 1}}

    
    ############################################################################
    # INITIALISATION
    ############################################################################
    # initialise all counts
    status = "alive"
    time = 0.0
    rls = 0
    divisionTimes = [0.0]
    it = 1
    phaseSplit = [[0.0, 0.0, NaN] for i = 1 : 3]
    currentPhase = "phase 1"
    phaseShift = [false, false]
    perturbedEnzymeUsage = zeros(Float64, length(perturbedEnzymes))

    # extract relevant model parts that are used over and over again
    fba = cell.fba
    model = fba.model
    boolean = cell.boolean
    ode = cell.ode

    # extract important references to the models
    refs = cell.relevantRefs
    glucoseUptakeFlux = refs["Uptake of glucose"]
    ethanolProductionFlux = refs["Production of ethanol"]
    ethanolUptakeFlux = refs["Production of ethanol (reversible)"]
    growthFlux = refs["Growth"]
    GAMComponent = refs["Maintainance for growth [cytosol]"]
    NGAM = refs["NGAM"]
    damageFluxes = [refs["protein damage exchange (cytosol)"],
                    refs["protein damage exchange (mitochondria)"]]
    trxFluxes = [refs["draw_prot_P22217"],
                 refs["draw_prot_P22803"]]
    booleanGlucose = refs["exGlc"]
    booleanPeroxoide = refs["H2O2"]
    booleanTrx = refs["Trx1_2"]

    # save initial enzyme, growth and glucose bounds
    growthLowerBound = lower_bound(growthFlux)
    glucoseUpperBound = upper_bound(glucoseUptakeFlux)
    enzymeLowerBounds = lower_bound.(model[:enzymes])
    enzymeUpperBounds = upper_bound.(model[:enzymes])

    # get maximal growth rate and NGAM 
    maxGrowthRate = upper_bound(growthFlux)
    initialGrowth = maxGrowthRate
    maxNGAM = upper_bound(NGAM)


    ############################################################################
    # SIMULATE TIMESTEP BY TIMESTEP
    ############################################################################
    for n = 1 : maximalTimesteps
    
        ########################################################################
        # UPDATE ecFBA PARAMETERS

        # restrict the usage of the deleted enzyme to 0 (if wanted)
        if !isempty(perturbedEnzymes) && perturbationDirection == "deletion" &&
            (perturbationPhase == "complete" || perturbationPhase == currentPhase)
            set_upper_bound.(perturbedEnzymes, 0.0)
        end
    
        # limit the enzyme pool according to current ODE state
        P = ode.state[1]
        set_upper_bound(model[:pool], P * 0.1799 * 0.4592)

        # set NGAM linearly dependent on the current amount of damage in the cell
        # between 0 and a maximal value (set as upper bound, from Lu et al., 2019)
        D = ode.state[2]
        set_lower_bound(NGAM, D * maxNGAM / (P + D))

        ########################################################################
        # SOLVE ecFBA MODEL FOR MAXIMAL GROWTH
        pFBA!(model, "Max", growthFlux)

        # stop if the FBA cannot solve at all anymore (no possible growth rate)
        if termination_status(model) != MOI.OPTIMAL 
            status = "dead (regular)"
            break
        end

        # extract and fix growth
        growth = value(growthFlux)
        set_lower_bound(growthFlux, growth * (1 - growthFlexibility))

        # update the GAM to this growth rate
        updateGAMValue!(growthFlux, GAMComponent, growth)

        ########################################################################
        # FIND THE LOGICAL STEADY STATE OF THE BOOLEAN MODEL AND REGULATE THE
        # ecFBA ACCORDINGLY

        # update the Boolean model input with the ones from nDelay timesteps
        # ago, i.e. the first element in the latestBooleanInputs array
        tmpBooleanInputs = cell.latestBooleanInputs[1]
        refs["exGlc"].present = tmpBooleanInputs[1]
        refs["H2O2"].present = tmpBooleanInputs[2]
        refs["Trx1_2"].present = tmpBooleanInputs[3]

        # run Boolean model 
        runBooleanModel!(boolean)

        # calculate ranks for gene expression given the transcription factor 
        # activity
        ranks = getExpressionRank(boolean.components, fba.enzymeNames, 
                                  cell.tfTargets)

        # regulate FBA according to the ranks from the Boolean model
        regulateFBA!(model, ranks, regulationFactor)

        ########################################################################
        # SOLVE THE REGULATED ecFBA MODEL
        pFBA!(model, "Max", growthFlux)

        # stop if the FBA cannot solve anymore
        if termination_status(model) != MOI.OPTIMAL 
            status = "dead (regulation)"
            break
        end

        # in case of an enzyme overexpression
        if !isempty(perturbedEnzymes) && perturbationDirection == "overexpression" &&
            (perturbationPhase == "complete" || perturbationPhase == currentPhase)

            # perturb the respective usage by some percentage
            perturbation = 0.5
            set_lower_bound.(perturbedEnzymes, value.(perturbedEnzymes) * 
                                               (1 + perturbation - precision))
            set_upper_bound.(perturbedEnzymes, value.(perturbedEnzymes) *  
                                               (1 + perturbation + precision))
            
            # increase pool to compensate for the overexpression
            proteinPool = findfirst(x -> x == "prot_pool", fba.componentNames)
            molecularWeights = normalized_coefficient.(model[:stochiometry][proteinPool],
                                                       perturbedEnzymes) .* (-1)
            set_upper_bound.(model[:pool], upper_bound(model[:pool]) + 
                             (molecularWeights' * value.(perturbedEnzymes) .* 
                             (perturbation + precision))) 

            # and solve once more
            pFBA!(model, "Max", growthFlux)

            if termination_status(model) != MOI.OPTIMAL 

                # if infeasible, allow to go down in the growth rate, to avoid 
                # that the model becomes infeasible only because of the overexpression
                set_lower_bound(growthFlux, 0.0)

                # and solve once more
                pFBA!(model, "Max", growthFlux)

                # break if it still cannot solve
                if termination_status(model) != MOI.OPTIMAL 
                    status = "dead (overexpression)"
                    break
                end

            end
        end

        ########################################################################
        # ADD ENZYME USAGES TO TOTAL USAGE OF THE PERTURBES ENZYMES
        perturbedEnzymeUsage .+= value.(perturbedEnzymes)

        ########################################################################
        # SOLVE ODE MODEL WITH PARAMETERS TAKEN FROM THE OPTIMAL ecFBA SOLUTION

        # extract the growth and damage fluxes of the optimal solution
        damage = sum(value.(damageFluxes))
        repair = 0.0
        growth = value(growthFlux)

        # set variable rate parameters for the ODE model accordingly
        ode.variableParameters = [damage, repair, growth]

        # update ODE for one timestep
        solveODE!(ode, timestep)
    
        ########################################################################
        # CELL DIVISION
        # Update cell cycle model parameters based on current state
        # solve around timestep
        cell.cellcycle.tspan=(0.9*time,time+1.5*time)
        updateCellCycleParameters!(cell.cellcycle, cell.fba, cell.boolean)
        # Solve cell cycle ODE for one timestep
        cp1 = getCCParamValue(cell.cellcycle, "compartment_1")
        sol = solve(cell.cellcycle.cell_cycle_prob, Tsit5())
        cell.cellcycle.u0 = cp1 .* sol(time) 

        # Check division based on cell cycle state
        if shouldDivide(cell.cellcycle)
            #println("SPLIT")
            division!(ode, sizeProportion, retention)
            rls += 1
            push!(divisionTimes, time)
            # Reset cycle after division
            cell.cellcycle = initialiseCellCycleModel()
        end

        ########################################################################
        # UPDATE THE BOOLEAN MODEL ACCORDING TO THE OPTIMAL ecFBA SOLUTION
        # update cell's signalling according to solution
        triggerSignalling!([glucoseUptakeFlux], cell.signallingThresholds[1], 
                booleanGlucose)
        triggerSignalling!(damageFluxes, cell.signallingThresholds[2], 
                        booleanPeroxoide)
        triggerSignalling!(trxFluxes, cell.signallingThresholds[3], booleanTrx)

        # update current Boolean model to be used in nDelay time steps,
        # i.e. the last element in the latestBooleanInputs array
        cell.latestBooleanInputs[end] = [refs["exGlc"].present, 
                                         refs["H2O2"].present,
                                         refs["Trx1_2"].present]

        # and shift all elements by one for the next time step
        cell.latestBooleanInputs[1:end-1] = cell.latestBooleanInputs[2:end]

        ########################################################################
        # DECIDE WHICH PHASE THE CELL IS IN
        # distinguish between phases (respect precision of solver)
        it == 1 ? initialGrowth = growth : nothing
        highestGrowth = growth - initialGrowth * 0.95 >= -precision
        ethanolUptake = (value(ethanolUptakeFlux) - value(ethanolProductionFlux)) > 
                         precision

        # phase 1
        # update time, rls and damage if the growth is still above the maximal 
        # values with some flexibility (and phase 2 has not been reached before)
        if highestGrowth && !ethanolUptake && phaseShift[1] == false
            currentPhase = "phase 1"
            phaseSplit[1][1] += timestep
            phaseSplit[1][2] = rls
            phaseSplit[1][3] = cell.ode.state[2]
        # phase 2
        # update time, rls and damage if the cell is switching to respiration
        # (and phase 3 has not been reached before)
        elseif !highestGrowth && !ethanolUptake && phaseShift[2] == false
            currentPhase = "phase 2"
            phaseShift[1] = true
            phaseSplit[2][1] += timestep
            phaseSplit[2][2] = rls - phaseSplit[1][2] 
            phaseSplit[2][3] = cell.ode.state[2]
        # phase 3
        # update time, rls and damage if there is ethanol respiration
        elseif ethanolUptake
            currentPhase = "phase 3"
            phaseShift[2] = true
            phaseSplit[3][1] += timestep
            phaseSplit[3][2] = rls - phaseSplit[1][2] - phaseSplit[2][2] 
            phaseSplit[3][3] = cell.ode.state[2]
        end

        ########################################################################
        # UPDATE ecFBA PARAMETERS FOR NEXT TIMESTEP
        # unconstrain enzymes, growth and glucose again
        # (those are the ones that are changed in each simulation step)
        set_lower_bound(growthFlux, growthLowerBound)
        set_upper_bound(glucoseUptakeFlux, glucoseUpperBound)
        set_lower_bound.(model[:enzymes], enzymeLowerBounds)
        set_upper_bound.(model[:enzymes], enzymeUpperBounds)

        # update time and iteration
        time += timestep
        it += 1
        
    end 
    
    ########################################################################
    # GET AVERAGE GENERATION TIME AND RETURN RESULTS
    if rls > 0
        averageGenerationTime = mean(divisionTimes[2:end-1] .- 
                                     divisionTimes[1:end-2])
    else
        averageGenerationTime = 0.0
    end

    ########################################################################
    # RETURN
    output = (status, rls, averageGenerationTime, phaseSplit,
              perturbedEnzymeUsage)
    return output

end

function updateCellCycleParameters!(cellcycle::CellCycleModel, fba::FBAModel, boolean::BooleanModel)
    ########################################################################
    # Update cell cycle parameters based on FBA and Boolean states
    #println(getCCParamValue(cellcycle, "kp_Cln3"))
    #println(getCCParamValue(cellcycle, "nutrition_factor"))
 
    # Growth-dependent parameters
    growth_rate = value(fba.model[:fluxes][findfirst(x -> x == "Growth", fba.reactionNames)])
    # Update parameters based on FBA and Boolean states
    glucose_flux = value(fba.model[:fluxes][findfirst(x -> x == "Uptake of glucose", fba.reactionNames)])
    setCCParamValue(cellcycle, "kp_Cln3", getCCParamOrigValue(cellcycle, "kp_Cln3") * (1 + glucose_flux * growth_rate))
    
    # What should this be?
    reference_mass = 1 
    
    protein_mass = value(fba.model[:pool])
    setCCParamValue(cellcycle, "nutrition_factor", protein_mass / reference_mass)

    # Stress response affects multiple parameters
    #stress_active = boolean.components[findfirst(x -> x.name == "H202", boolean.components)].present
    #stress_factor=0.9
    #if stress_active
    #    setCCParamValue(cellcycle, "kp_Far1", getCCParamOrigValue(cellcycle, "kp_Far1") * stress_factor )
    #    setCCParamValue(cellcycle, "kp_Cln2", getCCParamOrigValue(cellcycle, "kp_Cln2") / stress_factor )
    #else
    #    setCCParamValue(cellcycle, "kp_Far1", getCCParamOrigValue(cellcycle, "kp_Far1"))
    #    setCCParamValue(cellcycle, "kp_Cln2", getCCParamOrigValue(cellcycle, "kp_Cln2"))
    #end
     
    # Update the cell cycle ode problem
    updateCCProblem(cellcycle)
end

function shouldDivide(cellcycle::CellCycleModel)
    ########################################################################
    # Check if certain states are above / below threshold values
    return getCCStateValue(cellcycle,"Clb2") < 0.24 && getCCStateValue(cellcycle,"Cdc14") > 3400 && getCCStateValue(cellcycle,"Sic1") > 300
end
