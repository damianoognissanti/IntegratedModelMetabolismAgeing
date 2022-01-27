################################################################################
################################################################################
# FBA.jl includes useful functions for performing a (parsimonious) flux balance
# analysis:
#   - updateGAMValue!(growthReaction, GAMComponent, growthRate)
#   - pFBA!(model, extreme, objectiveFlux)
#   - minFlux, maxFlux = calculateFVA(model, fluxes)
#   - minUsage, maxUsage = calculateEVA(model, enzymes)
#   - enzymesToComponents, componentsToFluxes = mapEnzymeIdx(fba)
#   - eccMatrix, kcatChange = calculateECC(fba, controlledFluxes, perturbation, 
#                             extreme, objectiveFlux)
################################################################################
################################################################################


################################################################################
# calculate and update GAM depending on the growth rate 
# according to Österberg et al., PLOS Computational Biology 2021 
#
# input parameters:
#   - growthReaction : reference to the variable in the JuMPmodel that 
#                      corresponds to the growth
#   - GAMComponent : reference to the mass balance constraint in the JuMPmodel  
#                    that corresponds to the GAM
#   - growthRate : current growth rate (obs: formula only valid for growth rates
#                  between 0.0 and 0.4)
#
################################################################################
function updateGAMValue!(growthReaction::VariableRef, GAMComponent::ConstraintRef, 
                         growthRate::Float64)::Nothing

    # linear change from 18 to 30 mmol / gDW h
    # for growth rates between 0.0 and 0.285
    if 0.0 <= growthRate <= 0.285
        GAM = 18 + growthRate * (30 - 18) / 0.285
        
        # update the reaction stochiometry accordingly
        set_normalized_coefficient(GAMComponent, growthReaction, - GAM)

    # linear change from 30 to 25 mmol / gDW h
    # for growth rates between 0.285 and 0.4
    elseif 0.285 < growthRate <= 0.4
        GAM = 30 + (growthRate - 0.285) * (25 - 30) / (0.4 - 0.285)

        # update the reaction stochiometry accordingly
        set_normalized_coefficient(GAMComponent, growthReaction, - GAM)
    end

    # return
    return nothing
end


################################################################################
# solve parsimonious FBA for a JuMP model:
# two successive optimisations
#   (1) optimise a given flux (from input) 
#   (2) given the optimal solution of (1), minimze the enzyme usage and the
#       sum of all fluxes 
#
# input parameters:
#   - model : JuMP optimisation model of the metabolic network
#   - extreme : "Min" or "Max" for minimisation or maximisation problem 
#               respectively
#   - objectiveFlux : reference to the variable in the JuMP model that should 
#                     be optimised
#
################################################################################
function pFBA!(model::Model, extreme::String, 
               objectiveFlux::VariableRef)::Nothing

    # set objective for first optimisation
    if extreme == "Max"
        @objective(model, Max, objectiveFlux)
    elseif extreme == "Min"
        @objective(model, Min, objectiveFlux)
    end

    # perform first optimisation of parsimonious FBA
    optimize!(model)
    
    # only if the first optimisation was successful, do the second
    if termination_status(model) == MOI.OPTIMAL
    
        # get the optimal value to set as a constraint in the second
        # optimisation
        optimalValue = value(objectiveFlux)
        
        # save bounds for setting model back to input state
        tmpLowerBound = lower_bound(objectiveFlux)
        tmpUpperBound = upper_bound(objectiveFlux)
        
        # keep optimal glucose uptake for the second optimisation
        # we assume that in an irreversible model all fluxes are positive
        set_lower_bound(objectiveFlux, optimalValue * (1 - precision))
        set_upper_bound(objectiveFlux, optimalValue * (1 + precision))
        
        # set objective value for second optimisation
        parsimoniousObjective = model[:pool] + sum(model[:fluxes])
        @objective(model, Min, parsimoniousObjective)
        
        # perform second optimisation of parsimonious FBA
        optimize!(model)

        #= error message if not feasible
        if termination_status(model) != MOI.OPTIMAL
            @printf("ERROR: parsimonious FBA with status %s.\n",check 
                    termination_status(model))
        end =#
        
        # reset bounds
        set_lower_bound(objectiveFlux, tmpLowerBound)
        set_upper_bound(objectiveFlux, tmpUpperBound)

    #= error message if not feasible
    else
        @printf("ERROR: FBA with status %s.\n", termination_status(model)) =#
    end

    # return
    return nothing
end


################################################################################
# perform an flux variablility analysis for the FBA model by minimising and 
# maximising given fluxes of interest to see how sensitive a previously constraint 
# system is to changes the respective fluxes
# (for further details see Mahadevan et al., Metabolic Engeneering 2003)
#
# input parameters:
#   - model : JuMP optimisation model of the metabolic network
#   - fluxes : references to the fluxes that should be investigated with FVA
#
# output parameters:
#   - minFlux : minimal fluxes under the given conditions
#   - maxFlux : maximal fluxes under the given conditions
#
################################################################################
function calculateFVA(model::Model, fluxes::Array{VariableRef, 1} 
                      = Array{VariableRef, 1}(undef, 0))::Tuple{Array{Float64, 1}, 
                      Array{Float64, 1}}

    # if no indices are given check all fluxes
    isempty(fluxes) ? fluxes = model[:fluxes] : nothing

    # initialise output
    nFluxes = length(fluxes)
    minFlux = [NaN for i = 1 : nFluxes]
    maxFlux = [NaN for i = 1 : nFluxes]

    # go through all enzymes and minimise and maximise the respective usage given
    # the input model (that should potentially have constraints on other fluxes)
    for i = 1 : nFluxes

        # solve for minimal usage
        @objective(model, Min, fluxes[i])
        optimize!(model)

        # save if optimisation was successful
        if termination_status(model) == MOI.OPTIMAL
            minFlux[i] = value(fluxes[i])
        end

        # solve for maximal usage
        @objective(model, Max, fluxes[i])
        optimize!(model)

        # save if optimisation was successful
        if termination_status(model) == MOI.OPTIMAL
            maxFlux[i] = value(fluxes[i])
        end
    end

    # return
    return minFlux, maxFlux
end


################################################################################
# perform an enzyme variablility analysis (EVA) for the FBA model by minimising
# and maximising given enzyme usages of interest to see how sensitive a 
# previously constraint system is to changes the respective usages
# (for further details see Mahadevan et al., Metabolic Engeneering 2003, 
# Österberg et al., PLOS Computational Biology 2021)
#
# input parameters:
#   - model : JuMP optimisation model of the metabolic network
#   - enzymes : references to the enzymes that should be investigated with EVA
#
# output parameters:
#   - minUsage : minimal usages of given enzymes under the given conditions
#   - maxUsage : maximal usage of given enzymes under the given conditions
#
################################################################################
function calculateEVA(model::Model, enzymes::Array{VariableRef, 1} 
                      = Array{VariableRef, 1}(undef, 0))::Tuple{Array{Float64, 1}, 
                      Array{Float64, 1}}
    
    # if no indices are given check all fluxes
    isempty(enzymes) ? enzymes = model[:enzymes] : nothing

    # initialise output
    nEnzymes = length(enzymes)
    minUsage = [NaN for i = 1 : nEnzymes]
    maxUsage= [NaN for i = 1 : nEnzymes]
    
    # go through all enzymes and minimise and maximise the respective usage given
    # the input model (that should potentially have constraints on other fluxes)
    for i = 1 : nEnzymes

        # solve for minimal usage
        @objective(model, Min, enzymes[i])
        optimize!(model)
        
        # save if optimisation was successful
        if termination_status(model) == MOI.OPTIMAL
            minUsage[i] = value(enzymes[i])
        end
        
        # solve for maximal usage
        @objective(model, Max, enzymes[i])
        optimize!(model)
        
        # save if optimisation was successful
        if termination_status(model) == MOI.OPTIMAL
            maxUsage[i] = value(enzymes[i])
        end
    end
    
    # return
    return minUsage, maxUsage
end


################################################################################
# maps each enzyme to the components and reactions that it is involved 
#
# input parameters:
#   - fba : FBA model struct with all important information
#
# output parameters:
#   - enzymesToComponents : dictionary to map between the enzyme idx in the 
#                           enzymes variable/names and its idx in the components
#   - componentsToFluxes : dictionary to map between the enzyme idx in the 
#                          components to the fluxes/reactions that it is 
#                          involved in
#
################################################################################
function mapEnzymeIdx(fba::FBAModel)::Tuple{Dict{Int, Int}, Dict{Int, 
                                      Array{Int, 1}}}

    # extract relevant information from the FBA model
    model = fba.model
    fluxes = model[:fluxes]
    stochiometry = model[:stochiometry]
    components = fba.componentNames
    enzymes = fba.enzymeNames

    # initialise output dictionary
    enzymesToComponents = Dict{Int, Int}()
    componentsToFluxes = Dict{Int, Array{Int, 1}}()

    # go through enzymes and save the corresponding reaction indices 
    # (only in fluxes)
    for i = 1:size(enzymes, 1)
        nameInModel = enzymes[i, 1]
        idx = findfirst(x -> x == "prot_" * nameInModel, components)
        coefficients = normalized_coefficient.(stochiometry[idx], fluxes)
        fluxIndices = findall(x -> x != 0.0, coefficients)

        # add to dicitionary
        enzymesToComponents[i] = idx
        componentsToFluxes[idx] = fluxIndices
    end

    # return 
    return enzymesToComponents, componentsToFluxes
end


################################################################################
# calculate enzyme control coefficients (ECC) in the FBA model
# for further details see Kacser et al., Biochem Soc Trans 1995, Sanchez et al.,
# Molecular Systems Biology 2017, Österberg et al., PLOS Computational Biology
# 2021 (and many more)
#
# input parameters:
#   - fba : FBA model struct with all important information
#   - controlledFluxes : references to the variables in the JuMP model the ECC should 
#                        be calculated for
#   - perturbation : perturbation to the enzymes
#   - extreme : "Min" or "Max" for minimisation or maximisation problem 
#               respectively
#   - objectiveFlux : reference to the variable in the JuMP model that should 
#                     be optimised
#
# output parameters:
#   - eccs : matrix of approximate enzyme control coefficients for all enzymes
#            for the given reactions in fluxIdx under the given optimisation
#            problem, its columns are :
#               1) perturbed enzyme (index) (2:end) flux control 
#               coefficients for the given reactions
#   - kcatChange : binary vector that indicates if it was enough to adapt the  
#                  enzyme usage (0) or if the kcat values had to be adapted (1)
#
# NOTE: protein complexes are not respected
################################################################################
function calculateECC(fba::FBAModel, controlledFluxes::Array{VariableRef, 1}, 
                      perturbation::Float64, extreme::String, 
                      objectiveFlux::VariableRef)::Tuple{Array{Float64, 2}, 
                                                   Array{Bool, 1}}

    # extract relevant information from FBA model
    model = fba.model
    fluxes =  model[:fluxes]
    enzymes = model[:enzymes]
    stochiometry = model[:stochiometry]

    # map also the enzyme indices to the component indices for all enzymes
    # and to the reactions it is involved in
    enzymesToComponents, componentsToFluxes = mapEnzymeIdx(fba)

    # perform first optimisation and take current solution
    # as a baseline (all enzymes and the fluxes of interest)
    pFBA!(model, extreme, objectiveFlux)
    if termination_status(model) != MOI.OPTIMAL
        return nothing, nothing
    end
    initFluxes = value.(controlledFluxes)
    initEnzymes = value.(enzymes)

    # get the number of non-zero enzyme usages
    relevantIdx = findall(x -> x > precision, initEnzymes)
    relevantEnzymes = enzymes[relevantIdx]
    nEnzymes = length(relevantEnzymes)
    initEnzymes = initEnzymes[relevantIdx]

    # set objective for further optimisation
    if extreme == "Max"
        @objective(model, Max, objectiveFlux)
    elseif extreme == "Min"
        @objective(model, Min, objectiveFlux)
    end

    # prepare output matrix
    nControlledFluxes = length(controlledFluxes)
    eccMatrix = [NaN for i = 1 : nEnzymes, j = 1 : nControlledFluxes+1]

    # initialise boolean vector to indicate for which enzyme the k_cat
    # values were perturbed
    kcatChange = falses(nEnzymes)

    # go through all non-zero enzymes and perturb the respective usage to
    # investigate the influence on the given fluxes
    for i = 1 : nEnzymes

        # save bounds for setting model back to input state
        tmpLowerBound = lower_bound(relevantEnzymes[i])
        tmpUpperBound = upper_bound(relevantEnzymes[i])

        # perturb enzyme usage
        set_lower_bound(relevantEnzymes[i], initEnzymes[i] * (1 + perturbation) *
                        (1 - precision))
        set_upper_bound(relevantEnzymes[i], initEnzymes[i] * (1 + perturbation) *
                        (1 + precision))

        # re-optimise the system
        optimize!(model)

        # get the component index for the enzyme from mapping
        component = enzymesToComponents[relevantIdx[i]]

        # if the optimisation worked use the result directly for the
        # enzyme control coefficients
        if termination_status(model) == MOI.OPTIMAL
            # calculate the control coefficients (in a simplified
            # equation without neither k_cat nor e_i) according to the papers
            # mentioned above
            ecc = (value.(controlledFluxes) .- initFluxes) ./ initFluxes ./ 
                  perturbation

        # if it is not feasible perturb the turnover rate k_cat for all
        # reaction where the enzyme is involved
        else
            # mark that it was needed to change the k_cat values
            kcatChange[i] = true

            # get which reactions the enzyme is involved in
            reactionIdx = componentsToFluxes[component]
            reactions = fluxes[reactionIdx]
            nReactions = length(reactionIdx)

            # get corresponding coefficients
            tmpCoefficients = normalized_coefficient.(stochiometry[component], 
                                                      reactions)

            # set constant enzyme usage
            set_lower_bound(relevantEnzymes[i], initEnzymes[i] * (1 - precision))
            set_upper_bound(relevantEnzymes[i], initEnzymes[i] * (1 + precision))

            # perturb all turnover rates (k_cat is in the denominator of the
            # coefficient, so we divide by the perturbation)
            set_normalized_coefficient.(stochiometry[component], reactions, 
                                        tmpCoefficients ./ (1 + perturbation))
                        
            # re-optimise the system
            optimize!(model)

            if termination_status(model) == MOI.OPTIMAL
                # get the new fluxes and calculate the control
                # coefficient (in a simplified equation without neither
                # k_cat nor e_i) according to the papers mentioned above
                ecc = (value.(controlledFluxes) .- initFluxes) ./ initFluxes ./ 
                      perturbation
            else
                ecc = [NaN for n = 1 : nControlledFluxes]
            end

            # set back the coefficient
            set_normalized_coefficient.(stochiometry[component], reactions,
                                        tmpCoefficients)
        end

        # sort out the control coefficients that are below a certain precision
        # by setting their value to NaN
        tooSmallPrecision = abs.(ecc) .<= (precision / perturbation)
        ecc[tooSmallPrecision] .= NaN

        # save results
        eccMatrix[i, :] = vcat(relevantIdx[i], ecc)

        # set back the bounds to initial value for next enzyme
        set_lower_bound(relevantEnzymes[i], tmpLowerBound)
        set_upper_bound(relevantEnzymes[i], tmpUpperBound)
    end

    # accept only the rows where all control coefficients are valid
    # i.e. they are not NaN
    valid = vec(sum(isnan.(eccMatrix[:, 2:end]), dims = 2) .< 
                nControlledFluxes)
    eccMatrix = eccMatrix[valid, :]
    kcatChange = kcatChange[valid]

    return eccMatrix, kcatChange

end

