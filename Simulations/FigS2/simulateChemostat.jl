################################################################################
################################################################################
# simulation of the chemostat experiment (van Hoek et al., Applied and 
# Environmental Microbiology 1998) with a regulated enzyme-constraint metabolic 
# model, extended by ROS/RNS producing reactions in the metabolic model 
# and a Boolean stress signalling layer for regulation
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
@printf("Initialise models and outputs ... \n")
################################################################################
# ecFBA
fba = initialiseJuMPModel(fbaPath)
prepareChemostatExperiment!(fba)

# Boolean
boolean = initialiseBooleanModel(booleanSpeciesPath, booleanRulesPath)
targets = readTranscriptionalTargets(TFpath, fba.enzymeNames)

# output
outputFile = "exchangeFluxes.txt"
output = [NaN for i = 1:8, j = 1:nSteps+1]
damageFile = "damage.txt"
outputDamage = [NaN for i = 1:9, j = 1:nSteps+1]


################################################################################
@printf("Extract relevant reactions and components from FBA model ... \n")
################################################################################
# extract the JuMP optimisation variables
model = fba.model
reactionNames = fba.reactionNames
fluxes = model[:fluxes]
enzymes = model[:enzymes]
pool = model[:pool]

# extract growth reaction, GAM and glucose uptake
growthFlux = fluxes[findfirst(x -> x == "Growth", reactionNames)]
GAMConstraint = model[:stochiometry][findfirst(x -> 
                      x == "Maintainance for growth [cytosol]", fba.componentNames)]
glucoseUptakeFlux = fluxes[findfirst(x -> x == "Uptake of glucose", reactionNames)]

# get indices for exchange reactions
# (reversible reactions are 0 so the netto flux is given by those reactions only)
exchangeReactions = ["Uptake of O2", "Production of CO2", "Uptake of glucose", 
                     "Production of ethanol", "Production of acetate", 
                     "Production of glycerol"]
exchangeFluxes = fluxes[map(x -> findfirst(y -> y == x, reactionNames), 
                        exchangeReactions)]

# get other reactions responsible for damage production
damageReactions = ["superoxide production (mitochondria)", 
                   "superoxide oxidoreductase (arm, mitochondria)", 
                   "glutathione peroxidase (GPX3, mitochondria)", 
                   "thioredoxin peroxidase (TRX3, mitochondria)", 
                   "hydrogen peroxide catalase (CTA1, mitochondria)", 
                   "protein damage exchange (mitochondria)", 
                   "superoxide oxidoreductase (SOD1, cytosol)", 
                   "hydrogen peroxide transport (reversible)",
                   "hydrogen peroxide transport", 
                   "glutathione peroxidase (arm, cytosol)",
                   "thioredoxin peroxidase (arm, cytosol)", 
                   "hydrogen peroxide catalase (CTT1, cytosol)", 
                   "protein damage exchange (cytosol)"]
damageFluxes = fluxes[map(x -> findfirst(y -> y == x, reactionNames), 
                      damageReactions)]
damageProductionFluxes = damageFluxes[[6, 13]]

# get indices for exchange reactions
trxReactions = ["draw_prot_P22217", "draw_prot_P22803"]
trxFluxes = enzymes[map(x -> findfirst(y -> y == x, reactionNames), 
                   trxReactions) .- length(fluxes)]
    
# get also the Boolean components, in particular the ones needed for the 
# signalling
booleanGlucose = boolean.components[findfirst(y -> y.name == "exGlc",
                                    boolean.components)]
booleanPeroxoide= boolean.components[findfirst(y -> y.name == "H2O2",
                                     boolean.components)]
booleanTrx = boolean.components[findfirst(y -> y.name == "Trx1_2",
                                boolean.components)]

                                
################################################################################
@printf("Save initial lower and upper bounds ... \n")
################################################################################
glucoseUpperBound = upper_bound(glucoseUptakeFlux)
enzymeLowerBounds = lower_bound.(enzymes)
enzymeUpperBounds = upper_bound.(enzymes)


################################################################################
@printf("Perform optimisation ... \n")
################################################################################
# go through all growth rates and minimise the glucose uptake rate
# with parsimonious FBA
delta = (maxGrowth - minGrowth) / nSteps
for i = 1 : nSteps+1

    # unconstrain all enzymes and the glucose uptake
    set_upper_bound(glucoseUptakeFlux, glucoseUpperBound)
    set_lower_bound.(enzymes, enzymeLowerBounds)
    set_upper_bound.(enzymes, enzymeUpperBounds)

    # set temporary growth rate
    tmpGrowth = minGrowth + (i - 1) * delta
    set_lower_bound(growthFlux, tmpGrowth * (1 - precision))
    set_upper_bound(growthFlux, tmpGrowth * (1 + precision))
    
    # adapt the growth-associated maintenance (GAM) depending on the growth rate
    updateGAMValue!(growthFlux, GAMConstraint, tmpGrowth)

    # solve the parsimonious FBA
    pFBA!(model, "Min", glucoseUptakeFlux)
    
    if termination_status(model) == MOI.OPTIMAL
        # set nutrients in Boolean model according to solution of FBA
        triggerSignalling!([glucoseUptakeFlux], glucoseThreshold, booleanGlucose)

        # set stress in Boolean model according to solution of FBA
        triggerSignalling!(damageProductionFluxes, damageThreshold, 
                           booleanPeroxoide)
        triggerSignalling!(trxFluxes, trxThreshold, booleanTrx)

        # run Boolean model to obtain logical steady state for the regulation
        runBooleanModel!(boolean)

        # calculate ranks for gene expression given the transcription factor 
        # activity in the Boolean model
        ranks = getExpressionRank(boolean.components, fba.enzymeNames, targets)
        
        # and constrain it with some flexibility 
        # (factor defined in ../../Functions/Intergation.jl)
        set_upper_bound(glucoseUptakeFlux, value(glucoseUptakeFlux) * 
                        (1 + glucoseFlexibility))

        # regulate FBA according to the ranks from the Boolean model
        regulateFBA!(model, ranks, regulationFactor)

        # solve the regulated parsimonious FBA
        pFBA!(model, "Min", glucoseUptakeFlux)

        # save exchange fluxes and enzyme usage in output matrix
        output[:, i] = [tmpGrowth; value.(exchangeFluxes); value(pool)]

        # save fluxes related to ROS production and removal,
        # add them together if they belong together
        damage = value.(damageFluxes)
        damage[3] += damage[4] + damage[5]
        damage[8] -= damage[9]
        damage[10] += damage[11] + damage[12]
        outputDamage[:, i] = [tmpGrowth; damage[[1, 2, 3, 6, 7, 8, 10, 13]]]
    end
end


################################################################################
@printf("Save in output file ...\n")
################################################################################
pretext = 
"# EXCHANGE FLUXES IN REGULATED ecFBA MODEL
# FBA model $fbaPath with objective minGlc
# Boolean model $booleanRulesPath, $booleanSpeciesPath
# TF layer $TFpath
growth rate\tO2 uptake\tCO2 production\tglucose uptake\t" *
"ethanol production\tacetate production\tglycerol production\tenzyme pool
"

open(outputFile, "w") do file
    write(file, pretext)
    writedlm(file, transpose(output), "\t")
end

pretext = 
"# DAMAGE IN REGULATED ecFBA MODEL
# FBA model $fbaPath with objective minGlc
# Boolean model $booleanRulesPath, $booleanSpeciesPath
# TF layer $TFpath
growth rate\tsuperoxide production\tperoxide production (m)\t" *
"peroxide removal (m)\tdamage (m)\tperoxide production (c)\t" *
"peroxide transport\tperoxide removal (c)\tdamage (c)
"

open(damageFile, "w") do file
    write(file, pretext)
    writedlm(file, transpose(outputDamage), "\t")
end
