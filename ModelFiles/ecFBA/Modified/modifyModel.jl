################################################################################
################################################################################
# slightly adapts the published reduced ecYeast model from Österberg et al., 
# PLOS Computational Biology 2021, to be ordered and have the compartments 
# included in the component names. Just as in the published model, an exchange 
# reaction for pyruvate was added.
# Only a minimal version that is needed for the simulations is saved.
################################################################################
################################################################################


using Printf
using DelimitedFiles
using MAT
using SparseArrays


################################################################################
@printf("Read original ecYeast model ... \n")
################################################################################
filepath = "../Original/reducedEcYeast.mat"
file = matread(filepath)
modelName = collect(keys(file))[1]
model = file[modelName]


################################################################################
@printf("Extract model information ... \n")
################################################################################
componentNames = model["metNames"][:]
reactionNames = model["rxnNames"][:]
reactionPathways = model["subSystems"][:]
stochiometry = model["S"]
changes = model["b"]
nComponents, nReactions = size(stochiometry)
lowerBounds = model["lb"][:]
upperBounds = model["ub"][:]
enzymesModel = model["enzymes"][:]
enzymesSystematic = model["enzGenes"][:]
enzymesStandard = model["enzNames"]
enzymePathways = model["pathways"][:]


################################################################################
@printf("Add compartment to each component ... \n")
################################################################################
# add compartment to the component names to make them unique,
# unless it is an enzyme or pseudo component for arm reactions
compartmentNames = model["compNames"][:] 
compartmentIdx = Int.(model["metComps"][:])
for i = 1 : nComponents
    if occursin("prot_", componentNames[i]) 
        nothing
    else
        componentNames[i] *= " [" .* compartmentNames[compartmentIdx[i]] .* "]"
    end
end


################################################################################
@printf("Add two reactions for pyruvate exchange ... \n")
################################################################################
reactionNames = [reactionNames; "Production of pyruvate"; 
                 "Production of pyruvate (reversible)"]
reactionPathways = [reactionPathways; [["Other"], ["Other"]]]

# construct the new reaction
pyruvate = findfirst(x -> x == "pyruvate [cytosol]", componentNames)
productionReaction = zeros(Float64, nComponents)
productionReaction[pyruvate] = -1.0
uptakeReaction = zeros(Float64, nComponents)
uptakeReaction[pyruvate] = +1.0

# update the model
# (but don't allow uptake for now, to keep it consistent with 
#  Österberg et al., 2021)
stochiometry = [stochiometry productionReaction uptakeReaction]
lowerBounds = [lowerBounds; 0.0; 0.0]
upperBounds = [upperBounds; 1000.0; 0.0]


################################################################################
@printf("Remove some reactions that are double ... \n")
################################################################################
reactionsToDelete = []

# rename uptake of acetate and delete its reversible reaction
idx = findfirst(x -> x == "Uptake of acetate", reactionNames)
reactionNames[idx] = "Production of acetate (reversible)"
idx = findfirst(x -> x == "Uptake of acetate (reversible)", reactionNames)
reactionsToDelete = [reactionsToDelete; idx]

# rename uptake of ethanol and delete its reversible reaction
idx = findfirst(x -> x == "Uptake of ethanol", reactionNames)
reactionNames[idx] = "Production of ethanol (reversible)"
idx = findfirst(x -> x == "Uptake of ethanol (reversible)", reactionNames)
reactionsToDelete = [reactionsToDelete; idx]

# rename uptake of H2O and delete its reversible reaction
idx = findfirst(x -> x == "Uptake of f H2O", reactionNames)
reactionNames[idx] = "Production of H2O (reversible)"
idx = findfirst(x -> x == "Uptake of f H2O (reversible)", reactionNames)
reactionsToDelete = [reactionsToDelete; idx]

# rename production of Phosphate and delete its reversible reaction
idx = findfirst(x -> x == "Production of Phosphate", reactionNames)
reactionNames[idx] = "Uptake of Phosphate (reversible)"
idx = findfirst(x -> x == "Production of Phosphate (reversible)", reactionNames)
reactionsToDelete = [reactionsToDelete; idx]

# delete reactions 
reactionsToKeep = setdiff(1:size(stochiometry, 2), reactionsToDelete)    
stochiometry = stochiometry[:, reactionsToKeep] 
reactionNames = reactionNames[reactionsToKeep]
lowerBounds = lowerBounds[reactionsToKeep] 
upperBounds = upperBounds[reactionsToKeep]
reactionPathways = reactionPathways[reactionsToKeep]


################################################################################
@printf("Order ... \n")
################################################################################
# get the right order
enzymeIdx = contains.(reactionNames, "draw_prot")
poolIdx = contains.(reactionNames, "prot_pool")
fluxIdx = (enzymeIdx + poolIdx) .< 1

order = [findall(x -> x == 1, fluxIdx); findall(x -> x == 1, enzymeIdx);
         findall(x -> x == 1, poolIdx)]
         
# update orders
reactionNames = reactionNames[order]
stochiometry = stochiometry[:, order]
lowerBounds = lowerBounds[order]
upperBounds = upperBounds[order]

# update pathways 
reactionPathways = join.(reactionPathways, ", ")
reactionPathways[findall(x -> x == 1, enzymeIdx)] .= "Protein usage"
reactionPathways[findall(x -> x == 1, poolIdx)] .= "Protein pool"
reactionPathways = reactionPathways[order]

# update empty entries in pathways to other
reactionPathways[reactionPathways .== ""] .= "Other"

# and update some pathways to have capital letter
reactionPathways[reactionPathways .== "pentose phosphate"] .= "Pentose phosphate"


################################################################################
@printf("Prepare string matrix for enzymes ... \n")
################################################################################
# exchange one enzyme name in the enzNames list which were filled with 
# systematic names in the original
# (note: enzyme 119 does not have a standard name and there the systematic 
# name remains for both)
enzymesStandard[98] = "SDH9"

# make sure that the order in the enzymes names is the same as the order of the
# enzyme reactions in reactionNames
nFluxes = sum(fluxIdx)
namesInReactions = map(x -> x[11:end], reactionNames[nFluxes+1:end-1])
orderInEnzymes = map(x -> findfirst(y -> y == x, enzymesModel), 
                     namesInReactions)

# update order
enzymesModel = enzymesModel[orderInEnzymes]
enzymesStandard= enzymesStandard[orderInEnzymes]
enzymesSystematic = enzymesSystematic[orderInEnzymes]

enzymePathways = enzymePathways[orderInEnzymes]

# put the names together in a matrix
enzymes = [enzymesModel enzymesStandard enzymesSystematic]


################################################################################
@printf("Save new minimal MAT model ... \n")
################################################################################
newModel = Dict("metNames" => componentNames, "rxnNames" => reactionNames,
                "S" => stochiometry, "b" => changes, "lb" => lowerBounds,
                "ub" => upperBounds, "enzymes" => enzymes, 
                "enzymePathways" => enzymePathways, 
                "rxnPathways" => reactionPathways)
newFile = Dict("reducedEcYeast" => newModel)
matwrite("reducedEcYeast_modified.mat", newFile)