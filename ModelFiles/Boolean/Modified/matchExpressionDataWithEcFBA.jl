################################################################################
################################################################################
# goes through all files obtained from Yeastract and Lee et al., 2019, for the
# transcription factors and creates a compressed file that only includes the 
# targets that actually are included in the ecFBA model
################################################################################
################################################################################


using Printf
using DelimitedFiles
include("../../../Functions/IntegratedModel.jl")


################################################################################
@printf("Get all filenames from Yeastract...\n")
################################################################################
files = readdir("../../../ExperimentalData/YeastractLists/")
nFiles = length(files)
nTF = Int(nFiles / 2)


################################################################################
@printf("Read ecFBA model file to compare with ...\n")
################################################################################
model = "../../ecFBA/Modified/reducedEcYeast_modified_withDamage.mat"
fba = initialiseJuMPModel(model)
enzymeNames = fba.enzymeNames


################################################################################
@printf("Initialise output ...\n")
################################################################################
output = ["" for i = 1:nTF, j = 1:6]
outputFile = "TFtargets.txt"


################################################################################
@printf("Go through all files and match the target genes ...\n")
################################################################################
row = 1
for i = 1 : nFiles

    global row

    # get transcription factor and if it is positiv or negative regulation
    tf, regulation = split(files[i], "_")
    regulation = regulation[1:end-4]

    # read file
    TFtargets = readdlm("../../../ExperimentalData/YeastractLists/" * files[i], 
                        '\t', String, '\n')
    nTargets = size(TFtargets, 1)

    # go through targets and compare to the ecFBA model
    targetList = ""
    for j = 1 : nTargets

        # get target systematic name
        systematicName = TFtargets[j, 2]

        # search for the protein in the enzymes of the ecFBA model
        idx = findfirst(x -> x == systematicName, enzymeNames[:, 3])

        # add standard name to list if it was found in the model
        if !isnothing(idx)
            targetList *= enzymeNames[idx, 2] * ","
        end
    end

    # define the column in the output depending on if it is activation 
    # or inhibition
    if regulation == "activation"
        col = 2
    else
        col = 3
    end

    # fill the output
    idx = findfirst(x -> x == tf, output[:, 5])
    if isnothing(idx)

        # add the actual transcription factor name (can be more specific
        # if proteins are grouped together in the Boolean model)
        output[row, 5] = tf

        # add reference
        output[row, 6] = "Yeastract database"

        # add the correct name in the Boolean model as first column
        # (some are combined there)
        if tf == "Msn2" ||  tf == "Msn4" 
            output[row, 1] = "Msn2_4"
        elseif tf == "Rtg1" ||  tf == "Rtg3" 
            output[row, 1] = "Rtg1_3"   
        else
            output[row, 1] = tf
        end

        # fill in activation or inhibition targets as second or third column
        output[row, col] = targetList

        row += 1
    
    # if the transcription factor was already there fill only
    # in the target list as second or third column
    else
        output[idx, col] = targetList
    end  
end


################################################################################
@printf("Find the genes that occur both in activation and inhibition ...\n")
################################################################################
for i = 1 : nTF

    global output

    # get all individual genes that are upregulated by the transcription factor
    activatedGenes = split(output[i, 2], ",")
    activatedGenes = activatedGenes[1:end-1]

    # check for all those genes if they also occur in the inhibition column
    for j = 1 : length(activatedGenes)
        tmpGene = activatedGenes[j] * ","
        doubleRegulation = contains(output[i, 3], tmpGene)

        # if the gene occurs in both
        if doubleRegulation == true

            # remove the gene from both lists
            output[i, 2] = replace(output[i, 2], tmpGene => "")
            output[i, 3] = replace(output[i, 3], tmpGene => "")

            # add it to the column for both regulations
            output[i, 4] *= tmpGene
        end
    end

    # delete unnecessary commas
    !isempty(output[i, 2]) ? output[i, 2] = output[i, 2][1:end-1] : nothing
    !isempty(output[i, 3]) ? output[i, 3] = output[i, 3][1:end-1] : nothing
    !isempty(output[i, 4]) ? output[i, 4] = output[i, 4][1:end-1] : nothing

end


################################################################################
@printf("Add data from Lee et al., 2019 ...\n")
################################################################################
# read file 
TFtargets = readdlm("../../../ExperimentalData/Lee1999_Yap1_Skn7_RegulationTargets.txt", 
                    '\t', String, '\n', skipstart = 1)
nTF = size(TFtargets, 2)

# go through all transcription factors and add them to the output
for i = 1 : nTF

    global output

    # go trough all genes to specifically check if they are there in the model
    # (the second column indicated the up regulation, there is no down regulation)
    geneTargets = split(TFtargets[i, 2], ",")
    nTargets = length(geneTargets)

    # create updated target list
    targetList = ""
    for j = 1 : nTargets

        # search for the target in the enzyme list 
        target = geneTargets[j]
        idx = findfirst(x -> x == target, enzymeNames[:, 2])

        # add standard name to list if it was found in the model
        if !isnothing(idx)
            targetList *= enzymeNames[idx, 2] * ","
        end
    end

    # delete unnecessary commas
    !isempty(targetList) ? targetList = targetList[1:end-1] : nothing

    # add the new list to the output
    output = [output; [TFtargets[i, 1] targetList "" "" "" "Lee et al. 2019"]]

end


################################################################################
@printf("Save in output file ...\n")
################################################################################
open(outputFile, "w") do file
    write(file, "# activeBooleanTf\tonly activation/positive regulation")
    write(file, "\tonly inhibition/negative regulation\tboth regulations")
    write(file, "\tactualTf\treference\n")
    writedlm(file, output, "\t")
end

