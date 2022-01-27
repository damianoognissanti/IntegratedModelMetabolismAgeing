################################################################################
################################################################################
# estimates the fraction of the proteome that is included in the model
################################################################################
################################################################################


using Printf
using DelimitedFiles


################################################################################
@printf("Read in data from the integrated PaxDB database ...\n")
################################################################################
filenamePaxDBData = "../../../../ExperimentalData/paxDBdatabase.txt"

# read data file without header (13 lines)
# there are both strings and floats in the data, so we use the type Any
paxDBData = readdlm(filenamePaxDBData, '\t', Any, '\n', skipstart = 13)

# remove numbers in front of protein name in original data
paxDBData[:, 2] = map(x -> x[6:end], paxDBData[:, 2])


################################################################################
@printf("Read list of genes that are included in the model ...\n")
################################################################################
filenameEnzymes = "geneList.txt"

# read data file without header (4 lines)
enzymes = readdlm(filenameEnzymes, '\t', Any, '\n', skipstart = 4)


################################################################################
@printf("Extract abundances of included enzymes ...\n")
################################################################################
nEnzymes = size(enzymes, 1)
abundances = zeros(nEnzymes)
count = 0

for i = 1:nEnzymes

    global count 

    name = enzymes[i, 1]
    idx = findfirst(x -> x == name, paxDBData[:, 2])

    if isnothing(idx)
        @printf("Enzyme %s not found. \n", name) 
        count += 1
    else
        abundances[i] = float(paxDBData[idx, 3])
    end
end


################################################################################
@printf("Calculate fractional abundances ...\n")
################################################################################
abundanceOld = sum(abundances[enzymes[:, 2] .== 0])
abundanceNew = sum(abundances[enzymes[:, 2] .== 1])

# assume that the abundances cover 96% of the total protein abundance
totalAbundance = sum(paxDBData[:, 3]) / 0.96 


################################################################################
@printf("Save ...\n")
################################################################################
open("fractionalAbundances.txt", "w") do file
    write(file, "Original files:\n")
    write(file, "\"$filenamePaxDBData\" \n")
    write(file, "\"$filenameEnzymes\" \n\n")
    write(file, "$count enzymes were not found in the data base \n")
    write(file, "The fractional abundance of the enymes used previously ")
    write(file, @sprintf("(marked 0) is %.4f, that of the new enzymes (marked 1) is ", 
          abundanceOld / totalAbundance))
    write(file, @sprintf("%.4f, and consequently the total new fractional abundance ", 
          abundanceNew / totalAbundance))
    write(file, @sprintf("is %.4f. \n", (abundanceOld + abundanceNew) / totalAbundance))
end
