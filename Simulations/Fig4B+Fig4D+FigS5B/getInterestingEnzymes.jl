################################################################################
################################################################################
# extract interesting enzymes from deletion experiment
################################################################################
################################################################################

using Printf
using DelimitedFiles
using Statistics
include("../../Functions/IntegratedModel.jl")
         

################################################################################
@printf("Initialise the output ... \n")
################################################################################
outputFile = "interestingEnzymes.txt"


################################################################################
@printf("Read file ... \n")
################################################################################
filename = "enzymeDeletions.txt"
data = readdlm(filename, '\t', Any, '\n', skipstart = 7)


################################################################################
@printf("Extract enzymes ... \n")
################################################################################
# remove the wildtype
data = data[data[:, 1] .!= "wildtype", :]

# get the ones that are used more even though deleted at some point
subset1 = data[data[:, 7] .> 1.001, :]
reason1 = "more used despite deletion"

# get the ones that increase lifespan
subset2 = data[data[:, 6] .> 1, :]
reason2 = "increased lifespan"

# get the ones that increase lifespan but only in a specific phase
subset3 = data[(data[:, 6] .> 1) .& (data[:, 4] .!= "complete"), :]
idx = [] 
for k = 1 : size(subset3, 1)

    global idx

    tmpData = data[(data[:, 1] .== subset3[k, 1]) .& (data[:, 4] .== "complete"), :]
    tmpData[6] <= 1 ? idx = [idx; k] : nothing
end
subset3 = subset3[idx, :]
reason3 = "increased lifespan but only in specific phase"

# get the ones that decreased lifespan but only in a specific phase
subset4 = data[(data[:, 6] .< 1) .& (data[:, 4] .!= "complete"), :]
idx = [] 
for k = 1 : size(subset4, 1)

    global idx

    tmpData = data[(data[:, 1] .== subset4[k, 1]) .& (data[:, 4] .== "complete"), :]
    tmpData[6] >= 1 ? idx = [idx; k] : nothing
end
subset4 = subset4[idx, :]
reason4 = "decreased lifespan but only in specific phase"

# save in output
output = [subset1 [reason1 for i = 1:size(subset1, 1)];
          subset2 [reason2 for i = 1:size(subset2, 1)];
          subset3 [reason3 for i = 1:size(subset3, 1)];
          subset4 [reason4 for i = 1:size(subset4, 1)]]

# sort
output = sortslices(output, dims = 1, by = x -> (x[end], x[4], x[6]))


################################################################################
@printf("Save in output file ...\n")
################################################################################
pretext =
"# INTERESTING ENZYMES
# from deletion experiment 
enzyme(s)\tstandard name\tpathway\tdeletion span\trls\trls change\tenzyme change\t" *
"lcc\ttime phase 1\trls phase 1\tdamage phase 1\ttime phase 2\trls phase 2\t" *
"damage phase 2\ttime phase 3\trls phase 3\tdamage phase 3\tstatus\twhy
"

open(outputFile, "w") do file
    write(file, pretext)
    writedlm(file, output, "\t")
end 