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
filename = "enzymeOverexpressions.txt"
data = readdlm(filename, '\t', Any, '\n', skipstart = 7)


################################################################################
@printf("Extract enzymes ... \n")
################################################################################
# remove the wildtype
data = data[data[:, 1] .!= "wildtype", :]

# get the ones that increase lifespan but only in a specific phase
subset1 = data[(data[:, 6] .> 1) .& (data[:, 4] .!= "complete"), :]
idx = [] 
for k = 1 : size(subset1, 1)

    global idx

    tmpData = data[(data[:, 1] .== subset1[k, 1]) .& (data[:, 4] .== "complete"), :]
    tmpData[6] <= 1 ? idx = [idx; k] : nothing
end
subset1 = subset1[idx, :]
reason1 = "increased lifespan but only in specific phase"

# get the ones that decreased lifespan but only in a specific phase
subset2 = data[(data[:, 6] .< 1) .& (data[:, 4] .!= "complete"), :]
idx = [] 
for k = 1 : size(subset2, 1)

    global idx

    tmpData = data[(data[:, 1] .== subset2[k, 1]) .& (data[:, 4] .== "complete"), :]
    tmpData[6] >= 1 ? idx = [idx; k] : nothing
end
subset2 = subset2[idx, :]
reason2 = "decreased lifespan but only in specific phase"

# get the ones that increase lifespan even though used less
subset3 = data[(data[:, 7] .< 0.999) .& (data[:, 6] .> 1), :]
reason3 = "increased lifespan and used less"

# get the ones that increase lifespan because used more
subset4 = data[(data[:, 7] .> 1.001) .& (data[:, 6] .> 1), :]
reason4 = "increased lifespan and used more"

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
# from overexpression experiment 
enzyme(s)\tstandard name\tpathway\toverexpression span\trls\trls change\tenzyme change\t" *
"lcc\ttime phase 1\trls phase 1\tdamage phase 1\ttime phase 2\trls phase 2\t" *
"damage phase 2\ttime phase 3\trls phase 3\tdamage phase 3\tstatus\twhy
"

open(outputFile, "w") do file
    write(file, pretext)
    writedlm(file, output, "\t")
end 