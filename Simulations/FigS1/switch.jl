################################################################################
################################################################################
# simulation of a Boolean model of nutrient and stress signalling,
# to test if it switches correctly between the nutrient availabilities
# and the stress conditions
################################################################################
################################################################################


using Printf
include("../../Functions/IntegratedModel.jl")


################################################################################
@printf("Define parameters ... \n")
################################################################################
include("parameters.jl")


################################################################################
@printf("Initialise Boolean model ... \n")
################################################################################
boolean = initialiseBooleanModel(booleanSpecies, booleanRules)
names = map(x -> x.name, boolean.components)
pathways = map(x -> x.pathway, boolean.components)
glucoseIdx = findfirst(x -> x == "exGlc", names)
nitrogenIdx = findfirst(x -> x == "NH3", names)
h2o2Idx = findfirst(x -> x == "H2O2", names)
trxIdx = findfirst(x -> x == "Trx1_2", names)


################################################################################
@printf("Go through initial conditions and find steady state after switching ... \n")
################################################################################
for n = 1 : length(glucose)

    ############################################################################
    # Initialise output file
    ############################################################################
    folder = @sprintf("glc%i_nit%i_trx%i", glucose[n], nitrogen[n], 
                      trx[n])
    mkpath(folder)
    outputFile = folder * "/activity.txt"

    ############################################################################
    # Run until steady state of the initial nutrient availability reached ... \n")
    ############################################################################
    # set glucose and nitrogen availability and stress level
    boolean.components[glucoseIdx].present = glucose[n]
    boolean.components[nitrogenIdx].present = nitrogen[n]
    boolean.components[trxIdx].present = trx[n]

    # set initial state of switch
    boolean.components[h2o2Idx].present = h2o2_switch[1]

    # iterate over the Boolean model until the steady state is reached
    runBooleanModel!(boolean)

    # initialise the output
    output = [names pathways map(x -> x.active, boolean.components)]
    iteration = 0
    switch = zeros(Int64, length(h2o2_switch) - 1)

    ############################################################################
    # Switch nutrient availability and run again
    ############################################################################
    for k = 2 : length(h2o2_switch)

        # switch h2o2 state
        boolean.components[h2o2Idx].present = h2o2_switch[k]

        # iterate over the Boolean model step by step until the steady state 
        # is reached and save the state in each iteration
        steadyState = false
        while !steadyState
            iteration += 1
            steadyState, ~ = runBooleanModel!(boolean, 1)

            # extract activities and save in output
            activities = map(x -> x.active, boolean.components)
            output = [output activities]
        end

        switch[k-1] =  iteration
    
    end

    switch = switch[1:end-1]

    ############################################################################
    # Save in output file
    ############################################################################
    open(outputFile, "w") do file
        write(file, "# BOOLEAN MODEL ITERATIONS WITH SWITCH\n")
        write(file, "# model defined in $booleanSpecies and $booleanRules\n")
        write(file, "# glucose $(glucose[n]), nitrogen $(nitrogen[n]), ")
        write(file, "trx $(trx[n])\n")
        write(file, "# h2o2 switch: $h2o2_switch\n")
        write(file, "# switch at iteration $switch\n")
        write(file, "# name\tpathway\tinitActivity\tactivity($iteration iterations)\n")
        writedlm(file, output, "\t")
    end
end