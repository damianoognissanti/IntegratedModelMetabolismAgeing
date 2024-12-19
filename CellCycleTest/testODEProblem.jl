using ModelingToolkit, OrdinaryDiffEq, Plots
include("../Functions/CellCycle.jl")
include("../Functions/InitialiseCellCycleModel.jl")

Omph = initialiseCellCycleModel()
# Get compartment value
cp1 = getCCParamValue(Omph,"compartment_1")

# Change some initial values, update the system and solve
setCCStateValue(Omph,"Clb2",520)
setCCStateValue(Omph,"Cdc14",500)
updateCCProblem(Omph)
sol = solve(Omph.cell_cycle_prob, Tsit5())
# Make solution into array + plots in article are scaled up
sol_as_array = cp1 .* convert(Array,sol)
sol_SBF = sol_as_array[getCCStateIndex(Omph,"SBF"),:]
sol_Cln2 = sol_as_array[getCCStateIndex(Omph,"Cln2"),:]
index=10000  
plot(sol_SBF[index:end], sol_Cln2[index:end], linecolor=2)
  
# Reset system and try some other initial values
resetCCProblem(Omph)  
setCCStateValue(Omph,"SBF",100)   
setCCStateValue(Omph,"Cln2",200)  
updateCCProblem(Omph)
sol = solve(Omph.cell_cycle_prob, Tsit5())
sol_as_array = cp1 .* convert(Array,sol)
sol_SBF = sol_as_array[getCCStateIndex(Omph,"SBF"),:]
sol_Cln2 = sol_as_array[getCCStateIndex(Omph,"Cln2"),:]
index=10000
plot!(sol_SBF[index:end], sol_Cln2[index:end], linecolor=1)

savefig("10b.png")

plot(sol.t, sol_SBF)
plot!(sol.t, sol_Cln2)
savefig("20b.png")

sol_Clb2 = sol_as_array[getCCStateIndex(Omph,"Clb2"),:]
sol_Cdc14_p = sol_as_array[getCCStateIndex(Omph,"Cdc14_p"),:]
plot(sol_Clb2, sol_Cdc14_p)
savefig("30b.png")
plot(sol.t, sol_Clb2)
plot!(sol.t, sol_Cdc14_p)
savefig("40b.png")

