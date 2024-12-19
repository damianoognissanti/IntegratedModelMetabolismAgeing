using ModelingToolkit, OrdinaryDiffEq, Plots
include("../Functions/CellCycle.jl")
include("../Functions/InitialiseCellCycleModel.jl")

Omph = initialiseCellCycleModel()
# Get indices of states and parameters that we wish to check
index_Clb2 = findfirst(x->x=="Clb2",Omph.u0_names)
index_Cdc14 = findfirst(x->x=="Cdc14",Omph.u0_names)
index_SBF = findfirst(x->x=="SBF",Omph.u0_names)
index_Cln2 = findfirst(x->x=="Cln2",Omph.u0_names)
index_Clb2 = findfirst(x->x=="Clb2",Omph.u0_names)
index_Cdc14_p = findfirst(x->x=="Cdc14_p",Omph.u0_names)
index_compartment_1 = findfirst(x->x=="compartment_1", Omph.p_names)
# Get compartment value
cp1 = Omph.p[index_compartment_1]

# Change some initial values, update the system and solve
Omph.u0[index_Clb2] = 520
Omph.u0[index_Cdc14] = 500
Omph.cell_cycle_prob = createCellCycleProblem(Omph.cell_cycle_sys, Omph.u0, Omph.tspan, Omph.p)
sol = solve(Omph.cell_cycle_prob, Tsit5())
# Make solution into array + plots in article are scaled up
sol_as_array = cp1 .* convert(Array,sol)
sol_SBF = sol_as_array[index_SBF,:]
sol_Cln2 = sol_as_array[index_Cln2,:]
index=10000
plot(sol_SBF[index:end], sol_Cln2[index:end], linecolor=2)

# Reset system and try some other initial values
Omph = initialiseCellCycleModel()
Omph.u0[index_SBF] = 100
Omph.u0[index_Cln2] = 200
Omph.cell_cycle_prob = createCellCycleProblem(Omph.cell_cycle_sys, Omph.u0, Omph.tspan, Omph.p)
sol = solve(Omph.cell_cycle_prob, Tsit5())
sol_as_array = cp1 .* convert(Array,sol)
sol_SBF = sol_as_array[index_SBF,:]
sol_Cln2 = sol_as_array[index_Cln2,:]
index=10000
plot!(sol_SBF[index:end], sol_Cln2[index:end], linecolor=1)

savefig("10b.png")

plot(sol.t, sol_SBF)
plot!(sol.t, sol_Cln2)
savefig("20b.png")

sol_Clb2 = sol_as_array[index_Clb2,:]
sol_Cdc14_p = sol_as_array[index_Cdc14_p,:]
plot(sol_Clb2, sol_Cdc14_p)
savefig("30b.png")
plot(sol.t, sol_Clb2)
plot!(sol.t, sol_Cdc14_p)
savefig("40b.png")

