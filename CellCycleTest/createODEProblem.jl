using ModelingToolkit, OrdinaryDiffEq, Plots
include("cellCycleSBFCln2.jl")
equations_as_string = string.(equations(osys)) 
variable_definition_elements = replace.(equations_as_string, r"\s*~.*"=>"")
variable_definition_elements = replace.(variable_definition_elements, r"Differential\(t\)\((.*)\)"=>s"\1")


cp1_index = findall( x -> occursin("compartment_1", x), string.(p))[1]
cp1 = p[cp1_index].second

tspan=(1,350)
oprob = ODEProblem(structural_simplify(osys), u0, tspan, p)
sol = solve(oprob, Tsit5())
sol_as_array = convert(Array,sol)
SBF_index = findall(x->x=="SBF(t)",variable_definition_elements)[1]
Cln2_index = findall(x->x=="Cln2(t)",variable_definition_elements)[1]
# plots in article are scaled up
sol_SBF = cp1.*sol_as_array[SBF_index,:]
sol_Cln2 = cp1.*sol_as_array[Cln2_index,:]
index=10000
plot(sol_SBF[index:end], sol_Cln2[index:end])

include("cellCycleClb2Cdc14.jl")
oprob = ODEProblem(structural_simplify(osys), u0, tspan, p)
sol = solve(oprob, Tsit5())
sol_as_array = convert(Array,sol)
SBF_index = findall(x->x=="SBF(t)",variable_definition_elements)[1]
Cln2_index = findall(x->x=="Cln2(t)",variable_definition_elements)[1]
# plots in article are scaled up
sol_SBF = cp1.*sol_as_array[SBF_index,:]
sol_Cln2 = cp1.*sol_as_array[Cln2_index,:]
plot!(sol_SBF[index:end], sol_Cln2[index:end])

savefig("10.png")

plot(sol.t, sol_SBF)
plot!(sol.t, sol_Cln2)
savefig("20.png")

include("cellCycleSBFCln2.jl")
oprob = ODEProblem(structural_simplify(osys), u0, tspan, p)
sol = solve(oprob, Tsit5())
sol_as_array = convert(Array,sol)
Clb2_index = findall(x->x=="Clb2(t)",variable_definition_elements)[1]
Cdc14_p_index = findall(x->x=="Cdc14_p(t)",variable_definition_elements)[1]
# plots in article are scaled up
sol_Clb2 = cp1.*sol_as_array[Clb2_index,:]
sol_Cdc14_p = cp1.*sol_as_array[Cdc14_p_index,:]
plot(sol_Clb2, sol_Cdc14_p)
savefig("30.png")
plot(sol.t, sol_Clb2)
plot!(sol.t, sol_Cdc14_p)
savefig("40.png")

