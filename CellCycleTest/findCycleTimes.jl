using ModelingToolkit, OrdinaryDiffEq, Plots, Statistics
include("cellCycle.jl")
equations_as_string = string.(equations(osys)) 
variable_definition_elements = replace.(equations_as_string, r"\s*~.*"=>"")
variable_definition_elements = replace.(variable_definition_elements, r"Differential\(t\)\((.*)\)"=>s"\1")

cp1_index = findall( x -> occursin("compartment_1", x), string.(p))[1]
cp1 = p[cp1_index].second

tspan=(1,350)
oprob = ODEProblem(structural_simplify(osys), u0, tspan, p)
sol = solve(oprob, Tsit5())
sol_as_array = convert(Array,sol)
Cln2_index = findall(x->x=="Cln2(t)",variable_definition_elements)[1]
sol_Cln2 = cp1.*sol_as_array[Cln2_index,:]
cyc = sol_Cln2 .> quantile(sol_Cln2)[2]
cycle_start_index = findall(diff([0; cyc]) .== 1)
cycle_start = sol.t[cycle_start_index]

