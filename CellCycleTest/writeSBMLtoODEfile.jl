using SBMLImporter, ModelingToolkit

pathSBML = "Yeast_CellCycle_model_HU_Berlin_Klipp_group_2022_revision.xml"

prn, cb = load_SBML(pathSBML; massaction=true)
sys = convert(ODESystem, prn.rn)

# Extract states, parameters and equations
states = [string(state.first) for state in prn.u0]
params = [string(param.first) for param in prn.p]
eqs = equations(sys)
# Extract state and parameter values
statesval = [string(state.second) for state in prn.u0]
paramsval = [string(param.second) for param in prn.p]

open("CellCycle.jl", "w") do io
    println(io, "function createCellCycleModel()")
    println(io, "")

    println(io, "   ModelingToolkit.@variables t " * join(states," "))
    println(io,"")
    println(io, "   cell_cycle_state_array = [" * join(replace.(states, "(t)"=>""),", ") * "]")
    println(io,"")
    println(io, "   ModelingToolkit.@parameters " * join(params," "))
    println(io,"")
    println(io, "   cell_cycle_parameter_array = [" * join(params,", ") * "]")
    println(io,"")
    println(io,"    cellCycleEq = [")
    for (i,eq) in enumerate(eqs)
        println(io, "       ",replace(replace(string(eq.lhs),"(t)"=>""),"Differential"=>"Differential(t)")," ~ 1/(compartment_1) * (",replace(string(eq.rhs),"(t)"=>""),")", i != length(eqs) ? "," : "")
    end
    println(io,"    ]")

    println(io, "   @named cell_cycle_sys = ODESystem(cellCycleEq, t, cell_cycle_state_array, cell_cycle_parameter_array)")
    println(io,"")
    println(io, "   return cell_cycle_sys")
    println(io, "   end")
    println(io,"")
    println(io,"function createCellCycleProblem(cell_cycle_sys, state, tspan, param)")
    println(io, "   ModelingToolkit.@variables t " * join(states," "))
    println(io,"")
    println(io, "   ModelingToolkit.@parameters " * join(parameters(sys)," "))
    println(io,"")
    println(io, "   p = [")
    for (i,p) in enumerate(params)
        println(io, "       ",p," => param[",i,"]",i != length(params) ? "," : "")
    end
    println(io, "   ]")
    println(io,"")
    println(io, "   u0 = [")
    for (i,u0) in enumerate(states)
        println(io, "       ",replace(u0,"(t)"=>"")," => 1/(compartment_1) * state[",i,"]",i != length(states) ? "," : "")
    end
    println(io, "]")
    println(io,"")
    println(io,"    return ODEProblem(structural_simplify(cell_cycle_sys), u0, tspan, p)")
    println(io,"")
    println(io,"end")
end

function writeVector(io,vector,name,isstr)
    println(io, "   ",name," = [")
    for (ind,val) in enumerate(vector)
        val = isstr ? "\"" * replace(val,"(t)"=>"") * "\"" : val
        println(io, "       ",val,ind != length(vector) ? "," : "")
    end
    println(io,"    ]")
    println(io,"")
end

open("InitialiseCellCycleModel.jl", "w") do io
    println(io,"mutable struct CellCycleModel")
    println(io,"    cell_cycle_sys::ODESystem")
    println(io,"    cell_cycle_prob::ODEProblem")
    println(io,"    p::Array")
    println(io,"    u0::Array")
    println(io,"    p_names::Array")
    println(io,"    u0_names::Array")
    println(io,"    tspan::Tuple")
    println(io,"    porig::Array")
    println(io,"    u0orig::Array")
    println(io,"end")
    println(io,"")

    println(io,"function initialiseCellCycleModel()")
    println(io, "")
    writeVector(io,paramsval,"porig",false)
    writeVector(io,statesval,"u0orig",false)
    writeVector(io,params,"p_names",true)
    writeVector(io,states,"u0_names",true)
    println(io,"cell_cycle_sys = createCellCycleModel()")
    println(io,"    tspan=(1,350)")
    println(io,"    cell_cycle_prob = createCellCycleProblem(cell_cycle_sys, u0orig, tspan, porig)  ")
    println(io,"")
    println(io,"    return CellCycleModel(cell_cycle_sys, cell_cycle_prob, porig, u0orig, p_names, u0_names, tspan, porig, u0orig)")
    println(io,"end")
end
