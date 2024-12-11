using SBMLImporter, ModelingToolkit

pathSBML = "Yeast_CellCycle_model_HU_Berlin_Klipp_group_2022_revision.xml"

prn, cb = load_SBML(pathSBML; massaction=true)
sys = convert(ODESystem, prn.rn)

function writeVectorToFile(vector,vectorName,equations,io)
    println(io, vectorName * " = [")
    for (line_number,line) in enumerate(vector)
        linesuffix = line_number != length(vector) ? "," : ""
        prnline = equations ? replace(replace(string(line),"(t)"=>""),"Differential"=>"Differential(t)") : string(line)
        prnline = equations ? replace(prnline,r"(.*~ |.*=> )(.*)"=>s"\1 1/(compartment_1) * (\2)") : string(line) 
        println(io, prnline * linesuffix)
    end
    println(io, "]")
    println(io,"")
    
end

open("cellCycle.jl", "w") do io
    # Extract function names (states) from equations
    equations_as_string = string.(equations(sys)) 
    variable_definition_elements = replace.(equations_as_string, r"\s*~.*"=>"")
    variable_definition_elements = replace.(variable_definition_elements, r"Differential\(t\)\((.*)\)"=>s"\1")

    println(io, "ModelingToolkit.@variables t " * join(variable_definition_elements," "))
    println(io,"")
    variable_definition_elements = replace.(variable_definition_elements, "(t)"=>"")
    println(io, "state_array = [" * join(variable_definition_elements,", ") * "]")
    println(io,"")
    
    
    println(io, "ModelingToolkit.@parameters " * join(parameters(sys)," "))
    println(io,"")
    println(io, "parameter_array = [" * join(parameters(sys),", ") * "]")
    println(io,"")

    writeVectorToFile(equations(sys),"Eq",true,io)
    println(io, "@named osys = ODESystem(Eq, t, state_array, parameter_array)")
    println(io,"")
    writeVectorToFile(prn.p,"p",false,io)
    writeVectorToFile(prn.u0,"u0",true,io)

end

# check ! prn.u0 spara som den är
# check ! prn.p spara som den är
# check ! sys spara som den nu är efter massa bök! 

