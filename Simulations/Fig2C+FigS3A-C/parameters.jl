# FBA parameters
const fbaPath = "../../ModelFiles/ecFBA/Modified/reducedEcYeast_modified_withDamage.mat"

# Boolean parameters
const booleanSpeciesPath = "../../ModelFiles/Boolean/Modified/species.txt"
const booleanRulesPath = "../../ModelFiles/Boolean/Modified/rules.txt"
const TFpath = "../../ModelFiles/Boolean/Modified/TFtargets.txt"

# ODE parameters
const formationRate0 = collect(0.0:0.0001:0.002)
const repairRate0 = collect(0.0:0.002:0.02)
const maxGrowthRate = 0.35

# set initial conditions for ODE
# assume 46% of the drymass are functional proteins (Famili et al., PNAS 2003), 
# there is no damage and the initial drymass is 1 since we are only interested
# in its relativ change, not the absolute 
const P = 0.46
const D = 0.0
const mass = 1.0
const timestep = 0.1
const maxTimesteps = 1000

# parameters for division
const sizeProportion = 0.64
const retention = 0.3

# signalling parameters
const nDelay = 5
const glucoseThreshold = 3.2914
const damageThreshold = 0.001
const trxThreshold = 0.000000002

# regulation parameters
const regulationFactor = [0.0, 0.02, 0.04]
const growthFlexibility = 0.5

