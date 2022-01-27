# FBA parameters
const fbaPath = "../../ModelFiles/ecFBA/Modified/reducedEcYeast_modified_withDamage.mat"

# Boolean parameters
const booleanSpeciesPath = "../../ModelFiles/Boolean/Modified/species.txt"
const booleanRulesPath = "../../ModelFiles/Boolean/Modified/rules.txt"
const TFpath = "../../ModelFiles/Boolean/Modified/TFtargets.txt"

# knockout signalling proteins
const knockouts = ["Snf1", "Tor1+Tor2", "Msn2_4", "Msn2_4+Rim15", "Sch9", 
                   "Tpk1+Tpk2+Tpk3", "Yap1", "Skn7", "Yap1+Skn7", 
                   "Yap1+Skn7+Msn2_4"]

# ODE parameters
const parameterFile = "../Fig2A/f0r0grid.txt"
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
const regulationFactor = 0.04
const growthFlexibility = 0.5


