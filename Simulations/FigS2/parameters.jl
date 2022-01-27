# define where the FBA model is stored
const fbaPath = "../../ModelFiles/ecFBA/Modified/reducedEcYeast_modified_withDamage.mat"

# define where the Boolean model is stored
const booleanSpeciesPath = "../../ModelFiles/Boolean/Modified/species.txt"
const booleanRulesPath = "../../ModelFiles/Boolean/Modified/rules.txt"
const TFpath = "../../ModelFiles/Boolean/Modified/TFtargets.txt"

# set simulation steps
const minGrowth = 0.0
const maxGrowth = 0.4
const nSteps = 100

# signalling parameters
const glucoseThreshold = 3.2914
const damageThreshold = 0.001
const trxThreshold = 0.000000002

# regulation parameter
const regulationFactor = 0.04
const glucoseFlexibility = 0.15
