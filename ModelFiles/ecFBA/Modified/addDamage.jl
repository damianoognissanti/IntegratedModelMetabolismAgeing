################################################################################
################################################################################
# creates a new MAT model from the modified reduced ecYeast model
# that now also includes the damage reactions
################################################################################
################################################################################


using Printf
using MAT
using SparseArrays


################################################################################
@printf("Read minimal reduced ecYeast model ... \n")
################################################################################
filepath = "reducedEcYeast_modified.mat"
file = matread(filepath)
modelName = collect(keys(file))[1]
model = file[modelName]


################################################################################
@printf("Extract model information ... \n")
################################################################################
componentNames = model["metNames"][:]
reactionNames = model["rxnNames"][:]
reactionPathways = model["rxnPathways"][:]
stochiometry = model["S"]
changes = model["b"]
nComponents, nReactions = size(stochiometry)
lowerBounds = model["lb"][:]
upperBounds = model["ub"][:]
enzymes = model["enzymes"]
enzymePathways = model["enzymePathways"][:]


################################################################################
@printf("Update component names ... \n")
################################################################################
newComponents = ["O2 [mitochondria]",
                 "superoxide [mitochondria]",
                 "superoxide leakage minimum [mitochondria]", 
                 "superoxide leakage maximum [mitochondria]",
                 "H2O2 [mitochondria]",
                 "pmet_superoxide [mitochondria]",
                 "prot_P00445",
                 "prot_P00447",
                 "hydroxyl radical [mitochondria]",
                 "glutathione [mitochondria]",
                 "glutathione disulfide [mitochondria]",
                 "thioredoxin [mitochondria]",
                 "thioredoxin disulfide [mitochondria]",
                 "prot_P40581",
                 "prot_P25372",
                 "prot_P41921",
                 "prot_P38816",
                 "prot_P15202",
                 "iron(2+) [mitochondria]",
                 "iron(3+) [mitochondria]",
                 "protein damage [mitochondria]",
                 "H2O2 [cytosol]",
                 "superoxide [cytosol]",
                 "hydroxyl radical [cytosol]",
                 "h+ [cytosol]",
                 "glutathione [cytosol]",
                 "glutathione disulfide [cytosol]",
                 "pmet_glutathione [cytosol]",
                 "thioredoxin [cytosol]",
                 "thioredoxin disulfide [cytosol]",
                 "pmet_thioredoxin [cytosol]",
                 "prot_P36014",
                 "prot_P38143",
                 "prot_P22217",
                 "prot_P22803",
                 "prot_P29509",
                 "prot_P06115",
                 "iron(2+) [cytosol]",
                 "iron(3+) [cytosol]",
                 "protein damage [cytosol]"]
nNewComponents = length(newComponents)
componentNames = [componentNames; newComponents]


################################################################################
@printf("Update enzyme names ... \n")
################################################################################
newEnzymes = ["P00445" "SOD1" "YJR104C";
              "P00447" "SOD2" "YHR008C";
              "P40581" "GPX3" "YIR037W";
              "P25372" "TRX3" "YCR083W";
              "P41921" "GLR1" "YPL091W";
              "P38816" "TRR2" "YHR106W";
              "P15202" "CTA1" "YDR256C";
              "P36014" "GPX1" "YKL026C";
              "P38143" "GPX2" "YBR244W";   
              "P22217" "TRX1" "YLR043C";
              "P22803" "TRX2" "YGR209C";
              "P29509" "TRR1" "YDR353W";
              "P06115" "CTT1" "YGR088W"]
enzymes = [enzymes; newEnzymes]
nNewEnzymes= size(newEnzymes, 1)
kegg = ["sce00NNN  Damage production/removal pathway",
        "sce04146  Peroxisome",
        "sce04213  Longevity regulating pathway - multiple species",
        "sce00480  Glutathione metabolism",
        "sce01100  Metabolic pathways",
        "sce00450  Selenocompound metabolism",
        "sce00380  Tryptophan metabolism",
        "sce00630  Glyoxylate and dicarboxylate metabolism",
        "sce01110  Biosynthesis of secondary metabolites",
        "sce01200  Carbon metabolism",
        "sce04011  MAPK signaling pathway - yeast"]
newPathways = [join(kegg[[1, 2, 3]]," "),
               join(kegg[[1, 2, 3]]," "),
               kegg[1], kegg[1],
               join(kegg[[1, 4, 5]]," "),
               join(kegg[[1, 6]]," "),
               join(kegg[[1, 2, 3, 5, 7, 8, 9, 10, 11]]," "),
               kegg[1], kegg[1], kegg[1], kegg[1],
               join(kegg[[1, 6]]," "),
               join(kegg[[1, 2, 3, 5, 7, 8, 9, 10, 11]]," ")]
enzymePathways = [enzymePathways; newPathways]


################################################################################
@printf("Update right hand side (b) ... \n")
################################################################################
newChanges = zeros(nNewComponents)
changes = [changes; newChanges]


################################################################################
@printf("Update reaction names ... \n")
################################################################################
newReactions = ["oxygen transport",
                "oxygen transport (reversible)",
                "superoxide production (mitochondria)", 
                "pseudo superoxide leakage minimum",
                "pseudo superoxide leakage maximum", 
                "superoxide oxidoreductase (arm, mitochondria)",
                "superoxide oxidoreductase (SOD1, mitochondria)",
                "draw_prot_P00445",
                "superoxide oxidoreductase (SOD2, mitochondria)",
                "draw_prot_P00447",
                "hydroxyl production via peroxynitrite (mitochondria)",
                "glutathione peroxidase (GPX3, mitochondria)",
                "glutathione reductase (GLR1, mitochondria)",
                "thioredoxin peroxidase (TRX3, mitochondria)",
                "thioredoxin reductase (TRR2, mitochondria)",
                "draw_prot_P40581",
                "draw_prot_P41921",
                "draw_prot_P25372",
                "draw_prot_P38816",
                "hydrogen peroxide catalase (CTA1, mitochondria)",
                "draw_prot_P15202",
                "Fenton (mitochondria)",
                "Haber Weiss (mitochondria)",
                "protein damage production (mitochondria)",
                "protein damage exchange (mitochondria)",
                "proton transport",
                "proton transport (reversible)",
                "hydrogen peroxide transport",
                "hydrogen peroxide transport (reversible)",
                "superoxide production (cytosol)", 
                "superoxide oxidoreductase (SOD1, cytosol)",
                "hydroxyl production via peroxynitrite (cytosol)",
                "glutathione peroxidase (arm, cytosol)",
                "glutathione peroxidase (GPX1, cytosol)",
                "glutathione peroxidase (GPX2, cytosol)",
                "glutathione reductase (GLR1, cytosol)",
                "thioredoxin peroxidase (arm, cytosol)",
                "thioredoxin peroxidase (TRX1, cytosol)",
                "thioredoxin peroxidase (TRX2, cytosol)",
                "thioredoxin reductase (TRR1, cytosol)",
                "draw_prot_P36014",
                "draw_prot_P38143",
                "draw_prot_P22217",
                "draw_prot_P22803",
                "draw_prot_P29509",
                "hydrogen peroxide catalase (CTT1, cytosol)",
                "draw_prot_P06115",
                "Fenton (cytosol)",
                "Haber Weiss (cytosol)",
                "protein damage production (cytosol)",
                "protein damage exchange (cytosol)",
                "proton production (cytosol)"]
nNewReactions = length(newReactions)
reactionNames = [reactionNames; newReactions]
newPathways = ["Oxidative stress" for i = 1 : nNewReactions]
newPathways[contains.(newReactions, "draw_prot")] .= "Protein usage"
newPathways[end] = "Other"
reactionPathways = [reactionPathways; newPathways]


################################################################################
@printf("Update lower and upper bounds ... \n")
################################################################################
newLowerBounds = zeros(nNewReactions)
lowerBounds = [lowerBounds; newLowerBounds]

newUpperBounds = ones(nNewReactions) .* 1000
upperBounds = [upperBounds; newUpperBounds]


################################################################################
@printf("Resize stochiometry matrix ... \n")
################################################################################
elements = findnz(stochiometry)
stochiometry = sparse(elements[1], elements[2], elements[3], 
                      nComponents + nNewComponents, 
                      nReactions + nNewReactions)


################################################################################
@printf("Fill stochiometry matrix with new damage related reactions ... \n")
################################################################################
# all reactions that were added are based on the ecYeast model 8 in Lu et al., 
# Nature Communications 2019, and literature as Kanti Das et al., Archives of 
# Neuroscience 2014, Temple et al., Trends In Cell Biology 2005, Zhao et al., 
# International Journal of Molecular Medicine 2019, Cobley, Antioxidants 2020, 
# and many more
#
# in all variable names: C -> in cytosol, M -> in mitochondria
# see all details in supplementary file damageReactions.csv

# new reaction 1+2 #############################################################
# O2 transport from cytosol to mitochondria (and reversible)
oxygenC = findfirst(x -> x == "O2 [cytosol]", componentNames)
oxygenM = findfirst(x -> x == "O2 [mitochondria]", componentNames)
idx = [oxygenC, oxygenM]
coefficients = [-1.0, +1.0]
stochiometry[idx, nReactions + 1] = coefficients
stochiometry[idx, nReactions + 2] = -coefficients

# new reaction 3 ###############################################################
# mitochondrial superoxide production
superoxideM = findfirst(x -> x == "superoxide [mitochondria]", componentNames)
idx = [oxygenM, superoxideM]
coefficients = [-1.0, +1.0]
stochiometry[idx, nReactions + 3] = coefficients

# new reaction 4 ###############################################################
# pseudo constraint: flux through superoxide production [mitochondria] 
# bigger that 0.2% of flux through complex 1 (weight 70%) and 
# complex 3 (weight 30%)
complex1 = findfirst(x -> x == "NADH:ubiquinone oxidoreductase (No1)", 
                     reactionNames)
complex3 = findfirst(x -> x == "ubiquinol:ferricytochrome c reductase (No1)",
                     reactionNames)
component = findfirst(x -> x == "superoxide leakage minimum [mitochondria]",
                      componentNames)
idx = [nReactions + 3, complex1, complex3, nReactions + 4]
coefficients = [-1.0, +0.002 * 0.7, +0.002 * 0.3, +1.0]
stochiometry[component, idx] = coefficients 

# new reaction 5 ###############################################################
# pseudo constraint: flux through superoxide production [mitochondria] 
# smaller that 2% of flux through complex 1 (weight 70%) and 
# complex 3 (weight 30%)
component = findfirst(x -> x == "superoxide leakage maximum [mitochondria]",
                      componentNames)
idx = [nReactions + 3, complex1, complex3, nReactions + 5] 
coefficients = [+1.0, -0.02 * 0.7, -0.02 * 0.3, +1.0]
stochiometry[component, idx] = coefficients

# new reaction 6 ###############################################################
# superoxide dismutase (arm, mitochondria)
protonM = findfirst(x -> x == "h+ [mitochondria]", componentNames)
pseudoSuperoxideM = findfirst(x -> x == "pmet_superoxide [mitochondria]", 
                              componentNames)                            
idx = [superoxideM, protonM, pseudoSuperoxideM]
coefficients = [-2.0, -2.0, +1.0]
stochiometry[idx, nReactions + 6] = coefficients

# new reaction 7 ###############################################################
# superoxide dismutase (SOD1, mitochondria)
sod1 = findfirst(x -> x == "prot_P00445", componentNames)
peroxideM = findfirst(x -> x == "H2O2 [mitochondria]", componentNames)
idx = [pseudoSuperoxideM, sod1, peroxideM, oxygenM]
coefficients = [-1.0, -2.7778e-9, +1.0, +1.0]
stochiometry[idx, nReactions + 7] = coefficients

# new reaction 8 ###############################################################
# draw SOD1 from pool
pool = findfirst(x -> x == "prot_pool", componentNames)
idx = [pool, sod1]
coefficients = [-15.8544, +1.0]
stochiometry[idx, nReactions + 8] = coefficients

# new reaction 9 ###############################################################
# superoxide dismutase (SOD2, mitochondria)
sod2 = findfirst(x -> x == "prot_P00447", componentNames)
idx = [pseudoSuperoxideM, sod2, peroxideM, oxygenM]
coefficients = [-1.0, -2.7778e-9, +1.0, +1.0]
stochiometry[idx, nReactions + 9] = coefficients

# new reaction 10 ##############################################################
# draw SOD2 from pool
idx = [pool, sod2]
coefficients = [-25.7738, +1.0]
stochiometry[idx, nReactions + 10] = coefficients 

# new reaction 11 ##############################################################
# hydroxyl production via peroxynitrite
hydroxylM = findfirst(x -> x == "hydroxyl radical [mitochondria]", componentNames)
idx = [superoxideM, protonM, hydroxylM]
coefficients = [-1.0, -1.0, +1.0]
stochiometry[idx, nReactions + 11] = coefficients

# new reaction 12 ##############################################################
# removal of H2O2 via glutathione peroxidase (GPX3, mitochondria)
glutathioneM = findfirst(x -> x == "glutathione [mitochondria]", componentNames)
glutathioneDidulfideM = findfirst(x -> x == "glutathione disulfide [mitochondria]",
                                  componentNames)
gpx3 = findfirst(x -> x == "prot_P40581", componentNames)
waterM = findfirst(x -> x == "H2O [mitochondria]", componentNames)
idx = [glutathioneM, peroxideM, gpx3, waterM, glutathioneDidulfideM]
coefficients = [-2.0, -1.0, -2.9026e-7, +2.0, +1.0]
stochiometry[idx, nReactions + 12] = coefficients

# new reaction 13 ##############################################################
# back reaction to create glutathione again via glutathione 
# reductase (GLR1, mitochondria)
glr1 = findfirst(x -> x == "prot_P41921", componentNames)
nadphM = findfirst(x -> x == "NADPH [mitochondria]", componentNames)
nadpM = findfirst(x -> x == "nadp+ [mitochondria]", componentNames)
idx = [glutathioneDidulfideM, protonM, nadphM, glr1, glutathioneM, nadpM]
coefficients = [-1.0, -1.0, -1.0, -3.0864e-7, +2.0, +1.0]
stochiometry[idx, nReactions + 13] = coefficients

# new reaction 14 ##############################################################
# removal of H2O2 via thioredoxin peroxidase (TRX3, mitochondria)
thioredoxinM = findfirst(x -> x == "thioredoxin [mitochondria]", componentNames)
thioredoxinDidulfideM = findfirst(x -> x == "thioredoxin disulfide [mitochondria]",
                                  componentNames)
trx3 = findfirst(x -> x == "prot_P25372", componentNames)
idx = [thioredoxinM, peroxideM, trx3, waterM, thioredoxinDidulfideM]
coefficients = [-1.0, -1.0, -2.9026e-7, +2.0, +1.0]
stochiometry[idx, nReactions + 14] = coefficients

# new reaction 15 ##############################################################
# back reaction to create thioredoxin again via thioredoxin 
# reductase (TRR2, mitochondria)
trr2 = findfirst(x -> x == "prot_P38816", componentNames)
idx = [thioredoxinDidulfideM, protonM, nadphM, trr2, thioredoxinM, nadpM]
coefficients = [-1.0, -1.0, -1.0, -8.3417e-6, +1.0, +1.0]
stochiometry[idx, nReactions + 15] = coefficients

# new reaction 16 ##############################################################
# draw GPX3 from pool
idx = [pool, gpx3]
coefficients = [-18.6412, +1.0]
stochiometry[idx, nReactions + 16] = coefficients 

# new reaction 17 ##############################################################
# draw GLR1 from pool
idx = [pool, glr1]
coefficients = [-53.4402, +1.0]
stochiometry[idx, nReactions + 17] = coefficients 

# new reaction 18 ##############################################################
# draw TRX3 from pool
idx = [pool, trx3]
coefficients = [-14.4320, +1.0]
stochiometry[idx, nReactions + 18] = coefficients 

# new reaction 19 ##############################################################
# draw TRR2 from pool
idx = [pool, trr2]
coefficients = [-37.0870, +1.0]
stochiometry[idx, nReactions + 19] = coefficients 

# new reaction 20 ##############################################################
# removal of H2O2 via hydrogen peroxide catalase (CTA1, mitochondria)
cta1 = findfirst(x -> x == "prot_P15202", componentNames)
idx = [peroxideM, cta1, waterM, oxygenM]
coefficients = [-1.0, -3.2680e-10, +1.0, +1.0]
stochiometry[idx, nReactions + 20] = coefficients

# new reaction 21 ##############################################################
# draw CTA1 from pool
idx = [pool, cta1]
coefficients = [-58.5548, +1.0]
stochiometry[idx, nReactions + 21] = coefficients 

# new reaction 22 ##############################################################
# Fenton reaction to make Haber Weiss reaction possible and
iron2M = findfirst(x -> x == "iron(2+) [mitochondria]", componentNames)
iron3M = findfirst(x -> x == "iron(3+) [mitochondria]", componentNames)
idx = [iron3M, superoxideM, iron2M, oxygenM]
coefficients = [-1.0, -1.0, +1.0, +1.0]
stochiometry[idx, nReactions + 22] = coefficients 

# new reaction 23 ##############################################################
# Haber Weiss reaction to create hydroxyl radicals
idx = [iron2M, peroxideM, iron3M, hydroxylM]
coefficients = [-1.0, -1.0, +1.0, +1.0]
stochiometry[idx, nReactions + 23] = coefficients 

# new reaction 24 ##############################################################
# protein damage production 
damageM = findfirst(x -> x == "protein damage [mitochondria]", componentNames)
idx = [hydroxylM, damageM]
coefficients = [-1.0, 1.0]
stochiometry[idx, nReactions + 24] = coefficients

# new reaction 25 ##############################################################
# damage exchange (mitochondria)
idx = [damageM]
coefficients = [-1.0]
stochiometry[idx, nReactions + 25] = coefficients

# new reaction 26+27 ###########################################################
# proton transport and (reversible)
# otherwise infeasible
protonC = findfirst(x -> x == "h+ [cytosol]", componentNames)
idx = [protonC, protonM]
coefficients = [-1.0, +1.0]
stochiometry[idx, nReactions + 26] = coefficients
stochiometry[idx, nReactions + 27] = -coefficients

# new reaction 28+29 ###########################################################
# hydrogen peroxide transport and (reversible)
peroxideC = findfirst(x -> x == "H2O2 [cytosol]", componentNames)
idx = [peroxideC, peroxideM]
coefficients = [-1.0, +1.0]
stochiometry[idx, nReactions + 28] = coefficients
stochiometry[idx, nReactions + 29] = -coefficients

# new reaction 30 ##############################################################
# cytosolic superoxide production 
superoxideC = findfirst(x -> x == "superoxide [cytosol]", componentNames)
idx = [oxygenC, superoxideC]
coefficients = [-1.0, +1.0]
stochiometry[idx, nReactions + 30] = coefficients

# new reaction 31 ##############################################################
# superoxide dismutase (SOD1, cytosol)
idx = [superoxideC, protonC, sod1, peroxideC, oxygenC]
coefficients = [-2.0, -2.0, -2.7778e-9, +1.0, +1.0]
stochiometry[idx, nReactions + 31] = coefficients

# new reaction 32 ##############################################################
# hydroxyl production via peroxynitrite (cytosol)
hydroxylC = findfirst(x -> x == "hydroxyl radical [cytosol]", componentNames)
idx = [superoxideC, protonC, hydroxylC]
coefficients = [-1.0, -1.0, +1.0]
stochiometry[idx, nReactions + 32] = coefficients

# new reaction 33 ##############################################################
# glutathione peroxidase (arm, cytosol)
pseudoGlutathioneC = findfirst(x -> x == "pmet_glutathione [cytosol]", 
                         componentNames)
glutathioneC = findfirst(x -> x == "glutathione [cytosol]", componentNames)
idx = [glutathioneC, peroxideC, pseudoGlutathioneC]
coefficients = [-2.0, -1.0, +1.0]
stochiometry[idx, nReactions + 33] = coefficients

# new reaction 34 ##############################################################
# removal of H2O2 via glutathione peroxidase (GPX1, cytosol)
glutathioneDidulfideC = findfirst(x -> x == "glutathione disulfide [cytosol]",
                                  componentNames)
gpx1 = findfirst(x -> x == "prot_P36014", componentNames)
waterC = findfirst(x -> x == "H2O [cytosol]", componentNames)
idx = [pseudoGlutathioneC, gpx1, waterC, glutathioneDidulfideC]
coefficients = [-1.0, -2.9026e-7, +2.0, +1.0]
stochiometry[idx, nReactions + 34] = coefficients

# new reaction 35 ##############################################################
# removal of H2O2 via glutathione peroxidase (GPX2, cytosol)
gpx2 = findfirst(x -> x == "prot_P38143", componentNames)
idx = [pseudoGlutathioneC, gpx2, waterC, glutathioneDidulfideC]
coefficients = [-1.0, -2.9026e-7, +2.0, +1.0]
stochiometry[idx, nReactions + 35] = coefficients

# new reaction 36 ##############################################################
# back reaction to create glutathione again via glutathione reductase 
# (GLR1, cytosol)
glr1 = findfirst(x -> x == "prot_P41921", componentNames)
nadphC = findfirst(x -> x == "NADPH [cytosol]", componentNames)
nadpC = findfirst(x -> x == "nadp+ [cytosol]", componentNames)
idx = [glutathioneDidulfideC, protonC, nadphC, glr1, glutathioneC, nadpC]
coefficients = [-1.0, -1.0, -1.0, -3.0864e-7, +2.0, +1.0]
stochiometry[idx, nReactions + 36] = coefficients

# new reaction 37 ##############################################################
# thioredoxin peroxidase (arm, cytosol)
pseudoThioredoxinC = findfirst(x -> x == "pmet_thioredoxin [cytosol]", 
                         componentNames)
thioredoxinC = findfirst(x -> x == "thioredoxin [cytosol]", componentNames)
idx = [thioredoxinC, peroxideC, pseudoThioredoxinC]
coefficients = [-1.0, -1.0, +1.0]
stochiometry[idx, nReactions + 37] = coefficients

# new reaction 38 ##############################################################
# removal of H2O2 via thioredoxin peroxidase (TRX1, cytosol)
thioredoxinDidulfideC = findfirst(x -> x == "thioredoxin disulfide [cytosol]",
                                  componentNames)
trx1 = findfirst(x -> x == "prot_P22217", componentNames)
idx = [pseudoThioredoxinC, trx1, waterC, thioredoxinDidulfideC]
coefficients = [-1.0, -2.9026e-7, +2.0, +1.0]
stochiometry[idx, nReactions + 38] = coefficients

# new reaction 39 ##############################################################
# removal of H2O2 via thioredoxin peroxidase (TRX2, cytosol)
trx2 = findfirst(x -> x == "prot_P22803", componentNames)
idx = [pseudoThioredoxinC, trx2, waterC, thioredoxinDidulfideC]
coefficients = [-1.0, -2.9026e-7, +2.0, +1.0]
stochiometry[idx, nReactions + 39] = coefficients

# new reaction 40 ##############################################################
# back reaction to create thioredoxin again via thioredoxin reductase 
# (TRR1, cytosol)
trr1 = findfirst(x -> x == "prot_P29509", componentNames)
idx = [thioredoxinDidulfideC, protonC, nadphC, trr1, thioredoxinC, nadpC]
coefficients = [-1.0, -1.0, -1.0, -8.3417e-6, +1.0, +1.0]
stochiometry[idx, nReactions + 40] = coefficients

# new reaction 41 ##############################################################
# draw gpx1 from pool
idx = [pool, gpx1]
coefficients = [-19.4842, +1.0]
stochiometry[idx, nReactions + 41] = coefficients 

# new reaction 42 ##############################################################
# draw gpx2 from pool
idx = [pool, gpx2]
coefficients = [-18.4058, +1.0]
stochiometry[idx, nReactions + 42] = coefficients 

# new reaction 43 ##############################################################
# draw trx1 from pool
idx = [pool, trx1]
coefficients = [-11.2348, +1.0]
stochiometry[idx, nReactions + 43] = coefficients 

# new reaction 44 ##############################################################
# draw trx2 from pool
idx = [pool, trx2]
coefficients = [-11.2038, +1.0]
stochiometry[idx, nReactions + 44] = coefficients 

# new reaction 45 ##############################################################
# draw trr1 from pool
idx = [pool, trr1]
coefficients = [-34.2377, +1.0]
stochiometry[idx, nReactions + 45] = coefficients 

# new reaction 46 ##############################################################
# removal of H2O2 via hydrogen peroxide catalase (CTT1, mitochondria)
ctt1 = findfirst(x -> x == "prot_P06115", componentNames)
idx = [peroxideC, ctt1, waterC, oxygenC]
coefficients = [-1.0, -3.2680e-10, +1.0, +1.0]
stochiometry[idx, nReactions + 46] = coefficients

# new reaction 47 ##############################################################
# draw ctt1 from pool
idx = [pool, ctt1]
coefficients = [-64.5825, +1.0]
stochiometry[idx, nReactions + 47] = coefficients 

# new reaction 48 ##############################################################
# Fenton (cytosol) and (reversible)
iron2C = findfirst(x -> x == "iron(2+) [cytosol]", componentNames)
iron3C = findfirst(x -> x == "iron(3+) [cytosol]", componentNames)
idx = [iron3C, superoxideC, iron2C, oxygenC]
coefficients = [-1.0, -1.0, +1.0, +1.0]
stochiometry[idx, nReactions + 48] = coefficients 

# new reaction 49 ##############################################################
# Haber Weiss (cytosol) and (reversible)
idx = [iron2C, peroxideC, iron3C, hydroxylC]
coefficients = [-1.0, -1.0, +1.0, +1.0]
stochiometry[idx, nReactions + 49] = coefficients 

# new reaction 50 ##############################################################
# protein damage production (cytosol)
damageC = findfirst(x -> x == "protein damage [cytosol]", componentNames)
idx = [hydroxylC, damageC]
coefficients = [-1.0, 1.0]
stochiometry[idx, nReactions + 50] = coefficients

# new reaction 51 ##############################################################
# damage exchange (cytosol)
idx = [damageC]
coefficients = [-1.0]
stochiometry[idx, nReactions + 51] = coefficients

# new reaction 52 ##############################################################
# proton production (cytosol)
idx = [protonC]
coefficients = [+1.0]
stochiometry[idx, nReactions + 52] = coefficients


################################################################################
@printf("Add reaction for non-growth associated maintenance ... \n")
################################################################################
# from ecYeast model v8 in Lu et al., Nature Communications 2019
reactionNames = [reactionNames; "NGAM"]

# construct the new reaction
actors = ["ADP [cytosol]", "ATP [cytosol]", "h+ [cytosol]", 
          "phosphate [cytosol]", "H2O [cytosol]"] 
idx = map(y -> findfirst(x -> x == y, componentNames), actors)
newReaction = zeros(Float64, nComponents + nNewComponents)
newReaction[idx] = [1.0, -1.0, 1.0, 1.0, -1.0]

# update the model
stochiometry = [stochiometry newReaction]
lowerBounds = [lowerBounds; 0.0]
upperBounds = [upperBounds; 1000.0] 
reactionPathways = [reactionPathways; "Other"]


################################################################################
@printf("Order ... \n")
################################################################################
# get the right order
enzymeIdx = contains.(reactionNames, "draw_prot")
poolIdx = contains.(reactionNames, "prot_pool")
fluxIdx = (enzymeIdx + poolIdx) .< 1

order = [findall(x -> x == 1, fluxIdx); findall(x -> x == 1, enzymeIdx);
         findall(x -> x == 1, poolIdx)]
         
# update orders
reactionNames = reactionNames[order]
reactionPathways = reactionPathways[order]
stochiometry = stochiometry[:, order]
lowerBounds = lowerBounds[order]
upperBounds = upperBounds[order]

# make sure that the enzymes names are updated accordingly again
# to have them in the same order as the enymes occur in the reactionNames
# make sure that the order in the enzymes names is the same as the order of the
# enzyme reactions in reactionNames
nFluxes = sum(fluxIdx)
namesInReactions = map(x -> x[11:end], reactionNames[nFluxes+1:end-1])
orderInEnzymes = map(x -> findfirst(y -> y == x, enzymes[:, 1]), 
                     namesInReactions)

# update order
enzymes = enzymes[orderInEnzymes, :]
enzymePathways = enzymePathways[orderInEnzymes]


################################################################################
@printf("Save new MAT model ... \n")
################################################################################
newModel = Dict("metNames" => componentNames, "rxnNames" => reactionNames,
                "S" => stochiometry, "b" => changes, "lb" => lowerBounds,
                "ub" => upperBounds, "enzymes" => enzymes, 
                "enzymePathways" => enzymePathways, 
                "rxnPathways" => reactionPathways)
newFile = Dict("reducedEcYeastWithDamage" => newModel)
matwrite("reducedEcYeast_modified_withDamage.mat", newFile)


