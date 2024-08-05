import cobra
from cobra.flux_analysis.loopless import loopless_solution
from cobra import Model, Reaction, Metabolite

# Load
model=cobra.io.read_sbml_model('iLC858.sbml')
model.solver = 'glpk'

######################
# PHB metabolism
######################

model.remove_reactions([model.reactions.get_by_id("rxn27735_c")])

rxn01453 = Reaction("rxn01453_c")
rxn01453.add_metabolites({
    model.metabolites.get_by_id("cpd00006_c"): -1,
    model.metabolites.get_by_id("cpd02234_c"): -1,
    model.metabolites.get_by_id("cpd00005_c"): 1,
    model.metabolites.get_by_id("cpd00067_c"): 1,
    model.metabolites.get_by_id("cpd00279_c"): 1
})
rxn01453.bounds = (-1000, 1000)
rxn01453.name = "(R)-3-Hydroxybutanoyl-CoA:NADP+ oxidoreductase"
rxn01453.gene_reaction_rule = "PN96_18045"
model.add_reactions([rxn01453])

######################
# Anaerobic fumarate reductase (menaquinol electron donor)
######################

model.reactions.get_by_id("rxn08527_c").name = "Fumarate reductase (menaquinone8) (anaerobic)"
model.reactions.get_by_id("rxn08527_c").gene_reaction_rule = "PN96_14370 and PN96_14375 and PN96_14380 and PN96_14385"
model.reactions.get_by_id("rxn08527_c").bounds = (0,1000)

model.reactions.get_by_id("rxn08528_c").name = "Fumarate reductase (2-Demethylmenaquinone8) (anaerobic)"
model.reactions.get_by_id("rxn08528_c").gene_reaction_rule = "PN96_14370 and PN96_14375 and PN96_14380 and PN96_14385"
model.reactions.get_by_id("rxn08528_c").bounds = (0,1000)

# incorrect as ubiquinone can not be used under anaerobic conditions (can not be synthesised)
model.remove_reactions([model.reactions.get_by_id("rxn09272_c")]) 

######################
# Aerobic succinate dehydrogenase (ubiquinone electron acceptor)
######################

# FAD is covalently bound to this enzyme, so this reaction does not need to be split between first succinate -> FADH then FADH -> ubiquinone
model.reactions.get_by_id("rxn00288_c").subtract_metabolites(model.reactions.get_by_id("rxn00288_c").metabolites) # remove current metabolites
model.reactions.get_by_id("rxn00288_c").add_metabolites({
    model.metabolites.get_by_id("cpd00036_c"): -1, # succinate
    model.metabolites.get_by_id("cpd15560_c"): -1, # ubiquinone
    model.metabolites.get_by_id("cpd00106_c"): 1, # fumarate
    model.metabolites.get_by_id("cpd15561_c"): 1 # ubiquinol
})
model.reactions.get_by_id("rxn00288_c").name = "Succinate dehydrogenase (ubiquinone) (aerobic)"
model.reactions.get_by_id("rxn00288_c").gene_reaction_rule = "PN96_09285 and PN96_09290 and PN96_09295 and PN96_09300"
model.reactions.get_by_id("rxn00288_c").bounds = (0,1000)

######################
# Cytochrome bo3
# Mostly active under aerobic conditions -> oxidises ubiquinol with the pumping out of protons
######################

model.reactions.get_by_id("rxn10113_c").gene_reaction_rule = "PN96_21415 and PN96_21420 and PN96_21425 and PN96_21430"
# pumping 4 protons rather than 2.5
model.reactions.get_by_id("rxn10113_c").add_metabolites({
    model.metabolites.get_by_id("cpd00067_c"): -1.5,
    model.metabolites.get_by_id("cpd00067_e"): 1.5
})

######################
# Cytochrome bd
# Mostly active under low aerobic conditions -> oxidises menaquinol. 
# Oxygen reduction costs two cytoplasmatic protons, while menquinol oxidation releases 2 protons in the periplasm, 
# so without pumping still conributes two protons to proton motive force
######################

model.reactions.get_by_id("rxn10806_c").gene_reaction_rule = "PN96_08165 and PN96_08170"

######################
# Cytochrome c oxidase cbb3
######################

rxn14426 = Reaction("rxn14426_c")
rxn14426.add_metabolites({
    model.metabolites.get_by_id("cpd00007_c"): -0.5, # oxygen
    model.metabolites.get_by_id("cpd00067_c"): -4, # protons (2 for H2O, 2 for pumping)
    model.metabolites.get_by_id("cpd00110_c"): -2, # reduced ferrocytochrome
    model.metabolites.get_by_id("cpd00067_e"): 2, # pumped protons
    model.metabolites.get_by_id("cpd00109_c"): 2, # oxidised ferrocytochrome
    model.metabolites.get_by_id("cpd00001_c"): 1 # water
})
rxn14426.bounds = (0, 1000)
rxn14426.name = "cytochrome c oxidase cbb3"
rxn14426.gene_reaction_rule = "PN96_05720 and PN96_05730 and PN96_05735 and PN96_05740"
model.add_reactions([rxn14426])

######################
# Ubiquinol-cytochrome c oxidoreductase bc1
######################

rxn35348 = Reaction("rxn35348_c")
rxn35348.add_metabolites({
    model.metabolites.get_by_id("cpd00067_c"): -2, # intracellular protons
    model.metabolites.get_by_id("cpd00109_c"): -2, # oxidised ferrocytochrome
    model.metabolites.get_by_id("cpd15561_c"): -1, # ubiquinol
    model.metabolites.get_by_id("cpd00067_e"): 4, # extracellular protons
    model.metabolites.get_by_id("cpd00110_c"): 2, # reduced ferrocytochrome
    model.metabolites.get_by_id("cpd15560_c"): 1 # ubiquinone
})
rxn35348.bounds = (0, 1000)
rxn35348.name = "ubiquinol-cytochrome c oxidoreductase cb1 (complex III)"
rxn35348.gene_reaction_rule = "PN96_11285 and PN96_11290 and PN96_11295"
model.add_reactions([rxn35348])

######################
# Cytochrome c oxidase aa3
# This oxidase is a terminal oxidase with the same reaction as cbb3
# The difference in bacteria is that aa3 has lower affinity for oxygen and is therefore more used in high aerobic conditions
######################

rxn19357 = Reaction("rxn19357_c")
rxn19357.add_metabolites({
    model.metabolites.get_by_id("cpd00007_c"): -0.5, # oxygen
    model.metabolites.get_by_id("cpd00067_c"): -4, # protons (2 for H2O, 2 for pumping)
    model.metabolites.get_by_id("cpd00110_c"): -2, # reduced ferrocytochrome
    model.metabolites.get_by_id("cpd00067_e"): 2, # pumped protons
    model.metabolites.get_by_id("cpd00109_c"): 2, # oxidised ferrocytochrome
    model.metabolites.get_by_id("cpd00001_c"): 1 # water
})
rxn19357.bounds = (0, 1000)
rxn19357.name = "cytochrome c oxidase aa3"
rxn19357.gene_reaction_rule = "PN96_22635 and PN96_22640 and PN96_22625"
model.add_reactions([rxn19357])

######################
# NADPH:quinone oxidoreductase
######################

rxn08977 = Reaction("rxn08977_c")
rxn08977.add_metabolites({
    model.metabolites.get_by_id("cpd00005_c"): -1, # NADPH
    model.metabolites.get_by_id("cpd00067_c"): -1, # proton
    model.metabolites.get_by_id("cpd15560_c"): -1, # ubiquinone
    model.metabolites.get_by_id("cpd00006_c"): 1, # NADP
    model.metabolites.get_by_id("cpd15561_c"): 1, # ubiquinol
})
rxn08977.bounds = (0, 1000)
rxn08977.name = "NADPH Quinone Reductase (Ubiquinone-8)"
rxn08977.gene_reaction_rule = "PN96_01220"

rxn08978 = Reaction("rxn08978_c")
rxn08978.add_metabolites({
    model.metabolites.get_by_id("cpd00005_c"): -1, # NADPH
    model.metabolites.get_by_id("cpd00067_c"): -1, # proton
    model.metabolites.get_by_id("cpd15500_c"): -1, # menquinone
    model.metabolites.get_by_id("cpd00006_c"): 1, # NADP
    model.metabolites.get_by_id("cpd15499_c"): 1, # menquinol
})
rxn08978.bounds = (0, 1000)
rxn08978.name = "NADPH Quinone Reductase (Menaquinone-8)"
rxn08978.gene_reaction_rule = "PN96_01220"

rxn08979 = Reaction("rxn08979_c")
rxn08979.add_metabolites({
    model.metabolites.get_by_id("cpd00005_c"): -1, # NADPH
    model.metabolites.get_by_id("cpd00067_c"): -1, # proton
    model.metabolites.get_by_id("cpd15352_c"): -1, # 2-Demethyl menquinone
    model.metabolites.get_by_id("cpd00006_c"): 1, # NADP
    model.metabolites.get_by_id("cpd15353_c"): 1, # 2-Demethyl menquinol
})
rxn08979.bounds = (0, 1000)
rxn08979.name = "NADPH Quinone Reductase (2-Demethylmenaquinone-8)"
rxn08979.gene_reaction_rule = "PN96_01220"

model.add_reactions([rxn08977, rxn08978, rxn08979])

######################
# MKH2 synthesis
######################

model.reactions.get_by_id("rxn08333_c").gene_reaction_rule = "PN96_12040"

# UbiE acts on MKH2 rather MK
# Sources: 10.1021/bi700810x & 10.1128/jb.179.5.1748-1754.1997
model.reactions.get_by_id("rxn10094_c").add_metabolites({
    model.metabolites.get_by_id("cpd15352_c"): 1, # DMK
    model.metabolites.get_by_id("cpd15353_c"): -1, # DMKH2
    model.metabolites.get_by_id("cpd15500_c"): -1, # MK
    model.metabolites.get_by_id("cpd15499_c"): 1, # MKH2
})
model.reactions.get_by_id("rxn10094_c").gene_reaction_rule = "PN96_12880"

######################
# TCA
######################

model.reactions.get_by_id("rxn08094_c").gene_reaction_rule = "PN96_01345 and PN96_09280 and PN96_09275"

######################
# biomass
######################

# Replace MK by MKH2 and DMK by DMKH2 to enable biosynthesis of the reduced quinols rather than the oxidised quinones (which would require thermodynamically unfavourable NAD:quinol oxidoreductive reactions)
model.reactions.get_by_id("bio1").add_metabolites({
    model.metabolites.get_by_id("cpd15499_c"): -0.00309646685192537,
    model.metabolites.get_by_id("cpd15500_c"): 0.00309646685192537,
    model.metabolites.get_by_id("cpd15353_c"): -0.00309646685192537,
    model.metabolites.get_by_id("cpd15352_c"): 0.00309646685192537,
})

######################
# Thermodynamic feasability of lactate - cytochrome oxidoreductive step. 
# This modification is crucial otherwise NADH can be gained from lactate dehydrogenase and then used to start over the electron transport chain (yielding much more than the 2 cytochromes required to reduce pyruvate in rxn00145)
######################

model.reactions.get_by_id("rxn00145_c").bounds = (0,1000)
model.reactions.get_by_id("rxn00146_c").bounds = (0,1000)

######################
# ATP maintenance
# Example bound from E. coli iJO1366: 3.15
#######################

rxn00062 = Reaction("rxn00062_c")
rxn00062.add_metabolites({
    model.metabolites.get_by_id("cpd00001_c"): -1, # H2O
    model.metabolites.get_by_id("cpd00002_c"): -1, # ATP
    model.metabolites.get_by_id("cpd00008_c"): 1, # ADP
    model.metabolites.get_by_id("cpd00009_c"): 1, # phosphate
    model.metabolites.get_by_id("cpd00067_c"): 1, # proton
})
rxn00062.bounds = (3.15,1000)
rxn00062.name = "ATP maintenance"
model.add_reactions([rxn00062])

######################
# THF pathway
#######################

model.reactions.get_by_id("rxn01211_c").gene_reaction_rule = "PN96_09125"

######################
# Iron transport
# Current mechanism: succinate secretion -> succinate - citrate antiport -> dicitrate-ferric iron (Fe3+) import
# PN96_16000 (fecB, iron-dicitrate transporter substrate-binding subunit) was detected in proteomics on glucose aerobic growth even though citrate wasnt added. So current best guess
#######################

model.reactions.get_by_id("rxn05557_c").gene_reaction_rule = "PN96_16000"
model.reactions.get_by_id("rxn00068_c").gene_reaction_rule = "PN96_13805"

######################
# Peptidoglycan synthesis
# This reaction was mostly added for compatability with RBA modelling
# #######################

cpd15511 = Metabolite(
    'cpd15511_c',
    formula='C80H124N16O42',
    name='two linked disacharide pentapeptide murein units (uncrosslinked, middle of chain), 4',
    compartment='c')
rxn08957 = Reaction("rxn08957_c")
rxn08957.add_metabolites({
    model.metabolites.get_by_id("cpd03495_c"): -2, # undecaprenyl diphosphate carrier + murein pentapeptide
    model.metabolites.get_by_id("cpd02229_c"): 2, # undecaprenyl diphosphate carrier
    cpd15511: 1, # two linked murein units
})
rxn08957.bounds = (0, 1000)
rxn08957.name = "murein polymerizing transglycosylase"
rxn08957.gene_reaction_rule = "PN96_00770 or PN96_02030 or PN96_02370 or PN96_10195"
model.add_reactions([rxn08957])


# Write
cobra.io.write_sbml_model(model, 'iLC858_v1.1.sbml')
