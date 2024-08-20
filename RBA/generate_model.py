from __future__ import absolute_import, division, print_function

import os
import sys
sys.path.append("/Users/lucascoppens/Documents/Phd/Active/Vnat modelling/Vnat_v5/RBA/RBApy")
import rba

# package imports
import re
import cobra

def main():
    import_gem()

    # Inital run of model generation creates helper files
    vnat_rba = rba.RbaModel.from_data('params.in')

    # Set a growth medium
    vnat_rba.set_medium('data/medium.tsv')

    update_processes(vnat_rba)

    set_efficiencies(vnat_rba)

    set_maintenance_reaction(vnat_rba)

    set_protein_params(vnat_rba)

    vnat_rba.write()
    
# This GEM modification ensures protons and sodiums are not perceived as medium components
# Otherwise it would complicate the way RBA handles exchange reactions
def import_gem():
    gem = cobra.io.read_sbml_model('/Users/lucascoppens/Documents/Phd/Active/Vnat modelling/Vnat_v5/github_repo/GSMM/iLC858_v1.1.sbml')
    cobra.io.write_sbml_model(gem, 'data/sbml.sbml')

# customize some process parameters
def update_processes(model):
    # set dna concentration to be equal to E. coli RBA model
    fn = model.parameters.functions.get_by_id('dna_concentration')
    fn.parameters.get_by_id('CONSTANT').value = 0.0807
    
    # Remove test processes from default model that are not neccessary
    for pr in ['test_process_0', 'test_process_1', 'test_process_2']:
        pr = model.processes.processes.get_by_id(pr)
        model.processes.processes.remove(pr)

# set k_app default efficiencies analogous to E.coli
def set_efficiencies(model):
    # set default kapp to equal the E coli RBA model
    fn = model.parameters.functions.get_by_id('default_efficiency')
    fn.parameters.get_by_id('CONSTANT').value = 45000
    
    # set kapp respiratory efficiencies
    model.set_enzyme_efficiencies('data/kapp_resp.tsv')

# implement correct ATP maintenance reaction
def set_maintenance_reaction(model):
    model.metabolism.reactions.get_by_id('R_maintenance_atp').reactants.append(rba.xml.SpeciesReference('M_cpd00001_c', 1))
    model.metabolism.reactions.get_by_id('R_maintenance_atp').reactants.append(rba.xml.SpeciesReference('M_cpd00002_c', 1))
    model.metabolism.reactions.get_by_id('R_maintenance_atp').products.append(rba.xml.SpeciesReference('M_cpd00008_c', 1))
    model.metabolism.reactions.get_by_id('R_maintenance_atp').products.append(rba.xml.SpeciesReference('M_cpd00009_c', 1))
    model.metabolism.reactions.get_by_id('R_maintenance_atp').products.append(rba.xml.SpeciesReference('M_cpd00067_c', 1))

# set protein constraints
def set_protein_params(model):
    # Long et al (2017) found an amino acid concentration of 0.47 g/gDW.
    # Assuming an average molecular weight per amino acid of of 110 mg / mmol, 
    # 0.47g/gDCW would translate to 4.2727 mmol/gDCW.
    fn = model.parameters.functions.get_by_id('amino_acid_concentration')
    fn.parameters.get_by_id('LINEAR_CONSTANT').value = 4.2727
    fn.parameters.get_by_id('LINEAR_COEF').value = 0
    fn.parameters.get_by_id('X_MIN').value = 0
    fn.parameters.get_by_id('X_MAX').value = 1
    
    # set amount of secreted protein to 0
    fn = model.parameters.functions.get_by_id('fraction_protein_Secreted')
    fn.parameters.get_by_id('CONSTANT').value = 0

if __name__ == "__main__":
    main()
