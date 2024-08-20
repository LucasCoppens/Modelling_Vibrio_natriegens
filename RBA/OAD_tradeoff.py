from __future__ import absolute_import, division, print_function

import os
import sys
sys.path.append("/Users/lucascoppens/Documents/Phd/Active/Vnat modelling/Vnat_v5/RBA/RBApy")
import rba

# package imports
import re
import cobra
import matplotlib.pyplot as plt 
from cobra.flux_analysis.loopless import loopless_solution
import numpy as np


def main():
    RBA = rba.RbaModel.from_xml("model")
    GEM = cobra.io.read_sbml_model('/Users/lucascoppens/Documents/Phd/Active/Vnat modelling/Vnat_v5/github_repo/GSMM/iLC858_v1.1.sbml')
    GEM.solver = 'glpk'

    # Get FBA WT growth
    GEM.reactions.get_by_id("rxn05488_c").bounds=(-1000, -23.3) # ensure acetate secretion
    GEM_WT_growth = loopless_solution(GEM).objective_value

    # Get FBA Na-OAD KO growth
    GEM.reactions.get_by_id("rxn30509_c").bounds=(0, 0)
    GEM_NaOad_growth = loopless_solution(GEM).objective_value

    # Get RBA WT growth
    RBA_WT_growth = RBA.solve().mu_opt

    # Get RBA Na-OAD KO growth
    RBA.set_enzyme_efficiencies('kapp_NaOAD_KO.tsv')

    RBA_NaOad_growth = RBA.solve().mu_opt

    plot_bar(GEM_WT_growth, GEM_NaOad_growth, RBA_WT_growth, RBA_NaOad_growth)

def plot_bar(GEM_WT_growth, GEM_NaOad_growth, RBA_WT_growth, RBA_NaOad_growth):
    fig, ax = plt.subplots(figsize=(5, 6))

    bars1 = [GEM_WT_growth, RBA_WT_growth]
    bars2 = [GEM_NaOad_growth, RBA_NaOad_growth]

    barWidth = 0.35
    r1 = np.array([0, 1])
    r2 = r1 + barWidth

    plt.bar(r1, bars1, color='#A92549', width=barWidth, label='WT')
    plt.bar(r2, bars2, color='#224F70', width=barWidth, label='$Na^{+}$-OAD KO')

    plt.xticks([r + barWidth/2 for r in range(2)], ['WT', '$Na^{+}$-OAD KO'])  # Fixed line
    plt.ylabel('Growth rate ($h^{-1}$)', fontsize=16.5)
    plt.ylim(0, 2.1)

    plt.yticks(fontsize=16.5)
    plt.xticks(fontsize=16.5)

    plt.tight_layout()
    plt.legend(fontsize = 14.5)
    plt.savefig("NaOAD_tradeoff.png")
    plt.show()

if __name__ == '__main__':
    main()