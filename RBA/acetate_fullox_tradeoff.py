from __future__ import absolute_import, division, print_function

import os
import sys
sys.path.append("/Users/lucascoppens/Documents/Phd/Active/Vnat modelling/Vnat_v5/RBA/RBApy")
import rba

# package imports
import re
import cobra
import matplotlib.pyplot as plt 

tca_resp_reactions = [
        "rxn00256_c",
        "rxn00973_c",
        "rxn00198_c",
        "rxn08094_c",
        "rxn00285_c",
        "rxn00288_c",
        "rxn00799_c",
        "rxn00935_c",
        "rxn10042_c", # ATP synthase
        # "rxn30509_c", # Na-OAD
        "rxn37569_c", # Na-NQR
        "rxn10113_c", # cytochrome bo3
        "rxn10806_c", # cytochrome bd
        "rxn14426_c", # cytochrome cbb3
        "rxn35348_c", # cytochrome bc1
        "rxn19357_c", # cytochrome aa3
    ]


def main():
    model = rba.RbaModel.from_xml("model")

    efficiencies = list(range(100000, 1000000, 100000))
    acetate_secretions = []
    ATP_synthases = []
    growth_rates = []
    for efficiency in efficiencies:
        for enz in model.enzymes.enzymes:
            if enz.id[2:12] in tca_resp_reactions:
                model.parameters.functions.get_by_id(enz.id+"_forward_efficiency").parameters.get_by_id('CONSTANT').value = efficiency
                model.parameters.functions.get_by_id(enz.id+"_backward_efficiency").parameters.get_by_id('CONSTANT').value = efficiency

        sol = model.solve(bissection_tol = 0.001)
        rf = sol.reaction_fluxes()
        acetate_secretions.append(abs(rf["R_rxn05488_c"]))
        ATP_synthases.append(abs(rf["R_rxn10042_c"]))
        growth_rates.append(sol.mu_opt)
    
    plot_results(efficiencies, acetate_secretions, ATP_synthases, growth_rates)


def plot_results(efficiencies, acetate_secretions, ATP_synthases, growth_rates):
    
    fig, ax1 = plt.subplots(figsize=(10, 6))
    legend_labels = {}

    color = '#A92549'
    ax1.set_xlabel('TCA & respiratory complexes catalytic rate (*100000$h^{-1}$)', fontsize=16.5)
    ax1.set_ylabel('Growth rate ($h^{-1}$)', fontsize=16.5)
    line1, = ax1.plot(efficiencies, growth_rates, color=color, linewidth=2, label = "Growth rate")
    ax1.set_xticklabels([f"{x/100000:.1f}" for x in efficiencies])

    ax2 = ax1.twinx()

    color = '#224F70'
    line2, = ax2.plot(efficiencies, acetate_secretions, color=color, linewidth=2, label = "Acetate secretion")

    color = '#FFA500'
    ax2.set_ylabel('Acetate secretion / ATP synthase flux\n($mmol * gDW^{-1} * h^{-1}$)', fontsize=16.5)
    line3, = ax2.plot(efficiencies, ATP_synthases, color=color, linewidth=2, label = "ATP synthase")

    plt.yticks(fontsize=16.5)
    plt.xticks(fontsize=16.5)

    ax1.tick_params(axis='x', labelsize=16.5)
    ax1.tick_params(axis='y', labelsize=16.5)
    ax2.tick_params(axis='x', labelsize=16.5)
    ax2.tick_params(axis='y', labelsize=16.5)

    fig.tight_layout()

    lines = [line1, line2, line3]
    labels = [l.get_label() for l in lines]

    plt.legend(lines, labels, loc='upper left', fontsize=14.5)
    plt.savefig("acetate_fullox_tradeoff.png")
    plt.show()

if __name__ == '__main__':
    main()
