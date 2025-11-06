# Modelling *Vibrio natriegens*

This repository contains genome-scale metabolic models (GSMMs) and resource balance analysis (RBA) models for *Vibrio natriegens*, the fastest-growing organism known.

## Table of Contents
- [Overview](#overview)
- [GSMM (Genome-Scale Metabolic Model)](#gsmm-genome-scale-metabolic-model)
  - [Model Files](#model-files)
  - [Installation](#installation)
  - [Loading and Using the Model](#loading-and-using-the-model)
  - [Model Modifications in v1.1](#model-modifications-in-v11)
  - [Editing the Model](#editing-the-model)
- [RBA (Resource Balance Analysis)](#rba-resource-balance-analysis)
- [Citation](#citation)

---

## Overview

*Vibrio natriegens* is a fast-growing marine bacterium with significant potential for biotechnology applications. This repository provides:

1. **iLC858** - A genome-scale metabolic model with 858 genes
2. **iLC858_v1.1** - An updated version with improved annotations and corrected metabolic pathways
3. **RBA models** - Resource balance analysis models for more detailed constraint-based modeling

---

## GSMM (Genome-Scale Metabolic Model)

### Model Files

The GSMM models are located in the `GSMM/` directory:

- **`iLC858_v1.1.sbml`** - **Recommended** - Latest version with improved annotations and pathway corrections
- **`iLC858.sbml`** - Original version (legacy)
- **`updates/update_v1.1.py`** - Python script documenting all modifications made to create v1.1

### Installation

To work with the GSMM, you'll need Python 3.7+ and the COBRApy library:

```bash
# Create a virtual environment (recommended)
python3 -m venv vnat-env
source vnat-env/bin/activate  # On macOS/Linux
# or
vnat-env\Scripts\activate  # On Windows

# Install COBRApy and dependencies
pip install cobra
pip install python-libsbml  # For SBML support
```

Optional but useful packages:
```bash
pip install matplotlib  # For flux visualization
pip install pandas      # For data analysis
pip install jupyter     # For interactive notebooks
```

### Loading and Using the Model

#### Basic Usage

```python
import cobra

# Load the model
model = cobra.io.read_sbml_model('GSMM/iLC858_v1.1.sbml')

# Set the solver (options: 'glpk', 'cplex', 'gurobi')
model.solver = 'glpk'

# Optimize growth
solution = model.optimize()
print(f"Growth rate: {solution.objective_value:.3f} h⁻¹")

# View model statistics
print(f"Genes: {len(model.genes)}")
print(f"Reactions: {len(model.reactions)}")
print(f"Metabolites: {len(model.metabolites)}")
```

#### Simulate Growth on Different Carbon Sources

```python
# Example: Growth on glucose
model.reactions.EX_cpd00027_e.lower_bound = -10  # Glucose uptake
solution = model.optimize()
print(f"Growth on glucose: {solution.objective_value:.3f} h⁻¹")

# Example: Growth on acetate
model.reactions.EX_cpd00027_e.lower_bound = 0    # Block glucose
model.reactions.EX_cpd00029_e.lower_bound = -10  # Acetate uptake
solution = model.optimize()
print(f"Growth on acetate: {solution.objective_value:.3f} h⁻¹")
```

#### Flux Balance Analysis (FBA)

```python
# Perform FBA and examine fluxes
solution = model.optimize()

# View top active reactions
flux_dict = solution.fluxes.to_dict()
sorted_fluxes = sorted(flux_dict.items(), key=lambda x: abs(x[1]), reverse=True)

print("Top 10 reactions by flux:")
for rxn_id, flux in sorted_fluxes[:10]:
    rxn = model.reactions.get_by_id(rxn_id)
    print(f"{rxn_id}: {flux:.3f} | {rxn.name}")
```

#### Flux Variability Analysis (FVA)

```python
from cobra.flux_analysis import flux_variability_analysis

# Analyze flux variability for all reactions
fva_result = flux_variability_analysis(model, fraction_of_optimum=0.9)
print(fva_result.head())

# For specific reactions
reactions_of_interest = ['rxn00288_c', 'rxn10113_c', 'rxn14426_c']
fva_subset = flux_variability_analysis(model, reaction_list=reactions_of_interest)
print(fva_subset)
```

#### Gene Knockout Analysis

```python
from cobra.manipulation import delete_model_genes

# Single gene knockout
with model:  # Context manager ensures changes are reverted
    gene = model.genes.get_by_id("PN96_09285")
    gene.knock_out()
    solution = model.optimize()
    print(f"Growth after knocking out {gene.id}: {solution.objective_value:.3f}")

# Systematic gene knockout study
from cobra.flux_analysis import single_gene_deletion

deletion_results = single_gene_deletion(model)
essential_genes = deletion_results[deletion_results['growth'] < 0.01]
print(f"Number of essential genes: {len(essential_genes)}")
```

#### Export Results

```python
# Export fluxes to CSV
solution = model.optimize()
solution.fluxes.to_csv('flux_results.csv')

# Export model to JSON (easier to read than SBML)
cobra.io.save_json_model(model, 'iLC858_v1.1.json')
```

### Model Modifications in v1.1

The v1.1 update includes significant improvements to biological accuracy. Key modifications:

#### 1. **Electron Transport Chain Corrections**
- **Cytochrome bo3** (`rxn10113_c`): Corrected proton pumping stoichiometry (4H⁺ instead of 2.5H⁺)
  - Gene: `PN96_21415`, `PN96_21420`, `PN96_21425`, `PN96_21430`
- **Cytochrome bd** (`rxn10806_c`): Updated for low-oxygen conditions
  - Gene: `PN96_08165`, `PN96_08170`
- **Cytochrome c oxidase cbb3** (`rxn14426_c`): Added terminal oxidase
  - Gene: `PN96_05720`, `PN96_05730`, `PN96_05735`, `PN96_05740`
- **Cytochrome c oxidase aa3** (`rxn19357_c`): Added high-oxygen terminal oxidase
  - Gene: `PN96_22635`, `PN96_22640`, `PN96_22625`
- **Ubiquinol-cytochrome c oxidoreductase bc1** (`rxn35348_c`): Added complex III
  - Gene: `PN96_11285`, `PN96_11290`, `PN96_11295`

#### 2. **Quinone Metabolism**
- Corrected succinate dehydrogenase to use ubiquinone under aerobic conditions
- Fixed fumarate reductase to use menaquinone under anaerobic conditions
- Removed thermodynamically incorrect ubiquinone usage under anaerobic conditions
- Added NADPH:quinone oxidoreductase reactions for ubiquinone, menaquinone, and demethylmenaquinone

#### 3. **PHB (Polyhydroxybutyrate) Metabolism**
- Fixed PHB synthesis pathway with correct gene annotations
- Reaction: `rxn01453_c` - (R)-3-Hydroxybutanoyl-CoA:NADP+ oxidoreductase
  - Gene: `PN96_18045`

#### 4. **Biomass Composition**
- Updated quinone requirements to use reduced forms (MKH2 and DMKH2) instead of oxidized forms
- This prevents thermodynamically unfavorable reactions

#### 5. **ATP Maintenance**
- Added ATP maintenance requirement: 3.15 mmol ATP/gDW/h (based on *E. coli* iJO1366)
- Reaction: `rxn00062_c`

#### 6. **Thermodynamic Corrections**
- Made lactate-cytochrome oxidoreductive steps irreversible to prevent thermodynamically impossible ATP generation loops

#### 7. **Iron Transport**
- Updated iron-dicitrate transport mechanism
- Gene: `PN96_16000` (fecB)

#### 8. **Peptidoglycan Synthesis**
- Added murein polymerization for compatibility with RBA modeling
- Reaction: `rxn08957_c`
- Genes: `PN96_00770`, `PN96_02030`, `PN96_02370`, `PN96_10195`

#### 9. **Other Pathway Corrections**
- THF pathway: Gene annotation for `rxn01211_c` → `PN96_09125`
- TCA cycle: Gene annotation for `rxn08094_c` → `PN96_01345`, `PN96_09280`, `PN96_09275`
- ATP synthase: Complete gene annotation for `rxn10042_c`

### Editing the Model

#### Making Modifications to the Model

To create your own model variants:

```python
import cobra
from cobra import Reaction, Metabolite

# Load the model
model = cobra.io.read_sbml_model('GSMM/iLC858_v1.1.sbml')

# Example 1: Add a new reaction
new_rxn = Reaction('my_new_reaction')
new_rxn.name = 'My new enzymatic reaction'
new_rxn.subsystem = 'Custom metabolism'
new_rxn.lower_bound = 0  # Irreversible
new_rxn.upper_bound = 1000

# Add metabolites to the reaction
new_rxn.add_metabolites({
    model.metabolites.get_by_id("cpd00002_c"): -1,  # ATP consumed
    model.metabolites.get_by_id("cpd00008_c"): 1,   # ADP produced
})

# Add gene association
new_rxn.gene_reaction_rule = "PN96_00000 or PN96_00001"

# Add to model
model.add_reactions([new_rxn])

# Example 2: Modify an existing reaction
rxn = model.reactions.get_by_id("rxn00288_c")
rxn.lower_bound = 0  # Make irreversible
rxn.upper_bound = 500  # Reduce maximum flux

# Example 3: Add a new metabolite
new_met = Metabolite(
    'cpd99999_c',
    formula='C6H12O6',
    name='My custom metabolite',
    compartment='c'
)
model.add_metabolites([new_met])

# Example 4: Remove a reaction
model.remove_reactions([model.reactions.get_by_id("rxn27735_c")])

# Validate the model
print(f"Model validates: {model.medium is not None}")

# Test the modified model
solution = model.optimize()
print(f"Growth rate after modifications: {solution.objective_value:.3f}")

# Save the modified model
cobra.io.write_sbml_model(model, 'GSMM/iLC858_custom.sbml')
```

#### Systematic Model Update Script

Use the `update_v1.1.py` script as a template for documenting changes:

```python
import cobra
from cobra import Model, Reaction, Metabolite

# Load base model
model = cobra.io.read_sbml_model('iLC858.sbml')
model.solver = 'glpk'

# Make modifications with clear comments
######################
# Your modification category
######################

# Modify specific reactions
model.reactions.get_by_id("rxn_id").gene_reaction_rule = "gene_id"

# Add new reactions
new_rxn = Reaction("rxn_id")
# ... (reaction details)
model.add_reactions([new_rxn])

# Save updated model
cobra.io.write_sbml_model(model, 'iLC858_v1.2.sbml')
```

#### Best Practices for Model Editing

1. **Always document changes**: Use clear comments explaining the biological rationale
2. **Test after modifications**: Run `model.optimize()` to ensure the model still solves
3. **Validate stoichiometry**: Check mass and charge balance for new reactions
4. **Version control**: Save modified models with new version numbers
5. **Gene annotations**: Use proper locus tags (PN96_xxxxx format for *V. natriegens*)
6. **Literature support**: Document sources for new reactions or modified parameters
7. **Reversibility**: Set appropriate bounds based on thermodynamics
8. **Use context managers**: When testing changes, use `with model:` to avoid permanent modifications

#### Common Model Queries

```python
# Find all reactions associated with a gene
gene = model.genes.get_by_id("PN96_09285")
print(f"Reactions catalyzed by {gene.id}:")
for rxn in gene.reactions:
    print(f"  {rxn.id}: {rxn.name}")

# Find reactions in a specific subsystem
tca_reactions = [r for r in model.reactions if 'TCA' in r.subsystem]

# Find all exchange reactions
exchange_rxns = model.exchanges
print(f"Number of exchange reactions: {len(exchange_rxns)}")

# Check metabolite participation
met = model.metabolites.get_by_id("cpd00002_c")  # ATP
print(f"Reactions using {met.name}:")
for rxn in met.reactions:
    print(f"  {rxn.id}: {rxn.reaction}")

# Get model summary
print(model.summary())

# Get metabolite summary
print(model.metabolites.cpd00067_c.summary())  # Proton summary
```

---

## RBA (Resource Balance Analysis)

The RBA models provide more detailed constraint-based modeling that accounts for:
- Proteome allocation
- Ribosome capacity
- tRNA availability
- Macromolecular composition

RBA models are located in the `RBA/` directory. See RBA-specific documentation for usage.

**Note**: RBA models use the updated iLC858_v1.1.sbml as their metabolic network base.

---

## Citation

If you use these models in your research, please cite the appropriate publications:

[Add relevant publication citations here]

---

## License

[Add license information]

---

## Contact

For questions or issues, please open a GitHub issue or contact the repository maintainers.

---

## Additional Resources

- [COBRApy Documentation](https://cobrapy.readthedocs.io/)
- [SBML Format Specification](http://sbml.org/)
- [BiGG Models Database](http://bigg.ucsd.edu/)
- [ModelSEED](https://modelseed.org/)