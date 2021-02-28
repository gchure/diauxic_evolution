---
class: dry_run
status: accepted
keywords: [growth_curve, RNA_extraction]
person: gchure 
notes: >
    This was a pilot experiment to see how well I could extract total RNA from 
    our wild-type strain grown in a permissive growth condition. Status is 
    'rejected' as it is a dry run experiment and produces no useable data 
---

# 2021-02-25 Total RNA Extraction Dry Run

## Purpose
This experiment was a dry-run pilot experiment to test our extraction method
of total RNA from the wild-type strain NCM3722 in permissive glucose growth
medium. To assess the quality of the resulting total RNA, we looked at the
fragment resolution on an Agligent Bioanalyzer 2100.

## Results 
### Summary 
This was a successful pilot experiment demonstrating that I can successfully
extract total RNA from WT lab-strain *E. coli* with high integrity. 

The observed growth rate was λ ≈ 0.90 / hr, notably slower than the
"standard" growth rate of λ ≈ 0.94 / hr. This may be due to culturing 10 mL
in 30 mm diameter test tubes. Due to space constraints, the cultures were
grown in a "vertical" position rather than tilted as typically done for
smaller tube sizes. Future experiments should see if angled positioning of
the tubes is sufficient to make the growth rates more reproducible.

### Plots

### Growth Curve 
![](./output/2021-02-25_NCM3722_growth_statistics.png)

### Bioanalyzer Calibration
![](./output/2021-02-25_Bioanalyzer_calibration.png)

### Bioanalyzer Results
![](./output/2021-02-25_Bioanalyzer_NCM3722_total_RNA..png)


## Materials 
### Strains
In this experiment, we considered only one strain, the WT NC3722 *E. coli* used 
in myriad other experiments.

| **Strain** | **Genotype** | **Plasmid** | **Lab Strain ID** | **Location** |
| :--: | :--: | :--: | :--: | :--: |
| NCM3722 *E. coli* | -- | -- | `GE046` | `GE Box 1, Pos 46` |

### Growth Medium
This experiment was conducted in our standard `N-C-` Base Buffer supplemented
with Glucose, Nitrogen, and and Micronutrients. See `/miscellaneous/media_recipes.md` 
for information about the `N-C-` medium composition.

|**Chemical Species** | **Function** | **Concentration** |
|:--|:--|:--|
|N-C- Buffer Base | pH Buffering and osmolarity | 1X |
|Micronutrients | Trace metals to support growth | 1X |
| Glucose | Carbon Source | 10 mM |
| Ammonium Chloride | Nitrogen Source | 10 mM|

### RNA Extraction Reagents
The Following reagents were used in the extraction of total RNA in preparation 
for analysis using the Agligent Bioanalyzer.  See `/miscellaneous/reagents_buggeres.md`
for more information 

| **Reagent** | **Purpose** | **Concentration or Amount** |
|:--: | :--: | :--: |
| TMN Cell Resuspension Buffer | RNA Extraction & Stabilization | 1X |
|TRI Reagent (Phenol + Guanidinium isothiocyanate) | RNA Extraction & Stabilization | 250 µL |
| Zymo RNA Clean & Concentrator Kit | RNA Extraction and Stabilization | -- |
| Chloroform | RNA Extraction and Stabilization |  25 µL |
| 200 Proof Ethanol | RNA Extraction and Stabilization | 500 µL |
| RNA-ase free water | RNA Extraction and Stabilization | 15 µL |

### Major Equipment 

| **Instrument** | **Location** | **Step** |
|:--: |:--: |:--: |
| Genesys 30 Spectrophotometer (Thermo Scientific) | Clark W303 | OD<sub>600nm</sub> Measurement |
| Nanodrop ND-1000 Spectrophotometer (Thermo Scientific) | Clark W303 | RNA Concentration Determination |
| Agligent Bioanalyzer 2100 | PAN | RNA integrity determination |

## Experimental Protocol
Bolded regions represent values for this particular experiment while non-bold
text corresponds to the standard protocol used for all similar experiments.

### Bacterial Growth

1. The previous day, a seed culture was prepared by picking 1 colony from an
LB Miller agarose plate and resuspending in 3 mL of LB Miller liquid medium.
The seed culture was allowed to grow for 5 to 6 hours at 37° C with shaking
until saturated.
2. From the saturated culture, a dilution series was performed in the
experimental growth medium to use as precultures for the actual experiment. A
dilution series of 1:10<sup>3</sup>, 1:10<sup>4</sup> and 1:10<sup>5</sup>
was prepared via serial dilution. Precultures were placed in a 37° C water
bath with aeration at 250 rpm at **≈ 10:30 pm** on 2021-02-24.
3. Around 11.5 hours later (**at ≈ 9:00 AM** on 2021-02-25), the
**1:10<sup>4</sup>** preculture reached an OD<sub>600nm</sub> of **≈ 0.4**.
This preculture was diluted **≈ 1:30** into a **10 mL** of growth medium pre
warmed to 37° C. The resulting OD<sub>600nm</sub> of this culture was
**≈0.015**.
4. After dilution, measurements of optical density at 600 nm were made every
10 to 20 minutes. The density was recorded by hand into a paper notebook and
later transferred to a CSV file.
5. Once the OD<sub>600nm</sub> was **≈ 0.365**, the culture was harvested and
total RNA was isolated as described below.
6. Growth curve data was analyzed (see `growth_processing.py`) and a growth
rate was estimated via MCMC implemented by the stan code
(`code/stan/growth_rate_glm.stan`).

### Total RNA Extraction

1. Once cells were at an OD<sub>600nm</sub> of **≈0.365**, the cell volume
was transferred to a 15 mL Falcon tube. Cells were then pelleted by
centrifuging at 5000xg for 10 min at 4° C.
2. Supernatant was discarded and cells were resuspended in 250 µL of 1x TMN
Buffer and was transferred to a sterile 1.5 mL Eppendorf tube.
3. Resuspended pellet was mixed with 250 µL of Trizol and vortexed for 30
seconds. Solution was left to incubate at room temperature for ≈ 5 minutes.
4. 25 µL of chloroform was added to the solution and vortexed for 30 seconds
until thoroughly mixed. Solution was allowed to rest at room temperature for
another 5 minutes.
5. Solution was centrifuged for 2 min at 16000 x g at 4° C.
6. Clear phase of the supernatant (≈ 500 µL) was carefully removed and
transferred to a clean 1.5 mL Eppendorf tube.
7. An equal volume (≈ 500 µL) of 200 proof ethanol was mixed with the
supernatant.
8. Mixture was transferred to a Zymo QI RNA clean & concentrator column.
Column was spun at 16000xg for 1 minute and flow through was discarded.
9. 400 µL of Zymo RNA prep buffer was added to the column and was centrifuged
for 1 min at 16000xg. Flow through was discarded
10. 700 µL of Zymo RNA wash buffer was added to the column and wad
centrifuged for 1 min at 16000 xg. Flow through was discarded.
11. 400 µL of Zymo RNA wash buffer was added to the column and was
centrifuged for 2 min at 16000 xg. Flow through was discarded.
12. 15 µL of RNAase free water (**here I used filter-sterilized and
autoclaved water in lieu of commercial RNA-ase free water) was added to the
column. Column was transferred to a clean 1.5 mL Eppendorf tube and was
centrifuged for 1 min at 16000 xg. Column was discarded.
13. The concentration of the extracted RNA was determined using a Nanodrop.
Concentration was determined to be **≈ 263 ng/µL**. The absorbance ratios were
(260 nm / 280 nm) **≈ 2.08** and (260 nm/ 230 nm) **≈ 2.26**. 
14. A 3 µL aliquot of the extracted RNA was transferred to clean and labeled
PCR tube and was submitted for analysis on an Agligent Bioanalyzer 2100
through the [Protein and Nucleic Acid (PAN) Research
Center](http://panserver11.stanford.edu/).
