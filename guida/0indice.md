# FDG-PET in ALS — Indice della Pipeline

---

## Parte 1 — Preprocessing e Calcolo degli Z-Score ROI

| File | Contenuto | Step |
|------|-----------|------|
| [Preprocessing](1preprocessing) | Selezione ADNI, MATLAB/SPM, NIfTI, smoothing, intensity normalization | 1–4 |
| [ROI & ComBat](parte1_roi_combat.md) | Atlante AAL3, estrazione ROI da uptake, armonizzazione inter-scanner | 5–6 |
| [Normativo & Z-Score](parte1_normativo_zscore.md) | Modello normativo OLS, calcolo z-score ROI, controllo qualità pre-SuStaIn | 7–8 |
| [Appendice](parte1_appendice.md) | Script completi: MATLAB (preprocessing, intensity norm, modello voxel-wise) e Python (percorso visivo) | — |

---

## Parte 2 — SuStaIn

| File | Contenuto |
|------|-----------|
| [SuStaIn](parte2_sustain.md) | Teoria, preparazione dati, fitting, CVIC, PVD, staging |
| [Analisi post-hoc](parte2_analisi_posthoc.md) | Configurazione ALS, correlazioni cliniche, SOD1, EBM |
| [Workflow tecnico](parte2_workflow.md) | Server, tmux, struttura directory, checkpoint pickle |

---

## Flusso dati

```
DICOM (ADNIe/o med nucleare + ALS)
    │
    ▼ Step 1–4 (MATLAB/SPM)
i*_w.nii  (uptake intensità-normalizzato)
    │
    ▼ Step 5 (Python/nibabel)
U_ctrl, U_ALS  (matrici ROI uptake)
    │
    ▼ Step 6 (neuroCombat)
Ũ_ctrl, Ũ_ALS  (matrici armonizzate)
    │
    ├──▶ Step 7 (OLS su Ũ_ctrl)
    │         normative_model_roi.npz
    │
    ▼ Step 8 (Python)
Z  (N_ALS × R, z-score invertiti)
    │
    ▼ Parte 2
SuStaIn → sottotipi + stadi
```

---

## Riferimenti metodologici chiave

1. Young AL et al. (2018). Uncovering the heterogeneity and temporal complexity of neurodegenerative diseases with Subtype and Stage Inference. *Nature Communications*, 9:4273.
2. Cabras S et al. (2025). Role of 2-[18F]FDG-PET as a biomarker of upper motor neuron involvement in amyotrophic lateral sclerosis. *Journal of Neurology*, 272:766.
3. Tan HHG et al. (2022). MRI Clustering Reveals Three ALS Subtypes With Unique Neurodegeneration Patterns. *Annals of Neurology*, 92:1030–1045.
4. Della Rosa et al. (2014). A standardized [18F]-FDG-PET template for spatial normalization in statistical parametric mapping of dementia. *Neuroinformatics*, 12:575–593.
5. Rolls ET et al. (2020). Automated anatomical labelling atlas 3. *NeuroImage*, 206:116189.
6. Johnson WE et al. (2007). Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics*, 8:118–127.
7. Fortin, Jean-Philippe, Drew Parker, Birkan Tunç, et al. «Harmonization of Multi-Site Diffusion Tensor Imaging Data». NeuroImage 161 (novembre 2017): 149–70. https://doi.org/10.1016/j.neuroimage.2017.08.047.
8. Pomponio R et al. (2020). Harmonization of large MRI datasets for the analysis of brain imaging patterns throughout the lifespan. *NeuroImage*, 208:116450.
9. Aksman LM et al. (2021). pySuStaIn: A Python Implementation of the Subtype and Stage Inference Algorithm. *SoftwareX*, 16:100811.
10. https://github.com/PasqualeDellaRosa/Dementia-Specific-18F-FDG-PET-template
11. https://github.com/ucl-pond/pySuStaIn
