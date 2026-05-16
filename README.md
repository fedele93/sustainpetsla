# FDG-PET in ALS: Subtypes and Staging with SuStaIn

Pipeline metodologica per l'applicazione dell'algoritmo SuStaIn (*Subtype and Stage Inference*) a dati FDG-PET in pazienti con Sclerosi Laterale Amiotrofica (ALS), con l'obiettivo di identificare sottotipi metabolici e i loro correlati clinici, genetici e prognostici.

---

## Obiettivo

La pipeline implementa un approccio a due stadi:

1. **SuStaIn** sul dataset ALS completo e sui sottogruppi genetici (C9orf72, SOD1) per identificare pattern di progressione metabolica distinti.
2. **EBM** (*Event-Based Model*, `kde_ebm`) su coorti geneticamente omogenee per costruire modelli di staging applicabili a soggetti presintomatici nei trial clinici.

---

## Dataset

| Gruppo | N | Fonte | Ruolo |
|--------|---|-------|-------|
| Controlli CN | 231 | ADNI (*Co-registered, Averaged*) | Popolazione normativa |
| Pazienti ALS | in corso | Centri italiani di neurologia | Target dell'analisi |
| — di cui C9orf72 | ~70 | — | Analisi isolata (Cmax = 3) |
| — di cui SOD1 | ~20 | — | Controllo biologico positivo |

---

## Struttura della pipeline

```
DICOM (ADNI CN + pazienti ALS)
    │
    ▼  Parte 1 — Step 1–4  (MATLAB/SPM25)
i*_w.nii  [uptake FDG intensità-normalizzato, spazio MNI 2mm]
    │
    ▼  Step 5  (Python / nibabel)
U_ctrl, U_ALS  [matrici ROI uptake — atlante AAL3, 24 regioni]
    │
    ▼  Step 6  (neuroCombat)
Ũ_ctrl, Ũ_ALS  [matrici armonizzate inter-scanner]
    │
    ├──▶  Step 7  (OLS su Ũ_ctrl)
    │    normative_model_roi.npz  [β, σ per ROI]
    │
    ▼  Step 8  (Python)
Z  [N_ALS × 24, z-score invertiti — positivo = ipometabolismo]
    │
    ▼  Parte 2
SuStaIn  →  sottotipi + stadi
    │
    ▼
EBM (kde_ebm)  →  staging mutation-specific
```

**Percorso visivo secondario** (non alimenta SuStaIn): il modello normativo voxel-wise e le mappe z-score `z_*.nii` sono prodotti opzionalmente per visualizzazione in MRIcroGL e analisi SPM second-level. Descritti in [Appendice](guida/4appendicepart1.md).

---

## Documentazione — Wiki

| File | Contenuto |
|------|-----------|
| [Indice](guida/0indice.md) | Mappa completa della pipeline e riferimenti |
| **Parte 1** | |
| [Preprocessing](guida/1parte1preprocessing.md) | Step 1–4: ADNI, SPM25, normalizzazione spaziale, smoothing, intensity normalization |
| [ROI & ComBat](guida/2parte1roicombat.md) | Step 5–6: atlante AAL3, estrazione uptake ROI, armonizzazione inter-scanner |
| [Normativo & Z-Score](guida/3parte1normativozscore.md) | Step 7–8: modello OLS, z-score ROI, controllo qualità pre-SuStaIn |
| [Appendice](guida/4appendicepart1.md) | Script completi MATLAB e Python (percorso visivo) |
| **Parte 2** | |
| [SuStaIn](guida/5part2sustain.md) | Teoria, fitting, CVIC, PVD, staging |
| [Analisi post-hoc](guida/6part2posthoc.md) | Correlazioni cliniche, validazione SOD1, EBM |
| [Workflow tecnico](guida/7tecnico.md) | Server, tmux, struttura directory, checkpoint pickle |

---

## Dipendenze

### MATLAB
- MATLAB R2025b
- SPM25
- Template FDG-PET: [Della Rosa et al. 2014](https://github.com/PasqualeDellaRosa/Dementia-Specific-18F-FDG-PET-template) (`TEMPLATE_FDGPET_100.nii`)

### Python
```
nibabel
numpy
pandas
scipy
neuroCombat          # pip install neuroCombat
pySuStaIn            # https://github.com/ucl-pond/pySuStaIn
kde_ebm
matplotlib
scikit-learn
```

Ambiente consigliato: `conda activate sustain_tutorial_env`

### Atlante
- AAL3 v1 (`ROI_MNI_V7.nii`) — Rolls et al. 2020, *NeuroImage* 206:116189

---

## Infrastruttura computazionale

| Macchina | Ruolo | Accesso |
|----------|-------|---------|
| Fedora (laptop) | Sviluppo, preprocessing MATLAB | locale |
| GEEKOM Ubuntu 24.04 (`geekfed`) | Calcoli SuStaIn pesanti | `ssh -L 8888:localhost:8888 fidel@geekfed` |

Sincronizzazione file: `rsync` via Tailscale (`100.126.238.112`).  
Sessioni lunghe: `tmux new -s sustain_als`.

---

## Stato del progetto

- [x] Selezione e preprocessing controlli ADNI (Step 1–4)
- [x] Modello normativo voxel-wise (percorso visivo — Appendice)
- [x] Estrazione ROI da uptake intensità-normalizzato (Step 5)
- [ ] Armonizzazione ComBat — in attesa dei dati scanner/sito ALS (Step 6)
- [ ] Modello normativo ROI + z-score ROI (Step 7–8)
- [ ] SuStaIn dataset completo (Analisi 2a)
- [ ] SuStaIn C9orf72 isolato (Analisi 2b)
- [ ] Analisi post-hoc e validazione SOD1
- [ ] EBM mutation-specific

---

## Riferimenti

1. Young AL et al. (2018). Uncovering the heterogeneity and temporal complexity of neurodegenerative diseases with Subtype and Stage Inference. *Nature Communications*, 9:4273.
2. Cabras S et al. (2025). Role of 2-[18F]FDG-PET as a biomarker of upper motor neuron involvement in amyotrophic lateral sclerosis. *Journal of Neurology*, 272:766.
3. Tan HHG et al. (2022). MRI Clustering Reveals Three ALS Subtypes With Unique Neurodegeneration Patterns. *Annals of Neurology*, 92:1030–1045.
4. Della Rosa et al. (2014). A standardized [18F]-FDG-PET template for spatial normalization in statistical parametric mapping of dementia. *Neuroinformatics*, 12:575–593.
5. Rolls ET et al. (2020). Automated anatomical labelling atlas 3. *NeuroImage*, 206:116189.
6. Johnson WE et al. (2007). Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics*, 8:118–127.
7. Fortin JP et al. (2018). Harmonization of multi-site diffusion tensor imaging data. *NeuroImage*, 161:149–170.
8. Pomponio R et al. (2020). Harmonization of large MRI datasets for the analysis of brain imaging patterns throughout the lifespan. *NeuroImage*, 208:116450.
9. Aksman LM et al. (2021). pySuStaIn: A Python Implementation of the Subtype and Stage Inference Algorithm. *SoftwareX*, 16:100811.

---
 
*Progetto di tesi di specializzazione di neurologia*
