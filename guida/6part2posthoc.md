# Parte 2 — Analisi Post-Hoc e Workflow Tecnico

Questo file raccoglie due sezioni distinte:
- [Analisi post-hoc](#analisi-post-hoc) — correlazioni cliniche, SOD1, EBM
- [Workflow tecnico](#workflow-tecnico) — server, tmux, struttura directory

---

## Analisi Post-Hoc

**Input**: `ml_subtype`, `ml_stage` (output SuStaIn) + dati clinici (ALSFRS-R, NfL, ECAS, mutazione genetica)

---

### Configurazione specifica per il dataset ALS

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import pySuStaIn
from scipy.stats import spearmanr, mannwhitneyu

# Carica i dati clinici
df_clinical = pd.read_csv('clinical_data.csv')
# Colonne attese: SubjectID, alsfrs_total, alsfrs_velocity, onset_site,
#                 mutation, nfl_plasma, ecas_total, ptau217, disease_duration

# Aggiungi subtype e stage al dataframe clinico
df_clinical['sustain_subtype'] = ml_subtype
df_clinical['sustain_stage']   = ml_stage
```

### Analisi 2a — Dataset completo (Cmax = 5)

```python
mask_all = np.ones(len(df_clinical), dtype=bool)   # tutti i pazienti ALS

sustain_full = pySuStaIn.ZscoreSustain(
    data               = zdata[mask_all],
    Z_vals             = Z_vals,
    Z_max              = Z_max,
    biomarker_labels   = roi_labels,
    N_startpoints      = 25,
    N_S_max            = 5,
    N_iterations_MCMC  = int(1e6),
    output_folder      = 'output_full',
    dataset_name       = 'ALS_full'
)
```

### Analisi 2b — Sottogruppo C9orf72 (Cmax = 3)

```python
mask_c9 = (df_clinical['mutation'] == 'C9orf72').values

sustain_c9 = pySuStaIn.ZscoreSustain(
    data               = zdata[mask_c9],
    Z_vals             = Z_vals,
    Z_max              = Z_max,
    biomarker_labels   = roi_labels,
    N_startpoints      = 25,
    N_S_max            = 3,   # Cmax ridotto per n~70
    N_iterations_MCMC  = int(1e6),
    output_folder      = 'output_c9',
    dataset_name       = 'C9orf72'
)
```

---

### Correlazioni cliniche

```python
# Sede di insorgenza per sottotipo
print(pd.crosstab(df_clinical['sustain_subtype'],
                  df_clinical['onset_site'],
                  normalize='index').round(2))

# ALSFRS-R totale vs. stadio SuStaIn (Spearman)
r, p = spearmanr(df_clinical['sustain_stage'], df_clinical['alsfrs_total'])
print(f'ALSFRS-R vs. stadio: Spearman r={r:.3f}, p={p:.4f}')

# ALSFRS-R velocity vs. stadio
r, p = spearmanr(df_clinical['sustain_stage'], df_clinical['alsfrs_velocity'])
print(f'ALSFRS-R velocity vs. stadio: Spearman r={r:.3f}, p={p:.4f}')

# ECAS vs. stadio
r, p = spearmanr(df_clinical['sustain_stage'].dropna(),
                 df_clinical['ecas_total'].dropna())
print(f'ECAS vs. stadio: Spearman r={r:.3f}, p={p:.4f}')
```

### Confronto sottotipi su variabili continue

```python
# Mann-Whitney se 2 sottotipi, Kruskal-Wallis se >2
from scipy.stats import kruskal

for var in ['nfl_plasma', 'ecas_total', 'alsfrs_velocity', 'ptau217']:
    groups = [df_clinical.loc[df_clinical['sustain_subtype'] == k, var].dropna().values
              for k in range(N_ottimale)]
    if N_ottimale == 2:
        stat, p = mannwhitneyu(groups[0], groups[1])
        test = 'Mann-Whitney U'
    else:
        stat, p = kruskal(*groups)
        test = 'Kruskal-Wallis H'
    print(f'{var}: {test}, stat={stat:.1f}, p={p:.4f}')
```

### Distribuzione delle mutazioni per sottotipo

```python
print(pd.crosstab(df_clinical['sustain_subtype'],
                  df_clinical['mutation'],
                  normalize='index').round(2))
```

---

### Validazione biologica — SOD1 come controllo positivo

I SOD1 non entrano nel training — vengono stadiati con il modello cross-validato del dataset completo. L'aspettativa biologica è che si aggreghino nel sottotipo a pattern motorio puro.

```python
mask_sod1  = (df_clinical['mutation'] == 'SOD1').values
zdata_sod1 = zdata[mask_sod1]

ml_subtype_sod1, _, ml_stage_sod1, _, _ = \
    sustain_full.subtype_and_stage_individuals_newdata(
        zdata_sod1,
        samples_sequence_cval,
        samples_f_cval,
        N_samples=1000
    )

print('Distribuzione SOD1 per sottotipo:')
for k in range(N_ottimale):
    n = np.sum(ml_subtype_sod1 == k)
    pct = 100 * n / mask_sod1.sum()
    print(f'  Sottotipo {k+1}: {n}/{mask_sod1.sum()} pazienti ({pct:.1f}%)')
```

---

### Event-Based Model (EBM) per staging mutation-specific

L'EBM (`kde_ebm`) viene applicato successivamente a SuStaIn per costruire modelli di staging specifici per mutazione, applicabili a soggetti presintomatici nei trial clinici.

```python
from kde_ebm.mixture_model import fit_all_gmm_models, get_prob_mat
from kde_ebm.mcmc import mcmc, BootstrapSamples
import numpy as np

# Input: matrice z-score del sottogruppo (es. C9orf72 sintomatici)
# e z-score dei controlli armonizzati (per il modello GMM)
zdata_c9   = zdata[mask_c9]
zdata_ctrl = np.load('zscores_controls_roi.npy')   # da calcolare analogamente a Step 8
