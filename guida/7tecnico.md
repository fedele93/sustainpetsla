# Parte 2 — Workflow Tecnico

Setup server, tmux, struttura directory, checkpoint pickle.


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

# Fit del modello GMM misto (malato/sano) per ogni biomarcatore
mixtures = fit_all_gmm_models(zdata_c9, zdata_ctrl)

# Matrice di probabilità evento: P(evento avvenuto | soggetto i, biomarcatore j)
prob_mat = get_prob_mat(zdata_c9, mixtures)

# Campionamento MCMC per la sequenza degli eventi
mcmc_samples = mcmc(prob_mat)

# Bootstrap per l'incertezza
bs = BootstrapSamples(zdata_c9, zdata_ctrl, n_samples=100)
```

> **Nota**: l'EBM è applicato *dopo* SuStaIn, non in alternativa. SuStaIn identifica i sottotipi; EBM costruisce un modello di staging sequenziale all'interno di un sottogruppo genetico omogeneo, più adatto per lo staging presintomatico.

---

---

## Workflow Tecnico

---

### Setup — dal laptop (Fedora) al server (Ubuntu)

I calcoli SuStaIn con `1e6` iterazioni richiedono ore o giorni. Usare il server per i calcoli pesanti.

```bash
# Dal laptop (Fedora): copia dati e notebook sul server
rsync -avz ~/sustain/dati/   fidel@geekfed:~/sustain/dati/
rsync -avz ~/sustain/notebooks/ fidel@geekfed:~/sustain/notebooks/

# Connettiti al server e avvia una sessione tmux
ssh fidel@geekfed
tmux new -s sustain_als

# Attiva l'environment e avvia Jupyter
conda activate sustain_tutorial_env
cd ~/sustain/notebooks
jupyter notebook --no-browser --port=8888
```

```bash
# In un altro terminale sul laptop: forward della porta
ssh -L 8888:localhost:8888 fidel@geekfed
# Poi nel browser: http://localhost:8888
```

### Gestione delle sessioni tmux

```bash
# Detach senza chiudere (il calcolo continua in background)
# Ctrl+b, poi d

# Riattacca dopo disconnessione
tmux attach -t sustain_als

# Lista sessioni attive
tmux ls

# Crea una seconda sessione per monitoraggio
tmux new -s monitor
```

### Checkpoint e pickle

pySuStaIn salva i risultati in file pickle dopo ogni numero di sottotipi. Se il kernel crasha, il calcolo riprende dai pickle esistenti.

```bash
# Verifica stato dei pickle
ls output_full/pickle_files/
# ALS_full_subtype0.pickle   → modello 1 sottotipo
# ALS_full_subtype1.pickle   → modello 2 sottotipi
# ...

# Dimensione dei file (crescono con Cmax)
du -sh output_full/pickle_files/*.pickle
```

```python
# Carica un pickle per verificare il contenuto
import pickle
pk = pickle.load(open('output_full/pickle_files/ALS_full_subtype1.pickle', 'rb'))
print(pk.keys())
# dict_keys(['samples_sequence', 'samples_f', ...])
```

### Sincronizzazione risultati dal server al laptop

```bash
# Copia i pickle e le figure sul laptop per l'analisi
rsync -avz fidel@geekfed:~/sustain/output_full/ ~/sustain/output_full/
rsync -avz fidel@geekfed:~/sustain/figures/      ~/sustain/figures/
```

### Struttura di directory consigliata

```
sustain_als/
├── data/
│   ├── zscores_roi_als_final.csv       # output Step 8 — input SuStaIn
│   ├── roi_uptake_controls_adni_harmonized.csv
│   ├── normative_model_roi.npz
│   └── clinical_data.csv
├── notebooks/
│   ├── 01_roi_extraction.ipynb         # Step 5
│   ├── 02_combat.ipynb                 # Step 6
│   ├── 03_normative_zscore.ipynb       # Step 7–8
│   ├── 04_sustain_full.ipynb           # Analisi 2a
│   ├── 05_sustain_c9.ipynb             # Analisi 2b
│   └── 06_posthoc.ipynb                # Correlazioni, SOD1, EBM
├── output_full/
│   └── pickle_files/
├── output_c9/
│   └── pickle_files/
└── figures/
```

### Monitoraggio dell'uso delle risorse

```bash
# Sul server: controlla CPU e RAM durante il fitting
htop

# Stima del tempo rimanente (approssimativa)
# pySuStaIn non mostra progress bar nativa, ma scrive nei log
tail -f output_full/ALS_full_log.txt
```
