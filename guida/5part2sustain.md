# Parte 2 — SuStaIn: Teoria, Fitting e Interpretazione

**Input**: `zscores_roi_als_final.csv` — matrice $\mathbf{Z} \in \mathbb{R}^{N_{\text{ALS}} \times R}$ (output Step 8)  
**Strumenti**: Python 3, pySuStaIn, ambiente conda `sustain_tutorial_env`

---

## 1. Cos'è SuStaIn — intuizione di base

SuStaIn (*Subtype and Stage Inference*) è un algoritmo di apprendimento non supervisionato che risponde a questa domanda:

> *Dati un gruppo di pazienti con profili biomarker diversi, esistono sottogruppi con sequenze di progressione distinte? E a che punto della sequenza si trova ciascun paziente?*

Le malattie neurodegenerative sono eterogenee. Nei pazienti ALS, i profili FDG-PET mostrano variabilità sistematica: alcuni hanno ipometabolismo prevalentemente motorio, altri coinvolgono il frontale, altri ancora il cervelletto. Questa variabilità può derivare da soggetti che seguono **sequenze di progressione biologicamente distinte** — cioè sottotipi — oppure da soggetti allo stesso punto di una progressione unica ma con variabilità individuale. SuStaIn distingue queste due ipotesi simultaneamente, senza richiedere etichette cliniche a priori.

**Output per ogni paziente**: un *sottotipo* (quale sequenza segue) e uno *stadio* (a che punto di quella sequenza si trova).

---

## 2. Il modello lineare z-score

La variante utilizzata è **ZscoreSuStaIn**, che descrive la progressione in termini di deviazioni standard dalla norma.

**Z-score events**: ogni biomarcatore può attraversare soglie di patologia definite (z = 1, 2, 3), producendo 3 eventi per ROI. Con 24 ROI si hanno 72 eventi totali per sottotipo.

**Stadi totali**: il numero di stadi in un modello è $N_{\text{ROI}} \times N_{\text{z-events}}$. Uno stadio non è una misura temporale — è una misura di quantità di patologia accumulata secondo la sequenza del sottotipo.

Esempio: se un sottotipo inizia con Precentral → SMA → DLPFC:
- Stadio 1: Precentral_L supera z = 1
- Stadio 2: Precentral_L supera z = 2
- Stadio 3: Precentral_L supera z = 3
- Stadio 4: SMA_L supera z = 1
- ...

---

## 3. Come funziona l'algoritmo

**Approccio gerarchico**: SuStaIn costruisce i sottotipi in modo incrementale — prima cerca il miglior modello con 1 sottotipo, poi il 2°, fino a Cmax. Questo è più stabile rispetto a cercare tutti i sottotipi contemporaneamente.

**MCMC**: la soluzione non è una singola sequenza ottimale, ma una distribuzione di sequenze plausibili stimata con Markov Chain Monte Carlo. Il risultato cattura l'incertezza nei dati rumorosi.

**Startpoints**: per evitare ottimi locali, l'algoritmo riparte da `N_startpoints = 25` punti iniziali diversi e confronta le soluzioni.

**Log-likelihood**: misura quanto un modello spiega i dati. Più sottotipi → log-likelihood sempre migliore, ma con rischio di overfitting — da qui la necessità del CVIC.

---

## 4. Preparazione del dataset

### Input di SuStaIn

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import pySuStaIn

# Carica la matrice z-score ROI (output Step 8)
df_z = pd.read_csv('zscores_roi_als_final.csv')

roi_labels = [c for c in df_z.columns if c != 'SubjectID']
zdata = df_z[roi_labels].values   # shape: (N_als, R)
N_subjects, N_biomarkers = zdata.shape
print(f'Pazienti: {N_subjects} | ROI: {N_biomarkers}')

# Verifica: nessun NaN (SuStaIn non li gestisce)
assert not np.isnan(zdata).any(), 'Presenza di NaN nella matrice — rimuovere o imputare prima'
```

### Definizione degli eventi z-score

```python
# 3 soglie per ogni ROI: z = 1 (patologia lieve), 2 (moderata), 3 (severa)
Z_vals = np.array([[1, 2, 3]] * N_biomarkers)

# Z_max: valore massimo atteso per la visualizzazione (~95° percentile)
Z_max = np.array([np.percentile(zdata[:, r], 95) for r in range(N_biomarkers)])
Z_max = np.maximum(Z_max, 5.0)   # minimo 5 per la scala del PVD
```

### Cosa NON serve a SuStaIn

SuStaIn non usa diagnosi cliniche (sono post-hoc), tempo dalla diagnosi, o variabili demografiche (già rimosse nel modello normativo al Step 7).

---

## 5. Configurazione e avvio

### Analisi 2a — Dataset completo (Cmax = 5)

```python
sustain_full = pySuStaIn.ZscoreSustain(
    data               = zdata,
    Z_vals             = Z_vals,
    Z_max              = Z_max,
    biomarker_labels   = roi_labels,
    N_startpoints      = 25,
    N_S_max            = 5,
    N_iterations_MCMC  = int(1e6),
    output_folder      = 'output_full',
    dataset_name       = 'ALS_full',
    use_parallel_startpoints = False
)

# Avvio — può richiedere ore/giorni, usare tmux sul server
samples_sequence, samples_f, ml_subtype, prob_ml_subtype, \
    ml_stage, prob_ml_stage, prob_subtype_stage = sustain_full.run_sustain_algorithm()
```

> **Nota su Cmax**: per il dataset completo ALS usare Cmax = 5. Per il sottogruppo C9orf72 (~70 soggetti) usare Cmax = 3. La regola pratica è che ogni sottotipo deve avere almeno 10 pazienti per evento z nel sottotipo più piccolo.

> **Checkpoint**: pySuStaIn salva automaticamente i risultati in file pickle in `output_folder/pickle_files/`. Se il calcolo viene interrotto, riprende dai pickle esistenti.

---

## 6. Valutazione della convergenza MCMC

Prima di interpretare i risultati, verificare che il MCMC abbia convergito.

```python
# Carica il pickle del modello a N sottotipi
s = 1   # es. modello a 2 sottotipi (0-indexed)
pk = pickle.load(open(f'output_full/pickle_files/ALS_full_subtype{s}.pickle', 'rb'))
samples_sequence = pk["samples_sequence"]
samples_f        = pk["samples_f"]
```

Un trace che **converge** mostra log-likelihood stabile nella seconda metà del campionamento. Oscillazioni ampie fino alla fine indicano stazionarietà non raggiunta — aumentare `N_iterations_MCMC`.

`samples_f` contiene le proporzioni di pazienti per sottotipo. Distribuzioni molto sovrapposte tra due sottotipi suggeriscono che non sono separabili — il numero ottimale di sottotipi probabilmente è inferiore.

---

## 7. Selezione del numero ottimale di sottotipi — CVIC

Più sottotipi → log-likelihood migliore, ma con rischio di overfitting. Il CVIC (*Cross-Validated Information Criterion*) bilancia fit e complessità.

```python
# Cross-validation a 10 fold
cv_fold_assignments = sustain_full.cross_validate_sustain_model(n_folds=10)
```

pySuStaIn produce automaticamente il grafico del CVIC per ogni numero di sottotipi. Cercare il punto in cui la curva raggiunge un plateau o un chiaro massimo. Se la curva sale fino a Cmax, considerare di aumentarlo (verificando che il dataset lo supporti).

---

## 8. Interpretazione — Positional Variance Diagrams

I PVD sono il principale strumento di interpretazione. Ogni PVD descrive un sottotipo.

### Come si legge un PVD

- **Righe**: eventi (biomarcatore × z-score). Con 24 ROI e 3 z-score = 72 righe
- **Colonne**: stadi (da 1 a N_stadi totali)
- **Colore**: probabilità che quell'evento avvenga in quello stadio
  - Rosso = z = 1 (patologia lieve)
  - Magenta = z = 2 (moderata)
  - Blu = z = 3 (severa)

Un PVD **ben definito** mostra struttura diagonale chiara: ogni evento ha alta probabilità in uno stadio specifico. Un PVD **disperso** indica alta incertezza nella sequenza — il sottotipo potrebbe non essere ben supportato dai dati.

```python
# Plot dei PVD per il modello ottimale
N_ottimale = 2   # da determinare con CVIC
s = N_ottimale - 1   # 0-indexed
pk = pickle.load(open(f'output_full/pickle_files/ALS_full_subtype{s}.pickle', 'rb'))

M = len(zdata)
pySuStaIn.ZscoreSustain._plot_sustain_model(
    sustain_full,
    pk["samples_sequence"],
    pk["samples_f"],
    M,
    subtype_order=tuple(range(N_ottimale))
)
```

### Pattern attesi nel contesto ALS-PET

- Un sottotipo con ipometabolismo prevalentemente **motorio** (Precentral, SMA) — atteso come il più comune e dominante nei SOD1
- Un sottotipo con coinvolgimento **frontale e cingolato** — correlato cognitivo, frequente in C9orf72
- Eventualmente un sottotipo con pattern **cerebellare diffuso** — frequente in C9orf72 e forme rare

---

## 9. Staging e subtyping dei pazienti

### Modello completo vs. cross-validato

**Modello completo** (`ml_subtype`, `ml_stage`): migliore per descrivere il dataset, ma non per predire su nuovi soggetti.

**Modello cross-validato**: da usare per stadiare soggetti nuovi (es. presintomatici C9orf72).

```python
# Distribuzione dei sottotipi
for k in range(N_ottimale):
    n = np.sum(ml_subtype == k)
    print(f'Sottotipo {k+1}: {n} pazienti ({100*n/N_subjects:.1f}%)')

# Stadio medio per sottotipo
for k in range(N_ottimale):
    mask_k = ml_subtype == k
    print(f'Sottotipo {k+1}: stadio medio = {np.mean(ml_stage[mask_k]):.1f} ± {np.std(ml_stage[mask_k]):.1f}')
```

```python
# Distribuzione degli stadi per sottotipo
fig, axes = plt.subplots(1, N_ottimale, figsize=(5*N_ottimale, 4))
for k in range(N_ottimale):
    mask_k = ml_subtype == k
    axes[k].hist(ml_stage[mask_k], bins=20, color=f'C{k}', edgecolor='white')
    axes[k].set_title(f'Sottotipo {k+1} (n={mask_k.sum()})')
    axes[k].set_xlabel('Stadio SuStaIn')
    axes[k].set_ylabel('N pazienti')
plt.tight_layout()
plt.savefig('distribuzione_stadi.png', dpi=150)
```

### Staging di soggetti nuovi (es. presintomatici C9orf72)

```python
# new_data: matrice (N_nuovi, N_ROI) con z-score dei nuovi soggetti
ml_subtype_new, prob_subtype_new, ml_stage_new, prob_stage_new, _ = \
    sustain_full.subtype_and_stage_individuals_newdata(
        new_data,
        samples_sequence_cval,
        samples_f_cval,
        N_samples=1000
    )
```

---

## Note metodologiche

### Perché si inverte il segno degli z-score

FDG-PET misura metabolismo glucidico. La patologia in ALS causa **ipometabolismo** → z-score negativi rispetto ai controlli. SuStaIn è costruito per biomarker che *aumentano* con la malattia. L'inversione di segno trasforma "quanto è bassa l'attività metabolica" in "quanto è avanzata la patologia metabolica".

### Perché Cmax = 3 per C9orf72 (e non 5)

Con ~70 soggetti e 24 ROI × 3 z-score = 72 eventi, Young et al. raccomandano almeno 10 pazienti per evento z nel sottotipo più piccolo. Con Cmax = 5, il sottotipo minimo potrebbe avere ~14 pazienti — insufficiente. Cmax = 3 garantisce ~23 pazienti per il sottotipo più piccolo.

### Validazione

In un dataset reale non esiste ground truth. La validazione si basa su: coerenza biologica dei sottotipi con la letteratura; replicabilità della cross-validation tra fold; validazione esterna (i SOD1 si aggregano nel sottotipo motorio puro?); correlazioni cliniche con variabili indipendenti (sede di esordio, ALSFRS-R velocity, ECAS).

---

**Continua in**: [Analisi post-hoc](parte2_analisi_posthoc.md)
