# Guida a SuStaIn per il Progetto FDG-PET in SLA

*Guida di riferimento — basata sul tutorial pySuStaIn e adattata al progetto ALS*

---

## Indice

 1. Cos'è SuStaIn — intuizione di base
 2. Il modello lineare z-score
 3. Come funziona l'algoritmo
 4. Preparazione del dataset
 5. Configurazione e avvio di SuStaIn
 6. Valutazione dell'output — MCMC e convergenza
 7. Selezione del numero ottimale di sottotipi — CVIC
 8. Interpretazione — Positional Variance Diagrams
 9. Staging e subtyping dei pazienti
10. Applicazione al progetto ALS-PET
11. Ambiente di lavoro e workflow pratico (valido per il mio pc linux + mio server linux)

---

## 1. Cos'è SuStaIn — intuizione di base

SuStaIn (*Subtype and Stage Inference*) è un algoritmo di apprendimento non supervisionato che risponde a questa domanda:

> *Dati un gruppo di pazienti con profili biomarker diversi, esistono sottogruppi con sequenze di progressione distinte? E a che punto della sequenza si trova ciascun paziente?*

### Il problema che SuStaIn risolve

Le malattie neurodegenerative sono eterogenee. Se prendi 200 pazienti ALS e confronti i loro profili FDG-PET, non vedi un quadro uniforme: alcuni hanno ipometabolismo prevalentemente motorio, altri coinvolgono il frontale, altri ancora mostrano pattern cerebellari. Questa variabilità, semplificando, può avere due cause:

- I pazienti sono **allo stesso stadio** di una progressione unica ma con variabilità individuale
- I pazienti seguono **sequenze di progressione diverse** — cioè esistono sottotipi biologici distinti

SuStaIn distingue queste due ipotesi simultaneamente. Non richiede etichette cliniche a priori (non gli dici "questi sono i motori, questi i cognitivi") — le trova lui dai dati.

### L'output fondamentale

Per ogni paziente, SuStaIn assegna:

- Un **sottotipo** (quale sequenza di progressione segue)
- Uno **stadio** (a che punto di quella sequenza si trova)

E per ogni sottotipo, produce una **sequenza** che descrive in quale ordine i biomarcatori diventano patologici.

---

## 2. Il modello lineare z-score

SuStaIn ha diverse implementazioni interne. Quella che utilizzeremo per  il nostro caso è il **modello lineare z-score**, che descrive la progressione di ciascun biomarcatore in termini di sigma (deviazioni standard dalla normalità).

### Concetti fondamentali

**Z-score events**: ogni biomarcatore può attraversare soglie di patologia definite. Per esempio si possono usare z = 1, 2, 3. Per ogni biomarcatore hai quindi 3 eventi possibili: "ha superato 1σ", "ha superato 2σ", "ha superato 3σ".

**Z_vals**: la matrice che definisce quali z-score sono possibili per ciascun biomarcatore. Tipicamente:

```python
Z_vals = np.array([[1, 2, 3]] * N_biomarkers)  # forma: (N_biomarkers, N_z_events)
```

> In questo caso ha usato 3 eventi possibili, per tutti, ma potrebbe non essere sempre utile avere lo stesso numero per tutti

**Z_max**: il valore massimo atteso per ciascun biomarcatore, usato come riferimento per la visualizzazione. In pratica si usa il 95° percentile osservato nei pazienti — per esempio un valore ragionevole potrebbe essere 5, ovvero si stima che a piena patologia un biomarcatore raggiunga z = 5.

**Stadi totali**: il numero di stadi in un modello con N biomarcatori e K z-score per biomarcatore è N × K. Con 18 ROI e 3 z-score hai 54 stadi. Ogni stadio corrisponde all'attraversamento di una singola soglia in un singolo biomarcatore.

> sono tanti, penso

### Cosa significa uno stadio

Uno stadio **non è una misura temporale**. Non si può dire "stadio 10 = 3 anni dalla diagnosi". È una misura di *quantità di patologia accumulata* secondo la sequenza del sottotipo.

Esempio: se la sequenza di un sottotipo inizia con motore_sinistro → motore_destro → DLPFC, allora:

- Stadio 1: motore_sinistro ha superato z = 1
- Stadio 2: motore_sinistro ha superato z = 2
- Stadio 3: motore_sinistro ha superato z = 3
- Stadio 4: motore_destro ha superato z = 1
- …

---

## 3. Come funziona l'algoritmo

### Approccio gerarchico

SuStaIn costruisce i sottotipi in modo incrementale:

1. Prima cerca il miglior modello con **1 sottotipo** (un singolo ordine di progressione)
2. Poi, tenendo fisso quello, cerca il migliore **2° sottotipo** aggiuntivo
3. Continua fino a Cmax sottotipi

Questo approccio gerarchico è più stabile computazionalmente rispetto a cercare tutti i sottotipi contemporaneamente.

### MCMC — come funziona

La soluzione di SuStaIn non è un singolo ordine di progressione, ma una **distribuzione di ordini plausibili**. Per stimare questa distribuzione si usa il MCMC (Markov Chain Monte Carlo).

Il processo:

1. Parti da una sequenza casuale
2. Proponi una piccola modifica (scambia due elementi nella sequenza)
3. Se la nuova sequenza spiega meglio i dati → accettala sempre
4. Se la spiega peggio → accettala con una certa probabilità (per non restare bloccato in un ottimo locale)
5. Ripeti milioni di volte

Il risultato è un campione di sequenze, distribuite secondo la loro plausibilità. Le sequenze frequenti nel campione sono quelle più supportate dai dati.

**Perché non basta trovare la sequenza migliore?** Perché con dati rumorosi (come quelli PET) ci possono essere molte sequenze quasi ugualmente buone. La distribuzione MCMC cattura questa incertezza.

**Numero di iterazioni consigliato**: `1e6`. Per esplorazioni preliminari si può usare `1e5`, ma per i risultati finali usa `1e6`.

### Startpoints — perché esistono

L'algoritmo può restare intrappolato in ottimi locali (soluzioni buone ma non ottime). Per evitarlo, SuStaIn riparte da punti iniziali diversi e confronta le soluzioni. Il parametro `N_startpoints` (tipicamente 25) definisce quante partenze diverse vengono usate.

### Log-likelihood

La qualità di una soluzione si misura con la log-verosimiglianza: quanto è probabile osservare i dati dato quel modello di progressione. SuStaIn massimizza questa quantità. Ogni volta che aggiungi un sottotipo la log-likelihood migliora (più parametri = migliore fit), ma il miglioramento deve essere giustificato — vedi CVIC.

---

## 4. Preparazione del dataset

### Input di SuStaIn

SuStaIn si aspetta una matrice di forma `(N_soggetti, N_biomarkers)` dove ogni valore è lo **z-score** di quel soggetto in quel biomarcatore rispetto ai controlli sani.

### Z-scoring rispetto ai controlli

Lo step 5 (modello normativo OLS) e step 6 (calcolo z-score) producono esattamente questo. Per ogni paziente ALS e ogni ROI hai già uno z-score calcolato come:

```
z = (valore_paziente - predizione_modello(età, sesso)) / sigma_residui
```

Con il segno invertito: `z_final = -z` perché in ALS l'ipometabolismo FDG dà z negativi, ma SuStaIn si aspetta che i valori **crescano** con la patologia (z positivi = più patologico).

### Struttura della matrice ROI

Dopo lo step 7 (estrazione AAL), avrai una matrice di dimensioni circa `(N_pazienti_ALS, 18-20)`. Ogni colonna è una ROI, ogni riga è un paziente.

### Cosa NON serve a SuStaIn

SuStaIn non usa:

- Le diagnosi cliniche (sono post-hoc)
- Il tempo dalla diagnosi
- Le variabili demografiche (sono già state rimosse nel modello normativo)

---

## 5. Configurazione e avvio di SuStaIn

### Installazione e ambiente

L'environment conda (già configurato per le esercitazioni) è `sustain_tutorial_env`. Per i dati reali potresti voler creare un environment dedicato con le stesse dipendenze.

> Qui potrebbe essere utile descrivere come configurare le variabili d’ambiente

```bash
conda activate sustain_tutorial_env
cd ~/path/al/tuo/notebook
jupyter notebook
```

### Import fondamentali

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import pySuStaIn
from pathlib import Path
```

### Caricamento dei dati

```python
# Carica la matrice z-score ROI
# Forma attesa: (N_pazienti, N_ROI)
data = pd.read_csv('zscores_roi_als.csv')

# Il DataFrame deve avere solo le colonne delle ROI (non ID soggetto ecc.)
roi_labels = ['PreCentrale_sin', 'PreCentrale_dx', 'SMA', 
              'Frontale_medio_sin', 'Frontale_medio_dx', ...]  # i tuoi nomi ROI

zdata = data[roi_labels].values  # array numpy (N, 18)
N_subjects, N_biomarkers = zdata.shape
print(f"Pazienti: {N_subjects}, ROI: {N_biomarkers}")
```

### Definizione degli eventi z-score

```python
# 3 eventi per ogni biomarcatore: z=1, z=2, z=3
Z_vals = np.array([[1, 2, 3]] * N_biomarkers)

# Valore massimo atteso per la visualizzazione (~95° percentile)
Z_max = np.array([5] * N_biomarkers)
```

> Bisogna scegliere il z_max che racchiuda il 95 % dei soggetti

### Creazione dell'oggetto SuStaIn

```python
sustain_input = pySuStaIn.ZscoreSustain(
    data            = zdata,
    Z_vals          = Z_vals,
    Z_max           = Z_max,
    biomarker_labels= roi_labels,
    N_startpoints   = 25,
    N_S_max         = 5,          # Cmax: numero massimo di sottotipi da cercare
    N_iterations_MCMC = int(1e6), # 1.000.000 iterazioni MCMC
    output_folder   = 'output_sustain',
    dataset_name    = 'ALS_FDG',
    use_parallel_startpoints = False
)
```

> **Nota su Cmax**: per il dataset completo ALS usa Cmax = 5. Per il sottogruppo C9orf72 (~70 soggetti) usa Cmax = 3. La regola pratica di Young et al. è che servono almeno 10 pazienti per ogni evento z nel sottotipo più piccolo. 
>
> (non sono certo che sia scritto da qualche parte nell'articolo, forse nei tutorial)

### Avvio del fitting

```python
# Questo può richiedere ore/giorni — usa tmux per avviare il lavoro sul server!
samples_sequence, samples_f, ml_subtype, prob_ml_subtype, \
    ml_stage, prob_ml_stage, prob_subtype_stage = sustain_input.run_sustain_algorithm()
```

**Importante**: SuStaIn salva automaticamente i risultati in file pickle nella cartella `output_folder`. Se il calcolo viene interrotto, può riprendere dal punto in cui si era fermato grazie a questi file.

---

## 6. Valutazione dell'output — MCMC e convergenza

### MCMC trace

Prima di interpretare i risultati, verifica che il MCMC abbia effettivamente convergito (N.B. Pare esista veramente questa parola). Il trace mostra l'andamento della log-likelihood durante le iterazioni MCMC.

```python
# Il trace viene salvato automaticamente da pySuStaIn
# Puoi caricarlo dal pickle e visualizzarlo:
s = 1  # 1 split = 2 sottotipi
pickle_filename = f'output_sustain/pickle_files/ALS_FDG_subtype{s}.pickle'
pk = pickle.load(open(pickle_filename, 'rb'))
samples_sequence = pk["samples_sequence"]
samples_f        = pk["samples_f"]
```

Un trace che **converge** mostra valori stabili nella seconda metà del campionamento. Se il trace continua a oscillare ampiamente fino alla fine, il MCMC non ha raggiunto stazionarietà — considera di aumentare `N_iterations_MCMC`.

> Magari questo lo si capisce meglio provando nello jupyter nootebook

### Distribuzione della frazione di pazienti per sottotipo

`samples_f` contiene i campioni MCMC della proporzione di pazienti in ciascun sottotipo. Se due sottotipi hanno distribuzioni molto sovrapposte, probabilmente non sono sottotipi separabili — questo suggerisce un numero ottimale di sottotipi più basso.

---

## 7. Selezione del numero ottimale di sottotipi — CVIC 

### Perché non basta la log-likelihood

Più sottotipi aggiungi, più la log-likelihood migliora — ma rischi l'overfitting: il modello impara a memoria i dati invece di catturare struttura biologica reale.

### Cross-validation e CVIC

La soluzione è la cross-validation. I dati vengono divisi in K fold (tipicamente 10). Per ogni numero di sottotipi da 1 a Cmax:

- Si addestra il modello su K-1 fold
- Si valuta la log-likelihood sull'1 fold escluso

Il **CVIC** (Cross-Validated Information Criterion) è la somma delle log-likelihood di test su tutti i fold. Il numero di sottotipi ottimale è quello con il **CVIC massimo**.

```python
# La cross-validation usa StratifiedKFold di sklearn
# pySuStaIn la esegue automaticamente con:
cv_fold_assignments = sustain_input.cross_validate_sustain_model(n_folds=10)
```

### Come leggere il risultato

pySuStaIn produce un grafico del CVIC per ogni numero di sottotipi. Stai cercando il punto in cui la curva raggiunge un plateau o mostra un chiaro massimo. Se la curva continua a salire fino a Cmax, considera di aumentare Cmax (ma prima verifica che il dataset lo supporti).

> Prima dicevamo che il gruppo più piccolo deve avere almeno 10 pazienti per evento, mi sembra

---

## 8. Interpretazione — Positional Variance Diagrams

I PVD sono il principale strumento di interpretazione di SuStaIn. Ogni PVD descrive un sottotipo.

### Come si legge un PVD

Il grafico ha:

- **Righe**: eventi (biomarcatore × z-score). Per esempio con 18 ROI e 3 z-score hai 54 righe
- **Colonne**: stadi (da 1 a N_stadi)
- **Colore**: probabilità che quell'evento avvenga in quello stadio
  - Rosso = z = 1 (patologia lieve)
  - Magenta = z = 2 (moderata)
  - Blu = z = 3 (severa)

Un PVD **ben definito** mostra una struttura diagonale chiara: ogni evento ha alta probabilità in uno stadio specifico e bassa negli altri. Indica che il sottotipo è identificato con certezza.

Un PVD **disperso** — con colori distribuiti su molti stadi per lo stesso evento — indica alta incertezza nella sequenza. Può significare che il sottotipo non è ben supportato dai dati.

### Analogia per ricordarlo

È come chiedere a 1000 esperti di mettere in ordine 54 carte. Se tutti concordano, ogni carta sta sempre nella stessa posizione → PVD con colori concentrati. Se ci sono disaccordi, la stessa carta appare in posizioni diverse → colori dispersi.

```python
# Per plottare i PVD del modello a N sottotipi ottimale:
s = N_ottimale - 1  # 0-indexed
pk = pickle.load(open(f'output_sustain/pickle_files/ALS_FDG_subtype{s}.pickle', 'rb'))
samples_sequence = pk["samples_sequence"]
samples_f        = pk["samples_f"]
M = len(zdata)
tmp = pySuStaIn.ZscoreSustain._plot_sustain_model(
    sustain_input, samples_sequence, samples_f, M,
    subtype_order=tuple(range(N_ottimale))
)
```

### Cosa cercare nel contesto ALS-PET

Nel  dataset potremmo trovare pattern coerenti con la letteratura:

- Un sottotipo con ipometabolismo prevalentemente **motorio** (precentrale, SMA) — atteso come il più comune e come dominante nei SOD1
- Un sottotipo con coinvolgimento **frontale e cingolato** — correlato cognitivo, frequente in C9orf72
- Eventualmente un sottotipo con pattern più diffuso incluso **cerebellare** — più frequente in C9orf72 e forme rare

---

## 9. Staging e subtyping dei pazienti

### Modello completo vs. cross-validato

Dopo il fitting hai due modelli:

**Modello completo** (`ml_subtype`, `ml_stage`): addestrato su tutti i dati. È il miglior modello per descrivere il dataset, ma può essere sovraadattato — non usarlo per fare predizioni su nuovi pazienti.

**Modello cross-validato** (`ml_subtype_cval`, `ml_stage_cval`): derivato dalla cross-validation. È il modello da usare per stadiare soggetti nuovi (es. presintomatici C9orf72).

### Recupero dei risultati di staging

```python
# Subtyping/staging dal modello completo
# ml_subtype: array (N_pazienti,) — sottotipo più probabile per ciascuno
# ml_stage:   array (N_pazienti,) — stadio più probabile
# prob_ml_subtype: array (N_pazienti,) — probabilità del sottotipo assegnato
# prob_ml_stage:   array (N_pazienti,) — probabilità dello stadio assegnato
# prob_subtype_stage: array (N_pazienti, N_sottotipi, N_stadi) — distribuzione completa

print("Distribuzione dei sottotipi:")
for k in range(N_sottotipi):
    n = np.sum(ml_subtype == k)
    print(f"  Sottotipo {k+1}: {n} pazienti ({100*n/N_subjects:.1f}%)")

print("\nStadio medio per sottotipo:")
for k in range(N_sottotipi):
    mask = ml_subtype == k
    print(f"  Sottotipo {k+1}: {np.mean(ml_stage[mask]):.1f} ± {np.std(ml_stage[mask]):.1f}")
```

### Visualizzazione della distribuzione degli stadi

```python
fig, axes = plt.subplots(1, N_sottotipi, figsize=(5*N_sottotipi, 4))
for k in range(N_sottotipi):
    mask = ml_subtype == k
    axes[k].hist(ml_stage[mask], bins=20, color=f'C{k}', edgecolor='white')
    axes[k].set_title(f'Sottotipo {k+1} (n={mask.sum()})')
    axes[k].set_xlabel('Stadio SuStaIn')
    axes[k].set_ylabel('N pazienti')
plt.tight_layout()
plt.savefig('distribuzione_stadi.png', dpi=150)
plt.show()
```

### Staging di soggetti nuovi (es. presintomatici)

Per applicare il modello a soggetti non inclusi nel training (come i presintomatici C9orf72), usa il modello cross-validato:

```python
# new_data: matrice (N_nuovi, N_ROI) con z-score dei nuovi soggetti
# Carica le sequenze cross-validate e applica
ml_subtype_new, prob_subtype_new, ml_stage_new, prob_stage_new, _ = \
    sustain_input.subtype_and_stage_individuals_newdata(
        new_data,
        samples_sequence_cval,
        samples_f_cval,
        N_samples=1000
    )
```

---

## 10. Applicazione al progetto ALS-PET

### Schema del flusso dati

```
PARTE 1 (già completata o in corso)
Steps 1-7 MATLAB/SPM25
      ↓
matrice z-score ROI: (N_ALS, 18-20)
      ↓
PARTE 2 — questo notebook
      ↓
Analisi 2a: SuStaIn intero dataset (Cmax=5)
      ↓
Analisi 2b: SuStaIn C9orf72 isolato (Cmax=3)
      ↓
Analisi 2c: confronto profili genetici
      ↓
Analisi 2d: validazione SOD1
      ↓
Analisi 2e/2f: correlazioni cliniche e biomarcatori
```

### Configurazione specifica per il dataset ALS

```python
# Analisi 2a — dataset completo
sustain_full = pySuStaIn.ZscoreSustain(
    data              = zdata_all,   # tutti i pazienti ALS
    Z_vals            = Z_vals,
    Z_max             = Z_max,
    biomarker_labels  = roi_labels,
    N_startpoints     = 25,
    N_S_max           = 5,
    N_iterations_MCMC = int(1e6),
    output_folder     = 'output_full',
    dataset_name      = 'ALS_full'
)

# Analisi 2b — solo C9orf72
mask_c9 = (df_clinical['mutation'] == 'C9orf72')
sustain_c9 = pySuStaIn.ZscoreSustain(
    data              = zdata_all[mask_c9],
    Z_vals            = Z_vals,
    Z_max             = Z_max,
    biomarker_labels  = roi_labels,
    N_startpoints     = 25,
    N_S_max           = 3,           # Cmax ridotto per n~70
    N_iterations_MCMC = int(1e6),
    output_folder     = 'output_c9',
    dataset_name      = 'C9orf72'
)
```

### Analisi post-hoc: correlazioni cliniche

Una volta ottenuti `ml_subtype` e `ml_stage`, puoi aggiungerli al dataframe clinico e analizzare le associazioni:

```python
df_clinical['sustain_subtype'] = ml_subtype
df_clinical['sustain_stage']   = ml_stage

# Sede di insorgenza per sottotipo
pd.crosstab(df_clinical['sustain_subtype'], df_clinical['onset_site'])

# ALSFRS-R per stadio (correlazione di Spearman)
from scipy.stats import spearmanr
r, p = spearmanr(df_clinical['sustain_stage'], df_clinical['alsfrs_total'])
print(f"Spearman r={r:.3f}, p={p:.4f}")

# Confronto sottotipi su variabili continue (Mann-Whitney se 2 sottotipi)
from scipy.stats import mannwhitneyu
sub0 = df_clinical[df_clinical['sustain_subtype']==0]['nfl_plasma']
sub1 = df_clinical[df_clinical['sustain_subtype']==1]['nfl_plasma']
stat, p = mannwhitneyu(sub0.dropna(), sub1.dropna())
print(f"NfL: sottotipo 0 vs 1, U={stat:.0f}, p={p:.4f}")
```

### Validazione biologica SOD1

```python
# I SOD1 non entrano nel training — li stadi con il modello cross-validato del dataset completo
mask_sod1 = (df_clinical['mutation'] == 'SOD1')
zdata_sod1 = zdata_all[mask_sod1]

ml_subtype_sod1, _, ml_stage_sod1, _, _ = \
    sustain_input.subtype_and_stage_individuals_newdata(
        zdata_sod1,
        samples_sequence_cval,
        samples_f_cval,
        N_samples=1000
    )

# Aspettativa: i SOD1 si aggregano nel sottotipo a pattern motorio puro
print("Distribuzione SOD1 per sottotipo:")
for k in range(N_ottimale):
    n = np.sum(ml_subtype_sod1 == k)
    print(f"  Sottotipo {k+1}: {n}/{mask_sod1.sum()} pazienti ({100*n/mask_sod1.sum():.1f}%)")
```

---

## 11. Ambiente di lavoro e workflow pratico

### Setup — dal laptop (Fedora) al server (ubuntu server)

I calcoli SuStaIn con `1e6` iterazioni richiedono ore o giorni. Usa il server per i calcoli pesanti:

```bash
# Dal laptop (Fedora): copia i dati sul server (ubuntu)
rsync -avz ~/path/dati/ fidel@geekfed:~/sustain/dati/

# Connettiti al server (ubuntu server) e avvia una sessione tmux
ssh fidel@geekfed
tmux new -s sustain_als

# Attiva l'environment e avvia jupyter
conda activate sustain_tutorial_env
cd ~/sustain/notebooks
jupyter notebook --no-browser --port=8888
```

```bash
# In un altro terminale sul laptop (Fedora KDE): forward della porta
ssh -L 8888:localhost:8888 fidel@geekfed

# Poi apri nel browser del laptop (Fedora KDE):
# http://localhost:8888
```

### Gestione delle sessioni tmux

```bash
# Riattacca una sessione esistente dopo essersi disconnessi
tmux attach -t sustain_als

# Detach senza chiudere (il calcolo continua)
# Ctrl+b, poi d

# Lista sessioni attive
tmux ls
```

### Checkpoint e pickle

pySuStaIn salva i risultati in file pickle dopo ogni numero di sottotipi. Se il kernel crasha o devi riavviare, il calcolo riprende dai pickle esistenti. Assicurati che `output_folder` punti alla stessa directory.

```bash
# Controlla i pickle esistenti
ls output_sustain/pickle_files/
# ALS_FDG_subtype0.pickle
# ALS_FDG_subtype1.pickle
# ...
```

### Struttura di directory consigliata

```
sustain_als/
├── notebooks/
│   ├── 01_setup_e_dati.ipynb
│   ├── 02_sustain_full_dataset.ipynb
│   ├── 03_sustain_c9.ipynb
│   └── 04_analisi_postehoc.ipynb
├── data/
│   ├── zscores_roi_als.csv        # output step 7 MATLAB
│   └── clinical_data.csv
├── output_full/
│   └── pickle_files/
├── output_c9/
│   └── pickle_files/
└── figures/
```

---

## Note metodologiche per la tesi

### Perché Cmax = 3 per C9orf72 (e non 5)

Con ~70 soggetti e 18 ROI × 3 z-score = 54 eventi, Young et al. raccomandano almeno 10 pazienti per evento z nel sottotipo più piccolo. Con Cmax = 5 il sottotipo minimo potrebbe avere ~14 pazienti — insufficiente per stimare sequenze di 54 eventi. Cmax = 3 garantisce che anche il sottotipo più piccolo abbia ~23 pazienti.

### Perché si inverte il segno degli z-score

FDG-PET misura metabolismo glucidico. La patologia in ALS causa **ipometabolismo** → z-score negativi rispetto ai controlli. SuStaIn è costruito per biomarker che *aumentano* con la malattia. L'inversione di segno è quindi concettualmente corretta: trasforma "quanto è bassa l'attività metabolica" in "quanto è avanzata la patologia metabolica".

### Ground truth e validazione

In un dataset reale non hai la ground truth (non sai quali sono i "veri" sottotipi). La validazione si basa su:

- **Coerenza biologica**: i sottotipi identificati hanno senso rispetto alla letteratura?
- **Replicabilità**: la cross-validation dà risultati stabili tra fold?
- **Validazione esterna**: i SOD1 si aggregano nel sottotipo motorio puro come atteso?
- **Correlazioni cliniche**: i sottotipi differiscono per variabili cliniche indipendenti (onset site, ALSFRS-R velocity, ECAS)?

---

*Guida basata sull'esercitazione pySuStaIn (Young et al. 2018, Aksman et al. 2021) e adattata al progetto FDG-PET in SLA — v1.1*
