# Parte 1 — Modello Normativo e Z-Score ROI (Step 7–8)

**Input**: `roi_uptake_controls_adni_harmonized.csv`, `roi_uptake_patients_als_harmonized.csv`  
**Output**: `normative_model_roi.npz`, `zscores_roi_als_final.csv` → **input a SuStaIn**  
**Strumenti**: Python 3, NumPy, pandas

---

## Step 7 — Modello Normativo ROI

### 7.1 Razionale

Il modello normativo descrive, per ogni regione anatomica di interesse $r$, il metabolismo FDG relativo atteso in un soggetto sano in funzione della sua età e del suo sesso. Stimato esclusivamente sui controlli ADNI armonizzati, esso fornisce il riferimento rispetto al quale calcolare la deviazione individuale di ciascun paziente ALS (Step 8).

A differenza del modello normativo voxel-wise (descritto nell'[Appendice](parte1_appendice.md) per il percorso visivo), che stima ∼400.000 equazioni OLS indipendenti, il modello normativo ROI stima $R = 24$ equazioni OLS. La ridotta dimensionalità lo rende computazionalmente trascurabile e statisticamente più stabile: con 231 controlli e 24 equazioni indipendenti, il rapporto soggetti/parametri è ∼77, ampiamente superiore alla soglia convenzionale di robustezza.

### 7.2 Formulazione

Per ogni ROI $r$, il modello OLS stimato sui $N_{\text{ctrl}}$ controlli armonizzati è:

$$\tilde{u}_{ir} \sim \beta_{0r} + \beta_{1r} \cdot \text{age}_{z,i} + \beta_{2r} \cdot \text{sex}_{c,i} + \varepsilon_{ir}$$

dove $\tilde{u}_{ir}$ è il valore di uptake armonizzato del controllo $i$ nella ROI $r$; $\text{age}_{z,i}$ è l'età standardizzata sulla media e deviazione standard del campione di controllo ($\mu_{\text{age}}$, $\sigma_{\text{age}}$); $\text{sex}_{c,i}$ è il sesso centrato sulla proporzione del campione ($\mu_{\text{sex}}$); $\varepsilon_{ir} \sim \mathcal{N}(0, \sigma_r^2)$.

Per ogni ROI $r$ vengono stimati e salvati:
- $\hat{\beta}_{0r}$: uptake relativo atteso alla media di età e sesso
- $\hat{\beta}_{1r}$: effetto dell'età standardizzata
- $\hat{\beta}_{2r}$: effetto del sesso centrato
- $\hat{\sigma}_r$: deviazione standard dei residui ($df = N_{\text{ctrl}} - 3$)

### 7.3 Implementazione

```python
import numpy as np
import pandas as pd

# -----------------------------------------------------------------------
# Carica la matrice dei controlli armonizzati (output Step 6)
# -----------------------------------------------------------------------
df_ctrl = pd.read_csv('roi_uptake_controls_adni_harmonized.csv')

roi_labels = [c for c in df_ctrl.columns
              if c not in ['SubjectID', 'age', 'sex', 'scanner_id', 'group']]
R = len(roi_labels)

# -----------------------------------------------------------------------
# Standardizza le covariabili — parametri da salvare per Step 8
# -----------------------------------------------------------------------
age_mean = df_ctrl['age'].mean()
age_std  = df_ctrl['age'].std()
sex_mean = df_ctrl['sex'].mean()   # proporzione M nel campione

age_z = (df_ctrl['age'].values - age_mean) / age_std
sex_c = df_ctrl['sex'].values - sex_mean

X = np.column_stack([np.ones(len(df_ctrl)), age_z, sex_c])
XtXinv = np.linalg.inv(X.T @ X)

# -----------------------------------------------------------------------
# Stima OLS per ogni ROI
# -----------------------------------------------------------------------
U_ctrl = df_ctrl[roi_labels].values   # shape: (N_ctrl, R)
beta_roi  = np.zeros((R, 3))
sigma_roi = np.zeros(R)

for r in range(R):
    y = U_ctrl[:, r]
    beta_roi[r, :] = XtXinv @ X.T @ y
    residui = y - X @ beta_roi[r, :]
    sigma_roi[r] = np.sqrt(np.sum(residui**2) / (len(y) - 3))

# -----------------------------------------------------------------------
# Salva i parametri del modello normativo ROI
# -----------------------------------------------------------------------
np.savez('normative_model_roi.npz',
         roi_labels = roi_labels,
         beta_roi   = beta_roi,
         sigma_roi  = sigma_roi,
         age_mean   = age_mean,
         age_std    = age_std,
         sex_mean   = sex_mean,
         N_ctrl     = len(df_ctrl),
         df         = len(df_ctrl) - 3)

print(f'Modello normativo ROI stimato su {len(df_ctrl)} controlli.')
print(f'  Età:   media={age_mean:.1f} anni, SD={age_std:.1f} anni')
print(f'  Sesso: {sex_mean*100:.1f}% M')
```

### 7.4 Output

File `normative_model_roi.npz` contenente:
- `beta_roi` — matrice $(R \times 3)$: coefficienti $[\hat{\beta}_0, \hat{\beta}_1, \hat{\beta}_2]$ per ogni ROI
- `sigma_roi` — vettore $(R)$: deviazione standard dei residui per ogni ROI
- `age_mean`, `age_std`, `sex_mean` — parametri di standardizzazione da applicare identici ai pazienti ALS

---

## Step 8 — Calcolo degli Z-Score ROI (Pazienti ALS)

### 8.1 Principio

Per ogni paziente ALS $i$ e ogni ROI $r$, lo z-score quantifica di quante deviazioni standard il metabolismo FDG relativo armonizzato si discosta da quello atteso per un soggetto sano della stessa età e dello stesso sesso:

$$z_{ir} = \frac{\tilde{u}_{ir} - \left(\hat{\beta}_{0r} + \hat{\beta}_{1r} \cdot \text{age}_{z,i} + \hat{\beta}_{2r} \cdot \text{sex}_{c,i}\right)}{\hat{\sigma}_r}$$

I parametri di standardizzazione ($\mu_{\text{age}}$, $\sigma_{\text{age}}$, $\mu_{\text{sex}}$) sono quelli stimati sul campione di controllo e salvati in `normative_model_roi.npz` — **non vengono ricalcolati sui pazienti**. Questo garantisce che lo z-score misuri la distanza dalla norma sana, non dalla distribuzione patologica.

### 8.2 Inversione del segno

FDG-PET misura metabolismo glucidico: la patologia in ALS causa *ipometabolismo*, quindi le regioni colpite presentano z-score **negativi** rispetto ai controlli. ZscoreSuStaIn è costruito per valori positivi crescenti con la patologia. Si applica pertanto l'inversione del segno:

$$z_{ir}^{\text{final}} = -z_{ir}$$

Dopo questa trasformazione, $z^{\text{final}} > 0$ indica ipometabolismo (patologia), $z^{\text{final}} = 0$ indica metabolismo nella norma. Nelle didascalie delle figure va sempre specificato esplicitamente il segno utilizzato.

### 8.3 Implementazione

```python
import numpy as np
import pandas as pd

# -----------------------------------------------------------------------
# Carica parametri del modello normativo (output Step 7)
# -----------------------------------------------------------------------
params     = np.load('normative_model_roi.npz', allow_pickle=True)
roi_labels = list(params['roi_labels'])
beta_roi   = params['beta_roi']    # shape: (R, 3)
sigma_roi  = params['sigma_roi']   # shape: (R,)
age_mean   = float(params['age_mean'])
age_std    = float(params['age_std'])
sex_mean   = float(params['sex_mean'])
R          = len(roi_labels)

# -----------------------------------------------------------------------
# Carica la matrice dei pazienti ALS armonizzati (output Step 6)
# -----------------------------------------------------------------------
df_als    = pd.read_csv('roi_uptake_patients_als_harmonized.csv')
U_als     = df_als[roi_labels].values   # shape: (N_als, R)
N_als     = len(df_als)

# Standardizza le covariabili dei pazienti con i parametri dei controlli
age_z_als = (df_als['age'].values - age_mean) / age_std
sex_c_als = df_als['sex'].values - sex_mean

# -----------------------------------------------------------------------
# Calcolo z-score per ogni paziente e ogni ROI
# -----------------------------------------------------------------------
Z_als = np.full((N_als, R), np.nan)
for i in range(N_als):
    x_i    = np.array([1.0, age_z_als[i], sex_c_als[i]])
    pred_i = beta_roi @ x_i          # predizione per il paziente i
    z_i    = (U_als[i, :] - pred_i) / sigma_roi
    Z_als[i, :] = -z_i               # inversione: positivo = ipometabolismo

# -----------------------------------------------------------------------
# Verifica: z-score dei controlli deve essere ≈ 0
# -----------------------------------------------------------------------
df_ctrl_harm = pd.read_csv('roi_uptake_controls_adni_harmonized.csv')
U_ctrl       = df_ctrl_harm[roi_labels].values
age_z_ctrl   = (df_ctrl_harm['age'].values - age_mean) / age_std
sex_c_ctrl   = df_ctrl_harm['sex'].values - sex_mean

Z_ctrl = np.full((len(df_ctrl_harm), R), np.nan)
for i in range(len(df_ctrl_harm)):
    x_i = np.array([1.0, age_z_ctrl[i], sex_c_ctrl[i]])
    pred_i = beta_roi @ x_i
    Z_ctrl[i, :] = -(U_ctrl[i, :] - pred_i) / sigma_roi

print('Verifica controlli (devono essere ≈ 0):')
print('  media z:', np.round(np.mean(Z_ctrl, axis=0), 3))
print('  SD z:   ', np.round(np.std(Z_ctrl,  axis=0), 3))
print()
print('Verifica pazienti ALS (media > 0 = ipometabolismo atteso):')
print('  media z:', np.round(np.mean(Z_als, axis=0), 3))

# -----------------------------------------------------------------------
# Salva la matrice z-score finale — input diretto a SuStaIn
# -----------------------------------------------------------------------
df_z = pd.DataFrame(Z_als, columns=roi_labels)
df_z.insert(0, 'SubjectID', df_als['SubjectID'].values)
df_z.to_csv('zscores_roi_als_final.csv', index=False)

print(f'\nMatrice z-score salvata: {N_als} × {R}')
print('File: zscores_roi_als_final.csv  →  input per SuStaIn (Parte 2)')
```

### 8.4 Controllo qualità pre-SuStaIn

Prima di passare a SuStaIn verificare che la matrice rispetti le assunzioni del modello:

1. **Media z nei controlli ≈ 0 per ogni ROI** — confermato dal codice sopra. Deviazioni > 0.1 suggeriscono errori nelle covariabili o nell'armonizzazione.
2. **SD z nei controlli ≈ 1 per ogni ROI** — atteso per costruzione. SD > 1.2 può indicare che ComBat ha rimosso meno varianza di batch del previsto.
3. **Media z nei pazienti ALS > 0 per le ROI di interesse** — conferma che l'ipometabolismo è correttamente codificato come valori positivi. ROI con media ≈ 0 o negativa nei pazienti potrebbero non portare informazione utile a SuStaIn.
4. **Assenza di NaN** — ROI con NaN per > 10% dei soggetti vanno escluse o imputate prima di SuStaIn, che non gestisce valori mancanti.

**Visualizzazione consigliata**:

- *Violin plot / boxplot* per ogni ROI sull'intera coorte ALS: le regioni motorie (Precentral, SMA) devono avere distribuzione spostata verso valori positivi.
- *Heatmap* soggetti × ROI: identifica soggetti outlier (righe con z-score uniformemente estremi) e ROI con NaN sistematici.
- *Matrice di correlazione inter-ROI*: regioni adiacenti (Precentral_L e Precentral_R) devono mostrare correlazioni positive moderate (r = 0.3–0.7). Correlazioni negative inattese segnalano errori sistematici nel preprocessing di un sottogruppo.

### 8.5 Output — Input a SuStaIn

**`zscores_roi_als_final.csv`** — matrice $\mathbf{Z} \in \mathbb{R}^{N_{\text{ALS}} \times R}$ con z-score ROI invertiti. Questa è la matrice che viene caricata all'inizio della Parte 2 come `data` per l'istanziazione di `pySuStaIn.ZscoreSustain`.

---

**Continua in**: [SuStaIn (Parte 2)](parte2_sustain.md)  
**Script completi**: [Appendice](parte1_appendice.md)
