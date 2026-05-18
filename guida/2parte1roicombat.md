# Parte 1 — Estrazione ROI e Armonizzazione ComBat (Step 5–6)

**Input**: file `i*_w.nii` — volumi FDG-PET intensità-normalizzati (controlli ADNI + pazienti ALS)  
**Output**: `roi_uptake_controls_adni_harmonized.csv`, `roi_uptake_patients_als_harmonized.csv`  
**Strumenti**: Python 3, nibabel, neuroCombat

---

## Step 5 — Estrazione dei Valori ROI

### 5.1 Razionale

L'estrazione ROI riduce ciascun volume intensità-normalizzato (`i*_w.nii`) a un vettore di $R$ valori scalari — uno per regione anatomica di interesse — mediante l'atlante AAL3. Questa riduzione dimensionale è il primo passo del percorso quantitativo primario e prerequisito all'armonizzazione ComBat (Step 6): ComBat opera su matrici numeriche compatte (soggetti × ROI), non su volumi voxel-wise.

L'operazione viene applicata a **tutti** i soggetti — controlli ADNI e pazienti ALS — producendo due matrici di uptake FDG relativo (intensità-normalizzato):

$$\mathbf{U}_{\text{ctrl}} \in \mathbb{R}^{N_{\text{ctrl}} \times R}, \qquad \mathbf{U}_{\text{ALS}} \in \mathbb{R}^{N_{\text{ALS}} \times R}$$

con $N_{\text{ctrl}} = 231$ controlli ADNI e $R = 24$ ROI selezionate. Le due matrici vengono concatenate in:

$$\mathbf{U} = \begin{pmatrix} \mathbf{U}_{\text{ctrl}} \\ \mathbf{U}_{\text{ALS}} \end{pmatrix} \in \mathbb{R}^{(N_{\text{ctrl}} + N_{\text{ALS}}) \times R}$$

che costituisce l'input al passo successivo (ComBat, Step 6).

La riduzione dimensionale persegue tre obiettivi complementari:

1. **Riduzione della dimensionalità**: da $V \approx 300.000$ voxel a $R = 24$ regioni, rendendo il modello SuStaIn computazionalmente trattabile.
2. **Riduzione del rumore di misura**: la media di migliaia di voxel all'interno di una regione riduce la varianza da rumore casuale proporzionalmente alla radice quadrata del numero di voxel.
3. **Interpretabilità**: i biomarker che entrano in SuStaIn corrispondono a strutture neuroanatomiche nominabili con significato clinico diretto.

### 5.2 Atlante AAL3

Come schema di parcellizzazione viene utilizzato l'atlante **AAL3** (*Automated Anatomical Labeling 3*, Rolls et al., 2020, *NeuroImage* 206:116189), versione `AAL3v1`, file `ROI_MNI_V7.nii` a risoluzione **2×2×2 mm³** nello spazio MNI152 — identica alla risoluzione dei volumi `i*_w.nii`. Non è quindi richiesto alcun ricampionamento dell'atlante prima dell'estrazione.

AAL3 aggiunge rispetto alle versioni precedenti la suddivisione del cingolato anteriore nelle sue componenti subgenuale, pregenuale e sopracallosale, e la parcellizzazione del talamo nei suoi nuclei specifici — entrambe rilevanti per i pattern ALS/ALS-FTD.

> **Nota — retrocompatibilità con AAL1**: i codici numerici per le regioni corticali e subcorticali classiche sono identici tra AAL1 e AAL3. I codici 35/36 (cingolato anteriore in AAL1) sono **vuoti** in AAL3: al loro posto, le componenti pregenuale e sopracallosale hanno codici 153–158.

### 5.3 Formalizzazione matematica

Sia $F_i(\mathbf{v})$ il valore di uptake FDG relativo del soggetto $i$ nel voxel $\mathbf{v}$ (intensità-normalizzato). Sia $\mathcal{M}_r$ l'insieme dei voxel appartenenti alla ROI 

$r$

secondo la maschera binaria AAL3, e 

$\mathcal{M}_{\text{FOV},i}$

la maschera del campo visivo del soggetto $i$ (i voxel non-NaN nel volume).

Il valore ROI del soggetto $i$ per la regione $r$ è la media dell'uptake sui voxel nell'intersezione tra maschera anatomica e FOV:

$$\bar{u}_{i,r} = \frac{1}{|\mathcal{M}_r \cap \mathcal{M}_{\text{FOV},i}|} \sum_{\mathbf{v} \in \mathcal{M}_r \cap \mathcal{M}_{\text{FOV},i}} F_i(\mathbf{v})$$

L'intersezione con $\mathcal{M}_{\text{FOV},i}$ è necessaria perché alcune scansioni PET hanno un FOV assiale ridotto che non copre le porzioni inferiori del cervelletto. Il denominatore si adatta automaticamente ai voxel effettivamente disponibili; la ROI risulta NaN solo se la copertura è zero.

> **Nota sulla scelta della media**: la media è lo stimatore standard nella letteratura SuStaIn (Young et al., 2018; Tan et al., 2022) e viene mantenuta per comparabilità. Come misura di robustezza può essere applicata una winsorizzazione preliminare (clipping dei valori oltre ±5 unità di WBM) prima dell'averaging.

> **Nota sul FOV**: scanner dedicati al neuroimaging con FOV assiale ≥ 15 cm coprono in genere l'intero cervello incluso il cervelletto. Scanner più datati o protocolli con FOV ridotto possono escludere i lobuli cerebellari inferiori. I soggetti con > 10% di ROI mancanti vengono esclusi dall'analisi.

### 5.4 Set di ROI selezionate

Il set comprende **24 ROI** (12 coppie bilaterali), selezionate sulla base della stadiazione neuropatologica TDP-43 di Brettschneider et al. (2013), dei cluster metabolici PET identificati da Tan et al. (2022), e del fenotipo clinico ALS/ALS-FTD nelle sue varianti genetiche.

Le regioni bilaterali vengono estratte separatamente (sinistra e destra) per consentire a SuStaIn di rilevare eventuali asimmetrie di progressione.

| Dominio | Nome AAL3 | Codice | Vol. (voxel) | Giustificazione |
|---|---|---|---|---|
| **Motorio primario** | Precentral_L / _R | 1 / 2 | 3526 / 3381 | Coinvolgimento UMN universale; ipometabolismo precoce in tutti i sottotipi |
| **Motorio supplementare** | Supp_Motor_Area_L / _R | 15 / 16 | 2147 / 2371 | Stadio TDP-43 I–II; presente in forme bulbari e spinali |
| **Frontale superiore mediale** | Frontal_Sup_Medial_L / _R | 19 / 20 | 2992 / 2134 | Componente cognitivo-motoria; pattern C9orf72 |
| **Prefrontale dorsolaterale** | Frontal_Mid_2_L / _R | 5 / 6 | 4507 / 4860 | Cluster FT (Tan et al.); compromissione esecutiva; C9orf72 |
| **Orbitofrontale mediale** | Frontal_Med_Orb_L / _R | 21 / 22 | 719 / 856 | ALS-FTD comportamentale; correlato di disinibizione |
| **Cingolato anteriore (pregenuale)** | ACC_pre_L / _R | 153 / 154 | 627 / 648 | Cluster CPT (Tan et al.); marcatore cognitivo-comportamentale |
| **Cingolato posteriore** | Cingulate_Post_L / _R | 39 / 40 | 463 / 335 | Default mode network; deterioramento cognitivo in stadio avanzato |
| **Temporale superiore** | Temporal_Sup_L / _R | 85 / 86 | 2296 / 3141 | ALS-FTD variante semantica/logopenica; C9orf72 |
| **Temporale medio** | Temporal_Mid_L / _R | 89 / 90 | 4942 / 4409 | Compromissione linguistica; memoria episodica |
| **Parietale superiore** | Parietal_Sup_L / _R | 63 / 64 | 2065 / 2222 | Cluster CPT; correlato visuospaziale ed esecutivo |
| **Cervelletto lobulo IV–V** | Cerebelum_4_5_L / _R | 101 / 102 | 1125 / 861 | Cluster CPT; C9orf72; forme atassiche |
| **Cervelletto lobulo VIII** | Cerebelum_8_L / _R | 107 / 108 | 1887 / 2308 | Pattern C9orf72; progressione tardiva |

**Regioni escluse**: il **talamo** (nuclei AAL3 con 10–265 voxel ciascuno) è troppo piccolo per produrre un biomarker stabile dopo smoothing 10 mm. La corteccia **occipitale** è tipicamente risparmiata nell'ALS. Il complesso **ippocampale** è rilevante nell'ALS-FTD ma si sovrappone ai pattern AD. **Putamen e caudato** sono considerati per un'analisi di sottogruppo nelle varianti parkinsoniane.

> **Da verificare prima dell'analisi finale**: le giustificazioni cliniche per ciascuna ROI vanno riconciliate con i risultati delle analisi esplorative preliminari. Il numero e la selezione delle ROI possono essere aggiornati dopo la prima esecuzione esplorativa.

### 5.5 Implementazione Python

```python
import numpy as np
import nibabel as nib
import pandas as pd
from pathlib import Path

# =====================================================================
# CONFIGURAZIONE
# =====================================================================
aal_path       = Path('/home/fluisi/Dati ricerca/AAL3/ROI_MNI_V7.nii')
uptake_dir_ctrl = Path('/home/fluisi/Dati ricerca/ADNI_nii_v3_intensity')
uptake_dir_als  = Path('/home/fluisi/Dati ricerca/ALS_nii_intensity')
output_csv_ctrl = Path('/home/fluisi/Dati ricerca/roi_uptake_controls_adni.csv')
output_csv_als  = Path('/home/fluisi/Dati ricerca/roi_uptake_patients_als.csv')

# =====================================================================
# TABELLA ROI — codici verificati su ROI_MNI_V7_vol.txt (AAL3, Rolls 2020)
# =====================================================================
roi_table = [
    (1,   'Precentral_L'),      (2,   'Precentral_R'),
    (15,  'Supp_Motor_Area_L'), (16,  'Supp_Motor_Area_R'),
    (19,  'Frontal_Sup_Medial_L'), (20, 'Frontal_Sup_Medial_R'),
    (5,   'Frontal_Mid_2_L'),   (6,   'Frontal_Mid_2_R'),
    (21,  'Frontal_Med_Orb_L'), (22,  'Frontal_Med_Orb_R'),
    (153, 'ACC_pre_L'),         (154, 'ACC_pre_R'),
    (39,  'Cingulate_Post_L'),  (40,  'Cingulate_Post_R'),
    (85,  'Temporal_Sup_L'),    (86,  'Temporal_Sup_R'),
    (89,  'Temporal_Mid_L'),    (90,  'Temporal_Mid_R'),
    (63,  'Parietal_Sup_L'),    (64,  'Parietal_Sup_R'),
    (101, 'Cerebelum_4_5_L'),   (102, 'Cerebelum_4_5_R'),
    (107, 'Cerebelum_8_L'),     (108, 'Cerebelum_8_R'),
]

# =====================================================================
# CARICA ATLANTE E VERIFICA CODICI
# =====================================================================
aal_img = nib.load(aal_path)
aal_vol = np.round(aal_img.get_fdata()).astype(int)
codici_presenti = set(np.unique(aal_vol[aal_vol > 0]))

mancanti = [code for code, _ in roi_table if code not in codici_presenti]
if mancanti:
    raise ValueError(f'Codici ROI non trovati nel volume AAL3: {mancanti}')
print(f'Verifica OK: tutti i {len(roi_table)} codici ROI presenti nel volume.')

roi_codes = [r[0] for r in roi_table]
roi_names = [r[1] for r in roi_table]
R = len(roi_table)

# =====================================================================
# FUNZIONE DI ESTRAZIONE (riutilizzata per ctrl e ALS)
# Input: directory con file i*_w.nii (o equivalente per ALS)
# pattern: pattern glob per trovare i file
# =====================================================================
def extract_roi_matrix(vol_dir, pattern, subj_id_fn):
    """
    vol_dir:    Path della cartella contenente i volumi NIfTI
    pattern:    stringa glob (es. 'i*_w.nii')
    subj_id_fn: funzione che estrae l'ID soggetto dal nome file
    Restituisce: DataFrame (N_soggetti × R+1) con SubjectID + valori ROI
    """
    files = sorted(vol_dir.glob(pattern))
    N = len(files)
    if N == 0:
        raise FileNotFoundError(f'Nessun file trovato: {vol_dir}/{pattern}')
    print(f'File trovati: {N}')

    matrix = np.full((N, R), np.nan)
    subj_ids = []

    for i, fpath in enumerate(files):
        vol = nib.load(fpath).get_fdata()
        subj_ids.append(subj_id_fn(fpath.stem))

        for r, (code, name) in enumerate(roi_table):
            mask_roi = (aal_vol == code)
            voxels   = vol[mask_roi & np.isfinite(vol)]
            if len(voxels) > 0:
                matrix[i, r] = np.mean(voxels)
            # altrimenti rimane NaN (FOV incompleto)

        if (i + 1) % 20 == 0 or (i + 1) == N:
            n_nan = np.sum(np.isnan(matrix[i]))
            print(f'  [{i+1}/{N}] {subj_ids[-1]} — ROI mancanti: {n_nan}/{R}')

    # Controllo qualità: soggetti con > 10% ROI mancanti
    nan_per_subj = np.sum(np.isnan(matrix), axis=1)
    soglia = int(0.10 * R)
    da_escludere = [subj_ids[j] for j in range(N) if nan_per_subj[j] > soglia]
    if da_escludere:
        print(f'\nATTENZIONE: {len(da_escludere)} soggetti con > 10% ROI mancanti:')
        for s in da_escludere:
            print(f'  {s}')

    df = pd.DataFrame(matrix, columns=roi_names)
    df.insert(0, 'SubjectID', subj_ids)
    return df

# =====================================================================
# ESTRAZIONE CONTROLLI ADNI
# Pattern: i<SUBJECT_ID>_w.nii
# =====================================================================
df_ctrl = extract_roi_matrix(
    uptake_dir_ctrl,
    pattern     = 'i*_w.nii',
    subj_id_fn  = lambda stem: stem.lstrip('i').rstrip('_w')
)

# =====================================================================
# ESTRAZIONE PAZIENTI ALS
# Adattare pattern e subj_id_fn alla nomenclatura del proprio dataset
# =====================================================================
df_als = extract_roi_matrix(
    uptake_dir_als,
    pattern    = 'i*_w.nii',
    subj_id_fn = lambda stem: stem.lstrip('i').rstrip('_w')
)

# =====================================================================
# AGGIUNGI COLONNE DEMOGRAFICHE (necessarie per ComBat)
# Le colonne 'age', 'sex', 'scanner_id', 'group' vanno caricate
# da un CSV demografico separato e unite per SubjectID
# =====================================================================
# df_demo_ctrl = pd.read_csv('demographics_adni.csv')
# df_ctrl = df_ctrl.merge(df_demo_ctrl[['SubjectID','age','sex','scanner_id']], on='SubjectID')
# df_ctrl['group'] = 0  # 0 = controllo

# df_demo_als = pd.read_csv('demographics_als.csv')
# df_als = df_als.merge(df_demo_als[['SubjectID','age','sex','scanner_id']], on='SubjectID')
# df_als['group'] = 1  # 1 = paziente ALS

# =====================================================================
# SALVA
# =====================================================================
df_ctrl.to_csv(output_csv_ctrl, index=False)
df_als.to_csv(output_csv_als,  index=False)
print(f'\nControlli salvati: {df_ctrl.shape}  →  {output_csv_ctrl}')
print(f'Pazienti ALS salvati: {df_als.shape}  →  {output_csv_als}')
```

### 5.6 Output

- `roi_uptake_controls_adni.csv` — matrice $\mathbf{U}_{\text{ctrl}}$ ($N_{\text{ctrl}} \times R$), una riga per soggetto ADNI
- `roi_uptake_patients_als.csv` — matrice $\mathbf{U}_{\text{ALS}}$ ($N_{\text{ALS}} \times R$), una riga per paziente ALS

Entrambi i file devono contenere le colonne `SubjectID`, `age`, `sex`, `scanner_id`, `group` — necessarie per ComBat al passo successivo.

---

## Step 6 — Armonizzazione Inter-Scanner con ComBat

### 6.1 Razionale

Le immagini FDG-PET che alimentano questa pipeline provengono da due fonti eterogenee: i controlli ADNI acquisiti su scanner diversi distribuiti su scala nordamericana, e i pazienti ALS acquisiti su scanner della rete italiana di medicina nucleare. Anche dopo la normalizzazione spaziale (Step 2) e la normalizzazione dell'intensità individuale (Step 4), permangono differenze sistematiche di origine tecnica attribuibili allo scanner: variazioni nella PSF, nei kernel di ricostruzione iterativa (OSEM), nella correzione dell'attenuazione, e nella calibrazione fotopicca. Queste costituiscono un *effetto batch* statistico che, se non rimosso, rischia di introdurre pattern artificiosi negli z-score e di condizionare i sottotipi identificati da SuStaIn.

### 6.2 Il modello ComBat

ComBat (Johnson et al., 2007; neuroCombat: Fortin et al., 2018) fattorizza ogni misura osservata come:

$$y_{ijk} = \alpha_k + \mathbf{x}_j^{\top}\boldsymbol{\beta}_k + \gamma_{ik} + \delta_{ik}\,\varepsilon_{ijk}$$

dove $\gamma_{ik}$ è l'effetto additivo del batch $i$ (spostamento di media) e $\delta_{ik}$ l'effetto moltiplicativo (differenza di scala). I parametri di batch sono stimati con approccio **Bayesiano empirico** (*borrowing of strength* tra ROI), il che rende il metodo robusto anche con batch di piccole dimensioni — situazione tipica dei centri ALS italiani.

Il dato armonizzato è:

$$\hat{y}_{ijk} = \frac{y_{ijk} - \hat{\alpha}_k - \mathbf{x}_j^{\top}\hat{\boldsymbol{\beta}}_k - \hat{\gamma}_{ik}}{\hat{\delta}_{ik}} + \hat{\alpha}_k + \mathbf{x}_j^{\top}\hat{\boldsymbol{\beta}}_k$$

### 6.3 Covariabili da proteggere

L'inclusione di `group` (controllo vs. paziente ALS) come covariabile è **critica**: senza di essa ComBat interpreterebbe le differenze biologiche malattia–controllo come effetti di batch e le rimuoverebbe. Vengono quindi protette:

- `group`: 0 = controllo CN, 1 = paziente ALS (categorica)
- `age`: età alla data di acquisizione (continua)
- `sex`: sesso biologico (categorica, 0 = F, 1 = M)

### 6.4 Definizione del batch

La variabile di batch identifica lo scanner o il sito per ogni soggetto. Per i controlli ADNI: metadati del portale (tabella `PETQC`, campo `SCANNERTYPE`). Per i pazienti ALS: centro di acquisizione. Se alcuni batch hanno < 5–10 soggetti, aggregare per modello di scanner o per sito.

### 6.5 Implementazione

```python
import numpy as np
import pandas as pd
from neuroCombat import neuroCombat

# -----------------------------------------------------------------------
# Carica le matrici ROI con demografici (output Step 5)
# -----------------------------------------------------------------------
df_ctrl = pd.read_csv('roi_uptake_controls_adni.csv')
df_als  = pd.read_csv('roi_uptake_patients_als.csv')

roi_labels = [c for c in df_ctrl.columns
              if c not in ['SubjectID', 'group', 'age', 'sex', 'scanner_id']]

# Concatena: prima i controlli, poi i pazienti (l'ordine è importante)
df_all = pd.concat([df_ctrl, df_als], axis=0, ignore_index=True)
n_ctrl = len(df_ctrl)

# neuroCombat richiede forma (N_features, N_subjects)
data_matrix = df_all[roi_labels].values.T

covars = pd.DataFrame({
    'batch': df_all['scanner_id'].values,
    'age':   df_all['age'].values,
    'sex':   df_all['sex'].values,
    'group': df_all['group'].values
})

# -----------------------------------------------------------------------
# Esegui ComBat
# -----------------------------------------------------------------------
combat_out = neuroCombat(
    dat              = data_matrix,
    covars           = covars,
    batch_col        = 'batch',
    categorical_cols = ['sex', 'group'],
    continuous_cols  = ['age']
)

# Ritrasposizione: (N_soggetti, N_ROI)
data_harmonized = combat_out['data'].T

# -----------------------------------------------------------------------
# Separa controlli e pazienti armonizzati
# -----------------------------------------------------------------------
roi_ctrl_harm = data_harmonized[:n_ctrl, :]
roi_als_harm  = data_harmonized[n_ctrl:, :]

# -----------------------------------------------------------------------
# Salva con colonne demografiche (necessarie per Step 7)
# -----------------------------------------------------------------------
df_ctrl_harm = pd.DataFrame(roi_ctrl_harm, columns=roi_labels)
df_ctrl_harm.insert(0, 'SubjectID', df_ctrl['SubjectID'].values)
for col in ['age', 'sex', 'scanner_id']:
    df_ctrl_harm[col] = df_ctrl[col].values
df_ctrl_harm.to_csv('roi_uptake_controls_adni_harmonized.csv', index=False)

df_als_harm = pd.DataFrame(roi_als_harm, columns=roi_labels)
df_als_harm.insert(0, 'SubjectID', df_als['SubjectID'].values)
for col in ['age', 'sex', 'scanner_id']:
    df_als_harm[col] = df_als[col].values
df_als_harm.to_csv('roi_uptake_patients_als_harmonized.csv', index=False)

print(f'ComBat completato.')
print(f'  Controlli armonizzati: {roi_ctrl_harm.shape}')
print(f'  Pazienti ALS armonizzati: {roi_als_harm.shape}')
```

> **Nota — batch con pochi soggetti**: se alcuni scanner hanno < 5–10 soggetti, usare `eb=False` in neuroCombat (OLS puro, senza regolarizzazione Bayesiana).

### 6.6 Controllo qualità post-armonizzazione

**Verifica visiva (PCA biplot)**: costruire una PCA della matrice ROI prima e dopo l'armonizzazione, colorando per scanner e per gruppo. Dopo ComBat i cluster per scanner devono dissolversi; la separazione biologica controlli/ALS deve mantenersi o rafforzarsi.

**Verifica quantitativa (ANOVA F-test)**:

```python
from scipy import stats

print(f"{'ROI':<30} {'F_batch_pre':>12} {'F_batch_post':>13} {'F_group_post':>13}")
print("-" * 72)

for idx, roi in enumerate(roi_labels):
    batches = covars['batch'].values
    groups  = covars['group'].values

    # F del batch prima dell'armonizzazione
    grp_pre = [data_matrix[idx, batches == b] for b in np.unique(batches)]
    F_pre, _ = stats.f_oneway(*grp_pre)

    # F del batch dopo l'armonizzazione
    grp_post = [data_harmonized[:, idx][batches == b] for b in np.unique(batches)]
    F_post, _ = stats.f_oneway(*grp_post)

    # F del gruppo (ctrl vs ALS) dopo l'armonizzazione
    F_grp, _ = stats.f_oneway(
        data_harmonized[:n_ctrl, idx],
        data_harmonized[n_ctrl:, idx]
    )
    print(f"{roi:<30} {F_pre:>12.2f} {F_post:>13.2f} {F_grp:>13.2f}")
```

Dopo ComBat: `F_batch_post` deve ridursi sostanzialmente rispetto a `F_batch_pre`; `F_group_post` deve mantenersi stabile o aumentare.

### 6.7 Output

- `roi_uptake_controls_adni_harmonized.csv` — $\tilde{\mathbf{U}}_{\text{ctrl}}$ ($N_{\text{ctrl}} \times R$) + colonne demografiche
- `roi_uptake_patients_als_harmonized.csv` — $\tilde{\mathbf{U}}_{\text{ALS}}$ ($N_{\text{ALS}} \times R$) + colonne demografiche

---

**Continua in**: [Normativo & Z-Score (Step 7–8)](parte1_normativo_zscore.md)
