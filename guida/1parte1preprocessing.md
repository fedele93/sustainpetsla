# Parte 1 — Preprocessing (Step 1–4)

**Input**: immagini DICOM FDG-PET (ADNI controlli CN + pazienti ALS)  
**Output**: file `i*_w.nii` — volumi FDG-PET intensità-normalizzati in spazio MNI  
**Strumenti**: MATLAB R2025b + SPM25

---

## Panoramica della Parte 1

La Parte 1 segue due percorsi paralleli con finalità distinte:

- **Percorso quantitativo primario** (step 1→8): produce la matrice 
$\mathbf{Z} \in \mathbb{R}^{N_{\text{ALS}} \times R}$ di z-score ROI armonizzati, input diretto a SuStaIn.
- **Percorso visivo secondario** (step V): produce le mappe z-score voxel-wise `z_*.nii` per visualizzazione su MRIcroGL e analisi esplorative SPM second-level. Non alimenta SuStaIn.

| # | Step | Input → Output |
|---|------|----------------|
| 1 | Selezione controlli ADNI | — → lista soggetti CN baseline |
| 2 | Preprocessing e normalizzazione spaziale | DICOM → `s*_w.nii` |
| 3 | Smoothing 10 mm FWHM | `s*_w.nii` → `s*_w.nii` (smoothati) |
| 4 | Intensity normalization (whole brain mean) | `s*_w.nii` → `i*_w.nii` |
| 5 | Estrazione ROI (ctrl ADNI + pazienti ALS) | `i*_w.nii` tutti → $\mathbf{U}$ |
| 6 | Armonizzazione ComBat inter-scanner | $\mathbf{U}$ → $\tilde{\mathbf{U}}$ |
| 7 | Modello normativo ROI (OLS su ctrl armonizzati) | $\tilde{\mathbf{U}}_{\text{ctrl}}$ → $\hat{\boldsymbol{\beta}}$, $\hat{\boldsymbol{\sigma}}$ |
| 8 | Z-score ROI (pazienti ALS vs. norma) | $\tilde{\mathbf{U}}_{\text{ALS}}$ + modello → $\mathbf{Z}$ → **SuStaIn** |
| V | *(Opzionale)* Z-score voxel-wise per mappe visive | `i*_w.nii` ALS + modello voxel-wise → `z_*.nii` |

> **Nota sul percorso visivo (Step V)**: il modello normativo voxel-wise e il calcolo degli z-score voxel-wise sono descritti nell'[Appendice](parte1_appendice.md). La loro esecuzione è indipendente dal percorso primario e produce mappe ad uso esclusivamente visivo: non sono armonizzate con ComBat e non alimentano SuStaIn.

---

## Step 1 — Selezione della Popolazione di Controllo (ADNI)

I soggetti di controllo sono scaricati dall'Alzheimer's Disease Neuroimaging Initiative (ADNI). Criteri di selezione:

- **Stato cognitivo**: cognitivamente normale (CN) alla valutazione baseline
- **Formato immagine**: *Co-registered, Averaged* — il livello di preprocessing ADNI che include co-registrazione e media dei frame temporali, senza smoothing aggiuntivo. Scaricato in formato DICOM multi-slice.
- **Scansioni multiple**: in presenza di scansioni multiple per lo stesso soggetto, è stata selezionata la scansione con la data di acquisizione più antica (baseline), identificata tramite l'Image Data ID fornito dal portale ADNI.
- **Range di età**: 56–94 anni
- **N**: 231 soggetti

Sono stati esclusi i file provenienti da acquisizione ECAT o HRRT, e i soggetti con MMSE < 27 nei dati demografici.

> **Nota — varianza da scanner**: le immagini ADNI provengono da scanner PET di modelli e produttori diversi. Differenze nella risoluzione spaziale, nella PSF e nei protocolli di ricostruzione costituiscono una potenziale fonte di varianza non biologica. Anche le PET ALS sono state acquisite con scanner diversi rispetto a quelli dell'ADNI. Per tale motivo viene implementato ComBat (Step 6).

---

## Step 2 — Preprocessing e Normalizzazione Spaziale

### 2.1 — Conversione DICOM → NIfTI

Le immagini ADNI nel formato *Co-registered, Averaged* sono distribuite come serie DICOM multi-slice (tipicamente 47 file `.dcm` per soggetto). La conversione utilizza il modulo SPM `DICOM Import` con le seguenti impostazioni: output in formato NIfTI-1 `.nii` (non `.nii.gz`); struttura directory piatta (`root = 'flat'`); nessuna rotazione dell'orientamento durante l'import.

**Nota sul formato Analyze per i pazienti ALS**: i file provenienti da scanner di medicina nucleare italiani possono essere in formato Analyze 7.5 (coppia `.hdr`/`.img`). La conversione avviene tramite nibabel in Python (non tramite SPM DICOM Import). Dopo la conversione va eseguita una verifica visiva per escludere flip dell'asse x, che può verificarsi a causa dell'ambiguità nel campo `pixdim[0]` del formato Analyze.

### 2.2 — Correzione dell'origine dell'header

Le immagini pre-processate ADNI presentano frequentemente un'origine dell'header molto distante dallo spazio MNI (tipicamente y ≈ −240 mm), che impedisce la corretta inizializzazione della normalizzazione spaziale automatica. Prima della normalizzazione, l'origine di ciascun volume viene reimpostata al centro geometrico del volume modificando la matrice di trasformazione nell'header NIfTI, senza alterare i dati immagine. Questo passaggio consente all'algoritmo di normalizzazione di trovare una soluzione convergente.

```matlab
% Correzione origine per file NIfTI (controlli ADNI e pazienti ALS)
function reset_origin(fpath)
    V          = spm_vol(fpath);
    center_vox = V.dim / 2;
    T          = V.mat;
    T(1:3, 4)  = -T(1:3, 1:3) * center_vox(:);
    spm_get_space(fpath, T);
end
```

### 2.3 — Normalizzazione spaziale al template FDG-PET

La normalizzazione spaziale registra ogni volume nello spazio MNI. Il template usato è il **template FDG-PET di Della Rosa et al. (2014)** (file `TEMPLATE_FDGPET_100.nii`), costruito mediando 100 immagini FDG-PET di soggetti sani e pazienti con demenza. L'uso di un template FDG-PET — invece del template MNI T1 standard di SPM — riduce gli errori sistematici di registrazione dovuti alla differenza di contrasto tra modalità.

Il template viene impiegato esclusivamente come riferimento per la registrazione spaziale, non come riferimento metabolico; il riferimento metabolico è fornito dai controlli ADNI.

Viene utilizzato il modulo **Old Normalise: Estimate & Write** di SPM25 con i seguenti parametri:

- Trasformazione affine a 12 parametri seguita da deformazioni non lineari
- Basi DCT: 7×8×7 nelle dimensioni x, y, z; 16 iterazioni non lineari
- Cutoff: 25 mm; regolarizzazione: media (*medium*)
- Modalità *Preserve Concentrations*; interpolazione trilineare
- Dimensioni delle immagini normalizzate: voxel isotropici 2×2×2 mm³; FOV: x=−78:78, y=−112:76, z=−70:85 mm

> **Nota — versione SPM**: Della Rosa et al. utilizzavano SPM5. I parametri *Old Normalise* sono stati mantenuti identici in SPM12 e SPM25 per garantire la retrocompatibilità; il modulo è disponibile in SPM25 sotto *Tools → Old Normalise*.

Al termine della normalizzazione, ogni immagine viene ispezionata visivamente per verificare la convergenza della registrazione ed escludere immagini con field of view parziale o normalizzazione fallita.

---

## Step 3 — Smoothing Spaziale

Dopo la normalizzazione spaziale, ogni immagine viene convolta con un kernel gaussiano isotropico a **10 mm FWHM**, come riportato in Cabras et al. (2025). Lo smoothing migliora il rapporto segnale-rumore, riduce gli effetti di volume parziale residui e rende le immagini conformi all'assunzione di gaussianità richiesta per le analisi statistiche parametriche.

Lo smoothing è applicato dopo la normalizzazione spaziale — non prima — in modo da essere omogeneo tra controlli ADNI e pazienti ALS, indipendentemente dall'unità di misura originale delle immagini.

Gli script MATLAB completi per i Step 2 e 3 sono in [Appendice A.1](parte1_appendice.md#a1--preprocessing-e-normalizzazione-spaziale).

---

## Step 4 — Intensity Normalization Individuale

Per ogni immagine smoothata, ogni voxel viene diviso per la **media dell'intero cervello** (*whole brain mean*, WBM) del soggetto corrispondente. La WBM viene calcolata sui voxel intra-cerebrali, definiti come quelli con intensità superiore a zero dopo la normalizzazione spaziale.

La normalizzazione in intensità standardizza la magnitudine dei valori tra soggetti, correggendo per la variabilità nella quantità di radioattività iniettata e nei fattori di scala specifici del centro di acquisizione.

La scelta del whole brain mean come reference region è preferita rispetto al cervelletto o al ponte perché in ALS il cervelletto può essere coinvolto nella patologia — in particolare nelle forme C9orf72 (Tan et al., 2022). L'uso di una reference region potenzialmente patologica introdurrebbe un bias sistematico negli z-score delle regioni coinvolte.

**Naming convention output**: per ogni soggetto con file di input `s<SUBJECT_ID>_w.nii`, lo script produce `i<SUBJECT_ID>_w.nii` in una cartella separata.

Lo script MATLAB completo è in [Appendice A.2](parte1_appendice.md#a2--intensity-normalization).

---

## Applicazione della pipeline ai pazienti ALS

I pazienti ALS seguono esattamente la stessa sequenza di operazioni (Step 2–4) con gli stessi parametri. Le uniche differenze operative riguardano il formato del file di input:

- Se il file è già in formato NIfTI: si procede direttamente dalla correzione dell'origine (Step 2.2).
- Se il file è in formato Analyze 7.5 (`.hdr`/`.img`): va prima convertito in NIfTI con nibabel, con verifica visiva del flip sull'asse x, poi si procede dalla correzione dell'origine.

I dettagli della conversione Analyze → NIfTI per i pazienti ALS sono in [Appendice A.1b](parte1_appendice.md#a1b--conversione-analyze-nifti-per-pazienti-als).

---

**Continua in**: [ROI & ComBat (Step 5–6)](parte1_roi_combat.md)
