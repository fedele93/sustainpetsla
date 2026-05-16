# FDG PET in ALS: Pipeline Metodologica


---

## Panoramica

La pipeline si articola in due parti principali.

**Parte 1** è comune a tutto il dataset ALS e ai controlli ADNI: produce le mappe z-score voxel-wise e i valori ROI per ogni soggetto ALS.

**Parte 2** utilizza questi valori ROI come input per SuStaIn, con analisi sull'intero dataset ALS e su sottogruppi selezionati.

| PARTE 1 — Preprocessing & Z-score |
|---|
| 1. Selezione controlli ADNI |
| 2. Preprocessing e normalizzazione spaziale |
| 3. Smoothing 10 mm FWHM |
| 4. Intensity normalization (whole brain mean) |
| 5a. Estrazione ROI preliminare — controlli ADNI + pazienti ALS |
| 5b. Armonizzazione ComBat inter-scanner** |
| 6. Modello normativo OLS (età + sesso) — sui valori ROI armonizzati |
| 7. Calcolo z-score ROI — pazienti ALS vs. modello normativo |



---

# PARTE 1 — Preprocessing e Calcolo degli Z-Score

## Step 1 — Selezione della Popolazione di Controllo (ADNI)

I soggetti di controllo sono scaricati dall'Alzheimer's Disease Neuroimaging Initiative (ADNI). Criteri di selezione:

- **Stato cognitivo**: cognitivamente normale (CN) alla valutazione baseline
- **Formato immagine**: *Co-registered, Averaged* — il livello di preprocessing ADNI che include co-registrazione e media dei frame temporali, senza smoothing aggiuntivo. Scaricato in formato DICOM multi-slice.
- **Scansioni multiple**: in presenza di scansioni multiple per lo stesso soggetto, è stata selezionata la scansione con la data di acquisizione più antica (baseline), identificata tramite l'Image Data ID fornito dal portale ADNI.
- ~~**Range di età**: 56–94 anni~~
- ~~**N**: 365 soggetti dopo applicazione dei criteri~~ 

> ~~Le immagini ADNI in formato *Co-registered, Averaged* provengono da scanner PET di modelli e produttori diversi e, a seconda del sito di acquisizione, possono essere in unità diverse: conteggi grezzi (scanner con output in conteggi fotone, tipicamente DICOM), oppure unità calibrate SUV o Bq/ml (scanner con output in formato ECAT o HRRT). Questa eterogeneità di scala viene corretta dalla normalizzazione in intensità (Step 4), che riporta tutti i volumi alla stessa media unitaria.~~

Vienne deciso di escludere i file provenienti da acquisizione ECAT o HRRT. Inoltre ho escluso dei pazienti che nei demografic avevano dei mmse < 27. Per cui: 
- **Range di età**: 56-94 anni
- **N**: 231

> **Nota — varianza da scanner**: le immagini ADNI provengono da scanner PET di modelli e produttori diversi. Differenze nella risoluzione spaziale, nella PSF e nei protocolli di ricostruzione costituiscono una potenziale fonte di varianza non biologica. Anche le PET ALS sono state acquisite con scanner diversi ripetto a quelli dell' ADNI . Per tale motivo verrà implementato ComBat. 

> **N.B.** Nei demografics ADNI, sebbene siano stati selezionati esclusivamente pazienti CN (cognitive normal), erano presenti mmse < 27, che ora sono stati eliminati. Il dato dell'mmse era però presente per solo per una piccola quota di pazienti, pertanto non è possibile escludere che le altre immagini provengano da pazienti con un deficit cognitivo lieve. La presenza di questi soggetti potrebbe comunque non essere un problema: la popolazione dei pazienti ALS potrebbe avere copatologia Alzheimer o inziali deficit cognitivi; l'utilizzo di questa popolazione di riferimento per l'assegnazione dei z score nelle PET ALS potrebbe essere utile per evidenziare specificatamente "la variazione in senso SLA). 
>   Se si avranno a disposizione anchele pet dei controlli sani della società italiana di medicina nucleare si potrebbe ripetere la pipeline anche con quella popolazione di controllo per assegnarei z score.

---

## Step 2 — Preprocessing e Normalizzazione Spaziale

Il preprocessing delle immagini FDG-PET è eseguito in MATLAB R2025b con SPM25. Le operazioni sono applicate nell'ordine seguente a tutte le immagini (controlli ADNI e pazienti ALS (quest'ultima affermazione è da verificare una volta avute le pet SLA)).

### 2.1 — Conversione DICOM → NIfTI

I file DICOM di ciascun soggetto vengono convertiti in un unico volume NIfTI tridimensionale utilizzando il modulo *DICOM Import* di SPM25. Per le immagini *Co-registered, Averaged*, i DICOM rappresentano i singoli slice assiali del volume già mediato sui frame temporali; la conversione produce un unico file `.nii` per soggetto.

> I file dei controlli della società di Italiana di medicina nucleare sono nel formato precedente al NIfTI che veniva utilizzato da SPM (Analyze). 

> SPM utilizzava il formato Analyze 7.5 (.img / .hdr), un formato a due file separati: .hdr — file header con le informazioni sull'immagine (dimensioni, voxel size, tipo di dato, ecc.) .img — file binario con i dati volumetrici veri e propri. Era il formato standard di SPM nelle versioni fino a SPM2/SPM5 circa. Il NIfTI (.nii o coppia .nii.hdr/.nii.img) è stato introdotto nel 2004 come evoluzione dell'Analyze, unificando header e dati in un unico file e risolvendo alcune ambiguità del formato precedente — in particolare l'assenza di informazioni sull'orientamento spaziale (che era una fonte notoria di errori di flip left/right).

### 2.1b Contesto e motivazione

I dati FDG-PET dei pazienti ALS sono archiviati nel formato **Analyze 7.5** (coppia di file `.img` + `.hdr`). A differenza dei DICOM ADNI, che richiedono una conversione esplicita in NIfTI, i file Analyze sono già volumi tridimensionali ricostruiti pronti per l'elaborazione in SPM.

Questo non significa però che i file Analyze possano essere passati direttamente alla normalizzazione spaziale senza un'ispezione preliminare. L'Analyze 7.5 aveva una notoria limitazione nell'header: **non codifica informazioni sull'orientamento spaziale** (la matrice di rotazione e il verso degli assi è assente o ambiguo). In NIfTI questa informazione è esplicita nelle trasformazioni `sform`/`qform`. L'assenza dell'orientamento nell'Analyze significa che lo stesso file può essere interpretato con assi L→R oppure R→L a seconda del software, generando flip left/right invisibili all'ispezione visiva ma che invalidano l'analisi di lateralizzazione.

### Procedura di integrazione nella pipeline

**Passo A — Verifica dell'integrità della coppia `.img`/`.hdr`**

Prima di qualunque operazione, verificare che per ogni paziente esista una coppia di file con lo stesso nome base e che le dimensioni dichiarate nell'header siano consistenti con la dimensione del file `.img`. In MATLAB, `spm_vol` legge correttamente i file Analyze 7.5 e restituisce una struttura `V` con i campi `dim`, `mat` e `fname`; se il file è corrotto o la coppia è incompleta, `spm_vol` genererà un errore. In Python, `nibabel.load` accetta file `.img`/`.hdr` e restituisce un oggetto `Spm2AnalyzeImage`.

**Passo B — Conversione Analyze → NIfTI**

La conversione è necessaria per uniformare il formato con le immagini dei controlli ADNI e per codificare esplicitamente l'orientamento spaziale nell'header. Si raccomanda di convertire prima della correzione dell'origine (Passo C), così da operare sempre su file `.nii` da quel momento in poi.

In MATLAB con SPM25:

```matlab
% === Conversione Analyze → NIfTI con SPM25 ===
% Per ogni soggetto ALS con file .img/.hdr

als_dir = '/path/to/ALS_analyze/';
out_dir = '/path/to/ALS_nifti/';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

subj_list = dir(fullfile(als_dir, '*.hdr'));

for s = 1:length(subj_list)
    hdr_path = fullfile(als_dir, subj_list(s).name);
    img_path = strrep(hdr_path, '.hdr', '.img');

    % Leggi il volume Analyze
    V = spm_vol(hdr_path);
    vol_data = spm_read_vols(V);

    % Prepara header NIfTI
    V_out         = V;
    [~, fname, ~] = fileparts(subj_list(s).name);
    V_out.fname   = fullfile(out_dir, [fname '.nii']);
    V_out.descrip = sprintf('Converted from Analyze: %s', subj_list(s).name);

    % Scrivi come NIfTI
    spm_write_vol(V_out, vol_data);
    fprintf('Convertito: %s\n', fname);
end
```

In alternativa, la conversione può essere eseguita con **dcm2niix** (opzione `-f %p_%s`, accetta file `.img`/`.hdr` come input) o con nibabel in Python:

```python
import nibabel as nib
from pathlib import Path

als_dir = Path('/path/to/ALS_analyze')
out_dir = Path('/path/to/ALS_nifti')
out_dir.mkdir(exist_ok=True)

for hdr_file in sorted(als_dir.glob('*.hdr')):
    # nibabel carica la coppia .img/.hdr come Analyze
    img = nib.load(hdr_file)
    # Converti in NIfTI1 (preserva dati e matrice affine)
    nii = nib.Nifti1Image(img.get_fdata(), img.affine, img.header)
    out_path = out_dir / hdr_file.with_suffix('.nii').name
    nib.save(nii, out_path)
    print(f'Convertito: {hdr_file.name}')
```

> **Attenzione al flip left/right**: dopo la conversione, aprire un campione di file NIfTI con FSLeyes o MRIcroGL e verificare che la corteccia motoria destra sia anatomicamente a sinistra nell'orientamento radiologico (o a destra in quello neurologico), come atteso. Se si osserva un flip, il problema era nell'orientamento dell'Analyze originale: la correzione va fatta prima della conversione modificando la matrice `mat` del volume Analyze (scambiando il segno della prima colonna se gli assi x sono invertiti).

**Passo C — Correzione dell'origine dell'header (identica ai controlli ADNI)**

Una volta convertiti in NIfTI, i file dei pazienti ALS subiscono esattamente la stessa correzione dell'origine descritta per i controlli ADNI nello Step 2.2: l'origine viene reimpostata al centro geometrico del volume modificando la matrice di trasformazione nell'header. Questo garantisce che l'algoritmo Old Normalise di SPM trovi una sovrapposizione iniziale sufficiente con il template FDG-PET di Della Rosa.

```matlab
% === Correzione origine per file NIfTI dei pazienti ALS ===
% (identica alla funzione usata per i controlli ADNI)

nii_dir    = '/path/to/ALS_nifti/';
nii_files  = dir(fullfile(nii_dir, '*.nii'));

for s = 1:length(nii_files)
    fpath = fullfile(nii_dir, nii_files(s).name);
    V     = spm_vol(fpath);

    % Calcola centro geometrico e imposta come origine
    center_vox = V.dim / 2;
    T          = V.mat;
    T(1:3, 4)  = -T(1:3, 1:3) * center_vox(:);  % nuova traslazione
    spm_get_space(fpath, T);

    fprintf('Origine corretta: %s\n', nii_files(s).name);
end
```

**Passo D — Normalizzazione spaziale, smoothing e intensity normalization**

Da questo punto in poi, i file NIfTI dei pazienti ALS seguono esattamente la stessa pipeline dei controlli ADNI: normalizzazione spaziale al template Della Rosa con Old Normalise (Step 2.3), smoothing 10 mm FWHM (Step 3), e normalizzazione in intensità whole-brain (Step 4). I parametri sono identici; non è necessario alcun adattamento specifico per l'origine Analyze.

### Considerazioni sull'orientamento e il formato Analyze

Vale la pena spendere qualche parola sull'architettura del formato per evitare possibili errori. Il formato Analyze 7.5 codifica i dati in un file binario `.img` e le informazioni geometriche in un file `.hdr` di 348 byte. Il campo `pixdim` dell'header contiene le dimensioni dei voxel, ma il campo relativo all'orientamento (`hist.orient`) ha una semantica volutamente lasciata ambigua dalla specifica originale, tant'è che SPM nelle versioni storiche lo ignorava e assumeva sempre un orientamento implicito basato sulla convenzione radiologica (R→L per l'asse x). NIfTI 1.0 (2004) ha risolto questo problema introducendo le matrici `sform` e `qform`, che codificano l'orientamento in modo non ambiguo.

Quando nibabel legge un file Analyze e lo converte in NIfTI, costruisce la matrice affine a partire dai `pixdim` e da una convenzione di default: l'orientamento corretto dell'immagine risultante dipende quindi da come il software di acquisizione/ricostruzione del laboratorio di medicina nucleare scriveva il campo `pixdim[0]` (il segno della prima dimensione voxel, che alcune implementazioni usano per codificare la lateralità). È per questo che la verifica visiva del flip post-conversione (Passo B) è un passaggio non opzionale.

---


### 2.2 — Correzione dell'origine dell'header

Le immagini pre-processate ADNI presentano frequentemente un'origine dell'header molto distante dallo spazio MNI (tipicamente y ≈ −240 mm), che impedisce la corretta inizializzazione della normalizzazione spaziale automatica. Prima della normalizzazione, l'origine di ciascun volume viene reimpostata al centro geometrico del volume modificando la matrice di trasformazione nell'header NIfTI, senza alterare i dati immagine. Questo passaggio consente all'algoritmo di normalizzazione di trovare una soluzione convergente.

### 2.3 — Normalizzazione spaziale al template FDG-PET

La normalizzazione spaziale è il passaggio più critico della pipeline. La registrazione funziona in modo ottimale quando il template ha lo stesso contrasto dell'immagine da normalizzare: l'uso del template MNI T1 standard fornito con SPM introduce errori sistematici di registrazione per immagini FDG-PET.

Le immagini vengono pertanto normalizzate al **template FDG-PET di Della Rosa et al. (2014)** (file `TEMPLATE_FDGPET_100.nii`, disponibile pubblicamente su GitHub con licenza MIT), utilizzando il modulo *Old Normalise: Estimate & Write* di SPM25 con i seguenti parametri:

- Trasformazione affine a 12 parametri seguita da deformazioni non lineari
- Basi DCT: 7×8×7 nelle dimensioni x, y, z; 16 iterazioni non lineari
- Cutoff: 25 mm; regolarizzazione: media (*medium*)
- Modalità *Preserve Concentrations*; interpolazione trilineare
- Dimensioni delle immagini normalizzate: voxel isotropici 2×2×2 mm³; FOV: x=−78:78, y=−112:76, z=−70:85 mm

> Il template di Della Rosa et al. è stato costruito mediando 100 immagini FDG-PET normalizzate in spazio MNI, provenienti da soggetti sani e da pazienti con diverse forme di demenza, e validato su coorti indipendenti di AD, MCI, bvFTD e DLB. Viene qui impiegato esclusivamente come riferimento per la registrazione spaziale, non come riferimento metabolico; il riferimento metabolico è fornito dai controlli ADNI (oppure dai controlli della società italiana di medicina nucleare).

Al termine della normalizzazione, ogni immagine viene ispezionata visivamente per verificare la convergenza della registrazione ed escludere immagini con field of view parziale o normalizzazione fallita.

> **Nota — versione SPM**: Della Rosa et al. utilizzavano SPM5. I parametri *Old Normalise* sono stati mantenuti identici in SPM12 e SPM25 per garantire la retrocompatibilità; il modulo *Old Normalise* è disponibile in SPM25 sotto *Tools → Old Normalise*.

---

## Step 3 — Smoothing Spaziale

Dopo la normalizzazione spaziale, ogni immagine viene convolta con un kernel gaussiano isotropico a **10 mm FWHM**, come fatto in Cabras et al. (2025). Lo smoothing migliora il rapporto segnale-rumore, riduce gli effetti di volume parziale residui e rende le immagini conformi all'assunzione di gaussianità richiesta per le analisi statistiche parametriche.

Lo smoothing è applicato dopo la normalizzazione spaziale — non prima — in modo da essere omogeneo tra controlli ADNI e pazienti ALS, indipendentemente dall'unità di misura originale delle immagini.

---

## Step 4 — Intensity Normalization Individuale

Per ogni immagine smoothata, ogni voxel viene diviso per la **media dell'intero cervello** (*whole brain mean*, WBM) del soggetto corrispondente. La WBM viene calcolata sui voxel intra-cerebrali, definiti come quelli con intensità superiore a zero dopo la normalizzazione spaziale.

La normalizzazione in intensità standardizza la magnitudine dei valori tra soggetti, correggendo per la variabilità nella quantità di radioattività iniettata, nei fattori di scala specifici del centro di acquisizione. ~~e nelle differenze di unità di misura tra scanner (conteggi grezzi vs. SUV vs. Bq/ml- [*che credo sia relativo al formato del file nella versione grezza; quest'ultima considerazione non vale se utilizzo ADNI solo acquisiti in DICOM o solo i controlli ADNI*]).~~

La scelta del whole brain mean come reference region è preferita rispetto al cervelletto o al ponte perché in ALS il cervelletto può essere coinvolto nella patologia — vedi ad esempio le forme C9orf72 (Tan et al., 2022). L'uso di una reference region potenzialmente patologica introdurrebbe un bias sistematico negli z-score delle regioni coinvolte.

---
# Step 5 — Estrazione ROI Preliminare (Controlli + Pazienti ALS)

Prima di procedere al modello normativo, è necessario ridurre ciascun volume
intensità-normalizzato (file `i*_w.nii`, prodotto dallo Step 4) a un vettore di
valori scalari per regione anatomica. Questa operazione — descritta in dettaglio
nel §7 per le sue implicazioni metodologiche — viene qui anticipata rispetto
all'ordine narrativo della pipeline perché è **logicamente prerequisita**
all'armonizzazione ComBat (Step 4c): ComBat opera su una matrice numerica di
dimensioni ridotte (soggetti × ROI), non sui volumi voxel-wise.

L'estrazione ROI viene applicata a **tutti** i soggetti — controlli ADNI e
pazienti ALS — usando l'atlante AAL3 (`ROI_MNI_V7.nii`), con la procedura
descritta nel §7. Si ottengono due matrici:

- $\mathbf{U}_{\text{ctrl}} \in \mathbb{R}^{N_{\text{ctrl}} \times R}$: valori
  ROI dei controlli ADNI ($N_{\text{ctrl}} = 231$, $R =$ numero di ROI
  selezionate)
- $\mathbf{U}_{\text{ALS}} \in \mathbb{R}^{N_{\text{ALS}} \times R}$: valori ROI
  dei pazienti ALS

Le due matrici vengono concatenate verticalmente in un'unica matrice combinata
$\mathbf{U} \in \mathbb{R}^{(N_{\text{ctrl}} + N_{\text{ALS}}) \times R}$, che
costituisce l'input al ComBat.

---

# Step 5b — Armonizzazione Inter-Scanner con ComBat

## 5b.1 Razionale

Le immagini FDG-PET che alimentano questa pipeline provengono da due fonti
eterogenee: i controlli ADNI acquisiti su scanner diversi distribuiti su scala
nordamericana, e i pazienti ALS acquisiti su scanner della rete italiana di
medicina nucleare. Anche dopo la normalizzazione spaziale al template comune
(Step 2) e la normalizzazione dell'intensità individuale (Step 4), permangono
differenze sistematiche di origine tecnica attribuibili allo scanner: variazioni
nella *point spread function* (PSF) e nel relativo effetto di volume parziale,
differenze nei kernel di ricostruzione iterativa (OSEM) e nei parametri di
correzione dell'attenuazione, nonché eterogeneità nella calibrazione fotopicca
e nell'efficienza di rilevazione tra scanner di diversi produttori e generazioni.
Queste differenze non riflettono la biologia del soggetto ma costituiscono una
fonte di varianza non biologica — un *effetto batch* nel senso statistico del
termine — che, se non rimossa, rischia di introdurre pattern artificiosi nelle
mappe di z-score e, a cascata, di condizionare la stima del modello normativo
e i sottotipi identificati da SuStaIn.

Il problema è qualitativamente analogo a quello incontrato in genomica quando
si combinano dati di espressione genica da piattaforme o laboratori diversi.
L'armonizzazione con ComBat (*Combat-harmonization*) è la soluzione
metodologica sviluppata per questo contesto e successivamente adattata
all'imaging cerebrale: rimuove l'effetto additivo e moltiplicativo del batch
(scanner/sito) mentre preserva la variabilità biologica di interesse (età, sesso,
differenza malattia–controllo).

## 4b.2 Il modello ComBat

ComBat è stato originalmente proposto da Johnson et al. (2007, *Biostatistics*)
per l'armonizzazione di dati di microarray. La sua estensione all'imaging
cerebrale — denominata *neuroCombat* — è stata formalizzata da Fortin et al.
(2017, *NeuroImage*, per dati di spessore corticale da MRI strutturale) e da
Fortin et al. (2018, *NeuroImage*, per dati PET e dMRI).

Il modello assume che il valore osservato del biomarcatore $k$ nel soggetto $j$
appartenente al batch (scanner) $i$ possa essere decomposto come:

$$y_{ijk} = \alpha_k + \mathbf{x}_j^{\top}\boldsymbol{\beta}_k + \gamma_{ik} + \delta_{ik}\,\varepsilon_{ijk}$$

dove $\alpha_k$ è l'intercetta globale per il biomarcatore $k$; $\mathbf{x}_j$
è il vettore delle covariabili biologiche del soggetto $j$ (età, sesso,
appartenenza al gruppo) e $\boldsymbol{\beta}_k$ il vettore dei relativi
coefficienti; $\gamma_{ik}$ è l'**effetto additivo** del batch $i$ sul
biomarcatore $k$ (spostamento sistematico della media); $\delta_{ik}$ è
l'**effetto moltiplicativo** del batch $i$ sul biomarcatore $k$ (differenza di
scala o varianza); $\varepsilon_{ijk}$ è il residuo, assunto normale con media
zero e varianza unitaria.

I parametri di batch $\gamma_{ik}$ e $\delta_{ik}$ sono stimati con un
approccio **Bayesiano empirico** (*Empirical Bayes*, EB): invece di stimarli
indipendentemente per ogni biomarcatore — il che sarebbe instabile con pochi
soggetti per batch — la stima di ciascun $\gamma_{ik}$ viene regolarizzata
verso la distribuzione *a posteriori* costruita sfruttando l'informazione
condivisa tra tutti i biomarcatori. Questa proprietà di *borrowing of strength*
tra ROI rende ComBat particolarmente robusto in presenza di batch con un numero
ridotto di soggetti, che è esattamente la situazione di questo dataset (i
pazienti ALS sono distribuiti su più centri con numerosità variabile).

Il dato armonizzato è:

$$\hat{y}_{ijk} = \frac{y_{ijk} - \hat{\alpha}_k - \mathbf{x}_j^{\top}\hat{\boldsymbol{\beta}}_k - \hat{\gamma}_{ik}}{\hat{\delta}_{ik}} + \hat{\alpha}_k + \mathbf{x}_j^{\top}\hat{\boldsymbol{\beta}}_k$$

In parole: si rimuovono gli effetti additivi e moltiplicativi del batch, e si
ricolloca il dato sulla scala dell'effetto biologico globale e individuale.

## 5b.3 Covariabili da proteggere

L'inclusione corretta delle covariabili biologiche nel modello ComBat è
metodologicamente critica. Se la variabile `gruppo` (controllo *vs.* paziente
ALS) non viene inclusa come covariabile da preservare, ComBat interpreterà le
differenze metaboliche malattia–controllo come un effetto di batch — ovvero come
differenze tecniche — e le rimuoverà, distruggendo il segnale biologico
principale. Lo stesso vale per età e sesso, le cui associazioni con il
metabolismo cerebrale devono essere protette per non compromettere la stima del
modello normativo al passo successivo.

Le covariabili incluse nel modello ComBat come variabili di *protezione* sono
pertanto:

- `group`: variabile categorica (0 = controllo CN, 1 = paziente ALS)
- `age`: età alla data di acquisizione della scansione PET, trattata come
  variabile continua
- `sex`: sesso biologico, trattato come variabile categorica (0 = F, 1 = M)

## 5b.4 Definizione del batch

La variabile di batch identifica, per ogni soggetto, lo scanner (o il sito di
acquisizione) su cui è stata acquisita la PET. Per i controlli ADNI, questa
informazione è disponibile nei metadati del portale (tabella `PETQC` o campo
`SCANNERTYPE`). Per i pazienti ALS, la variabile di batch corrisponde al centro
di acquisizione. Qualora la granularità degli scanner ADNI risulti eccessiva
(molti scanner con pochissimi soggetti ciascuno, rendendo l'Empirical Bayes
instabile), è appropriato aggregare per modello di scanner (es. tutte le unità
Siemens Biograph mMR dello stesso sito) o per sito di acquisizione.

## 4c.5 Posizionamento nell'ordine della pipeline

La scelta di applicare ComBat **dopo** l'estrazione ROI — piuttosto che
direttamente sui volumi voxel-wise — è dettata da considerazioni pratiche e
metodologiche. Dal punto di vista pratico, applicare ComBat su volumi 3D con
circa 400.000 voxel per soggetto richiederebbe risorse computazionali ingenti
(memoria RAM proporzionale a $N \times V$, con $V \sim 4 \times 10^5$) senza
portare un vantaggio sostanziale: gli effetti di scanner si manifestano come
differenze sistematiche tra regioni e non come pattern rumorosi voxel-wise
casuali, e sono quindi catturati adeguatamente a livello ROI. Dal punto di vista
metodologico, applicare ComBat sulla matrice ROI è coerente con l'approccio
adottato in letteratura per pipeline analoghe che combinano popolazioni di
controllo di riferimento con coorti di pazienti acquisiti in centri diversi
(Pomponio et al., 2020, *NeuroImage*; Beer et al., 2020, *Human Brain Mapping*).

Di conseguenza, l'ordine logico della pipeline è:

1. Estrazione ROI su controlli e pazienti → matrice $\mathbf{U}$ (Step 4b)
2. Armonizzazione ComBat su $\mathbf{U}$ → matrice armonizzata $\tilde{\mathbf{U}}$ (Step 4c)
3. Separazione di $\tilde{\mathbf{U}}$ in controlli e pazienti
4. Stima del modello normativo sui valori ROI dei controlli armonizzati (Step 5, adattato)
5. Calcolo degli z-score ROI per i pazienti ALS rispetto al modello normativo armonizzato (Step 6, adattato)

**Nota sul rapporto con il modello normativo voxel-wise**: il documento descrive
anche una versione voxel-wise del modello normativo e degli z-score (Step 5–6
originali), utile per la visualizzazione di mappe z-score cerebrali complete.
Quella pipeline rimane valida per le analisi esplorative visive e per le
eventuali analisi second-level con SPM. La matrice ROI armonizzata $\tilde{\mathbf{U}}$
è invece l'input diretto a SuStaIn, e rappresenta il percorso quantitativo
primario.

## 5b.6 Implementazione computazionale

L'armonizzazione è eseguita con la libreria Python `neuroCombat` (Fortin et al.,
2018), installabile via `pip install neuroCombat`. Il codice opera sulla matrice
combinata controlli + pazienti; la separazione dei due gruppi avviene dopo
l'armonizzazione, non prima.

```python
import numpy as np
import pandas as pd
from neuroCombat import neuroCombat

# -----------------------------------------------------------------------
# 1. Carica le matrici ROI prodotte allo Step 4b
#    Ogni riga è un soggetto, ogni colonna una ROI
#    Le colonne non-ROI ('SubjectID', 'group', 'age', 'sex', 'scanner_id')
#    devono essere presenti nel DataFrame ma non passate come dati a ComBat
# -----------------------------------------------------------------------
df_ctrl = pd.read_csv('roi_controls_adni.csv')   # shape: (N_ctrl, N_ROI + meta)
df_als  = pd.read_csv('roi_patients_als.csv')    # shape: (N_als,  N_ROI + meta)

roi_labels = [c for c in df_ctrl.columns
              if c not in ['SubjectID', 'group', 'age', 'sex', 'scanner_id']]

# Concatena in ordine: prima i controlli, poi i pazienti
# (l'ordine è importante per la separazione successiva)
df_all = pd.concat([df_ctrl, df_als], axis=0, ignore_index=True)
n_ctrl = len(df_ctrl)

# -----------------------------------------------------------------------
# 2. Matrice dati in forma (N_features, N_subjects) — richiesta da neuroCombat
# -----------------------------------------------------------------------
data_matrix = df_all[roi_labels].values.T   # shape: (N_ROI, N_ctrl + N_als)

# -----------------------------------------------------------------------
# 3. DataFrame delle covariabili (una riga per soggetto)
#    'batch'  → scanner/sito (variabile da rimuovere)
#    'age'    → età (covariabile biologica da preservare, continua)
#    'sex'    → sesso (covariabile biologica da preservare, categorica)
#    'group'  → 0=controllo, 1=ALS (CRITICO: deve essere incluso per
#               proteggere le differenze malattia–controllo dall'essere
#               interpretate come effetto di batch)
# -----------------------------------------------------------------------
covars = pd.DataFrame({
    'batch': df_all['scanner_id'].values,
    'age':   df_all['age'].values,
    'sex':   df_all['sex'].values,
    'group': df_all['group'].values
})

# -----------------------------------------------------------------------
# 4. Esegui ComBat
#    categorical_cols: variabili trattate come fattori (dummies interne)
#    continuous_cols:  variabili trattate come covariate lineari
# -----------------------------------------------------------------------
combat_out = neuroCombat(
    dat              = data_matrix,
    covars           = covars,
    batch_col        = 'batch',
    categorical_cols = ['sex', 'group'],
    continuous_cols  = ['age']
)

# Output armonizzato: forma (N_ROI, N_soggetti) → ritrasponim a (N_soggetti, N_ROI)
data_harmonized = combat_out['data'].T

# -----------------------------------------------------------------------
# 5. Separa controlli e pazienti armonizzati
# -----------------------------------------------------------------------
roi_ctrl_harm = data_harmonized[:n_ctrl, :]   # shape: (N_ctrl, N_ROI)
roi_als_harm  = data_harmonized[n_ctrl:, :]   # shape: (N_als,  N_ROI)

# -----------------------------------------------------------------------
# 6. Salva le matrici armonizzate come CSV per i passi successivi
# -----------------------------------------------------------------------
df_ctrl_harm = pd.DataFrame(roi_ctrl_harm, columns=roi_labels)
df_ctrl_harm.insert(0, 'SubjectID', df_ctrl['SubjectID'].values)
df_ctrl_harm.to_csv('roi_controls_adni_harmonized.csv', index=False)

df_als_harm = pd.DataFrame(roi_als_harm, columns=roi_labels)
df_als_harm.insert(0, 'SubjectID', df_als['SubjectID'].values)
df_als_harm.to_csv('roi_patients_als_harmonized.csv', index=False)

print(f'ComBat completato. Matrici armonizzate salvate.')
print(f'  Controlli: {roi_ctrl_harm.shape}')
print(f'  Pazienti ALS: {roi_als_harm.shape}')
```

> **Nota implementativa — batch con pochi soggetti**: se alcuni scanner
> contengono un numero molto ridotto di soggetti (< 5–10), la stima
> dell'Empirical Bayes per quel batch può essere instabile. In questo caso
> è opportuno aggregare scanner dello stesso sito o dello stesso modello in
> un unico batch, oppure adottare la variante `neuroCombat` con l'opzione
> `eb=False` (OLS puro senza regolarizzazione Bayesiana), che è meno efficiente
> ma non richiede assunzioni distributive sui parametri di batch.

## 4c.7 Controllo qualità post-armonizzazione

Dopo l'armonizzazione è essenziale verificare che ComBat abbia rimosso l'effetto
di scanner senza alterare la struttura biologica del dataset. Il controllo
qualità comprende due verifiche complementari.

La prima è di tipo **visivo**: si costruisce una PCA biplot della matrice ROI
prima e dopo l'armonizzazione, colorando i punti per batch (scanner) e per
gruppo (controllo/ALS). Prima dell'armonizzazione è atteso che i soggetti si
raggruppino parzialmente per scanner; dopo l'armonizzazione i cluster per
scanner devono dissolversi, mentre la separazione tra controlli e pazienti —
se presente — deve essere preservata o rafforzata.

La seconda è **quantitativa**: si calcola, per ogni ROI, la statistica F di un
ANOVA a effetti fissi con il batch come fattore, prima e dopo l'armonizzazione,
e se ne confronta il valore p. Dopo ComBat, i valori F associati al batch
devono ridursi sostanzialmente; i valori F associati al gruppo (controllo/ALS),
calcolati separatamente, devono invece mantenersi stabili o aumentare, a
conferma che il segnale biologico è stato preservato.

```python
from scipy import stats

print(f"{'ROI':<30} {'F_batch_pre':>12} {'F_batch_post':>13} {'F_group_post':>13}")
print("-" * 72)

for i, roi in enumerate(roi_labels):
    batches  = covars['batch'].values
    groups   = covars['group'].values

    # F del batch prima dell'armonizzazione
    groups_pre  = [data_matrix[i, batches == b] for b in np.unique(batches)]
    F_pre, _    = stats.f_oneway(*groups_pre)

    # F del batch dopo l'armonizzazione
    groups_post = [data_harmonized[:, i][batches == b] for b in np.unique(batches)]
    F_post, _   = stats.f_oneway(*groups_post)

    # F del gruppo (ctrl vs ALS) dopo l'armonizzazione
    ctrl_vals   = data_harmonized[:n_ctrl, i]
    als_vals    = data_harmonized[n_ctrl:, i]
    F_grp, _    = stats.f_oneway(ctrl_vals, als_vals)

    print(f"{roi:<30} {F_pre:>12.2f} {F_post:>13.2f} {F_grp:>13.2f}")
```


## Step 6 — Modello Normativo Voxel-Wise

### Fase 6.1 — Stima del modello

Per ogni voxel $j$, viene stimato un modello di regressione lineare OLS (minimi quadrati ordinari) sui 231 controlli ADNI (*o il numero dei controlli  della società italiana di medicina nucleare*):

$$\text{uptake}_j \sim \beta_0 + \beta_1 \cdot \text{age}_z + \beta_2 \cdot \text{sex}_c$$ dove $$\text{age}_z$$  è l'età del soggetto alla data di acquisizione della scansione PET, centrata e standardizzata sulla media e deviazione standard del campione di controllo

 ( $$\mu_{\text{age}} = 75.1$$  anni,  $$\sigma_{\text{age}} = 6.4$$  anni); $$\text{sex}_c$$ è il sesso codificato come variabile binaria (M=1, F=0) e centrato sulla proporzione del campione (49.0% M).

Il sesso è incluso come covariata per controllare le differenze di metabolismo cerebrale tra uomini e donne che possono influenzare il modello normativo indipendentemente dall'età.

Per ogni voxel vengono stimati e salvati:
- $\hat{\beta}_{0j}$: intercetta (metabolismo predetto alla media di età e sesso)
- $\hat{\beta}_{1j}$: coefficiente dell'età standardizzata
- $\hat{\beta}_{2j}$: coefficiente del sesso centrato
- $\hat{\sigma}_j$: deviazione standard dei residui ($df = 362$)

La stima è eseguita in forma vettorizzata sull'intera matrice voxel × soggetti in un'unica operazione matriciale, previa applicazione di una maschera cerebrale (voxel con media > 0.1 nel campione di controllo; 403.117 voxel su 585.390 totali, 68.9%).

### Fase 6.2 — Output del modello

Il modello normativo produce tre file NIfTI:

- `normative_beta.nii` — volume 4D con i tre coefficienti β₀, β₁, β₂
- `normative_sigma.nii` — mappa della deviazione standard dei residui
- `normative_mean.nii` — mappa del metabolismo predetto alla media di età e sesso (β₀)

Un file `normative_params.mat` contiene i parametri di standardizzazione ($\mu_{\text{age}}$, $\sigma_{\text{age}}$, $\mu_{\text{sex}}$) necessari per applicare il modello ai pazienti ALS.
**(sono arrivato qui!!)**

---

## Step 7 — Calcolo degli Z-Score Individuali 

### 7.1 Principio generale e razionale

L'obiettivo di questo step è trasformare le immagini FDG-PET dei pazienti ALS da valori di uptake metabolico assoluto (o relativo, dopo la normalizzazione in intensità) in misure di *deviazione dalla norma* espressa in unità di deviazione standard. La logica è la stessa di qualsiasi z-score: si vuole rispondere alla domanda "di quante deviazioni standard il metabolismo di questo paziente in questo voxel si discosta da quello atteso per un soggetto sano della stessa età e dello stesso sesso?".

Il confronto avviene rispetto al **modello normativo voxel-wise** stimato al passo precedente (Step 5) sull'intera popolazione di controllo ADNI. Il modello normativo formalizza, per ogni voxel $j$ dello spazio MNI, quale sia il metabolismo atteso in funzione dell'età standardizzata

$\text{age}_{z}$ 

 e del sesso centrato

 $\text{sex}_{c}$ ,

 e quanto sia la variabilità fisiologica attorno a quella predizione (la deviazione standard dei residui $\hat{\sigma}_j$).

Il passo è concettualmente separato dalla stima del modello normativo perché i parametri del modello vengono stimati *una volta sola* sui controlli, poi applicati a ciascun paziente ALS individualmente. Questa separazione è importante: garantisce che le mappe z-score dei pazienti siano indipendenti l'una dall'altra e non "contaminino" la distribuzione di riferimento.

### 6.2 Formula del calcolo

Per ogni paziente ALS $i$ e ogni voxel $j$ all'interno della maschera cerebrale, lo z-score è calcolato come:

$$z_{ij} = \frac{x_{ij} - \left(\hat{\beta}_{0j} + \hat{\beta}_{1j} \cdot \text{age}_{z,i} + \hat{\beta}_{2j} \cdot \text{sex}_{c,i}\right)}{\hat{\sigma}_j}$$

dove $x_{ij}$ è il valore di uptake FDG (normalizzato in intensità e smoothato) del paziente $i$ nel voxel $j$; il numeratore è la differenza tra il valore osservato e il valore *predetto* dal modello normativo alla specifica età e sesso del paziente; il denominatore è la deviazione standard dei residui stimata sul campione di controllo.

I parametri di standardizzazione dell'età e del sesso ($\mu_{\text{age}}$, $\sigma_{\text{age}}$, $\mu_{\text{sex}}$ ) sono quelli salvati nel file `normative_params.mat` e sono applicati identici per tutti i pazienti:

$$\text{age}_{z,i} = \frac{\text{age}_i - \mu_{\text{age}}}{\sigma_{\text{age}}} \qquad \text{sex}_{c,i} = \text{sex}_i - \mu_{\text{sex}}$$

Questo è essenziale: usare i parametri del campione di controllo — e non ricalcolarli sui pazienti — assicura che lo z-score rifletta la distanza dalla norma sana, non dalla distribuzione patologica. Un paziente di 80 anni non viene confrontato con la media dei pazienti ALS, ma con la predizione del modello normativo per un ipotetico soggetto sano di 80 anni.

### 7.3 Inversione del segno e convenzione SuStaIn

L'FDG-PET è un marcatore di **metabolismo glucidico**: la patologia neurodegenerativa causa riduzione del metabolismo, quindi le regioni colpite presentano z-score *negativi* rispetto ai controlli sani. SuStaIn, nella sua implementazione a z-score (ZscoreSustain), è concepito per lavorare con valori positivi crescenti, dove valori più alti indicano maggiore patologia.

Per adottare questa convenzione, si applica l'inversione del segno:

$$z_{ij}^{\text{final}} = -z_{ij}$$

Dopo questa trasformazione, uno z-score finale di +2 in una regione frontale significa che il paziente ha un metabolismo 2 deviazioni standard *sotto* la norma attesa per la sua età, ovvero un ipometabolismo frontale significativo. Lo z-score finale di 0 significa metabolismo perfettamente nella norma. L'inversione non altera la struttura matematica del modello; è una convenzione di segno necessaria per la coerenza con ZscoreSustain, dove le soglie di "evento" sono definite come $z = +1$, $z = +2$, $z = +3$ (crescenti con la patologia).

> **Nota critica — implicazione per l'interpretazione**: i valori scalari che entrano in SuStaIn sono già con il segno invertito. Nelle analisi esplorative (distribuzione per ROI, heatmap), i valori positivi corrispondono all'ipometabolismo. È buona pratica indicare esplicitamente nelle didascalie delle figure se si utilizza la convenzione di segno originale (z negativo = ipometabolismo) o invertita (z positivo = ipometabolismo).

### 7.4 Output e organizzazione dei file

Per ogni paziente ALS $i$, viene prodotto un file NIfTI `z_<SUBJECT_ID>.nii` contenente la mappa z-score finale (con segno invertito) nello spazio MNI (voxel 2×2×2 mm³, dimensioni identiche alle immagini normalizzate dei controlli). I voxel fuori dalla maschera cerebrale vengono impostati a NaN (o 0) per chiarezza nella visualizzazione.

L'insieme di tutti i file `z_*.nii` costituisce la *libreria di mappe z-score* che alimenta sia lo Step 7 (estrazione ROI scalare per SuStaIn) sia eventuali analisi voxel-wise esplorative di gruppo (SPM second-level, se future analisi le richiedessero).

### 7.5 Implementazione computazionale

Il calcolo è eseguito in Python (ambiente Conda, notebook Jupyter), sfruttando la vettorizzazione NumPy per operare su tutti i pazienti in parallelo. Di seguito la struttura logica del codice; la versione completa è inclusa in Appendice A.4.

```python
import numpy as np
import nibabel as nib
import pandas as pd
from pathlib import Path

# === Carica parametri normativo ===
# (equivalente di normative_params.mat, salvato in formato .npz o .pkl)
params = np.load('normative_params.npz', allow_pickle=True)
age_mean   = float(params['age_mean'])
age_std    = float(params['age_std'])
sex_mean   = float(params['sex_mean'])

beta_img   = nib.load('normative_beta.nii')    # shape: (X, Y, Z, 3)
sigma_img  = nib.load('normative_sigma.nii')   # shape: (X, Y, Z)
beta_vol   = beta_img.get_fdata()              # β0, β1, β2 per ogni voxel
sigma_vol  = sigma_img.get_fdata()
affine     = beta_img.affine

# === Maschera cerebrale ===
mask = sigma_vol > 0   # voxel stimati nel modello normativo

# === Dati pazienti ALS ===
# df_als: DataFrame con colonne ['subject_id', 'age_at_scan', 'sex', 'filepath_nii']
df_als = pd.read_csv('als_subjects.csv')

output_dir = Path('zmaps_als')
output_dir.mkdir(exist_ok=True)

for _, row in df_als.iterrows():
    # --- Standardizza covariabili del paziente ---
    age_z_i  = (row['age_at_scan'] - age_mean) / age_std
    sex_c_i  = (1.0 if row['sex'] == 'M' else 0.0) - sex_mean

    # --- Carica immagine del paziente ---
    pet_img  = nib.load(row['filepath_nii'])
    x_vol    = pet_img.get_fdata()

    # --- Calcolo predizione normativa voxel-wise ---
    pred_vol = (beta_vol[..., 0]               # β0 (intercetta)
              + beta_vol[..., 1] * age_z_i     # β1 * età_z
              + beta_vol[..., 2] * sex_c_i)    # β2 * sesso_c

    # --- Z-score e inversione segno ---
    z_vol = np.full(x_vol.shape, np.nan)
    z_vol[mask] = (x_vol[mask] - pred_vol[mask]) / sigma_vol[mask]
    z_vol[mask] = -z_vol[mask]   # inversione: positivo = ipometabolismo

    # --- Salva NIfTI ---
    out_path = output_dir / f"z_{row['subject_id']}.nii"
    nib.save(nib.Nifti1Image(z_vol.astype(np.float32), affine), out_path)

print(f"Z-score calcolati per {len(df_als)} pazienti.")
```

> **Nota implementativa**: il ciclo `for` su ogni paziente è necessario perché le covariabili (età, sesso) differiscono tra pazienti; la predizione normativa è quindi calcolata su misura per ciascun soggetto. La parte computazionalmente onerosa — caricare `beta_vol` e `sigma_vol` — avviene una sola volta prima del ciclo.

### 7.6 Controllo qualità

Prima di procedere all'estrazione ROI (Step 7), è raccomandato eseguire un controllo qualità delle mappe z-score individuali. Per ogni paziente va verificata la distribuzione dei valori: in un soggetto senza grave patologia cerebrale diffusa, la distribuzione degli z-score cerebrali dovrebbe essere approssimativamente gaussiana centrata intorno a 0, con coda positiva nelle regioni ipometaboliche. Distribuzioni con media molto spostata (e.g. media > 1.5 su tutto il cervello) suggeriscono un problema di preprocessing (normalizzazione in intensità fallita, FOV incompleto, o normalizzazione spaziale divergente).

Si raccomanda inoltre di sovrapporre visivamente alcune mappe z-score individuali sul template MNI, per verificare che le regioni di ipometabolismo atteso (corteccia motoria, aree frontali) siano anatomicamente plausibili e che non vi siano artefatti di registrazione (iperintensità nei lobuli parietali posteriori, pattern a virgola nei giri parahippocampali, che sono segni tipici di registrazione fallita).

---

# Step 8 — Estrazione dei Valori ROI (Atlante AAL3)

## 8.1 Razionale della riduzione dimensionale

Le mappe z-score voxel-wise prodotte allo Step 6 descrivono ciascun paziente con un vettore di circa 200.000–400.000 valori scalari — uno per ogni voxel della maschera cerebrale nello spazio MNI a risoluzione 2 mm isotropica. Applicare SuStaIn direttamente a questa rappresentazione sarebbe computazionalmente proibitivo e statisticamente instabile: il modello ZscoreSuStaIn deve stimare una sequenza ordinata di "eventi" (ciascuno definito da una regione e da una soglia z), e con centinaia di migliaia di potenziali biomarker il problema diventerebbe non identificabile con campioni dell'ordine di qualche centinaio di soggetti.

La riduzione a un singolo valore scalare per regione anatomica di interesse (ROI) risolve questo problema perseguendo tre obiettivi complementari. Il primo è la **riduzione della dimensionalità**: da $V \approx 300.000$ voxel a $R \approx 20\text{–}25$ regioni, rendendo il modello SuStaIn computazionalmente trattabile. Il secondo è la **riduzione del rumore di misura**: la media di migliaia di voxel all'interno di una regione riduce la varianza da rumore casuale (termico, statistico, residuo di misallineamento) proporzionalmente alla radice quadrata del numero di voxel — è essenzialmente un averaging spaziale vincolato dall'anatomia. Il terzo è l'**interpretabilità**: i biomarker che entrano in SuStaIn corrispondono ora a strutture neuroanatomiche nominabili, e gli eventi di progressione del modello avranno direttamente un significato clinico e biologico leggibile.

---

## 8.2 Atlante utilizzato: AAL3

Come schema di parcellizzazione viene utilizzato l'atlante **AAL3** (*Automated Anatomical Labeling 3*, Rolls et al., 2020, *NeuroImage* 206:116189), nella versione `AAL3v1` distribuita dal Neurofunctional Imaging Group (GIN, UMR5296, Bordeaux) con licenza GNU GPL.

AAL3 è la terza generazione della famiglia AAL (Tzourio-Mazoyer et al., 2002), e incorpora tutte le regioni delle versioni precedenti aggiungendo strutture di particolare rilevanza per le malattie neurodegenerative: la suddivisione del cingolato anteriore nelle sue componenti subgenuale, pregenuale e sopracallosale; la parcellizzazione del talamo nei suoi nuclei specifici; e strutture sottocorticali come nucleo accumbens, substantia nigra, area tegmentale ventrale e nuclei del rafe. In totale comprende **166 regioni** con codici numerici da 1 a 170 (con alcuni salti dove le vecchie etichette AAL2, ad esempio i codici 35/36 per il cingolato anteriore e 81/82 per il talamo complessivo, sono stati svuotati e rimpiazzati dalle nuove suddivisioni).

Il file atlas utilizzato è `ROI_MNI_V7.nii`, con risoluzione **2×2×2 mm³** nello spazio MNI152, identica alla risoluzione delle mappe z-score prodotte allo Step 6. Non è quindi richiesto alcun ricampionamento dell'atlante prima dell'estrazione.

> **Nota — perché AAL3 e non AAL1/AAL2**: la letteratura PET-ALS più datata usa spesso AAL1 (116 regioni, codici 1–116). AAL3 è pienamente retrocompatibile con AAL1 per le regioni corticali e subcorticali classiche — i codici numerici sono identici — ma aggiunge la suddivisione del cingolato anteriore e del talamo che sono rilevanti per questo progetto. La scelta di AAL3 è giustificata dalla necessità di distinguere le sottoregioni del cingolato anteriore (in particolare la componente pregenuale, coinvolta nei pattern cognitivo-motori dell'ALS) e dalla sua adozione crescente come standard de facto nei lavori di neuroimaging recenti. Il riferimento da citare è Rolls et al. (2020).

---

## 8.3 Formalizzazione matematica

Sia

 $Z_i(\mathbf{v})$ 

la mappa z-score del paziente  $i$  nello spazio MNI, dove

 $\mathbf{v}$ 

indica le coordinate voxel. Sia

 $\mathcal{M}_r$ 

l'insieme dei voxel appartenenti alla ROI $$r$$ secondo la maschera binaria AAL3, e $\mathcal{M}_{\text{FOV},i}$ la maschera del campo visivo del paziente $i$ (i voxel non-NaN nella sua mappa z-score).

Il valore ROI del paziente $i$ per la regione $r$ è la media degli z-score sui voxel nell'intersezione tra maschera anatomica e FOV:

$$\bar{z}_{i,r} = \frac{1}{|\mathcal{M}_r \cap \mathcal{M}_{\text{FOV},i}|} \sum_{\mathbf{v} \in \mathcal{M}_r \cap \mathcal{M}_{\text{FOV},i}} Z_i(\mathbf{v})$$

L'intersezione con $\mathcal{M}_{\text{FOV},i}$ è necessaria perché alcune scansioni PET hanno un FOV assiale ridotto (tipicamente < 15 cm) che non copre le porzioni inferiori del cervelletto: il denominatore si adatta automaticamente ai voxel effettivamente disponibili, e la ROI risulta NaN solo se la copertura è zero.

> **Nota sul FOV**: il Field of View (campo visivo) è il volume fisico acquisito dallo scanner PET. Scanner dedicati al neuroimaging con FOV assiale ≥ 15 cm coprono in genere l'intero cervello incluso il cervelletto. Scanner più datati o protocolli con FOV ridotto possono escludere i lobuli cerebellari inferiori (VIII–IX). I soggetti con > 10% di ROI mancanti per FOV incompleto vengono esclusi dall'analisi principale (vedi §7.6).

L'insieme delle estrazioni per tutti i pazienti produce la **matrice dei biomarker**:

$$\mathbf{X} \in \mathbb{R}^{N \times R}$$

dove $N$ è il numero di pazienti ALS e $R$ il numero di ROI selezionate. Questa matrice $\mathbf{X}$ è l'input diretto al modello SuStaIn.

> **Nota sulla scelta della media**: alternative come la mediana o il 25° percentile sono state proposte per ridurre l'influenza di voxel con z-score estremi da artefatti locali. La media è tuttavia lo stimatore standard nella letteratura SuStaIn (Young et al., 2018; Tan et al., 2022) e viene qui mantenuta per garantire la comparabilità. Come misura di robustezza, può essere applicata una winsorizzazione preliminare (clipping dei valori oltre ±5 z) prima dell'averaging.

---

## 8.4 Set di ROI selezionate

Il set comprende **24 ROI** (11 coppie bilaterali sinistra/destra + 1 coppia di cingolato anteriore pregenuale), selezionate sulla base di tre livelli di evidenza convergenti: la stadiazione neuropatologica TDP-43 di Brettschneider et al. (2013), i cluster metabolici PET identificati da Tan et al. (2022) nella coorte ALS multicentrica, e il fenotipo clinico ALS/ALS-FTD nelle sue varianti genetiche.

Le regioni bilaterali vengono estratte separatamente (sinistra e destra) per consentire a SuStaIn di rilevare eventuali asimmetrie di progressione — documentate nelle forme bulbari e nelle varianti flail-arm/flail-leg, dove l'ipometabolismo tende a essere lateralizzato al lato di esordio.

| Dominio | Nome AAL3 | Codice | Vol. (voxel) | Giustificazione |
|---|---|---|---|---|
| **Motorio primario** | Precentral_L / _R | 1 / 2 | 3526 / 3381 | Coinvolgimento UMN universale; ipometabolismo precoce in tutti i sottotipi |
| **Motorio supplementare** | Supp_Motor_Area_L / _R | 15 / 16 | 2147 / 2371 | Stadio TDP-43 I–II; presente in forme bulbari e spinali |
| **Frontale superiore mediale** | Frontal_Sup_Medial_L / _R | 19 / 20 | 2992 / 2134 | Componente cognitivo-motoria del lobo frontale; pattern C9orf72 |
| **Prefrontale dorsolaterale** | Frontal_Mid_2_L / _R | 5 / 6 | 4507 / 4860 | Cluster FT (Tan et al.); compromissione esecutiva; C9orf72 |
| **Orbitofrontale mediale** | Frontal_Med_Orb_L / _R | 21 / 22 | 719 / 856 | ALS-FTD comportamentale; correlato di disinibizione |
| **Cingolato anteriore (pregenuale)** | ACC_pre_L / _R | 153 / 154 | 627 / 648 | Cluster CPT (Tan et al.); marcatore cognitivo-comportamentale |
| **Cingolato posteriore** | Cingulate_Post_L / _R | 39 / 40 | 463 / 335 | Default mode network; deterioramento cognitivo in stadio avanzato |
| **Temporale superiore** | Temporal_Sup_L / _R | 85 / 86 | 2296 / 3141 | ALS-FTD variante semantica/logopenica; C9orf72 |
| **Temporale medio** | Temporal_Mid_L / _R | 89 / 90 | 4942 / 4409 | Compromissione linguistica; memoria episodica |
| **Parietale superiore** | Parietal_Sup_L / _R | 63 / 64 | 2065 / 2222 | Cluster CPT; correlato visuospaziale ed esecutivo |
| **Cervelletto lobulo IV–V** | Cerebelum_4_5_L / _R | 101 / 102 | 1125 / 861 | Cluster CPT; C9orf72; forme atassiche |
| **Cervelletto lobulo VIII** | Cerebelum_8_L / _R | 107 / 108 | 1887 / 2308 | Pattern C9orf72; progressione tardiva |

I volumi in voxel sono riportati dalla colonna `vol_vox` del file `ROI_MNI_V7_vol.txt` incluso nella distribuzione AAL3 e verificati con nibabel sul file `ROI_MNI_V7.nii`.

**Regioni escluse e motivazione.** Il **talamo** è presente in AAL3 come 15 paia di nuclei distinti (codici 121–150), ma la maggior parte di essi ha dimensioni di 10–265 voxel — troppo piccole per produrre un biomarker stabile dopo smoothing a 10 mm FWHM, come avvertito dalla stessa documentazione AAL3. Il talamo viene quindi escluso dall'analisi primaria e potrà essere considerato come ROI aggregata (unione di tutti i nuclei) in un'analisi di sensibilità. La corteccia **occipitale** è tipicamente risparmiata nell'ALS e aggiunge solo rumore. Il complesso **ippocampale** è rilevante nell'ALS-FTD ma si sovrappone ai pattern dell'AD, con rischio di confounding in coorti miste. **Putamen e caudato** vengono considerati per un'analisi di sottogruppo nelle varianti parkinsoniane. L'**insula**, pur rilevante nei pattern FTD, ha in AAL3 una parcellizzazione anatomicamente discutibile per questo scopo.

> **Da verificare prima dell'analisi finale!**: le giustificazioni cliniche per ciascuna ROI vanno riconciliate con i risultati delle analisi esplorative preliminari. Il numero e la selezione delle ROI possono essere aggiornati dopo la prima esecuzione esplorativa.

---

## 8.5 Implementazione tecnica in Python (nibabel)

L'estrazione viene eseguita interamente in Python, senza passaggio per MATLAB/SPM, garantendo continuità con l'ambiente in cui gira SuStaIn e piena riproducibilità tramite il codice versionato su GitHub.

**File necessari** (tutti inclusi nella distribuzione AAL3):
- `ROI_MNI_V7.nii` — volume atlas 2×2×2 mm³ (file da caricare con nibabel)
- `ROI_MNI_V7_vol.txt` — tabella delle etichette con codici numerici verificati

```python
import numpy as np
import nibabel as nib
import pandas as pd
from pathlib import Path

# =====================================================================
# CONFIGURAZIONE — adattare i percorsi all'ambiente in uso
# =====================================================================
aal_path   = Path('/home/fluisi/Dati ricerca/AAL3/ROI_MNI_V7.nii')
zmap_dir   = Path('/home/fluisi/Dati ricerca/zmaps_als')   # output Step 6
output_csv = Path('/home/fluisi/Dati ricerca/roi_matrix_als.csv')

# =====================================================================
# TABELLA ROI — codici verificati su ROI_MNI_V7_vol.txt (AAL3, Rolls 2020)
# Formato: (codice_numerico_AAL3, 'nome_descrittivo')
# =====================================================================
roi_table = [
    # Motorio primario
    (1,   'Precentral_L'),
    (2,   'Precentral_R'),
    # Area motoria supplementare
    (15,  'Supp_Motor_Area_L'),
    (16,  'Supp_Motor_Area_R'),
    # Frontale superiore mediale
    (19,  'Frontal_Sup_Medial_L'),
    (20,  'Frontal_Sup_Medial_R'),
    # Prefrontale dorsolaterale
    (5,   'Frontal_Mid_2_L'),
    (6,   'Frontal_Mid_2_R'),
    # Orbitofrontale mediale
    (21,  'Frontal_Med_Orb_L'),
    (22,  'Frontal_Med_Orb_R'),
    # Cingolato anteriore pregenuale (AAL3: ACC_pre; i vecchi codici 35/36 sono vuoti in AAL3)
    (153, 'ACC_pre_L'),
    (154, 'ACC_pre_R'),
    # Cingolato posteriore
    (39,  'Cingulate_Post_L'),
    (40,  'Cingulate_Post_R'),
    # Temporale superiore
    (85,  'Temporal_Sup_L'),
    (86,  'Temporal_Sup_R'),
    # Temporale medio
    (89,  'Temporal_Mid_L'),
    (90,  'Temporal_Mid_R'),
    # Parietale superiore
    (63,  'Parietal_Sup_L'),
    (64,  'Parietal_Sup_R'),
    # Cervelletto lobulo IV-V
    (101, 'Cerebelum_4_5_L'),
    (102, 'Cerebelum_4_5_R'),
    # Cervelletto lobulo VIII
    (107, 'Cerebelum_8_L'),
    (108, 'Cerebelum_8_R'),
]

# =====================================================================
# VERIFICA PRELIMINARE — eseguire una volta prima dell'analisi reale
# Controlla che tutti i codici ROI esistano nel volume atlas
# =====================================================================
aal_img = nib.load(aal_path)
aal_vol = np.round(aal_img.get_fdata()).astype(int)
codici_nel_volume = set(np.unique(aal_vol[aal_vol > 0]))

codici_mancanti = [code for code, _ in roi_table if code not in codici_nel_volume]
if codici_mancanti:
    raise ValueError(f'Codici ROI non trovati nel volume AAL3: {codici_mancanti}')
else:
    print(f'Verifica OK: tutti i {len(roi_table)} codici ROI presenti nel volume.')

roi_codes = [r[0] for r in roi_table]
roi_names = [r[1] for r in roi_table]

# =====================================================================
# ESTRAZIONE — per ogni paziente e per ogni ROI
# =====================================================================
zfiles = sorted(zmap_dir.glob('z_*.nii'))
N      = len(zfiles)
R      = len(roi_table)

if N == 0:
    raise FileNotFoundError(f'Nessuna mappa z-score trovata in: {zmap_dir}')

print(f'Pazienti trovati: {N} | ROI: {R}')

roi_matrix = np.full((N, R), np.nan)   # NaN di default (FOV incompleto)
subj_ids   = []

for i, zpath in enumerate(zfiles):
    z_vol = nib.load(zpath).get_fdata()
    subj_ids.append(zpath.stem)        # es. 'z_ALS001'

    for r, (code, name) in enumerate(roi_table):
        mask_roi = (aal_vol == code)                        # maschera binaria della ROI
        voxels   = z_vol[mask_roi & np.isfinite(z_vol)]    # interseca con FOV (esclude NaN)

        if len(voxels) == 0:
            # FOV incompleto: la ROI non è disponibile per questo soggetto
            # roi_matrix[i, r] rimane NaN
            pass
        else:
            roi_matrix[i, r] = np.mean(voxels)

    # Log di avanzamento ogni 10 soggetti
    if (i + 1) % 10 == 0 or (i + 1) == N:
        n_nan = np.sum(np.isnan(roi_matrix[i]))
        print(f'  [{i+1}/{N}] {subj_ids[-1]} — ROI mancanti: {n_nan}/{R}')

# =====================================================================
# CONTROLLO QUALITÀ RAPIDO — soggetti con troppi NaN
# =====================================================================
nan_per_subj = np.sum(np.isnan(roi_matrix), axis=1)
soglia_esclusione = int(0.10 * R)   # > 10% di ROI mancanti
sogg_da_escludere = [subj_ids[i] for i in range(N) if nan_per_subj[i] > soglia_esclusione]

if sogg_da_escludere:
    print(f'\nATTENZIONE: {len(sogg_da_escludere)} soggetti con > 10% ROI mancanti:')
    for s in sogg_da_escludere:
        print(f'  {s}')
    print('Questi soggetti vanno esclusi o verificati manualmente.')

# =====================================================================
# EXPORT CSV
# =====================================================================
df_out = pd.DataFrame(roi_matrix, columns=roi_names)
df_out.insert(0, 'SubjectID', subj_ids)
df_out.to_csv(output_csv, index=False)
print(f'\nMatrice {N}×{R} salvata in: {output_csv}')
```

---

## 8.6 Utilizzo dell'atlante AAL3 in MRIcroGL (per le figure)

MRIcroGL supporta nativamente AAL3 tramite i file distribuiti con l'archivio. L'installazione consente di usare lo stesso atlante sia per l'estrazione quantitativa (Python) sia per la visualizzazione qualitativa e la produzione delle figure di pubblicazione, garantendo coerenza anatomica end-to-end.

**Installazione dei file atlante in MRIcroGL (Fedora Linux)**

Individua la cartella `templates` di MRIcroGL sul tuo sistema. Su Fedora, a seconda di come è stato installato, può trovarsi in `~/.local/share/mricron/templates/` oppure nella directory di installazione del programma (verificabile con `which MRIcroGL` e navigando nella cartella superiore). Una volta individuata:

```bash
# Copia i file AAL3 per MRIcroGL dalla cartella di download
cp ~/Scaricati/AAL3/AAL3v1.nii.gz    ~/.local/share/mricron/templates/
cp ~/Scaricati/AAL3/AAL3v1.nii.txt   ~/.local/share/mricron/templates/
```

Il file `.nii.gz` è il volume atlas compresso, il file `.txt` è la lista di etichette che MRIcroGL usa per mostrare i nomi delle regioni al click. Entrambi sono necessari.

**Dopo aver copiato i file**, riavvia MRIcroGL. L'atlante sarà disponibile nel menu `Atlases` (o `Atlas` a seconda della versione) con il nome `AAL3v1`.

**Utilizzo per le figure**

Il workflow tipico per una figura di pubblicazione è il seguente. Si carica prima il template anatomico di riferimento (ad esempio il template MNI152 incluso in MRIcroGL, oppure il template FDG-PET di Della Rosa se si vuole un background PET). Si sovrappone poi la mappa z-score media del gruppo (o la mappa dei risultati SuStaIn) come overlay colorato. Si attiva infine l'atlante AAL3 per verificare l'identificazione anatomica delle regioni: cliccando su un punto dell'immagine, MRIcroGL mostra il nome della regione AAL3 corrispondente. Questo permette di verificare che le regioni con gli eventi SuStaIn più precoci corrispondano effettivamente alle strutture elencate nella tabella del §7.4.

Per evidenziare specifiche ROI, è possibile generare maschere binarie da Python (estraendo i voxel con un codice AAL3 specifico dal volume `ROI_MNI_V7.nii`) e caricarle come overlay separati in MRIcroGL per colorare selettivamente le regioni di interesse. Questo è il metodo raccomandato per le figure che mostrano la localizzazione anatomica delle ROI incluse nell'analisi.

---

## 7.7 Controllo qualità della matrice $\mathbf{X}$

Prima di passare al modello SuStaIn è necessario ispezionare la struttura della matrice $\mathbf{X}$ sotto tre profili.

**Distribuzione per ROI.** Violin plot o boxplot degli z-score per ciascuna delle 24 regioni sull'intera coorte ALS. Ci si attende che la corteccia motoria primaria (Precentral) e la SMA siano tra le regioni con distribuzione più spostata verso valori positivi (ipometabolismo). Regioni come il cingolato posteriore e il temporale medio dovrebbero mostrare maggiore variabilità, riflettendo l'eterogeneità dei sottotipi. Una ROI con distribuzione centrata su zero suggerirebbe assenza di coinvolgimento sistematico in questa coorte.

**Heatmap $N \times R$.** La visualizzazione soggetti × regioni è lo strumento diagnostico più efficace: permette di identificare in un'unica immagine i soggetti *outlier* (righe con z-score uniformemente estremi, suggestive di fallimento del preprocessing individuale) e le ROI con dati mancanti sistematici (colonne con molti NaN, tipicamente le cerebellari inferiori in scansioni con FOV ridotto).

**Correlazione inter-ROI.** Regioni anatomicamente adiacenti o funzionalmente connesse — ad esempio Precentral_L e Precentral_R, oppure SMA e Precentral — dovrebbero mostrare correlazioni positive moderate (r = 0.3–0.7). Correlazioni negative inattese tra regioni che ci si aspetta co-varino positivamente sono un segnale di errore sistematico nel preprocessing di un sottogruppo di soggetti.

**Soglia di missing data.** Soggetti con più del 10% di ROI mancanti (ovvero più di 2–3 regioni NaN su 24) vengono esclusi dall'analisi principale. Questi corrispondono tipicamente a scansioni con FOV assiale < 13 cm che esclude completamente il cervelletto. Un'analisi di sensibilità senza le ROI cerebellari può essere condotta includendo anche questi soggetti.

---
---

[vai a  PARTE 2](https://github.com/fedele93/sustainpetsla/wiki/Funzionamento-di-SuStain)
---



---

# Riferimenti Metodologici Chiave

1. **Young AL et al. (2018)**. Uncovering the heterogeneity and temporal complexity of neurodegenerative diseases with Subtype and Stage Inference. *Nature Communications*, 9:4273.

2. **Cabras S et al. (2025)**. Role of 2-[18F]FDG-PET as a biomarker of upper motor neuron involvement in amyotrophic lateral sclerosis. *Journal of Neurology*, 272:766.

3. **Tan HHG et al. (2022)**. MRI Clustering Reveals Three ALS Subtypes With Unique Neurodegeneration Patterns. *Annals of Neurology*, 92:1030–1045.

4. **Della Rosa et al. (2014)**. A standardized [18F]-FDG-PET template for spatial normalization in statistical parametric mapping of dementia. *Neuroinformatics*, 12:575–593.

5. **Aksman LM et al. (2021)**. pySuStaIn: A Python Implementation of the Subtype and Stage Inference Algorithm. *SoftwareX*, 16:100811.

6. https://github.com/PasqualeDellaRosa/Dementia-Specific-18F-FDG-PET-template — Template FDG-PET (licenza MIT)

7. https://github.com/ucl-pond/pySuStaIn — Repository pySuStaIn (licenza MIT)

8. https://adni.loni.usc.edu/data-samples/adni-data/#AccessData — Portale download immagini ADNI
9. **Johnson WE, Li C, Rabinovic A (2007)**. Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics*, 8(1):118–127. (da controllare)
10. **Fortin JP et al. (2017)**. Harmonization of cortical thickness measurements across scanners and sites. *NeuroImage*, 167:104–120. (da controllare)
11. **Fortin JP et al. (2018)**. Harmonization of multi-site diffusion tensor imaging data. *NeuroImage*, 161:149–170. (da controllare)
12. **Pomponio R et al. (2020)**. Harmonization of large MRI datasets for the analysis of brain imaging patterns throughout the lifespan. *NeuroImage*, 208:116450. (da conrollare)

---

---

# Appendice A — Implementazione Computazionale

Tutti gli script sono scritti in MATLAB R2025b con SPM25 (Step 1–3) e Python 3 (script ausiliari). I percorsi nella sezione `CONFIGURAZIONE` di ciascuno script devono essere adattati al sistema locale prima dell'esecuzione. Gli script sono eseguibili sequenzialmente nell'ordine indicato.

Le immagini ADNI sono organizzate nella cartella `ADNI_CoregAvg_baseline/` con struttura `SUBJECT_ID/IMAGE_ID/*.dcm`, contenente esclusivamente le scansioni baseline selezionate tramite lo script `extract_baseline_dicoms.py` (Appendice A.0).

---

## A.0 — Selezione e Organizzazione delle Scansioni Baseline
**File**: `extract_baseline_dicoms.py`

Script Python che, a partire dal CSV baseline (`ADNI_FDG_CN_baseline_CoregAvg_withDates.csv`), estrae dalle cartelle ADNI scaricate le sole scansioni baseline per ciascun soggetto, copiandole in una struttura pulita e ordinata. Garantisce che la pipeline MATLAB elabori unicamente le immagini selezionate.

```python
"""
extract_baseline_dicoms.py

Estrae le cartelle DICOM delle sole scansioni baseline (prima acquisizione)
dalla directory ADNI scaricata, copiandole in una struttura pulita e ordinata.

Struttura output:
    output_dir/
    ├── SUBJECT_ID/
    │   └── IMAGE_ID/
    │       ├── ADNI_..._.dcm
    │       └── ...
    └── ...

Input richiesto:
    - CSV baseline con colonne: image_data_id, subject_id, acq_date
      (image_data_id include già il prefisso 'I', es. 'I240524')
    - Directory ADNI con struttura:
        adni_root/SUBJECT_ID/Co-registered__Averaged/DATE/IMAGE_ID/*.dcm
"""

import shutil
import pandas as pd
from pathlib import Path

# ============================================================
# CONFIGURAZIONE
# ============================================================
adni_root    = '/home/fluisi/Dati ricerca/ADNI'
baseline_csv = '/home/fluisi/Dati ricerca/ADNI_FDG_CN_baseline_CoregAvg_withDates.csv'
output_dir   = '/home/fluisi/Dati ricerca/ADNI_CoregAvg_baseline'
# ============================================================

def find_image_dir(subj_dir: Path, image_id: str) -> Path | None:
    coreg_candidates = [
        subj_dir / 'Co-registered__Averaged',
        subj_dir / 'Co-registered, Averaged',
    ]
    coreg_dir = None
    for candidate in coreg_candidates:
        if candidate.exists():
            coreg_dir = candidate
            break
    if coreg_dir is None:
        return None
    matches = list(coreg_dir.rglob(image_id))
    matches = [m for m in matches if m.is_dir()]
    if matches:
        return matches[0]
    return None

def main():
    df = pd.read_csv(baseline_csv)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    log_lines = [f'Estrazione baseline DICOM — {pd.Timestamp.now()}\n']
    n_ok = n_skip = n_missing = 0

    for _, row in df.iterrows():
        subj_id  = str(row['subject_id']).strip()
        image_id = str(row['image_data_id']).strip()
        acq_date = str(row['acq_date']).strip()
        dest_dir = output_path / subj_id / image_id

        if dest_dir.exists() and any(dest_dir.iterdir()):
            log_lines.append(f'SKIP: {subj_id} | {image_id}\n')
            n_skip += 1
            continue

        subj_dir = Path(adni_root) / subj_id
        if not subj_dir.exists():
            log_lines.append(f'MISSING_SUBJ: {subj_id}\n')
            n_missing += 1
            continue

        img_dir = find_image_dir(subj_dir, image_id)
        if img_dir is None:
            log_lines.append(f'MISSING_IMGID: {subj_id} | {image_id}\n')
            n_missing += 1
            continue

        dcm_files = list(img_dir.glob('*.dcm')) + list(img_dir.glob('*.DCM'))
        if not dcm_files:
            log_lines.append(f'NO_DICOM: {subj_id} | {image_id}\n')
            n_missing += 1
            continue

        dest_dir.mkdir(parents=True, exist_ok=True)
        for dcm in dcm_files:
            shutil.copy2(dcm, dest_dir / dcm.name)

        log_lines.append(f'OK: {subj_id} | {image_id} | {len(dcm_files)} file\n')
        n_ok += 1

    summary = f'\nSuccessi: {n_ok} | Skip: {n_skip} | Errori: {n_missing}\n'
    log_lines.append(summary)
    with open(output_path / 'extract_log.txt', 'w') as f:
        f.writelines(log_lines)

if __name__ == '__main__':
    main()
```

---

## A.1 — Preprocessing e Normalizzazione Spaziale
**File**: `pipeline_step1_CoregAvg_v2.m`

Corrisponde agli Step 2 e 3 del Methods. Per ogni soggetto: individua la cartella DICOM della scansione baseline tramite l'image_id, converte i DICOM in NIfTI, corregge l'origine dell'header, normalizza al template FDG-PET di Della Rosa et al. con Old Normalise, applica smoothing 10mm FWHM. Produce un file `s<SUBJECT_ID>_w.nii` per soggetto.

```matlab
% =========================================================================
% pipeline_step1_CoregAvg_v2.m
%
% Pipeline preprocessing FDG-PET ADNI — Controlli CN
% Formato input: Co-registered, Averaged (DICOM multi-slice)
%                struttura: ADNI_CoregAvg_baseline/SUBJECT_ID/IMAGE_ID/*.dcm
%
% Step eseguiti per ogni soggetto:
%   1. DICOM Import → unico volume NIfTI 3D
%   2. Reset origine al centro geometrico
%   3. Normalizzazione spaziale (Old Normalise) → template Della Rosa 2014
%   4. Smoothing isotropico 10mm FWHM
%
% Output: file s<SUBJECT_ID>_w.nii per ogni soggetto in output_dir
% =========================================================================

%% CONFIGURAZIONE
adni_root    = '/home/fluisi/Dati ricerca/ADNI_CoregAvg_baseline';
baseline_csv = '/home/fluisi/Dati ricerca/ADNI_FDG_CN_baseline_CoregAvg_withDates.csv';
output_dir   = '/home/fluisi/Dati ricerca/ADNI_nii_v3';
template     = '/home/fluisi/Dati ricerca/Template normalizzazione/TEMPLATE_FDGPET_100.nii';

%% INIZIALIZZAZIONE
spm('defaults', 'PET');
spm_jobman('initcfg');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% LEGGI CSV BASELINE
fid = fopen(baseline_csv, 'r');
header = fgetl(fid);
baseline = struct('image_data_id', {}, 'subject_id', {});
while ~feof(fid)
    line = fgetl(fid);
    if ischar(line) && ~isempty(strtrim(line))
        parts = strsplit(line, ',');
        if numel(parts) >= 2
            baseline(end+1).image_data_id = strtrim(parts{1});
            baseline(end).subject_id      = strtrim(parts{2});
        end
    end
end
fclose(fid);
fprintf('Soggetti nel CSV baseline: %d\n\n', numel(baseline));

%% LOG FILE
log_path = fullfile(output_dir, 'pipeline_log.txt');
log_fid  = fopen(log_path, 'w');
fprintf(log_fid, 'Pipeline log — %s\n\n', datestr(now));

%% LOOP SUI SOGGETTI
n_ok = 0;
for s = 1:numel(baseline)
    subj_id  = baseline(s).subject_id;
    image_id = baseline(s).image_data_id;
    fprintf('[%d/%d] Soggetto: %s | Image ID: %s\n', s, numel(baseline), subj_id, image_id);

    out_nii = fullfile(output_dir, sprintf('s%s_w.nii', subj_id));
    if exist(out_nii, 'file')
        fprintf('  -> Già processato, skip.\n');
        fprintf(log_fid, 'SKIP: %s\n', subj_id);
        continue;
    end

    dicom_dir = fullfile(adni_root, subj_id, image_id);
    if ~exist(dicom_dir, 'dir')
        warning('  Cartella DICOM non trovata: %s', dicom_dir);
        fprintf(log_fid, 'MISSING_DIR: %s | %s\n', subj_id, image_id);
        continue;
    end

    dcm_files = dir(fullfile(dicom_dir, '*.dcm'));
    if isempty(dcm_files), dcm_files = dir(fullfile(dicom_dir, '*.DCM')); end
    if isempty(dcm_files)
        warning('  Nessun DICOM in %s', dicom_dir);
        fprintf(log_fid, 'NO_DICOM: %s\n', subj_id);
        continue;
    end
    fprintf('  DICOM trovati: %d slice\n', numel(dcm_files));

    dicom_list = cell(numel(dcm_files), 1);
    for d = 1:numel(dcm_files)
        dicom_list{d} = fullfile(dicom_dir, dcm_files(d).name);
    end

    tmp_dir = fullfile(output_dir, ['tmp_' subj_id]);
    if ~exist(tmp_dir, 'dir'), mkdir(tmp_dir); end

    % STEP 1 — DICOM Import
    try
        matlabbatch = {};
        matlabbatch{1}.spm.util.dicom.data             = dicom_list;
        matlabbatch{1}.spm.util.dicom.root              = 'flat';
        matlabbatch{1}.spm.util.dicom.outdir            = {tmp_dir};
        matlabbatch{1}.spm.util.dicom.protfilter        = '.*';
        matlabbatch{1}.spm.util.dicom.convopts.format   = 'nii';
        matlabbatch{1}.spm.util.dicom.convopts.icedims  = 0;
        spm_jobman('run', matlabbatch);
    catch ME
        warning('  DICOM Import fallito: %s', ME.message);
        fprintf(log_fid, 'DICOM_FAIL: %s | %s\n', subj_id, ME.message);
        rmdir(tmp_dir, 's'); continue;
    end

    nii_files = dir(fullfile(tmp_dir, '*.nii'));
    if isempty(nii_files)
        fprintf(log_fid, 'NO_NII: %s\n', subj_id);
        rmdir(tmp_dir, 's'); continue;
    end
    nii_path = fullfile(tmp_dir, nii_files(1).name);

    % STEP 2 — Reset origine
    try
        V = spm_vol(nii_path);
        center_vox = (V.dim / 2)';
        M = V.mat;
        M(1:3, 4) = -M(1:3, 1:3) * center_vox;
        spm_get_space(nii_path, M);
    catch ME
        warning('  Reset origine fallito: %s', ME.message);
        fprintf(log_fid, 'ORIGIN_FAIL: %s | %s\n', subj_id, ME.message);
        rmdir(tmp_dir, 's'); continue;
    end

    % STEP 3 — Old Normalise
    try
        matlabbatch = {};
        matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.source   = {[nii_path ',1']};
        matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.wtsrc    = '';
        matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.resample = {[nii_path ',1']};
        matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.template = {[template ',1']};
        matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.weight   = '';
        matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smosrc   = 8;
        matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smoref   = 0;
        matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.regtype  = 'mni';
        matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.cutoff   = 25;
        matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.nits     = 16;
        matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.reg      = 1;
        matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
        matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.bb       = [-78 -112 -70; 78 76 85];
        matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.vox      = [2 2 2];
        matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.interp   = 1;
        matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.wrap     = [0 0 0];
        matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.prefix   = 'w';
        spm_jobman('run', matlabbatch);
    catch ME
        warning('  Old Normalise fallito: %s', ME.message);
        fprintf(log_fid, 'NORM_FAIL: %s | %s\n', subj_id, ME.message);
        rmdir(tmp_dir, 's'); continue;
    end

    w_files = dir(fullfile(tmp_dir, 'w*.nii'));
    if isempty(w_files)
        fprintf(log_fid, 'NO_WNII: %s\n', subj_id);
        rmdir(tmp_dir, 's'); continue;
    end
    w_path = fullfile(tmp_dir, w_files(1).name);

    % STEP 4 — Smoothing 10mm FWHM
    try
        matlabbatch = {};
        matlabbatch{1}.spm.spatial.smooth.data   = {[w_path ',1']};
        matlabbatch{1}.spm.spatial.smooth.fwhm   = [10 10 10];
        matlabbatch{1}.spm.spatial.smooth.dtype  = 0;
        matlabbatch{1}.spm.spatial.smooth.im     = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        spm_jobman('run', matlabbatch);
    catch ME
        warning('  Smoothing fallito: %s', ME.message);
        fprintf(log_fid, 'SMOOTH_FAIL: %s | %s\n', subj_id, ME.message);
        rmdir(tmp_dir, 's'); continue;
    end

    sw_files = dir(fullfile(tmp_dir, 'sw*.nii'));
    if isempty(sw_files)
        fprintf(log_fid, 'NO_SWNII: %s\n', subj_id);
        rmdir(tmp_dir, 's'); continue;
    end

    copyfile(fullfile(tmp_dir, sw_files(1).name), out_nii);
    rmdir(tmp_dir, 's');
    fprintf(log_fid, 'OK: %s\n', subj_id);
    fprintf('  Output: %s  [OK]\n\n', out_nii);
    n_ok = n_ok + 1;
end

fclose(log_fid);
fprintf('\nPipeline completata. Log: %s\n', log_path);
```

---

## A.2 — Intensity Normalization (Whole Brain Mean)
**File**: `pipeline_step2_intensity_norm_v1.m`

Corrisponde allo Step 4 del Methods. Per ogni file `s<SUBJECT_ID>_w.nii`, divide ogni voxel per la media del segnale intra-cerebrale (voxel con intensità > 0). Produce un file `i<SUBJECT_ID>_w.nii` nella cartella output separata.

```matlab
% =========================================================================
% pipeline_step2_intensity_norm_v1.m
%
% Intensity normalization individuale per ogni immagine FDG-PET.
% Per ogni soggetto: divide ogni voxel per la whole brain mean.
%
% Input:  file s*_w.nii in input_dir (normalizzati + smoothati)
% Output: file i*_w.nii in output_dir (cartella separata)
% =========================================================================

%% CONFIGURAZIONE
input_dir  = '/home/fluisi/Dati ricerca/ADNI_nii_v3';
output_dir = '/home/fluisi/Dati ricerca/ADNI_nii_v3_intensity';

if ~exist(output_dir, 'dir'), mkdir(output_dir); end

nii_files = dir(fullfile(input_dir, 's*_w.nii'));
fprintf('File trovati: %d\n\n', numel(nii_files));

log_path = fullfile(output_dir, 'intensity_norm_log.txt');
log_fid  = fopen(log_path, 'w');
fprintf(log_fid, 'Intensity normalization log — %s\n\n', datestr(now));

n_ok = 0; n_fail = 0; n_skip = 0;

for s = 1:numel(nii_files)
    fname   = nii_files(s).name;
    subj_id = regexprep(fname, '^s(.+)_w\.nii$', '$1');
    in_path  = fullfile(input_dir,  fname);
    out_path = fullfile(output_dir, ['i' subj_id '_w.nii']);

    fprintf('[%d/%d] %s\n', s, numel(nii_files), subj_id);

    if exist(out_path, 'file')
        fprintf('  -> Già presente, skip.\n');
        fprintf(log_fid, 'SKIP: %s\n', subj_id);
        n_skip = n_skip + 1; continue;
    end

    try
        V    = spm_vol(in_path);
        data = spm_read_vols(V);
        brain_voxels = data(data > 0);
        if isempty(brain_voxels)
            fprintf(log_fid, 'EMPTY: %s\n', subj_id);
            n_fail = n_fail + 1; continue;
        end
        wb_mean   = mean(brain_voxels);
        data_norm = data / wb_mean;
        V_out = V;
        V_out.fname   = out_path;
        V_out.descrip = sprintf('Intensity normalized (WB mean=%.4f)', wb_mean);
        spm_write_vol(V_out, data_norm);
        fprintf(log_fid, 'OK: %s | WB_mean=%.4f\n', subj_id, wb_mean);
        n_ok = n_ok + 1;
    catch ME
        fprintf(log_fid, 'FAIL: %s | %s\n', subj_id, ME.message);
        n_fail = n_fail + 1;
    end
end

fprintf(log_fid, '\nOK: %d | Skip: %d | Errori: %d\n', n_ok, n_skip, n_fail);
fclose(log_fid);
fprintf('Log: %s\n', log_path);
```

---

## A.3 — Modello Normativo Voxel-Wise
**File**: `pipeline_step3_normative_model_v2.m`

Corrisponde allo Step 5 del Methods. Carica tutti i 365 volumi `i*_w.nii` in memoria, stima per ogni voxel un modello OLS con età standardizzata e sesso centrato come covariate, e salva le mappe dei coefficienti e della deviazione standard dei residui.

```matlab
% =========================================================================
% pipeline_step3_normative_model_v2.m
%
% Modello normativo OLS voxel-wise sui controlli CN ADNI.
% Per ogni voxel: uptake ~ β0 + β1*età_z + β2*sesso_c
%
% Output in output_dir:
%   normative_beta.nii   — coefficienti (4D: β0, β1, β2)
%   normative_sigma.nii  — SD dei residui
%   normative_mean.nii   — metabolismo predetto alla media (β0)
%   normative_params.mat — parametri per applicazione ai pazienti
% =========================================================================

%% CONFIGURAZIONE
input_dir    = '/home/fluisi/Dati ricerca/ADNI_nii_v3_intensity';
baseline_csv = '/home/fluisi/Dati ricerca/ADNI_FDG_CN_baseline_CoregAvg_withDates.csv';
output_dir   = '/home/fluisi/Dati ricerca/normative_model_v2';

if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% LEGGI CSV E ABBINA FILE
opts = detectImportOptions(baseline_csv);
opts = setvartype(opts, 'sex', 'char');
T = readtable(baseline_csv, opts);

subj_list = {}; age_list = []; sex_list = [];
for i = 1:height(T)
    subj_id = strtrim(T.subject_id{i});
    fpath   = fullfile(input_dir, sprintf('i%s_w.nii', subj_id));
    if ~exist(fpath, 'file'), continue; end
    sex_str = strtrim(T.sex{i});
    if strcmpi(sex_str, 'M'), sex_val = 1;
    elseif strcmpi(sex_str, 'F'), sex_val = 0;
    else, continue; end
    subj_list{end+1} = fpath;
    age_list(end+1)  = T.age_at_scan(i);
    sex_list(end+1)  = sex_val;
end
N = numel(subj_list);
fprintf('Soggetti abbinati: %d\n', N);

%% STANDARDIZZA COVARIATE
age_mean = mean(age_list); age_std = std(age_list);
age_z    = (age_list - age_mean) / age_std;
sex_mean = mean(sex_list);
sex_c    = sex_list - sex_mean;
fprintf('Età: media=%.1f, SD=%.1f | Sesso: %.1f%% M\n\n', age_mean, age_std, sex_mean*100);

%% MATRICE DEL DESIGN
X = [ones(N,1), age_z(:), sex_c(:)];
XtXinv_Xt = (X' * X) \ X';

%% CARICA VOLUMI
V_ref = spm_vol(subj_list{1}); dim = V_ref.dim; n_vox = prod(dim);
Y = zeros(n_vox, N, 'single');
for s = 1:N
    V = spm_vol(subj_list{s});
    data = spm_read_vols(V);
    Y(:,s) = single(data(:));
    if mod(s,50)==0, fprintf('  Caricati %d/%d\n', s, N); end
end

%% MASCHERA E STIMA OLS
mean_Y = mean(Y, 2);
mask   = mean_Y > 0.1;
fprintf('Voxel nella maschera: %d (%.1f%%)\n', sum(mask), 100*sum(mask)/n_vox);

beta_map  = zeros(n_vox, 3, 'single');
sigma_map = zeros(n_vox, 1, 'single');
Y_masked  = Y(mask, :);
beta_masked = (XtXinv_Xt * Y_masked')';
Y_hat       = (X * beta_masked')';
residui     = Y_masked - Y_hat;
df          = N - 3;
sigma_masked = sqrt(sum(residui.^2, 2) / df);
beta_map(mask, :) = beta_masked;
sigma_map(mask)   = sigma_masked;

%% SALVA OUTPUT
V_out = V_ref; V_out.dt = [16 0];
beta_names = {'intercetta', 'beta_eta', 'beta_sesso'};
for b = 1:3
    V_out.fname = fullfile(output_dir, 'normative_beta.nii');
    V_out.n = [b 1];
    V_out.descrip = sprintf('Normative model: %s', beta_names{b});
    spm_write_vol(V_out, reshape(beta_map(:,b), dim));
end

V_out.n = [1 1];
V_out.fname = fullfile(output_dir, 'normative_sigma.nii');
V_out.descrip = 'Normative model: sigma';
spm_write_vol(V_out, reshape(sigma_map, dim));

V_out.fname = fullfile(output_dir, 'normative_mean.nii');
V_out.descrip = 'Normative model: mean predicted';
spm_write_vol(V_out, reshape(beta_map(:,1), dim));

save(fullfile(output_dir, 'normative_params.mat'), ...
    'age_mean', 'age_std', 'sex_mean', 'dim', 'mask', ...
    'subj_list', 'age_list', 'sex_list', 'N', 'df');

fprintf('\n=== MODELLO NORMATIVO COMPLETATO ===\n');
fprintf('N=%d | df=%d | Voxel stimati=%d\n', N, df, sum(mask));
```

---


