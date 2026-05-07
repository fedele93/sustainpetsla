# FDG-PET SuStaIn nella SLA — Pipeline di analisi

> **Applicazione dell'algoritmo Subtype and Stage Inference (SuStaIn) all'imaging metabolico FDG-PET nella Sclerosi Laterale Amiotrofica**

---

## Panoramica

Questa repository contiene la pipeline di analisi sviluppata nell'ambito di un progetto di tesi di specializzazione che studia i sottotipi metabolici nella SLA tramite neuroimaging FDG-PET. La domanda centrale è se esistano pattern distinti di ipometabolismo cerebrale tra i pazienti — e se questi pattern si raggruppino in sottotipi biologicamente significativi che correlano con il background genetico, il fenotipo clinico e la prognosi.

L'approccio applica l'**algoritmo SuStaIn** (Subtype and Stage Inference; Young et al., *Nature Communications* 2018) a mappe z-score voxel-wise derivate da dati FDG-PET. SuStaIn è particolarmente adatto alla ricerca sulle malattie neurodegenerative perché inferisce simultaneamente i sottotipi di malattia e colloca ciascun paziente lungo una traiettoria di progressione specifica per sottotipo, senza richiedere etichette diagnostiche a priori.

La coorte SLA comprende forme sporadiche e genetiche (C9orf72, SOD1, TARDBP, FUS e mutazioni più rare) con dati clinici associati (ALSFRS-R, sede di esordio, durata di malattia, ECAS, neurofilamenti, ptau217/181). Soggetti cognitivamente normali del database **ADNI** costituiscono la popolazione di controllo normativa.

---

## Contesto scientifico

La FDG-PET misura il metabolismo glucidico cerebrale regionale. Nella SLA, la neurodegenerazione causa **ipometabolismo** nelle regioni coinvolte — in modo più evidente nella corteccia motoria primaria, nell'area motoria supplementare e nelle regioni frontali — sebbene il pattern spaziale vari considerevolmente tra pazienti e background genetici.

L'ipotesi centrale è che questa eterogeneità non sia rumore casuale, ma rifletta sottotipi biologici distinti. SuStaIn la modella apprendendo un insieme di **sequenze z-score**: progressioni ordinate di ipometabolismo regionale che definiscono ciascun sottotipo. I pazienti vengono poi assegnati probabilisticamente a un sottotipo e posizionati lungo il suo asse di staging.

Le analisi sono pianificate a tre livelli:

**Coorte SLA completa** — sottotipizzazione esplorativa (fino a 5 sottotipi) per caratterizzare l'eterogeneità globale.

**C9orf72 isolato** — il sottogruppo genetico più numeroso (~70 soggetti), analizzato separatamente data la sua ampia variabilità fenotipica, dalla forma motoria pura alla ALS-FTD.

**SOD1 come controllo biologico** — la SLA-SOD1 presenta un coinvolgimento prevalentemente motorio in assenza di patologia TDP-43; funge da validazione interna: il modello dovrebbe recuperare un pattern circoscritto alle regioni motorie. (dubbio, potrei decidere di inserirli nella coorte completa)



---

## Struttura della pipeline

La pipeline è divisa in due parti sequenziali. La **Parte 1** (implementata in MATLAB/SPM) produce le mappe z-score per soggetto e i valori ROI a partire dalle immagini PET grezze. La **Parte 2** (implementata in Python) applica SuStaIn a quei valori ROI ed esegue tutte le analisi a valle.

### Parte 1 — Preprocessing e calcolo degli z-score (MATLAB / SPM25)

| Step | Descrizione |
|------|-------------|
| 1 | Selezione dei controlli cognitivamente normali ADNI (scansione baseline per soggetto) |
| 2 | Conversione DICOM→NIfTI, correzione dell'origine dell'header, normalizzazione spaziale al template FDG-PET di Della Rosa et al. (2014) tramite SPM Old Normalise |
| 3 | Smoothing spaziale — kernel gaussiano isotropico 10 mm FWHM |
| 4 | Normalizzazione in intensità — divisione per la media dell'intero cervello individuale |
| 5 | Modello normativo OLS voxel-wise stimato sui controlli ADNI (covariate: età, sesso) |
| 6 | Mappe z-score per paziente — deviazione dal modello normativo in ogni voxel |
| 7 | Estrazione dei valori ROI con atlante AAL tramite MarsBaR (~18–20 regioni) |

Gli script MATLAB che implementano gli Step 1–5 sono documentati nella **Wiki** e riprodotti integralmente nell'Appendice A della stessa .

### Parte 2 — Modellazione SuStaIn e analisi (Python / pySuStaIn)

| Analisi | Descrizione |
|---------|-------------|
| 2a | SuStaIn sull'intero dataset SLA (Cmax = 5) |
| 2b | SuStaIn sul sottogruppo C9orf72 isolato (Cmax = 3) |
| 2c | SuStaIn sugli altri sottogruppi genetici |
| 2d | SOD1 — validazione biologica |
| 2e | Correlazione tra sottotipi/stadi e variabili cliniche e biomarcatori |
| 2f | Analisi prognostica per sottotipo |

---

## Contenuto della repository

```
├── Vesione_0_3_SuStaIn_ALS_FDG_template.ipynb   ← Notebook principale di analisi (Parte 2)
└── README.md                                      ← Questo file
```

La **Wiki** (scheda Wiki in alto) contiene:

*Preparazione dei dati* — descrizione completa della pipeline di preprocessing (Parte 1), con la normalizzazione spaziale, lo smoothing, la normalizzazione in intensità e la costruzione del modello normativo, inclusi gli script MATLAB completi.

*Passaggi SuStaIn* — guida dettagliata a ciascuno step della Parte 2, con i controlli di convergenza MCMC, la selezione del numero di sottotipi tramite CVIC, i Positional Variance Diagrams e l'output dello staging.

---

## Notebook principale: `Vesione_0_3_SuStaIn_ALS_FDG_template.ipynb`

È il notebook centrale del progetto. Attualmente gira su **dati sintetici** con una struttura progettata per rispecchiare il dataset SLA reale — in modo da poter sviluppare, testare e revisionare l'intera pipeline prima che i dati FDG-PET reali siano disponibili. Quando i valori z-score ROI reali saranno pronti, sarà sufficiente sostituire la cella di caricamento dati; il resto del notebook non cambia.

Il notebook copre: esplorazione del dataset e distribuzioni z-score per ROI; configurazione e fitting del modello SuStaIn; diagnostica di convergenza MCMC; selezione del numero di sottotipi (CVIC e confronto tra modelli); Positional Variance Diagrams come visualizzazione principale dell'output SuStaIn; staging e assegnazione al sottotipo per ogni paziente; correlazioni cliniche e con biomarcatori post-hoc; analisi C9orf72 isolata; validazione biologica SOD1.

---

## Dipendenze

L'ambiente Python richiede:

```bash
pip install pySuStaIn numpy pandas matplotlib seaborn scipy scikit-learn
```

Per la pipeline di preprocessing MATLAB sono necessari **MATLAB R2024+** con **SPM12 o SPM25** e il [template FDG-PET di Della Rosa et al.](https://github.com/PasqualeDellaRosa/Dementia-Specific-18F-FDG-PET-template).

---

## Riferimenti principali

- Young AL, et al. (2018). Uncovering the heterogeneity and temporal complexity of neurodegenerative diseases with Subtype and Stage Inference. *Nature Communications*, 9, 4273.
- Della Rosa PA, et al. (2014). A standardized [18F]-FDG-PET template for spatial normalization in statistical parametric mapping of dementia. *Neuroinformatics*, 12(4), 575–593.


---

## Stato del progetto

Il progetto è in sviluppo attivo. La Parte 1 (preprocessing fino al modello normativo) è largamente implementata. Il notebook principale gira dall'inizio alla fine sui dati sintetici. L'estrazione delle ROI (Step 7) e l'applicazione della pipeline z-score alla coorte SLA reale sono in corso.

*Progetto di tesi — Università di Torino/Bari*
