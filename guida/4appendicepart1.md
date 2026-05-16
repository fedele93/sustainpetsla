# Parte 1 — Appendice: Script Completi

Tutti gli script sono numerati nell'ordine di esecuzione. I percorsi nella sezione `CONFIGURAZIONE` di ciascuno vanno adattati al sistema locale.

---

## A.0 — Selezione e Organizzazione delle Scansioni Baseline
**File**: `extract_baseline_dicoms.py`

Estrae le cartelle DICOM delle sole scansioni baseline da una directory ADNI scaricata, copiandole in una struttura pulita `SUBJECT_ID/IMAGE_ID/*.dcm`.

```python
"""
extract_baseline_dicoms.py
Struttura output: output_dir/SUBJECT_ID/IMAGE_ID/*.dcm
Input CSV: image_data_id, subject_id, acq_date
"""
import shutil
import pandas as pd
from pathlib import Path

adni_root    = '/home/fluisi/Dati ricerca/ADNI'
baseline_csv = '/home/fluisi/Dati ricerca/ADNI_FDG_CN_baseline_CoregAvg_withDates.csv'
output_dir   = '/home/fluisi/Dati ricerca/ADNI_CoregAvg_baseline'

def find_image_dir(subj_dir, image_id):
    for candidate in [subj_dir / 'Co-registered__Averaged', subj_dir / 'Co-registered, Averaged']:
        if candidate.exists():
            matches = [m for m in candidate.rglob(image_id) if m.is_dir()]
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
        dest_dir = output_path / subj_id / image_id

        if dest_dir.exists() and any(dest_dir.iterdir()):
            log_lines.append(f'SKIP: {subj_id} | {image_id}\n'); n_skip += 1; continue

        subj_dir = Path(adni_root) / subj_id
        if not subj_dir.exists():
            log_lines.append(f'MISSING_SUBJ: {subj_id}\n'); n_missing += 1; continue

        img_dir = find_image_dir(subj_dir, image_id)
        if img_dir is None:
            log_lines.append(f'MISSING_IMGID: {subj_id} | {image_id}\n'); n_missing += 1; continue

        dcm_files = list(img_dir.glob('*.dcm')) + list(img_dir.glob('*.DCM'))
        if not dcm_files:
            log_lines.append(f'NO_DICOM: {subj_id} | {image_id}\n'); n_missing += 1; continue

        dest_dir.mkdir(parents=True, exist_ok=True)
        for dcm in dcm_files:
            shutil.copy2(dcm, dest_dir / dcm.name)
        log_lines.append(f'OK: {subj_id} | {image_id} | {len(dcm_files)} file\n')
        n_ok += 1

    log_lines.append(f'\nSuccessi: {n_ok} | Skip: {n_skip} | Errori: {n_missing}\n')
    with open(output_path / 'extract_log.txt', 'w') as f:
        f.writelines(log_lines)

if __name__ == '__main__':
    main()
```

---

## A.1 — Preprocessing e Normalizzazione Spaziale
**File**: `pipeline_step1_CoregAvg_v2.m`

Corrisponde agli Step 2 e 3. Per ogni soggetto: DICOM Import, reset origine, Old Normalise al template Della Rosa, smoothing 10 mm FWHM. Produce `s<SUBJECT_ID>_w.nii`.

```matlab
% =========================================================================
% pipeline_step1_CoregAvg_v2.m
% Step eseguiti: DICOM Import → reset origine → Old Normalise → Smoothing
% Output: s<SUBJECT_ID>_w.nii
% =========================================================================

%% CONFIGURAZIONE
adni_root    = '/home/fluisi/Dati ricerca/ADNI_CoregAvg_baseline';
baseline_csv = '/home/fluisi/Dati ricerca/ADNI_FDG_CN_baseline_CoregAvg_withDates.csv';
output_dir   = '/home/fluisi/Dati ricerca/ADNI_nii_v3';
template     = '/home/fluisi/Dati ricerca/Template normalizzazione/TEMPLATE_FDGPET_100.nii';

spm('defaults', 'PET');
spm_jobman('initcfg');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% LEGGI CSV BASELINE
fid = fopen(baseline_csv, 'r'); fgetl(fid);
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

log_fid = fopen(fullfile(output_dir, 'pipeline_log.txt'), 'w');
fprintf(log_fid, 'Pipeline log — %s\n\n', datestr(now));

%% LOOP
for s = 1:numel(baseline)
    subj_id  = baseline(s).subject_id;
    image_id = baseline(s).image_data_id;
    out_nii  = fullfile(output_dir, sprintf('s%s_w.nii', subj_id));

    if exist(out_nii, 'file')
        fprintf(log_fid, 'SKIP: %s\n', subj_id); continue
    end

    dicom_dir = fullfile(adni_root, subj_id, image_id);
    if ~exist(dicom_dir, 'dir')
        fprintf(log_fid, 'MISSING_DIR: %s\n', subj_id); continue
    end

    dcm_files = dir(fullfile(dicom_dir, '*.dcm'));
    if isempty(dcm_files), dcm_files = dir(fullfile(dicom_dir, '*.DCM')); end
    if isempty(dcm_files)
        fprintf(log_fid, 'NO_DICOM: %s\n', subj_id); continue
    end

    dicom_list = cell(numel(dcm_files), 1);
    for d = 1:numel(dcm_files)
        dicom_list{d} = fullfile(dicom_dir, dcm_files(d).name);
    end

    tmp_dir = fullfile(output_dir, ['tmp_' subj_id]);
    if ~exist(tmp_dir, 'dir'), mkdir(tmp_dir); end

    try
        % DICOM Import
        matlabbatch = {};
        matlabbatch{1}.spm.util.dicom.data             = dicom_list;
        matlabbatch{1}.spm.util.dicom.root              = 'flat';
        matlabbatch{1}.spm.util.dicom.outdir            = {tmp_dir};
        matlabbatch{1}.spm.util.dicom.protfilter        = '.*';
        matlabbatch{1}.spm.util.dicom.ignorefilters     = 0;
        matlabbatch{1}.spm.util.dicom.icedims           = 0;
        matlabbatch{1}.spm.util.dicom.convopts.format   = 'nii';
        matlabbatch{1}.spm.util.dicom.convopts.meta     = 0;
        matlabbatch{1}.spm.util.dicom.convopts.rot90    = 0;
        spm_jobman('run', matlabbatch);

        nii_files = dir(fullfile(tmp_dir, '*.nii'));
        if isempty(nii_files), error('Nessun NIfTI prodotto'), end
        nii_path = fullfile(tmp_dir, nii_files(1).name);

        % Reset origine
        V = spm_vol(nii_path);
        center_vox = V.dim / 2;
        T = V.mat;
        T(1:3, 4) = -T(1:3, 1:3) * center_vox(:);
        spm_get_space(nii_path, T);

        % Old Normalise
        matlabbatch = {};
        matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.source = {[nii_path ',1']};
        matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.wtsrc  = '';
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

        w_nii = fullfile(tmp_dir, ['w' nii_files(1).name]);
        if ~exist(w_nii, 'file'), error('File warpato non trovato'), end

        % Smoothing 10 mm FWHM
        matlabbatch = {};
        matlabbatch{1}.spm.spatial.smooth.data   = {[w_nii ',1']};
        matlabbatch{1}.spm.spatial.smooth.fwhm   = [10 10 10];
        matlabbatch{1}.spm.spatial.smooth.dtype  = 0;
        matlabbatch{1}.spm.spatial.smooth.im     = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        spm_jobman('run', matlabbatch);

        sw_nii = fullfile(tmp_dir, ['sw' nii_files(1).name]);
        movefile(sw_nii, out_nii);
        rmdir(tmp_dir, 's');
        fprintf(log_fid, 'OK: %s\n', subj_id);

    catch ME
        fprintf(log_fid, 'FAIL: %s | %s\n', subj_id, ME.message);
        if exist(tmp_dir, 'dir'), rmdir(tmp_dir, 's'); end
    end
end
fclose(log_fid);
```

---

## A.1b — Conversione Analyze → NIfTI per Pazienti ALS

I file dei centri italiani in formato Analyze 7.5 vanno convertiti in NIfTI prima di procedere con la pipeline. Dopo la conversione è necessaria una verifica visiva per escludere flip sull'asse x.

```python
"""
convert_analyze_to_nifti.py
Converte file Analyze (.hdr/.img) in NIfTI (.nii) per i pazienti ALS.
Dopo la conversione: verificare visivamente l'orientamento in MRIcroGL.
"""
import nibabel as nib
import numpy as np
from pathlib import Path

input_dir  = Path('/path/to/ALS_analyze')
output_dir = Path('/path/to/ALS_nifti')
output_dir.mkdir(exist_ok=True)

for hdr_file in sorted(input_dir.glob('*.hdr')):
    img = nib.load(hdr_file)
    data = img.get_fdata()

    # Verifica flip: il valore del voxel (0,0,0) deve corrispondere all'angolo
    # posteriore-inferiore-destro in convenzione radiologica.
    # Se il cervello appare flippato sull'asse x in MRIcroGL, decommenta:
    # data = data[::-1, :, :]

    nii_out = nib.Nifti1Image(data.astype(np.float32), img.affine)
    out_path = output_dir / (hdr_file.stem + '.nii')
    nib.save(nii_out, out_path)
    print(f'Convertito: {out_path.name}')
```

Dopo la conversione, procedere con la **correzione dell'origine** (Step 2.2) usando il codice MATLAB già illustrato, poi normalizzazione spaziale, smoothing e intensity normalization con gli stessi parametri dei controlli ADNI.

---

## A.2 — Intensity Normalization
**File**: `pipeline_step2_intensity_norm_v1.m`

Corrisponde allo Step 4. Divide ogni voxel per la whole brain mean. Produce `i<SUBJECT_ID>_w.nii`.

```matlab
% =========================================================================
% pipeline_step2_intensity_norm_v1.m
% Input:  s*_w.nii in input_dir
% Output: i*_w.nii in output_dir
% =========================================================================

input_dir  = '/home/fluisi/Dati ricerca/ADNI_nii_v3';
output_dir = '/home/fluisi/Dati ricerca/ADNI_nii_v3_intensity';

if ~exist(output_dir, 'dir'), mkdir(output_dir); end
nii_files = dir(fullfile(input_dir, 's*_w.nii'));
fprintf('File trovati: %d\n\n', numel(nii_files));

log_fid = fopen(fullfile(output_dir, 'intensity_norm_log.txt'), 'w');
fprintf(log_fid, 'Intensity normalization log — %s\n\n', datestr(now));
n_ok = 0; n_fail = 0; n_skip = 0;

for s = 1:numel(nii_files)
    fname    = nii_files(s).name;
    subj_id  = regexprep(fname, '^s(.+)_w\.nii$', '$1');
    in_path  = fullfile(input_dir,  fname);
    out_path = fullfile(output_dir, ['i' subj_id '_w.nii']);

    if exist(out_path, 'file')
        fprintf(log_fid, 'SKIP: %s\n', subj_id); n_skip = n_skip + 1; continue
    end

    try
        V    = spm_vol(in_path);
        data = spm_read_vols(V);
        brain_voxels = data(data > 0);
        if isempty(brain_voxels)
            fprintf(log_fid, 'EMPTY: %s\n', subj_id); n_fail = n_fail + 1; continue
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
```

---

## A.3 — Percorso Visivo: Modello Normativo Voxel-Wise
**File**: `pipeline_step3_normative_model_v2.m`

> **Nota**: questo script fa parte del percorso visivo secondario (Step V) e **non alimenta SuStaIn**. Produce le mappe `normative_beta.nii` e `normative_sigma.nii` necessarie per calcolare le mappe z-score voxel-wise individuali (A.4). I parametri demografici di standardizzazione sono identici a quelli del percorso primario, ma l'armonizzazione ComBat **non** viene applicata ai volumi voxel-wise.

```matlab
% =========================================================================
% pipeline_step3_normative_model_v2.m
% Modello normativo OLS voxel-wise — percorso visivo (non alimenta SuStaIn)
% Output: normative_beta.nii, normative_sigma.nii, normative_mean.nii,
%         normative_params.mat
% =========================================================================

input_dir    = '/home/fluisi/Dati ricerca/ADNI_nii_v3_intensity';
baseline_csv = '/home/fluisi/Dati ricerca/ADNI_FDG_CN_baseline_CoregAvg_withDates.csv';
output_dir   = '/home/fluisi/Dati ricerca/normative_model_v2';

if ~exist(output_dir, 'dir'), mkdir(output_dir); end

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

age_mean = mean(age_list); age_std = std(age_list);
age_z    = (age_list - age_mean) / age_std;
sex_mean = mean(sex_list);
sex_c    = sex_list - sex_mean;

X = [ones(N,1), age_z(:), sex_c(:)];
XtXinv_Xt = (X' * X) \ X';

V_ref = spm_vol(subj_list{1}); dim = V_ref.dim; n_vox = prod(dim);
Y = zeros(n_vox, N, 'single');
for s = 1:N
    V = spm_vol(subj_list{s});
    data = spm_read_vols(V);
    Y(:,s) = single(data(:));
    if mod(s,50)==0, fprintf('  Caricati %d/%d\n', s, N); end
end

mean_Y = mean(Y, 2);
mask   = mean_Y > 0.1;
fprintf('Voxel nella maschera: %d (%.1f%%)\n', sum(mask), 100*sum(mask)/n_vox);

beta_map  = zeros(n_vox, 3, 'single');
sigma_map = zeros(n_vox, 1, 'single');
Y_masked  = Y(mask, :);
beta_masked  = (XtXinv_Xt * Y_masked')';
Y_hat        = (X * beta_masked')';
residui      = Y_masked - Y_hat;
sigma_masked = sqrt(sum(residui.^2, 2) / (N - 3));
beta_map(mask, :) = beta_masked;
sigma_map(mask)   = sigma_masked;

V_out = V_ref; V_out.dt = [16 0];
for b = 1:3
    V_out.fname   = fullfile(output_dir, 'normative_beta.nii');
    V_out.n       = [b 1];
    V_out.descrip = sprintf('Normative model: beta%d', b-1);
    spm_write_vol(V_out, reshape(beta_map(:,b), dim));
end
V_out.n = [1 1];
V_out.fname   = fullfile(output_dir, 'normative_sigma.nii');
V_out.descrip = 'Normative model: sigma';
spm_write_vol(V_out, reshape(sigma_map, dim));
V_out.fname   = fullfile(output_dir, 'normative_mean.nii');
V_out.descrip = 'Normative model: mean predicted';
spm_write_vol(V_out, reshape(beta_map(:,1), dim));

save(fullfile(output_dir, 'normative_params.mat'), ...
    'age_mean', 'age_std', 'sex_mean', 'dim', 'mask', 'N');
fprintf('\nModello normativo completato. N=%d | df=%d\n', N, N-3);
```

---

## A.4 — Percorso Visivo: Z-Score Voxel-Wise

> **Nota**: questo script fa parte del percorso visivo secondario e **non alimenta SuStaIn**. Produce le mappe `z_<SUBJECT_ID>.nii` per visualizzazione su MRIcroGL e analisi SPM second-level.

```python
"""
zscore_voxelwise.py — Percorso visivo (non alimenta SuStaIn)
Produce mappe z-score voxel-wise per ogni paziente ALS.
"""
import numpy as np
import nibabel as nib
import pandas as pd
from pathlib import Path

# Parametri del modello normativo voxel-wise (output A.3)
params     = np.load('/home/fluisi/Dati ricerca/normative_model_v2/normative_params.npz',
                     allow_pickle=True)
age_mean   = float(params['age_mean'])
age_std    = float(params['age_std'])
sex_mean   = float(params['sex_mean'])

beta_img   = nib.load('/home/fluisi/Dati ricerca/normative_model_v2/normative_beta.nii')
sigma_img  = nib.load('/home/fluisi/Dati ricerca/normative_model_v2/normative_sigma.nii')
beta_vol   = beta_img.get_fdata()   # shape: (X, Y, Z, 3)
sigma_vol  = sigma_img.get_fdata()
affine     = beta_img.affine
mask       = sigma_vol > 0

df_als     = pd.read_csv('als_subjects.csv')
# Colonne attese: subject_id, age_at_scan, sex, filepath_nii

output_dir = Path('/home/fluisi/Dati ricerca/zmaps_als')
output_dir.mkdir(exist_ok=True)

for _, row in df_als.iterrows():
    age_z_i  = (row['age_at_scan'] - age_mean) / age_std
    sex_c_i  = (1.0 if row['sex'] == 'M' else 0.0) - sex_mean

    pet_img  = nib.load(row['filepath_nii'])
    x_vol    = pet_img.get_fdata()

    pred_vol = (beta_vol[..., 0]
              + beta_vol[..., 1] * age_z_i
              + beta_vol[..., 2] * sex_c_i)

    z_vol = np.full(x_vol.shape, np.nan)
    z_vol[mask] = (x_vol[mask] - pred_vol[mask]) / sigma_vol[mask]
    z_vol[mask] = -z_vol[mask]   # positivo = ipometabolismo

    out_path = output_dir / f"z_{row['subject_id']}.nii"
    nib.save(nib.Nifti1Image(z_vol.astype(np.float32), affine), out_path)

print(f'Z-score voxel-wise calcolati per {len(df_als)} pazienti.')
print(f'Output: {output_dir}')
print('Uso: visualizzazione in MRIcroGL, figure di pubblicazione.')
print('NOTA: queste mappe NON alimentano SuStaIn.')
```

**Controllo qualità delle mappe z-score voxel-wise**: in un soggetto senza grave patologia diffusa, la distribuzione degli z-score cerebrali deve essere approssimativamente gaussiana centrata intorno a 0, con coda positiva nelle regioni ipometaboliche. Media > 1.5 su tutto il cervello suggerisce un problema di preprocessing. Sovrapporre visivamente alcune mappe sul template MNI per verificare che le regioni di ipometabolismo atteso (corteccia motoria, aree frontali) siano anatomicamente plausibili.
