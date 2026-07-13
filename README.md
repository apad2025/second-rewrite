# Changes
- If on Windows, manually add each to the MATLAB path
    - GRMD
    - hernando
    - functions
    - bipolar_fat_water_separation
- Change line 1559 `mainpth` variable in `dogexplorer.m`
- Rerun `dogexplorer.m` to regenerage `DD.mat`
- Change line 8 in `processeddogdata.m`
- In `DogAnalysis.m`:
    - Set `pth_data` to the same path as `mainpth` in `dogexplorer.m`
    - Set `pth_code` to the GRMD folder

Comment out the following in `DogAnalysis.m`:
- `pth_main` since it is unused
- **Bipolar correction** lines 269-315
- **Unwrap data** lines 316-357
- Lines 376-406
    - `dfat = -3.5*D.F0;` to `D.CSCorrection.FatWaterSwapChecked.Automatic = true;`
- **Check for fat/water swaps - automatic** lines 411-439
- **Check for fat/water swaps - manual** lines 440-575
- Lines 578 and 579:
```matlab
fields2rm = {'Image', 'WeightedMagnitude'};
if isfield(D.Data, 'Image'), D.Data = rmfield(D.Data, 'Image'); end
```

# Linux only changes
## dogexplorer.m
Change line 1638 from
```matlab
save(join([mainpth, "DD.mat"], '\'), "DD")
```
to
```matlab
save(fullfile(mainpth, "DD.mat"), "DD")
```

All other lines:

Preprocess / spectra section
- 233 — _noiseCov.mat save
- 258, 259 — path_save_raw, path_save_preproc_twix

Fit section
- 422, 424, 426 — _preproc.mat, _preproc_phased.mat, _preproc_twix.mat
- 442 — _proc.mat
- 459 — _OXSA_opts.mat
- 536 — _OXSA_results.mat

Analyze section
- 556, 557, 558 — path_proc, path_proc_twix, path_oxsa
- 569 — load(...) twix

FixPPM0
- 1155, 1171 — the two Ps2change literal-backslash concatenations

GenPath
- 1291, 1292 — GRE raw/save
- 1296, 1301, 1306, 1311 — the four identical Path.Main + Hydrogen.Main lines (done via replace_all)
- 1297, 1302, 1307, 1312 — Images.Main / Noise / Unsuppressed / Suppressed
- 1316 — Phosphorus

parseInputs
- 1413, 1464 — the two path_tmp lookups

Initialize
- 1614 — Path.Main (date, Subject)
- 1616, 1618 — Path.DICOM (−1 / −2 suffix)
- 1621, 1622 — Path.Main / Path.DICOM (single-date)

## DogInitialize.m
Change line 31
```matlab
load(append(pth_DD,"\DD"), 'DD')
```
to
```matlab
load(fullfile(pth_DD,"DD"), 'DD')
```

## DogAnalysis.m
Comment out the unused Windows-path variables (lines 6-7). They are dead code — only `pth_code`/`pth_data` are used.

Reconcile the `DD.mat` location so it does not need to be copied by hand. `dogexplorer.m` saves `DD.mat` to `mainpth` (the `DICOM_Files` root), but `DogInitialize` loads it from whatever directory it is passed. Point that at the same root.

Change lines 6-8 from
```matlab
pth = "C:\Users\apad2\Desktop\Fat_water_separation\DICOM_Files";
pth_main = pth + "_Project DMDiv\Dog Data";
pth_code = "/scratch/user/apad/GRMD";
```
to
```matlab
% pth_main = pth + "_Project DMDiv\Dog Data";                       % legacy, unused on Linux
pth_data = "/scratch/user/apad/Fat_water_separation/DICOM_Files";   % location of DD.mat (MUST match dogexplorer's mainpth)
pth_code = "/scratch/user/apad/GRMD";                               % code root (kept for reference)
```

Change the `DogInitialize` call (line 53) from
```matlab
[snames, DD, HDR, pth_1H, plrange, vout, filelocs] = DogInitialize(pth_code, dog, date, flags, cscFLAG);
```
to
```matlab
[snames, DD, HDR, pth_1H, plrange, vout, filelocs] = DogInitialize(pth_data, dog, date, flags, cscFLAG);
```

Note: `pth_data` here MUST match `mainpth` in `dogexplorer.m` (line 1559). If the data root moves, update both.

Add four `addpath` lines (after the path variables) so the code and its dependencies are on the MATLAB path:
```matlab
addpath(genpath(pth_code))                                 % GRMD
addpath(genpath("/scratch/user/apad/CREAM_PDFF/hernando")) % hernando
addpath(genpath("/scratch/user/apad/functions"))           % external utilities (loadima, get_echoMIP, GenPlotRange, plotmygraph, ...)
addpath(genpath("/scratch/user/apad/bipolar_fat_water_separation")) % Function_Bipolar_GC (BipolarIGC slice separation)
```

## loadima.m
DICOM file paths were built with a hardcoded backslash separator, which fails on Linux (only the raw-DICOM preprocessing path hits this; it is skipped when preprocessed `.mat` files already exist). Change lines 206, 207, 209, 210, 226, 227, 229, 230 from
```matlab
data(:,:,:,j)   = dicomread([dcm '\' filelist(j).name]);
hdr(j)          = dicominfo([dcm '\' filelist(j).name]);
data_p(:,:,:,j) = dicomread([dcm_p '\' filelist_p(j).name]);
hdr_tmp         = dicominfo([dcm_p '\' filelist_p(j).name]);
data(:,:,i,j)   = dicomread([dcm '\' filelist(idx).name]);
hdr(i,j)        = dicominfo([dcm '\' filelist(idx).name]);
data_p(:,:,i,j) = dicomread([dcm_p '\' filelist_p(idx).name]);
hdr_tmp         = dicominfo([dcm_p '\' filelist_p(idx).name]);
```
to use `fullfile`, e.g.
```matlab
data(:,:,:,j)   = dicomread(fullfile(dcm, filelist(j).name));
hdr(j)          = dicominfo(fullfile(dcm, filelist(j).name));
data_p(:,:,:,j) = dicomread(fullfile(dcm_p, filelist_p(j).name));
hdr_tmp         = dicominfo(fullfile(dcm_p, filelist_p(j).name));
data(:,:,i,j)   = dicomread(fullfile(dcm, filelist(idx).name));
hdr(i,j)        = dicominfo(fullfile(dcm, filelist(idx).name));
data_p(:,:,i,j) = dicomread(fullfile(dcm_p, filelist_p(idx).name));
hdr_tmp         = dicominfo(fullfile(dcm_p, filelist_p(idx).name));
```
Leave line 64 (`strsplit(hdr.ImageType,'\')`) unchanged — that backslash is the DICOM `ImageType` field delimiter, not a filesystem path.