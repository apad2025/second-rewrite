function D = hernando_fat_water(D, opts)
%HERNANDO_FAT_WATER Water/fat separation via Hernando's graph-cut algorithm.
%
%   D = hernando_fat_water(D)
%   D = hernando_fat_water(D, opts)
%
%   Runs the regularized-fieldmap graph-cut fat/water separation
%   (fw_i2cm1i_3pluspoint_hernando_graphcut) slice-by-slice on the complex
%   multi-echo volume in D.Data.Image, filling in:
%       D.Data.Water      - complex water image      [nx ny nz]
%       D.Data.Fat        - complex fat image        [nx ny nz]
%       D.Data.TotalField - B0 field map (Hz)        [nx ny nz]
%       D.Data.R2StarMap  - R2* map (s^-1)           [nx ny nz]
%
%   opts fields (all optional):
%       .subsample  : spatial subsampling of the field map for speed (default 1)
%       .numWorkers : parallel workers for the per-slice parfor (default 2)
%       .verbose    : print per-slice progress (default true)
%       .species    : override the water/fat spectral model (struct array)
%       .range_fm   : field-map search range in Hz (default [-600 600])
%
%   This is the "IGC" (iterative graph-cut) path condensed from
%   chemicalshiftcorrection / GetCSCParams / IGC in the original DogAnalysis.m.
%
%   Reference: Hernando D, Kellman P, Haldar JP, Liang ZP. Robust water/fat
%   separation in the presence of large field inhomogeneities using a graph
%   cut algorithm. Magn Reson Med. 2010;63(1):79-90.

    if nargin < 2, opts = struct(); end
    if ~isfield(opts, 'subsample'),  opts.subsample  = 1;    end
    if ~isfield(opts, 'numWorkers'), opts.numWorkers = 2;    end
    if ~isfield(opts, 'verbose'),    opts.verbose    = true; end
    if ~isfield(opts, 'range_fm'),   opts.range_fm   = [-600 600]; end

    % ---- Multi-peak fat model + graph-cut algorithm parameters ----
    % (same values as the original DogAnalysis 'IGC' pipeline)
    if isfield(opts, 'species')
        algoParams.species = opts.species;
    else
        algoParams.species = struct('name', {'water', 'fat'}, ...
                               'frequency', {0, [-3.8,  -3.4,  -2.6,  -1.94, -0.39,  0.6]}, ...
                                 'relAmps', {1, [ 0.087, 0.693, 0.128, 0.004, 0.039, 0.048]});
    end
    algoParams.size_clique           = 1;               % MRF neighborhood size (8-neighborhood)
    algoParams.range_r2star          = [0 100];         % R2* range (s^-1)
    algoParams.NUM_R2STARS           = 101;             % R2* quantization levels
    algoParams.range_fm              = opts.range_fm;   % field map range (Hz)
    algoParams.NUM_FMS               = 1001;            % field map discretization
    algoParams.NUM_ITERS             = 40;              % graph cut iterations
    algoParams.SUBSAMPLE             = opts.subsample;  % spatial subsampling for speed
    algoParams.DO_OT                 = 1;               % optimization-transfer refinement
    algoParams.LMAP_POWER            = 2;               % spatially-varying regularization
    algoParams.lambda                = 0.05;            % regularization parameter
    algoParams.LMAP_EXTRA            = 0.05;            % extra smoothing for low-signal regions
    algoParams.TRY_PERIODIC_RESIDUAL = 0;               % periodic residual (uniform TEs)

    % ---- Format data for the graph-cut routine ----
    dataParams = struct('FieldStrength', D.B0, ...
                                   'TE', D.TE, ...
                'PrecessionIsClockwise', 1);

    [nx, ny, SL, nTE] = size(D.Data.Image);
    images = zeros(nx, ny, SL, 1, nTE);
    images(:,:,:,1,:) = D.Data.Image;

    % Per-slice cells (parfor-friendly)
    Wc  = cell(SL,1);  Fc  = cell(SL,1);
    FMc = cell(SL,1);  R2c = cell(SL,1);
    dp  = cell(SL,1);
    for sl = 1:SL
        dp{sl} = dataParams;
        dp{sl}.images = images(:,:,sl,1,:);
    end

    verboseFLAG = opts.verbose;

    % ---- Ensure a parallel pool (guarded) ----
    if isempty(gcp('nocreate'))
        try
            parpool(opts.numWorkers);
        catch
            % fall back to serial if the Parallel Computing Toolbox is unavailable
        end
    end

    fprintf('\nPerforming Hernando graph-cut fat/water separation (%d slices)...\n', SL);
    parfor (sl = 1:SL, opts.numWorkers)
        if verboseFLAG, fprintf('  slice %d...\n', sl); end
        % Pass verbose as Hernando's DEBUG flag so its per-stage progress
        % ("Estimating field map...", "R2*...", "optimization transfer...")
        % is printed -- lets you confirm a slow slice is working, not hung.
        out = fw_i2cm1i_3pluspoint_hernando_graphcut(dp{sl}, algoParams, verboseFLAG);
        Wc{sl}  = out.species(1).amps;
        Fc{sl}  = out.species(2).amps;
        FMc{sl} = out.fieldmap;
        R2c{sl} = out.r2starmap;
    end
    delete(gcp('nocreate'))

    % ---- Reassemble volumes ----
    D.Data.Water      = zeros(nx, ny, SL);
    D.Data.Fat        = zeros(nx, ny, SL);
    D.Data.TotalField = zeros(nx, ny, SL);
    D.Data.R2StarMap  = zeros(nx, ny, SL);
    for sl = 1:SL
        D.Data.Water(:,:,sl)      = Wc{sl};
        D.Data.Fat(:,:,sl)        = Fc{sl};
        D.Data.TotalField(:,:,sl) = FMc{sl};
        D.Data.R2StarMap(:,:,sl)  = R2c{sl};
    end
    D.FieldMap = D.Data.TotalField;

    % ---- Bookkeeping ----
    D.Size = size(D.Data.TotalField);
    D.Flags.CorrectedChemicalShift = true;
    D.CSCorrection = struct('Method', 'IGC (Hernando graph-cut)', ...
                   'SubsampleFactor', algoParams.SUBSAMPLE, ...
                           'Species', algoParams.species, ...
                     'FieldMapRange', algoParams.range_fm, ...
                   'FieldMapNumbers', algoParams.NUM_FMS, ...
                'GraphCutIterations', algoParams.NUM_ITERS, ...
                            'Lambda', algoParams.lambda);
end
