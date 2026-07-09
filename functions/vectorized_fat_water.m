function D = vectorized_fat_water(D, opts)
%VECTORIZED_FAT_WATER Water/fat separation via the vectorized Bipolar_GC toolbox.
%
%   D = vectorized_fat_water(D)
%   D = vectorized_fat_water(D, opts)
%
%   Drop-in alternative to HERNANDO_FAT_WATER that uses Jorge Campos' vectorized
%   graph-cut toolbox (toolboxes/vectorized, "Bipolar_GC"). It runs
%   Function_Bipolar_GC on the complex multi-echo volume in D.Data.Image and
%   fills in the same fields as hernando_fat_water:
%       D.Data.Water      - complex water image      [nx ny nz]
%       D.Data.Fat        - complex fat image        [nx ny nz]
%       D.Data.TotalField - B0 field map (Hz)        [nx ny nz]
%       D.Data.R2StarMap  - R2* map (s^-1)           [nx ny nz]
%   plus the bipolar-specific correction maps:
%       D.Data.PhiMap     - phase modulation (rad)   [nx ny nz]
%       D.Data.EpsMap     - amplitude modulation     [nx ny nz]
%
%   IMPORTANT: this backend performs its OWN bipolar readout correction
%   internally by splitting the echo train into odd/even sub-trains. It must be
%   run on RAW (uncorrected) bipolar data -- do NOT run bipolar_correct first,
%   or the readout phase will be corrected twice. It therefore requires an EVEN
%   number of echoes.
%
%   opts fields (all optional):
%       .subsample  : spatial subsampling of the field map for speed (default 1)
%       .verbose    : print progress (default true)
%       .species    : override the water/fat spectral model (struct array)
%       .range_fm   : field-map search range in Hz (default [-500 500])
%       .plot_debug : show the toolbox's intermediate figures (default false)
%       .weight     : odd/even correction weight for the phi/eps fit (default 0.5)
%       .tik_reg    : Tikhonov regularization in the phi/eps fit (default 0)
%       .fm_init    : use an initial field-map guess (default 0)
%       .trim_last_echo : if the echo count is odd, drop the last echo to make
%                         it even instead of erroring (default true). The
%                         trimmed D.Data.Image / D.TE are returned, so the saved
%                         result reflects the echoes actually used.
%       .gc_scale_target : target magnitude for the max-flow min-cut so its
%                         int32 edge capacities neither overflow nor lose
%                         precision (default 1e8). Raise/lower only if you see
%                         the "minimum cut ~= max-flow" overflow warning at the
%                         default. Does not change the result, only the scaling.
%
%   Requires the Optimization Toolbox (lsqlin) and, on the path, both
%   toolboxes/vectorized and toolboxes/hernando (Function_Bipolar_GC calls the
%   plain Hernando graph-cut on the corrected "synthetic unipolar" signal).
%
%   Reference: Campos J. A flexible approach for fat-water separation with
%   bipolar readouts and correction of gradient-induced phase and amplitude
%   effects. Based on Hernando D et al. Magn Reson Med. 2010;63(1):79-90.

    if nargin < 2, opts = struct(); end
    if ~isfield(opts, 'subsample'),  opts.subsample  = 1;          end
    if ~isfield(opts, 'verbose'),    opts.verbose    = true;       end
    if ~isfield(opts, 'range_fm'),   opts.range_fm   = [-500 500]; end
    if ~isfield(opts, 'plot_debug'), opts.plot_debug = false;      end
    if ~isfield(opts, 'weight'),     opts.weight     = 0.5;        end
    if ~isfield(opts, 'tik_reg'),    opts.tik_reg    = 0;          end
    if ~isfield(opts, 'fm_init'),    opts.fm_init    = 0;          end
    if ~isfield(opts, 'trim_last_echo'), opts.trim_last_echo = true; end
    % Target magnitude for the max-flow cut so its int32 edge capacities neither
    % overflow (INT32_MAX = 2.147e9) nor lose precision under floor(). See the
    % normalization block below.
    if ~isfield(opts, 'gc_scale_target'), opts.gc_scale_target = 1e8; end

    [nx, ny, SL, nTE] = size(D.Data.Image);

    % ---- Bipolar_GC needs paired odd/even echoes ----
    % The algorithm splits the echo train into odd/even sub-trains, so nTE must
    % be even. If it is odd, optionally drop the final echo (trimming D so the
    % saved result reflects the echoes actually used); otherwise error.
    if mod(nTE, 2) ~= 0
        if opts.trim_last_echo
            warning('vectorized_fat_water:trimLastEcho', ...
                    ['Odd echo count (%d); dropping the last echo (TE = %g) to ' ...
                     'make it even for Bipolar_GC.'], nTE, D.TE(end));
            D.Data.Image = D.Data.Image(:,:,:,1:end-1);
            D.TE         = D.TE(1:end-1);
            nTE          = nTE - 1;
        else
            error('vectorized_fat_water:oddEchoes', ...
                  ['The Bipolar_GC backend splits the echo train into odd/even ' ...
                   'sub-trains and requires an even number of echoes (got %d). ' ...
                   'Set opts.trim_last_echo = true to drop the last echo.'], nTE);
        end
    end

    % ---- Multi-peak fat model + graph-cut algorithm parameters ----
    % Same 6-peak fat spectrum and core graph-cut settings as hernando_fat_water,
    % so results are directly comparable. Function_Bipolar_GC forms the fat/water
    % chemical shifts as gyro*(species(2).frequency - species(1).frequency(1))*B0,
    % so the water-at-0 / fat-negative-ppm convention matches the Hernando path.
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
    algoParams.DO_OT                 = 0;               % OT is forced off for the final unipolar fit anyway
    algoParams.LMAP_POWER            = 2;               % spatially-varying regularization
    algoParams.lambda                = 0.05;            % regularization parameter
    algoParams.LMAP_EXTRA            = 0.05;            % extra smoothing for low-signal regions
    algoParams.TRY_PERIODIC_RESIDUAL = 0;               % periodic residual (uniform TEs)
    algoParams.THRESHOLD             = 0.01;            % low-signal threshold

    % ---- Bipolar-specific parameters ----
    algoParams.slice_image      = max(1, round(SL/2)); % slice used for the debug figures
    algoParams.plot_debug       = opts.plot_debug;     % show intermediate figures
    algoParams.crameri_colormap = 0;                   % use MATLAB colormaps (no crameri dependency)
    algoParams.tik_reg          = opts.tik_reg;        % Tikhonov reg in the phi/eps fit
    algoParams.weight           = opts.weight;         % odd/even correction weight
    algoParams.fm_init          = opts.fm_init;        % initial field-map guess flag

    % ---- Format data for the graph-cut routine ----
    if isfield(D.Data, 'Mask')
        mask = D.Data.Mask;
    else
        mask = true(nx, ny, SL);   % separate everywhere if no mask available
    end

    imDataParams.FieldStrength         = D.B0;   % Tesla (D.B0 = ImagingFrequency / gyro)
    imDataParams.TE                    = D.TE;   % echo times (s)
    imDataParams.PrecessionIsClockwise = 1;
    imDataParams.voxelSize             = D.VoxelSize;
    imDataParams.mask_fwseparation     = mask;

    % ---- Normalize signal magnitude for the integer max-flow graph cut ----
    % matlab_bgl's max_flow floors edge capacities to int32. The capacities are
    % built from the fit residual, which scales as |signal|^2, so raw scanner
    % magnitudes make the min-cut overflow int32 -- the symptom is the warning
    % "the rounded (unrounded) value of the minimum cut is -2147483648 (...),
    % but the value of the max-flow is ...". The min-cut *solution* is invariant
    % to a global scale of the signal, so we rescale the complex data into a
    % safe band (cut ~ gc_scale_target, comfortably under INT32_MAX and well
    % above 1 for floor() precision) and undo the scale on the linear water/fat
    % outputs afterwards. Field map / R2* / phi / eps are scale-invariant.
    energy = sum(abs(D.Data.Image(:)).^2);
    if energy > 0
        sigScale = sqrt(opts.gc_scale_target / energy);
    else
        sigScale = 1;
    end

    imDataParams.images = zeros(nx, ny, SL, 1, nTE);
    imDataParams.images(:,:,:,1,:) = sigScale * D.Data.Image;

    % ---- Run the vectorized Bipolar_GC separation (all slices) ----
    vec_slices = 1:SL;

    fprintf(['\nPerforming vectorized Bipolar_GC fat/water separation ' ...
             '(%d slices, signal scaled by %.3g for the graph cut)...\n'], SL, sigScale);
    out = Function_Bipolar_GC(imDataParams, algoParams, vec_slices, opts.verbose);

    % ---- Reassemble volumes (undo the graph-cut scaling on water/fat) ----
    D.Data.Water      = out.species(1).amps / sigScale;
    D.Data.Fat        = out.species(2).amps / sigScale;
    D.Data.TotalField = out.fieldmap;
    D.Data.R2StarMap  = out.r2starmap;
    D.Data.PhiMap     = out.phi_map;   % phase modulation from bipolar readout (rad)
    D.Data.EpsMap     = out.eps_map;   % amplitude modulation from bipolar readout
    D.FieldMap        = D.Data.TotalField;

    % ---- Bookkeeping ----
    D.Size = size(D.Data.TotalField);
    % Bipolar_GC folds the readout correction into the separation, so both flags
    % are satisfied by this single step.
    D.Flags.CorrectedBipolarPhase   = true;
    D.Flags.CorrectedChemicalShift  = true;
    D.BipolarCorrection = struct('Method', 'Bipolar_GC (vectorized, internal odd/even)');
    D.CSCorrection = struct('Method', 'Bipolar_GC (Campos vectorized graph-cut)', ...
                   'SubsampleFactor', algoParams.SUBSAMPLE, ...
                           'Species', algoParams.species, ...
                     'FieldMapRange', algoParams.range_fm, ...
                   'FieldMapNumbers', algoParams.NUM_FMS, ...
                'GraphCutIterations', algoParams.NUM_ITERS, ...
                            'Lambda', algoParams.lambda, ...
                            'Weight', algoParams.weight);
end
