% Grab chemical shift correction flags
function algoParams = grabPars(flags)
    % Isolate flags
    verboseFLAG = flags.verbose;
    flags = flags.cscorrection;

    % Ensure input is valid
    if ~any(strcmp(flags.method, {'IGC', 'MixedIGC', 'BipolarIGC', 'vlGC', 'IDEAL-CE', 'SPURS', 'Hierarchical IDEAL', 'Golden Section'}))
        error('The entered technique for chemical shift correction has not been added to this pipeline.')
    end
    algoParams = struct('useCUDA', true);
    
    % Algorithm-specific parameters
    switch flags.method
        case {'IGC', 'MixedIGC'}
            algoParams.species = struct('name', {'water', 'fat'}, ...
                                   'frequency', {0, [-3.8,  -3.4,  -2.6,  -1.94, -0.39,  0.6]}, ...
                                     'relAmps', {1, [ 0.087, 0.693, 0.128, 0.004, 0.039, 0.048]});
            algoParams.size_clique           = 1;               % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
            algoParams.range_r2star          = [0 300];         % Range of R2* values
            algoParams.NUM_R2STARS           = 101;             % Number of R2* values for quantization (default = 11)
            algoParams.range_fm              = [-500 500];      % Range of field map values
            algoParams.NUM_FMS               = 1001;            % Number of field map values to discretize (default = 300)
            algoParams.NUM_ITERS             = 40;              % Number of graph cut iterations
            algoParams.SUBSAMPLE             = flags.subsample; % Spatial subsampling for field map estimation (for speed)
            algoParams.DO_OT                 = 1;               % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
            algoParams.LMAP_POWER            = 2;               % Spatially-varying regularization (2 gives ~ uniformn resolution)
            algoParams.lambda                = 0.10;            % Regularization parameter
            algoParams.LMAP_EXTRA            = 0.05;            % More smoothing for low-signal regions
            algoParams.TRY_PERIODIC_RESIDUAL = 0;               % Take advantage of periodic residual if uniform TEs (will change range_fm)  
            
            if strcmp(flags.method, 'IGC-mixed')
                algoParams.mixed.NUM_MAGN        = 1;               % Number of potentially phase-corrupted echoes (default = 1)
                algoParams.mixed.THRESHOLD       = 0;               % Signal threshold for processing voxels (default = 0)
                algoParams.mixed.range_r2star    = [0 300];         % Range of R2* values
            end

        case 'BipolarIGC'
            algoParams.species = struct('name', {'water', 'fat'}, ...
                                   'frequency', {0, [-3.8,  -3.4,  -2.6,  -1.94, -0.39,  0.6]}, ...
                                     'relAmps', {1, [ 0.087, 0.693, 0.128, 0.004, 0.039, 0.048]});
            algoParams.size_clique           = 1;               % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
            algoParams.range_r2star          = [0 300];         % Range of R2* values
            algoParams.NUM_R2STARS           = 101;             % Number of R2* values for quantization (default = 11)
            algoParams.range_fm              = [-500 500];      % Range of field map values
            algoParams.NUM_FMS               = 1001;            % Number of field map values to discretize (default = 300)
            algoParams.NUM_ITERS             = 120;             % Number of graph cut iterations
            algoParams.SUBSAMPLE             = flags.subsample; % Spatial subsampling for field map estimation (for speed)
            algoParams.DO_OT                 = 0;               % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
            algoParams.LMAP_POWER            = 2;               % Spatially-varying regularization (2 gives ~ uniformn resolution)
            algoParams.lambda                = 0.10;            % Regularization parameter
            algoParams.LMAP_EXTRA            = 0.2;             % Settle on 0.2 for smoothing
            algoParams.TRY_PERIODIC_RESIDUAL = 0;               % Take advantage of periodic residual if uniform TEs (will change range_fm)  
            algoParams.tik_reg               = 0;               % Tikhonov regularization binary flag
            algoParams.plot_debug            = false;
            algoParams.parallel              = true;
            
            if exist("/scratch/user/apad/residuals/checkpoint.mat", "file")
                load("/scratch/user/apad/residuals/checkpoint.mat");
                algoParams.residual = residual;
            end

        case 'vlGC'
            % Algorithm-specific parameters
            algoParams.species = struct('name', {'water', 'fat'}, ...
                                   'frequency', {0, [-3.8,    -3.4,    -3.1,    -2.68,   -2.46,   -1.95,     -0.5,    0.49, 0.59]}, ...
                                     'relAmps', {1,  [0.0899,  0.5834,  0.0599,  0.0849,  0.0599,  0.0150,    0.04,   0.01, 0.0569]});
            algoParams.range_r2star               = [0 300];         % Range of R2* values
            algoParams.NUM_R2STARS                = 101;             % Number of R2* values for quantization (default = 11)
            algoParams.range_fm                   = [-500 500];      % Range of field map values
            algoParams.sampling_stepsize          = flags.subsample; % Spatial subsampling (for speed)
            algoParams.airSignalThreshold_percent = 0;               % threshold for masking (default = 5)

        case 'IDEAL-CE'
            % Algorithm-specific parameters along with data
            algoParams = struct('species', struct('name', {'water', 'fat'}, ...
                                             'frequency', {0, [-3.8,  -3.4,  -2.6,  -1.94, -0.39,  0.6]}, ...
                                               'relAmps', {1, [ 0.087, 0.693, 0.128, 0.004, 0.039, 0.048]}), ...
                          'ShiftedkSpace', false,...
                                    'muB', 0.00, ...
                                    'muR', 0.01, ...
                          'SmoothedPhase', false, ...
                          'SmoothedField', false, ...
                                 'Filter', ones(3), ...
                               'NoiseSTD', [], ...
                          'MaxIterations', [10 10 5]);

        case 'Hierarchical IDEAL'
            algoParams.species = struct('name', {'water', 'fat'}, ...
                                   'frequency', {0, [-3.8,  -3.4,  -2.6,  -1.94, -0.39,  0.6]}, ...
                                     'relAmps', {1, [ 0.087, 0.693, 0.128, 0.004, 0.039, 0.048]});
            algoParams.Verbose                          = verboseFLAG;      % Show info (default = 1)
            algoParams.AlwaysShowGUI                    = false;            % Always show GUI to verify input (default = 1)
            algoParams.Visualize                        = false;            % Show graphics (default = 1)
            algoParams.Visualize_FatWaterMapMultiplier  = 1.0;              % Fat fraction multiplier for visualization (default = 1.5)
            algoParams.MinFractSizeToDivide             = 0.05;             % Min size of the image to stop dividing into smaller regions (default = 0.05)
            algoParams.MaxNumDiv                        = 6;                % Max number of hierarchical sub-divisions (default = 6)
            algoParams.AssumesSinglePeakAsWater         = true;             % If there is only one peak, assume it is water (default = 1)
            algoParams.SnrToAssumeSinglePeak            = 2.5;              % Min SNR to assess whether or not it is noise (default = 2.5)
            algoParams.CorrectAmpForT2star              = 1;                % Produce water & fat maps that have been corrected for T2* decay (default = 1)
            algoParams.MaxR2star                        = 250;              % Max acceptable R2* in case of erroneously large values from bad fit (default = 250)
    
        case 'Golden Section'
            algoParams.species = struct('name', {'water', 'fat'}, ...
                                   'frequency', {0, [-3.8,  -3.4,  -2.6,  -1.94, -0.39,  0.6]}, ...
                                     'relAmps', {1, [ 0.087, 0.693, 0.128, 0.004, 0.039, 0.048]});
    end
end