function outParams = calculateMinimaDirect(imDataParams, VARPROparams, verbose)
    % Check for verbose flag
    if nargin < 3
        verbose = false;
    end

    %% Define and calculate key numbers of the input image
    try
        nMaxMinimizers = 2*size(VARPROparams.species(2).frequency, 2)*round(VARPROparams.nSamplingPeriods);
    catch
        nMaxMinimizers = round(VARPROparams.nSamplingPeriods);
    end
    
    if nMaxMinimizers < 10
        nMaxMinimizers = 10;
        if verbose
            disp(['nMaxMinimizers = ', nMaxMinimizers])
        end
    end
            
    [nVoxel_Y, nVoxel_X, nVoxel_Z, ~, ~] = size(imDataParams.images);

    costLocalMinimaRescale = zeros(nVoxel_Y, nVoxel_X, nVoxel_Z, nMaxMinimizers);
    indexLocalMinimaRescale = zeros(nVoxel_Y, nVoxel_X, nVoxel_Z, nMaxMinimizers);
    nMinimaPerVoxel = zeros(nVoxel_Y, nVoxel_X, nVoxel_Z);
    INFTY = 1e8;
    gyro = VARPROparams.gyro;
    if verbose
        fprintf('Calculating Cost function (VARPRO Hernando)\n');
    end
    options = struct();
    options.nMaxMinimizers = nMaxMinimizers;
    options.rescale = [1/INFTY, 1];
    options.minDistance = 1.2*min(abs(diff(VARPROparams.species(2).frequency)))*gyro*imDataParams.FieldStrength;

    if ~isscalar(options.minDistance)
        options.minDistance = VARPROparams.sampling_stepsize;
    end

    pool = parpool;
    opts = parforOptions(pool, 'MaxNumWorkers',4);
    parfor (Z = 1:nVoxel_Z, opts)
        if verbose
            fprintf('Slice %d of %d: ', Z, nVoxel_Z);
        end
        tmp_imDataParams = imDataParams;
        tmp_imDataParams.images = imDataParams.images(:, :, Z, :, :);
        t1 = tic;
        [residual, ~] = Boehm_computeResidual(tmp_imDataParams, VARPROparams, verbose);
        t1 = toc(t1);
        if verbose
            fprintf('residual calculation done! (%.2fs)', t1);
        end
        %% Transfer to array
        tmp_masksignal = get_tissueMask(tmp_imDataParams.images, VARPROparams.airSignalThreshold_percent);
        t2 = tic;
        [tmp_nMinimaPerVoxel, tmp_costLocalMinimaRescale, tmp_indexLocalMinimaRescale] = findLocalMinima_and_rescale(residual, tmp_masksignal, options);
        t2 = toc(t2);
        if verbose
            fprintf('minima extraction done! (%.2fs)\n', t2);
        end
        tmp_indexLocalMinimaRescale = VARPROparams.gridspacing .* tmp_indexLocalMinimaRescale;

        %% Transfer to the 4D - Arrays
        costLocalMinimaRescale(:, :, Z, :) = tmp_costLocalMinimaRescale;
        indexLocalMinimaRescale(:, :, Z, :) = tmp_indexLocalMinimaRescale;
        nMinimaPerVoxel(:, :, Z) = tmp_nMinimaPerVoxel;
    end
    outParams.costLocalMinimaRescale = costLocalMinimaRescale;
    outParams.indexLocalMinimaRescale = indexLocalMinimaRescale;
    outParams.masksignal = nMinimaPerVoxel > 0;
    outParams.nMinimaPerVoxel = nMinimaPerVoxel;
end