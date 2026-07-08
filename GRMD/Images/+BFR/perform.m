% Remove background field
function data = perform(data, flags, dims, verboseFLAG) 
% Inputs:
%     data: data structure
%    flags: processing structure
% 
% Outputs:
%     data: data structure

    % Ensure input is valid
    if ~any(strcmp(flags.method, {'PDF - MEDI', 'PDF - QSMmaster', 'SHARP', 'V-SHARP - SEPIA', 'V-SHARP - STISuite', 'LBV'}))
        if strcmp(flags.method, 'PDF')
            PDFmethod = input('Two PDF algorithms are available: (1) MEDI and (2) QSMmaster. Which would you prefer?','s');
            if strcmpi(PDFmethod, {'MEDI', '1', '(1)'})
                flags.method = 'PDF - MEDI';
            elseif strcmpi(PDFmethod, {'QSMmaster', '2', '(2)'})
                flags.method = 'PDF - QSMmaster';
            else
                error('The entered input could not be interpreted.')
            end
        elseif strcmp(flags.method, 'V-SHARP')
            VSHARPmethod = input('Two V-SHARP algorithms are available: (1) SEPIA and (2) STISuite. Which would you prefer?','s');
            if strcmpi(VSHARPmethod, {'SEPIA', '1', '(1)'})
                flags.method = 'V-SHARP - SEPIA';
            elseif strcmpi(VSHARPmethod, {'STISuite', '2', '(2)'})
                flags.method = 'V-SHARP - STISuite';
            else
                error('The entered input could not be interpreted.')
            end
        else
            error('The entered technique for background field removal has not been added to this pipeline.')
        end
    end

    % Ensure B0 direction is positive
    [data, flippedFLAG] = Operations.Flip(data);

    % Check slice direction
    if ~isequal(data.B0Direction,[0 0 1])
        disp('This is angled slicing');
        disp(data.B0Direction);
    end

    fprintf('\nPerforming %iD background field removal using %s...', dims, flags.method)

    % Set algorithm specific parameters
    switch flags.method
        case 'PDF - MEDI'
            params = struct('tol', 0.1, ...                                Default 0.1         , tolerance level
                           'iter', 100, ...                                Default 30          , number of conjugate gradient iterations
                          'space', 'imagespace', ...                       Default 'imagespace', domain for kernel creation
                          'psize', 40);                                  % Default 40          , padding size
        case 'PDF - QSMmaster'
            params = struct('iter', 200);                                % Default 200         , number of iterations

        case 'SHARP'
            params = struct('kerrad', 0.45, ...                            Default 0.5         , (5*pixel) convolution kernel size (mm)
                              'tsvd', 0.025);                            % Default 0.05        , truncated singular value decomposition

        case 'V-SHARP - SEPIA'
            params = struct('kerrads', 8:-1:1);                          % Default [5,4,3,2,1] , convolution kernel sizes (#pixels)

        case 'V-SHARP - STISuite'
            params = struct('psize', data.Size, ...                        Default [12,12,12]  , padding size
                          'smvsize', 12);                                % Default 12          , kernel size (#pixels)

        case 'LBV'
            params = struct('tol', 0.001, ...                              Default 0.01        , stopping criteria on coursest grid
                           'peel', 0, ...                                  Default 0           , number of boundary layers to be peeled (quick but dirty)
                          'depth', -1, ...                                 Default -1          , number of length scales
                             'n1', 30, ...                                 Default 30          , iterations on each depth before recursive call
                             'n2', 100, ...                                Default 100         , iterations on each depth after recursive call
                             'n3', 100);                                 % Default 100         , iterations on finest scale after FMG
    end

    % Remove background field
    data.Data.LocalField = zeros(data.Size, 'like', data.Data.TotalField);
    switch flags.method
        case 'PDF - MEDI'
            %%% Projection onto Dipole Fields (PDF) - MEDI Toolbox
            %   [rdf, shim] = PDF(iFreq, N_std, mask, matrix_size, voxel_size, B0_dir, tol, n_CG, space, n_pad)
            %
            % Inputs
            %   iFreq       : unwrapped field map
            %   N_std       : noise standard deviation on the field map. (1 over SNR for single echo)
            %   mask        : binary 3D matrix denoting the Region Of Interest
            %   matrix_size : size of the 3D matrix
            %   voxel_size  : size of the voxel in mm
            %   B0_dir      : direction of the B0 field
            %   tol         : tolerance level (optional)
            %   n_CG        : conjugate gradient iterations (optional)
            %   space       : time or fourier space selection (optional)
            %   n_pad       : pad size (optional)
            %
            % Outputs
            %   rdf         : relative difference field, or local field
            %   shim        : cropped background dipole distribution

            data.Data.BackgroundField = zeros(data.Size, 'like', data.Data.TotalField);

            % Check 2D or 3D
            if dims == 2
                for sl = 1:data.Size(3)
                    if verboseFLAG; fprintf('\nSlice %i...', sl); end
                    [data.Data.LocalField(:,:,sl), data.Data.BackgroundField(:,:,sl)] = PDF(data.Data.TotalField(:,:,sl), data.Data.NoiseSTD(:,:,sl), data.Data.Mask(:,:,sl), [data.Size(1), data.Size(2)], data.VoxelSize, data.B0Direction, params.tol, params.iter, params.space, params.psize);
                end
            else
                [data.Data.LocalField, data.Data.BackgroundField] = PDF(data.Data.TotalField, data.Data.NoiseSTD, data.Data.Mask, data.Size, data.VoxelSize, data.B0Direction, params.tol, params.iter, params.space, params.psize);
            end

        case 'PDF - QSMmaster'
            %%% PDF - QSM-master
            %   [lfs, bkg_sus, bkg_field] = projectionontodipolefields(tfs, mask, vox, weight, v_norm, num_iter)
            %
            % Inputs
            %   tfs         : input total field shift
            %   mask        : binary mask defining the brain ROI
            %   vox         : voxel size (mm)
            %   weight      : weighting
            %   v_norm      : imaging axis
            %   num_iter    : maximum number of iterations, e.g. 200
            %
            % Outputs
            %   lfs         : local field shift after background removal
            %   bkg_sus     : background susceptibility
            %   bkg_field   : background field

            data.Data.BackgroundField = zeros(data.Size, 'like', data.Data.TotalField);

            if dims == 2
                for sl = 1:data.Size(3)
                    if verboseFLAG; fprintf('\nSlice %i...', sl); end
                    [data.Data.LocalField(:,:,sl), ~, data.Data.BackgroundField(:,:,sl)] = projectionontodipolefields(data.Data.TotalField(:,:,sl), data.Data.Mask(:,:,sl), data.VoxelSize, data.Data.WeightedMagnitude(:,:,sl), data.B0Direction, params.iter);
                end
            else
                [data.Data.LocalField, ~, data.Data.BackgroundField] = projectionontodipolefields(data.Data.TotalField, data.Data.Mask, data.VoxelSize, data.Data.WeightedMagnitude, data.B0Direction, params.iter);
            end

        case 'SHARP'
            %%% SHARP - QSM-master
            %   [lfs, m_ero] = SHARP(tfs, mask, vox, ker_rad, tsvd)
            %
            % Inputs
            %   tfs         : input total field shift
            %   mask        : binary mask defining the ROI
            %   vox         : voxel size (mm)
            %   ker_rad     : radius of convolution kernel (mm), e.g. 5
            %   tsvd        : truncated singular value decomposition, e.g. 0.05
            %
            % Outputs
            %   lfs         : local field shift after background removal
            %   m_ero       : eroded mask after convolution
            
            if dims == 2
                for sl = 1:data.Size(3)
                    if verboseFLAG; fprintf('\nSlice %i...', sl); end
                    [data.Data.LocalField(:,:,sl), data.Data.Mask(:,:,sl)] = sharp(data.Data.TotalField(:,:,sl), data.Data.Mask(:,:,sl), data.VoxelSize, params.kerrad, params.tsvd);
                end
            else
                [data.Data.LocalField, data.Data.Mask] = sharp(data.Data.TotalField, data.Data.Mask, data.VoxelSize, params.kerrad, params.tsvd);
            end

        case 'V-SHARP - SEPIA'
            %%% V-SHARP - SEPIA
            %   [lfs, mask_ero] = V_SHARP(tfs, mask, 'voxelsize', vsize, 'smvsize', smvsize)
            %
            % Inputs
            %   tfs         : input total field shift
            %   mask        : binary mask defining the brain ROI
            %   'radius'    : radius of convolution kernels (#pixels), default [5,4,3,2,1]
            %
            % Outputs
            %   lfs         : local field shift after background removal
            %   mask_ero    : eroded mask after convolution

            if dims == 2
                reverseStr = '';
                for sl = 1:data.Size(3)
                    if verboseFLAG; reverseStr = UpdatePercent(100*sl/data.Size(3), reverseStr); end
                    [data.Data.LocalField(:,:,sl), data.Data.Mask(:,:,sl)] = BKGRemovalVSHARP_2D(data.Data.TotalField(:,:,sl), data.Data.Mask(:,:,sl), [data.Size(1), data.Size(2)], 'radius', params.kerrads);
                end
                fprintf('\tDone!\n');
            else
                [data.Data.LocalField, data.Data.Mask] = BKGRemovalVSHARP(data.Data.TotalField, data.Data.Mask, data.Size, 'radius', params.kerrads);
            end


        case 'V-SHARP - STISuite'
            %%% V-SHARP - STI Suite
            %   TissuePhase = V_SHARP_2d(Unwrapped_Phase,BrainMask,'voxelsize',voxelsize,'padsize',padsize,'smvsize',smvsize);
            %   [TissuePhase,NewMask] = V_SHARP(Unwrapped_Phase,BrainMask,'voxelsize',voxelsize,'smvsize',smvsize);
            %
            % Inputs
            %   Unwrapped_Phase : 3D raw phase
            %   BrainMask       : binary mask
            %   voxelsize       : spatial resolution, default = [1 1 1]
            %   padsize         : size for padarray to increase numerical accuracy, default = [12 12 12] (2D only)
            %   smvsize         : max kernel size, default = 12
            %
            % Outputs
            %   TissuePhase     : local field
            %   NewMask         : eroded mask

            % Remove background field field as algorithm does not output it
            data = rmfield(data, 'BackgroundField');

            if dims == 2
                data.Data.LocalField = V_SHARP_2d(data.Data.TotalField, data.Data.Mask, 'voxelsize', data.VoxelSize, 'padsize', params.psize, 'smvsize', params.smvsize);

                % For some reason, 2D V-SHARP doesn't output trimmed mask. So have to rerun V-SHARP for the sole purpose of obtaining eroded mask
                data.Data.Mask = data.Data.LocalField~=0;
            else
                [data.Data.LocalField, data.Data.Mask] = V_SHARP(data.Data.TotalField, data.Data.Mask, 'voxelsize', data.VoxelSize, 'smvsize', params.smvsize);
            end

        case 'LBV'
            %%% Laplacian Boundary Value (LBV) - QSM-master/MEDI Toolbox
            %   [fL] = LBV(iFreq, Mask, matrix_size, voxel_size, tol, depth, peel, N1, N2, N3)
            %
            % Inputs
            %   iFreq       : total magnetic field
            %   Mask        : ROI
            %   matrix_size : dimension of the 3D image stack
            %   voxel_size  : dimensions of the voxels in unit of mm
            %   tol         : iteration stopping criteria on the coarsest grid, e.g. 0.01
            %   peel        : number of boundary layers to be peeled off (quick & dirty way to get rid of bad voxels at ROI boundary
            %   depth       : number of length scales (largest length scale is 2^depth*voxelsize
            %   N1          : iterations on each depth before the recursive call
            %   N2          : iterations on each depth after the recursive call
            %   N3          : iterations on the finest scale after the FMG is finished.
            %
            % Outputs
            %   fL          : local field

            if dims == 2
                error(['The LBV background removal algorithm will likely crash if run for 2D data. ' ...
                    'Please verify your system can handle it beforehand. If so, remove this error.'])
            else
                data.Data.LocalField = LBV(data.Data.TotalField, data.Data.Mask, data.Size, data.VoxelSize, params.tol, params.peel, params.depth, params.n1, params.n2, params.n3);
            end
    end

    % Save parameters to output structure
    switch flags.method
        case 'PDF - MEDI'
            BackgroundRemoval = struct('Method', flags.method, ...
                                    'Tolerance', params.tol, ...
                                'MaxIterations', params.iter, ...
                               'SpaceSelection', params.space, ...
                                      'PadSize', params.psize);

        case 'PDF - QSMmaster'
            BackgroundRemoval = struct('Method', flags.method, ...
                                'MaxIterations', params.iter);

        case 'SHARP'
            BackgroundRemoval = struct('Method', flags.method, ...
                                 'KernelRadius', params.kerrad, ...
                                 'TruncatedSVD', params.tsvd);

        case 'V-SHARP - SEPIA'
            BackgroundRemoval = struct('Method', flags.method, ...
                                  'KernelRadii', params.kerrads);

        case 'V-SHARP - STISuite'
            if dims == 2
                BackgroundRemoval = struct('Method', flags.method, ...
                                          'PadSize', params.psize, ...
                                    'MaxKernelSize', params.smvsize./data.VoxelSize(1));
            end
            data.BackgroundRemoval3D = struct('Method', flags.method, ...
                                       'MaxKernelSize', params.smvsize./data.VoxelSize(1));

        case 'LBV'
            BackgroundRemoval = struct('Method', flags.method, ...
                                    'Tolerance', params.tol, ...
                         'PeeledBoundaryLayers', params.peel, ...
                                 'LengthScales', params.depth, ...
                       'PreRecursiveIterations', params.n1, ...
                      'PostRecursiveIterations', params.n2, ...
                             'FinestIterations', params.n3);
    end

    % Assign to structure
    if dims == 2
        % Check if LBV somehow got snuck in
        if strcmp(flags.method,'LBV')
            error('LBV can only be used in 3D data!')
        else
            data.BackgroundRemoval2D = BackgroundRemoval;
        end
    elseif ~strcmp(flags.method,'V-SHARP - STISuite') % skip for special case
        data.BackgroundRemoval3D = BackgroundRemoval;
    end

    % Update flags
    if dims == 2
        data.Flags.Removed2DBackgroundField = true;
        if strcmp(flags.method,'V-SHARP - STISuite')
            data.Flags.Removed3DBackgroundField = true;
        end
    else
        data.Flags.Removed3DBackgroundField = true;
    end

    % Correct for flipped data
    data = Operations.unFlip(data, flippedFLAG);

    % Trim empty slices due to erosion
    sl1 = 1; sl2 = data.Size(3);
    for sl = 1:data.Size(3)
        if sum(data.Data.Mask(:,:,sl), "all") == 0
            if sl < data.Size(3)/2
                sl1 = sl+1;
            else
                sl2 = sl-1;
            end
        end
    end
    datafields = fieldnames(data.Data);
    for j = 1:numel(datafields)
        data.Data.(datafields{j}) = data.Data.(datafields{j})(:,:,sl1:sl2);
    end
    data.Size = size(data.Data.LocalField);
end
