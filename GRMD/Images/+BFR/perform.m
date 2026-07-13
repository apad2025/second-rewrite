% Remove background field
function D = perform(D, flags, dims, verboseFLAG) 
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
    [D, flippedFLAG] = Operations.Flip(D);

    % Check slice direction
    if ~isequal(D.B0Direction,[0 0 1])
        disp('This is angled slicing');
        disp(D.B0Direction);
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
            params = struct('psize', D.Size, ...                        Default [12,12,12]  , padding size
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
    D.Data.LocalField = zeros(D.Size, 'like', D.Data.TotalField);
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

            D.Data.BackgroundField = zeros(D.Size, 'like', D.Data.TotalField);

            % Check 2D or 3D
            if dims == 2
                for sl = 1:D.Size(3)
                    if verboseFLAG; fprintf('\nSlice %i...', sl); end
                    [D.Data.LocalField(:,:,sl), D.Data.BackgroundField(:,:,sl)] = PDF(D.Data.TotalField(:,:,sl), D.Data.NoiseSTD(:,:,sl), D.Data.Mask(:,:,sl), [D.Size(1), D.Size(2)], D.VoxelSize, D.B0Direction, params.tol, params.iter, params.space, params.psize);
                end
            else
                [D.Data.LocalField, D.Data.BackgroundField] = PDF(D.Data.TotalField, D.Data.NoiseSTD, D.Data.Mask, D.Size, D.VoxelSize, D.B0Direction, params.tol, params.iter, params.space, params.psize);
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

            D.Data.BackgroundField = zeros(D.Size, 'like', D.Data.TotalField);

            if dims == 2
                for sl = 1:D.Size(3)
                    if verboseFLAG; fprintf('\nSlice %i...', sl); end
                    [D.Data.LocalField(:,:,sl), ~, D.Data.BackgroundField(:,:,sl)] = projectionontodipolefields(D.Data.TotalField(:,:,sl), D.Data.Mask(:,:,sl), D.VoxelSize, D.Data.WeightedMagnitude(:,:,sl), D.B0Direction, params.iter);
                end
            else
                [D.Data.LocalField, ~, D.Data.BackgroundField] = projectionontodipolefields(D.Data.TotalField, D.Data.Mask, D.VoxelSize, D.Data.WeightedMagnitude, D.B0Direction, params.iter);
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
                for sl = 1:D.Size(3)
                    if verboseFLAG; fprintf('\nSlice %i...', sl); end
                    [D.Data.LocalField(:,:,sl), D.Data.Mask(:,:,sl)] = sharp(D.Data.TotalField(:,:,sl), D.Data.Mask(:,:,sl), D.VoxelSize, params.kerrad, params.tsvd);
                end
            else
                [D.Data.LocalField, D.Data.Mask] = sharp(D.Data.TotalField, D.Data.Mask, D.VoxelSize, params.kerrad, params.tsvd);
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
                for sl = 1:D.Size(3)
                    if verboseFLAG; reverseStr = UpdatePercent(100*sl/D.Size(3), reverseStr); end
                    [D.Data.LocalField(:,:,sl), D.Data.Mask(:,:,sl)] = BKGRemovalVSHARP_2D(D.Data.TotalField(:,:,sl), D.Data.Mask(:,:,sl), [D.Size(1), D.Size(2)], 'radius', params.kerrads);
                end
                fprintf('\tDone!\n');
            else
                [D.Data.LocalField, D.Data.Mask] = BKGRemovalVSHARP(D.Data.TotalField, D.Data.Mask, D.Size, 'radius', params.kerrads);
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
            D = rmfield(D, 'BackgroundField');

            if dims == 2
                D.Data.LocalField = V_SHARP_2d(D.Data.TotalField, D.Data.Mask, 'voxelsize', D.VoxelSize, 'padsize', params.psize, 'smvsize', params.smvsize);

                % For some reason, 2D V-SHARP doesn't output trimmed mask. So have to rerun V-SHARP for the sole purpose of obtaining eroded mask
                D.Data.Mask = D.Data.LocalField~=0;
            else
                [D.Data.LocalField, D.Data.Mask] = V_SHARP(D.Data.TotalField, D.Data.Mask, 'voxelsize', D.VoxelSize, 'smvsize', params.smvsize);
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
                D.Data.LocalField = LBV(D.Data.TotalField, D.Data.Mask, D.Size, D.VoxelSize, params.tol, params.peel, params.depth, params.n1, params.n2, params.n3);
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
                                    'MaxKernelSize', params.smvsize./D.VoxelSize(1));
            end
            D.BackgroundRemoval3D = struct('Method', flags.method, ...
                                       'MaxKernelSize', params.smvsize./D.VoxelSize(1));

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
            D.BackgroundRemoval2D = BackgroundRemoval;
        end
    elseif ~strcmp(flags.method,'V-SHARP - STISuite') % skip for special case
        D.BackgroundRemoval3D = BackgroundRemoval;
    end

    % Update flags
    if dims == 2
        D.Flags.Removed2DBackgroundField = true;
        if strcmp(flags.method,'V-SHARP - STISuite')
            D.Flags.Removed3DBackgroundField = true;
        end
    else
        D.Flags.Removed3DBackgroundField = true;
    end

    % Correct for flipped data
    D = Operations.unFlip(D, flippedFLAG);

    % Trim empty slices due to erosion
    sl1 = 1; sl2 = D.Size(3);
    for sl = 1:D.Size(3)
        if sum(D.Data.Mask(:,:,sl), "all") == 0
            if sl < D.Size(3)/2
                sl1 = sl+1;
            else
                sl2 = sl-1;
            end
        end
    end
    datafields = fieldnames(D.Data);
    for j = 1:numel(datafields)
        D.Data.(datafields{j}) = D.Data.(datafields{j})(:,:,sl1:sl2);
    end
    D.Size = size(D.Data.LocalField);
end
