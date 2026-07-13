% Dipole inversion
function [D] = perform(D, flags, algoParams)
% Inputs:
%         data: data structure
%        flags: processing structure
%   algoParams: algorithm parameters (optional)
% 
% Outputs:
%         data: data structure

    if ~isfield(flags, 'verbose') || flags.verbose
        verbose = true;
    else
        verbose = false;
    end

    % Ensure input is valid
    if ~any(strcmp(flags.method, {'TKD - MRIscm', 'TKD - SEPIA', 'TVDI', 'MEDI', 'iLSQR', 'StarQSM', 'FANSI', 'Closed Form L2', 'Closed Form', 'NDI', 'Magnitude Weighted L1'}))
        if strcmp(flags.method, 'TKD')
            TKDmethod = input('Two TKD algorithms are available: (1) SEPIA and (2) MRI susceptibility calculation methods. Which would you prefer?','s');
            if strcmpi(TKDmethod, {'SEPIA', '1', '(1)'})
                flags.method = 'TKD - SEPIA';
            elseif strcmpi(TKDmethod, {'MRI susceptibility calculation methods', 'MRIscm', '2', '(2)'})
                flags.method = 'TKD - MRIscm';
            else
                error('The entered input could not be interpreted.')
            end
        else
            error('The entered technique for dipole inversion has not been added to this pipeline.')
        end
    else
        if verbose
            fprintf('\nPerforming dipole inversion using %s...\n', flags.method)
        end
    end

    % Ensure B0 direction is positive
    [D, flippedFLAG] = Operations.Flip(D);

    % Algorithm specific parameters
    switch flags.method
        case 'TKD - MRIscm' 
            try thresh = algoParams.thresh;                 catch; thresh          = 2/3;                        end % Default 2/3   , kernel threshold 
        case 'TKD - SEPIA'
            try thresh = algoParams.thresh;                 catch; thresh          = 5/20;                       end % Default 3/20  , kernel threshold 
        case 'TVDI'
            try lambda = algoParams.lambda;                 catch; lambda          = 5e-4;                       end % Default 5e-4  , total variation regularization paramter, e.g. 5e-4
            try weight = algoParams.weight;                 catch; weight          = D.Data.WeightedMagnitude;end % Default mag   , data consistency weighting (mask or magnitude)
            try iter = algoParams.iter;                     catch; iter            = 50;                         end % Default 500   , number of NLCG iterations
            try pnorm = algoParams.pnorm;                   catch; pnorm           = 1;                          end % Default 1     , L1 or L2 regularization
        case 'MEDI'
            try lambda = algoParams.lambda;                 catch; lambda          = 200;                        end % Default 1000  , regularization parameter
            try perc = algoParams.perc;                     catch; perc            = 0.9;                        end % Default 0.9   , number of voxels considered 'edge' in L1 regularization
            try smv = algoParams.smv;                       catch;                                               end % Default 5     , approximation for Laplacian operator
            try max_iter = algoParams.max_iter;             catch; max_iter        = 10;                         end % Default 10    , number of Gauss-Newton solver iterations
            try cg_max_iter = algoParams.cg_max_iter;       catch; cg_max_iter     = 100;                        end % Default 100   , number of conjugate gradient iterations within each GN solver iteration
            try data_weighting = algoParams.data_weighting; catch; data_weighting  = 1;                          end % Default 1     , data weighting mode (0=uniform weighting, 1=SNR weighting)
            try tol_norm_ratio = algoParams.tol_norm_ratio; catch; tol_norm_ratio  = 0.1;                        end % Default 0.1   , threshold value of relative update changes
            try gpuFLAG = algoParams.gpuFLAG;               catch; gpuFLAG         = true;                       end
        case {'iLSQR', 'StarQSM'}
            try psize = algoParams.psize;                   catch; psize           = round(D.Size/10);        end % Default 1/10  , pad size
        case 'FANSI'
            try padFLAG = algoParams.padFLAG;               catch; padFLAG         = true;                       end % Default true  , padding flag
            try tgvFLAG = algoParams.tgvFLAG;               catch; tgvFLAG         = false;                      end % Default false , TV or TGV regularization
            try nonlinearFLAG = algoParams.nonlinearFLAG;   catch; nonlinearFLAG   = true;                       end % Default true  , linear or nonlinear algorithm
            try alpha = algoParams.alpha;                   catch; alpha           = 1e-5;                       end % Default 1e-4  , regularization parameter
            try maxOuterInter = algoParams.maxOuterInter;   catch; maxOuterInter   = 1000;                       end % Default 150   , maximum number of iterations
            try tol_update = algoParams.tol_update;         catch; tol_update      = 0.1;                        end % Default 0.1   , convergence limit/update ratio of the solution
            % try gradientMode = algoParams.gradientMode;     catch; gradientMode    = 1;                          end % Default 0     , gradient field used as regularization weights with TV/TGV
            try gradientMode = algoParams.gradientMode;     catch;                                               end % Default 0     , gradient field used as regularization weights with TV/TGV
            try kernelMode = algoParams.kernelMode;         catch;                                               end % Default 0     , dipole kernel formulation (not used if B0 direction is angled)
        case 'Closed Form L2'
            try padFLAG = algoParams.padFLAG;               catch; padFLAG         = true;                       end % Default true  , padding flag
            try lambda = algoParams.lambda;                 catch; lambda          = 0.015;                      end % Default 0.015 , regularization parameter
        case 'Closed Form'
            try lambda = algoParams.lambda;                 catch; lambda          = 10.^(-linspace(3,10,100));  end % Default 0.015 , regularization parameter
        case 'NDI'
            try padFLAG = algoParams.padFLAG;               catch; padFLAG         = true;                       end % Default true  , padding flag
            try alpha = algoParams.alpha;                   catch; alpha           = 1e-6;                       end % Default 1e-6  , regularization parameter
            try maxOuterInter = algoParams.maxOuterInter;   catch; maxOuterInter   = 1000;                       end % Default 1000  , maximum number of iterations
            try weight = algoParams.weight;                 catch; weight          = D.Data.WeightedMagnitude;end % Default mag   , data consistency/fidelity weighting
            try tau = algoParams.tau;                       catch; tau             = 1;                          end % Default 1     , gradient descent rate
            try precond = algoParams.precond;               catch; precond         = false;                      end % Default false , preconditioned solution for stability (start QSM with 3*weight*lfs instead of array of zeros)
            try isShowIters = algoParams.isShowIters;       catch; isShowIters     = false;                      end % Default false , verbosity flag
        case 'Magnitude Weighted L1'
            try padFLAG = algoParams.padFLAG;               catch; padFLAG         = true;                       end % Default true  , padding flag
            try lambda_L1 = algoParams.lambda_L1;           catch; lambda_L1       = 1e-4;                       end % Default 1e-4  , L1 regularization parameter
            try lambda_L2 = algoParams.lambda_L2;           catch; lambda_L2       = 1e-4;                       end % Default 1e-4  , L2 regularization parameter
    end

    % Run algorithm
    switch flags.method   
        case 'TKD - MRIscm'
            %%% Truncated k-space Division (TKD) - MRI susceptiblity calculation methods
            %   susc = TKD(Params)
            %
            % Inputs
            %   Params      :
            %                 Params.FieldMap    : input image in ppm
            %                 Params.Mask        : binary tissue mask (default = 1)
            %                 Params.Threshold   : kernel threshold (default = 2/3)
            %                 Params.Resolution  : image resolution vector (dx,dy,dz) in mm (default = 1 mm isotropic)
            %                 Params.B0direction : 3-element unit vector aligned with B0 (default = [0,0,1])
            %
            % Ouputs
            %   susc        : output susceptibility map in ppm
    
            % Convert to ppm
            Parameters.FieldMap = D.Data.LocalField./(42.5774780505984*D.B0); % in ppm
    
            % Add other inputs to structure
            Parameters.Mask = D.Data.Mask; % default = 1
            Parameters.Threshold = thresh; % delfault = 2/3
            Parameters.Resolution = D.VoxelSize; % default = [1,1,1]
            Parameters.B0direction = D.B0Direction; % default = [0,0,1]
    
            D.Data.SusceptibilityMap = TKD(Parameters);

        case 'TKD - SEPIA'
            %%% Truncated k-space Division (TKD) - SEPIA
            %   susc = TKD(Params)
            %
            % Inputs
            %   localField      : local field perturbatios
            %   mask            : user-defined mask
            %   matrixSize      : image matrix size
            %   voxelSize       : spatial resolution of image 
            %   varargin        : flags with
            %       'threshold'     -   threshold for k-space inversion 
            %
            % Ouputs
            %   susc        : output susceptibility map in ppm
       
            % Ensure threshold is less than 2/3
            if thresh >= 2/3
                thresh = input('The threshold must be less than 2/3. Please enter in a new threshold:', 's');
                if isnan(str2double(thresh))
                    error('The given response is not numeric. Dipole inversion has been cancelled.')
                elseif thresh >= 2/3
                    error('The given response is still greater than or equal to 2/3. Dipole inversion has been cancelled.')
                end
            end

            D.Data.SusceptibilityMap = qsmTKD(D.Data.LocalField, D.Data.Mask, D.Size, D.VoxelSize, 'threshold', thresh, 'b0dir', D.B0Direction);

            % Convert to ppm
            D.Data.SusceptibilityMap = D.Data.SusceptibilityMap./(42.5774780505984*D.B0);

        case 'TVDI'
            %%% Total Variation Dipole Inversion (TVDI) - QSM-master
            %   [sus, res] = TVDI(lfs, mask, vox, tv_reg, weights, z_prjs, itnlim, pnorm)
            %
            % Inputs
            %   lfs         : local field shift (field perturbation map)
            %   mask        : binary mask defining ROI
            %   vox         : voxel size, e.g. [1 1 1] for isotropic resolution
            %   tv_reg      : Total Variation regularization paramter, e.g. 5e-4
            %   weights     : weights for the data consistancy term, e.g. mask or magnitude
            %  *z_prjs      : normal vector of the imaging plane, default = [0,0,1]
            %  *itnlim      : interation numbers of nlcg, default = 500
            %  *pnorm       : L1 or L2 norm regularization, default = 1
            %
            % Outputs
            %   sus         : susceptibility distribution after dipole inversion
            %   res         : residual field after QSM fitting
            %   kernel      : dipole kernel
    
            % Convert to ppm
            lfs = D.Data.LocalField./(42.5774780505984*D.B0);

            [D.Data.SusceptibilityMap, ~, kernel] = tvdi(lfs, D.Data.Mask, D.VoxelSize, lambda, weight, D.B0Direction, iter, pnorm);

            % Calculate cost
            sus = padarray(D.Data.SusceptibilityMap,[0 0 20]);
            lfs = padarray(lfs,[0 0 20]);
            mask = padarray(D.Data.Mask,[0 0 20]);
            [data_cost, reg_cost] = compute_costs(sus.*mask, lfs.*mask, kernel);

        case 'MEDI'
            %%% Morphology Enabled Dipole Inversion (MEDI) - MEDI toolbox
            %   [sus, cost_reg, cost_data, resultsfile] = MEDI_L1(filename, lambda, data_weighting, merit, smv, DEBUG, lambda_CSF, percentage)
            %
            % Inputs (name-value pairs)
            %   filename           : name of file containing inputs, default = "RFD.mat"
            %   lambda             : regularization parameter, default = 1000
            %   data_weighting     : data weighting mode (0=uniform weighting, 1=SNR weighting), default = 1
            %   merit              : turn on model error reduction through iterative tuning (not name-value pair)
            %   smv                : turn on smv operation & set smv radius, default = 5 (not including this input turns off smv operation)
            %   DEBUG              : turn on debug mode (not name-value pair)
            %   lambda_CSF         : automatic zero reference (MEDI+0) also require Mask_CSF in RDF.mat
            %   percentage         : number of voxels regarded as "edge" in ROI during L1 regularization, default = 0.9
            %   *Note: also requires "RDF.mat" saved to path & including the following
            %      iFreq           : total field shift (unwrapped phase), units = rad/echo
            %      RDF             : local field shift, units = rad/echo
            %      N_std           : estimated noise standard deviation of raw phase
            %      iMag            : raw magnitude
            %      Mask            : binary mask
            %      matrix_size     : matrix size
            %      voxel_size      : voxel size, units = mm
            %      delta_TE        : echo spacing/TE for multi/single-echo, units = sec
            %      CF              : center frequency, units = Hz
            %      B0_dir          : magnet direction
            %
            % Outputs
            %   sus                : susceptibility distribution, units = ppm
            %   cost_reg           : cost of regularization term
            %   cost_data          : cost of data fidelity term
            %   results            : results file, containing iMag, Mask, QSM, RDF, and summary
    
            % Pad data
            pval = 10;
            tfs = padarray(D.Data.TotalField, [0 0 pval]);
            nstd = padarray(D.Data.NoiseSTD, [0 0 pval]);
            mag = padarray(D.Data.WeightedMagnitude, [0 0 pval]);
            lfs = padarray(D.Data.LocalField, [0 0 pval]);
            msk = padarray(D.Data.Mask, [0 0 pval]);

            % Rename variables for algorithm (and transfer to GPU if requested)
            data_struct = struct('iFreq', tfs.*(2*pi*D.deltaTE), ...
                                   'RDF', lfs.*(2*pi*D.deltaTE), ...
                                 'N_std', nstd, ...
                                  'iMag', mag, ...
                                  'Mask', msk, ...
                           'matrix_size', size(lfs), ...
                            'voxel_size', D.VoxelSize, ...
                              'delta_TE', D.deltaTE, ...
                                    'CF', D.F0*1e6, ...
                                'B0_dir', D.B0Direction);
            if gpuFLAG
                nfields = fieldnames(data_struct);
                for i = 1:numel(nfields)
                    data_struct.(nfields{i}) = gpuArray(data_struct.(nfields{i}));
                end

                lambda = gpuArray(lambda);
                if exist('smv', 'var')
                    smv = gpuArray(smv);
                end
                perc = gpuArray(perc);
                max_iter = gpuArray(max_iter);
                cg_max_iter = gpuArray(cg_max_iter);
                tol_norm_ratio = gpuArray(tol_norm_ratio);
            end
    
            if exist('smv', 'var')
                [D.Data.SusceptibilityMap, reg_cost, data_cost] = MEDI_L1('data_struct', data_struct, 'lambda', lambda, 'data_weighting', data_weighting, 'merit', 'smv', smv, 'percentage', perc, 'max_iter', max_iter, 'cg_max_iter', cg_max_iter, 'tol_norm_ratio', tol_norm_ratio, 'verbose', verbose);
            else
                [D.Data.SusceptibilityMap, reg_cost, data_cost] = MEDI_L1('data_struct', data_struct, 'lambda', lambda, 'data_weighting', data_weighting, 'merit',             'percentage', perc, 'max_iter', max_iter, 'cg_max_iter', cg_max_iter, 'tol_norm_ratio', tol_norm_ratio, 'verbose', verbose);
            end

            % Remove padding
            D.Data.SusceptibilityMap = D.Data.SusceptibilityMap(:,:,pval+1:end-pval);

            if gpuFLAG
                reg_cost = gather(reg_cost);
                data_cost = gather(data_cost);
            end
            reg_cost = real(reg_cost); 
            data_cost = real(data_cost);
            
            % Remove zeros from cost data
            reg_cost = reg_cost(reg_cost~=0);
            data_cost = data_cost(data_cost~=0);

            % Isolate the last iteration
            reg_cost = reg_cost(end);
            data_cost = data_cost(end);

        case 'iLSQR'
            %%% Initial susceptibility estimation using Sparse Linear Equation & Least-Squares (iLSQR) - STI Suite
            %   sus = QSM_iLSQR(lfs, m, 'TE', TE, 'B0', B0, 'H', H, 'padsize', psize, 'voxelsize', vsize);
            %
            % Inputs
            %   lfs          : tissue phase
            %   m            : binary mask of ROI
            %   TE           : echo time(s) in s
            %   B0           : magnet field strength in T
            %   H            : the field direction, e.g. H=[0 0 1];
            %   vsize        : spatial resolution
            %   psize        : size for padarray to increase numerical accuracy
            %
            % Outputs
            %       sus      : QSM images
    
            sus = QSM_iLSQR(D.Data.LocalField, D.Data.Mask, 'TE', D.TE, 'B0', D.B0, 'H', D.B0Direction, 'padsize', psize, 'voxelsize', D.VoxelSize);
    
            % Convert from ppb to ppm
            D.Data.SusceptibilityMap = sus./1000;
    
        case 'StarQSM'
            %%% Streaking Artifact Reduction for QSM (StarQSM) - STI Suite
            %  sus = QSM_star(lfs, m, 'TE', TE, 'B0', B0, 'H', H, 'padsize', psize, 'voxelsize', vsize);
            %
            % Inputs
            %   lfs          : tissue phase
            %   m            : binary mask of ROI
            %   TE           : echo time(s) in s
            %   B0           : magnet field strength in T
            %   H            : the field direction, e.g. H=[0 0 1];
            %   vsize        : spatial resolution
            %   psize        : size for padarray to increase numerical accuracy
            %
            % Outputs
            %       sus      : QSM images
    
            sus = QSM_star(D.Data.LocalField, D.Data.Mask, 'TE', D.TE, 'B0', D.B0, 'H', D.B0Direction, 'padsize', psize, 'voxelsize', D.VoxelSize);

            % Convert from ppb to ppm
            D.Data.SusceptibilityMap = sus./1000;

        case 'FANSI'
            %%% FAst Nonlinear Susceptibility Inversion (FANSI) - FANSI Toolbox
            %%% Nonlinear QSM and Total Variation regularization (nlTV)
            %%% Nonlinear QSM and Total Generalized Variation regularization (nlGTV)
            %%% Linear QSM and Total Variation regularization (wTV)
            %%% Linear QSM and Total Generalized Variation regularization (wTGV)
            %
            % Inputs
            %   phase                   : local field map data
            %   magn                    : magnitude data
            %   alpha                   : gradient L1 penalty, regularization weight
            %   options.isNonlinear     : linear or nonlinear algorithm? (default = true, i.e nonlinear method)
            %   options.isTGV           : TV or TGV regularization? (default = false, i.e. TV regularization)
            %   options.mu              : ADMM Lagrange multiplier (recommended = 100*alpha1)
            %   options.iterations      : maximum number of iterations (recommended = 150)
            %   options.update          : convergence limit, update ratio of the solution (recommended = 0.1)
            %   options.weight          : data fidelity spatially variable weight (recommended = magnitude_data).
            %   options.isPrecond       : preconditionate solution by smart initialization (default = true)
            %   options.isGPU           : GPU acceleration (default = true)
            %   options.gradientMode    : this creates a gradient field to be used as regularization weights for TV/TGV, with:
            %                               0 to use the vector field.
            %                               1 for the L1 norm, and
            %                               2 for the L2 norm
            %   options.noise           : noise standard deviation in the complex signal, required for the regularization weight based on the gradients.
            %   options.voxelSize       : Spatial resolution vector, in mm, or normalized (mean = 1).
            %   options.B0_dir          : main field direction, e.g. [0 0 1] (only for continuous kernel)
            %   options.kernelMode      : dipole kernel formulation, with:
            %                               0 for the continuous kernel proposed by Salomir, et al. 2003.
            %                               1 for the discrete kernel proposed by Milovic, et al. 2017.
            %                               2 for the Integrated Green function proposed by Jenkinson, et al. 2004
            %   options.mu2             : fidelity consistency weight (ADMM weight, recommended value = 1.0)
            %   options.alpha0          : curvature L1 penalty, regularization weight (recommended = 2*alpha1)
            %   options.mu0             : curvature consistency weight (ADMM weight, recommended = 2*mu1)
            %   options.regweight       : regularization spatially variable weight.
            %
            % Outputs
            %   out.x                   : calculated susceptibility map
            %   out.iter                : number of iterations needed
            %   out.time                : elapsed time (excluding pre-calculations)
            %   out.totalTime           : total elapsed time (including pre-calculations)

            if padFLAG
                % Set padding amount
                p = 4;

                % Pad data
                msk_pad = padarray(D.Data.Mask, floor(D.Size/p));
                lfs_pad = padarray(D.Data.LocalField, floor(D.Size/p));
                mag_pad = padarray(D.Data.WeightedMagnitude, floor(D.Size/p));
            else
                msk_pad = D.Data.Mask;
                lfs_pad = D.Data.LocalField;
                mag_pad = D.Data.WeightedMagnitude;
            end
            n = D.Size;
            N = size(msk_pad);

            % Create kernel (from FANSI to account for angled slices)
            kernel = dipole_kernel_angulated(N, D.VoxelSize, D.B0Direction);

            % Run main algorithm
            mu1 = 100*alpha;
            options = struct('isTGV', tgvFLAG, ...
                       'isNonlinear', nonlinearFLAG, ...
                         'voxelSize', D.VoxelSize, ...
                            'B0_dir', D.B0Direction, ...
                                'mu', mu1, ...
                        'iterations', maxOuterInter, ...
                            'update', tol_update);
            if exist('gradientMode', 'var')
                options.gradientMode = gradientMode;
            end
            out = FANSI(lfs_pad, mag_pad, alpha, options);
            D.Data.SusceptibilityMap = out.x;

            % Calculate cost
            [data_cost, reg_cost] = compute_costs(out.x.*msk_pad, lfs_pad.*msk_pad, kernel);

            % Crop output
            if padFLAG
                D.Data.SusceptibilityMap = out.x(1+floor(n(1)/p):floor((1+1/p)*n(1)), 1+floor(n(2)/p):floor((1+1/p)*n(2)), 1+floor(n(3)/p):floor((1+1/p)*n(3))).*D.Data.Mask;
            else
                D.Data.SusceptibilityMap = out.x;
            end

        case 'Closed Form L2'
            if padFLAG
                % Pad data
                msk_pad = padarray(D.Data.Mask, D.Size/2);
                nfm_pad = padarray(D.Data.LocalField, D.Size/2);
            end
            n = D.Size;
            N = size(msk_pad);

            % Create kernel (from FANSI to account for angled slices)
            kernel = dipole_kernel_angulated(N, D.VoxelSize, D.B0Direction);

            % Closed form QSM solution
            [k2,k1,k3] = meshgrid(0:N(2)-1, 0:N(1)-1, 0:N(3)-1);
            fdx = (1 - exp(-2*pi*1i*k1/N(1)));
            fdy = (1 - exp(-2*pi*1i*k2/N(2)));
            fdz = (1 - exp(-2*pi*1i*k3/N(3)));

            D_reg = kernel./(eps + abs(kernel).^2 + lambda*(abs(fdx).^2 + abs(fdy).^2 + abs(fdz).^2));
            D_regx = real(ifftn(D_reg .* fftn(nfm_pad)));

            % Calculate cost
            [data_cost, reg_cost] = compute_costs(D_regx.*msk_pad, nfm_pad.*msk_pad, kernel);

            % Crop data
            D.Data.SusceptibilityMap = D_regx(1+n(1)/2:1.5*n(1), 1+n(2)/2:1.5*n(2), 1+n(3)/2:1.5*n(3)).*D.Data.Mask;

        case 'Closed Form'
            % Convert to ppm
            lfs = D.Data.LocalField./(42.5774780505984*D.B0);

            [D.Data.SusceptibilityMap, lambda] = reconstructSusceptibility(lfs, D.Data.Mask, D.VoxelSize, lambda, D.B0Direction);

        case 'NDI'
            % Convert to rad from Hz (rad = Hz*2pi*sec)
            lfs = D.Data.LocalField.*(2*pi*D.deltaTE);

            % Pad data
            if padFLAG
                xpad = 4;
                ypad = 4;
                zpad = 60;
                lfs = padarray(lfs, [xpad, ypad, zpad], 0, 'both');
                weight = padarray(weight, [xpad, ypad, zpad], 0, 'both');
                msk = padarray(D.Data.Mask, [xpad, ypad, zpad], 0, 'both');
            end

            % Create kernel
            kernel = dipole_kernel_angulated(size(lfs), D.VoxelSize, D.B0Direction);

            % Run main algorithm
            options = struct('input', lfs, ...
                                 'K', kernel, ...
                             'alpha', alpha, ...
                     'maxOuterInter', maxOuterInter, ...
                         'voxelSize', D.VoxelSize, ...
                            'weight', weight, ...
                               'tau', tau, ...
                           'precond', precond, ...
                       'isShowIters', isShowIters, ...
                             'isGPU', true);
            out = ndi_auto(options);

            % Demean data
            out.x = fieldmap_demean(out.x, msk).*msk;

            % Calculate cost
            [data_cost, reg_cost] = compute_costs(out.x.*msk, lfs.*msk, kernel);

            % Transfer to output structure
            D.Data.SusceptibilityMap = out.x(xpad+1:end-xpad,ypad+1:end-ypad,zpad+1:end-zpad);
            weight = weight(xpad+1:end-xpad,ypad+1:end-ypad,zpad+1:end-zpad);

        case 'Magnitude Weighted L1'
            % Convert to rad from Hz (rad = Hz*2pi*sec)
            lfs = D.Data.LocalField.*(2*pi*D.deltaTE);

            % Pad data
            if padFLAG
                xpad = 32;
                ypad = 32;
                zpad = 8;
                lfs = padarray(lfs, [xpad, ypad, zpad], 0, 'both');
                msk = padarray(D.Data.Mask, [xpad, ypad, zpad], 0, 'both');
                mag = padarray(D.Data.WeightedMagnitude, [xpad, ypad, zpad], 0, 'both');
            end

            % Run algorithm
            D.Data.SusceptibilityMap = L1_magnWeighted_reconstruction(lfs, mag, msk, D.VoxelSize, lambda_L1, lambda_L2, D.B0Direction);

            % Transfer to output structure
            D.Data.SusceptibilityMap = D.Data.SusceptibilityMap(xpad+1:end-xpad,ypad+1:end-ypad,zpad+1:end-zpad);
    end

    % Save parameters to output structure
    switch flags.method
        case {'TKD - MRIscm', 'TKD - SEPIA'}
            D.DipoleInversion = struct('Method', flags.method, ...
                                       'Threshold', thresh);
        case 'TVDI'
            D.DipoleInversion = struct('Method', flags.method, ...
                         'RegularizationParameter', lambda, ...
                                   'MaxIterations', iter, ...
                             'LRegularizationType', pnorm, ...
                             'DataConsistencyCost', data_cost, ...
                              'RegularizationCost', reg_cost);
            if all(weight==D.Data.Mask,"all")
                D.DipoleInversion.DataConsistencyWeighting = 'Mask';
            else
                D.DipoleInversion.DataConsistencyWeighting = 'Magnitude';
            end
        case 'MEDI'
            if gpuFLAG
                D.Data.SusceptibilityMap = gather(D.Data.SusceptibilityMap);
                lambda = gather(lambda);
                perc = gather(perc);
                max_iter = gather(max_iter);
                cg_max_iter = gather(cg_max_iter);
                tol_norm_ratio = gather(tol_norm_ratio);
                if exist('smv', 'var')
                    smv = gather(smv);
                end
            end

            D.DipoleInversion = struct('Method', flags.method, ...
                         'RegularizationParameter', lambda, ...
                                     'SNRWeighted', data_weighting, ...
                                           'MERIT', true, ...
                                        'EdgeSize', perc, ...
                               'AutoZeroReference', false, ...
                                   'MaxIterations', max_iter, ...
                                 'MaxCGIterations', cg_max_iter, ...
                                       'Threshold', tol_norm_ratio, ...
                             'DataConsistencyCost', data_cost, ...
                              'RegularizationCost', reg_cost);
            if exist('smv', 'var')
                D.DipoleInversion.SMV = true;
                D.DipoleInversion.SMVRadius = smv;
            else
                D.DipoleInversion.SMV = false;
            end
        case {'iLSQR',  'Star'}
            D.DipoleInversion = struct('Method', flags.method, ...
                                         'PadSize', psize);
        case 'FANSI'
            D.DipoleInversion = struct('Method', flags.method, ...
                                       'PaddedFOV', padFLAG, ...
                                     'Generalized', tgvFLAG, ...
                                       'Nonlinear', nonlinearFLAG, ...
                                   'MaxIterations', maxOuterInter, ...
                                'LimitUpdateRatio', tol_update, ...
                         'RegularizationParameter', alpha, ...
                                             'mu1', mu1, ...
                             'DataConsistencyCost', data_cost, ...
                              'RegularizationCost', reg_cost);
            if exist('gradientMode', 'var')
                D.DipoleInversion.GradientMode = gradientMode;
            end
            if exist('kernelMode', 'var')
                D.DipoleInversion.KernelMode = kernelMode;
            end
        case 'Closed Form L2'
            D.DipoleInversion = struct('Method', flags.method, ...
                                       'PaddedFOV', padFLAG, ...
                         'RegularizationParameter', lambda, ...
                             'DataConsistencyCost', data_cost, ...
                              'RegularizationCost', reg_cost);
        case 'NDI'
            D.DipoleInversion = struct('Method', flags.method, ...
                                       'PaddedFOV', padFLAG, ...
                                   'MaxIterations', maxOuterInter, ...
                         'RegularizationParameter', alpha, ...
                             'GradientDescentRate', tau, ...
                             'DataConsistencyCost', data_cost, ...
                              'RegularizationCost', reg_cost);

            if all(weight==D.Data.Mask,"all")
                D.DipoleInversion.DataConsistencyWeighting = 'Mask';
            else
                D.DipoleInversion.DataConsistencyWeighting = 'Magnitude';
            end
    end

    % Correct for flipped data
    D = Operations.unFlip(D, flippedFLAG);

    % Update flags
    D.Flags.InvertedDipole = true;
end