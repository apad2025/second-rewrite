%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code: Function to perform fat-water separation using a modified version
% of the graph-cut method (Hernando et al. doi: 10.1002/mrm.22177)
% Copyright Jorge Campos 2025 - MIT License
%
% This code uses tools from the ISMRM fat-water toolbox. For this code to
% work, make sure to include the toolbox in matlab's directory.
% Phantom data was acquired in a 3T Philips scanner. Data was acquired
% enabling 'Delayed Reconstruction'

% This code uses an phase unwrapping (used only to display some results
% without phase wraps). Performing phase unwrapping is optional, but if
% used, include in the code directory the function
% 'qualityGuidedUnwrapping' from ortier and Levesque, DOI 10.1002/mrm.26989, 2017, https://gitlab.com/veronique_fortier/Quality_guided_unwrapping

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outParams = Function_Bipolar_GC(imDataParams, algoParams , vec_slices, VERBOSE)

% Check validity of params, and set default algorithm parameters if not provided
[validParams,algoParams] = checkParamsAndSetDefaults_graphcut_Bipolar_GC(imDataParams,algoParams,vec_slices);
if validParams==0
    disp('Exiting -- data not processed');
    outParams = [];
    return
end

% if nargin < 4 || algoParams.parallel
if nargin < 4
    VERBOSE = false;
end

% The size of the input signal needs to be [3D 1 #echoes], the ones is to
% represent that the signal corresponds to a single coil (or coil combined
% image) or [3D #coils #echoes] if coil combination is performed before the
% fat-water separation
imDataParams.deltaF = [0; algoParams.gyro*(algoParams.species(2).frequency(:) - algoParams.species(1).frequency(1))*(imDataParams.FieldStrength)];

%% Signal for the acquisition with bipolar gradients and without any correction

% If PrecessionIsClockwise=1 the data is loaded without any modification.
% If PrecessionIsClockwise=-1 data is changed to follow a consistent
% convention for the sign of the phase:
% Complex conjugat of data is used and code set imDataParams.PrecessionIsClockwise = 1
if imDataParams.PrecessionIsClockwise < 0
    imDataParams.images = conj(imDataParams.images);
    imDataParams.PrecessionIsClockwise = 1;
end

input_signal_bipolar_RO = imDataParams.images;

%% Matrix size and number of voxels

matrix_size = size(input_signal_bipolar_RO);
nSlices = numel(vec_slices);
numvox = prod([matrix_size(1:2) nSlices]); % only the slices in vec_slices are corrected

%% Number of echoes

numte = matrix_size(5);

%% Echo times for the full bipolar dataset

TEs_bipolar = imDataParams.TE;

%% Binary mask for fat water separation

mask = imDataParams.mask_fwseparation;

% Generate logical mask same size as data
if isscalar(mask) && mask==1
    mask = true(matrix_size(1:3));
    imDataParams.mask_fwseparation = mask;
end

%% Memory allocation

% Field maps and R2 star maps
Water_GC_odd = zeros(matrix_size(1:3));
% residual = zeros([algoParams.NUM_FMS, matrix_size(1:3)]);
if algoParams.parallel
    Water_GC_odd = cell(matrix_size(3),1);
    for kk = vec_slices
        Water_GC_odd{kk} = zeros(matrix_size(1:2));
    % residual = squeeze(mat2cell(residual, algoParams.NUM_FMS, matrix_size(1), matrix_size(2), ones([matrix_size(3),1])));
    end
end

Fat_GC_odd = Water_GC_odd;

Water_GC_even = Water_GC_odd;
Fat_GC_even = Water_GC_odd;

FieldMap_DualGC = Water_GC_odd;
R2_DualGC = Water_GC_odd;

Water_bipolar = Water_GC_odd;
Fat_bipolar = Water_GC_odd;
FieldMap_bipolar = Water_GC_odd;
R2_bipolar = Water_GC_odd;

%% Parallel pool sized to the SLURM allocation
if algoParams.parallel
    % Use one fewer worker than allocated cores in SLURM script
    nWorkers = 25;

    delete(gcp('nocreate'));
    c = parcluster('local');

    % fix queue issue
    jobid = getenv('SLURM_JOB_ID');
    if ~isempty(jobid)
        storage = fullfile(tempdir, ['matlab_pool_' jobid]);
        if ~exist(storage, 'dir')
            mkdir(storage);
        end
        c.JobStorageLocation = storage;
    end
    c.NumWorkers = max(c.NumWorkers, nWorkers);
    % Must be process pool because graph-cut solver calls max_flow MEX
    % which is unsupported by threads
    parpool(c, nWorkers);
end

%% Graph-cut fat-water separation for odd and even echoes
if VERBOSE
    fprintf('\nFat-water separation for odd and even echo datasets slice ');
end

if algoParams.parallel
    % assign only the relevant indexed parameters to each tmp worker
    imDataParams_cell = Water_GC_odd;
    for kk = vec_slices
        imDataParams_cell{kk} = imDataParams;
        imDataParams_cell{kk}.images = imDataParams.images(:,:,kk,:,:);
        imDataParams_cell{kk}.mask = mask(:,:,kk);
        imDataParams_cell{kk}.sliceofint = kk;
    end

    parfor kk = vec_slices
        outParams_GC = Function_i2cm1i_3pluspoint_hernando_Bipolar_GC(imDataParams_cell{kk}, algoParams, VERBOSE);
    
        Water_GC_odd{kk} = outParams_GC.species(1).amps;
        Fat_GC_odd{kk} = outParams_GC.species(2).amps;
    
        Water_GC_even{kk} = (outParams_GC.species(3).amps);
        Fat_GC_even{kk} = (outParams_GC.species(4).amps);
    
        FieldMap_DualGC{kk} = outParams_GC.fieldmap;
        R2_DualGC{kk} = outParams_GC.r2starmap;
        % residual{kk} = outParams_GC.residual;
    end

    % Extract from cell arrays
    Water_GC_odd_cell = Water_GC_odd; Water_GC_odd = zeros(matrix_size(1:3));
    Fat_GC_odd_cell = Fat_GC_odd; Fat_GC_odd = Water_GC_odd;
    Water_GC_even_cell = Water_GC_even; Water_GC_even = Water_GC_odd;
    Fat_GC_even_cell = Fat_GC_even; Fat_GC_even = Water_GC_odd;
    FieldMap_DualGC_cell = FieldMap_DualGC; FieldMap_DualGC = Water_GC_odd;
    R2_DualGC_cell = R2_DualGC; R2_DualGC = Water_GC_odd;
    for kk = vec_slices
        Water_GC_odd(:,:,kk) = Water_GC_odd_cell{kk};
        Fat_GC_odd(:,:,kk) = Fat_GC_odd_cell{kk};
        Water_GC_even(:,:,kk) = Water_GC_even_cell{kk};
        Fat_GC_even(:,:,kk) = Fat_GC_even_cell{kk};
        FieldMap_DualGC(:,:,kk) = FieldMap_DualGC_cell{kk};
        R2_DualGC(:,:,kk) = R2_DualGC_cell{kk};
    end
    clear Water_GC_odd_cell Fat_GC_odd_cell Water_GC_even_cell Fat_GC_even_cell FieldMap_DualGC_cell R2_DualGC_cell

else
    for kk = vec_slices

        if VERBOSE
            fprintf([num2str(kk) ', ']);
        end

        imDataParams.sliceofint = kk;
        outParams_GC = Function_i2cm1i_3pluspoint_hernando_Bipolar_GC(imDataParams, algoParams, VERBOSE);

        Water_GC_odd(:,:,kk) = outParams_GC.species(1).amps;
        Fat_GC_odd(:,:,kk) = outParams_GC.species(2).amps;

        Water_GC_even(:,:,kk) = (outParams_GC.species(3).amps);
        Fat_GC_even(:,:,kk) = (outParams_GC.species(4).amps);

        FieldMap_DualGC(:,:,kk) = outParams_GC.fieldmap;
        R2_DualGC(:,:,kk) = outParams_GC.r2starmap;
        % residual(:,:,:,kk) = outParams_GC.residual;
    end
end

%% Fat and water fraction masks generation

slice_image = algoParams.slice_image;

[c_ff, c_wf] = Function_Fat_Quantification_Bipolar_GC( Fat_GC_odd, Water_GC_odd);

c_ff = c_ff.*mask(:,:,:);
c_ff(c_ff>=1) = 1;
c_ff(c_ff<=0) = 0;

c_wf = c_wf.*mask(:,:,:);
c_wf(c_wf>=1) = 1;
c_wf(c_wf<=0) = 0;

if algoParams.plot_debug
% Plot to show the initial estimation of the fat fraction
    
    figure(101)
   
    if algoParams.crameri_colormap 
        colormap_plot = crameri('-lajolla');
    else
        colormap_plot = colormap("hot");
    end

    fp(1) = subplot(1,2,1);
    fig_imag = c_wf(:,:,slice_image);
    imagesc(fig_imag)
    axis image
    axis off
    clim([0 1])
    colormap(fp(1),colormap_plot);
    title('Initial Water Map')

    fp(2) = subplot(1,2,2);
    fig_imag = c_ff(:,:,slice_image);
    imagesc(fig_imag)
    axis image
    axis off
    clim([0 1])
    colormap(fp(2),colormap_plot);
    title('Initial Fat Map')
end

%% Fat and water mask (after thresholding using a specific values included in structure params)

ff = zeros(size(c_ff));
wf = ff;

% Set SNR threshold
mag = squeeze(sqrt(sum(abs(input_signal_bipolar_RO).^2, 5)));
if isfield(algoParams, 'snr_thresh')
    % Threshold supplied by the caller, so that slices corrected in separate
    % jobs all share the same volume-wide value
    snr_thresh = algoParams.snr_thresh;
else
    mask_mag = mag > 0.1*max(mag(:));
    snr_thresh = prctile(mag(mask_mag),25);
end

for kk = vec_slices
    for xx = 1:matrix_size(1)
        for yy = 1:matrix_size(2)
            % Evaluate pixel against SNR threshold
            if mag(xx,yy,kk) > snr_thresh
                ff(xx,yy,kk) = c_ff(xx,yy,kk);
                wf(xx,yy,kk) = c_wf(xx,yy,kk);
            else
                % Binary assignment
                if c_ff(xx,yy,kk) >= algoParams.weight
                    ff(xx,yy,kk) = 1;
                    wf(xx,yy,kk) = 0;
                else % c_ff(xx,yy,kk) < algoParams.weight
                    ff(xx,yy,kk) = 0;
                    wf(xx,yy,kk) = 1;
                end
            end
        end
    end
end

ff = ff.*mask(:,:,:);
wf = wf.*mask(:,:,:);

clear c_ff c_wf mag mask_mag

if algoParams.plot_debug
% Binary mask derived from the initial fat fraction

    figure(102)

    if algoParams.crameri_colormap 
        colormap_plot = crameri('-lajolla');
    else
        colormap_plot = colormap("hot");
    end

    fp(1) = subplot(1,2,1);
    imagesc(ff(:,:,slice_image))
    axis image
    axis off
    colormap(fp(1),colormap_plot);
    title('Fat Mask')

    fp(2) = subplot(1,2,2);
    imagesc(wf(:,:,slice_image))
    axis image
    axis off
    colormap(fp(2),colormap_plot);
    title('Water Mask')

end

%% Estimating the error maps using the water and fat signals (these maps are used as initial guesses to determine phase and amplitude effects through an optimization algorithm)

complex_map1_water = (Water_GC_odd.*conj(Water_GC_even))./(Water_GC_odd.*conj(Water_GC_odd));
complex_map1_fat = (Fat_GC_odd.*conj(Fat_GC_even))./(Fat_GC_odd.*conj(Fat_GC_odd));

complex_map1_combined = (complex_map1_water.^wf .* complex_map1_fat.^ff);

%% Calculation of initial phi and eps maps

phi_map_init = -angle(complex_map1_combined)/2;

eps_map_init = log(abs(complex_map1_combined))/2;

%% Least square solution

options = optimoptions('lsqlin','Algorithm','trust-region-reflective','Display','off');

tik_reg = algoParams.tik_reg;

% Preallocate
b1 = zeros(matrix_size(1:3));
b = zeros([2 matrix_size(1:3)]);
b2 = b1;
correction_map_unwrapped = b1;
correction_map = b1;

% Calculate b matrix
b1(:,:,vec_slices) = 0.5.*(log(Water_GC_even(:,:,vec_slices)) - log(Water_GC_odd(:,:,vec_slices)));
b2(:,:,vec_slices) = 0.5.*(log(Fat_GC_even(:,:,vec_slices)) - log(Fat_GC_odd(:,:,vec_slices)));
b(:,:,:,vec_slices) = cat(1, reshape(wf(:,:,vec_slices).*b1(:,:,vec_slices), [1 size(wf(:,:,vec_slices))]), reshape(ff(:,:,vec_slices).*b2(:,:,vec_slices), [1 size(wf(:,:,vec_slices))]));

if VERBOSE
    tic
end

if algoParams.parallel
    % Cell allocation
    b_cell = cell(nSlices,1);
    mask_cell = b_cell;
    wf_cell = b_cell;
    ff_cell = b_cell;
    phi_map_init_cell = b_cell;
    correction_map_unwrapped_cell = b_cell;
    correction_map_cell = b_cell;
    for kk = vec_slices
        b_cell{kk} = b(:,:,:,kk);
        mask_cell{kk} = mask(:,:,kk);
        wf_cell{kk} = wf(:,:,kk);
        ff_cell{kk} = ff(:,:,kk);
        phi_map_init_cell{kk} = phi_map_init(:,:,kk);
        correction_map_unwrapped_cell{kk} = correction_map_unwrapped(:,:,kk);
        correction_map_cell{kk} = correction_map(:,:,kk);
    end

    parfor kk = vec_slices
        % Initialize/clear solution array
        solution_reg_real = zeros([2 matrix_size(1:2)]); % 2 x nx x ny preallocation
        solution_reg_imag = solution_reg_real;

        for xx = 1:matrix_size(1)
            for yy = 1:matrix_size(2)
                if mask_cell{kk}(xx,yy) == 1
    
                    A = [1i*wf_cell{kk}(xx,yy),wf_cell{kk}(xx,yy);1i*ff_cell{kk}(xx,yy),ff_cell{kk}(xx,yy)];
    
                    A_tik = ctranspose(A)*A + tik_reg*eye(2);
                    b_tik = ctranspose(A) * b_cell{kk}(:,xx,yy);
    
                    if tik_reg == 0
                        solution_reg_real(:,xx,yy) = lsqminnorm(real(A),real(b_cell{kk}(:,xx,yy)));
                        solution_reg_imag(:,xx,yy) = lsqlin(imag(A),imag(b_cell{kk}(:,xx,yy)),[],[],[],[],[-pi;-pi],[pi;pi],[phi_map_init_cell{kk}(xx,yy);0],options);
                    else
                        solution_reg_real(:,xx,yy) = lsqminnorm(real(A_tik),real(b_tik));
                        solution_reg_imag(:,xx,yy) = lsqlin(imag(A_tik),imag(b_tik),[],[],[],[],[-pi;-pi],[pi;pi],[phi_map_init_cell{kk}(xx,yy);0],options);
                    end
                end
            end
        end

        correction_map_unwrapped_cell{kk} = squeeze(1i.*solution_reg_imag(1,:,:) + solution_reg_real(2,:,:));

        % Apply artificial wrap
        solution_reg_imag2 = squeeze(solution_reg_imag(1,:,:));
        solution_reg_imag2(solution_reg_imag2>pi/2) = solution_reg_imag2(solution_reg_imag2>pi/2) - pi;
        solution_reg_imag2(solution_reg_imag2<-pi/2) = solution_reg_imag2(solution_reg_imag2<-pi/2) + pi;
        correction_map_cell{kk} = 1i.*solution_reg_imag2 + squeeze(solution_reg_real(2,:,:));
    end
       
    % Extract from cell array
    for kk = vec_slices
        correction_map_unwrapped(:,:,kk) = correction_map_unwrapped_cell{kk};
        correction_map(:,:,kk) = correction_map_cell{kk};
    end

    clear b_cell mask_cell wf_cell ff_cell phi_map_init_cell correction_map_unwrapped_cell correction_map_cell

else
    % Initialize progress display
    perc = 0;
    maxPerc = nSlices*matrix_size(1);
    reverseStr = '';

    for kk = vec_slices
        % Initialize/clear solution array
        solution_reg_real = zeros([2 matrix_size(1:2)]);
        solution_reg_imag = solution_reg_real;

        for xx = 1:matrix_size(1)
            if VERBOSE
                perc = perc + 1;
                reverseStr = UpdatePercent(perc/maxPerc*100, reverseStr);
            end

            for yy = 1:matrix_size(2)
                if mask(xx,yy,kk) == 1

                    A = [1i*wf(xx,yy,kk),wf(xx,yy,kk);1i*ff(xx,yy,kk),ff(xx,yy,kk)];

                    A_tik = ctranspose(A)*A + tik_reg*eye(2);
                    b_tik = ctranspose(A) * b(:,xx,yy,kk);

                    if tik_reg == 0
                        solution_reg_real(:,xx,yy) = lsqminnorm(real(A),real(b(:,xx,yy,kk)));
                        solution_reg_imag(:,xx,yy) = lsqlin(imag(A),imag(b(:,xx,yy,kk)),[],[],[],[],[-pi;-pi],[pi;pi],[phi_map_init(xx,yy,kk);0],options);
                    else
                        solution_reg_real(:,xx,yy) = lsqminnorm(real(A_tik),real(b_tik));
                        solution_reg_imag(:,xx,yy) = lsqlin(imag(A_tik),imag(b_tik),[],[],[],[],[-pi;-pi],[pi;pi],[phi_map_init(xx,yy,kk);0],options);
                    end
                end
            end
        end

        correction_map_unwrapped(:,:,kk) = squeeze(1i.*solution_reg_imag(1,:,:) + solution_reg_real(2,:,:));

        % Apply artificial wrap
        solution_reg_imag2 = squeeze(solution_reg_imag(1,:,:));
        solution_reg_imag2(solution_reg_imag2>pi/2) = solution_reg_imag2(solution_reg_imag2>pi/2) - pi;
        solution_reg_imag2(solution_reg_imag2<-pi/2) = solution_reg_imag2(solution_reg_imag2<-pi/2) + pi;
        correction_map(:,:,kk) = 1i.*solution_reg_imag2 + squeeze(solution_reg_real(2,:,:));
    end
end

if VERBOSE
    tttime = toc;
    fprintf('Done (%.2f sec)\n', tttime)
end

% Compute residual for all voxels and all field values
% Note: the residual is computed in a vectorized way, for increased speed
phi_map = imag(correction_map);
eps_map = real(correction_map);

%% Comparison of initial guess and final solution for phi and eps maps

% This plot shows water and fat specific phase error and amplitude
% modulation maps and the corresponding maps combined for fat and water
% mixtures. The final row presents the maps that are used to correct for
% bipolar induced effects

if algoParams.plot_debug

    figure(103)
    tiledlayout(2,4,"TileSpacing","compact","Padding","compact")

    if algoParams.crameri_colormap 
        colormap_paramater = crameri('-grayC');
    else
        colormap_paramater = colormap("gray");
    end

    phi_water = -squeeze(angle(complex_map1_water(:,:,:)))/2;

    phi_fat = -squeeze(angle(complex_map1_fat(:,:,:)))/2;

    eps_water = log(abs(complex_map1_water))/2;

    eps_fat = log(abs(complex_map1_fat))/2;

    nexttile
    imagesc(phi_water(:,:,slice_image).*mask(:,:,slice_image))
    axis image
    axis off
    title("Water \phi_{initial}", "Interpreter",'tex')

    nexttile
    imagesc(phi_fat(:,:,slice_image).*mask(:,:,slice_image))
    axis image
    axis off
    title("Fat \phi_{initial}", "Interpreter",'tex')

    nexttile
    imagesc(phi_map_init(:,:,slice_image).*mask(:,:,slice_image))
    axis image
    axis off
    title("Combined \phi_{initial}", "Interpreter",'tex')

    nexttile
    fig_imag = phi_map(:,:,slice_image).*mask(:,:,slice_image);
    fig_imag(mask(:,:,slice_image)==0) = NaN;
    imagesc(fig_imag)
    axis image
    axis off
    title("\phi_{final}", "Interpreter",'tex')
    hc = colorbar;
    clim([-1.5 1.5])
    ylabel(hc,'[rad]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile
    imagesc(eps_water(:,:,slice_image).*mask(:,:,slice_image))
    axis image
    axis off
    clim([-1 1])
    title("Water \epsilon_{initial}", "Interpreter",'tex')

    nexttile
    imagesc(eps_fat(:,:,slice_image).*mask(:,:,slice_image))
    axis image
    axis off
    clim([-1 1])
    title("Fat \epsilon_{initial}", "Interpreter",'tex')

    nexttile
    imagesc(eps_map_init(:,:,slice_image).*mask(:,:,slice_image))
    axis image
    axis off
    clim([-1 1])
    title("Combined \epsilon_{initial}", "Interpreter",'tex')

    nexttile
    imagesc(eps_map(:,:,slice_image).*mask(:,:,slice_image))
    axis image
    axis off
    clim([-1 1])
    title("\epsilon_{final}", "Interpreter",'tex')
    hc = colorbar;
    clim([-1 1])
    ylabel(hc,'[rad]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);
    colormap(colormap_paramater)
end

%% Plot to check results

% These plots show the magnitude and phase results for the fat-water separation

if algoParams.plot_debug

    if algoParams.crameri_colormap 
        colormap_paramater = crameri('-grayC');
    else
        colormap_paramater = colormap("gray");
    end

    FW_plot = figure(104);
    FW_tiles = tiledlayout(3,6,"TileSpacing","compact","Padding","compact");

    nexttile(FW_tiles, 1)
    imagesc(sqrt(Water_GC_odd(:,:,slice_image).*conj(Water_GC_odd(:,:,slice_image))).*mask(:,:,slice_image))
    axis image
    axis off
    title("Water Image (TE_{odd})","Interpreter","tex")
    hc = colorbar;
    % clim([0 3000000])
    ylabel(hc,'[a.u.]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile(FW_tiles, 7)
    imagesc(sqrt(Water_GC_even(:,:,slice_image).*conj(Water_GC_even(:,:,slice_image))).*mask(:,:,slice_image))
    axis image
    axis off
    title("Water Image (TE_{even})","Interpreter","tex")
    hc = colorbar;
    % clim([0 3000000])
    ylabel(hc,'[a.u.]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile(FW_tiles, 2)
    imagesc(sqrt(Fat_GC_odd(:,:,slice_image).*conj(Fat_GC_odd(:,:,slice_image))).*mask(:,:,slice_image))
    axis image
    axis off
    title("Fat Image (TE_{odd})","Interpreter","tex")
    hc = colorbar;
    % clim([0 3000000])
    ylabel(hc,'[a.u.]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile(FW_tiles, 8)
    imagesc(sqrt(Fat_GC_even(:,:,slice_image).*conj(Fat_GC_even(:,:,slice_image))).*mask(:,:,slice_image))
    axis image
    axis off
    title("Fat Image (TE_{even})","Interpreter","tex")
    hc = colorbar;
    % clim([0 3000000])
    ylabel(hc,'[a.u.]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile(FW_tiles, 3)
    imagesc(FieldMap_DualGC(:,:,slice_image));
    axis image
    axis off
    title("Field Map Dual GC")
    hc = colorbar;
    clim(algoParams.range_fm)
    ylabel(hc,'[Hz]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile(FW_tiles, 4)
    imagesc(R2_DualGC(:,:,slice_image).*mask(:,:,slice_image));
    axis image
    axis off
    title("R_2^* Map Dual GC","Interpreter","tex")
    hc = colorbar;
    clim([0 algoParams.range_r2star(2)])
    ylabel(hc,'[s^{-1}]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile(FW_tiles, 5)
    imagesc(phi_map(:,:,slice_image).*mask(:,:,slice_image))
    axis image
    axis off
    title("Phase Modulation (TE_2)","Interpreter","tex")
    hc = colorbar;
    %clim([0 1])
    ylabel(hc,'[rad]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile(FW_tiles, 6)
    % fig_imag = exp(-1.*eps_map(:,:,slice_image));
    imagesc(eps_map(:,:,slice_image).*mask(:,:,slice_image))
    axis image
    axis off
    title("Amplitude Modulation (TE_2)","Interpreter","tex")
    hc = colorbar;
    %clim([0 1.5])
    ylabel(hc,'[a.u.]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);
    colormap(colormap_paramater);

    FW_plot2 = figure(105);
    FW_tiles2 = tiledlayout(3,6,"TileSpacing","compact","Padding","compact");

    nexttile(FW_tiles2, 1)
    imagesc(angle(Water_GC_odd(:,:,slice_image)).*mask(:,:,slice_image))
    axis image
    axis off
    title("Water Phase (TE_{odd})","Interpreter","tex")
    hc = colorbar;
    %clim([0 1500])
    ylabel(hc,'[a.u.]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile(FW_tiles2, 7)
    imagesc(angle(Water_GC_even(:,:,slice_image)).*mask(:,:,slice_image))
    axis image
    axis off
    title("Water Phase (TE_{even})","Interpreter","tex")
    hc = colorbar;
    %clim([0 1500])
    ylabel(hc,'[a.u.]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile(FW_tiles2, 2)
    imagesc(angle(Fat_GC_odd(:,:,slice_image)).*mask(:,:,slice_image))
    axis image
    axis off
    title("Fat Phase (TE_{odd})","Interpreter","tex")
    hc = colorbar;
    %clim([0 1500])
    ylabel(hc,'[a.u.]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile(FW_tiles2, 8)
    imagesc(angle(Fat_GC_even(:,:,slice_image)).*mask(:,:,slice_image))
    axis image
    axis off
    title("Fat Phase (TE_{even})","Interpreter","tex")
    hc = colorbar;
    %clim([0 1500])
    ylabel(hc,'[a.u.]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile(FW_tiles2, 3)
    imagesc(FieldMap_DualGC(:,:,slice_image));
    axis image
    axis off
    title("Field Map Dual GC")
    hc = colorbar;
    clim(algoParams.range_fm)
    ylabel(hc,'[Hz]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile(FW_tiles2, 4)
    imagesc(R2_DualGC(:,:,slice_image).*mask(:,:,slice_image));
    axis image
    axis off
    title("R_2^* Map Dual GC","Interpreter","tex")
    hc = colorbar;
    clim([0 algoParams.range_r2star(2)])
    ylabel(hc,'[s^{-1}]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile(FW_tiles2, 5)
    imagesc(phi_map(:,:,slice_image).*mask(:,:,slice_image))
    axis image
    axis off
    title("Phase Modulation (TE_2)","Interpreter","tex")
    hc = colorbar;
    %clim([0 1])
    ylabel(hc,'[rad]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile(FW_tiles2, 6)
    imagesc(eps_map(:,:,slice_image).*mask(:,:,slice_image))
    axis image
    axis off
    title("Amplitude Modulation (TE_2)","Interpreter","tex")
    hc = colorbar;
    %clim([0 1.5])
    ylabel(hc,'[a.u.]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);
    colormap(colormap_paramater);
end


%% Correction of phi map and unwrapping

% Unwrapping is included in case that corrections maps present phase wraps.
% Phase wraps will not affect the fat-water separation, but unwrapping was
% included in case that user wants to display correction maps without phase
% wraps. In case of wanting to use phase unwrapping, add to the directory
% the code from Fortier and Levesque, DOI 10.1002/mrm.26989, 2017, https://gitlab.com/veronique_fortier/Quality_guided_unwrapping

unwrapping = 0;

if unwrapping

    % Allocation of memory for unwrapped map
    unwrapped_map = zeros(size(eps_map));

    % Unwrapping
    unwrapped_map(:,:,vec_slices) = qualityGuidedUnwrapping(squeeze(imag(correction_map_unwrapped(:,:,vec_slices))), squeeze(mask(:,:,vec_slices)), squeeze(abs(complex_map1_combined(:,:,vec_slices))));

    % Calculation of phi map
    phi_map(:,:,vec_slices) = unwrapped_map(:,:,vec_slices);

else % Keep phi map as least square result
    phi_map(isnan(phi_map))=0;
end

%% Total exponential term for correction

total_correction = exp(1i*(phi_map - 1i*eps_map)); %exp(2*correction_map).^(1/2);

%% Plot to check results

if algoParams.plot_debug
    nexttile(FW_tiles, 11)
    imagesc(angle(total_correction(:,:,slice_image)).*mask(:,:,slice_image))
    axis image
    axis off
    title("Phase Modulation","Interpreter","tex")
    hc = colorbar;
    %clim([0 1])
    ylabel(hc,'[rad]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile(FW_tiles, 12)
    imagesc(log(abs(total_correction(:,:,slice_image))).*mask(:,:,slice_image))
    axis image
    axis off
    title("Amplitude Modulation","Interpreter","tex")
    hc = colorbar;
    %clim([0 1.5])
    ylabel(hc,'[a.u.]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile(FW_tiles2, 11)
    imagesc(angle(total_correction(:,:,slice_image)).*mask(:,:,slice_image))
    axis image
    axis off
    title("Phase Modulation","Interpreter","tex")
    hc = colorbar;
    %clim([0 1])
    ylabel(hc,'[rad]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile(FW_tiles2, 12)
    imagesc(log(abs(total_correction(:,:,slice_image))).*mask(:,:,slice_image))
    axis image
    axis off
    title("Amplitude Modulation","Interpreter","tex")
    hc = colorbar;
    %clim([0 1.5])
    ylabel(hc,'[a.u.]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);
end

%% Deconvolving the effect of the errors introduce by the readout

s0 = permute(input_signal_bipolar_RO(:,:,vec_slices,:,:), [length(matrix_size) 1:length(matrix_size)-1]);
s0 = reshape(s0,[numte numvox]);

counter_t = reshape(1:numte,[numte 1]);
counter_t = repmat(counter_t,[1 numvox]);

total_correction_reshaped = reshape(total_correction(:,:,vec_slices),[1 numvox]);

E = repmat(total_correction_reshaped,[numte 1]) .^ ((-1).^(counter_t));

corrected_signal = s0 ./ (E);

corrected_signal = permute(corrected_signal,[2 1]);

corrected_bipolar_signal = zeros([matrix_size(1:3) 1 numte], 'like', input_signal_bipolar_RO);
corrected_bipolar_signal(:,:,vec_slices,:,:) = reshape(corrected_signal,matrix_size(1),matrix_size(2),nSlices,1,numte);

outParams.bipolar_error_map_theta = phi_map - 1i*eps_map;

%% Fat water separation of the corrected signal

imDataParams.TE = TEs_bipolar;

if VERBOSE
    tic
    fprintf('\nFat-water separation in synthetic unipolar dataset slice ');
end


if algoParams.parallel
    % Allocate updated image data
    for kk = 1:matrix_size(3)
        imDataParams_cell{kk}.images = corrected_bipolar_signal(:,:,kk,:,:);
    end

    parfor kk = vec_slices        
        outParams_bipolar = fw_i2cm1i_3pluspoint_hernando_graphcut(imDataParams_cell{kk}, algoParams, false);
    
        Water_bipolar{kk} = outParams_bipolar.species(1).amps;
        Fat_bipolar{kk} = outParams_bipolar.species(2).amps;
        FieldMap_bipolar{kk} = outParams_bipolar.fieldmap;
        R2_bipolar{kk} = outParams_bipolar.r2starmap;
    end

    % Extract from cell arrays
    Water_bipolar_cell = Water_bipolar; Water_bipolar = zeros(matrix_size(1:3));
    Fat_bipolar_cell = Fat_bipolar; Fat_bipolar = Water_bipolar;
    FieldMap_bipolar_cell = FieldMap_bipolar; FieldMap_bipolar = Water_bipolar;
    R2_bipolar_cell = R2_bipolar; R2_bipolar = Water_bipolar;
    for kk = vec_slices
        Water_bipolar(:,:,kk) = Water_bipolar_cell{kk};
        Fat_bipolar(:,:,kk) = Fat_bipolar_cell{kk};
        FieldMap_bipolar(:,:,kk) = FieldMap_bipolar_cell{kk};
        R2_bipolar(:,:,kk) = R2_bipolar_cell{kk};
    end
    clear Water_bipolar_cell Fat_bipolar_cell FieldMap_bipolar_cell R2_bipolar_cell

else
    for kk = vec_slices
        if VERBOSE
            fprintf([num2str(kk) ', ']);
        end
        imDataParams.images = corrected_bipolar_signal(:,:,kk,:,:); % Originaly size[nx,ny,nz,ncoils,nechoes] Note: This specific order is required for GC algorithm

        outParams_bipolar = fw_i2cm1i_3pluspoint_hernando_graphcut(imDataParams, algoParams, false);

        Water_bipolar(:,:,kk) = outParams_bipolar.species(1).amps;
        Fat_bipolar(:,:,kk) = outParams_bipolar.species(2).amps;
        FieldMap_bipolar(:,:,kk) = outParams_bipolar.fieldmap;
        R2_bipolar(:,:,kk) = outParams_bipolar.r2starmap;
    end
end

if VERBOSE
    tttime = toc;
    fprintf('. Done (%.2f sec)\n', tttime)
end

%% Figure to check the results

if algoParams.plot_debug

    nexttile(FW_tiles, 13)
    imagesc(abs(Water_bipolar(:,:,slice_image)).*mask(:,:,slice_image));
    axis image
    axis off
    title("Water Image")
    hc = colorbar;
    % clim([0 3000000])
    ylabel(hc,'[a.u.]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile(FW_tiles, 14)
    imagesc(abs(Fat_bipolar(:,:,slice_image)).*mask(:,:,slice_image));
    axis image
    axis off
    title("Fat Image")
    hc = colorbar;
    % clim([0 3000000])
    ylabel(hc,'[a.u.]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile(FW_tiles, 15)
    imagesc(FieldMap_bipolar(:,:,slice_image));
    axis image
    axis off
    title("Field Map")
    hc = colorbar;
    clim(algoParams.range_fm)
    ylabel(hc,'[Hz]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile(FW_tiles, 16)
    imagesc(R2_bipolar(:,:,slice_image).*mask(:,:,slice_image));
    axis image
    axis off
    title("R_2^* Map","Interpreter","tex")
    hc = colorbar;
    clim([0 algoParams.range_r2star(2)])
    ylabel(hc,'[s^{-1}]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile(FW_tiles2, 13)
    imagesc(angle(Water_bipolar(:,:,slice_image)).*mask(:,:,slice_image));
    axis image
    axis off
    title("Water Phase")
    hc = colorbar;
    %clim([0 1500])
    ylabel(hc,'[a.u.]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile(FW_tiles2, 14)
    imagesc(angle(Fat_bipolar(:,:,slice_image)).*mask(:,:,slice_image));
    axis image
    axis off
    title("Fat Phase")
    hc = colorbar;
    %clim([0 1500])
    ylabel(hc,'[a.u.]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile(FW_tiles2, 15)
    imagesc(FieldMap_bipolar(:,:,slice_image));
    axis image
    axis off
    title("Field Map")
    hc = colorbar;
    % clim([-30 30])
    ylabel(hc,'[Hz]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);

    nexttile(FW_tiles2, 16)
    imagesc(R2_bipolar(:,:,slice_image).*mask(:,:,slice_image));
    axis image
    axis off
    title("R_2^* Map","Interpreter","tex")
    hc = colorbar;
    %clim([0 100])
    ylabel(hc,'[s^{-1}]','Units', 'normalized', 'Position', [2.5, 0.5],'Rotation',90);
end

%% Creating the output for the function

% Water signal for fat-water separation using all echoes
outParams.species(1).amps = Water_bipolar;

% Fat signal for fat-water separation using all echoes
outParams.species(2).amps = Fat_bipolar;

% R2 star map using all echoes
outParams.r2starmap = R2_bipolar;

% Field map using all echoes
outParams.fieldmap = FieldMap_bipolar;

% Phi map (related to phase modulation due to bipolar readout)
outParams.phi_map = phi_map;

% Epsilon map (related to amplitude modulation due to bipolar readout)
outParams.eps_map = eps_map;

% Corrected signal: bipolar signal transformed into the unipolar equivalent
outParams.corrected_bipolar_signal = corrected_bipolar_signal;

% Correction to remove bipolar induced effects
outParams.total_correction = total_correction;

% Results from dualGC odd and even echoes
outParams.Water_GC_odd = Water_GC_odd;
outParams.Fat_GC_odd = Fat_GC_odd;
outParams.Water_GC_even = Water_GC_even;
outParams.Fat_GC_even = Fat_GC_even;
outParams.FieldMap_DualGC = FieldMap_DualGC;
outParams.R2_DualGC = R2_DualGC;
% Preserve prior behavior: keep the residual of the last processed slice
% (outParams_GC no longer survives the parfor loop above).
% outParams.residual = residuals{vec_slices(end)};
end

% Update display percentage
function revstr = UpdatePercent(perc, revstr)
    msg = sprintf('%.2f percent. ', perc);
    fprintf([revstr, msg]);
    revstr = repmat(sprintf('\b'), 1, length(msg));
end