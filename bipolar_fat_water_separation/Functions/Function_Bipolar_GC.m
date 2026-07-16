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

%% Signal for the acquistion with bipolar gradients and without any correction

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
numvox = prod(matrix_size(1:3));

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
residual = zeros([algoParams.NUM_FMS, matrix_size(1:3)]);
if algoParams.parallel
    Water_GC_odd = squeeze(mat2cell(Water_GC_odd, matrix_size(1), matrix_size(2), ones([matrix_size(3),1])));
    residual = squeeze(mat2cell(residual, algoParams.NUM_FMS, matrix_size(1), matrix_size(2), ones([matrix_size(3),1])));
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
residual = Water_GC_odd;

%% Parallel pool sized to the SLURM allocation
if algoParams.parallel
    % Falls back to local core count if not sized to SLURM core count
    if isfield(algoParams, 'nWorkers')
        nWorkers = algoParams.nWorkers;
    else
        nWorkers = str2double(getenv('SLURM_CPUS_PER_TASK'));
        if isnan(nWorkers) || nWorkers < 1
            nWorkers = maxNumCompThreads;
        end
        nWorkers = max(1, min(nWorkers, numel(vec_slices)));
    end
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
    parpool(c, nWorkers);

    fprintf("Number of workers is set to: %i", nWorkers);
end

%% Graph-cut fat-water separation for odd and even echoes
if VERBOSE
    fprintf('\nFat-water separation for odd and even echo datasets slice ');
end

if algoParams.parallel
    % assign only the relevant indexed parameters to each tmp worker
    imDataParams_cell = Water_GC_odd;
    for kk = 1:matrix_size(3)
        imDataParams_cell{kk} = imDataParams;
        imDataParams_cell{kk}.images = imDataParams_cell{kk}.images(:,:,kk,:,:);
        imDataParams_cell{kk}.mask_fwseparation = imDataParams_cell{kk}.mask_fwseparation(:,:,kk);
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
        residual{kk} = outParams_GC.residual;
    end

    % Extract from cell arrays
    Water_GC_odd = cell2mat(reshape(Water_GC_odd,[1 1 matrix_size(3)]));
    Fat_GC_odd = cell2mat(reshape(Fat_GC_odd,[1 1 matrix_size(3)]));
    Water_GC_even = cell2mat(reshape(Water_GC_even,[1 1 matrix_size(3)]));
    Fat_GC_even = cell2mat(reshape(Fat_GC_even,[1 1 matrix_size(3)]));
    FieldMap_DualGC = cell2mat(reshape(FieldMap_DualGC,[1 1 matrix_size(3)]));
    R2_DualGC = cell2mat(reshape(R2_DualGC,[1 1 matrix_size(3)]));
    residual = cell2mat(reshape(residual,[1 1 1 matrix_size(3)]));

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
        residual(:,:,kk) = outParams_GC.residual;
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

%% Binary fat and water mask (after thresholding using a specific values included in structure params)

% ff = zeros(matrix_size(1:3));
% 
% ff (c_ff>=algoParams.weight) = 1;
% ff (c_ff<algoParams.weight) = 0;
% 
% ff = ff.*mask(:,:,:);
% wf = (1 - ff).*mask(:,:,:);

% option 1
ff = c_ff.*mask(:,:,:);
wf = c_wf.*mask(:,:,:); % c_wf == 1 - c_ff, already clamped to [0,1]

% option 2
% soft threshold: sharpness k controls how gradual the transition is
% k = 10; % larger k means closer to the old binary behavior
% ff = 1 ./ (1 + exp(-k*(c_ff - algoParams.weight)));
% ff = ff.*mask(:,:,:);
% wf = (1 - ff).*mask(:,:,:);

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

if algoParams.parallel
    % Cell allocation
    b = squeeze(mat2cell(b1, matrix_size(1), matrix_size(2), ones([matrix_size(3),1])));
    mask = squeeze(mat2cell(mask, matrix_size(1), matrix_size(2), ones([matrix_size(3),1])));
    wf = squeeze(mat2cell(wf, matrix_size(1), matrix_size(2), ones([matrix_size(3),1])));
    ff = squeeze(mat2cell(ff, matrix_size(1), matrix_size(2), ones([matrix_size(3),1])));
    phi_map_init = squeeze(mat2cell(phi_map_init, matrix_size(1), matrix_size(2), ones([matrix_size(3),1])));
    correction_map_unwrapped = squeeze(mat2cell(zeros(matrix_size(1:3)), matrix_size(1), matrix_size(2), ones([matrix_size(3),1])));
    correction_map = squeeze(mat2cell(zeros(matrix_size(1:3)), matrix_size(1), matrix_size(2), ones([matrix_size(3),1])));

    parfor kk = vec_slices
        % Initialize/clear solution array
        solution_reg_real = zeros([2 matrix_size(1:2)]);
        solution_reg_imag = solution_reg_real;

        for xx = 1:matrix_size(1)
            for yy = 1:matrix_size(2)
                if mask{kk} == 1
    
                    A = [1i*wf{kk},wf{kk};1i*ff{kk},ff{kk}];
    
                    A_tik = ctranspose(A)*A + tik_reg*eye(2);
                    b_tik = ctranspose(A) * b{kk};
    
                    if tik_reg == 0
                        solution_reg_real = lsqminnorm(real(A),real(b{kk}));
                        solution_reg_imag = lsqlin(imag(A),imag(b{kk}),[],[],[],[],[-pi;-pi],[pi;pi],[phi_map_init{kk};0],options);
                    else
                        solution_reg_real = lsqminnorm(real(A_tik),real(b_tik));
                        solution_reg_imag = lsqlin(imag(A_tik),imag(b_tik),[],[],[],[],[-pi;-pi],[pi;pi],[phi_map_init{kk};0],options);
                    end
                end
            end
        end

        correction_map_unwrapped{kk} = squeeze(1i.*solution_reg_imag(1,:,:) + solution_reg_real(2,:,:));

        % Apply artificial wrap
        solution_reg_imag2 = squeeze(solution_reg_imag(1,:,:));
        solution_reg_imag2(solution_reg_imag2>pi/2) = solution_reg_imag2(solution_reg_imag2>pi/2) - pi;
        solution_reg_imag2(solution_reg_imag2<-pi/2) = solution_reg_imag2(solution_reg_imag2<-pi/2) + pi;
        correction_map{kk} = 1i.*solution_reg_imag2 + squeeze(solution_reg_real(2,:,:));
    end

    correction_map_unwrapped = cell2mat(reshape(correction_map_unwrapped,[1 1 matrix_size(3)]));
    correction_map = cell2mat(reshape(correction_map,[1 1 matrix_size(3)]));

else
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
                        solution_reg_real(:,xx,yy,kk) = lsqminnorm(real(A),real(b(:,xx,yy,kk)));
                        solution_reg_imag(:,xx,yy,kk) = lsqlin(imag(A),imag(b(:,xx,yy,kk)),[],[],[],[],[-pi;-pi],[pi;pi],[phi_map_init(xx,yy,kk);0],options);
                    else
                        solution_reg_real(:,xx,yy,kk) = lsqminnorm(real(A_tik),real(b_tik));
                        solution_reg_imag(:,xx,yy,kk) = lsqlin(imag(A_tik),imag(b_tik),[],[],[],[],[-pi;-pi],[pi;pi],[phi_map_init(xx,yy,kk);0],options);
                    end
                end
            end
        end

        correction_map_unwrapped(:,:,kk) = squeeze(1i.*solution_reg_imag(1,:,:,kk) + solution_reg_real(2,:,:,kk));

        % Apply artificial wrap
        solution_reg_imag2 = squeeze(solution_reg_imag(1,:,:,kk));
        solution_reg_imag2(solution_reg_imag2>pi/2) = solution_reg_imag2(solution_reg_imag2>pi/2) - pi;
        solution_reg_imag2(solution_reg_imag2<-pi/2) = solution_reg_imag2(solution_reg_imag2<-pi/2) + pi;
        correction_map(:,:,kk) = 1i.*solution_reg_imag2 + squeeze(solution_reg_real(2,:,:,kk));
    end
end

if VERBOSE
    tttime = toc;
    fprintf('Done (%.2f sec)', tttime)
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

s0 = permute(input_signal_bipolar_RO, [length(matrix_size) 1:length(matrix_size)-1]);
s0 = reshape(s0,[numte numvox]);

counter_t = reshape(1:numte,[numte 1]);
counter_t = repmat(counter_t,[1 numvox]);

total_correction_reshaped = reshape(total_correction,[1 numvox]);

E = repmat(total_correction_reshaped,[numte 1]) .^ ((-1).^(counter_t));

corrected_signal = s0 ./ (E);

corrected_signal = permute(corrected_signal,[2 1]);

corrected_bipolar_signal = reshape(corrected_signal,matrix_size(1),matrix_size(2),matrix_size(3),1,numte);

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
    Water_bipolar = cell2mat(reshape(Water_bipolar,[1 1 matrix_size(3)]));
    Fat_bipolar = cell2mat(reshape(Fat_bipolar,[1 1 matrix_size(3)]));
    FieldMap_bipolar = cell2mat(reshape(FieldMap_bipolar,[1 1 matrix_size(3)]));
    R2_bipolar = cell2mat(reshape(R2_bipolar,[1 1 matrix_size(3)]));

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
    fprintf('. Done (%.2f sec)', tttime)
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
outParams.residual = residuals{vec_slices(end)};
end

% Update display percentage
function revstr = UpdatePercent(perc, revstr)
    msg = sprintf('%.2f percent. ', perc);
    fprintf([revstr, msg]);
    revstr = repmat(sprintf('\b'), 1, length(msg));
end