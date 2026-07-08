% Function name: Function_i2cm1i_3pluspoint_hernando_Bipolar_GC
%
% Description: Fat-water separation using regularized fieldmap formulation and graph cut solution.
% Modified to enable fat-water separation with bipolar readout gradient pulses
%
% Hernando D, Kellman P, Haldar JP, Liang ZP. Robust water/fat separation in the presence of large
% field inhomogeneities using a graph cut algorithm. Magn Reson Med. 2010 Jan;63(1):79-90.
%
% Some properties:
%   - Image-space
%   - 2 species (water-fat)
%   - Complex-fitting
%   - Multi-peak fat (pre-calibrated)
%   - Single-R2*
%   - Independent water/fat phase
%   - Requires 3+ echoes at arbitrary echo times (some choices are much better than others! see NSA...)
%
% Input: structures imDataParams and algoParams
%   - imDataParams.images: acquired images, array of size[nx,ny,1,ncoils,nTE]
%   - imDataParams.TE: echo times (in seconds)
%   - imDataParams.FieldStrength: (in Tesla)
%
%   - algoParams.species(ii).name = name of species ii (string)
%   - algoParams.species(ii).frequency = frequency shift in ppm of each peak within species ii
%   - algoParams.species(ii).relAmps = relative amplitude (sum normalized to 1) of each peak within species ii
%   Example
%      - algoParams.species(1).name = 'water' % Water
%      - algoParams.species(1).frequency = [0]
%      - algoParams.species(1).relAmps = [1]
%      - algoParams.species(2).name = 'fat' % Fat
%      - algoParams.species(2).frequency = [3.80, 3.40, 2.60, 1.94, 0.39, -0.60]
%      - algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048]
%
%   - algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
%   - algoParams.range_r2star = [0 0]; % Range of R2* values
%   - algoParams.NUM_R2STARS = 1; % Numbre of R2* values for quantization
%   - algoParams.range_fm = [-400 400]; % Range of field map values
%   - algoParams.NUM_FMS = 301; % Number of field map values to discretize
%   - algoParams.NUM_ITERS = 40; % Number of graph cut iterations
%   - algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
%   - algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
%   - algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
%   - algoParams.lambda = 0.05; % Regularization parameter
%   - algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
%   - algoParams.TRY_PERIODIC_RESIDUAL = 0; % Take advantage of periodic residual if uniform TEs (will change range_fm)
%   - algoParams.residual: in case we pre-computed the fit residual (mostly for testing)
% Extra parameters added for bipolar readout gradient pulses 

% % imDataParams.mask_fwseparation; % Enables execution of fat-water separation only on the voxels within the
% % % binary mask
% % 
% % algoParams.slice_image; % Value of slice to display the images to check intermediate results (mostly for debugging purposes)
% % algoParams.plot_debug; % Binary flag (1->debugging 0->no debugging) to display images that show intermediate results (useful for debugging code)
% % algoParams.tik_reg; % Binary flag to include Tikhonov regularization in the inverse problem to claculate phi and eps (1->yes regularization 0->no regularization). Regularization was not included for this paper, but it can enable future research about the the potential need and impact of regularization in these parameters 
% % algoParams.weight; % Weight for the inverse problem to determine correction factors phi and eps. In this paper we tested 0.5
% % algoParams.fm_init;% Binary flag (1->initial guess 0->no initial guess) to enable the use of an initial guess for the field inhomogeneity term 

% Output: structure outParams
%   - outParams.species(ii).name: name of the species (taken from algoParams)
%   - outParams.species(ii).amps: estimated water/fat images, size [nx,ny,ncoils]
%   - outParams.r2starmap: R2* map (in s^{-1}, size [nx,ny])
%   - outParams.fieldmap: field map (in Hz, size [nx,ny])
%
%
% Author: Diego Hernando
% Date created: August 5, 2011
% Date last modified: November 10, 2011
% Modified by Jorge Campos
% Date: May 05, 2025

function outParams = Function_i2cm1i_3pluspoint_hernando_Bipolar_GC(imDataParams, algoParams, VERBOSE)

if nargin < 3
    VERBOSE = 0;
end

% Check validity of params, and set default algorithm parameters if not provided
[validParams,algoParams] = checkParamsAndSetDefaults_graphcut_Bipolar_GC(imDataParams,algoParams);
if validParams==0
    disp('Exiting -- data not processed');
    outParams = [];
    return
end

% Get data dimensions
[sx,sy,sz,C,N] = size(imDataParams.images);

% If more than one slice, pick central slice
if sz > 1
    disp('Multi-slice data: processing slice of interest');
    imDataParams.images = imDataParams.images(:,:,imDataParams.sliceofint,:,:);
end

% If more than one channel, coil combine
if C > 1
    disp('Multi-coil data: coil-combining');
    imDataParams.images = coilCombine(imDataParams.images);
end

% If precession is clockwise (positive fat frequency) simply conjugate data
if imDataParams.PrecessionIsClockwise <= 0
    imDataParams.images = conj(imDataParams.images);
    imDataParams.PrecessionIsClockwise = 1;
end

% Check spatial subsampling option (speedup ~ quadratic SUBSAMPLE parameter)
SUBSAMPLE = algoParams.SUBSAMPLE;
if SUBSAMPLE > 1
    images0 = imDataParams.images;
    START = round(SUBSAMPLE/2);
    [sx,sy] = size(images0(:,:,1,1,1));
    allX = 1:sx;
    allY = 1:sy;
    subX = START:SUBSAMPLE:sx;
    subY = START:SUBSAMPLE:sy;
    imDataParams.images = images0(subX,subY,:,:,:);
end

% Regularization parameter
lambda = algoParams.lambda;

% Spatially-varying regularization.  The LMAP_POWER applies to the
% sqrt of the curvature of the residual, and LMAP_POWER=2 yields
% approximately uniform resolution.
LMAP_POWER = algoParams.LMAP_POWER;

% LMAP_EXTRA: Extra flexibility for including prior knowledge into
% regularization. For instance, it can be used to add more smoothing
% to noise regions (by adding, eg a constant LMAP_EXTRA), or even to
% add spatially-varying smoothing as a function of distance to
% isocenter...
LMAP_EXTRA = algoParams.LMAP_EXTRA;

% Finish off with some optimization transfer -- to remove discretization
DO_OT = algoParams.DO_OT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reshape data in odd and even echoes 

% Total number of echoes

numte = size(imDataParams.images,5);

% Dividing data in odd and even echoes

imDataParams.odd_echoes = imDataParams.images(:,:,:,:,1:2:numte);

imDataParams.even_echoes = imDataParams.images(:,:,:,:,2:2:numte);

imDataParams.TEs_odd = imDataParams.TE(1:2:numte);

imDataParams.TEs_even = imDataParams.TE(2:2:numte);

% Combining the data in a new order

imDataParams.images = cat(5,imDataParams.odd_echoes,imDataParams.even_echoes);

imDataParams.TE = [imDataParams.TEs_odd imDataParams.TEs_even];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Let's get the residual. If it's not already in the params, compute it
if isfield(algoParams, 'residual')
    % Grab the residual from the params structure
    residual = algoParams.residual;
else
    % Check for uniform TE spacings
    dTE = diff(imDataParams.TE);

    if algoParams.TRY_PERIODIC_RESIDUAL && sum(abs(dTE - dTE(1)))<1e-6 % If we have uniform TE spacing
        UNIFORM_TEs = 1;
    else
        UNIFORM_TEs = 0;
    end


    % if VERBOSE
    %     UNIFORM_TEs
    % end

    % Compute the residual
    if UNIFORM_TEs == 1 % TEST DH* 090801
        
    else
        % If not uniformly spaced TEs, get the residual for the whole range
        residual = computeResidual_Bipolar_GC(imDataParams, algoParams, VERBOSE);
    end

end

%save tempres.mat residual

% Setup the estimation, get the lambdamap,...
fms = linspace(algoParams.range_fm(1),algoParams.range_fm(2),algoParams.NUM_FMS);
dfm = fms(2)-fms(1);
lmap = getQuadraticApprox( residual, dfm );
lmap = (sqrt(lmap)).^LMAP_POWER;
lmap = lmap + mean(lmap(:))*LMAP_EXTRA;

% Initialize the field map indices
if algoParams.fm_init == 1
    % Grad the indices from the from the params structure if they were
    % already calculated
    cur_ind = algoParams.index_fm_init;
else
    cur_ind = ceil(length(fms)/2)*ones(size(imDataParams.images(:,:,1,1,1)));
end

% This is the core of the algorithm
if VERBOSE
    tic
    fprintf('\nEstimating field map through iterative graph-cuts...')
end
%fm = graphCutIterations_OLD(imDataParams,algoParams,residual,lmap,cur_ind);
[fm,index_fm] = graphCutIterations_Bipolar_GC(imDataParams,algoParams,residual,lmap,cur_ind);
outParams.index_fm = index_fm;

if VERBOSE
    tttime = toc;
    fprintf('Done (%.2f sec)', tttime)
end

if algoParams.plot_debug
    figure(100)
    subplot(1,2,2)
    imagesc(fm(:,:))
    axis image
    axis off
    clim(algoParams.range_fm)
    title('Final B_0 map','Interpreter','tex');
    colormap('gray')
    colorbar
end

% If we have subsampled (for speed), let's interpolate the field map
if SUBSAMPLE > 1
    fmlowres = fm;
    [SUBX,SUBY] = meshgrid(subY(:),subX(:));
    [ALLX,ALLY] = meshgrid(allY(:),allX(:));
    fm = interp2(SUBX,SUBY,fmlowres,ALLX,ALLY,'*spline');
    lmap = interp2(SUBX,SUBY,lmap,ALLX,ALLY,'*spline');
    fm(isnan(fm)) = 0;
    imDataParams.images = images0;
end
algoParams.lmap = lmap;

% Preallocate data
r2starmap = zeros(size(imDataParams.images,1), size(imDataParams.images,2), size(imDataParams.images,6));
w_odd = zeros(size(imDataParams.images,1), size(imDataParams.images,2), 1, size(imDataParams.images,6));
f_odd = zeros(size(imDataParams.images,1), size(imDataParams.images,2), 1, size(imDataParams.images,6));
w_even = w_odd;
f_even = f_odd;

if VERBOSE
    tic
    fprintf('\nEstimating R2star map...')
end

% Now take the field map fm and get the rest of the estimates
for ka = 1:size(imDataParams.images,6)

    curParams = imDataParams;
    curParams.images = imDataParams.images(:,:,:,:,:,ka);


    if algoParams.range_r2star(2) > 0
        % DH* 100422 use fine R2* discretization at this point
        %algoParams.NUM_R2STARS = round(algoParams.range_r2star(3)/2)+1;
        r2starmap(:,:,ka) = estimateR2starGivenFieldmap_Bipolar_GC(curParams, algoParams, fm);
    else
        r2starmap(:,:,ka) = zeros(size(fm));
    end

    if DO_OT ~= 1
        % If no Optimization Transfer, just get the water/fat images
        [amps,error] = decomposeGivenFieldMapAndDampings_Bipolar_GC(curParams,algoParams, fm,r2starmap(:,:,ka),r2starmap(:,:,ka));

        waterimage_odd = squeeze(amps(:,:,1,:));
        fatimage_odd = squeeze(amps(:,:,2,:));

        waterimage_even = squeeze(amps(:,:,3,:));
        fatimage_even = squeeze(amps(:,:,4,:));

        w_odd(:,:,:,ka) = waterimage_odd;
        f_odd(:,:,:,ka) = fatimage_odd;

        w_even(:,:,:,ka) = waterimage_even;
        f_even(:,:,:,ka) = fatimage_even;

    end
end
if VERBOSE
    tttime = toc;
    fprintf('Done (%.2f sec)', tttime)
end

% Put results in outParams structure
try
    outParams.species(1).name = algoParams.species(1).name;
    outParams.species(2).name = algoParams.species(2).name;
catch
    outParams.species(1).name = 'water';
    outParams.species(2).name = 'fat';
end

outParams.species(1).amps = w_odd;
outParams.species(2).amps = f_odd;

outParams.species(3).amps = w_even;
outParams.species(4).amps = f_even;

outParams.r2starmap = r2starmap;
outParams.fieldmap = fm;