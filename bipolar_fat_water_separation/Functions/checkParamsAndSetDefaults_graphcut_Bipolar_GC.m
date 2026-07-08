% Function: checkParamsAndSetDefaults_graphcut_Bipolar_GC
%
% Description: Check validity of input parameters and set defaults for unspecified parameters
%
% Input:
%   - imDataParams: TEs, images and field strength
%   - algoParams: algorithm parameters
%
% Output:
%   - validParams: binary variable (0 if parameters are not valid for this algorithm)
%   - algoParams2: "completed" algorithm parameter structure (after inserting defaults for unspecified parameters)
% 
%
% Author: Diego Hernando
% Date created: August 19, 2011
% Date last modified: November 10, 2011
% Modified by Jorge Campos
% Date: May 05, 2025

function [validParams,algoParams2] = checkParamsAndSetDefaults_graphcut_Bipolar_GC(imDataParams,algoParams, vec_slices)

imDataParams2 = imDataParams;
algoParams2 = algoParams;
validParams = 1;

% Start by checking validity of provided data and recon parameters
if size(imDataParams,3) > 1
    disp('ERROR: 2D recon -- please format input data as array of size SX x SY x 1 X nCoils X nTE')
    validParams = 0;
end

if length(algoParams.species) > 2
    disp('ERROR: Water=fat recon -- use a multi-species function to separate more than 2 chemical species')
    validParams = 0;
end

if length(imDataParams.TE) < 3
    disp('ERROR: 3+ point recon -- please use a different recon for acquisitions with fewer than 3 TEs')
    validParams = 0;
end

%%   - algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
if isfield(algoParams, 'size_clique')
    algoParams2.size_clique = algoParams.size_clique;
else 
    algoParams2.size_clique = 1;
end

%%   - algoParams.range_r2star = [0 0]; % Range of R2* values
if isfield(algoParams, 'range_r2star')
    algoParams2.range_r2star = algoParams.range_r2star;
else 
    algoParams2.range_r2star = [0 0];
end

%%   - algoParams.NUM_R2STARS = 1; % Numbre of R2* values for quantization
if isfield(algoParams, 'NUM_R2STARS')
    algoParams2.NUM_R2STARS = algoParams.NUM_R2STARS;
else 
    algoParams2.NUM_R2STARS = 1;
end

%%   - algoParams.range_fm = [-400 400]; % Range of field map values
if isfield(algoParams, 'range_fm')
    algoParams2.range_fm = algoParams.range_fm;
else 
    algoParams2.range_fm = [-400 400]; 
end

%%   - algoParams.NUM_FMS = 301; % Number of field map values to discretize
if isfield(algoParams, 'NUM_FMS')
    algoParams2.NUM_FMS = algoParams.NUM_FMS;
else 
    algoParams2.NUM_FMS = 301;
end

%%   - algoParams.NUM_ITERS = 40; % Number of graph cut iterations
if isfield(algoParams, 'NUM_ITERS')
    algoParams2.NUM_ITERS = algoParams.NUM_ITERS;
else 
    algoParams2.NUM_ITERS = 40;
end

%%   - algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
if isfield(algoParams, 'SUBSAMPLE')
    algoParams2.SUBSAMPLE = algoParams.SUBSAMPLE;
else 
    algoParams2.SUBSAMPLE = 1;
end

%%   - algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
if isfield(algoParams, 'DO_OT')
    algoParams2.DO_OT = algoParams.DO_OT;
else 
    algoParams2.DO_OT = 0;
end

%%   - algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
if isfield(algoParams, 'LMAP_POWER')
    algoParams2.LMAP_POWER = algoParams.LMAP_POWER;
else 
    algoParams2.LMAP_POWER = 2; 
end

%%   - algoParams.lambda = 0.05; % Regularization parameter
if isfield(algoParams, 'lambda')
    algoParams2.lambda = algoParams.lambda;
else 
    algoParams2.lambda = 0.05;
end

%%   - algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
if isfield(algoParams, 'LMAP_EXTRA')
    algoParams2.LMAP_EXTRA = algoParams.LMAP_EXTRA;
else 
    algoParams2.LMAP_EXTRA = zeros(size(imDataParams.images(:,:,1,1,1)));
end

%%   - algoParams.TRY_PERIODIC_RESIDUAL = 1; % Take advantage of periodic residual if uniform TEs (will change range_fm)  
if isfield(algoParams, 'TRY_PERIODIC_RESIDUAL')
    algoParams2.TRY_PERIODIC_RESIDUAL = algoParams.TRY_PERIODIC_RESIDUAL;
else 
    algoParams2.TRY_PERIODIC_RESIDUAL = 1;
end

%%   - imDataParams.PrecessionIsClockwise (1 = fat has positive frequency; -1 = fat has negative frequency)
if isfield(algoParams, 'PrecessionIsClockwise')
    imDataParams2.PrecessionIsClockwise = imDataParams.PrecessionIsClockwise;
    % If precession is clockwise (positive fat frequency) simply conjugate data
    if imDataParams2.PrecessionIsClockwise <= 0
        imDataParams2.images = conj(imDataParams.images);
        imDataParams2.PrecessionIsClockwise = -1;
    end
else
    imDataParams2.PrecessionIsClockwise = -1;
end

%%   - imDataParams.gyro = 42.5774780505984
if isfield(algoParams, 'gyro')
    algoParams2.gyro = algoParams.gyro;
else
    algoParams2.gyro = 42.5774780505984;
end

%% Extra parameters to enable fat-water separation with bipolar readouts

%% - imDataParams.mask_fwseparation (binary mask to perform fat-water separation in voxels within the mask)
if isfield(algoParams, 'mask_fwseparation')
    imDataParams2.mask_fwseparation = imDataParams.mask_fwseparation;
else
    imDataParams2.mask_fwseparation = 1;
end

%% - algoParams.slice_image (slice to display images to preview results of fat-water separation)
if isfield(algoParams, 'slice_image')
    algoParams2.slice_image = algoParams.slice_image;
elseif nargin > 2 && isscalar(vec_slices)
    algoParams2.slice_image = vec_slices;
else
    algoParams2.slice_image = round(size(imDataParams.images,3)/2);
end

%% - algoParams.plot_debug (Binary flag (1->debugging 0->no debugging) to display images that show intermediate results (useful for debugging code))
if isfield(algoParams, 'plot_debug')
    algoParams2.plot_debug = algoParams.plot_debug;
else
    algoParams2.plot_debug = false;
end

%% - algoParams.crameri_colormap (Binary flag (1->present 0->not present))
if isfield(algoParams, 'crameri_colormap')
    algoParams2.crameri_colormap = algoParams.crameri_colormap;
else
    algoParams2.crameri_colormap = true;
end

%% - algoParams.tik_reg Binary flag to include Tikhonov regularization in the inverse problem to claculate phi and eps (1->yes regularization 0->no regularization). Regularization was not included for this paper, but it can enable future research about the the potential need and impact of regularization in these parameters 
if isfield(algoParams, 'tik_reg')
    algoParams2.tik_reg = algoParams.tik_reg;
else
    algoParams2.tik_reg = 0;
end

%% - algoParams.weight (Weight for the inverse problem to determine correction factors phi and eps. In this paper we tested 0.5)
if isfield(algoParams, 'weight')
    algoParams2.weight = algoParams.weight;
else
    algoParams2.weight = 0.5;
end

%% - algoParams.fm_init (Binary flag (1->initial guess 0->no initial guess) to enable the use of an initial guess for the field inhomogeneity term )
if isfield(algoParams, 'fm_init')
    algoParams2.fm_init = algoParams.fm_init;
else
    algoParams2.fm_init = 0;
end

if isfield(algoParams, 'dkg')
    algoParams2.dkg = algoParams.dkg;
else
    algoParams2.dkg = 15; % After dkg iterations, we may switch to a more homogeneous
                          % regularization, to achieve more smoothness in noise-only
                          % regions
end

if isfield(algoParams, 'SMOOTH_NOSIGNAL')
    algoParams2.SMOOTH_NOSIGNAL = algoParams.SMOOTH_NOSIGNAL;
else
    algoParams2.SMOOTH_NOSIGNAL = true; % Whether to "homogenize" the lambdamap after
                                        % some iterations, to get a smoother fieldmap in
                                        % low-signal regions
end

if isfield(algoParams, 'STARTBIG')
    algoParams2.STARTBIG = algoParams.STARTBIG;
else
    algoParams2.STARTBIG = true;
end