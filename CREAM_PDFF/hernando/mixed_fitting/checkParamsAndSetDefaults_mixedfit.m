%% Function: checkParamsAndSetDefaults_mixedfit
%%
%% Description: Check validity of input parameters and set defaults for unspecified parameters
%%
%% Input:
%%   - imDataParams: TEs, images and field strength
%%   - algoParams: algorithm parameters
%%
%% Output:
%%   - validParams: binary variable (0 if parameters are not valid for this algorithm)
%%   - algoParams2: "completed" algorithm parameter structure (after inserting defaults for unspecified parameters)
%% 
%%
%% Author: Diego Hernando
%% Date created: August 19, 2011
%% Date last modified: November 10, 2011
%%

function [validParams,algoParams2] = checkParamsAndSetDefaults_mixedfit( imDataParams,algoParams )

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

%% Whether to use constrained estimation with set bounds (e.g., on R2* estimates)
if isfield(algoParams, 'use_bounds')
    algoParams2.use_bounds = algoParams.use_bounds;
else
    algoParams2.use_bounds = 1;
end

%% Initial guess for R2* map
if isfield(algoParams, 'r2starmap')
    algoParams2.r2starmap = algoParams.r2starmap;
else
  algoParams2.r2starmap = zeros(size(imDataParams.images(:,:,1,1,1,1)));
  disp('No initial guess for R2* map provided. Initializing with zeroes');
end

%% Initial guess for field map
if isfield(algoParams, 'fieldmap')
    algoParams2.fieldmap = algoParams.fieldmap;
else
    algoParams2.fieldmap = zeros(size(imDataParams.images(:,:,1,1,1,1)));
  disp('No initial guess for B0 field map provided. Initializing with zeroes');
end

%% Number of echoes with potentially corrupted phase
if isfield(algoParams, 'NUM_MAGN')
    algoParams2.NUM_MAGN = algoParams.NUM_MAGN;
else
    algoParams2.NUM_MAGN = 1;
end

%% Signal threshold for processing voxels (by default process all)
if isfield(algoParams, 'THRESHOLD')
    algoParams2.THRESHOLD = algoParams.THRESHOLD;
else
    algoParams2.THRESHOLD = 0.0;
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

%%   - algoParams.lambda = 0.05; % Regularization parameter
if isfield(algoParams, 'lambdamap')
    algoParams2.lambdamap = algoParams.lambdamap;
else
    algoParams2.lambdamap = ones(size(imDataParams.images(:,:,1,1,1,1)));
end

%%   - algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
if isfield(algoParams, 'LMAP_EXTRA')
    algoParams2.LMAP_EXTRA = algoParams.LMAP_EXTRA;
else
    algoParams2.LMAP_EXTRA = zeros(size(imDataParams.images(:,:,1,1,1)));
end

%%   - algoParams.TRY_PERIODIC_RESIDUAL = 0; % Take advantage of periodic residual if uniform TEs (will change range_fm)  
if isfield(algoParams, 'TRY_PERIODIC_RESIDUAL')
    algoParams2.TRY_PERIODIC_RESIDUAL = algoParams.TRY_PERIODIC_RESIDUAL;
else
    algoParams2.TRY_PERIODIC_RESIDUAL = 0;
end

%%   - imDataParams.PrecessionIsClockwise (1 = fat has positive frequency; -1 = fat has negative frequency)
if isfield(algoParams, 'PrecessionIsClockwise')
    imDataParams2.PrecessionIsClockwise = imDataParams.PrecessionIsClockwise;
    if imDataParams2.PrecessionIsClockwise <= 0
        imDataParams2.PrecessionIsClockwise = -1;
    end
else
    algoParams2.PrecessionIsClockwise = -1;
end

