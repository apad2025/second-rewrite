% Function: computeResidual_Bipolar_GC
%
% Description: compute fit residual for water/fat imaging
% 
% Parameters:
% Input: structures imDataParams and algoParams
%   - imDataParams.images: acquired images, array of size[nx,ny,1,ncoils,nTE]
%   - imDataParams.TEs: echo times (in seconds)
%   - imDataParams.fieldStrength: (in Tesla)
%
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

% % imDataParams.mask_fwseparation; % Enables execution of fat-water separation only on the voxels within the
% % % binary mask
% % 
% % algoParams.slice_image; % Value of slice to display the images to check intermediate results (mostly for debugging purposes)
% % algoParams.plot_debug; % Binary flag (1->debugging 0->no debugging) to display images that show intermediate results (useful for debugging code)
% % algoParams.tik_reg; % Binary flag to include Tikhonov regularization in the inverse problem to claculate phi and eps (1->yes regularization 0->no regularization). Regularization was not included for this paper, but it can enable future research about the the potential need and impact of regularization in these parameters 
% % algoParams.weight; % Weight for the inverse problem to determine correction factors phi and eps. In this paper we tested 0.5
% % algoParams.fm_init;% Binary flag (1->initial guess 0->no initial guess) to enable the use of an initial guess for the field inhomogeneity term 


%
% Returns: 
%  - residual: the residual, of size NUM_FMS X sx X sy
%
% Author: Diego Hernando
% Date created: August 13, 2011
% Date last modified: December 8, 2011
% Modified by Jorge Campos
% Date: May 05, 2025

function residual = computeResidual_Bipolar_GC(imDataParams, algoParams, VERBOSE)

if nargin < 3
    VERBOSE = 0;
end

images = imDataParams.images;

if isfield(imDataParams, 'PrecessionIsClockwise')
    precessionIsClockwise = imDataParams.PrecessionIsClockwise;
    % If precession is clockwise (positive fat frequency) simply conjugate data
    if precessionIsClockwise <= 0
        imDataParams.images = conj(imDataParams.images);
        imDataParams.PrecessionIsClockwise = 1;
    end
else
    precessionIsClockwise = 1;
end

deltaF = imDataParams.deltaF;
relAmps = algoParams.species(2).relAmps;
range_fm = algoParams.range_fm;
t = imDataParams.TE;
NUM_FMS = algoParams.NUM_FMS;

range_r2star = algoParams.range_r2star;
NUM_R2STARS = algoParams.NUM_R2STARS;


% Images are size sx X sy, N echoes, C coils
sx = size(images,1);
sy = size(images,2);
N = size(images,5);
C = size(images,4);

% Number of acquisitions
num_acq = size(images,6);

% Get VARPRO-formulation matrices for given echo times and chemical shifts 
Phi = getPhiMatrixMultipeak_Bipolar_GC(deltaF,relAmps,t);

psis = linspace( range_fm(1),range_fm(2),NUM_FMS );

% Compute residual
if length(range_r2star) == 2
    r2s = linspace(range_r2star(1),range_r2star(2),NUM_R2STARS);
else
    extra_step_r2s = algoParams.extra_step_r2s;
    r2s = [linspace(range_r2star(1),range_r2star(2),NUM_R2STARS),range_r2star(2)+extra_step_r2s:extra_step_r2s:range_r2star(3)];
end

% Precompute all projector matrices (one per field value) for VARPRO
if VERBOSE
    tic
    fprintf('\nPrecomputing all projector matrices...')
end
P = zeros(N*NUM_FMS,N,length(r2s));
for kr = 1:length(r2s)%NUM_R2STARS
    P1 = [];
    for k = 1:NUM_FMS
        Psi = diag(exp(j*2*pi*psis(k)*t - abs(t)*r2s(kr)));
        P1 = [P1;(eye(N)-Psi*Phi*pinv(Psi*Phi))];
    end
    P(:,:,kr) = P1;
end
if VERBOSE
    tttime = toc;
    fprintf('Done (%.2f sec)', tttime)
end

% Compute residual for all voxels and all field values
% Note: the residual is computed in a vectorized way, for increased speed
if VERBOSE
    tic
    fprintf('\nComputing residual for all voxels & field values...')
    reverseStr = '';
end

% Go line-by-line in the image to avoid using too much memory, while
% still reducing the loops significantly
residual = zeros(NUM_FMS,sx,sy);
for ka = 1:num_acq
    for ky = 1:sy
        if VERBOSE
            reverseStr = UpdatePercent(ky/sy*100, reverseStr);
        end

        temp = reshape(squeeze(permute(images(:,ky,:,:,:,ka),[1 2 3 5 4])),[sx N*C]).';
        temp = reshape(temp,[N sx*C]);
        for kr=1:length(r2s)%NUM_R2STARS
            temp2(:,:,kr) = reshape(sum(abs(reshape(P(:,:,kr)*temp,[N C*NUM_FMS*sx])).^2,1),[NUM_FMS C*sx]).';
            temp3(:,kr) = sum(reshape(temp2(:,:,kr),[C NUM_FMS*sx]),1);
        end
        [mint3,imint3] = min(temp3,[],2);

        residual(:,:,ky) = squeeze(squeeze(residual(:,:,ky)).' + reshape(mint3,[sx NUM_FMS])).';
    end
end
if VERBOSE
    tttime = toc;
    fprintf('Done (%.2f sec)', tttime)
end

% if algoParams.plot_debug
% figure(29)
% hold on
% grid minor
% title('In vivo abdomen')
% plot(psis,residual(:,28,38),'--k')
% 
% figure(611)
% hold on
% grid minor
% title('Comparison single and dual')
% plot(psis,residual(:,28,38),'--k')
% 
% end

end

% Update display percentage
function revstr = UpdatePercent(perc, revstr)
    msg = sprintf('%.2f percent. ', perc);
    fprintf([revstr, msg]);
    revstr = repmat(sprintf('\b'), 1, length(msg));
end