% Function: estimateR2starGivenFieldmap_Bipolar_GC
%
% Description: estimate R2* map, given the fieldmap
% 
% Parameters:
% Input: structures imDataParams and algoParams
%   - imDataParams.images: acquired images, array of size[nx,ny,1,ncoils,nTE]
%   - imDataParams.TEs: echo times (in seconds)
%   - imDataParams.fieldStrength: (in Tesla)
%
%   - algoParams.species(ii).name = string containing the name of the species
%   - algoParams.species(ii).frequency = frequency shift in ppm of each peak within species ii
%   - algoParams.species(ii).relAmps = relative amplitude (sum normalized to 1) of each peak within species ii
%   Example
%      - algoParams.species(1).name = 'water' % Water
%      - algoParams.species(1).frequency = [0] 
%      - algoParams.species(1).relAmps = [1]   
%      - algoParams.species(2).name = 'fat' % Fat
%      - algoParams.species(2).frequency = [3.80, 3.40, 2.60, 1.94, 0.39, -0.60]
%      - algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048]
%   - algoParams.range_r2star = [0 0]; % Range of R2* values
%   - algoParams.NUM_R2STARS = 1; % Numbre of R2* values for quantization
% 

% % imDataParams.mask_fwseparation; % Enables execution of fat-water separation only on the voxels within the
% % % binary mask
% % 
% % algoParams.slice_image; % Value of slice to display the images to check intermediate results (mostly for debugging purposes)
% % algoParams.plot_debug; % Binary flag (1->debugging 0->no debugging) to display images that show intermediate results (useful for debugging code)
% % algoParams.tik_reg; % Binary flag to include Tikhonov regularization in the inverse problem to claculate phi and eps (1->yes regularization 0->no regularization). Regularization was not included for this paper, but it can enable future research about the the potential need and impact of regularization in these parameters 
% % algoParams.weight; % Weight for the inverse problem to determine correction factors phi and eps. In this paper we tested 0.5
% % algoParams.fm_init;% Binary flag (1->initial guess 0->no initial guess) to enable the use of an initial guess for the field inhomogeneity term 

%  - fm: the estimated B0 field map
%
% Returns: 
%  - r2starmap: the estimated R2* map
%  - residual: fit error residual
%
% Author: Diego Hernando
% Date created: 2009
% Date last modified: August 18, 2011
% Modified by Jorge Campos
% Date: May 05, 2025


function [r2starmap,residual] = estimateR2starGivenFieldmap_Bipolar_GC (imDataParams, algoParams, fm)


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

range_r2star = algoParams.range_r2star;
NUM_R2STARS = algoParams.NUM_R2STARS;
deltaF = imDataParams.deltaF;
relAmps = algoParams.species(2).relAmps;
images = imDataParams.images;
t = imDataParams.TE;
t = reshape(t,[],1);

sx = size(images,1);
sy = size(images,2);
C = size(images,4);
N = size(images,5);
num_acqs = size(images,6);

images = reshape(permute(images,[1 2 5 4 6 3]),[sx sy N C*num_acqs]);

% Undo effect of field map
images2 = zeros(size(images));
for kt = 1:N
    for kc = 1:C*num_acqs
      images2(:,:,kt,kc) = squeeze(images(:,:,kt,kc)).*exp(-1i*2*pi*fm*t(kt));
    end
end

% Compute residual as a function of r2
if length(range_r2star) == 2
    r2s = linspace(range_r2star(1),range_r2star(2),NUM_R2STARS);
else
    extra_step_r2s = algoParams.extra_step_r2s;
    r2s = [linspace(range_r2star(1),range_r2star(2),NUM_R2STARS),range_r2star(2)+extra_step_r2s:extra_step_r2s:range_r2star(3)];
    NUM_R2STARS = length(r2s);
end
%r2s = linspace(range_r2star(1),range_r2star(3),NUM_R2STARS);
Phi = getPhiMatrixMultipeak_Bipolar_GC(deltaF,relAmps,t);
P = [];
for k = 1:NUM_R2STARS
  Psi = diag(exp(-r2s(k)*t));
  P = [P;(eye(N)-Psi*Phi*pinv(Psi*Phi))];
end

% Compute residual for all voxels and all field values
% Note: the residual is computed in a vectorized way, for increased speed
residual = zeros(sx,sy,NUM_R2STARS);

% Go line-by-line in the image to avoid using too much memory, while
% still reducing the loops significantly
% $$$ disp('Calculating residual...')
for kx = 1:sx
  temp = reshape(squeeze(images2(kx,:,:,:)),[sy N*C*num_acqs]).';
  temp = reshape(temp,[N sy*C*num_acqs]);
  temp2 = reshape(sum(abs(reshape(P*temp,[N C*num_acqs*NUM_R2STARS*sy])).^2,1),[NUM_R2STARS C*num_acqs*sy]).';
  temp2 = sum(reshape(temp2,[C*num_acqs NUM_R2STARS*sy]),1);
  residual(kx,:,:) = reshape(temp2,[sy NUM_R2STARS]);
end
%  residual = shiftdim(residual,2);
% $$$ disp('done computing residual.');

[minres,iminres] = min(residual,[],3);
r2starmap = r2s(iminres);

% if algoParams.plot_debug
% figure(412)
% hold on
% grid minor
% pos_x = 28;
% pos_y = 38;
% title('R2* map residues')
% plot(r2s,squeeze(residual(pos_x,pos_y,:)))
% scatter(r2starmap(pos_x,pos_y),minres(pos_x,pos_y))
% 
% figure(52)
% imag_plot = r2starmap;
% imagesc(imag_plot)
% crameri oslo
% caxis([0 150])
% end