% Function: decomposeGivenFieldMapAndDampings_Bipolar_GC
%
% Description: estimate water/fat images given the nonlinear parameters
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

% % imDataParams.mask_fwseparation; % Enables execution of fat-water separation only on the voxels within the
% % % binary mask
% % 
% % algoParams.slice_image; % Value of slice to display the images to check intermediate results (mostly for debugging purposes)
% % algoParams.plot_debug; % Binary flag (1->debugging 0->no debugging) to display images that show intermediate results (useful for debugging code)
% % algoParams.tik_reg; % Binary flag to include Tikhonov regularization in the inverse problem to claculate phi and eps (1->yes regularization 0->no regularization). Regularization was not included for this paper, but it can enable future research about the the potential need and impact of regularization in these parameters 
% % algoParams.weight; % Weight for the inverse problem to determine correction factors phi and eps. In this paper we tested 0.5
% % algoParams.fm_init;% Binary flag (1->initial guess 0->no initial guess) to enable the use of an initial guess for the field inhomogeneity term 

%
%  - fieldmap: the estimated B0 field map
%  - r2starWater: the estimated water R2* map
%  - r2starFat: the estimated fat R2* map
%
% Returns: 
%  - amps: the amplitudes for all chemical species and coils
%  - remerror: fit error norm
%
% Author: Diego Hernando
% Date created: 
% Date last modified: August 18, 2011
% Modified by Jorge Campos
% Date: May 05, 2025

function [amps,remerror] = decomposeGivenFieldMapAndDampings_Bipolar_GC(imDataParams,algoParams,fieldmap,r2starWater,r2starFat)

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

try 
  ampW = algoParams.species(1).relAmps;
catch
  ampW = 1.0;
end

deltaF = imDataParams.deltaF;
relAmps = algoParams.species(2).relAmps;
images = imDataParams.images;
t = imDataParams.TE;

sx = size(images,1);
sy = size(images,2);
N = size(images,5);
C = size(images,4);

relAmps = reshape(relAmps,1,[]);

B1 = zeros(N,4);
B = zeros(N,4);

for n = 1:N/2
   B1(n,:) = [ampW*exp(1i*2*pi*deltaF(1)*t(n)),sum(relAmps(:).*exp(1i*2*pi*deltaF(2:end)*t(n))),0,0];
end

for n = N/2+1:N
    B1(n,:) = [0,0,ampW*exp(1i*2*pi*deltaF(1)*t(n)),sum(relAmps(:).*exp(1i*2*pi*deltaF(2:end)*t(n)))];
end

remerror = zeros(sx,sy);
for kx = 1:sx
  for ky= 1:sy
    s = reshape( squeeze(images(kx,ky,:,:,:)), [C N]).';

    B(:,1) = B1(:,1).*exp(1i*2*pi*fieldmap(kx,ky)*t(:) - r2starWater(kx,ky)*t(:));
    B(:,2) = B1(:,2).*exp(1i*2*pi*fieldmap(kx,ky)*t(:) - r2starFat(kx,ky)*t(:));

    B(:,3) = B1(:,3).*exp(1i*2*pi*fieldmap(kx,ky)*t(:) - r2starWater(kx,ky)*t(:));
    B(:,4) = B1(:,4).*exp(1i*2*pi*fieldmap(kx,ky)*t(:) - r2starFat(kx,ky)*t(:));

    amps(kx,ky,:,:) = B\s;

    if nargout > 1
      remerror(kx,ky) = norm(s - B*squeeze(amps(kx,ky,:,:)),'fro');
    end
  end
end

if nargout > 1

% % %     figure(7)
% % %     imagesc(remerror)
% % %     %caxis([0 2*10^(-7)])
% % %     
% % %     figure(8)
% % %     imshow3D(abs(amps(:,:,1)))
% % %     figure(9)
% % %     imshow3D(abs(amps(:,:,2)))
% % % 
% % %     figure(10)
% % %     imshow3D(abs(amps(:,:,3)))
% % %     figure(11)
% % %     imshow3D(abs(amps(:,:,4)))
end