% Function name: fw_i2xm1c_3pluspoint_hernando_mixedfit.m
%
% Description: Fat-water separation using mixed magnitude/complex fitting.
%
% Hernando D, Hines CDG, Yu H, Reeder SB. Addressing phase errors in fat-water imaging 
% using a mixed magnitude/complex fitting method. Magn Reson Med; 2011.
% 
% Some properties:
%   - Image-space
%   - 2 species (water-fat)
%   - Mixed-fitting
%   - Multi-peak fat (pre-calibrated)
%   - Single-R2*
%   - Common water/fat phase
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
%   - algoParams.range_r2star = [0 0]; % Range of R2* values
%   - algoParams.range_fm = [-400 400]; % Range of field map values
%   - algoParams.NUM_ITERS = 40; % Number of descent iterations
%
% Output: structure outParams
%   - outParams.species(ii).name: name of the species (taken from algoParams)
%   - outParams.species(ii).amps: estimated water/fat images, size [nx,ny,ncoils] 
%   - outParams.r2star: R2* map (in s^{-1}, size [nx,ny])
%   - outParams.fieldmap: field map (in Hz, size [nx,ny])
%
% Author: Diego Hernando
% Date created: August 18, 2011
% Date last modified: February 10, 2012

function outParams = fw_i2xm1c_3pluspoint_hernando_mixedfit(imDataParams, algoParams, DEBUG)

if nargin < 3
    DEBUG = 0;
end

% Check validity of params, and set default algorithm parameters if not provided
[validParams,algoParams] = checkParamsAndSetDefaults_mixedfit(imDataParams,algoParams);
if validParams == 0
    disp('Exiting -- data not processed');
    outParams = [];
    return;
end

% Get data dimensions
[sx,sy,sz,C,N] = size(imDataParams.images);

% If more than one slice, pick central slice
if sz > 1
    disp('Multi-slice data: processing central slice');
    imDataParams.images = imDataParams.images(:,:,ceil(end/2),:,:);
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

%% Get images and recon parameters
im = double(imDataParams.images);
[sx,sy,sz,C,N] = size(im);
if isfield(algoParams, 'gyro')
    gyro = algoParams.gyro;
else
    gyro = 42.5774780505984;
end
freqs = [0; gyro*(algoParams.species(2).frequency(:) - algoParams.species(1).frequency(1))*(imDataParams.FieldStrength)];
M = length(freqs);

% Number of potentially phase-corrupted echoes
NUM_MAGN = algoParams.NUM_MAGN;

% Signal threshold for processing voxels (by default process all)
THRESHOLD = algoParams.THRESHOLD;

% Fat relative amplitudes
relAmps = algoParams.species(2).relAmps(:);

% Initial guess for R2* map
r2init = algoParams.r2starmap;

% Initial guess for B0 fieldmap
fminit = algoParams.fieldmap;

% Take the maximum allowed R2* value
MAXR2 = algoParams.range_r2star(2);

% Check whether we should use bounds at all
USE_BOUNDS = algoParams.use_bounds;
if MAXR2 < 0
    USE_BOUNDS = 0;
end

% Get the signal norm at each voxel (used for thresholding)
smagn = sqrt(sum(sum(abs(im).^2,5),4));
smagn = smagn/max(smagn(:));

% Create figure
if DEBUG
    figure;
end

% Initialize maps and iterate over voxels
ampsmap = zeros(sx,sy,2);
r2starmap = zeros(sx,sy);
fm = zeros(sx,sy);
phi = zeros(sx,sy);
for kx = 1:sx
    for ky = 1:sy
        % if kx==33333 && ky==3
        %   DEBUG = 1;
        % else
        %   DEBUG = 0;
        % end

        if smagn(kx,ky) > THRESHOLD
            [ampsmap(kx,ky,:),r2starmap(kx,ky),fm(kx,ky),phi(kx,ky)] = fit_Mixed_1R2_MP_SingleVoxel(reshape(squeeze(im(kx,ky,:,:,:)),[N C]), imDataParams.TE, freqs, relAmps,fminit(kx,ky),r2init(kx,ky), NUM_MAGN, DEBUG, MAXR2, USE_BOUNDS);
        end
    end
end

% CLose figure
if DEBUG
    close
end

% Put results in outParams structure
try
    outParams.species(1).name = algoParams.species(1).name;
    outParams.species(2).name = algoParams.species(2).name;
catch
    outParams.species(1).name = 'water';
    outParams.species(2).name = 'fat';
end

outParams.species(1).amps = ampsmap(:,:,1).*exp(1i*(phi));
outParams.species(2).amps = ampsmap(:,:,2).*exp(1i*(phi));
outParams.r2starmap = r2starmap;
outParams.fieldmap = fm;

%% Perform mixed fitting at each voxel
function [amps,r2,fm,phi] = fit_Mixed_1R2_MP_SingleVoxel(s, t, fs, relAmps, fminit, r0, NUM_MAGN, DEBUG, MAXR2, USE_BOUNDS)
% s = voxel signal, t = TEs, fs = frequencies, r0 = initial r2* map
% Initialize amplitudes
N = length(t); % # of echoes
M = length(fs); % # of fat & water frequencies
C = size(s,2); % # of coils
B = zeros(N,M);
t = t(:); % echo times

r0b = r0*ones(M,1);
for n = 1:N
    B(n,:) = exp(1i*2*pi*fs(:)*t(n) + 1i*2*pi*fminit*t(n) - r0b(:)*t(n)); % eq 1 (most of it)
end
B = [B(:,1) , sum(B(:,2:end)*diag(relAmps),2)]; % rest of eq 1

%amps0 = B(NUM_MAGN+1:end,:)\s(NUM_MAGN+1:end);
amps0 = B\s(:);

magn0 = abs(amps0);
[maxval,maxi] = max(magn0);
phi0 = angle(amps0(maxi));

magn0 = real(amps0*exp(-1i*phi0));

x0 = [magn0;r0;fminit;phi0];

if DEBUG
    x0 %#ok<NOPRT>
end

% Set some bounds
MAXS = max(abs(s(:)));
lb = [-2*MAXS*exp(MAXR2*t(1))*ones(size(amps0(:)));zeros(size(r0(:)));-1e3;-2*pi ];
ub = [2*MAXS*exp(MAXR2*t(1))*ones(size(amps0(:)));MAXR2*ones(size(r0(:)));1e3; 2*pi];

if ~USE_BOUNDS
    lb = [];
    ub = [];
end

% Run it!
options = optimset('Display','off','Jacobian','on','DerivativeCheck','off','MaxIter',15);
x = lsqnonlin( @(y)fitError_Mixed_1R2_MP(y,s,t,N,M,C,fs,relAmps,NUM_MAGN,DEBUG),x0,lb,ub,options);
amps = x(1:2);
r2 = x(3);
fm = x(4);
phi = x(5);

%% Compute fit error and Jacobian
function [fiter, J] = fitError_Mixed_1R2_MP(x,s,t,N,M,C,freqs,relAmps,NUM_MAGN,DEBUG)

% Grab the parameters
amps = x(1:2); % magn0
r2 = x(3); % r0
fm = x(4); % fminit
phi = x(5); % phi0

% t1 are the magnitude fitting TEs, t2 are the complex fitting TEs
t1 = t(1:NUM_MAGN);
t2 = t(NUM_MAGN+1:end);

% Find the error
r2b = r2*ones(length(freqs),1);
B = zeros(N,M);
for n = 1:N
    B(n,:) = exp(1i*2*pi*(freqs(:)+fm)*t(n) - r2b(:)*t(n)); % same as B in fit_Mixed_1R2_MP_SingleVoxel, just slightly simplified visually
    cf(n,1) = sum(relAmps.*exp(1i*2*pi*freqs(2:end)*t(n)));
end
B = [B(:,1), sum(B(:,2:end)*diag(relAmps),2)];

shat = B*amps*exp(1i*phi);

% Compute mixed fitting fit error
shat1 = abs(shat(1:NUM_MAGN));
s1 = abs(s(1:NUM_MAGN));
shat2 = shat(NUM_MAGN+1:end);
s2 = s(NUM_MAGN+1:end);

fiter1 = shat1(:) - s1(:); % abs(shat(1:NUM_MAGN)) - abs(s(1:NUM_MAGN))
fiter2 = shat2(:) - s2(:); % shat(NUM_MAGN+1:end) - s(NUM_MAGN+1:end)
fiter2 = [real(fiter2); imag(fiter2)];

cf1 = cf(1:NUM_MAGN);
cf2 = cf(NUM_MAGN+1:end);

fiter = [fiter1;fiter2];

% If in DEBUG mode, display some stuff
if DEBUG
    disp('')
    disp(['fieldmap = ' num2str(fm)]);
    disp(['r2 = ' num2str(r2.')]);
    disp(['amps = ' num2str(abs(amps).')])
    disp(['fit error = ' num2str(norm(shat(:) - s(:)))])
    plot(t*1000,abs(s),'k','DisplayName','s');hold on;
    plot(t*1000,abs(shat),'g','DisplayName','$$\hat{s}$$');hold off;
    xlabel('TE (msec)');
    legend('Interpreter','latex')
    drawnow;pause(0.02);
end

% Get the Jacobian, if needed
if nargout > 1
    % Magnitude portion
    denom = sqrt(amps(1)^2 + amps(2)^2*abs(cf1).^2 + 2*amps(1)*amps(2)*real(cf1));
    J1(:,1) = exp(-r2*t1).*(amps(1) + amps(2)*real(cf1))./denom;
    J1(:,2) = exp(-r2*t1).*(amps(2).*abs(cf1).^2 + amps(1)*real(cf1))./denom;
    J1(:,3) = -t1.*exp(-r2*t1).*denom;
    J1(:,4) = 0;
    J1(:,5) = 0;

    % Complex portion
    J2(:,1) = exp(1i*2*pi*fm*t2).*exp(-r2*t2)*exp(1i*phi);
    J2(:,2) = cf2.*exp(1i*2*pi*fm*t2).*exp(-r2*t2)*exp(1i*phi);
    J2(:,3) = -t2.*shat2;
    J2(:,4) = 1i*2*pi*t2.*shat2;
    J2(:,5) = 1i*shat2;
    J2 = [real(J2);imag(J2)];

    % Combine all data
    J = [J1;J2];
end