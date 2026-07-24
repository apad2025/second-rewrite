%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code: Example for performing fat-water separation with unipolar and
% bipolar readout gradient pulses using a graph-cut fat-water separation
% approach based on Hernando et al. technique doi: 10.1002/mrm.22177

% Copyright Jorge Campos 2025 - MIT License
%
% This code uses tools from the ISMRM fat-water toolbox. For this code to
% work, make sure to include the toolbox in matlab's directory.
% Phantom data was acquired in a 3T Philips scanner. Data was acquired
% enabling 'Delayed Reconstruction

% This code uses a phase unwrapping (used only to display some results
% without phase wraps). Performing phase unwrapping is optional, but if
% used, include in the code directory the function
% 'qualityGuidedUnwrapping' from ortier and Levesque, DOI 10.1002/mrm.26989, 2017, https://gitlab.com/veronique_fortier/Quality_guided_unwrapping

% This code uses the crameri perceptually uniform scientific colormaps to
% display some images. Add this toolbox to your matlab directory if needed https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps 
% If the crameri colormpas are not available in the directory, the code
% uses matlab colormaps


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except dataset
close all
clc

%% Include all directories in current folder

addpath(genpath(pwd));

%% Include ISMRM toolbox path

addpath(genpath('.\ISMRM_toolbox'));

%% Include crameri toolbox path

addpath(genpath('.\Crameri_color_maps'));

contents = dir('.\Crameri_color_maps');

if numel(contents) <= 3
    disp('The directory for the crameri colormaps is empty, script uses matlab colormaps');
    algoParams.crameri_colormap = 0;
else
    addpath(genpath('.\Crameri_color_maps'));
    algoParams.crameri_colormap = 1;
end

%% Name of dataset for postprocessing

name_experiment = '20240212_Experiment1_Bipolar_O_6_TE_15';

%% Selection of fat water separation approach
% Choose 'GC' to use graph-cut fat-water separation (Hernando et al. doi: 10.1002/mrm.22177)
% Choose 'Bipolar_GC' to use a modified version of the graph-cut fat-water
% separation approach that works with bipolar gradients

fat_water_separation_method = 'Bipolar_GC';

%% Input folder location (folder stores mat file with MR information for postprocessing)

input_folder = fullfile('.','Data/Input');

%% Output folder location (folder to save output)

output_folder = fullfile('.','Data/Output');

%% Selection of slices for postprocessing

sliceNb_ROI = 13; % Main slice for the postprocessing (fat-water separation is done in 2D over slice# = 'sliceNb_ROI')
extra_slices = 0; % Number of slices below and above sliceNb_ROI for postprocessing. If specific slices (that, for example, are not next to each other) need to be postprocessed, an array can be passed in the parameter imDataParams.sliceofint

%% Load the data from a mat file

data = load([input_folder '/' name_experiment  '.mat']); % Complex

MagnData(:,:,:,:)=double(data.dataScan.magn_rescaleFP);  %magnitude data without Philips scaling
realIm(:,:,:,:)=double(data.dataScan.real); %real MR signal
imIm(:,:,:,:)=double(data.dataScan.im); %imaginary MR signal
PhaseData(:,:,:,:)=double(data.dataScan.phaseReconstructed); %Phase data reconstructed from real and imaginary

voxelSize=data.dataScan.voxelSize; %voxel size [mm]
matrixSize=data.dataScan.matrixSize; %matrix size
gamma=267.513*10^6; % constant [rad/sT]
CF=data.dataScan.CF; % center frequency [Hz]
mainField=2*pi*CF/gamma; % main magnetic field [T]
TE=data.dataScan.TE; % timing of echo train length [s];
if length(TE)>1
    delta_TE=TE(2)-TE(1);
end

if isempty(PhaseData)==0
    iField(:,:,:,:) = MagnData(:,:,:,:).*exp(1i*PhaseData(:,:,:,:)); %--------------------------
end

%% Mask for phantom
% This piece of code calculates a binary mask for the big cylindrical
% enclosure

image_for_mask = MagnData(:,:,:,1);

% 'mask_bacground' analyzes a small portion of the large phantom
% compartement
mask_background = zeros(size(MagnData(:,:,sliceNb_ROI,1)));

mask_background(40:end-40,85:end-13) = 1;

%figure,imshowpair(mask_background,MagnData(:,:,sliceNb_ROI,1),"montage")

slice_background = MagnData(:,:,sliceNb_ROI,1);

rough_mask_filled1 = mean(slice_background(mask_background==1))-2.5*std(slice_background(mask_background==1))<MagnData(:,:,:,1);

BW2 = activecontour(image_for_mask,rough_mask_filled1,300);

se = strel('disk',10);

for kk = 1:size(BW2,3)
    mask_filled(:,:,kk) = imclose(BW2(:,:,kk),se);
end

%figure,imshowpair(mask_filled(:,:,sliceNb_ROI),MagnData(:,:,sliceNb_ROI,1),"montage")

% Saving mask into 'data' dataset
data.mask_filled = mask_filled;

%% Set of slices for postprocessing ()

data.sliceNb_ROI = sliceNb_ROI;

imDataParams.sliceofint = [sliceNb_ROI-extra_slices:sliceNb_ROI+extra_slices]; 
vec_slices = imDataParams.sliceofint;
slice_image = sliceNb_ROI;

%% Initialization for fat-water separation

% three-point Dixon - ISMRM toolbox
algoParams.species(1).name = 'water';
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
% Fat spectrum (6 resonance model for peanut oil) 
algoParams.species(1).frequency = 4.7;
algoParams.species(2).frequency = [0.80 1.20 2.00 2.66 4.21 5.20];
algoParams.species(2).relAmps = [0.087 0.694 0.128 0.004 0.039 0.048];

imDataParams.voxelSize=voxelSize;
imDataParams.FieldStrength=mainField;
imDataParams.PrecessionIsClockwise=1;
imDataParams.images=reshape((iField(:,:,:,1:end)),[matrixSize(1),matrixSize(2),matrixSize(3),1,length(TE(1:end))]);
imDataParams.TE=TE;

% Algorithm-specific parameters
algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
algoParams.range_r2star = [0 100]; % Range of R2* values
algoParams.NUM_R2STARS = 26;%11; % Numbre of R2* values for quantization
algoParams.range_fm = [-500 500];%[-1 1]; % Range of field map values
algoParams.NUM_FMS = 501; %50; % Number of field map values to discretize
algoParams.MAX_ITERS = 80; % Number of graph cut iterations
algoParams.SUBSAMPLE = 0; % Spatial subsampling for field map estimation (for speed)
algoParams.DO_OT = 0; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
algoParams.lambda = 0.05; % Regularization parameter
algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
algoParams.TRY_PERIODIC_RESIDUAL = 0;
algoParams.THRESHOLD = 0.01;

% Extra parameters created for the algorithm to do fat-water separation
% with bipolar readout gradients

imDataParams.mask_fwseparation = mask_filled; % Enables execution of fat-water separation only on the voxels within the
% binary mask

algoParams.slice_image = slice_image; % Value of slice to display the images to check intermediate results (mostly for debugging purposes)
algoParams.plot_debug = 1; % Binary flag (1->debugging 0->no debugging) to display images that show intermediate results (useful for debugging code)
algoParams.tik_reg = 0; % Binary flag to include Tikhonov regularization in the inverse problem to claculate phi and eps (1->yes regularization 0->no regularization). Regularization was not included for this paper, but it can enable future research about the the potential need and impact of regularization in these parameters 
algoParams.weight = 0.5; % Weight for the inverse problem to determine correction factors phi and eps. In this paper we tested 0.5
algoParams.fm_init = 0;% Binary flag (1->initial guess 0->no initial guess) to enable the use of an initial guess for the field inhomogeneity term 

if strcmp(fat_water_separation_method,'GC')

    % Memory allocation
    Fat_img = zeros([size(iField,1),size(iField,2),size(iField,3)]);
    Water_img = zeros([size(iField,1),size(iField,2),size(iField,3)]);
    FieldMap = zeros([size(iField,1),size(iField,2),size(iField,3)]);
    R2star = zeros([size(iField,1),size(iField,2),size(iField,3)]);

    % Fat-water separation
    for kk=1:length(vec_slices)

        imDataParams.sliceofint = vec_slices(kk);
        outParams = fw_i2cm1i_3pluspoint_hernando_graphcut( imDataParams, algoParams );

        Fat_img(:,:,vec_slices(kk)) = outParams.species(2).amps;
        Water_img(:,:,vec_slices(kk)) = outParams.species(1).amps;
        FieldMap(:,:,vec_slices(kk)) = outParams.fieldmap;
        R2star(:,:,vec_slices(kk)) = outParams.r2starmap;
    end


elseif strcmp(fat_water_separation_method,'Bipolar_GC')

    % Memory allocation
    Fat_img = zeros([size(iField,1),size(iField,2),size(iField,3)]);
    Water_img = zeros([size(iField,1),size(iField,2),size(iField,3)]);
    FieldMap = zeros([size(iField,1),size(iField,2),size(iField,3)]);
    R2star = zeros([size(iField,1),size(iField,2),size(iField,3)]);

    real_phase_error = zeros([size(iField,1),size(iField,2),size(iField,3)]);
    imag_phase_error = zeros([size(iField,1),size(iField,2),size(iField,3)]);

    %

    Fat_img_odd = zeros([size(iField,1),size(iField,2),size(iField,3)]);
    Water_img_odd = zeros([size(iField,1),size(iField,2),size(iField,3)]);

    Fat_img_even = zeros([size(iField,1),size(iField,2),size(iField,3)]);
    Water_img_even = zeros([size(iField,1),size(iField,2),size(iField,3)]);

    FieldMap_odd_even = zeros([size(iField,1),size(iField,2),size(iField,3)]);
    R2star_odd_even = zeros([size(iField,1),size(iField,2),size(iField,3)]);

    % Fat-water separation with bipolar readout gradient pulses

    outParams_DualGC = Function_Bipolar_GC( imDataParams, algoParams ,vec_slices);

    Fat_img(:,:,vec_slices) = outParams_DualGC.species(2).amps(:,:,vec_slices);
    Water_img(:,:,vec_slices) = outParams_DualGC.species(1).amps(:,:,vec_slices);
    FieldMap(:,:,vec_slices) = outParams_DualGC.fieldmap(:,:,vec_slices);
    R2star(:,:,vec_slices) = outParams_DualGC.r2starmap(:,:,vec_slices);

    real_phase_error(:,:,vec_slices) = outParams_DualGC.phi_map(:,:,vec_slices);
    imag_phase_error(:,:,vec_slices) = outParams_DualGC.eps_map(:,:,vec_slices);

    %

    Fat_img_odd(:,:,vec_slices) = outParams_DualGC.Fat_GC_odd(:,:,vec_slices);
    Water_img_odd(:,:,vec_slices) = outParams_DualGC.Water_GC_odd(:,:,vec_slices);

    Fat_img_even(:,:,vec_slices) = outParams_DualGC.Fat_GC_even(:,:,vec_slices);
    Water_img_even(:,:,vec_slices) = outParams_DualGC.Water_GC_even(:,:,vec_slices);

    FieldMap_odd_even(:,:,vec_slices) = outParams_DualGC.FieldMap_DualGC(:,:,vec_slices);
    R2star_odd_even(:,:,vec_slices) = outParams_DualGC.R2_DualGC(:,:,vec_slices);

end

%% Storing resulting images

data.Fat_img = Fat_img;
data.Water_img = Water_img;
data.FieldMap = FieldMap;
data.R2star = R2star;
if strcmp(fat_water_separation_method,'GC') | strcmp(fat_water_separation_method,'R2starIDEAL')
    data.real_phase_error = [];
    data.imag_phase_error = [];
elseif strcmp(fat_water_separation_method,'BipolarR2starIDEAL')
    data.real_phase_error = real_phase_error;
    data.imag_phase_error = imag_phase_error;
elseif strcmp(fat_water_separation_method,'Bipolar_GC')
    data.real_phase_error = real_phase_error;
    data.imag_phase_error = imag_phase_error;
    data.Fat_img_odd = Fat_img_odd;
    data.Water_img_odd = Water_img_odd;
    data.Fat_img_even = Fat_img_even;
    data.Water_img_even = Water_img_even;
    data.FieldMap_odd_even = FieldMap_odd_even;
    data.R2star_odd_even = R2star_odd_even;
end

%% Save results

save([output_folder '/' name_experiment '_' fat_water_separation_method '.mat'], ...
        'data','-v7.3');


