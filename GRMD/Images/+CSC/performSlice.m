% Correct chemical shift for a single slice
function [D, Data] = performSlice(D, flags, sl, snr_thresh)
% Inputs:
%             D: data structure
%         flags: processing structure
%            sl: slice to correct
%    snr_thresh: volume-wide SNR threshold (see CSC.snrThreshold)
%
% Outputs:
%             D: data structure, with the echo times updated to match
%          Data: corrected data for the requested slice
%
% Note: this is the BipolarIGC branch of CSC.perform restricted to one slice,
% so that the slices can be run as independent jobs and stacked afterwards by
% CSC.assemble.

    fprintf('\nPerforming chemical shift correction using %s on slice %i...\n', flags.cscorrection.method, sl)

    if ~strcmp(flags.cscorrection.method, 'BipolarIGC')
        error('Slice-wise chemical shift correction has only been added for BipolarIGC.')
    end

    if sl < 1 || sl > D.Size(3)
        error('Slice %i falls outside of the dataset (%i slices).', sl, D.Size(3))
    end

    % Algorithm-specific parameters
    algoParams = CSC.grabPars(flags);
    algoParams.parallel   = false;      % one slice per job, so no pool is needed
    algoParams.snr_thresh = snr_thresh; % keep the threshold identical across jobs

    % Isolate flags
    verboseFLAG = flags.verbose;

    % Format data to correct shape
    images = zeros(D.Size(1), D.Size(2), D.Size(3), 1, D.Size(4));
    images(:,:,:,1,:) = D.Data.Image;

    % Remove last echo, if necessary
    if rem(size(images,5),2)==1
        images = images(:,:,:,:,1:end-1);
        D.TE = D.TE(1:end-1);
        D.deltaTE = mean(diff(D.TE));
        D.Size(4) = size(images,5);
    end

    % Create data structure
    dataParams = struct('FieldStrength', D.B0, ...
                                   'TE', D.TE, ...
                'PrecessionIsClockwise', 1);

    % Grab extra parameters
    dataParams.voxelSize = D.VoxelSize;
    dataParams.images = images;
    dataParams.mask_fwseparation = 1;

    outParams = Function_Bipolar_GC(dataParams, algoParams, sl, verboseFLAG);

    % Combine final data
    % Note: kept as [nx ny 1 ...] so that CSC.assemble can stack along the third dimension
    Data.Image = reshape(outParams.corrected_bipolar_signal(:,:,sl,:,:), [D.Size(1), D.Size(2), 1, size(images,5)]); % bipolar signal transformed into unipolar equivalent
    Data.Water = outParams.species(1).amps(:,:,sl);
    Data.Fat = outParams.species(2).amps(:,:,sl);
    Data.TotalField = outParams.fieldmap(:,:,sl);
    Data.R2Star = outParams.r2starmap(:,:,sl);
    Data.Phi = outParams.phi_map(:,:,sl); % related to phase modulation due to bipolar readout
    Data.Epsilon = outParams.eps_map(:,:,sl); % related to amplitude modulation due to bipolar readout
    Data.BipolarError = outParams.bipolar_error_map_theta(:,:,sl); % phi - i*eps;
    Data.Correction = outParams.total_correction(:,:,sl); % correction to remove bipolar induced effects, e^(i*BipolarError)

    Data.WaterOdd = outParams.Water_GC_odd(:,:,sl);
    Data.FatOdd = outParams.Fat_GC_odd(:,:,sl);
    Data.WaterEven = outParams.Water_GC_even(:,:,sl);
    Data.FatEven = outParams.Fat_GC_even(:,:,sl);
    Data.TotalFieldDualGC = outParams.FieldMap_DualGC(:,:,sl);
    Data.R2StarDualGC = outParams.R2_DualGC(:,:,sl);
end
