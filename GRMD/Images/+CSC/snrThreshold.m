% Calculate the volume-wide SNR threshold used to binarize the fat/water masks
function snr_thresh = snrThreshold(D, flags)
% Inputs:
%     D: data structure
% flags: processing structure
%
% Outputs:
%    snr_thresh: threshold to hand to Function_Bipolar_GC via algoParams
%
% Note: the threshold is a reduction over the whole volume, so it cannot be
% recalculated by a job that only holds a single slice. Cache it here and pass
% it in, otherwise slice-wise runs will not match a whole-volume run.

    % Format data to correct shape
    images = zeros(D.Size(1), D.Size(2), D.Size(3), 1, D.Size(4));
    images(:,:,:,1,:) = D.Data.Image;

    % Remove last echo, if necessary
    if any(strcmp(flags.cscorrection.method, {'IGC','MixedIGC','BipolarIGC'})) && rem(size(images,5),2)==1
        images = images(:,:,:,:,1:end-1);
    end

    % Set SNR threshold
    % Note: sum over coils (dim 4) and echoes (dim 5) to match the magnitude
    % definition in Function_Bipolar_GC's mask generation.
    mag = squeeze(sqrt(sum(abs(images).^2, [4 5])));
    mask_mag = mag > 0.1*max(mag(:));
    snr_thresh = prctile(mag(mask_mag),25);
end
