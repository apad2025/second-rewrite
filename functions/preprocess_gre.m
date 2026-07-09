function D = preprocess_gre(D, flags)
%PREPROCESS_GRE Prepare a multi-echo GRE volume for fat/water separation.
%
%   D = preprocess_gre(D, flags)
%
%   Performs, in order:
%     1. Optional k-space zero-filling / interpolation (flags.zipped)
%     2. Weighted magnitude (echo MIP, normalised to [0,1])
%     3. Binary tissue mask (GenMask) + z-orientation check
%     4. Field-of-view crop to a padded square ROI (ROIisol)
%
%   flags fields (all optional):
%       .zipped  : logical, zero-fill in k-space to reduce ghosting (default false)
%       .verbose : logical, print progress (default true)
%       .nobone  : logical, remove bone/marrow from the mask (default false)
%
%   Condensed from the "Preprocessing" section of the original DogAnalysis.m.
%   The interactive through-plane slice-trimming dialog has been removed so the
%   pipeline can run unattended; in-plane FOV cropping (ROIisol) is retained.

    if ~isfield(flags, 'zipped'),  flags.zipped  = false; end
    if ~isfield(flags, 'verbose'), flags.verbose = true;  end
    if ~isfield(flags, 'nobone'),  flags.nobone  = false; end

    % ---- 1. Zero-fill (interpolate) in k-space to prevent ghosting ----
    if flags.zipped && ~D.Flags.Interpolated
        psize = [D.Size(1), D.Size(2), 0, 0]; % pad size
        img = padarray(D.Data.Image, psize/2, 0, 'both');

        % Forward to k-space
        kSpace = zeros(size(img), 'like', D.Data.Image);
        for j = 1:D.Size(4)
            for i = 1:D.Size(3)
                kSpace(:,:,i,j) = ifftshift(ifft2(ifftshift(img(:,:,i,j))));
            end
        end

        if flags.verbose, fprintf('\nZero-filling data...'); tic; end
        kSpace = padarray(kSpace, psize, 0, 'both');

        % Back to image space
        D.Data.Image = zeros(size(kSpace), 'like', kSpace);
        for j = 1:D.Size(4)
            for i = 1:D.Size(3)
                D.Data.Image(:,:,i,j) = fftshift(fft2(fftshift(kSpace(:,:,i,j))));
            end
        end
        clear kSpace

        % Remove excess padding
        D.Data.Image = D.Data.Image(D.Size(1)+1:end-D.Size(1), D.Size(2)+1:end-D.Size(2), :, :);

        % Update data information
        D.VoxelSize(1:2) = D.VoxelSize(1:2).*(D.Size(1:2)./size(D.Data.Image,[1 2]));
        D.Size = size(D.Data.Image);
        D.Flags.Interpolated = true;
        if flags.verbose, tm = toc; fprintf('Done (%0.2f sec)', tm); end
    end

    % ---- 2. Weighted magnitude (echo MIP) ----
    if ~isfield(D.Data, 'WeightedMagnitude')
        if flags.verbose, fprintf('\nCalculating weighted magnitude...'); tic; end
        mag = get_echoMIP(D.Data.Image);
        D.Data.WeightedMagnitude = (mag - min(mag,[],'all'))./(max(mag,[],'all') - min(mag,[],'all')); % [0, 1]
        clear mag
        if flags.verbose, tm = toc; fprintf('Done (%0.2f sec)', tm); end
    end

    % ---- 3. Tissue mask + orientation check ----
    if ~isfield(D.Data, 'Mask')
        D.Data.Mask = GenMask(D.Data.Image, D.Data.WeightedMagnitude, flags.verbose, 3, flags.nobone);

        % Flip so the larger cross-section is at slice 1 (consistent B0 sign)
        if sum(D.Data.Mask(:,:,1),'all') > sum(D.Data.Mask(:,:,end),'all')
            D.Data.Image = flip(D.Data.Image, 3);
            D.Data.Mask = flip(D.Data.Mask, 3);
            D.Data.WeightedMagnitude = flip(D.Data.WeightedMagnitude, 3);
            D.B0Direction(3) = -D.B0Direction(3);
        end
    end

    % ---- 4. Crop field of view to a padded square ROI ----
    if ~D.Flags.Trimmed.InPlane
        if flags.verbose, fprintf('\nTrimming field of view...'); tic; end
        D = ROIisol(D);
        if flags.verbose, tm = toc; fprintf('Done (%0.2f sec)', tm); end
    end
end
