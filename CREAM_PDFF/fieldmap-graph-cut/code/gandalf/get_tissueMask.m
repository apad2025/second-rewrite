function tissueMask = get_tissueMask(signal, airSignalThreshold_percent)
% tissueMask = get_tissueMask(signal, airSignalThreshold_percent)

    % Combine if multiple echoes
    echoMIP = get_echoMIP(signal);

    threshold = airSignalThreshold_percent/100 * max(echoMIP(:));

    tissueMask = echoMIP >= threshold .* ones(size(echoMIP));

    % Fill holes
    for z = 1:size(tissueMask,3)
        tissueMask(:,:,z) = imfill(tissueMask(:,:,z), "holes");
    end

    %% Remove large noise regions
    % Find connected components
    CC = bwconncomp(tissueMask); 

    % Only proceed if more than one region - JND 4/1/25
    if CC.NumObjects > 1
        % Determine size of each connected region
        numRegPixels = cellfun(@numel,CC.PixelIdxList); 
    
        % Determine total number of pixels
        % numPixels = numel(tissueMask);
        numPixels = sum(tissueMask,"all");
    
        % Include regions larger than 10% of total pixels
        includedRegs = numRegPixels > numPixels/10;
    
        % Combine regions
        tissueMask = false(size(tissueMask));
        for i = 1:length(numRegPixels)
            if includedRegs(i)
                tissueMask(CC.PixelIdxList{i}) = true;
            end
        end
    end
end