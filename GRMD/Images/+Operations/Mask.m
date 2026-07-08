% Generate Mask
function msk = Mask(img, mag, verboseFLAG, thresh, boneFLAG, noiseSTD)
    if verboseFLAG, fprintf('\nGenerating boolean mask...'); tic; end

    % Generate basic 3D mask
    if ndims(img) > ndims(mag)
        % Combine if multiple echoes
        mag = get_echoMIP(img);
    else
        mag = abs(img);
    end
    mask1 = mag >= (thresh/100*max(mag(:))).*ones(size(mag));

    % Fill holes
    for z = 1:size(mask1,3)
        mask1(:,:,z) = imfill(mask1(:,:,z), "holes");
    end

    %% Remove large noise regions
    % Find connected components
    CC = bwconncomp(mask1);

    % Only proceed if more than one region
    if CC.NumObjects > 1
        numRegPixels = cellfun(@numel,CC.PixelIdxList); % size of each connected region
        numPixels = sum(mask1,"all"); % total number of pixels

        % Include regions larger than 10% of total pixels
        includedRegs = numRegPixels > numPixels/10;

        % Combine regions
        mask1 = false(size(mask1));
        for i = 1:length(numRegPixels)
            if includedRegs(i)
                mask1(CC.PixelIdxList{i}) = true;
            end
        end
    end

    % Find initial contour
    mask_basic = mag > max(mag(:))./10;
    mask_basic = imfill(mask_basic, "holes");
    mask_basic = imdilate(mask_basic, strel('disk', 20));

    % Obtain refined mask
    mask2 = activecontour(mag, mask_basic, 'Chan-vese');
    mask2 = imfill(mask2, "holes");

    % Combine
    mask = mask1+mask2 > 0;

    % Remove bones
    if boneFLAG
        mask = MaskBone(mask, img);
    end

    % Incorporate Noise STD
    if nargin == 6
        maskNSTD = noiseSTD < quantile(noiseSTD(:),0.135); % Only keep points in lowest quartile
        mask = mask - ((mask - maskNSTD) > 0);
    end

    for i = 1:50
        mask(:,:,i) = imfill(mask(:,:,i),'holes');
    end

    % Remove any sections significantly smaller than main
    CC = bwconncomp(mask);
    NumPixels = cellfun(@numel, CC.PixelIdxList);
    msk = false(size(mask));
    for i = 1:length(NumPixels)
        if NumPixels(i) > max(NumPixels)*0.05
            msk(CC.PixelIdxList{i}) = true;
        end
    end

    % Final revisions
    msk = imfill(msk, "holes");
    msk = imerode(msk, strel('disk', 1));
    msk = imdilate(msk, strel('disk', 1));
    msk = msk > 0; % just to ensure data is logical

    if verboseFLAG, tm = toc; fprintf('Done (%0.2f sec)', tm); end

    function mask_new = MaskBone(mask, IM)
        % Generate erosion & dilation shapes
        sphere4 = strel('sphere', 4);
        sphere1 = strel('sphere', 1);
        cube1 = strel('cube', 1);

        % Erode current mask to remove low signal/partial volume artifacts around skin layer
        n = size(mask);
        mask_er = imerode(mask, sphere4);
        mask_bone_total = false(n);
        tot = zeros(n);

        for eco = size(IM,4):-1:1
            im = abs(IM(:,:,:,eco));

            % Generate mask with bone
            mask_inv = im < max(im,[],'all')*0.1;
            mask_inv_er = mask_er.*mask_inv == 1;

            % Remove any regions that cross from one side to the other
            cc = bwconncomp(mask_inv_er);
            numpixels = cellfun(@numel, cc.PixelIdxList);
            mask_lr = false(n);
            for j = 1:length(numpixels)
                mask_temp = false(n);
                mask_temp(cc.PixelIdxList{j}) = true;

                if sum(mask_temp(:,1:round(n(2)/2),:),'all')==sum(mask_temp,'all') || sum(mask_temp(:,round(n(2)/2):end,:),'all')==sum(mask_temp,'all') % only on left or right side
                    mask_lr(cc.PixelIdxList{j}) = true;
                end
            end

            % Add previous mask
            mask_lr = (mask_bone_total + mask_lr) > 0;

            % Fill holes
            mask_filled = imdilate(imdilate(mask_lr, sphere1), cube1);
            for j = 1:n(3)
                mask_filled(:,:,j) = imfill(mask_filled(:,:,j), "holes");
            end
            mask_filled = imerode(imerode(mask_filled, sphere1), cube1);

            % Add previous mask
            mask_filled = (mask_bone_total + mask_filled) > 0;

            % Remove any super small regions
            cc = bwconncomp(mask_filled);
            numpixels = cellfun(@numel, cc.PixelIdxList);
            mask_bone = false(n);
            for j = 1:length(numpixels)
                if numpixels(j) > round(max(numpixels)*0.02)
                    mask_bone(cc.PixelIdxList{j}) = true;
                end
            end

            % Determine which are the two bones (assume the two largest clusters)
            cc = bwconncomp(mask_bone);
            numpixels = cellfun(@numel, cc.PixelIdxList);
            [~,boneIdx1] = max(numpixels);
            NumPixelstmp = numpixels;
            NumPixelstmp(boneIdx1) = 1;
            [~,boneIdx2] = max(NumPixelstmp);
            clear NumPixelstmp

            % Mark all regions aside from bones
            regs2check = 1:length(numpixels);
            regs2check([boneIdx1,boneIdx2]) = NaN;
            regs2check = regs2check(~isnan(regs2check));

            % Determine how much overlap, if any, there is with main bone mask
            mask_bone_overlap = false(n);
            if eco < 7
                for j = regs2check
                    mask_temp = false(n);
                    mask_temp(cc.PixelIdxList{j}) = true;
                    overlap = sum(mask_temp.*mask_bone_total,'all')/numpixels(j);

                    if overlap > 0
                        mask_bone_overlap(cc.PixelIdxList{j}) = true;
                    end
                end
            else
                % Only save two largest bones
                mask_bone_overlap(cc.PixelIdxList{boneIdx1}) = true;
                mask_bone_overlap(cc.PixelIdxList{boneIdx2}) = true;
            end

            % Add to main mask
            mask_bone_total = (mask_bone_total + mask_bone_overlap) > 0;
            tot = tot + 1.*mask_bone_overlap;
        end
        mask_new = mask.*(1-mask_bone_total) > 0;
    end
end