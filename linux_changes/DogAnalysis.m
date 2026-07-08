%% Select data & analysis pipeline
clear; clc; close all
% function DogAnalysis(dog, date)
% error('min cut/=max-flow in max_flow')
% Main path
% pth_main = pth + "_Project DMDiv\Dog Data";                       % legacy, unused on Linux
pth_data = "/scratch/user/apad/Fat_water_separation/DICOM_Files";   % location of DD.mat (MUST match dogexplorer's mainpth)
pth_code = "/scratch/user/apad/GRMD";                               % code root (kept for reference)

% Select dog & date (high PDFF are 4,3 & 2,3) (3,3 isn't finished being processed)
% 1 - Waylon
% 2 - Sushi
% 3 - Selene
% 4 - Aphrodite
% 5 - EOS
dog = 2; 
date = 5;

% Select flags
flags = struct('trimmed', true, ...
                'zipped', true, ...
     'bipolarcorrection', struct('method', 'MEDI'), ...
               'verbose', true, ...
                  'plot', false, ...
                'nobone', false, ...
            'unwrapping', struct('method', 'GraphCuts', ...
                                 'subsample', 1, ...
                                 'corrected', true), ...
          'cscorrection', struct('method', 'BipolarIGC', ...
                              'subsample', 1), ...
             'bfremoval', struct('D2', struct('method','V-SHARP - SEPIA'), ...
                                 'D3', struct('method','V-SHARP - SEPIA')), ...
       'dipoleinversion', struct('method', 'MEDI'));
saveFLAG = true;
cscFLAG = true; % chemical shift correction instead of unwrap + echo combination
debugFLAG = false;

% Remove subsampling if not graph-cuts
if ~strcmp(flags.unwrapping.method, 'GraphCuts'), flags.unwrapping = rmfield(flags.unwrapping, 'subsample'); end

% Remove CSC correction method if applicable
if ~cscFLAG, flags.cscorrection = struct('method', false); end

% Force 3D background field removal to be V-SHARP (STISuite) if that is used for 2D
if strcmp(flags.bfremoval.D2.method, 'V-SHARP - STISuite')
    flags.bfremoval.D3.method = 'V-SHARP - STISuite';
    jointD23FLAG = true;
else
    jointD23FLAG = false;
end

%% Load data
[snames, DD, HDR, pth_1H, plrange, vout, filelocs] = DogInitialize(pth_data, dog, date, flags, cscFLAG);

% Extract individual saves
for i = 1:length(vout)
    D = vout{i};
    if filelocs.Raw == i
        D_raw = D;
    elseif filelocs.Preprocessed == i
        D_preproc = D; 
    elseif filelocs.BipolarCorrected == i
        D_bpc = D;
    elseif filelocs.UnwrappedPhase == i
        D_unwrap = D; 
    elseif filelocs.CorrectedCS == i
        D_cscorr = D;
    elseif filelocs.LocalFieldD2 == i
        D_2Dlocfield = D; 
    elseif filelocs.LocalFieldD3 == i || filelocs.LocalFieldD23 == i
        D_3Dlocfield = D; 
    elseif filelocs.SusceptibilityMap == i
        D_sus = D; 
    end
end
clear vout filelocs

% Remove data from raw if zipped flag but raw unzipped was loaded
if flags.zipped && ~D.Flags.Interpolated && isfield(D.Data, 'Mask')
    clear D_raw

    % Remove out-of-date fields
    D = rmfield(D.Data, 'Mask');
    if isfield(D.Data, 'WeightedMagnitude'), D = rmfield(D.Data, 'WeightedMagnitude'); end

    % Correct name
    snames.Preprocessed = [snames.Preprocessed 'zip'];
end

%% Preprocessing
ppFLAG = false;

% Pad array
if flags.zipped && ~D.Flags.Interpolated
    % Pad image space data to prevent ghosting
    psize = [D.Size(1), D.Size(2), 0, 0]; % pad size
    img = padarray(D.Data.Image, psize/2, 0, 'both');

    % Convert to kspace
    kSpace = zeros(size(img), 'like', D.Data.Image);
    for j = 1:D.Size(4)
        for i = 1:D.Size(3)
            kSpace(:,:,i,j) = ifftshift(ifft2(ifftshift(img(:,:,i,j))));
        end
    end

    if flags.verbose, fprintf('\nZero-filling data...'); tic; end
    kSpace = padarray(kSpace, psize, 0, 'both');

    % Convert back to image data
    D.Data.Image = zeros(size(kSpace), 'like', kSpace);
    for j = 1:D.Size(4)
        for i = 1:D.Size(3)
            D.Data.Image(:,:,i,j) = fftshift(fft2(fftshift(kSpace(:,:,i,j))));
        end
    end
    clear kSpace

    % Remove excess padding
    D.Data.Image = D.Data.Image(D.Size(1)+1:end-D.Size(1),D.Size(2)+1:end-D.Size(2),:,:);

    % Update data information
    D.VoxelSize(1:2) = D.VoxelSize(1:2).*(D.Size(1:2)./size(D.Data.Image,[1 2]));
    D.Size = size(D.Data.Image);

    % Update flag & plot range
    ppFLAG = true;
    D.Flags.Interpolated = true;
    plrange = {1:D.Size(1), 1:D.Size(2), 5:5:D.Size(3)-5};
    if flags.verbose, tm = toc; fprintf('Done (%0.2f sec)', tm); end
end

% Obtain magnitude weighting
if ~isfield(D.Data, 'WeightedMagnitude')
    if flags.verbose; fprintf('\nCalculating weighted magnitude...'); tic; end
    mag = get_echoMIP(D.Data.Image);
    D.Data.WeightedMagnitude = (mag - min(mag,[],'all'))./(max(mag,[],'all') - min(mag,[],'all')); % [0, 1]
    clear mag
    ppFLAG = true;
    if flags.verbose, tm = toc; fprintf('Done (%0.2f sec)', tm); end
end

% Create mask (from vlGC & MATLAB)
if ~isfield(D.Data, 'Mask')
    D.Data.Mask = Operations.Mask(D.Data.Image, D.Data.WeightedMagnitude, flags.verbose, 3, flags.nobone);

    ppFLAG = true;

    % Check if data needs to be flipped
    if sum(D.Data.Mask(:,:,1),'all') > sum(D.Data.Mask(:,:,end),'all')
        D.Data.Image = flip(D.Data.Image, 3);
        D.Data.Mask = flip(D.Data.Mask, 3);
        D.Data.WeightedMagnitude = flip(D.Data.WeightedMagnitude, 3);
        D.B0Direction(3) = -D.B0Direction(3);
    end
end

% Trim off bad/unecessary slices
if flags.trimmed && isscalar(D.Flags.Trimmed.ThroughPlane) && ~D.Flags.Trimmed.ThroughPlane
    % Combine magnitude & phase
    comb = abs(D.Data.Image(:,:,:,end));
    comb = -pi + ((comb - min(comb,[],'all'))./(max(comb,[],'all') - min(comb,[],'all'))).*(2*pi); % rescale amplitudes to range [-pi, pi]
    comb(:,1:D.Size(2)/2,:) = angle(D.Data.Image(:,1:D.Size(2)/2,:,end));

    % Trim sides of image temporarily
    plrange_tmp = plrange;
    plrange_tmp{2} = plrange_tmp{2}(ceil(D.Size(2)*0.15)+1:end-ceil(D.Size(2)*0.15));

    % Plot data & request input
    % fig = plotmygraph(comb, PlotTitle=['D' num2str(dog) 'D' num2str(date) ' Magnitude'], DataRange=plrange_tmp, IsotropicVoxel=true);
    % 
    % pause(1)
    % figure(fig)
    % trimFLAG = questdlg('Are there any slices that should be removed prior to data analysis?', ...
    %                     'Trim Check', ...
    %                     'Yes', 'No', 'No');
    % close
    
    trimFLAG = 'No';

    if strcmp(trimFLAG, 'Yes')
        plrange_tmp{3} = 1:D.Size(3);

        % Create cell array of indices
        possibleidx = cell(length(plrange_tmp{3}),1);
        for i = 1:length(possibleidx)
            possibleidx{i} = num2str(plrange_tmp{3}(i));
        end

        % Show all slices
        fig = plotmygraph(comb, PlotTitle=['D' num2str(dog) 'D' num2str(date) ' Magnitude'], DataRange=plrange_tmp, IsotropicVoxel=true);
        pause(1)
        figure(fig)

        % Request start index
        startidx = listdlg('PromptString', 'Select the new start index for your dataset', ...
                           'ListString', possibleidx, ...
                           'SelectionMode','single');

        if isempty(startidx)
            close
            error('No start index was selected.');
        end

        % Request end index
        endidx = listdlg('PromptString', 'Select the new end index for your dataset', ...
                         'ListString', possibleidx, ...
                         'SelectionMode','single');

        if isempty(endidx)
            close
            error('No end index was selected.');
        end

        close; clear fig

        % Round to nearest four
        while rem(endidx-startidx+1,4) ~= 0
            df = 4-rem(endidx-startidx+1,4); % Amount to be added
            if startidx > 1
                if startidx > df+1
                    startidx = startidx - df;
                else
                    startidx = 1;
                end
            else
                endidx = endidx - (4-df);
            end
        end

        % Trim data
        D.Data.Image = D.Data.Image(:,:,startidx:endidx,:);
        D.Data.Mask = D.Data.Mask(:,:,startidx:endidx);
        D.Data.WeightedMagnitude = D.Data.WeightedMagnitude(:,:,startidx:endidx);
        D.Size = size(D.Data.Image);

        ppFLAG = true;
        D.Flags.Trimmed.ThroughPlane = true;
        if ~D.Flags.Trimmed.InPlane
            data.TrimmedIndices = zeros([3 2]);
        end
        data.TrimmedIndices(3,:) = [startidx endidx];
    else
        disp('No slices will be removed.');
    end
end

% Trim FOV
if ~D.Flags.Trimmed.InPlane
    if flags.verbose, fprintf('\nTrimming field of view...'); tic; end
    D = Operations.roiIsol(D);
    ppFLAG = true;

    if flags.verbose, tm = toc; fprintf('Done (%0.2f sec)', tm); end
end

% Update changes & save data
if ppFLAG
    plrange = GenPlotRange(D.Size);
    D_raw = Operations.Save(D, saveFLAG, pth_1H, snames.Preprocessed, flags.verbose);
end

% Plot data
if flags.plot
    plotmygraph(D_raw.Data.WeightedMagnitude, 'PlotTitle', ['D' num2str(dog) 'D' num2str(date) ' Magnitude'], 'DataRange', plrange);
    plotmygraph(angle(D_raw.Data.Image), 'PlotTitle', ['D' num2str(dog) 'D' num2str(date) ' Raw Phase'], 'ColorbarTitle', 'Radians', 'DataRange', plrange);
    plotmygraph(D_raw.Data.Mask, 'PlotTitle', ['D' num2str(dog) 'D' num2str(date) ' Mask'], 'DataRange', plrange);
end
clear comb

%% Bipolar correction
% ppFLAG = false;
% 
% % Correct bipolar phase error
% if ~D.Flags.CorrectedBipolarPhase
%     D = D_raw;
% 
%     if ~strcmp(flags.cscorrection.method, 'BipolarIGC')
%         D = BC.perform(D, D_raw, flags);
% 
%         ppFLAG = true;
%     end
% end
% 
% % Calculate noise standard deviation (if not already found)
% if ~isfield(D.Data, 'NoiseSTD')
%     if flags.verbose; fprintf('\nCalculating Noise STD...'); end
%     opts.max_iter = 60;
%     [~, D.Data.NoiseSTD, relres] = Fit_ppm_complex(D.Data.Image, opts); % No need to include unequal spacing, since only first three echoes used to fit. Also, no difference when calculated 3D vs 2D
% 
%     ppFLAG = true;
%     if flags.verbose, tm = toc; fprintf('Done (%0.2f sec)', tm); end
% 
%     if flags.verbose; fprintf('\nRe-calculating binary mask with noise STD...'); end
%     D.Data.Mask = Operations.Mask(D.Data.Image, D.Data.WeightedMagnitude, flags.verbose, 3, flags.nobone, D.Data.NoiseSTD);
%     if flags.verbose, tm = toc; fprintf('Done (%0.2f sec)', tm); end
% end
% 
% % Save data
% if ppFLAG
%     D_bpc = Operations.Save(D, saveFLAG, pth_1H, snames.BipolarCorrected, flags.verbose);
% end
% 
% % Plot data
% if flags.plot && ~strcmp(flags.chemicalshiftcorrection.method, 'BipolarIGC')
%     if strcmp(flags.bipolarcorrection.method, 'hernando')
%         plotmygraph(sqrt((abs(D_bpc.Data.compWq) - abs(D_bpc.Data.magnWq)).^2 + (abs(D_bpc.Data.compFq) - abs(D_bpc.Data.magnFq)).^2).*D_bpc.Data.Mask, 'PlotTitle', ['D' num2str(dog) 'D' num2str(date) ' Bipolar Correction'], 'ColorbarTitle', 'Radians');
%     else
%         plotmygraph(squeeze(angle(D_bpc.Data.Image(:,:,round(D_bpc.Size(3)/2),:)) - angle(D_raw.Data.Image(:,:,round(D_bpc.Size(3)/2),:))), 'PlotTitle', ['D' num2str(dog) 'D' num2str(date) ' Bipolar Correction'], 'ColorbarTitle', 'Radians');
%         if strcmp(flags.bipolarcorrection.method, 'SEPIAslow')
% 
% 
%             plotmygraph(squeeze(D_bpc.Data.FIT3D(:,:,round(D_bpc.Size(3)/2),:)), 'PlotTitle', ['D' num2str(dog) 'D' num2str(date) ' FIT3D'], 'ColorbarTitle', 'Radians');
%         end
%     end
% end

%% Unwrap data
% Only case this is used is no chemical shift correction
% if ~D.Flags.UnwrappedPhase
%     if ~cscFLAG
%         if ~D.Flags.CombinedEchoes
%             D = PU.fitOptimumWeights(D, flags.unwrapping.method, flags.verbose);
%             D.Data = rmfield(D.Data, "Image");
%         else
%             D = PU.perform(D, flags);
%         end
%     % elseif strcmp(flags.bipolarcorrection, 'MEDI')
%     %     D = PU.perform(D, flags);
%     end
%     % Correct 2D phase unwrapping
%     if flags.unwrapping.corrected && D.Flags.UnwrappedPhase
%         for ec = 1:D.Size(4)
%             if flags.verbose; fprintf('\nCorrecting echo %i, slice %i...', ec, 1); end
%             for sl = 1:D.Size(3)
%                 if flags.verbose; fprintf([repmat('\b',[1 numel(num2str(sl-1))]) '\b\b\b%i...'], sl); end
%                 D.Data.UnwrappedPhase(:,:,sl,ec) = CorrectUnwrap(D.Data.UnwrappedPhase(:,:,sl,ec), angle(D.Data.Image(:,:,sl,ec)), D.Data.Mask(:,:,sl), 'DCOnly');
%             end
%             if flags.verbose; fprintf('\n'); end
%         end
% 
%         % Update flag
%         D.Unwrapping.SliceOffsetCorrection = true;
%     end
% 
%     % Save data
%     if D.Flags.UnwrappedPhase, D_unwrap = Operations.Save(D, saveFLAG, pth_1H, snames.UnwrappedPhase, flags.verbose); end
% end
% 
% % Plot data
% if ~cscFLAG
%     if flags.plot
%         plotmygraph(D_unwrap.Data.UnwrappedPhase, 'PlotTitle', ['D' num2str(dog) 'D' num2str(date) ' Unwrapped Phase'], 'ColorbarTitle', 'Radians', 'DataRange', plrange);
%         if isfield(D_unwrap.Data, 'TotalField')
%             plotmygraph(D_unwrap.Data.TotalField, 'PlotTitle', ['D' num2str(dog) 'D' num2str(date) ' Total Field'], 'ColorbarTitle', 'Frequency (Hz)', 'DataRange', plrange);
%         end
%     end
% end

%% CS Correction
ppFLAG = false;

% Check for previous data
if ~D.Flags.CorrectedChemicalShift && cscFLAG
    clear D
    if ~strcmp(flags.cscorrection.method,'BipolarIGC')
        D = D_bpc;
    else
        D = D_raw;
    end
    
    % Perform chemical shift correction
    delete(gcp('nocreate'))
    D = CSC.perform(D, flags);

    % % Calculate fat frequency shift
    % dfat = -3.5*D.F0;
    % 
    % % Calculate effective fat?
    % effect_fat_Hz = dfat + floor(0.5-dfat*D.deltaTE)/D.deltaTE;
    % 
    % % Iterate through slices
    % if flags.verbose
    %     disp('Initialize extended IDEAL fitting');
    % end
    % for sl = 1:D.Size(3)
    %     if flags.verbose
    %         fprintf('Slice #%i\n', sl);
    %     end
    %     % Fit field map to IDEAL for further refinement
    %     [wwater, wfat, wfreq] = fit_IDEAL(D_raw.Image(:,:,sl,:), D.TE, dfat, D.Data.TotalField(:,:,sl),[],100);
    % 
    %     % Check for fat water swap
    %     if sum(abs(wfat(:)).^2) > sum(abs(wwater(:)).^2)
    %         if flags.verbose
    %             disp('Potential water fat swap');
    %         end
    %         D.Data.TotalField(:,:,sl) = D.Data.TotalField(:,:,sl)+effect_fat_Hz;
    %         [D.Data.Water(:,:,sl), D.Data.Fat(:,:,sl), D.Data.TotalField(:,:,sl)] = fit_IDEAL(D_raw.Image(:,:,sl,:), D.TE, dfat, D.Data.TotalField(:,:,sl),[],100);
    %     else
    %         D.Data.Water(:,:,sl) = wwater;
    %         D.Data.Fat(:,:,sl) = wfat;
    %         D.Data.TotalField(:,:,sl) = wfreq;
    %     end
    % end
    % 
    % % Update flag
    % D.CSCorrection.FatWaterSwapChecked.Automatic = true;

    ppFLAG = true;
end

%% Check for fat/water swaps - automatic
% if ~D.CSCorrection.FatWaterSwapChecked.Automatic && ~D.CSCorrection.FatWaterSwapChecked.Manual
%     % Preallocate array
%     if ~isfield(D.CSCorrection, 'SwappedSlice'), D.CSCorrection.SwappedSlice = false(1, D.Size(3)); end
% 
%     % Iterate through slices
%     for sl = 1:D.Size(3)
%         % Check if more fat than water signal (but only if 50% greater, otherwise might just be fatty slice)
%         if sum(abs(D.Data.Fat(:,:,sl)).^2,'all') > 1.5*sum(abs(D.Data.Water(:,:,sl)).^2,'all')
%             if flags.verbose
%                 disp('Potential fat water swap')
%             end
% 
%             % Correct the data
%             truefat = D.Data.Water(:,:,sl);
%             truewater = D.Data.Fat(:,:,sl);
%             D.Data.Water(:,:,sl) = truewater;
%             D.Data.Fat(:,:,sl) = truefat;
% 
%             % Mark slice as swapped
%             D.CSCorrection.SwappedSlice(sl) = ~D.CSCorrection.SwappedSlice(sl);
%         end
%     end
% 
%     % Update flags
%     D.CSCorrection.FatWaterSwapChecked.Automatic = true;
%     ppFLAG = true;
% end

%% Check for fat/water swaps - manual
% if ~D.CSCorrection.FatWaterSwapChecked.Manual
%     % Preallocate array
%     if ~isfield(D.CSCorrection, 'SwappedSlice'), D.CSCorrection.SwappedSlice = false(1, D.Size(3)); end
% 
%     % Update Mask
%     D.Data.Mask = Operations.Mask(D.Data.Water+D.Data.Fat, D.Data.WeightedMagnitude, flags.verbose, 3, flags.nobone);
% 
%     % Calculate fat fraction
%     FF = CSC.PDFF(D.Data.Water, D.Data.Fat, D.Data.Mask);
% 
%     % Iterate through slices
%     Dispnum = 5;
%     for sl = 1:Dispnum:D.Size(3)
%         if sl+Dispnum-1 > D.Size(3)
%             possibleswaps = sl:D.Size(3);
%         else
%             possibleswaps = sl:(sl+Dispnum-1);
%         end
%         dispnum = length(possibleswaps);
%         % Plot data
%         Z = CSC.showSwaps(abs(D_raw.Data.Image(:,:,:,1)), FF, sl, dispnum);
% 
% 
%         % Request input
%         swapFLAG = questdlg('Are there any swapped fat and water maps?', ...
%                             'Swap Check', ...
%                             'Yes', 'No', 'Cancel', 'Cancel');
% 
%         % Correct for swap
%         if strcmp(swapFLAG, 'Yes')
%             % Set while loop flag
%             happy = false;
%             % Request further input
%             while ~happy
%                 possibleidx = cell(dispnum,1);
%                 for i = 1:length(possibleswaps), possibleidx{i} = num2str(possibleswaps(i)); end
%                 swapidx = listdlg('PromptString', 'Select which slices are swapped', ...
%                                   'ListString', possibleidx);
% 
%                 % Ensure input was given
%                 if isempty(swapidx)
%                     warning('No input was given. No changes will be made to the data.')
%                 else
%                     swaps = possibleswaps(swapidx);
% 
%                     % Check if there are multiple components in slice
%                     for sl2 = swaps
%                         CC = bwconncomp(FF(:,:,sl2));
%                         numRegPixels = cellfun(@numel,CC.PixelIdxList); 
%                         nPix = sum(cellfun(@numel,CC.PixelIdxList));
% 
%                         if sum(numRegPixels./sum(numRegPixels)*100 > 0.1) > 1
%                             compswapFLAG = questdlg('Is the whole image swapped, or just a single component?', ...
%                                                     ['Slice #' num2str(sl2) ' Swap Check'], ...
%                                                     'Entire Image', 'Single Component', 'Entire Image');
% 
%                             if strcmp(compswapFLAG, 'Entire Image')
%                                 compswapFLAG = false;
%                             elseif strcmp(compswapFLAG, 'Single Component')
%                                 compswapFLAG = true;
%                             elseif isempty(compswapFLAG)
%                                 warning('No input was given. It will be assumed that the entire slice is swapped.')
%                                 compswapFLAG = false;
%                             end
% 
%                             % Assume the smallest component has been swapped
%                             % Note: this is most likely, but also hard to code for user input on which component
%                             if compswapFLAG
%                                 idx = find(numRegPixels == min(numRegPixels(numRegPixels./sum(numRegPixels)*100 > 0.1)));
% 
%                                 % swap smallest component
%                                 tmpwater = D.Data.Water(:,:,sl2); truewater = D.Data.Water(:,:,sl2);
%                                 tmpfat = D.Data.Fat(:,:,sl2); truefat = D.Data.Fat(:,:,sl2);
% 
%                                 truewater(CC.PixelIdxList{idx}) = tmpfat(CC.PixelIdxList{idx});
%                                 truefat(CC.PixelIdxList{idx}) = tmpwater(CC.PixelIdxList{idx});
%                                 D.Data.Water(:,:,sl2) = truewater;
%                                 D.Data.Fat(:,:,sl2) = truefat;
%                             end
%                         else
%                             compswapFLAG = false;
%                         end
% 
%                         % Correct the data
%                         if ~compswapFLAG
%                             truefat = D.Data.Water(:,:,sl2);
%                             truewater = D.Data.Fat(:,:,sl2);
%                             D.Data.Water(:,:,sl2) = truewater;
%                             D.Data.Fat(:,:,sl2) = truefat;
%                         end
%                     end
% 
%                     FF(:,:,swaps) = CSC.PDFF(D.Data.Water(:,:,swaps), D.Data.Fat(:,:,swaps), D.Data.Mask(:,:,swaps));
% 
%                     % Mark slice as swapped
%                     D.CSCorrection.SwappedSlice(swaps) = ~D.CSCorrection.SwappedSlice(swaps);
%                 end
% 
%                 % Plot data
%                 close
%                 CSC.showSwaps(abs(D_raw.Data.Image), FF, sl, dispnum);
% 
%                 % Ask if any other indexes need to be changed
%                 swapFLAG = questdlg('Do any other slices need to be swapped?', ...
%                                     'Swap Check', ...
%                                     'Yes', 'No', 'Cancel', 'Cancel');
% 
%                 % Check new input
%                 if strcmp(swapFLAG, 'No')
%                     happy = true;
%                 elseif strcmp(swapFLAG, 'Cancel') || isempty(swapFLAG)
%                     close
%                     error('You have ended the current fat-water swap check loop before completion. No data will be saved.')
%                 end
%             end
%             close
%         elseif strcmp(swapFLAG, 'No')
%             close
%         elseif strcmp(swapFLAG, 'Cancel') || isempty(swapFLAG)
%             fprintf('\nYou have ended the fat-water swap checks before completion. No data will be saved.\n')
%             close
%             break
%         end
%     end
% 
%     % Update flags
%     if ~strcmp(swapFLAG, 'Cancel')
%         D.CSCorrection.FatWaterSwapChecked.Manual = true;
%         ppFLAG = true;
%         plotmygraph(FF);
%     else
%         ppFLAG = false;
%     end
% end

%% Save data
if ppFLAG
    % fields2rm = {'Image', 'WeightedMagnitude'};
    % if isfield(D.Data, 'Image'), D.Data = rmfield(D.Data, 'Image'); end
    D_cscorr = Operations.Save(D, saveFLAG, pth_1H, snames.CorrectedCS, flags.verbose);
end

% No need for relaxation & flip angle correction: https://doi.org/10.1002/mrm.21301
%   Since flip angle is small (8 deg) & TR is long (616msec), T1 & flip angle effects can be ignored
%   And since T2* >> TE, T2* effects can be ignored too
% % Correct for relaxation & flip angle
% if ~D.CSCorrection.Corrected
%     % Assumed values for T1
%     % Note: T2* >> TE & thus can be ignored
%     T1w = 1.3;
%     T1f = 0.25;
%     T2sw = 0.025;
%     T2sf = 0.01;
% 
%     % Correct
%     D.Data.Water = (D.Data.Water.*((1 - exp(-D.RepetitionTime/T1w))*exp(-TE/T2sw)*sin(D.FlipAngle*(pi/180))))/(1 - exp(-D.RepetitionTime/T1w)*cos(D.FlipAngle*(pi/180)));
%     D.Data.Fat   =   (D.Data.Fat.*((1 - exp(-D.RepetitionTime/T1f))*exp(-TE/T2sf)*sin(D.FlipAngle*(pi/180))))/(1 - exp(-D.RepetitionTime/T1f)*cos(D.FlipAngle*(pi/180)));
% 
%     % Update flag
%     D.CSCorrection.Corrected = true;
% end

% Plot data
if flags.plot
    plotmygraph(abs(D_cscorr.Data.Water), 'PlotTitle', ['D' num2str(dog) 'D' num2str(date) ' Water Map'], 'DataRange', plrange);
    % plotmygraph(angle(D_cscorr.Data.Water), 'PlotTitle', ['D' num2str(dog) 'D' num2str(date) ' Water Phase'], 'ColorbarTitle', 'Radians', 'DataRange', plrange);
    plotmygraph(abs(D_cscorr.Data.Fat), 'PlotTitle', ['D' num2str(dog) 'D' num2str(date) ' Fat Map'], 'DataRange', plrange);
    % plotmygraph(angle(D_cscorr.Data.Fat), 'PlotTitle', ['D' num2str(dog) 'D' num2str(date) ' Fat Phase'], 'ColorbarTitle', 'Radians', 'DataRange', plrange);
    plotmygraph(CSC.PDFF(D_cscorr.Data.Water, D_cscorr.Data.Fat, D_cscorr.Data.Mask), 'PlotTitle', ['D' num2str(dog) 'D' num2str(date) ' PDFF'], 'ColorbarTitle', 'Fat Fraction', 'DataRange', plrange);
    plotmygraph(D_cscorr.Data.TotalField, 'PlotTitle', ['D' num2str(dog) 'D' num2str(date) ' Total Field'], 'ColorbarTitle', 'Frequency (Hz)', 'DataRange', plrange);
end

%%
pdff = CSC.PDFF(D_cscorr.Data.Water, D_cscorr.Data.Fat, D_cscorr.Data.Mask);
%% Background field removal

% Check for previous data
if ~D.Flags.Removed2DBackgroundField
    % Remove unecessary fields
    if ~cscFLAG
        D = D_unwrap; %#ok<*UNRCH>
    else
        D.Data = rmfield(D_cscorr.Data, {'Water', 'Fat', 'R2StarMap'}); 
    end
    rmfields = {'Error', 'UnwrappedPhase', 'UncorrectedTotalField'};
    for i = length(rmfields)
        if isfield(D.Data, rmfields{i})
            D.Data = rmfield(D.Data, rmfields{i});
        end
    end

    % Remove background field via 2D process
    % if ~cscFLAG % starting with just unwrapped data
    %     for ec = 1:D.Size(4)
    %         D.Data.Mask = repmat(D.Data.Mask, [1 1 1 7]);
    % 
    %         D.Data.TotalField(:,:,:,ec) = D.Data.UnwrappedPhase(:,:,:,ec)./(2*pi*D.TE(ec));
    %         reverseStr = '';
    %         for sl = 1:D.Size(3)
    %             if flags.verbose; reverseStr = UpdatePercent(100*sl/D.Size(3), reverseStr); end
    %             [D.Data.LocalField(:,:,sl,ec), D.Data.Mask(:,:,sl,ec)] = BKGRemovalVSHARP_2D(data.Data.TotalField(:,:,sl,ec), D.Data.Mask(:,:,sl,ec), D.Size(1:2), 'radius', 12:-1:1);
    %         end
    %         fprintf('\tDone!\n');
    % 
    %         [D.Data.LocalField(:,:,:,ec), D.Data.Mask(:,:,:,ec)] = BKGRemovalVSHARP(D.Data.TotalField(:,:,:,ec), D.Data.Mask(:,:,:,ec), D.Size(1:3), 'radius', 12:-1:1);
    %     end
    % 
    %     % Combine echoes
    %     iF1 = abs(D_raw.Image).*exp(1i*D.Data.LocalField);
    %     [D.Data.LocalField, D.Data.NoiseSTD] = Fit_ppm_complex_bipolar(iF1);
    % else
        D = BFR.perform(D, flags.bfremoval.D2, 2, flags.verbose);
    % end

    % Save data
    if ~jointD23FLAG
        D_2Dlocfield = Operations.Save(D, saveFLAG, pth_1H, snames.LocalField.D2, flags.verbose);
    else
        D_3Dlocfield = Operations.Save(D, saveFLAG, pth_1H, snames.LocalField, flags.verbose);
    end
end

% Check for previous data
if ~D.Flags.Removed3DBackgroundField
    % Reset data structure
    D = D_2Dlocfield;

    % Rename field
    D.Data.TotalField = D.Data.LocalField;
    D.Data = rmfield(D.Data, 'LocalField');

    % Remove background field via 3D process
    D = BFR.perform(D, flags.bfremoval.D3, 3, flags.verbose);
    plrange = GenPlotRange(D.Size);

    % Save data
    D_3Dlocfield = Operations.Save(D, saveFLAG, pth_1H, snames.LocalField.D3, flags.verbose);
end

% Plot data
if flags.plot
    plotmygraph(D_3Dlocfield.Data.LocalField, 'PlotTitle', ['D' num2str(dog) 'D' num2str(date) ' 3D Local Field'], 'ColorbarTitle', 'Frequency (Hz)', 'DataRange', plrange);
    if isfield(D_3Dlocfield.Data, 'BackgroundField')
        plotmygraph(D_3Dlocfield.Data.BackgroundField, 'PlotTitle', ['D' num2str(dog) 'D' num2str(date) ' 3D Background Field'], 'ColorbarTitle', 'Frequency (Hz)', 'DataRange', plrange);
    end
    if ~jointD23FLAG
        plotmygraph(D_2Dlocfield.Data.LocalField, 'PlotTitle', ['D' num2str(dog) 'D' num2str(date) ' 2D Local Field'], 'ColorbarTitle', 'Frequency (Hz)', 'DataRange', plrange);
        if isfield(D_2Dlocfield.Data, 'BackgroundField')
            plotmygraph(D_2Dlocfield.Data.BackgroundField, 'PlotTitle', ['D' num2str(dog) 'D' num2str(date) ' 2D Background Field'], 'ColorbarTitle', 'Frequency (Hz)', 'DataRange', plrange);
        end
    end
end

%% Dipole inversion L-curve analyis
if debugFLAG
    switch flags.dipoleinversion.method
        case 'MEDI'
            params = struct('perc', 0.9, ...                                  Default 0.9                   , number of voxels considered 'edge' in L1 regularization
                        'max_iter', 10, ...                                   Default 10                    , number of Gauss-Newton solver iterations
                     'cg_max_iter', 100, ...                                  Default 100                   , number of conjugate gradient iterations within each GN solver iteration
                  'data_weighting', 1, ...                                    Default 1                     , data weighting mode (0=uniform weighting, 1=SNR weighting)
                  'tol_norm_ratio', 0.2, ...                                  Default 0.1                   , threshold value of relative update changes
                         'gpuFLAG', true);  
            params.lambda = 10.^(1.75:0.25:8.25);
        case 'NDI'
            params = struct('padFLAG', true, ...                             Default true                   , padding flag
                      'maxOuterInter', 1000, ...                             Default 1000                   , maximum number of iterations
                             'weight', D_3Dlocfield.WeightedMagnitude, ...   Default mag                    , data consistency/fidelity weighting
                                'tau', 1, ...                                Default 1                      , gradient descent rate
                            'precond', false, ...                            Default false                  , preconditioned solution for stability (start QSM with 3*weight*lfs instead of array of zeros)
                        'isShowIters', false);  
            params.alpha = logspace(-10, -3, 40);
        case 'Magnitude Weighted L1'
            params = struct('padFLAG', true, ...                              Default true                  , padding flag
                          'lambda_L1', logspace(-7, -1, 80), ...              Default logspace(-6, -2, 30)  , maximum number of iterations
                          'lambda_L2', logspace(-7, -1, 80));               % Default logspace(-6, -2, 30)  , maximum number of iterations
        case 'TVDI'
            params = struct('lambda', logspace(-6, -2, 40), ...               Default 5e-4                  , regularization paramter
                          'weight', D.Data.WeightedMagnitude, ...                  Default mag                   , data consistency weighting (mask or magnitude)
                          'iter', 30, ...                                     Default 500                   , number of NLCG iterations
                          'pnorm', 1);                                      % Default 1                     , regularization type
    
    end
    
    [lambda, kappa, cost_data, cost_reg] = DI.lcurve(D_3Dlocfield, flags.dipoleinversion.method, params, plrange);
end

%% Dipole inversion
% Check for previous data
if ~D.Flags.InvertedDipole
    % Reset data structure
    D = D_3Dlocfield;

    % Obtain susceptibility map
    D = DI.perform(D, flags.dipoleinversion);

    % Save data
    rmfields = {'WeightedMagnitude', 'NoiseSTD', 'TotalField', 'LocalField'};
    for i = 1:numel(rmfields)
        if isfield(D.Data, rmfields{i})
            D.Data = rmfield(D.Data, rmfields{i});
        end
    end
    D_sus = Operations.Save(D, saveFLAG, pth_1H, snames.SusceptibilityMap, flags.verbose);
end

% Plot data
if flags.plot
    plotmygraph(real(D_sus.Data.SusceptibilityMap), 'PlotTitle', ['D' num2str(dog) 'D' num2str(date) ' Susceptibility Map'], 'ColorbarTitle', 'Magnetic Susceptibility (ppm)', 'DataRange', plrange);
end

% end
