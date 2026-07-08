% Optimum weights fitting
 function data = fitOptimumWeights(data, method, verboseFLAG)

    % Preallocate data
    data.Data.UnwrappedPhase = zeros(data.Size(1), data.Size(2), data.Size(3), data.Size(4)-1);
    data.Data.TotalField = zeros(data.Size(1), data.Size(2), data.Size(3));

    % Update all time variables to second (if necessary)
    if data.TE > 1
        data.TR = data.TR/1000;
        data.TE = data.TE./1000;
    end
    
    % Calculate dims
    dims = size(data.Data.Image(:,:,1,:));
    dims(4) = dims(4) - 1;

    % Iterate through slices
    for sl = 1:data.Size(3)
        if verboseFLAG; fprintf('slice %i...\n', sl); end

        % Isolate slice & mask
        ph = angle(data.Data.Image(:,:,sl,:));
        msk = data.Data.Mask(:,:,sl);

        %%%%%%%%%%%%%%%%%%%%%%%% Step 1: unwrap echo phase %%%%%%%%%%%%%%%%%%%%%%%%
        % compute wrapped phase shift between successive echoes & unwrap each echo phase shift
        UnwrappedPhase = zeros(dims);
        for ec = 1:dims(4)
            fprintf('Unwrapping #%i echo shift...\n',ec);
            tmp	= angle(exp(1i*ph(:,:,:,ec+1))./exp(1i*ph(:,:,:,ec)));
            switch method
                case 'GraphCuts'
                    tmp = unwrapping_gc(tmp, data.Data.WeightedMagnitude(:,:,sl), data.VoxelSize, 'no', 1);
                case 'RegGrow'
                    tmp = unwrapPhase(squeeze(abs(data.Data.Image(:,:,sl,ec))), squeeze(tmp), size(tmp));
            end
            UnwrappedPhase(:,:,:,ec) = tmp-round(mean(tmp(msk == 1))/(2*pi))*2*pi;
        end
        % get phase accumulation over all echoes
        UnwrappedPhase = cumsum(UnwrappedPhase,4);
    
        % get the unwrapped phase accumulation across echoes
        % unwrap first echo
        switch method
            case 'GraphCuts'
                tmp = unwrapping_gc(ph(:,:,:,1), data.Data.WeightedMagnitude(:,:,sl), data.VoxelSize, 'yes', 1);
            case 'RegGrow'
                tmp = unwrapPhase(squeeze(abs(data.Data.Image(:,:,sl,1))), squeeze(ph(:,:,:,1)), size(ph(:,:,:,1)));
        end
        tmp = tmp-round(mean(tmp(msk == 1))/(2*pi))*2*pi;
        ph = cat(4,tmp,UnwrappedPhase + tmp);
    
        %%%%%%%%%%%%%%%%%%%%%%%% Step 2: Compute weights %%%%%%%%%%%%%%%%%%%%%%%%
        % Robinson et al. 2017 NMR Biomed Appendix A2
        [weight, N_std] = compute_optimum_weighting_combining_phase_difference(abs(data.Data.Image(:,:,sl,:)), data.TE);
    
        % standard deviation of field map from weighted avearging
        N_std                 = sqrt(sum(weight.^2 .* N_std.^2,4));    % sqrt(Weighted variance) = SD
        N_std(isnan(N_std))	  = 0;
        N_std(isinf(N_std))	  = 0;
        data.Data.NoiseSTD(:,:,sl) = N_std./norm(N_std(msk~=0));
    
        %%%%%%%%%%%%%%%%%%%%%%%% Step 3: Weighted average %%%%%%%%%%%%%%%%%%%%%%%% units of radHz??
        % 20210803: use for-loop to reduce memory
        TotalField = zeros(dims(1:3), 'like', ph);
        for ec = 1:dims(4)
            % Weighted average of unwrapped phase shift
            TotalField = TotalField + weight(:,:,:,ec).*((ph(:,:,:,ec+1) - ph(:,:,:,1))/(data.TE(ec+1)-data.TE(1)));
        end

        % Combine into structure
        data.Data.UnwrappedPhase(:,:,sl,:) = UnwrappedPhase;
        data.Data.TotalField(:,:,sl) = TotalField;
    end
    
    % Correct output
    data.Data.TotalField(isnan(data.Data.TotalField))=0;
    data.Data.TotalField(isinf(data.Data.TotalField))=0;
    
    % Convert to correct units
    data.Data.TotalField = data.Data.TotalField./(2*pi);
    
    % Update flag
    data.Flags.UnwrappedPhase = true;
    data.Unwrapping = struct('Method', method);
    if strcmp(method, 'GraphCuts'), data.Unwrapping.SubsampleFactor = 1; end

    % Update data size
    data.Size = size(data.Data.UnwrappedPhase);
end