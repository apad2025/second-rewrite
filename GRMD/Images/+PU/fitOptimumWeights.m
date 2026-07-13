% Optimum weights fitting
 function D = fitOptimumWeights(D, method, verboseFLAG)

    % Preallocate data
    D.Data.UnwrappedPhase = zeros(D.Size(1), D.Size(2), D.Size(3), D.Size(4)-1);
    D.Data.TotalField = zeros(D.Size(1), D.Size(2), D.Size(3));

    % Update all time variables to second (if necessary)
    if D.TE > 1
        D.TR = D.TR/1000;
        D.TE = D.TE./1000;
    end
    
    % Calculate dims
    dims = size(D.Data.Image(:,:,1,:));
    dims(4) = dims(4) - 1;

    % Iterate through slices
    for sl = 1:D.Size(3)
        if verboseFLAG; fprintf('slice %i...\n', sl); end

        % Isolate slice & mask
        ph = angle(D.Data.Image(:,:,sl,:));
        msk = D.Data.Mask(:,:,sl);

        %%%%%%%%%%%%%%%%%%%%%%%% Step 1: unwrap echo phase %%%%%%%%%%%%%%%%%%%%%%%%
        % compute wrapped phase shift between successive echoes & unwrap each echo phase shift
        UnwrappedPhase = zeros(dims);
        for ec = 1:dims(4)
            fprintf('Unwrapping #%i echo shift...\n',ec);
            tmp	= angle(exp(1i*ph(:,:,:,ec+1))./exp(1i*ph(:,:,:,ec)));
            switch method
                case 'GraphCuts'
                    tmp = unwrapping_gc(tmp, D.Data.WeightedMagnitude(:,:,sl), D.VoxelSize, 'no', 1);
                case 'RegGrow'
                    tmp = unwrapPhase(squeeze(abs(D.Data.Image(:,:,sl,ec))), squeeze(tmp), size(tmp));
            end
            UnwrappedPhase(:,:,:,ec) = tmp-round(mean(tmp(msk == 1))/(2*pi))*2*pi;
        end
        % get phase accumulation over all echoes
        UnwrappedPhase = cumsum(UnwrappedPhase,4);
    
        % get the unwrapped phase accumulation across echoes
        % unwrap first echo
        switch method
            case 'GraphCuts'
                tmp = unwrapping_gc(ph(:,:,:,1), D.Data.WeightedMagnitude(:,:,sl), D.VoxelSize, 'yes', 1);
            case 'RegGrow'
                tmp = unwrapPhase(squeeze(abs(D.Data.Image(:,:,sl,1))), squeeze(ph(:,:,:,1)), size(ph(:,:,:,1)));
        end
        tmp = tmp-round(mean(tmp(msk == 1))/(2*pi))*2*pi;
        ph = cat(4,tmp,UnwrappedPhase + tmp);
    
        %%%%%%%%%%%%%%%%%%%%%%%% Step 2: Compute weights %%%%%%%%%%%%%%%%%%%%%%%%
        % Robinson et al. 2017 NMR Biomed Appendix A2
        [weight, N_std] = compute_optimum_weighting_combining_phase_difference(abs(D.Data.Image(:,:,sl,:)), D.TE);
    
        % standard deviation of field map from weighted avearging
        N_std                 = sqrt(sum(weight.^2 .* N_std.^2,4));    % sqrt(Weighted variance) = SD
        N_std(isnan(N_std))	  = 0;
        N_std(isinf(N_std))	  = 0;
        D.Data.NoiseSTD(:,:,sl) = N_std./norm(N_std(msk~=0));
    
        %%%%%%%%%%%%%%%%%%%%%%%% Step 3: Weighted average %%%%%%%%%%%%%%%%%%%%%%%% units of radHz??
        % 20210803: use for-loop to reduce memory
        TotalField = zeros(dims(1:3), 'like', ph);
        for ec = 1:dims(4)
            % Weighted average of unwrapped phase shift
            TotalField = TotalField + weight(:,:,:,ec).*((ph(:,:,:,ec+1) - ph(:,:,:,1))/(D.TE(ec+1)-D.TE(1)));
        end

        % Combine into structure
        D.Data.UnwrappedPhase(:,:,sl,:) = UnwrappedPhase;
        D.Data.TotalField(:,:,sl) = TotalField;
    end
    
    % Correct output
    D.Data.TotalField(isnan(D.Data.TotalField))=0;
    D.Data.TotalField(isinf(D.Data.TotalField))=0;
    
    % Convert to correct units
    D.Data.TotalField = D.Data.TotalField./(2*pi);
    
    % Update flag
    D.Flags.UnwrappedPhase = true;
    D.Unwrapping = struct('Method', method);
    if strcmp(method, 'GraphCuts'), D.Unwrapping.SubsampleFactor = 1; end

    % Update data size
    D.Size = size(D.Data.UnwrappedPhase);
end