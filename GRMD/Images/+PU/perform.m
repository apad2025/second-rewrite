% Unwrap phase
function D = perform(D, flags)
    switch flags.unwrapping.method
        case 'GraphCuts'
            % Transfer data to temporary variables
            tmp_UnwrappedPhase = cell(D.Size(4));
            tmp_RawPhase = cell(D.Size(4));
            tmp_WeightedMagnitude = D.Data.WeightedMagnitude;
            for ec = 1:length(D.TE)
                tmp_RawPhase{ec} = angle(D.Data.Image(:,:,:,ec));
            end
            EC = D.Size(4);
            SL = D.Size(3);
            verbose = flags.verbose;
        
            % Create parallel pool
            pool = parpool;
            opts = parforOptions(pool, 'MaxNumWorkers', 4);
    
            % Iterate through echoes
            parfor (ec = 1:EC, opts)
                if verbose; fprintf('\nUnwrapping echo %i...', ec); end
    
                % Transfer echo to temporary variable & preallocate data
                tmp_ph = tmp_RawPhase{ec};
                tmp_unph = zeros(size(tmp_ph));
    
                % Iterate through slices
                for sl = 1:SL
                    if verbose; fprintf('\nUnwrapping slice %i...\n', sl); end
                    tmp_unph(:,:,sl) = unwrapping_gc(tmp_ph(:,:,sl), tmp_WeightedMagnitude(:,:,sl), VoxelSize, 'no', subsample); %#ok<PFBNS>
                end
                tmp_UnwrappedPhase{ec} = tmp_unph;
            end
    
            % Extract from cell array
            D.Data.UnwrappedPhase = zeros(D.Size);
            for ec = 1:EC
                D.Data.UnwrappedPhase(:,:,:,ec) = tmp_UnwrappedPhase{ec};
            end
    
            % Delete parallel pool
            delete(gcp('nocreate'))
    
            % Release memory
            clear EC SL VoxelSize subsample verbose tmp_UnwrappedPhase tmp_RawPhase tmp_WeightedMagnitude pool opts tmp_ph tmp_unph
        
        case 'RegGrow'
            D.Data.UnwrappedPhase = zeros(D.Size);
            for ec = 1:length(D.TE)
                if flags.verbose; fprintf('\nUnwrapping echo %i...', ec); end
                for sl = 1:D.Size(3)
                    if flags.verbose; fprintf('\nUnwrapping slice %i...\n', sl); end
                    D.Data.UnwrappedPhase(:,:,sl,ec) = unwrapPhase(squeeze(abs(D.Data.Image(:,:,sl,ec))), squeeze(angle(D.Data.Image(:,:,sl,ec))), size(D.Data.Image(:,:,sl,ec)));
                end
            end
    end

    % Update flags
    D.Flags.UnwrappedPhase = true;
    D.Unwrapping = struct('Method', flags.unwrapping.method);
    if isfield(flags.unwrapping, 'subsample'), D.Unwrapping.SubsampleFactor = flags.unwrapping.subsample; end
end