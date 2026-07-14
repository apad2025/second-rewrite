% Correct chemical shift
function D = perform(D, flags) 
% Inputs:
%     data: data structure
%    flags: processing structure
% 
% Outputs:
%     data: data structure

    fprintf('\nPerforming chemical shift correction using %s...\n', flags.cscorrection.method)

    % Algorithm-specific parameters
    algoParams = CSC.grabPars(flags);
    if flags.zipped, mnw = 4; else, mnw = 6; end

    % Isolate flags
    verboseFLAG = flags.verbose;
    flags = flags.cscorrection;

    % Format data to correct shape
    images = zeros(D.Size(1), D.Size(2), D.Size(3), 1, D.Size(4));
    images(:,:,:,1,:) = D.Data.Image;

    % Remove last echo, if necessary
    if any(strcmp(flags.method, {'IGC','MixedIGC','BipolarIGC'})) && rem(size(images,5),2)==1
        images = images(:,:,:,:,1:end-1);
        D.TE = D.TE(1:end-1);
        D.deltaTE = mean(diff(D.TE));
        D.Size(4) = size(images,5);
    end

    % Create data structure
    dataParams = struct('FieldStrength', D.B0, ...
                                   'TE', D.TE, ...
                'PrecessionIsClockwise', 1);

    % Run correction algorithm
    switch flags.method
        case 'IGC'
            res = CSC.IGC(images, dataParams, algoParams, verboseFLAG, mnw, false);

            % Combine final data
            D.Data.TotalField = res.FM;
            D.FieldMap = res.FM;
            D.Data.Water = res.W;
            D.Data.Fat = res.F;

        case 'MixedIGC'
            [magn, comp] = CSC.IGC(images, dataParams, algoParams, verboseFLAG, mnw, true);

            % Combine final data
            D.Data.TotalField = comp.FM;
            D.FieldMap = magn.FM;
            D.Data.Water = comp.W;
            D.Data.Fat = comp.F;

        case 'BipolarIGC'
            % Grab extra parameters
            dataParams.FieldStrength = D.B0;
            dataParams.voxelSize = D.VoxelSize;
            dataParams.images = images;
            dataParams.mask_fwseparation = 1;

            outParams = Function_Bipolar_GC(dataParams, algoParams, 12, verboseFLAG);

            % Combine final data
            D.Data.Image = squeeze(outParams.corrected_bipolar_signal); % bipolar signal transformed into unipolar equivalent
            D.Data.Water = outParams.species(2).amps;
            D.Data.Fat = outParams.species(1).amps;
            D.Data.TotalField = outParams.fieldmap;
            D.Data.R2Star = outParams.r2starmap;
            D.Data.Phi = outParams.phi_map; % related to phase modulation due to bipolar readout
            D.Data.Epsilon = outParams.eps_map; % related to amplitude modulation due to bipolar readout
            D.Data.BipolarError = outParams.bipolar_error_map_theta; % phi - i*eps;
            D.Data.Correction = outParams.total_correction; % correction to remove bipolar induced effects, e^(i*BipolarError)

        case 'vlGC'
            % Check parameters
            algoParams = checkParamsAndSetDefaults_GANDALF(dataParams, algoParams, verboseFLAG);
    
            % Calculate key data in input images
            dataParams.images = images;
            VARPROparams = calculateMinimaDirect(dataParams, algoParams, verboseFLAG);

            % Delete parallel pool
            delete(gcp('nocreate'))

            if verboseFLAG; fprintf('\nCorrecting slices...'); end
            outWaterFatParams = GANDALF(dataParams, algoParams, VARPROparams, verboseFLAG);
            D.Data.Water = outWaterFatParams.water;
            D.Data.Fat = outWaterFatParams.fat;
            D.Data.TotalField = outWaterFatParams.fieldmap;
            D.Data.R2Star = outWaterFatParams.r2starmap;

        case 'IDEAL-CE'
            % Data preallocation
            D.Data.Water = zeros(D.Size(1), D.Size(2), D.Size(3));
            D.Data.Fat = D.Data.Water;
            D.Data.TotalField = D.Data.Water;
            D.Data.UnwrappedPhase = D.Data.Water;
            D.Data.Error = D.Data.Water;

            % Iterate through slices
            for sl = 1:D.Size(3)
                if verboseFLAG; fprintf('\nCorrecting slice %i...\n', sl); end

                % Run algorithm
                [params, soserror] = presco(D.TE, squeeze(images(:,:,sl,:,:)), D.B0, 'mask', D.Data.Mask(:,:,sl), 'species', species, 'kspace_shift', algoParams.ShiftedkSpace, 'muB', algoParams.muB, 'muR',  algoParams.muR, 'smooth_phase', algoParams.SmoothedPhase, 'smooth_field', algoParams.SmoothedField,  'filter', algoParams.Filter,  'noise', algoParams.NoiseSTD, 'maxit', algoParams.MaxIterations, 'display', false);
                
                % Combine into main structure
                D.Data.TotalField(:,:,sl) = params.B0; 
                D.Data.Water(:,:,sl) = params.W; 
                D.Data.Fat(:,:,sl) = params.F;
                D.Data.UnwrappedPhase(:,:,sl) = params.PH;
                D.Data.R2Star(:,:,sl) = params.R2;
                D.Data.Error(:,:,sl) = soserror;
            end

        case 'SPURS'
            % output: 
            %     - wwater: the water map
            %     - wfat:   the fat map
            %     - wfreq:  the field map in rad after running IDEAL as fine tunning,input for QSM
            %     - wunwph_uf:  the field map after unwrapping and unfat,initial guess for IDEAL
            %     - unwphw: phase unwrapping result
            %
            % input:
            %      - iField : a multi-echo 4 dimentional data (Note that if the water
            %      fat map totally swap, try conj(iField) instead of iField as input)
            %      - how to choose iField or conj(iField) as input:
            %             if PrecessionIsClockwise = 1, [] = spurs_gc(conj(iField),TE,CF,voxel_size);
            %             if PrecessionIsClockwise = -1, [] = purcs_gc(iField,TE,CF,voxel_size);

            % Data preallocation
            D.Data.Water = zeros(D.Size(1), D.Size(2), D.Size(3));
            D.Data.Fat = D.Data.Water;
            D.Data.TotalField = D.Data.Water;
            D.Data.UnwrappedPhase = D.Data.Water;
            D.Data.UncorrectedTotalField = D.Data.Water;

            % Iterate through slices
            for sl = 1:D.Size(3)
                if verboseFLAG; fprintf('\nCorrecting and unwrapping slice %i...\n', sl); end
                
                % Run algorithm
                [D.Data.Water(:,:,sl), D.Data.Fat(:,:,sl), D.Data.TotalField(:,:,sl), D.Data.UncorrectedTotalField(:,:,sl), D.Data.UnwrappedPhase(:,:,sl), ~] = spurs_gc(D.Data.Image(:,:,sl,:), D.TE, D.F0*1e6, D.VoxelSize, D.Flags.CorrectedBipolarPhase, flags.subsample);
            end

        case 'Hierarchical IDEAL'
            % Data preallocation
            D.Data.Water = zeros(D.Size(1), D.Size(2), D.Size(3));
            D.Data.Fat = D.Data.Water;
            D.Data.TotalField = D.Data.Water;
            D.Data.R2Star = D.Data.Water;
            D.Residual = D.Data.Water;

            for sl = 1:D.Size(3)
                if verboseFLAG; fprintf('\nCorrecting slice %i...\n', sl); end

                % Isolate slice data
                dataParams.images = images(:,:,sl,:,:);

                % Run algorithm
                outParams = fw_i2cm0c_3pluspoint_tsaojiang(dataParams, algoParams);

                % Combine into main structure
                D.Data.Water(:,:,sl) = outParams.species(1).amps;
                D.Data.Fat(:,:,sl) = outParams.species(2).amps;
                D.Residual(:,:,sl) = outParams.fiterror;
                D.Data.TotalField(:,:,sl) = outParams.phasemap;
                D.Data.R2Star(:,:,sl) = outParams.r2starmap;
            end

            % Unwrap field map
            if verboseFLAG
                verboselevel = 'yes';
            else
                verboselevel = 'no';
            end
            for sl = 1:D.Size(3)
                D.Data.TotalField(:,:,sl) = unwrapping_gc(angle(D.Data.TotalField(:,:,sl)), D.Data.WeightedMagnitude(:,:,sl), D.VoxelSize, verboselevel, 1);
            end

            % Correct according to tsao jiang code
            D.Data.TotalField = (D.Data.TotalField/(2*pi) - 1/3)/(42.5774780505984*D.deltaTE);

        case 'Golden Section'
            % Data preallocation
            D.Data.Water = zeros(D.Size(1), D.Size(2), D.Size(3));
            D.Data.Fat = zeros(D.Size(1), D.Size(2), D.Size(3));
            D.Data.TotalField = zeros(D.Size(1), D.Size(2), D.Size(3));

            for sl = 1:D.Size(3)
                if verboseFLAG; fprintf('\nCorrecting slice %i...\n', sl); end

                % Isolate slice data
                dataParams.images = images(:,:,sl,:,:);

                % Run algorithm
                outParams = fw_3point_wm_goldSect(dataParams, algoParams);

                % Combine into main structure
                D.Data.Water(:,:,sl) = outParams.species(1).amps;
                D.Data.Fat(:,:,sl) = outParams.species(2).amps;
                D.Data.TotalField(:,:,sl) = outParams.fieldmap;
            end
    end

    % Add parameters to output structure
    switch flags.method
        case {'IGC','MixedIGC'}
            D.CSCorrection = struct('Method', flags.method, ...
                              'SubsampleFactor', algoParams.SUBSAMPLE, ...
                                      'Species', struct('Name', {algoParams.species(1).name, algoParams.species(2).name}, ...
                                                   'Frequency', {algoParams.species(1).frequency, algoParams.species(2).frequency}, ...
                                           'RelativeAmplitude', {algoParams.species(1).relAmps, algoParams.species(2).relAmps}), ...
                                   'CliqueSize', algoParams.size_clique, ...
                                'FieldMapRange', algoParams.range_fm, ...
                              'FieldMapNumbers', algoParams.NUM_FMS, ...
                           'GraphCutIterations', algoParams.NUM_ITERS, ...
                  'OptimizationTransferDescent', algoParams.DO_OT, ...
                                       'Lambda', algoParams.lambda, ...
                               'LambdaMapPower', algoParams.LMAP_POWER, ...
                     'LambdaMapSmoothingFactor', algoParams.LMAP_EXTRA);

            if strcmp(flags.method, 'IGC-mixed')
                D.CSCorrection.Threshold = algoParams.THRESHOLD;
                D.CSCorrection.PhaseCorruptedEchos = algoParams.NUM_MAGN;
                D.CSCorrection.R2StarRange = algoParams.range_r2star;
            end

        case 'vlGC'
            D.CSCorrection = struct('Method', flags.method, ...
                              'SubsampleFactor', algoParams.sampling_stepsize, ...
                                      'Species', struct('Name', {algoParams.species(1).name, algoParams.species(2).name}, ...
                                                   'Frequency', {algoParams.species(1).frequency, algoParams.species(2).frequency}, ...
                                           'RelativeAmplitude', {algoParams.species(1).relAmps, algoParams.species(2).relAmps}), ...
                                'FieldMapRange', algoParams.range_fm, ...
                              'FieldMapNumbers', algoParams.NUM_FMS, ...
                              'SamplingPeriods', algoParams.nSamplingPeriods, ...
                                    'Threshold', algoParams.airSignalThreshold_percent/100, ...
                                     'TEPeriod', algoParams.period, ...
                                  'GridSpacing', algoParams.gridspacing);

        case 'IDEAL-CE'
            D.CSCorrection = struct('Method', flags.method, ...
                                      'Species', struct('Name', {species(1).name, species(2).name}, ...
                                                   'Frequency', {species(1).frequency, species(2).frequency}, ...
                                           'RelativeAmplitude', {species(1).relAmps, species(2).relAmps}), ...
                                'ShiftedkSpace', algoParams.ShiftedkSpace,...
                                          'muB', algoParams.muB, ...
                                          'muR',  algoParams.muR, ...
                                'SmoothedPhase', algoParams.SmoothedPhase, ...
                                'SmoothedField', algoParams.SmoothedField, ...
                                       'Filter', algoParams.Filter, ...
                                     'NoiseSTD', algoParams.NoiseSTD, ...
                                'MaxIterations', algoParams.MaxIterations);

        case 'SPURS'
            D.CSCorrection = struct('Method', flags.method);

        case 'Hierarchical IDEAL'
            D.CSCorrection = struct('Method', flags.method, ...
                                      'Species', struct('Name', {algoParams.species(1).name, algoParams.species(2).name}, ...
                                                   'Frequency', {algoParams.species(1).frequency, algoParams.species(2).frequency}, ...
                                           'RelativeAmplitude', {algoParams.species(1).relAmps, algoParams.species(2).relAmps}), ...
                              'MinFractionSize', algoParams.MinFractSizeToDivide, ...
                  'MaxHierarchicalSubdivisions', algoParams.MaxNumDiv, ...
                                       'MinSNR', algoParams.SnrToAssumeSinglePeak, ...
                        'T2starCorrectedFWMaps', algoParams.CorrectAmpForT2star, ...
                                    'MaxR2star', algoParams.MaxR2star); 
        case 'Golden Section'
            D.CSCorrection = struct('Method', flags.method, ...
                                      'Species', struct('Name', {algoParams.species(1).name, algoParams.species(2).name}, ...
                                                   'Frequency', {algoParams.species(1).frequency, algoParams.species(2).frequency}, ...
                                           'RelativeAmplitude', {algoParams.species(1).relAmps, algoParams.species(2).relAmps}));
    end

    % Update data size
    D.Size = size(D.Data.TotalField);

    % Update flags
    D.Flags.CorrectedChemicalShift = true;
    D.Flags.UnwrappedPhase = true;
    D.CSCorrection.FatWaterSwapChecked.Automatic = false;
    D.CSCorrection.FatWaterSwapChecked.Manual = false;
end