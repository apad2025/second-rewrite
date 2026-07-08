% Correct chemical shift
function data = perform(data, flags) 
% Inputs:
%     data: data structure
%    flags: processing structure
% 
% Outputs:
%     data: data structure

    % Algorithm-specific parameters
    algoParams = CSC.grabPars(flags);
    if flags.zipped, mnw = 4; else, mnw = 6; end

    % Isolate flags
    verboseFLAG = flags.verbose;
    flags = flags.cscorrection;

    % Format data to correct shape
    images = zeros(data.Size(1), data.Size(2), data.Size(3), 1, data.Size(4));
    images(:,:,:,1,:) = data.Data.Image;

    % Remove last echo, if necessary
    if any(strcmp(flags.method, {'IGC','MixedIGC','BipolarIGC'})) && rem(size(images,5),2)==1
        images = images(:,:,:,:,1:end-1);
        data.TE = data.TE(1:end-1);
        data.deltaTE = mean(diff(data.TE));
        data.Size(4) = size(images,5);
    end

    % Create data structure
    dataParams = struct('FieldStrength', data.B0, ...
                                   'TE', data.TE, ...
                'PrecessionIsClockwise', 1);

    % Preallocate data
    data.Data.Water = zeros(data.Size(1), data.Size(2), data.Size(3));
    data.Data.Fat = zeros(data.Size(1), data.Size(2), data.Size(3));
    data.Data.TotalField = zeros(data.Size(1), data.Size(2), data.Size(3));

    % Run correction algorithm
    switch flags.method
        case 'IGC'
            res = CSC.IGC(images, dataParams, algoParams, verboseFLAG, mnw, false);

            % Combine final data
            data.Data.TotalField = res.FM;
            data.FieldMap = res.FM;
            data.Data.Water = res.W;
            data.Data.Fat = res.F;

        case 'MixedIGC'
            [magn, comp] = CSC.IGC(images, dataParams, algoParams, verboseFLAG, mnw, true);

            % Combine final data
            data.Data.TotalField = comp.FM;
            data.FieldMap = magn.FM;
            data.Data.Water = comp.W;
            data.Data.Fat = comp.F;

        case 'BipolarIGC'
            dataParams.FieldStrength = data.B0;
            dataParams.voxelSize = data.VoxelSize;
            dataParams.images = images;
            dataParams.mask_fwseparation = data.Data.Mask;

            outParams = Function_Bipolar_GC(dataParams, algoParams, 25, verboseFLAG);
            data.Data = outParams;

            % data.Data.Water = outParams.DualGC.species(2).amps;
            % data.Data.Fat = outParams.DualGC.species(1).amps;
            % data.Data.TotalField = outParams.FieldMap_DualGC;
            % data.Data.R2StarMap = outParams.R2_DualGC;
            % data.Data.Odd.Water = outParams.Water_GC_odd;
            % data.Data.Odd.Fat = outParams.Fat_GC_odd;
            % data.Data.Even.Water = outParams.Water_GC_even;
            % data.Data.Even.Fat = outParams.Fat_GC_even;

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
            data.Data.Water = outWaterFatParams.water;
            data.Data.Fat = outWaterFatParams.fat;
            data.Data.TotalField = outWaterFatParams.fieldmap;
            data.Data.R2StarMap = outWaterFatParams.r2starmap;

        case 'IDEAL-CE'
            % Additional data preallocation
            data.Data.UnwrappedPhase = zeros(data.Size(1), data.Size(2), data.Size(3));
            data.Data.Error = zeros(data.Size(1), data.Size(2), data.Size(3));

            % Iterate through slices
            for sl = 1:data.Size(3)
                if verboseFLAG; fprintf('\nCorrecting slice %i...\n', sl); end

                % Run algorithm
                [params, soserror] = presco(data.TE, squeeze(images(:,:,sl,:,:)), data.B0, 'mask', data.Data.Mask(:,:,sl), 'species', species, 'kspace_shift', algoParams.ShiftedkSpace, 'muB', algoParams.muB, 'muR',  algoParams.muR, 'smooth_phase', algoParams.SmoothedPhase, 'smooth_field', algoParams.SmoothedField,  'filter', algoParams.Filter,  'noise', algoParams.NoiseSTD, 'maxit', algoParams.MaxIterations, 'display', false);
                
                % Combine into main structure
                data.Data.TotalField(:,:,sl) = params.B0; 
                data.Data.Water(:,:,sl) = params.W; 
                data.Data.Fat(:,:,sl) = params.F;
                data.Data.UnwrappedPhase(:,:,sl) = params.PH;
                data.Data.R2StarMap(:,:,sl) = params.R2;
                data.Data.Error(:,:,sl) = soserror;
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

            % Additional data preallocation
            data.Data.UnwrappedPhase = zeros(data.Size(1), data.Size(2), data.Size(3));
            data.Data.UncorrectedTotalField = zeros(data.Size(1), data.Size(2), data.Size(3));

            % Iterate through slices
            for sl = 1:data.Size(3)
                if verboseFLAG; fprintf('\nCorrecting and unwrapping slice %i...\n', sl); end
                
                % Run algorithm
                [data.Data.Water(:,:,sl), data.Data.Fat(:,:,sl), data.Data.TotalField(:,:,sl), data.Data.UncorrectedTotalField(:,:,sl), data.Data.UnwrappedPhase(:,:,sl), ~] = spurs_gc(data.Data.Image(:,:,sl,:), data.TE, data.F0*1e6, data.VoxelSize, data.Flags.CorrectedBipolarPhase, flags.subsample);
            end

        case 'Hierarchical IDEAL'
            % Additional data preallocation
            data.Residual = zeros(data.Size(1), data.Size(2), data.Size(3));

            for sl = 1:data.Size(3)
                if verboseFLAG; fprintf('\nCorrecting slice %i...\n', sl); end

                % Isolate slice data
                dataParams.images = images(:,:,sl,:,:);

                % Run algorithm
                outParams = fw_i2cm0c_3pluspoint_tsaojiang(dataParams, algoParams);

                % Combine into main structure
                data.Data.Water(:,:,sl) = outParams.species(1).amps;
                data.Data.Fat(:,:,sl) = outParams.species(2).amps;
                data.Residual(:,:,sl) = outParams.fiterror;
                data.Data.TotalField(:,:,sl) = outParams.phasemap;
                data.Data.R2StarMap(:,:,sl) = outParams.r2starmap;
            end

            % Unwrap field map
            if verboseFLAG
                verboselevel = 'yes';
            else
                verboselevel = 'no';
            end
            for sl = 1:data.Size(3)
                data.Data.TotalField(:,:,sl) = unwrapping_gc(angle(data.Data.TotalField(:,:,sl)), data.Data.WeightedMagnitude(:,:,sl), data.VoxelSize, verboselevel, 1);
            end

            % Correct according to tsao jiang code
            data.Data.TotalField = (data.Data.TotalField/(2*pi) - 1/3)/(42.5774780505984*data.deltaTE);

        case 'Golden Section'
            for sl = 1:data.Size(3)
                if verboseFLAG; fprintf('\nCorrecting slice %i...\n', sl); end

                % Isolate slice data
                dataParams.images = images(:,:,sl,:,:);

                % Run algorithm
                outParams = fw_3point_wm_goldSect(dataParams, algoParams);

                % Combine into main structure
                data.Data.Water(:,:,sl) = outParams.species(1).amps;
                data.Data.Fat(:,:,sl) = outParams.species(2).amps;
                data.Data.TotalField(:,:,sl) = outParams.fieldmap;
            end
    end

    % Add parameters to output structure
    switch flags.method
        case {'IGC','MixedIGC'}
            data.CSCorrection = struct('Method', flags.method, ...
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
                data.CSCorrection.Threshold = algoParams.THRESHOLD;
                data.CSCorrection.PhaseCorruptedEchos = algoParams.NUM_MAGN;
                data.CSCorrection.R2StarRange = algoParams.range_r2star;
            end

        case 'vlGC'
            data.CSCorrection = struct('Method', flags.method, ...
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
            data.CSCorrection = struct('Method', flags.method, ...
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
            data.CSCorrection = struct('Method', flags.method);

        case 'Hierarchical IDEAL'
            data.CSCorrection = struct('Method', flags.method, ...
                                      'Species', struct('Name', {algoParams.species(1).name, algoParams.species(2).name}, ...
                                                   'Frequency', {algoParams.species(1).frequency, algoParams.species(2).frequency}, ...
                                           'RelativeAmplitude', {algoParams.species(1).relAmps, algoParams.species(2).relAmps}), ...
                              'MinFractionSize', algoParams.MinFractSizeToDivide, ...
                  'MaxHierarchicalSubdivisions', algoParams.MaxNumDiv, ...
                                       'MinSNR', algoParams.SnrToAssumeSinglePeak, ...
                        'T2starCorrectedFWMaps', algoParams.CorrectAmpForT2star, ...
                                    'MaxR2star', algoParams.MaxR2star); 
        case 'Golden Section'
            data.CSCorrection = struct('Method', flags.method, ...
                                      'Species', struct('Name', {algoParams.species(1).name, algoParams.species(2).name}, ...
                                                   'Frequency', {algoParams.species(1).frequency, algoParams.species(2).frequency}, ...
                                           'RelativeAmplitude', {algoParams.species(1).relAmps, algoParams.species(2).relAmps}));
    end

    % Update data size
    % data.Size = size(data.Data.TotalField);

    % Update flags
    data.Flags.CorrectedChemicalShift = true;
    data.Flags.UnwrappedPhase = true;
    data.CSCorrection.FatWaterSwapChecked.Automatic = false;
    data.CSCorrection.FatWaterSwapChecked.Manual = false;
end