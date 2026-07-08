% DOGEXPLORER explore, load, & preprocess dog data
%   DOGEXPLORER displays the currently available database, including MIDs
%   and GREs, across subjects and imaging sessions.
%
%   path = DOGEXPLORER(subject, session, MID) lists the folder contents 
%   associated with the specified subject & session. 'subject' can be the 
%   name in string format, or a numeric index corresponding to a subject 
%   in the DD database. 'session' can be the date in "YYMMDD"/"YYYYMMDD" 
%   format, or a numeric index corresponding to session in the DD database.
%   'MID' must be an integer corresponding to the measurement ID. The 
%   folder path will be provided as output.
%
%   path = DOGEXPLORER(subject, session, 'DataType', dtype) 
%   lists the data associated with the specified subject, session, and 
%   datatype. 'DataType' can be any of the following (case insensitive):
%       - 1H GRE
%           - 'gradientecho'
%           - 'gre'
%       - 1H localizer
%           - 'localizer'
%           - 'loc'
%       - 1H MRS w/o water suppression
%           - 'hydrogennonsupressed'
%           - 'hydrogenunsupressed'
%           - 'hydrogennonsup'
%           - 'hydrogenunsup'
%           - 'hydrogen'
%           - '1hnonsupressed'
%           - '1hunsupressed'
%           - '1hnonsup'
%           - '1hunsup'
%           - '1h'
%       - 1H MRS w/ water suppression
%           - 'hydrogenwatersupressed'
%           - 'hydrogensupressed'
%           - 'hydrogenwatsup'
%           - 'hydrogensup'
%           - '1hwatersupressed'
%           - '1hsupressed'
%           - '1hwatsup'
%           - '1hsup'
%       - 1H MRS noise scan
%           - 'hydrogennoise'
%           - '1hnoise'
%       - 31P MRS
%           - 'Phosphorus'
%           - '31P'
%
%   [path, varargout] = DOGEXPLORER(subject, session, MID, 'Load') 
%   imports all data related to the subject. If no datatype is provided,
%   then all data will be provided.
%
%   [path, raw_struct, proc_struct, varargout] = DOGEXPLORER(subject, session, MID, 'Preprocess') 
%   imports data and preprocesses the data. If the data is MRS, the 
%   'preprocFIDA' function will be ran. If the data is an image, then only
%   loadima or loaddat will be ran. The function will return the input data
%   path, the raw data structure, the preprocessed data structure, and the 
%   preprocessed water twix file (if applicable).
%
%   [path, fitresults] = DOGEXPLORER(subject, session, MID, 'Fit') 
%   imports preprocessed data and fits it using OXSA. If the data is not
%   MRS, the function will return an error. The function will return the
%   preprocessed data path, and the fitting results structure. MID can be
%   replaced with 'DataType' input.
%
%   [path, raw_struct, proc_struct, varargout] = DOGEXPLORER(subject, session, MID, 'Preprocess', 'Parameters', pars) 
%   imports data and preprocesses the data with the parameters indicated by 
%   the pars structure.
%
%   [path, raw_struct, proc_struct, varargout] = DOGEXPLORER(subject, session, MID, 'Preprocess', 'Override', ovr) 
%   imports data and preprocesses the data. If data has already been
%   preprocessed, the logical input indicated by 'Override' will determine
%   whether or not it should be overriden. The function will return the
%   input data path, the raw data structure, the preprocessed data 
%   structure, and the preprocessed water twix file (if applicable).
%
%   [path, fitresults] = DOGEXPLORER(session, 'Analyze', 'DataType', dtype) 
%   imports all fitted data for the given session and displays it. Session 
%   must be given as the numeric index corresponding to the DD structure.
%
%   [path, fitresults] = DOGEXPLORER(subject, 'Analyze', 'DataType', dtype) 
%   imports all fitted data for the given subject and displays it. Subject
%   must be given in the char/string form of the subject name.

function [path, varargout] = dogexplorer(subject, session, varargin)

if nargin == 0
    subject = [];
    session = [];
end

% Parse inputs
[DD, paths, subj, sess, dtype, MID, flags] = parseInputs(subject, session, varargin{:});

if ischar(flags) && strcmp(flags, 'return')
    path = paths;
    varargout{1} = DD;
    return
else
    clear subject session varargin
end

% Non-GRE but invalid measurement ID
if all(isnan(MID)) && ~strcmp(dtype, 'GRE')
    warning('The selected dataset is not available for the selected subject & imaging session');
    varargout{1} = []; varargout{2} = []; varargout{3} = [];
    return
end

switch dtype
    %% Gradient Echo
    case 'GRE'
        % Construct paths
        [path_raw, path] = GenPath(DD, paths, subj, sess, dtype);
        cd(path_raw(1))
        if flags.Load
            [snames, ~, HDR, ~, ~, vout, filelocs] = DogInitialize(paths.Main, subj, sess); 

        elseif flags.Preprocess
            % Load magnitude data
            [img, HDR, ~] = loadima(path_raw);
        
            % Determine echoes
            TEs = zeros(1,HDR(end).EchoNumbers);
            for te = 1:size(img,3):size(img,3)*length(TEs)
                TEs(HDR(te).EchoNumbers) = HDR(te).EchoTime;
            end
        
            % Determine encoding steps
            PEsteps = HDR(1).NumberOfPhaseEncodingSteps;
            FEsteps = PEsteps/(HDR(1).PercentSampling/100);
        
            % Determine field of view
            FEfov = round(FEsteps*HDR(1).PixelSpacing(1));
            PEfov = round(FEfov*(HDR(1).PercentPhaseFieldOfView/100));
        
            % Determine B0 direction
            if any(strcmp(HDR(1).PatientPosition(1:2), {'FF', 'HF'})) % Means that B0 is along transverse
                % Calculate B0 direction as the unit vector of the vector between two slices [x,y,z], where 
                %       x is right to left
                %       y is anterior to posterior
                %       z is inferior to superior
                B0Direction = (HDR(2).ImagePositionPatient - HDR(1).ImagePositionPatient)./norm(HDR(2).ImagePositionPatient - HDR(1).ImagePositionPatient);
            else
                error('B0 direction cannot be determined from patient position.')
            end
        
            % Create data structure
                D = struct('RawDataPath', path_raw, ... 
                               'Patient', struct('Name', DD(subj, sess).Name, ...
                                                  'Sex', HDR(1).PatientSex, ...
                                            'BirthDate', HDR(1).PatientBirthDate, ...
                                               'Weight', HDR(1).PatientWeight), ...
                                    'TR', HDR(1).RepetitionTime/1000, ...
                                    'TE', TEs/1000, ...
                               'deltaTE', diff(TEs(1:2))/1000, ...
                              'Averages', HDR(1).NumberOfAverages, ...
                                    'B0', HDR(1).ImagingFrequency/42.5774780505984, ...
                                    'F0', HDR(1).ImagingFrequency, ...
                             'FlipAngle', HDR(1).FlipAngle, ...
                     'BandwidthPerPixel', HDR(1).PixelBandwidth, ...
         'InPlanePhaseEncodingDirection', 'ROW', ...
                           'FieldOfView', [FEfov, PEfov], ...
                            'MatrixSize', [FEsteps, PEsteps], ...
                             'VoxelSize', [FEfov/FEsteps, PEfov/PEsteps, HDR(1).SliceThickness], ...
                              'SliceGap', HDR(1).SpacingBetweenSlices, ...
                                  'Size', size(img), ...
                        'TrimmedIndices', zeros([3 2]), ...
                           'B0Direction', B0Direction, ...
                                 'Flags', struct('Interpolated', false, ...
                                                      'Trimmed', struct('ThroughPlane', false, ...
                                                                             'InPlane', false), ...
                                               'UnwrappedPhase', false, ...
                                        'CorrectedBipolarPhase', false, ...
                                               'CombinedEchoes', false, ...
                                       'CorrectedChemicalShift', false, ...
                                     'Removed2DBackgroundField', false, ...
                                     'Removed3DBackgroundField', false, ...
                                               'InvertedDipole', false), ...
                                  'Data', struct('Image', img));
    
            varargout{1} = D;
            varargout{2} = HDR;
    
        elseif flags.Fit
            error('GRE data cannot be fitted.')
        elseif flags.Analyze
        elseif flags.verbose
            fprintf('\nData associated with given input:\n%s', ls)
        end

    case 'Localizer'
        %% Localizer
        if flags.Fit
            error('Localizer data cannot be fitted.')
        elseif flags.Analyze
            error('Localizer data cannot be analyzed.')
        end
    
        % Create path
        [path_raw, path] = GenPath(DD, paths, subj, sess, dtype);
        clear paths
        cd(path_raw(1))
    
        if flags.Preprocess
            varargout{1} = loaddat(MID, 'Save');
            cd(path)
        elseif flags.verbose
            fprintf('\nData associated with given input:\n%s', ls)
        end
        
    case 'HydrogenNoise'
        %% Noise Scans
        if flags.Fit
            error('Hydrogen Noise data cannot be fitted.')
        elseif flags.Analyze
            error('Hydrogen Noise data cannot be analyzed.')
        end
    
        % Create path
        [path_raw, path] = GenPath(DD, paths, subj, sess, dtype);
        clear paths
        cd(path_raw(1))
    
        if flags.Preprocess
            varargout{1} = loaddat(MID, 'Save');
            cd(path)

            % Generate covariance data
            twix = op_averaging(varargout{1});

            Covariance = cov(twix.specs);
            save(fullfile(path, append("MID", string(MID), "_noiseCov.mat")),"Covariance")
        elseif flags.verbose
            fprintf('\nData associated with given input:\n%s', ls)
        end

    otherwise
        %% Spectra that can be processed
        [path_raw, path_save] = GenPath(DD, paths, subj, sess, dtype);
        clear paths
    
        if flags.verbose
            for i = 1:numel(sess)
                fprintf('%s %i %s\n', DD(subj(i),sess(i)).Name, sess(i), dtype)
            end
        end

        if flags.Preprocess
            %% Preprocess
            cd(path_raw(1))
    
            % Create file names
            fname_raw = append("MID", string(MID), "_raw.mat");
            fname_preproc_twix = append("MID", string(MID), "_preproc_twix.mat");
    
            % Create file paths
            path_save_raw = fullfile(path_save, fname_raw);
            path_save_preproc_twix = fullfile(path_save, fname_preproc_twix);

            % Check if files exist
            if ~exist(path_save_preproc_twix, 'file')
                procFLAG = 'y';
            elseif islogical(flags.Override)
                if flags.Override
                    procFLAG = 'y';
                else
                    procFLAG = 'n';
                end
            else % does exist but not automatically overridden
                dat = dir(path_save_preproc_twix);
                fprintf('\b: last edited on %s\n', dat.date)
                procFLAG = input('The selected dataset has already been processed. Would you like to override? (''y'' or ''n'') ','s');
            end
    
            if strcmpi(procFLAG,'y')
                % Load data file
                twix = loaddat(MID, 'Save');
                varargout{1} = twix;
                cd(path_save)
    
                % Grab options
                if ~islogical(flags.Parameters) && ~isempty(flags.Parameters)
                    opts = flags.Parameters;
                    opts.auto = true;
                else
                    switch dtype
                        case 'Phosphorus'
                            opts = struct('AddRCVRS', false, ...
                                   'PhaseCorrection', true, ...
                                     'RMBadAverages', true, ...
                                    'FreqCorrection', true, ...
                            'SeparateFreqCorrection', true, ...
                                           'Average', false, ...
                                         'LeftShift', true, ...
                                      'ECCorrection', false, ...
                                      'WaterRemoval', false, ...
                                            'Plot3D', false, ...
                                        'SaveFormat', 'm', ...
                                              'auto', true);
                        case 'Hydrogen'
                            opts = struct('AddRCVRS', false, ...
                                   'PhaseCorrection', true, ...
                                     'RMBadAverages', true, ...
                                    'FreqCorrection', true, ...
                            'SeparateFreqCorrection', false, ...
                                     'AlignChannels', false, ...
                                           'Average', false, ...
                                         'LeftShift', true, ...
                                      'ECCorrection', true, ...
                                      'WaterRemoval', false, ...
                                            'Plot3D', false, ...
                                        'SaveFormat', 'm', ...
                                              'auto', true);
                        case 'HydrogenWatSup'
                            opts = struct('AddRCVRS', false, ...
                                   'PhaseCorrection', true, ...
                                     'RMBadAverages', true, ...
                                    'FreqCorrection', true, ...
                            'SeparateFreqCorrection', false, ...
                                     'AlignChannels', false, ...
                                           'Average', false, ...
                                         'LeftShift', true, ...
                                      'ECCorrection', true, ...
                                      'WaterRemoval', false, ...
                                            'Plot3D', false, ...
                                        'SaveFormat', 'm', ...
                                              'auto', true);
                        otherwise
                            opts = 'AutoRun';
                    end

                    if any(strcmp(dtype, {'Hydrogen','HydrogenWatSup'}))
                        % Manually identified bad averages
                        if     strcmp(DD(subj,sess).Name,'Aphrodite') && strcmp(DD(subj,sess).Date, "240125") && MID==87
                            opts.rmbadaverages.badaverages = [1:15 18:31 33:34 36 40 81 84]';
                        elseif strcmp(DD(subj,sess).Name,'Sushi')     && strcmp(DD(subj,sess).Date, "240923") && MID==187
                            opts.rmbadaverages.badaverages = [1:5 12:30 33 34]';
                        elseif strcmp(DD(subj,sess).Name,'Aphrodite') && strcmp(DD(subj,sess).Date, "240924") && MID==85
                            opts.rmbadaverages.badaverages = [1 3 7:10 12:25 27:29 29:42 53:56 58 59 61:63 65 86 89]';
                        elseif strcmp(DD(subj,sess).Name,'Waylon')    && strcmp(DD(subj,sess).Date, "240925") && MID==62
                            opts.rmbadaverages.badaverages = [1:8 16]';
                        elseif strcmp(DD(subj,sess).Name,'Selene')    && strcmp(DD(subj,sess).Date, "240926") && MID==109
                            opts.rmbadaverages.badaverages = (91:95)';
                        elseif strcmp(DD(subj,sess).Name,'Waylon')    && strcmp(DD(subj,sess).Date, "250113") 
                            if     MID==136
                                opts.rmbadaverages.badaverages = [23 24 27 35:40 60:66 74 87 92 117:121]';
                            elseif MID==140
                                opts.rmbadaverages.badaverages = [2:12 67 73 74 76:81 84:87 90:97 100:108]';
                            end
                        elseif strcmp(DD(subj,sess).Name,'Selene')    && strcmp(DD(subj,sess).Date, "250115")
                            if MID==99
                                opts.rmbadaverages.badaverages = [38:48 50:65 67:72 75 81 83 86:88 109:111 113]';
                            elseif MID==101
                                opts.rmbadaverages.badaverages = [28 33 42 44 54 67]';
                            end
                        elseif strcmp(DD(subj,sess).Name,'Aphrodite') && strcmp(DD(subj,sess).Date, "250127") && MID==103
                            opts.rmbadaverages.badaverages = [1:15 17:19 21 23 24 28 31 33 67 100 101]';
                        elseif strcmp(DD(subj,sess).Name,'Waylon')    && strcmp(DD(subj,sess).Date, "250501") && MID==144
                            if     MID==144
                                opts.rmbadaverages.badaverages = [1:60 62 63]';
                            elseif MID==142
                                opts.rmbadaverages.badaverages = [63:66 68:73 90 95]';
                            end
                        elseif strcmp(DD(subj,sess).Name,'Aphrodite') && strcmp(DD(subj,sess).Date, "250505") && MID==139
                            opts.rmbadaverages.badaverages = [2:74 76:80 85 87 89 101 104 105 112 113 128:128]';
                        elseif strcmp(DD(subj,sess).Name,'Sushi')     && strcmp(DD(subj,sess).Date, "250506") && MID==95
                            opts.rmbadaverages.badaverages = [1 4:6 11:13 19:40 71 74 75 78 80 84 86 90 91 94 102 113:116 118:128]';
                        end
                        if isfield(opts, 'rmbadaverages')
                            badavg = opts.rmbadaverages.badaverages;
                            opts = rmfield(opts, 'rmbadaverages');
                            opts.rmbadaverages = struct('domain', NaN, 'nsd', NaN, 'badaverages', badavg);
                        end

                        % Override residual water removal where applicable
                        if (strcmp(DD(subj,sess).Name,'Sushi')     && strcmp(DD(subj,sess).Date, "250129") && MID==111) || ...
                           (strcmp(DD(subj,sess).Name,'Aphrodite') && strcmp(DD(subj,sess).Date, "250127") && MID==103)
                            opts.WaterRemoval = true;
                        end
                    end
                end

                % Preprocess data
                if strcmp(dtype, 'Phosphorus')
                    varargout{2} = preprocFIDA(MID, varargout{1}, opts);
                else
                    [varargout{2}, proc_strct_w] = preprocFIDA(MID, varargout{1}, opts);
                    if ~isempty(proc_strct_w)
                        varargout{3} = proc_strct_w;
                    end
                    % Manually shift reference ppm points
                    if ~isfield(varargout{2}, 'ReferencePPM')
                        [~,~,~,~,varargout{2}] = FixPPM0(DD(subj,sess).Name,DD(subj,sess).Date,MID,path_save);
                    end
                end

                % Set name
                FIGS = findall(gcf);
                if numel(FIGS) > 1
                    for fg = 2:numel(FIGS)
                        if strcmp(get(FIGS(fg),'Type'),'tiledlayout')
                            title(FIGS(fg), sprintf('%s %i %s', DD(subj,sess).Name, sess, dtype))
                            break
                        end
                    end
                end

            else % Do not preprocess, just load files
                varargout{1} = load(path_save_raw);
                varargout{2} = load(path_save_preproc_twix);
            end
            path = path_save;
    
        elseif flags.Fit
            %% Fit
            cd(path_raw(1))
            opts = [];
    
            % Load data
            if strcmp(dtype, 'Phosphorus')
                path_raw = fullfile(path_save, append("MID", string(MID), "_preproc.mat"));
                load(path_raw, 'Starting')
                path_raw = fullfile(path_save, append("MID", string(MID), "_preproc_phased.mat"));
                load(path_raw, 'Processed')
                path_raw = fullfile(path_save, append("MID", string(MID), "_preproc_twix.mat"));
                load(path_raw, 'twix')

                twix.fids = squeeze(permute(Starting.Data, [3 1 2])); % Reshape starting data into twix format
                twix.specs = fftshift(ifft(twix.fids,[],twix.dims.t),twix.dims.t);
                if isfield(Processed.Operations, 'FrequencyShift')
                    twix = op_freqshift(twix, Processed.Operations.FrequencyShift/2); % Apply frequency & phase shifts
                end
                twix = op_addphase(twix, Processed.A0, 0);
                clear Starting Processsed

                if ~isfield(twix, 'refppm')
                    twix.refppm = 0;
                end

            elseif any(strcmp(dtype, {'Hydrogen', 'HydrogenWatSup'}))
                path_raw = fullfile(path_save, append("MID", string(MID), "_proc.mat"));
                twix = load(path_raw, 'McMRSData'); twix = twix.McMRSData;

                % Check if PPM0 corrected
                if ~isfield(twix, 'ReferencePPM') || twix.ReferencePPM == 4.65
                    [~,~,twix,~,~] = FixPPM0(DD(subj,sess).Name,DD(subj,sess).Date,MID,path_save);
                end

                % Add temporary flag for OXSA
                if strcmp(dtype,'HydrogenWatSup')
                    twix.flags.isWaterSuppressed = true;
                else
                    twix.flags.isWaterSuppressed = false;
                end
            end

            % Check for IV changes
            path_opts = fullfile(path_save, append("MID", string(MID), "_OXSA_opts.mat"));
            if exist(path_opts, "file")
                load(path_opts, 'opts')
            end
    
            % Fix outdated options files
            if ~isempty(opts)
                fields = {'IVs', 'Bs', 'PKs'};
                % if ~isfield(opts.KeepPeaks, 'Unknown')
                %     for j = 1:length(fields)
                %         if isfield(opts,fields{j}) && ~isempty(opts.(fields{j}))
                %             for i = 1:length(opts.(fields{j}))
                %                 if opts.(fields{j})(i).Index >= 12
                %                     opts.(fields{j})(i).Index = opts.(fields{j})(i).Index + 1;
                %                 end
                %             end
                %         end
                %     end
                %     opts.KeepPeaks.Unknown = false;
                %     save(path_opts, 'opts')
                % end
    
                if strcmp(dtype,'Phosphorus')
                    pks = ["bATP","UDPG","NAD","NADH","UDPG","aATP","gATP","PCr","MP","GPC","GPE","Unknown","Pia","Pib","DPG","G1P","DPG","PC","PE","G6P","Ref"];
                elseif any(strcmp(dtype,{'Hydrogen','HydrogenWatSup'}))
                    pks = ["IMCL_CH3","EMCL_CH3","IMCL_CH2","EMCL_CH2","Lip20","AcC","Lip22","Lip27","TMA","Cr3","Cho1","Tau11","Tau12","Tau21","Tau22","Cho3","Cr2","Cho2","Lip41","Water","Lip51","Lip52","Car_C4","Car_C2"];
                end

                for j = 1:length(fields)
                    if isfield(opts, fields{j})
                        if isfield(opts.(fields{j}), 'Index')
                            for i = 1:length(opts.(fields{j}))
                                opts.(fields{j})(i).Peak = pks(opts.(fields{j})(i).Index);
                            end
                            opts.(fields{j}) = rmfield(opts.(fields{j}), 'Index');
                        end
                    else
                        opts.(fields{j}) = [];
                    end
                    save(path_opts, 'opts')
                end
    
                if strcmp(dtype, 'Phosphorus')
                    if isfield(opts.KeepPeaks, 'Ref')
                        opts.KeepPeaks.PPA = opts.KeepPeaks.Ref;
                        opts.KeepPeaks = rmfield(opts.KeepPeaks, 'Ref');
                        for i = 1:numel(fields)
                            for j = 1:length(opts.(fields{i}))
                                if strcmp(opts.(fields{i})(j).Peak, 'Ref')
                                    opts.(fields{i})(j).Peak = 'PPA';
                                end
                            end
                        end
                        save(path_opts, 'opts')
                    end
                end
            end
    
            % Fit data
            [fitResults, fitSDs, fitStatus, FQN, FQN_global, FQN_SNR, ~, data, pk] = OXSA.run(twix, opts);
   
            % Save
            Results = struct(...
                'Fit', data, ...
                'Parameters', rmfield(fitResults,'covariance'), ...
                'ParameterSDs', fitSDs, ...
                'Covariance', fitResults.covariance, ...
                'FQN', FQN, ...
                'GlobalFQN', FQN_global, ...
                'FQN_SNR', FQN_SNR, ...
                'Status', fitStatus, ...
                'Priors', pk, ...
                'BoundFlag', fitResults.boundFlag);
   
            title(sprintf('FQN=%.2f, FQN/SNR=%.5f, CoV=%0.2f', Results.GlobalFQN, Results.FQN_SNR, mean([Results.ParameterSDs.amplitude Results.ParameterSDs.linewidth Results.ParameterSDs.sigma]./[Results.Parameters.amplitude Results.Parameters.linewidth Results.Parameters.sigma],"all","omitmissing")))

            path = path_save;
            path_save = fullfile(path_save, append("MID", string(MID), "_OXSA_results.mat"));
            save(path_save, "Results")
            varargout{1} = Results;
    
        elseif flags.Analyze
            %% Analyze 
            % Remove any missing datasets
            keeps = true(size(MID));
            xTicks = 1:numel(MID);
            for i = 1:numel(MID)
                if isnan(MID(i))
                    keeps(i) = false;
                end
            end
            MID = MID(keeps);
            path_save = path_save(keeps);
            xTicks = xTicks(keeps);
            clear i keeps
    
            % Generate full path
            path_proc = fullfile(string(path_save), append("MID", string(MID), "_preproc_phased.mat"));
            path_proc_twix = fullfile(path_save, append("MID", string(MID), "_preproc_phased_twix.mat"));
            path_oxsa = fullfile(path_save, append("MID", string(MID), "_OXSA_results.mat"));
            
            nd = numel(MID); % number in dataset

            % Load data
            SNRs = zeros(nd,1);
            A1 = zeros(nd,1);
            for diDx = 1:nd
                % Load data
                load(path_proc(diDx), "Processed")
                load(path_oxsa(diDx), 'Results')
                load(fullfile(path_save(diDx), append("MID", string(MID(diDx)), "_preproc_twix.mat")), 'twix');

                % Shift data to center PCr at 0ppm
                refPPM = Results.Parameters.chemShift(AMARES.getActualRefPeakDx(Results.Status.pkWithLinLsq));
                Results.Parameters.chemShift = Results.Parameters.chemShift - refPPM;
                Results.Status.xFit(Results.Status.constraintsCellArray.chemShift{1}{2}:Results.Status.constraintsCellArray.chemShift{end}{2}) = Results.Status.xFit(Results.Status.constraintsCellArray.chemShift{1}{2}:Results.Status.constraintsCellArray.chemShift{end}{2}) - refPPM;
                Results.Status.exptParams.ppmAxis = Results.Status.exptParams.ppmAxis - refPPM;
                Processed.FreqAxis = Processed.FreqAxis - refPPM;
                Processed.NoiseRegionLimits = Processed.NoiseRegionLimits - refPPM.*Processed.NoiseRegionLimits;
                Processed.PeakLocation = Processed.PeakLocation - refPPM;
                twix.ppm = twix.ppm - refPPM;

                % Insert into arrays
                TWIX(diDx) = twix; %#ok<AGROW>
                SNRs(diDx) = Processed.SNR;
                A1(diDx) = Processed.A1;
                if isfield(Results.Parameters,'covariance')
                    Results.Covariance = Results.Parameters.covariance;
                    Results.Parameters = rmfield(Results.Parameters,'covariance');
                elseif isfield(Results.ParameterSDs,'covariance')
                    Results.Covariance = Results.ParameterSDs.covariance;
                    Results.ParameterSDs = rmfield(Results.ParameterSDs,'covariance');
                end
                if ~isfield(Results,'BoundFlag')
                    Results.BoundFlag = [];
                end
                if ~isfield(Results,'Priors')
                    Results.Priors = Results.Status.pkWithLinLsq;
                end
                RESULTS(diDx) = Results; %#ok<AGROW>
    
                % Clear for next iteration
                clear twix Processed Results refPPM
            end

            % Set outputs
            path = struct('twix', convertStringsToChars(path_proc_twix), 'OXSA', convertStringsToChars(path_oxsa));
            varargout{1} = TWIX;
            clear path_proc_twix path_oxsa path_proc path_raw path_save diDx

            % Reformat results
            [modelResults, peakResults, oldNames, pkNames] = ReformatResults(RESULTS, DD(subj(1),1).Name, [DD(subj(1),:).Date], dtype);
            varargout{2} = modelResults;
            varargout{3} = peakResults;

            %% Combine peaks
            % Loop through datasets
            dInfo = struct('n',0,'pNames',{{""}},'multDx',{{[]}},'peakDx',{{[]}});
            pnames_tmp = fieldnames(oldNames);
            for diDx = 1:nd
                dInfo(diDx).n = 0;

                % Extract peak names
                for pkDx = 1:size(peakResults,1)
                    % Check if peak is present
                    if ~isscalar(oldNames(diDx).(pnames_tmp{pkDx})) || ~strcmp(oldNames(diDx).(pnames_tmp{pkDx}).name,"")
                        dInfo(diDx).n = dInfo(diDx).n + 1;
                        dInfo(diDx).pNames{dInfo(diDx).n} = [oldNames(diDx).(pnames_tmp{pkDx}).name];
                        dInfo(diDx).multDx{dInfo(diDx).n} = {oldNames(diDx).(pnames_tmp{pkDx}).multDx};
                        dInfo(diDx).peakDx{dInfo(diDx).n} = [oldNames(diDx).(pnames_tmp{pkDx}).peakDx];
                    end
                end
            end
            clear diDx pkDx oldNames
            
            % Plot data
            if flags.verbose
                % Grab indices for peak names
                switch dtype
                    case 'Phosphorus'
                        pkN = dictionary(["ATP","tNAD","PCr","MP","PDE","Pi","PME","PPA"], 1:8);
                    case 'Hydrogen'
                        error('Not added yet!')
                end
                ns = RESULTS(1).Status.exptParams.samples; % number of points
                lb = (1./SNRs)/min(1./SNRs)*3;
                ppmAxis = zeros(nd,ns);
                specs_true = zeros(ns,nd);
                specs_fit = zeros(ns,nd);
                for diDx = 1:nd
                    exptParams = modelResults(diDx).exptParams;
                    ppmAxis(diDx,:) = exptParams.ppmAxis.';

                    % Construct bandwidth axis
                    if ~isfield(exptParams,'freqAxis')
                        exptParams.freqAxis = exptParams.ppmAxis*exptParams.imagingFrequency;
                    end
    
                    % Grab zero-order phase
                    if strcmp(dtype, 'Phosphorus')
                        zeroOrderPhaseRad = peakResults(pkN("PCr"),1,diDx).phase*pi/180;
                    else
                        error('Hydrogen content not yet added!')
                    end

                    % Combine zero- & first-order phase corrections
                    phaseCorrection = exp(-1i*(zeroOrderPhaseRad + 2*pi*exptParams.freqAxis*exptParams.beginTime));

                    % Phase correct fitted phase
                    inputSpec_phcor = modelResults(diDx).inputSpec.*phaseCorrection;
                    modelSpec_phcor = sum(modelResults(diDx).modelSpecs,2).*phaseCorrection;

                    % Apodize data
                    specs_true(:,diDx) = real(specFft(FilterFid(specInvFft(inputSpec_phcor), exptParams.timeAxis, lb(diDx))));
                    specs_fit(:,diDx) = real(specFft(FilterFid(specInvFft(modelSpec_phcor), exptParams.timeAxis, lb(diDx)))); % do the same to fit so LW is more comparable
                    clear exptParams inputSpec_phcor modelSpec_phcor
                end
                figure(Theme='light'); waterfall(repmat(xTicks.', [1 ns]), ppmAxis, specs_true.'); ylim([-20 25]); zlim tight;
                xlabel('Session'); ylabel('Frequency (ppm)'); zlabel('Amplitude (A.u.)'); colormap jet; view(gca,269.7,17); set(gca,'XTick',xTicks); title('Data');

                figure(Theme='light'); waterfall(repmat(xTicks.', [1 ns]), ppmAxis, specs_fit.'); ylim([-20 25]); zlim tight;
                xlabel('Session'); ylabel('Frequency (ppm)'); zlabel('Amplitude (A.u.)'); colormap jet; view(gca,269.7,17); set(gca,'XTick',xTicks); title('Fit');
            end
            
            % Display results
            if flags.verbose
                for diDx = nd:-1:1
                    truePeaks = [dInfo(diDx).pNames{:}];
                    fprintf('\n\nDataset %i', xTicks(diDx))
                    for i = 1:6
                        switch i
                            case 1
                                str1 = 'Peak';
                            case 2
                                str1 = 'chemShift';
                                strf = '%6.2f+/-%-7.2f';
                            case 3
                                str1 = 'amplitude';
                                strf = '%6.2f+/-%-7.2f';
                            case 4
                                str1 = 'linewidth';
                                strf = '%6.2f+/-%-7.2f';
                            case 5
                                str1 = 'phase';
                                strf = '%6.2f+/-%-7.2f';
                            case 6
                                str1 = 'sigma';
                                strf = '%6.2f+/-%-7.2f';
                        end
    
                        if i > 1
                            str2 = [peakResults(:,1,diDx).(str1)]';
                            str2(:,2) = [peakResults(:,1,diDx).([str1 'CRB'])]';
                        end

                        fprintf(['\n' pad(str1,12)])
                        
                        if i == 1
                            for piDx = numel(truePeaks):-1:1
                                fprintf(' %-15s', truePeaks(piDx));
                            end
                        else
                            for piDx = numel(truePeaks):-1:1
                                fprintf(strf, mean(str2(piDx,:),1));
                            end
                        end
                    end
                end
                fprintf('\n')
            end
        else
            if flags.verbose
                fprintf('\nData associated with given input:\n%s', ls)
            end
            path = path_save;
        end
end
%% Functions
    function [modelResults, peakResults, oldNames, pkNames_all] = ReformatResults(RES, subj, date, nuc)
        %{ 
            Reorganize data
        Global model and model quality will be stored in "modelResults"
        Peak-specific information will be stored in "peakResults"
        %} 

        % Set static information
        params = ["chemShift", "linewidth", "amplitude", "phase", "sigma"];
        params_all = [params, "chemShiftDelta"];
        if strcmp(nuc, 'Phosphorus')
            renamePeaks = dictionary(["bATP",  "NAD", "NADH", "tNAD", "aATP", "gATP", "PCr", "MP", "GPC", "GPE", "PDE", "Unknown", "Pia", "Pib", "Pi",  "PC",  "PE", "PME", "G6P", "PPA", "Ref"], ...
                                     [ "ATP", "tNAD", "tNAD", "tNAD",  "ATP",  "ATP", "PCr", "MP", "PDE", "PDE", "PDE", "Unknown",  "Pi",  "Pi", "Pi", "PME", "PME", "PME", "PME", "PPA", "PPA"]);
            pknames = unique(values(renamePeaks),'stable');
        else
            error('Hydrogen content not yet added!')
        end

        % Preallocate structures
        modelResults = struct('Subject',"");
        modelResults(1,numel(RES)).Subject = "";
        peakResults = struct('Subject',"", 'Session', "", 'Date', "", 'peakName', "", 'subPeakNames', "", 'multDx', {{}}, 'peakDx', []);
        ends = ["", "CRB", "IV", "B", "PK"];
        for par = 1:numel(params_all)
            for subpar = 1:numel(ends)
                if subpar < numel(ends)
                    peakResults.(join([params_all(par) ends(subpar)],"")) = [];
                else
                    peakResults.(join([params_all(par) ends(subpar)],"")) = "";
                end
            end
        end
        peakResults = rmfield(peakResults, {'chemShiftPK', 'chemShiftDeltaPK'});
        % peakResults(numel(pknames),size(modelResults,1),size(modelResults,2)).Subject = "";

        %% Loop through date
        subjFLAG = true;
        oldNames = struct();
        pkNames_all = {strings()};
        for y = 1:size(modelResults,2)
            dateFLAG = true;

            % Account for missing dates
            if ismissing(date(1))
                dy = y + 1;
            else
                dy = y;
            end

            % Append model results
            modelResults(1,y).Subject = subj;
            modelResults(1,y).Session = y;
            modelResults(1,y).Date = date(dy);
            modelResults(1,y).inputSpec = RES(y).Fit.inputSpec;
            modelResults(1,y).modelSpecs = RES(y).Fit.modelSpecs;
            modelResults(1,y).xFit = RES(y).Status.xFit;
            modelResults(1,y).residual = RES(y).Fit.residual;
            modelResults(1,y).noise_var = RES(y).Status.noise_var;
            modelResults(1,y).FQN = RES(y).GlobalFQN;
            modelResults(1,y).covariance = RES(y).Covariance;
            modelResults(1,y).constraintsCellArray = RES(y).Status.constraintsCellArray;
            if ~isfield(RES(y), 'Priors')
                modelResults(1,y).pk = RES(y).Status.pkWithLinLsq;
            else
                modelResults(1,y).pk = RES(y).Priors;
            end
            modelResults(1,y).exptParams = RES(y).Status.exptParams;

            %% Determine indices, and names, for individual peaks
            for p = 1:numel(pknames)
                oldNames(y).(pknames(p)) = struct('name',"",'multDx',[],'peakDx',[]);
            end
            subPeaks = strings();
            for p = 1:numel(modelResults(1,y).pk.bounds)
                subPeaks_tmp = modelResults(1,y).pk.bounds(p).peakName;
                if iscell(subPeaks_tmp)
                    subPeaks_tmp = cellfun(@(x) x(1:end-1), subPeaks_tmp, 'UniformOutput', false);
                end
                subPeaks = [subPeaks convertCharsToStrings(subPeaks_tmp)]; %#ok<AGROW>

                % Append pk structure index and name
                pkNames_all{y}(p) = subPeaks(end);
                if isscalar(oldNames(y).(renamePeaks(subPeaks(end)))) && strcmp(oldNames(y).(renamePeaks(subPeaks(end))).name,"")
                    oldNames(y).(renamePeaks(subPeaks(end)))(1).name = subPeaks(end);
                    oldNames(y).(renamePeaks(subPeaks(end)))(1).multDx = (numel(subPeaks)-numel(convertCharsToStrings(subPeaks_tmp))+1:numel(subPeaks)) - 1;
                    oldNames(y).(renamePeaks(subPeaks(end)))(1).peakDx = p;
                else
                    oldNames(y).(renamePeaks(subPeaks(end)))(end+1).name = subPeaks(end);
                    oldNames(y).(renamePeaks(subPeaks(end)))(end).multDx = (numel(subPeaks)-numel(convertCharsToStrings(subPeaks_tmp))+1:numel(subPeaks)) - 1;
                    oldNames(y).(renamePeaks(subPeaks(end)))(end).peakDx = p;
                end
            end

            %% Grab peak data

            % First, identify which peaks are actually present
            pDx_subs = {};
            mDx_subs = {};
            name_subs = {};
            sDx = 0;
            acsdDx = 0; % index for amplitude & chemShiftDelta values
            acsdDx_subs = {}; % used to correlate above amplitude index from above with specific peaks (amplitude index will be value, chemShiftDelta index will be value+1)
            derivedStrs = {};
            for p = 1:numel(pknames)
                for subP = 1:numel(oldNames(y).(pknames(p)))
                    if ~isempty(oldNames(y).(pknames(p))(subP).peakDx)
                        sDx = sDx + 1;
                        pDx_subs{sDx} = oldNames(y).(pknames(p))(subP).peakDx;
                        mDx_subs{sDx} = oldNames(y).(pknames(p))(subP).multDx;
                        name_subs{sDx} = oldNames(y).(pknames(p))(subP).name;

                        % Make derived string to calculate amplitude & chemShiftDelta
                        if numel(mDx_subs{sDx}) > 1
                            acsdDx = acsdDx + 2;
                            derivedStrs{acsdDx-1} = makeDerivedStr(name_subs(sDx), mDx_subs(sDx), "amplitude");
                            derivedStrs{acsdDx} = makeDerivedStr(name_subs(sDx), mDx_subs(sDx), "chemShiftDelta");

                            acsdDx_subs{sDx} = acsdDx - 1;
                        else
                            acsdDx_subs{sDx} = NaN;
                        end
                    end
                end
            end

            % Obtain amplitude & chemShiftDelta values
            [valMN_subs, valCRB_subs] = AMARES.estimateDerivedParamAndCRB(modelResults(1,y).pk, ...
                                                                          modelResults(1,y).xFit, ...
                                                                          modelResults(1,y).constraintsCellArray, ...
                                                                          modelResults(1,y).covariance, ...
                                                                          derivedStrs, {}, [], false);

            % Loop to organize into structure
            sDx = 0;
            for p = 1:numel(pknames) 
                peakResults(p,1,y).Subject = subj;
                peakResults(p,1,y).Session = y;
                peakResults(p,1,y).peakName = pknames(p);

                % Loop through subpeaks
                for subP = 1:numel(oldNames(y).(pknames(p)))
                    peakFLAG = true;

                    % Only proceed if peak is actually present
                    if ~isempty(oldNames(y).(pknames(p))(subP).peakDx)
                        sDx = sDx + 1;

                        % Insert into output structure
                        peakResults(p,1,y).multDx{subP} = mDx_subs{sDx};
                        peakResults(p,1,y).peakDx(subP) = pDx_subs{sDx};
                        peakResults(p,1,y).subPeakNames(subP) = name_subs{sDx};

                        % Grab other parameters
                        for par = 1:numel(params_all)
                            Par = params_all(par);
                            % Special cases
                            switch Par
                                case "amplitude"
                                    if numel(mDx_subs{sDx}) > 1
                                        valMN = valMN_subs(acsdDx_subs{sDx});
                                        valCRB = valCRB_subs(acsdDx_subs{sDx});
                                    else
                                        valMN = RES(y).Parameters.(Par)(mDx_subs{sDx}(1));
                                        valCRB = RES(y).ParameterSDs.(Par)(mDx_subs{sDx}(1));
                                    end
                                    valIV = modelResults(1,y).pk.initialValues(pDx_subs{sDx}).(Par);
                                    valB = modelResults(1,y).pk.bounds(pDx_subs{sDx}).(Par).';

                                case "chemShift"
                                    % If multiplet, always grab from first
                                    valMN = RES(y).Parameters.(Par)(mDx_subs{sDx}(1));
                                    valCRB = RES(y).ParameterSDs.(Par)(mDx_subs{sDx}(1));
                                    valIV = modelResults(1,y).pk.initialValues(pDx_subs{sDx}).(Par);
                                    valB = modelResults(1,y).pk.bounds(pDx_subs{sDx}).(Par).';

                                case "chemShiftDelta"
                                    % If multiplet, grab from middle peak or average (depending on triplet or doublet)
                                    if numel(mDx_subs{sDx}) > 1
                                        valMN = valMN_subs(acsdDx_subs{sDx}+1)*modelResults(1,y).exptParams.imagingFrequency; % Convert to Hz
                                        valCRB = valCRB_subs(acsdDx_subs{sDx}+1)*modelResults(1,y).exptParams.imagingFrequency;
                                        valIV = modelResults(1,y).pk.priorKnowledge(pDx_subs{sDx}).chemShiftDelta*modelResults(1,y).exptParams.imagingFrequency; % chemShiftDelta uses prior knowledge as initial value
                                        valB = modelResults(1,y).pk.bounds(pDx_subs{sDx}).(Par).'.*modelResults(1,y).exptParams.imagingFrequency;
                                    else
                                        valMN = NaN;
                                        valCRB = NaN;
                                        valIV = NaN;
                                        valB = [];
                                    end

                                otherwise
                                    valMN = RES(y).Parameters.(Par)(mDx_subs{sDx}(1));
                                    valCRB = RES(y).ParameterSDs.(Par)(mDx_subs{sDx}(1));
                                    valIV = modelResults(1,y).pk.initialValues(pDx_subs{sDx}).(Par);
                                    valB = modelResults(1,y).pk.bounds(pDx_subs{sDx}).(Par).';

                                    % chemShift & chemShiftDelta have no prior knowledge
                                    if isfield(modelResults(1,y).pk.priorKnowledge, join(["G_" Par],"")) && ... case for G_sigma, which is sometimes absent
                                            ~isempty(modelResults(1,y).pk.priorKnowledge(pDx_subs{sDx}).(join(["G_" Par],""))) % peak is grouped
                                        peakResults(p,1,y).(join([Par "PK"],""))(subP) = pkNames_all{y}(modelResults(1,y).pk.priorKnowledge(pDx_subs{sDx}).(join(["G_" Par],"")));
                                    else
                                        peakResults(p,1,y).(join([Par "PK"],""))(subP) = "";
                                    end
                            end
                            peakResults(p,1,y).(Par)(subP) = valMN;
                            peakResults(p,1,y).(join([Par "CRB"],""))(subP) = valCRB;
                            peakResults(p,1,y).(join([Par "IV"],""))(subP) = valIV;

                            % Correct for empty bounds
                            if ~isempty(valB)
                                peakResults(p,1,y).(join([Par "B"],""))(:,subP) = valB;
                            else
                                peakResults(p,1,y).(join([Par "B"],""))(:,subP) = [NaN; NaN];
                            end
                        end
                    else
                        if subjFLAG
                            fprintf('\nSubject %s', subj)
                            subjFLAG = false;
                        end
                        if dateFLAG
                            fprintf('\n\tDate %i\n\t\t', y)
                            dateFLAG = false;
                        end
                        if peakFLAG
                            fprintf('%s, ', pknames(p))
                            peakFLAG = false; %#ok<NASGU>
                        end
                    end
                end
            end
            % Rename peaks
            pkNames_all{y} = unique(renamePeaks(pkNames_all{y}),"stable");
        end

        %% Functions
        function CRB = EP(type, crb, cv, data)
            % Error propogation
            switch lower(type)
                case 'addition'
                    CRB = sqrt(sum(crb.^2) + 2*sum(cv));
                case 'subtraction'
                    CRB = sqrt(sum(crb.^2) - 2*sum(cv));
                case 'multiplication'
                    CRB = abs(prod(data))*sqrt(sum((crb./data).^2) + 2*sum(cv./(data(2:end).*data(1:end-1))));
                case 'division'
                    CRB = abs(prod(data))*sqrt(sum((crb./data).^2) - 2*sum(cv./(data(2:end).*data(1:end-1))));
            end
        end

        % Construct derived string
        function derivedStr = makeDerivedStr(pkNames, multDx, param)
            
            % Grab total number of peaks
            mults_pk = cellfun(@numel, multDx);
            pkNames = arrayfun(@(x,y) repmat(x,1,y), pkNames, mults_pk, 'UniformOutput', false);
            for subPk = 1:numel(mults_pk)
                if mults_pk(subPk) > 1
                    pkNames{subPk} = pkNames{subPk} + string(1:mults_pk(subPk));
                end
            end
            
            % Extract from the added bottom cell layer
            pkNames = horzcat(pkNames{:});
            
            % Add symbol
            switch param
                case {"chemShift","chemShiftDelta"}
                    pkNames = pkNames + "_cs";
                case "linewidth"
                    pkNames = pkNames + "_lw";
                case "amplitude"
                    pkNames = pkNames + "_am";
                case "phase"
                    pkNames = pkNames + "_ph";
                case "sigma"
                    pkNames = pkNames + "_sg";
            end
            
            % Add summation symbol, where applicable
            if strcmp(param, "chemShiftDelta")
                derivedStr = join(pkNames([numel(pkNames) 1]), " - ");
            else
                if ~isscalar(pkNames)
                    derivedStr = join(pkNames," + ");
                end
            end

            % Add divide by two if triplet chemShiftDelta
            if strcmp(param, "chemShiftDelta") && numel(pkNames)>2
                derivedStr = "(" + derivedStr + ")/" + string(numel(pkNames)-1);
            end
        end
    end

    % function baseline_smooth = GenBaseline(results, pDxs, mDxs, cDxs)
    %     % Generate array with phase data
    %     exParams = results.Status.exptParams;
    %     fitPhase = results.Parameters.phase;
    %     nP = numel(results.Status.pkWithLinLsq.bounds);
    % 
    %     % Generate phase correction
    %     phAll = GenPhase(fitPhase, pDxs, mDxs, cDxs, [exParams.samples nP]);
    % 
    %     % Set peaks to NaN
    %     nanPeaks = results.Fit.modelSpec.*exp(1i*phAll*pi/180);
    %     for P = 1:nP
    %         nanPeaks(pDxs{P}) = NaN;
    %     end
    % 
    %     % Replace area where peak used to be with sum of non-peak fits
    %     Dx_stop_prev = 1;
    %     while Dx_stop_prev < numel(nanPeaks) && any(isnan(nanPeaks(Dx_stop_prev:end)))
    %         % Find peak start & stop
    %         Dx_startP = find(isnan(nanPeaks(Dx_stop_prev:end)),1) + Dx_stop_prev-1;
    %         Dx_stopP = find(~isnan(nanPeaks(Dx_startP:end)),1) + Dx_startP-1 - 1;
    %         if any(isnan(nanPeaks(Dx_stopP+1:end)))
    %             Dx_start_nextP = find(isnan(nanPeaks(Dx_stopP+1:end)),1) + Dx_stopP;
    %         else
    %             Dx_start_nextP = numel(nanPeaks);
    %         end
    % 
    %         % Save sum of peaks with chemical shift outside of region
    %         ops = false(1,nP);
    %         otherSpec = zeros(exParams.samples,1);
    %         for op = 1:nP
    %             if pDxs{op}(1) > Dx_stopP || pDxs{op}(end) < Dx_startP
    %                 otherSpec = otherSpec + sum(results.Fit.modelSpecs(:,mDxs{op}),2);
    %                 ops(op) = true;
    %             end
    %         end
    % 
    %         % Generate phase correction without current peak
    %         % error('Peaks are not being selected correctly')
    %         fitPhaseOther = fitPhase(cell2mat(mDxs(ops)));
    %         mDxCurrent = mDxs(~ops);
    %         mDxCurrent = mDxCurrent{end}(end);
    %         mDxsOther = mDxs(ops);
    %         mDxOther = 0;
    %         for m = 1:numel(mDxsOther)
    %             if mDxsOther{m}(1) > mDxCurrent
    %                 mDxsOther(m:end) = cellfun(@(x) x-mDxCurrent+mDxOther, mDxsOther(m:numel(mDxsOther)), UniformOutput=false);
    %                 break
    %             else
    %                 mDxOther = mDxsOther{m}(end);
    %             end
    %         end
    %         phOther = GenPhase(fitPhaseOther, pDxs(ops), mDxsOther, cDxs(ops), [exParams.samples sum(ops)]);
    %         otherSpec = otherSpec.*exp(1i*phOther*pi/180);
    % 
    %         % Set replacement start & stop bounds
    %         if Dx_stop_prev == 1
    %             Dx_start = Dx_stop_prev;
    %         else
    %             Dx_start = round((Dx_stop_prev + Dx_startP)/2);
    %         end
    %         if Dx_start_nextP == numel(nanPeaks)
    %             Dx_stop = Dx_start_nextP;
    %         else
    %             Dx_stop = round((Dx_stopP + Dx_start_nextP)/2);
    %         end
    % 
    %         % Create scale array to soften transition from otherSpec to baseline
    %         sa_start = linspace(0, 1, Dx_startP-Dx_start+1).';
    %         sa_stop = linspace(0, 1, Dx_stop-Dx_stopP+1).';
    % 
    %         % Replace
    %         nanPeaks(Dx_start:Dx_startP) = otherSpec(Dx_start:Dx_startP).*sa_start + nanPeaks(Dx_start:Dx_startP).*flip(sa_start);
    %         nanPeaks(Dx_startP:Dx_stopP) = otherSpec(Dx_startP:Dx_stopP);
    %         nanPeaks(Dx_stopP:Dx_stop) = otherSpec(Dx_stopP:Dx_stop).*flip(sa_stop) + nanPeaks(Dx_stopP:Dx_stop).*sa_stop;
    % 
    %         Dx_stop_prev = Dx_stop;
    %     end
    % 
    %     % Smooth & re-apply phase
    %     phAll_smooth = smooth(phAll,100);
    %     baseline_smooth = smooth(nanPeaks,100);
    %     phcor = exp(1i*phAll_smooth*pi/180);
    %     inputSpec_ph = results.Fit.inputSpec.*phcor;
    %     modelSpec_ph = results.Fit.modelSpec.*phcor;
    % 
    %     % Correct baseline
    %     inputSpec_ph_cor = inputSpec_ph - baseline_smooth;
    %     modelSpec_ph_cor = modelSpec_ph - baseline_smooth;
    % end

    % Manually shift reference ppm points
    function varargout = FixPPM0(Name,Date,MID,path_save)
        % By default they are shifted by 4.65 ppm, but water does not always
        % land at that point for some reason.
        if strcmp(Name,'Waylon')        && strcmp(Date, "240709")
            newRefPPM = 4.5747;
        elseif strcmp(Name,'Aphrodite') && strcmp(Date, "240710")
            newRefPPM = 4.5985;
        elseif strcmp(Name,'Sushi')     && strcmp(Date, "240125")
            newRefPPM = 4.5826;
        elseif strcmp(Name,'Selene')    && strcmp(Date, "240711")
            newRefPPM = 4.5588;
        elseif strcmp(Name,'Aphrodite') && strcmp(Date, "240711")
            newRefPPM = 4.5430;
        elseif strcmp(Name,'Sushi')     && strcmp(Date, "240923")
            newRefPPM = 4.6064;
        elseif strcmp(Name,'Aphrodite') && strcmp(Date, "240924")
            newRefPPM = 4.5667;
        elseif strcmp(Name,'Waylon')    && strcmp(Date, "240925")
            newRefPPM = 4.5588;
        elseif strcmp(Name,'Selene')    && strcmp(Date, "240926")
            newRefPPM = 4.5747;
        elseif strcmp(Name,'Aphrodite') && strcmp(Date, "250127")
            newRefPPM = 4.6381;
        elseif strcmp(Name,'Sushi') && strcmp(Date, "250129")
            newRefPPM = 4.5033;
        else
            error('no reference ppm set')
        end

        % error('This may be saving unsuppressed 1H data as suppressed 1H data')

        newRefPPM = 4.65*2 - newRefPPM;
        Ps2change = fullfile(path_save, append("MID", string(MID), ["_raw.mat";"_preproc.mat";"_proc.mat"]))';
        v = 0;
        for P = Ps2change
            load(P,'McMRSData')
            for f = ["PeakLocation","NoiseRegionLimits"]
                if all(McMRSData.(f) > McMRSData.PPMAxis(1)) && all(McMRSData.(f) < McMRSData.PPMAxis(end))
                    McMRSData.(f) = McMRSData.(f) - 4.65 + newRefPPM;
                end
            end
            McMRSData.PPMAxis = McMRSData.PPMAxis - 4.65 + newRefPPM;
            McMRSData.ReferencePPM = newRefPPM;
            save(P,'McMRSData')

            v = v + 1;
            varargout{v} = McMRSData; %#ok<AGROW>
        end
        Ps2change = fullfile(path_save, append("MID", string(MID), ["_raw_twix.mat";"_preproc_twix.mat"]))';
        for P = Ps2change
            load(P,'twix')
            twix.ppm = twix.ppm - 4.65 + newRefPPM;
            twix.refppm = newRefPPM;
            save(P,'twix')

            v = v + 1;
            varargout{v} = twix;
        end
    end

    function fit = GenFitNoPhase(exptParams, constraintsCellArray, xFit)
        % Generate fit with zero phase
        xFit(constraintsCellArray.phase{1}{2}:constraintsCellArray.phase{end}{2}) = 0;
        [fit,~,~] = AMARES.makeModelFidAndJacobianReIm(xFit,constraintsCellArray,exptParams.beginTime,exptParams.dwellTime,exptParams.imagingFrequency,exptParams.samples, 'complexOutput', true);
        fit = specFft(fit);
    end

    function phAll = GenPhase(fitPhase, pDxs, mDxs, cDxs, dims)
        % Correct for any discontinuities (differences greater than pi) (excluding start & stop)
        for Px = 2:dims(2)
            if mean(fitPhase(mDxs{Px})) - mean(fitPhase(mDxs{Px-1})) > 180
                fitPhase(mDxs{Px}(1):end) = fitPhase(mDxs{Px}(1):end) - 360;
            elseif mean(fitPhase(mDxs{Px})) - mean(fitPhase(mDxs{Px-1})) < -180
                fitPhase(mDxs{Px}(1):end) = fitPhase(mDxs{Px}(1):end) + 360;
            end
        end
        
        % Construct array of phase
        phAll = zeros(dims(1),1);
        prevph0 = 0;
        prevstop = 1;
        for Px = 1:dims(2) % peaks
            ph0 = -mean(arrayfun(@(x) fitPhase(x), mDxs{Px}));
    
            % Check for overlap
            if Px == 1 || pDxs{Px}(1) > pDxs{Px-1}(end)
                newstart = pDxs{Px}(1);
            else
                newstart = min(cDxs{Px});
            end
            if Px == dims(2) || pDxs{Px}(end) < pDxs{Px+1}(1)
                newstop = pDxs{Px}(end);
            else
                newstop = max(cDxs{Px});
            end
            phAll(newstart:newstop) = ph0;
    
            % Interpolate between peaks
            if Px > 1
                phAll(prevstop:newstart) = linspace(prevph0, ph0, newstart-prevstop+1);
            end
    
            % Save for next iteration
            prevph0 = ph0;
            prevstop = newstop;
        end
    
        % Extrapolate ends based on discontinuities (slope greater than pi/2)
        dph = diff(phAll);
        discon = find(abs(dph) > 90); % points where discontinuity occurs
        slope = mean(dph(discon(1)+1:discon(2)-1));
        ptsLeft = 1:discon(1); % interpolated points
        ptsRight = discon(2)+1:numel(phAll);
        phAll(ptsLeft) = (-numel(ptsLeft):-1).*slope + phAll(discon(1)+1);
        phAll(ptsRight) = (1:numel(ptsRight)).*slope + phAll(discon(2));
    end

    function fid = FilterFid(fid, t, lb)
        % Reshape
        if size(fid,1) < size(fid,2)
            fid = fid.';
            flipFLAG = true;
        else
            flipFLAG = false;
        end
        if size(t,1) > size(t,2)
            t = t.';
        end
        sz = size(fid);
        t2 = 1/(pi*lb);

        % Create an exponential decay (lorentzian) filter
        lor = exp(-t/t2);

        % Make a bunch of vectors of ones that are the same lengths as each of the dimensions of the data. Store them in a cell array for ease of use.
        p = cell(1,length(sz));
        for n = 1:length(sz)
            p{n} = ones(sz(n), 1);
        end

        % Populate a filter array with the lorentzian that has the same dimensions as the data
        switch length(sz)
            case 1
                filt = lor;
            case 2
                [filt,~] = ndgrid(lor,p{2});
            case 3
                [filt,~,~] = ndgrid(lor,p{2},p{3});
            case 4
                [filt,~,~,~] = ndgrid(lor,p{2},p{3},p{4});
            case 5
                [filt,~,~,~,~] = ndgrid(lor,p{2},p{3},p{4},p{5});
        end

        % Multiply the data by the filter array
        fid = fid.*filt;

        if flipFLAG
            fid = fid.';
        end
    end

    function [path_raw, path_save] = GenPath(DD, paths, subj, sess, dtype)
        path_raw = strings(size(sess));
        path_save = strings(size(sess));
        switch dtype
            case 'GRE'
                for pDx = 1:length(subj)
                    path_raw(pDx) = fullfile(DD(subj(pDx),sess(pDx)).Path.DICOM, DD(subj(pDx),sess(pDx)).GRE);
                    path_save(pDx) = fullfile(DD(subj(pDx),sess(pDx)).Path.Main, paths.Hydrogen.Main, paths.Hydrogen.Images.Main, paths.Hydrogen.Images.GRE);
                end
            case 'Localizer'
                for pDx = 1:length(subj)
                    path_raw(pDx) = fullfile(DD(subj(pDx),sess(pDx)).Path.Main, paths.Hydrogen.Main);
                    path_save(pDx) = fullfile(path_raw(pDx), paths.Hydrogen.Images.Main);
                end
            case 'HydrogenNoise'
                for pDx = 1:length(subj)
                    path_raw(pDx) = fullfile(DD(subj(pDx),sess(pDx)).Path.Main, paths.Hydrogen.Main);
                    path_save(pDx) = fullfile(path_raw(pDx), paths.Hydrogen.Noise);
                end
            case 'Hydrogen'
                for pDx = 1:length(subj)
                    path_raw(pDx) = fullfile(DD(subj(pDx),sess(pDx)).Path.Main, paths.Hydrogen.Main);
                    path_save(pDx) = fullfile(path_raw(pDx), paths.Hydrogen.Unsuppressed);
                end
            case 'HydrogenWatSup'
                for pDx = 1:length(subj)
                    path_raw(pDx) = fullfile(DD(subj(pDx),sess(pDx)).Path.Main, paths.Hydrogen.Main);
                    path_save(pDx) = fullfile(path_raw(pDx), paths.Hydrogen.Suppressed);
                end
            case 'Phosphorus'
                for pDx = 1:length(subj)
                    path_raw(pDx) = fullfile(DD(subj(pDx),sess(pDx)).Path.Main, paths.Phosphorus);
                    path_save(pDx) = path_raw(pDx);
                end
        end
    end

    function [DD, paths, Subj, Sess, dtype, MID, flags] = parseInputs(subj, sess, varargin)
        % Initialize data
        [paths,Subjects,Sessions,Phosphorus,HydrogenWatSup,Hydrogen,HydrogenNoise,DD] = Initialize;

        % Check for no inputs
        if isempty(subj) && isempty(sess)
            disp('Dataset includes data from five dogs over the course of three imaging sessions.')
            for ddx = 1:size(Sessions,2)
                fprintf('Session %i:\n', ddx)
                T = table(Subjects, Sessions(:,ddx), Phosphorus(:,ddx), HydrogenWatSup(:,ddx), Hydrogen(:,ddx), HydrogenNoise(:,ddx), 'VariableNames', {'Name', 'Date', 'Phosphorus', 'HydrogenWatSup', 'Hydrogen', 'HydrogenNoise'});
                disp(T)
            end

            Subj = [];
            Sess = [];
            dtype = [];
            MID = [];
            flags = 'return';
            return
        end

        Sessions2 = append("20", Sessions);

        % Validate subject
        if isnumeric(subj)
            if subj < 1 || subj > 5 && ~isinteger(subj)
                error('Incorrect subject input.')
            else
                Subj = subj;
            end
        else
            subj = string(subj);
            if any(strcmpi(subj, Subjects))
                Subj = find(Subjects == subj);
            else
                error('Incorrect subject input.')
            end
        end

        % Validate session
        analyzeFLAG = false;
        if isnumeric(sess)
            if sess < 1 || sess > 5 && ~isinteger(sess)
                error('Incorrect session input.')
            else
                Sess = sess;
            end
        else
            sess = string(sess);
            % Check for analyze instead of session
            if any(strcmpi(sess, {'Analyze', 'Analysis'}))
                analyzeFLAG = true;
                if isnumeric(subj) % subj is actually sess
                    Subj = (1:size(DD,1))';
                    Sess = ones(size(Subj)).*subj;
                elseif any(strcmpi(subj, Subjects))
                    Sess = (1:size(DD,2))';
                    Subj = find(Subjects == subj);
                    Subj = ones(size(Sess)).*Subj;
                else
                    error('Incorrect subject or session input.')
                end
            elseif any(strcmpi(sess, Sessions))
                Sess = find(Sessions == sess);
            elseif any(strcmpi(sess, Sessions2))
                Sess = find(Sessions2 == sess);
            else
                error('Incorrect sessions input.')
            end
        end
    
        % Determine what type of input to grab
        dtypes = dictionary(["gradientecho", "gre", "GRE"], ["GRE", "GRE", "GRE"]);
        dtypes(["localizer", "loc"]) = "Localizer";
        dtypes(["31p", "31P","phosphorus"]) = "Phosphorus";
        dtypes(["hydrogennoise", "1hnoise"]) = "HydrogenNoise";
        dtypes(["hydrogensup", "1hsup", "hydrogensuppressed", "hydrogenwatsup", "1hwatsup", "hydrogenwatersuppressed", "1hwatersuppressed", "1hsuppressed"]) = "HydrogenWatSup";
        dtypes(["hydrogen", "1h", "hydrogennonsup", "1hnonsup", "hydrogenunsup", "1hunsup", "hydrogennonsuppressed", "1hnonsuppressed", "hydrogenunsuppressed", "1hunsuppressed"]) = "Hydrogen";
        
        if isnumeric(varargin{1}) % input is MID
            if ~analyzeFLAG
                % Extract MID
                idx = 1;
                MID = varargin{idx};
                if MID < 100, MID = ['0' char(MID)]; else, MID = char(MID); end
        
                % Check file
                dtype = [];
                paths_tmp = [paths.Phosphorus, paths.Hydrogen.Main];
                dtypes_tmp = ["Phosphorus", "Hydrogen"];
                for m = 1:numel(dtypes_tmp)
                    path_tmp = fullfile(DD(Subj,Sess).Path.Main, paths_tmp(m));
                    ls_tmp = ls(path_tmp);
                    for n = 1:size(ls_tmp,1) % iterate through folder contents
                        if strcmp(ls_tmp(n,11:13),MID)
                            dtype = char(dtypes_tmp(m));
                            break
                        end
                    end
                    if n ~= size(ls_tmp,1)
                        break
                    end
                end
                MID = double(MID);
        
                % Ensure valid input
                if isempty(dtype)
                    error('MID could not be located');
                % Determine specific dtype of hydrogen (if necessary)
                elseif strcmp(dtype,'Hydrogen')
                    ls_tmp2 = strip(ls_tmp(n,:));
                    if strcmp(ls_tmp2(end-8:end-4),'noise')
                        dtype = 'HydrogenNoise';
                    else
                        error('Not yet finished')
                    end
                end
            else
                error('Datatype must be given for analysis.')
            end
        else
            if (strcmpi(varargin{1}, 'datatype') && nargin >= 2)
                idx = 2;
            elseif isKey(dtypes, lower(varargin{1}))
                idx = 1;
            else
                error('Incorrect data input.')
            end

            % Determine dtype
            dtype = dtypes(lower(varargin{idx}));

            if strcmp(dtype, "GRE") && ismissing(DD(Subj,Sess).GRE)
                error('No GRE for selected subject & imaging session')
            end

            % Determine MID
            switch dtype
                case 'GRE'
                    MID = NaN;

                case 'Localizer'
                    path_tmp = fullfile(DD(Subj,Sess).Path.Main, paths.Hydrogen.Main, paths.Hydrogen.Images.Main);
                    ls_tmp = ls(path_tmp);
                    locs = false(size(ls_tmp,1), 1);
                    for n = 1:size(ls_tmp,1)
                        if contains(ls_tmp(n,:),'localizer')
                            locs(n) = true;
                        end
                    end
                    if sum(locs) == 1
                        MID = ls_tmp(locs==true, 11:13);
                    else
                        n = find(locs, 1, 'last');
                        MID = ls_tmp(n, 11:13);
                    end

                case {'Phosphorus', 'HydrogenNoise', 'HydrogenWatSup', 'Hydrogen'}
                    MID = zeros(size(Sess));
                    for pDx = 1:length(Sess)
                        MID(pDx) = DD(Subj(pDx),Sess(pDx)).(dtype);
                    end
            end
        end
    
        % Check for missing outputs
        if ~exist('dtype', 'var')
            error('Requested input type could not be determined.')
        elseif ~exist('MID', 'var')
            error('MID could not be determined.')
        elseif all(isnan(MID)) && ~strcmp(dtype, 'GRE') % Non-GRE but invalid measurement ID
            warning('The selected dataset is not available for the selected subject & imaging session');
        else
            % Check for flags
            flags = struct( ...
                'Load', false, ...
                'Preprocess', false, ...
                'Fit', false, ...
                'Override', NaN, ...
                'Analyze', analyzeFLAG, ...
                'Parameters', false, ...
                'verbose', true);
    
            while idx < length(varargin)
                idx = idx + 1;
                if isstring(varargin{idx}) || ischar(varargin{idx})
                    switch lower(varargin{idx})
                        case 'load'
                            flags.Load = true;

                        case 'preprocess'
                            flags.Preprocess = true;
        
                        case 'fit'
                            flags.Fit = true;

                        case 'override'
                            idx = idx + 1;
        
                            if islogical(varargin{idx})
                                flags.Override = varargin{idx};
                            else
                                error('Override input must be logical.')
                            end

                        case 'parameters'
                            idx = idx + 1;
    
                            if isstruct(varargin{idx})
                                flags.Parameters = varargin{idx};
                            else
                                error('Parameters input must be a structure')
                            end

                        case 'verbose'
                            idx = idx + 1;

                            if islogical(varargin{idx})
                                flags.verbose = varargin{idx};
                            else
                                error('Verbose input must be logical.')
                            end
                    end
                end
            end
    
            % Ensure no double flags
            if flags.Analyze
                flags.Preprocess = false;
                flags.Fit = false;
            elseif flags.Fit
                flags.Preprocess = false;
            end
        end

        function [paths,Subject,Session,Phosphorus,HydrogenWatSup,Hydrogen,HydrogenNoise,DD] = Initialize
            % Initialize data
            mainpth = "/scratch/user/apad/Fat_water_separation/DICOM_Files/";
            paths = struct('Main', mainpth, ...
                     'Phosphorus', "31P", ...
                       'Hydrogen', struct('Main', "1H", ...
                                    'Suppressed', "Suppressed", ...
                                  'Unsuppressed', "Unsuppressed", ...
                                         'Noise', "Noise", ...
                                        'Images', struct('Main', "Images", ...
                                                          'GRE', "GRE")), ...
                          'DICOM', "DICOM");
            
            % Sort data
            Subject = ["Waylon"; "Sushi"; "Selene"; "Aphrodite"; "EOS"];
            Sex = ["M"; "F"; "F"; "F"; NaN];
            Session = ["240124", "240709", "240925", "250113", "250501"; % Waylon
                            NaN, "240710", "240923", "250129", "250506"; % Sushi
                            NaN, "240711", "240926", "250115", "250502"; % Selene
                       "240125", "240711", "240924", "250127", "250505"; % Aphrodite
                       "240125",      NaN,      NaN,      NaN,      NaN];% EOS
            Phosphorus = [ 88,  64,  76, 148, 134; % Waylon
                          NaN,  94, 197, 119, 109; % Sushi
                          NaN,  68, 118, 122, 158; % Selene
                           93, 160,  96, 120, 172; % Aphrodite
                          154, NaN, NaN, NaN, NaN];% EOS
            HydrogenWatSup = [104,  93,  62, 140, 144; % Waylon
                              NaN, 137, 186, 111,  97; % Sushi
                              NaN,  89, 108,  99, 147; % Selene
                               87, 182,  85, 103, 146; % Aphrodite
                              NaN, NaN, NaN, NaN, NaN];% EOS
            Hydrogen = [105,  94,  63, 136, 142; % Waylon
                        NaN, 138, 187, 112,  95; % Sushi
                        NaN,  90, 109, 101, 145; % Selene
                         88, 183,  86, 104, 139; % Aphrodite
                        NaN, NaN, NaN, NaN, NaN];% EOS
            HydrogenNoise = [106,  97,  65, 142, 125; % Waylon
                             NaN, 140, 189, 114, 100; % Sushi
                             NaN,  92, 111, 103, 149; % Selene
                              89, 185,  88, 106, 155; % Aphrodite
                             NaN, NaN, NaN, NaN, NaN];% EOS
            GRE = [false,  true,  true,  true,  true; % Waylon
                   false,  true,  true,  true,  true; % Sushi
                   false,  true,  true,  true,  true; % Selene
                   false,  true,  true,  true,  true; % Aphrodite
                   false, false, false, false, false];% EOS
            
            % Combine into structure
            for I = 1:size(Session,1)
                for J = 1:size(Session,2)
                    DD(I,J).Name = Subject(I); %#ok<AGROW>
                    DD(I,J).Sex = Sex(I); %#ok<AGROW>
                    DD(I,J).Date = Session(I,J); %#ok<AGROW>
    
                    date = append("20", Session(I,J));
                    if sum(Session == Session(I,J),"all") > 1 % Two dogs in same 
                        
                        DD(I,J).Path.Main = fullfile(paths.Main, date, Subject(I)); %#ok<AGROW>
                        if strcmp(Subject(I), "Aphrodite")
                            DD(I,J).Path.DICOM = fullfile(paths.Main, "DICOM", append(date,"-1")); %#ok<AGROW>
                        elseif strcmp(Subject(I), "Selene")
                            DD(I,J).Path.DICOM = fullfile(paths.Main, "DICOM", append(date,"-2")); %#ok<AGROW>
                        end
                    else
                        DD(I,J).Path.Main = fullfile(paths.Main, date); %#ok<AGROW>
                        DD(I,J).Path.DICOM = fullfile(paths.Main, "DICOM", date); %#ok<AGROW>
                    end
    
                    DD(I,J).Phosphorus = Phosphorus(I,J); %#ok<AGROW>
                    DD(I,J).HydrogenWatSup = HydrogenWatSup(I,J); %#ok<AGROW>
                    DD(I,J).Hydrogen = Hydrogen(I,J); %#ok<AGROW>
                    DD(I,J).HydrogenNoise = HydrogenNoise(I,J); %#ok<AGROW>
                    if GRE(I,J)
                        DD(I,J).GRE = append("GRE2D_FATWATER_", upper(Subject(I)), "_0012"); %#ok<AGROW>
                    else
                        DD(I,J).GRE = string(nan); %#ok<AGROW>
                    end
                end
            end
    
            % Save
            save(fullfile(mainpth, "DD.mat"), "DD")
        end
    end
end
