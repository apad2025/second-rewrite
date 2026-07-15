function [snames, DD, HDR, pth_1H, plrange, vout, filelocs] = DogInitialize(pth_DD, dog, date, flags, cscFLAG)
    % Set defaults
    if nargin < 5
        cscFLAG = true;
        if nargin < 4
            flags = struct('trimmed', true, ...
                            'zipped', true, ...
                 'bipolarcorrection', struct('method', 'MEDI'), ...
                           'verbose', true, ...
                              'plot', true, ...
                            'nobone', false, ...
                        'unwrapping', struct('method', 'RegGrow', ...
                                             'subsample', 1, ...
                                             'corrected', true), ...
                      'cscorrection', struct('method', 'vlGC', ...
                                          'subsample', 1), ...
                         'bfremoval', struct('D2', struct('method','V-SHARP - SEPIA'), ...
                                             'D3', struct('method','V-SHARP - SEPIA')), ...
                   'dipoleinversion', struct('method', 'MEDI'));
        end
    end

    % Check for joint background removal
    if all(isfield(flags.bfremoval, {'D2','D3'})) && strcmp(flags.bfremoval.D2.method,'V-SHARP - STISuite') && strcmp(flags.bfremoval.D3.method,'V-SHARP - STISuite')
        jointD23FLAG = true;
    else
        jointD23FLAG = false;
    end

    % Import data
    load(append(pth_DD,"\DD"), 'DD')

    % Create path for 1H
    pth_1H = char(fullfile(DD(dog,date).Path.Main, '1H', 'Images', 'GRE'));

    cd(pth_1H)
    
    snames = CreateNames(flags, jointD23FLAG, cscFLAG);

    [dataFLAG, rawunzipFLAG] = Check4Preload(snames, jointD23FLAG, pth_1H);

    % Swap raw zipped data to unzipped, if necessary
    if rawunzipFLAG
        snames.Preprocessed = snames.Preprocessed(1:3);
        dataFLAG.Preprocessed = true;
    end

    [plrange, HDR, vout, filelocs] = DataLoad(dog, date, DD, flags, dataFLAG, snames, jointD23FLAG, pth_1H);

    %% Fix outdated saves
    for ii = 1:length(vout)
        Dtmp = vout{ii};
        if isfield(Dtmp.Flags, 'TrimmedFieldOfView')
            % Save trimmed settings
            trimmedFOV = Dtmp.Flags.TrimmedFieldOfView; 
            trimmedslices = false;
            if Dtmp.Size(3) < 50 && ~Dtmp.Flags.Interpolated
                trimmedslices = true;
            end
        
            % Update structure shape
            Dtmp.Flags = rmfield(Dtmp.Flags, 'TrimmedFieldOfView');
            Dtmp.Flags.Trimmed = struct('ThroughPlane', false, 'InPlane', false);
        
            % Append up-to-date values
            if trimmedFOV, Dtmp.Flags.Trimmed.InPlane = true; end
            if trimmedslices, Dtmp.Flags.Trimmed.ThroughPlane = true; end
        end
        
        % Update all time variables to second (if necessary)
        if Dtmp.TE > 1
            Dtmp.TR = Dtmp.TR/1000;
            Dtmp.TE = Dtmp.TE./1000;
            if isfield(Dtmp, 'deltaTE')
                Dtmp.deltaTE = Dtmp.deltaTE/1000;
            else
                % Dtmp.deltaTE = diff(Dtmp.TE(1:2));
                Dtmp.deltaTE = mean(diff(Dtmp.TE));
            end
        end

        % Check for data field
        if ~isfield(Dtmp, 'Data')
            Dtmp.Data = struct();
            DataFields = {'Image', 'Mask', 'WeightedMagnitude', 'NoiseSTD', 'UnwrappedPhase', 'Water', 'Fat', 'TotalField', 'R2StarMap', 'LocalField', 'BackgroundField', 'SusceptibilityMap', 'Error', 'UncorrectedTotalField'};
            for f = 1:length(DataFields)
                if isfield(Dtmp, DataFields{f})
                    Dtmp.Data.(DataFields{f}) = Dtmp.(DataFields{f});
                    Dtmp = rmfield(Dtmp, DataFields{f});
                end
            end
        end

        vout{ii} = Dtmp;
    end

    %% Functions
    % Create data names
    function snames = CreateNames(flags, jointD23FLAG, cscFLAG)
        % Initialize dictionaries
        dictBipCor = dictionary(["MEDI", "SEPIAslow", "SEPIAfast", "hernando", "none"], ...
                                [ "plp",       "bec",       "fbc",      "igc","nobpc"]);
        dictUnwrap = dictionary(["GraphCuts", "RegGrow"], ...
                                [       "gc",      "rg"]);
        dictCSCor = dictionary(["IGC", "MixedIGC", "BipolarIGC", "vlGC", "IDEAL-CE", "SPURS", "Hierarchical IDEAL", "Golden Section"], ...
                               ["igc",     "migc",       "bigc", "vlgc",  "idealce", "spurs",             "hideal",             "gs"]);
        dictBFR = dictionary(["PDF - MEDI", "PDF - QSMmaster", "SHARP", "V-SHARP - SEPIA", "V-SHARP - STISuite", "LBV"], ...
                             [   "medipdf",         "qsmmpdf"     "sh",        "sepiavsh",             "stivsh", "lbv"]);
        dictDInv = dictionary(["TKD - MRIscm", "TKD - SEPIA", "MEDI", "iLSQR", "StarQSM", "FANSI", "Closed Form L2", "NDI", "Magnitude Weighted L1"], ...
                              [      "scmtkd",    "sepiatkd", "medi", "ilsqr",    "star", "fansi",           "cfl2", "ndi",                  "mwl1"]);

        jointBCSFLAG = false;
        if strcmp(flags.cscorrection.method, "BipolarIGC")
            jointBCSFLAG = true;
        end

        % Create save names structure
        snames = struct('Preprocessed', 'RAW', ...
                    'BipolarCorrected', 'BPC', ...
                      'UnwrappedPhase', 'TF', ...
                         'CorrectedCS', 'FW', ...
                          'LocalField', struct('D2', 'LF2D', ...
                                               'D3', 'LF3D'), ...
                   'SusceptibilityMap', 'QSM');
        if jointD23FLAG
            snames.LocalField = 'LF3D';
        end
        
        base = '';

        % Preprocessed file name
        if flags.zipped
            base = 'zip';
        end
        if flags.nobone
            base = [base, 'NB'];
        end
        snames.Preprocessed = [snames.Preprocessed, base];

        % Bipolar correction file name
        if ~jointBCSFLAG
            if ~isempty(base)
                base = [char(dictBipCor(flags.bipolarcorrection.method)), '_', base];
            else
                base = char(dictBipCor(flags.bipolarcorrection.method));
            end
            snames.BipolarCorrected = [snames.BipolarCorrected, base];
        end
        
        % Unwrapped phase file name
        baseUW = char(dictUnwrap(flags.unwrapping.method));
        if isfield(flags.unwrapping, 'subsample') && flags.unwrapping.subsample == 2
            baseUW = [baseUW, 'SS'];
        end
        if flags.unwrapping.corrected
            baseUW = [baseUW, 'Cor'];
        end
        if ~isempty(base)
            baseUW = [baseUW, '_', base];
        end
        snames.UnwrappedPhase = [snames.UnwrappedPhase, baseUW];

        % Chemical shift corrected file name
        if cscFLAG
            baseCS = char(dictCSCor(flags.cscorrection.method));
            if flags.cscorrection.subsample == 2
                baseCS = [baseCS, 'SS'];
            end
            if ~isempty(base)
                baseCS = [baseCS, '_', base];
            end
            snames.CorrectedCS = [snames.CorrectedCS, baseCS];
        end

        % Local field
        base = char(dictBFR(flags.bfremoval.D2.method));
        if strcmp(flags.bfremoval.D2.method, 'V-SHARP - STISuite')
            snames.LocalField = [snames.LocalField, char(dictBFR(flags.bfremoval.D2.method))];
        else
            snames.LocalField.D2 = [snames.LocalField.D2, char(dictBFR(flags.bfremoval.D2.method))];
        end
        if ~jointD23FLAG
            if cscFLAG
                base = [base, '_', baseCS];
            else
                base = [base, '_', baseUW];
            end
            snames.LocalField.D2 = [snames.LocalField.D2, base];
            base = [char(dictBFR(flags.bfremoval.D3.method)), '_', base];
            snames.LocalField.D3 = [snames.LocalField.D3, base];
        else
            snames.LocalField = [snames.LocalField, base];
        end
        
        % Susceptibility map
        base = [char(dictDInv(flags.dipoleinversion.method)), '_', base];
        snames.SusceptibilityMap = [snames.SusceptibilityMap, base];
    end

    % Check for data
    function [dataFLAG, rawunzipFLAG] = Check4Preload(snames, jointD23FLAG, pth_1H)
        % Preallocate preload data flags
        rawunzipFLAG = false; % only used if checking for raw zipped data, but only raw unzipped data is found
        dataFLAG = struct('Preprocessed', false, ...
                      'BipolarCorrected', false, ...
                        'UnwrappedPhase', false, ...
                           'CorrectedCS', false, ...
                            'LocalField', struct('D2', false, ...
                                                 'D3', false), ...
                     'SusceptibilityMap', false);
        if jointD23FLAG
            dataFLAG.LocalField = false;
        end
        
        % Now check for saved data
        nfields = fieldnames(snames);
        for i = 1:numel(nfields)
            if isstruct(snames.(nfields{i}))
                nfields2 = fieldnames(snames.(nfields{i}));
                for j = 1:numel(nfields2)
                    if isfile(fullfile(pth_1H, [snames.(nfields{i}).(nfields2{j}), '.mat']))
                        dataFLAG.(nfields{i}).(nfields2{j}) = true;
                    elseif isfile(fullfile(pth_1H, [New2OldNames(snames.(nfields{i}).(nfields2{j})), '.mat']))
                        movefile(fullfile(pth_1H, [New2OldNames(snames.(nfields{i}).(nfields2{j})), '.mat']), fullfile(pth_1H, [snames.(nfields{i}).(nfields2{j}), '.mat']))
                        dataFLAG.(nfields{i}).(nfields2{j}) = true;
                    end
                end
            else
                if isfile(fullfile(pth_1H, [snames.(nfields{i}), '.mat']))
                    dataFLAG.(nfields{i}) = true;
                elseif isfile(fullfile(pth_1H, [New2OldNames(snames.(nfields{i})), '.mat']))
                    movefile(fullfile(pth_1H, [New2OldNames(snames.(nfields{i})), '.mat']), fullfile(pth_1H, [snames.(nfields{i}), '.mat']))
                    dataFLAG.(nfields{i}) = true;
                end
            end
        end

        % Check for raw unzipped data, if necessary
        if ~dataFLAG.Preprocessed && strcmp(snames.Preprocessed, 'RAWzip')
            if isfile(char(fullfile(pth_1H, [snames.Preprocessed(1:3), '.mat'])))
                rawunzipFLAG = true;
            end
        end

        function sname = New2OldNames(sname)
            snames_old = {'bec',     'fbc',     'gs_',       'sh_'};
            snames_new = {'sepslow', 'sepfast', 'goldsect_', 'sharp_'};
            for s = 1:length(snames_old)
                sname = replace(sname, snames_old{s}, snames_new{s});
            end
        end
    end

    % Load data
    function [plrange, HDR, vout, filelocs] = DataLoad(dog, date, DD, flags, dataFLAG, snames, jointD23FLAG, pth_1H)
        filelocs = struct('Raw', 0, ...
                 'Preprocessed', 0, ...
             'BipolarCorrected', 0, ...
               'UnwrappedPhase', 0, ...
                  'CorrectedCS', 0, ...
                 'LocalFieldD2', 0, ...
                 'LocalFieldD3', 0, ...
                'LocalFieldD23', 0, ...
            'SusceptibilityMap', 0);

        % Load data
        if flags.verbose
            dogdate = char(DD(dog,date).Date);
            infodisp = ['Selected dataset: ' char(DD(dog,date).Name) ' ' dogdate(3:4) '/' dogdate(5:6) '/20' dogdate(1:2)];
            disp(infodisp)
            clear dogdate infodisp
        end

        % Check for previous data
        vout = {};
        idx = 1;
        if dataFLAG.Preprocessed
            if flags.verbose; fprintf('Loading saved data: %s\n', snames.Preprocessed); end
            load(char(fullfile(pth_1H, [snames.Preprocessed, '.mat'])), 'D');
            load(char(fullfile(pth_1H, 'RAWHeader.mat')), 'HDR');
            filelocs.Preprocessed = idx;
            filelocs.Raw = idx;
            vout{idx} = D;

            nfields = fieldnames(snames);
            nfields = nfields(2:end);
            for i = 1:numel(nfields)
                if isstruct(snames.(nfields{i}))
                    nfields2 = fieldnames(snames.(nfields{i}));
                    for j = 1:numel(nfields2)
                        if dataFLAG.(nfields{i}).(nfields2{j})
                            if flags.verbose; fprintf('Loading saved data: %s\n', snames.(nfields{i}).(nfields2{j})); end
                            load(char(fullfile(pth_1H, [snames.(nfields{i}).(nfields2{j}), '.mat'])), 'D');
                            idx = idx + 1;
                            vout{idx} = D;
                            filelocs.([nfields{i} nfields2{j}]) = idx;
                        end
                    end

                elseif dataFLAG.(nfields{i})
                    if flags.verbose; fprintf('Loading saved data: %s\n', snames.(nfields{i})); end
                    load(char(fullfile(pth_1H, [snames.(nfields{i}), '.mat'])), 'D');
                    idx = idx + 1;
                    vout{idx} = D;

                    if strcmp(nfields{i}, 'LocalField') && jointD23FLAG
                        filelocs.LocalFieldD23 = idx;
                    else
                        filelocs.(nfields{i}) = idx;
                    end
                end
            end
        else
            [~, D, HDR] = dogexplorer(dog, date, 'gre', 'Preprocess');
            filelocs.Raw = idx;
            vout{idx} = D;

            % Save data
            if flags.verbose; disp('Saving data'); end
            save(char(fullfile(pth_1H, 'RAW.mat')), "D");
            save(char(fullfile(pth_1H, 'RAWHeader.mat')), 'HDR');
        end
        
        % Set plot range
        plrange = GenPlotRange(D.Size);
    end
end