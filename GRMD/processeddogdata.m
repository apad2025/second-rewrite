%% Request inputs
clear; clc;

% Set autoflag
autoFLAG = true;

% Initialize path
path = 'C:\Users\apad2\Desktop\Fat_water_separation\DICOM_Files';

switch autoFLAG
    case true
        ntype = '31P';
        extype = 'excel_oxsa';
        cont = true;

    case false
        % Determine nuclei
        fig = uifigure;
        selection = uiconfirm(fig, 'What nulcei?', 'Nuclei Selection', 'Icon', 'question', 'Interpreter', 'tex', ...
                              'Options', ["^{1}H", "^{31}P", "Cancel"], 'DefaultOption', 1, 'CancelOption', 3);
        cont = true;
        ntype = "Cancel";
        switch selection
            case {"1H", "31P"}
                ntype = char(selection);
            otherwise
                cont = false;
        end
        
        % Determine 1H data type, if applicable
        if cont && strcmp(ntype, '1H')
            fig = uifigure;
            selection = uiconfirm(fig, 'What type of data?', '^{1}H Data Selection', 'Icon', 'question', 'Interpreter', 'tex', ...
                                  'Options', ["No Suppression", "Water Suppression", "Noise", "Cancel"], 'DefaultOption', 1, 'CancelOption', 4);
            cont = true;
            Htype = "Cancel";
            switch selection
                case {"No Suppression", "Water Suppression", "Noise"}
                    Htype = selection;
                otherwise
                    cont = false;
            end
        end
        
        % Determine what to do with data
        if cont
            fig = uifigure;
            selection = uiconfirm(fig, 'What would you like to do with the data?', 'Export Options', 'Icon', 'question', ...
                                  'Options', ["Export preprocessed data to Excel", "Convert preprocessed data to jMRUI", "Combine preprocessed data in MATLAB", "Export fitted data to Excel", "Cancel"], 'DefaultOption', 4, 'CancelOption', 4);
            cont = true;
            extype = "Cancel";
            switch selection
                case "Export preprocessed data to Excel"
                    extype = 'excel_preproc';
                    disp('Preprocessed .mat files from McMRS shall be imported, combined, and exported to the "DogData" Excel file.')
                case "Convert preprocessed data to jMRUI"
                    extype = 'jmrui';
                    disp('Preprocessed .mat files from McMRS shall be imported and converted to the jMRUI .txt files.')
                case "Combine preprocessed data in MATLAB"
                    disp('Processed .txt and .mrui files from jMRUI shall be imported and combined into a MATLAB structure.')
                case "Export fitted data to Excel"
                    extype = 'excel_oxsa';
                    disp('Fitted .mat files from OXSA shall be imported, combined, and exported to the "DogData" Excel file.')
                otherwise
                    cont = false;
            end
        end
end

if cont
    %% Extract file paths
    folds_all = dir(path);
    
    % Find data folders
    folds = '';
    idx = 0;
    for i = 1:length(folds_all)
        if length(folds_all(i).name) == 8 % ensure folder names are correct length
            if strcmp(folds_all(i).name(1:3), '202') && folds_all(i).isdir % ensure date & is directory
                idx = idx + 1;
                folds(idx,:) = folds_all(i).name;
            end
        end
    end
    
    % Preallocate path array
    paths = strings(size(folds,1),1);
    
    % Extract data paths
    i = 0;
    while i < size(folds,1)
        i = i + 1;
        mdFLAG = false; % multiple dogs in a single day
    
         % Isolate date path
        paths(i) = append(path,"\",folds(i,:));
    
        % Verify subfile is present
        subfolds = dir(paths(i));
        for j = 1:length(subfolds), if strcmp(subfolds(j).name, ntype), break; end; end
    
        if ~strcmp(subfolds(j).name, ntype) % multiple dogs in a single day
            mdFLAG = true;
    
            % Duplicate date
            folds = [folds(1:i,:); folds(i:end,:)];
    
            % Append paths
            paths(i+1) = paths(i);
            paths(i) = append(paths(i),"\",subfolds(end-1).name,"\",ntype);
            paths(i+1) = append(paths(i+1),"\",subfolds(end).name,"\",ntype);
        else
            paths(i) = append(paths(i),"\",ntype);
        end
    
        % Grab subfolders, if applicable
        switch ntype
            case '1H'
                paths(i) = append(paths(i),"\",Htype);
                if mdFLAG, paths(i+1) = append(paths(i+1),"\",Htype); end
        end
    
        % Extra shift if multiple dogs
        if mdFLAG, i = i + 1; end
    end
    
    % Create file end names
    if strcmp(ntype, '1H') && ~strcmp(extype, 'excel_oxsa')
        if strcmp(Htype, "Water Suppression")
            ends = ["_WSVD", "_SNRw", "_ECC_SNRw"];
            n = 3;
        else
            ends = ["_WSVD", "_SNRw"];
            n = 2;
        end
    else
        ends = "";
        n = 1;
    end
    
    % Extract ID
    ID = zeros(1,size(folds,1));
    for i = 1:length(paths)
        % Ensure path exists
        if exist(paths(i),"dir")
            % Determine ID
            files = dir(paths(i));
            for j = 1:length(files)
                if length(files(j).name) > 5 && ~files(j).isdir % Check if data file
                    if strcmp(files(j).name(6),'_') % ID is under 100
                        ID(i) = str2double(files(j).name(4:5));
                    else % ID is over 100
                        ID(i) = str2double(files(j).name(4:6));
                    end
                    break
                end
            end
        else
            ID(i) = NaN;
        end
    end
    
    %% Import & process data
    % Preallocate data
    in = [];
    m = size(folds,1);
    switch extype
        case 'excel_preproc'
            in.S = zeros(n,m,5);
            in.N = zeros(n,m,5);
            in.LW = zeros(n,m,5);
        case 'jmrui'
            in = struct('FName', cell(n,m,5), ...
                        'NPts', cell(n,m,5), ...
                        'SInt', cell(n,m,5), ...
                        'F0', cell(n,m,5), ...
                        'B0', cell(n,m,5), ...
                        'PName', cell(n,m,5), ...
                        'ExpDate', cell(n,m,5), ...
                        'SigNames', cell(n,m,5), ...
                        'Sig', cell(n,m,5), ...
                        'OutName', cell(n,m,5));
        case 'matlab'
            in = struct();
        case 'excel_oxsa'
            M = 7; % mean & stdev for LW, amp, & chemshift. Only mean for integral
            in = struct('FQN', zeros(n,m*M,5), ...
                        'FQN_SNR', zeros(n,m*M,5), ...
                        'bATP', zeros(n,m*M,5), ...
                        'tNAD', zeros(n,m*M,5), ...
                        'aATP', zeros(n,m*M,5), ...
                        'gATP', zeros(n,m*M,5), ...
                        'PCr', zeros(n,m*M,5), ...
                        'MP', zeros(n,m*M,5), ...
                        'PDE', zeros(n,m*M,5), ...
                        'Unknown', zeros(n,m*M,5), ...
                        'Pia', zeros(n,m*M,5), ...
                        'Pib', zeros(n,m*M,5), ...
                        'PME', zeros(n,m*M,5), ...
                        'Ref', zeros(n,m*M,5));
    end
    
    switch extype
        case {'excel_preproc','jmrui'}
            revstr = 'Loading preprocessed .mat files... ';
        case 'matlab'
            revstr = 'Loading processed .txt and .mrui files... ';
        case 'excel_oxsa'
            revstr = 'Loading fitted .mat files... ';
    end

    % Import data
    i2 = 1;
    for i = 1:length(paths)
        % Ensure path exists
        if exist(paths(i),"dir")
            % Loop through preprocessed files
            for j = 1:n
                % Create full path
                switch extype
                    case 'excel_preproc'
                        fullpath = append(paths(i),"\MID",num2str(ID(i)),"_preproc_phased",ends(j),".mat");
                    case 'jmrui'
                        in(j,i).OutName = append("MID",num2str(ID(i)),"_preproc_phased",ends(j));
                        fullpath = append(paths(i),"\",in(j,i).OutName,".mat");
                    case 'matlab'
                        fullpath = append(paths(i),"\MID",num2str(ID(i)),"_preproc_phased",ends(j),"_AMARES.txt");
                    case 'excel_oxsa'
                        fullpath = append(paths(i),"\MID",num2str(ID(i)),"_OXSA_results.mat");
                end

                    % Extract data
                switch extype
                    case 'excel_preproc'
                        in = import4excel_preproc(fullpath, i, j, in);
    
                    case 'jmrui'
                        in(j,i) = import4jmrui(paths(i), path, fullpath, in(j,i)); %#ok<*SAGROW>
    
                    case 'matlab'
                        if i == 1 && j == 1
                            in = import4matlab(fullpath);
                        else
                            in(j,i) = import4matlab(fullpath);
                        end

                    case 'excel_oxsa'
                        in = import4excel_oxsa(fullpath, i2, j, M, in);
                        i2 = i2 + M;
                end
            end
        elseif strcmp(extype, 'excel_preproc')
            in = set2nan(in, i, n);
        elseif strcmp(extype, 'excel_oxsa')
            in = set2nan(in, i2, n, M);
        end
        revstr = UpdatePercent(i/length(paths)*100, revstr);
    end

    % Process data
    switch extype
        case 'excel'
            export2excel_preproc(path, ID, ntype, Htype, S, N, LW);
    
        case 'jmrui'
            export2jmrui(paths, ntype, n, in)

        case 'matlab'
            [f, fSTD, a, aSTD, l, lSTD, p, pSTD, A] = combinejmrui(in);
    end
end

%% Functions
function in = set2nan(in, i, n, M)
    fields = fieldnames(in);

    if nargin == 4
        for f = 1:length(fields)
            for J = 1:n
                in.(fields{f})(J,i:i+M-1) = NaN;
            end
        end
    else
        for f = 1:length(fields)
            for J = 1:n
                in.(fields{f})(J,i) = NaN;
            end
        end
    end
end

function in = import4excel_preproc(fullpath, i, j, in)
    % Load & extract data
    D = load(fullpath);
    if isfield(D.Processed, 'Peak_Maxes')
        in.S(j,i) = D.Processed.Peak_Maxes;
        in.N(j,i) = D.Processed.Real_noise;
        in.LW(j,i) = D.Processed.Linewidth_Hz;
    else
        in.S(j,i) = D.Processed.PeakMaxes;
        in.N(j,i) = D.Processed.RealNoise;
        in.LW(j,i) = D.Processed.Linewidth;
    end
end

function in = import4jmrui(pth, path, fullpath, in)
    % Load & extract data
    D = load(fullpath);
    in.FName = D.Processed.Filename;
    in.B0 = D.Processed.B0;
    in.Sig = D.Processed.Spectrum;
    in.NPts = length(in.Sig);
    if isfield(D.Processed, 'Peak_Maxes') % Old name format
        in.SInt = 1/D.Processed.BW_Hz*1e3;
        in.F0 = D.Processed.F0_MHz*1e6;
        in.SigNames = join(D.Processed.Source_files, ';');
    else % New name format
        in.SInt = (1/D.Processed.BW)*1e3;
        in.F0 = D.Processed.F0*1e6;
        in.SigNames = join(D.Processed.SourceFiles, ';');
    end

    % Determine date
    pth = char(pth);
    date = pth(length(path)+2:length(path)+9);
    in.ExpDate = string([date(5:6) '/' date(7:8) '/' date(1:4)]);

    % Determine patient name
    if strcmp(date(3:end), '240711') || strcmp(date(3:end), '240125')
        path_split = strsplit(pth,'\');
        in.PName = path_split(10);
    else
        switch date(3:end)
            case {'240124','240709','240925','250113'}
                in.PName = "Waylon";
            case {'240710','240923','250129'}
                in.PName = "Sushi";
            case {'240926','250115'}
                in.PName = "Selene";
            case {'240924','250127'}
                in.PName = "Aphrodite";
        end
    end
end

function in = import4matlab(fullpath)
    [in, ~] = loadAMARESresults(fullpath, 'noplot');
end

function in = import4excel_oxsa(fullpath, i, j, M, in)
    % Load & extract data
    load(fullpath, 'Results');
    in.FQN(j,i) = Results.FQN;
    in.FQN_SNR(j,i) = Results.FQN_SNR;

    % Range each peak
    [PDxs, MDxs] = RangePeak(Results.Parameters, Results.Status.exptParams, Results.Status.pkWithLinLsq.bounds, Results.Fit);

    % Determine true number of peaks
    bounds = Results.Status.pkWithLinLsq.bounds;
    pk = struct('bATP', 1, 'tNAD', 2, 'aATP', 3, 'gATP', 4, 'PCr', 5, 'MP', 6, 'PDE', 7, 'Unknown', 8, 'Pia', 9, 'Pib', 10, 'PME', 11, 'Ref', 12);
    Npk = length(fieldnames(pk));
    
    % Preallocate
    if M == 7
        chemShift = cell(Npk,1);
        chemShift_sd = cell(Npk,1);
        linewidth = cell(Npk,1);
        linewidth_sd = cell(Npk,1);
        amplitude = cell(Npk,1);
        amplitude_sd = cell(Npk,1);
        integral = cell(Npk,1);
    else
        error('add more or less parameters')
    end

    % Extract by peaks
    for nDx = 1:length(bounds)
        % Grab peak name for structure
        if iscell(bounds(nDx).peakName)
            pkName = bounds(nDx).peakName{1};
            pkName = pkName(1:end-1);
        else
            pkName = bounds(nDx).peakName;
        end

        % Correct pkName for PDE & PME
        if any(strcmp(pkName, {'MP', 'GPC', 'GPE'}))
            pkName = 'PDE';
            PDEPMEflag = true;
        elseif any(strcmp(pkName, {'x2_3_DPG2', 'G1P', 'x2_3_DPG1', 'PC', 'PE', 'G6P'}))
            pkName = 'PME';
            PDEPMEflag = true;
        else
            PDEPMEflag = false;
        end

        % Calculate j coupling constant
        if MDxs(1,nDx) ~= MDxs(2,nDx)
            jcoup = (Results.Parameters.chemShift(MDxs(1,nDx)+1)-Results.Parameters.chemShift(MDxs(1,nDx)))*Results.Status.exptParams.imagingFrequency;
        else
            jcoup = 0;
        end

        % Append data to params
        if ~PDEPMEflag
            chemShift{pk.(pkName)} = mean(Results.Parameters.chemShift(MDxs(1,nDx):MDxs(2,nDx)));
            chemShift_sd{pk.(pkName)} = Results.ParameterSDs.chemShift(MDxs(1,nDx));
            linewidth{pk.(pkName)} = Results.Parameters.linewidth(MDxs(1,nDx)) + (MDxs(2,nDx)-MDxs(1,nDx))*jcoup;
            linewidth_sd{pk.(pkName)} = Results.ParameterSDs.linewidth(MDxs(1,nDx));
            amplitude{pk.(pkName)} = max(Results.Parameters.amplitude(MDxs(1,nDx):MDxs(2,nDx)));
            amplitude_sd{pk.(pkName)} = max(Results.ParameterSDs.amplitude(MDxs(1,nDx):MDxs(2,nDx)));
            integral{pk.(pkName)} = CalcIntegral(Results.Fit.modelSpecs, Results.Fit.phase(MDx(1))*pi/180, exptParams, MDxs(:,nDx)', PDxs(:,nDx)');
        else
            chemShift{pk.(pkName)} = [chemShift{pk.(pkName)}, {mean(Results.Parameters.chemShift(MDxs(1,nDx):MDxs(2,nDx)))}];
            chemShift_sd{pk.(pkName)} = [chemShift_sd{pk.(pkName)}, {Results.ParameterSDs.chemShift(MDxs(1,nDx))}];
            linewidth{pk.(pkName)} = [linewidth{pk.(pkName)}, {Results.Parameters.linewidth(MDxs(1,nDx)) + (MDxs(2,nDx)-MDxs(1,nDx))*jcoup}];
            linewidth_sd{pk.(pkName)} = [linewidth_sd{pk.(pkName)}, {Results.ParameterSDs.linewidth(MDxs(1,nDx))}];
            amplitude{pk.(pkName)} = [amplitude{pk.(pkName)}, {max(Results.Parameters.amplitude(MDxs(1,nDx):MDxs(2,nDx)))}];
            amplitude_sd{pk.(pkName)} = [amplitude_sd{pk.(pkName)}, {max(Results.ParameterSDs.amplitude(MDxs(1,nDx):MDxs(2,nDx)))}];
            integral{pk.(pkName)} = [integral{pk.(pkName)}, {CalcIntegral(Results.Fit.modelSpecs, Results.Fit.phase(MDx(1))*pi/180, exptParams, MDxs(:,nDx)', PDxs(:,nDx)')}];
            error('figure out better approach here')
        end
    end

    % Convert to array
    for nDx = 1:Npk
        % Set peaks that were not present to NaN
        if chemShift{nDx}==0 && linewidth{nDx}==0 && amplitude{nDx}==0 && integral{nDx}==0
            chemShift{nDx} = '=NA()';
            chemShift_sd{nDx} = '=NA()';
            linewidth{nDx} = '=NA()';
            linewidth_sd{nDx} = '=NA()';
            amplitude{nDx} = '=NA()';
            amplitude_sd{nDx} = '=NA()';
            integral{nDx} = '=NA()';
        end

        in.(pkName)(j,i:i+M-1) = [chemShift;
                                  chemShift_sd;
                                  linewidth;
                                  linewidth_sd;
                                  amplitude;
                                  amplitude_sd;
                                  integral];
    end

    function A = CalcIntegral(Specs, zeroOrderPhaseRad, exptParams, MDx, PDx)
        % First, correct the peak phase
        if ~isfield(exptParams,'freqAxis')
            exptParams.freqAxis = exptParams.ppmAxis*exptParams.imagingFrequency;
        end

        spec = zeros(size(Specs(:,1)), 'like', Specs(:,1));
        for Dx = MDx(1):MDx(2)
            spec = spec + Specs(:,Dx).*exp(-1i*(zeroOrderPhaseRad + 2*pi*exptParams.freqAxis*exptParams.beginTime));
        end

        % Integrate (approximately)
        A = trapz(exptParams.ppmAxis(PDx(1):PDx(2)), spec(PDx(1):PDx(2)));
    end

    function [PDxs, MDxs, combSpecs] = RangePeak(fitResults, exptParams, bounds, data)
        nc = length(bounds);
        hzPerPoint = (1/exptParams.dwellTime)/exptParams.samples;
        normAmp = ((fitResults.amplitude - min(fitResults.amplitude))/(max(fitResults.amplitude) - min(fitResults.amplitude)) + 1)/2; % rescale amplitudes to range [0.5, 1]

        % Preallocate
        combSpecs = zeros(size(data.inputFid,1), nc);
        PDxs = zeros(2, nc); % peak indexes (along x-axis)
        MDxs = zeros(2, nc); % multiplet indexes (along results structures/arrays)

        for pDx = 1:nc
            if pDx > 1
                MDxs(1,pDx) = MDxs(2,pDx-1) + 1;
            else
                MDxs(1,pDx) = 1;
            end

            % Check if multiplet
            if iscell(bounds(pDx).peakName)
                MDxs(2,pDx) = MDxs(1,pDx) + length(bounds(pDx).peakName) - 1;
            else
                MDxs(2,pDx) = MDxs(1,pDx);
            end

            if MDxs(2,pDx) == MDxs(1,pDx)
                combSpecs(:,pDx) = data.modelSpecs(:,MDxs(2,pDx));
            else
                for Dx = MDxs(1,pDx):MDxs(2,pDx)
                    combSpecs(:,pDx) = combSpecs(:,pDx) + data.modelSpecs(:,Dx);
                end
            end

            % Calculate the points to plot over so we don't end up with lots of baselines overlapping
            if MDxs(2,pDx) == MDxs(1,pDx) % singlet
                [~,peakCentreIndex] =  min(abs(exptParams.ppmAxis - fitResults.chemShift(MDxs(2,pDx)))); % Find the point in ppmAxis that is most similar to chemical shift
            else % doublet & triplet
                [~,peakCentreIndex] =  min(abs(exptParams.ppmAxis - (fitResults.chemShift(MDxs(1,pDx)) + fitResults.chemShift(MDxs(2,pDx)))/2));
            end


            if MDxs(2,pDx) == MDxs(1,pDx)
                Jcoup = 0;
                AmpMult = normAmp(MDxs(1,pDx)); % amplitude multiplier so shorter peaks are not weighted more than taller peaks
            else
                Jcoup = (fitResults.chemShift(MDxs(1,pDx)+1)-fitResults.chemShift(MDxs(1,pDx)))*exptParams.imagingFrequency; % in Hz
                AmpMult = max(normAmp(MDxs(1,pDx):MDxs(2,pDx)));
            end
            TrueLW = fitResults.linewidth(MDxs(2,pDx)) + Jcoup*(MDxs(2,pDx) - MDxs(1,pDx)); % takes into account j-coupling
            PeakWidth = 2.5*round((TrueLW+2*sqrt(2*log(2))*fitResults.sigma(MDxs(2,pDx))/hzPerPoint)*AmpMult); %%%%%%%% WHAT IS SIGMA %%%%%%%%
            PDxs(:,pDx) = [floor(peakCentreIndex-PeakWidth) ceil(peakCentreIndex+PeakWidth)]; % linewidth x 5

            % Ensure indices stay within bounds
            if PDxs(1,pDx) < 1
                PDxs(1,pDx) = 1;
            end
            if PDxs(2,pDx) > exptParams.samples
                PDxs(2,pDx) = exptParams.samples;
            end
        end
    end
end

function export2excel_preproc(path, ID, ntype, Htype, in)
    S = in.S;
    N = in.N;
    LW = in.LW;

    % Determine export parameters
    ROW.ID = 1;
    ROW.S = 5; % BSC & no BSC
    ROW.N = 6;
    ROW.LW = 8;
    switch ntype
        case '31P'
            sheet = '31P';
        case '1H'
            ROW.S(2) = 10; % SNR weighted data
            ROW.N(2) = 11;
            ROW.LW(2) = 13;
    
            switch Htype
                case "Water Suppression"
                    sheet = '1H Suppressed';
                case "No Suppression"
                    sheet = '1H Unsuppressed';
            end
    end
    
    % Check for NaNs
    ID = num2cell(ID);
    S = num2cell(S);
    N = num2cell(N);
    LW = num2cell(LW);
    if any(isnan(ID))
        nans = isnan(ID);
        ID{nans} = '=NA()';
        S{:,nans} = '=NA()';
        N{:,nans} = '=NA()';
        LW{:,nans} = '=NA()';
    end

    % Write data
    file = [path '\DogData.xlsx'];
    revstr = 'Exporting preprocessed data to "DogData" Excel file... ';
    for i = 1:length(ROW.S)
        SN = [S(i,:); N(i,:)];
        writecell(ID,file,'Sheet',sheet,'Range',['C' num2str(ROW.ID) ':Q' num2str(ROW.ID)],'AutoFitWidth',false)
        writecell(SN,file,'Sheet',sheet,'Range',['C' num2str(ROW.S(i)) ':Q' num2str(ROW.N(i))],'AutoFitWidth',false)
        writecell(LW(i,:),file,'Sheet',sheet,'Range',['C' num2str(ROW.LW(i)) ':Q' num2str(ROW.LW(i))],'AutoFitWidth',false)
        revstr = UpdatePercent(i/length(ROW.S)*100, revstr);
    end
end

function export2jmrui(paths, ntype, n, in)
    % Select nucleus output
    switch ntype
        case "1H"
            ntype = 1;
        case "31P"
            ntype = 2;
        case "13C"
            ntype = 3;
        case "19F"
            ntype = 4;
        case "23Na"
            ntype = 5;
    end

    revstr = 'Exporting preprocessed data as .mrui files... ';
    for i = 1:length(paths)
        if exist(paths(i),"dir")
            for j = 1:n
                fileID = fopen(fullfile(paths(i),strcat(in(j,i).OutName,'.txt')),'w');
                % Generate Header
                fprintf(fileID,'jMRUI Data Textfile\r\n');
                fprintf(fileID,'\r\n');
                fprintf(fileID,'Filename: %s\r\n',in(j,i).OutName);
                fprintf(fileID,'\r\n');
                fprintf(fileID,'PointsInDataset: %d\r\n',in(j,i).NPts);
                fprintf(fileID,'DatasetsInFile: %d\r\n',1);
                fprintf(fileID,'SamplingInterval: %.4E\r\n',in(j,i).SInt);
                fprintf(fileID,'ZeroOrderPhase: 0E0\r\n');
                fprintf(fileID,'BeginTime: 0E0\r\n');
                fprintf(fileID,'TransmitterFrequency: %.4E\r\n',in(j,i).F0);
                fprintf(fileID,'MagneticField: %.4E\r\n',in(j,i).B0);
                fprintf(fileID,'TypeOfNucleus: %.0E\r\n',ntype);
                fprintf(fileID,'NameOfPatient: %s\r\n',in(j,i).PName);
                fprintf(fileID,'DateOfExperiment: %s\r\n',in(j,i).ExpDate);
                fprintf(fileID,'Spectrometer: Siemens Verio 3T VB19B \r\n');
                fprintf(fileID,'AdditionalInfo: %s\r\n', in(j,i).FName);
                fprintf(fileID,'SignalNames: %s\r\n', in(j,i).SigNames);
                fprintf(fileID,'\r\n');
                fprintf(fileID,'\r\n');
                fprintf(fileID,'Signal and FFT\r\n');
                fprintf(fileID,'sig(real)\t sig(imag)\t fft(real)\t fft(imag)\r\n');

                fprintf(fileID,'Signal 1 out of 1 in file\r\n');
                specs = in(j,i).Sig;
                fids = ifft(ifftshift(specs,2),[],2);
                for row = 1:in(j,i).NPts
                    fprintf(fileID,'%.4E\t %.4E\t %.4E\t %.4E\r\n', real(fids(row)), imag(fids(row)), real(specs(row)), imag(specs(row)));
                end
            end
        end
        revstr = UpdatePercent(i/length(ROW.S)*100, revstr);
    end
end

function export2excel_oxsa(path, ID, ntype, Htype, in)
    % Determine export parameters
    pk = {'bATP', 'tNAD', 'aATP', 'gATP', 'PCr', 'MP', 'PDE', 'Unknown', 'Pia', 'Pib', 'PME', 'Ref'};
    ROW.FQN = 18;
    ROW.FQN_SNR = 19;
    for i = 1:length(pk)
        ROW.(pk{i}) = 21 + i;
    end
    switch ntype
        case '31P'
            sheet = '31P';
        case '1H'
            ROWs = fieldnames(ROW);
            for i = 1:length(ROWs)
                ROW.(ROWs{i}) = ROW.(ROWs{i}) + 10;
            end
    
            switch Htype
                case "Water Suppression"
                    sheet = '1H Suppressed';
                case "No Suppression"
                    sheet = '1H Unsuppressed';
            end
    end

    % Write data
    file = [path '\DogData.xlsx'];
    revstr = 'Exporting preprocessed data to "DogData" Excel file... ';

    FQNdata = [in.FQN; in.FQN_SNR];
    PKdata = [in.bATP; in.tNAD; in.aATP; in.gATP; in.PCr; in.MP; in.PDE; in.Unknown; in.Pia; in.Pib; in.Pib; in.PME; in.Ref];

    writecell(ID,file,'Sheet',sheet,'Range',['C' num2str(ROW.ID) ':Q' num2str(ROW.ID)],'AutoFitWidth',false)
    writecell(SN,file,'Sheet',sheet,'Range',['C' num2str(ROW.S(i)) ':Q' num2str(ROW.N(i))],'AutoFitWidth',false)
    writecell(LW(i,:),file,'Sheet',sheet,'Range',['C' num2str(ROW.LW(i)) ':Q' num2str(ROW.LW(i))],'AutoFitWidth',false)
    UpdatePercent(1/length(ROW.S)*100, revstr);
end

function [f, fSTD, a, aSTD, l, lSTD, p, pSTD, A] = combinejmrui(in)
    % Determine size for preallocation & loops
    K = length(in(1,1).estimates);
    mets = string(K,1);
    for i = 1:size(in,1)
        for j = 1:size(in,2)
            if length(in(i,j).estimates) > K
                K = length(in(i,j).estimates);
                for k = 1:K, mets(k) = string(in(i,j).estimates(k).Name); end
            end
        end
    end

    % Preallocate arrays for data
    f = nan(size(in), K);
    fSTD = nan(size(in), K);
    a = nan(size(in), K);
    aSTD = nan(size(in), K);
    l = nan(size(in), K);
    lSTD = nan(size(in), K);
    p = nan(size(in), K);
    pSTD = nan(size(in), K);
    A = nan(size(in), K);

    for i = 1:size(in,1)
        for j = 1:size(in,2)
            kk = 0;
            while kk < length(in(i,j).estimates)
                kk = kk + 1;

                for k = 1:K
                    if strcmpi(in(i,j).estimates(kk).Name, mets(k))
                        num = length(in(i,j).estimates(kk).Frequency);
                        if num == 1
                                f(i,j,k) = in(i,j).estimates(kk).Frequency;
                                fSTD(i,j,k) = in(i,j).estimates(kk).FrequencySTD;
                                a(i,j,k) = in(i,j).estimates(kk).Amplitude;
                                aSTD(i,j,k) = in(i,j).estimates(kk).AmplitudeSTD;
                                l(i,j,k) = in(i,j).estimates(kk).Linewidth;
                                lSTD(i,j,k) = in(i,j).estimates(kk).LinewidthSTD;
                                p(i,j,k) = in(i,j).estimates(kk).Phase;
                                pSTD(i,j,k) = in(i,j).estimates(kk).PhaseSTD;
                        else % num == 2 || num == 3
                            if num == 2
                                f(i,j,k) = mean(in(i,j).estimates(kk).Frequency);
                            else % num == 3
                                f(i,j,k) = in(i,j).estimates(kk).Frequency(2);
                            end
                            fSTD(i,j,k) = in(i,j).estimates(kk).FrequencySTD(2);
                            a(i,j,k) = in(i,j).estimates(kk).Amplitude(2);
                            aSTD(i,j,k) = in(i,j).estimates(kk).AmplitudeSTD(2);
                            l(i,j,k) = in(i,j).estimates(kk).Linewidth(2);
                            lSTD(i,j,k) = in(i,j).estimates(kk).LinewidthSTD(2);
                            p(i,j,k) = in(i,j).estimates(kk).Phase(2);
                            pSTD(i,j,k) = in(i,j).estimates(kk).PhaseSTD(2);
                        end
                        A(i,j,k) = in(i,j).estimates(kk).Area;
                    end
                end
            end
        end
    end
end