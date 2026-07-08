function [data, exptParams, pk, opts, A2P, P2A, pkInfo, Params, ParamsFit, SNR, changedBs, csbIndices, noisevar, Groups, Gpeaks] = prep(twix, opts)
% Check if fid-A or McMRSData format

if isfield(twix, 'fids')
    nuc = twix.nucleus;

    if twix.dims.averages ~= 0
        twix = op_averaging(twix);
    end
    
    data.inputFid = double(twix.fids);
    data.inputSpec = specFft(data.inputFid);
    
    exptParams = struct('samples', twix.sz(twix.dims.t), ...
                        'imagingFrequency', twix.txfrq/1e6, ...
                        'timeAxis', twix.t.', ...
                        'dwellTime', twix.dwelltime, ...
                        'ppmAxis', linspace(twix.ppm(end),twix.ppm(1),twix.sz(twix.dims.t)).', ...
                        'beginTime', 450e-6, ...
                        'offset', -twix.refppm);
else
    nuc = twix.Nucleus;

    twix = proc.CombAcqs(twix);

    data.inputFid = squeeze(twix.TimeDomain);
    data.inputSpec = squeeze(twix.Spectrum);

    exptParams = struct('samples', twix.Size.Points, ...
                        'imagingFrequency', twix.Frequency, ...
                        'timeAxis', (0:1/twix.SpectralWidth:twix.Size.Points/twix.SpectralWidth-1/twix.SpectralWidth).', ...
                        'dwellTime', 1/twix.SpectralWidth, ...
                        'ppmAxis', linspace(twix.PPMAxis(1),twix.PPMAxis(end),twix.Size.Points).', ...
                        'beginTime', 0, ... 
                        'offset', -twix.ReferencePPM);
end

if isfield(twix,'fids')
    if strcmp(nuc, '1H')
        ppm = struct('minS', 0.5, 'maxS', 2, 'minN', min(exptParams.ppmAxis), 'maxN', -0.5);
    elseif strcmp(nuc, '31P')
        ppm = struct('minS', -2, 'maxS', 2, 'minN', min(exptParams.ppmAxis), 'maxN', min(exptParams.ppmAxis)*5/6);
    end
    SNR = op_getSNR(twix, ppm.minS, ppm.maxS, ppm.minN, ppm.maxN, true);
    clear ppm
else
    SNR = meas.StoN(twix);
end

% Center shift ppm axis if 1H data
if strcmp(nuc,'1H') && abs(exptParams.ppmAxis(round(numel(exptParams.ppmAxis)/2))) > 1
    exptParams.ppmAxis = exptParams.ppmAxis + exptParams.offset;
end

% Possible Options
%{ 
    FIELD             ASSOCIATED FUNCTION     DEFAULT VALUE               DESCRIPTION
  - quiet             amaresFit               false                       suppresses output
  - MaxFunEvals       lsqcurvefit             500                         max number of function evaluations
  - TolFun            lsqcurvefit             1e-6*sqrt(max(abs(fid)))    termination tolerance on the function value
  - MaxIter           lsqcurvefit             100                         max number of iterations
  - apodization       amaresPlot              30                          apodizes input & fit
  - xUnits            amaresPlot              'PPM'                       sets x-axis units
  - plotIndividual    amaresPlot              true                        plot individual peaks
  - plotResidual      amaresPlot              true                        plot residual
  - plotInitial       amaresPlot              true                        plot input
  - offset            amaresPlot              0                           applies frequency offset
  - xlims             amaresPlot              [-25 15]                    x-axis limits
  - overideXAxis      amaresPlot              N/A                         uses given x-axis instead of one provided via exptParams
  - hFig              amaresPlot              N/A                         figure number
  - firstOrder        amaresPlot              true                        applies 0 & 1st order phasing to input & fit

Possible peaks
 1) β-ATP     -16.2 (triplet)
 2) UDPG (2)   -9.7 (doublet)
 3) NAD        -8.3
 4) NADH       -8.2
 5) UDPG (1)   -8.0 (doublet)
 6) α-ATP      -7.6 (doublet)
 7) γ-ATP      -2.6 (doublet)
 8) PCr         0.0
 9) MP          2.3
10) GPC         3.0 (triplet)
11) GPE         3.5
12) Unknown     4.4 (doublet
13) Pia         4.9
14) Pib         5.2
15) 2,3-DPG (2) 4.9 (doublet)
16) G1P         5.0 (doublet)
17) 2,3-DPG (1) 5.8 (triplet)
18) PC          6.3 (triplet)
19) PE          6.8 (triplet)
20) G6P         7.0 (triplet)
21) PPA        18.4
%}

if ~isfield(twix.flags, 'isWaterSuppressed')
    twix.flags.isWaterSuppressed = false;
end

% Sort through options
[pk, opts, changedBs] = Initialize(opts, nuc, -exptParams.offset, twix.flags.isWaterSuppressed);
opts.nucleus = nuc;

% Find number of peaks
pkInfo.Np.pk = numel(pk.initialValues); % number of peaks, multiplets combined only in pk file
A2P = {}; % transfer function from Np_all to pkInfo.Np.pk
P2A = 1:pkInfo.Np.pk; % transfer function from pkInfo.Np.pk to Np_all
pkInfo.Np.all = 0; % number of peaks, all multiplets combined
pkInfo.Names.all = {};
pkInfo.Names.pk = cell(1,pkInfo.Np.pk);
peak = 0;
while peak < pkInfo.Np.pk % step through peaks w/ combined multiplets
    peak = peak + 1;

    % Check if first test case is cell array
    pkInfo.Names.pk{peak} = pk.priorKnowledge(peak).peakName;
    if iscell(pkInfo.Names.pk{peak})
        pkInfo.Names.pk{peak} = pkInfo.Names.pk{peak}{1}(1:end-1);
    end

    % Ensure current peak hasn't been checked
    if all(peak ~= cell2mat(A2P)) 
        % Append
        pkInfo.Np.all = pkInfo.Np.all + 1;
        A2P{pkInfo.Np.all} = peak; %#ok<AGROW>
        P2A(peak) = pkInfo.Np.all;
        pkInfo.Names.all{pkInfo.Np.all} = pkInfo.Names.pk{peak};
        if peak ~= pkInfo.Np.pk % verify we haven't reached the very end
            for peak2 = peak+1:pkInfo.Np.pk % step through remaining peaks
                % Check if second test case is cell array
                name2 = pk.priorKnowledge(peak2).peakName;
                if iscell(name2)
                    name2 = name2{1}(1:end-1);
                end

                if (strcmp(pkInfo.Names.pk{peak},'Tau111') && numel(name2)>=3 && strcmp(name2(1:3),'Tau')) || ... if peaks are Tau
                        (~strcmp(pkInfo.Names.pk{peak}(1),'X') && numel(pkInfo.Names.pk{peak})>=3 && ~strcmp(pkInfo.Names.pk{peak}(1:3),'Lip') && ... if peak1 is not unknown/lipid
                        (strcmp(pkInfo.Names.pk{peak}, name2(1:end-1)) || strcmp(pkInfo.Names.pk{peak}(1:end-1), name2(1:end-1)))) && ... if peak1 ends with 1 & peak2 minus number equals peak1 minus number
                        ~isempty(pk.priorKnowledge(peak).G_linewidth) && ~isempty(pk.priorKnowledge(peak2).G_linewidth) && pk.priorKnowledge(peak).G_linewidth==pk.priorKnowledge(peak2).G_linewidth && ... peaks are in linewidth group
                        ~isempty(pk.priorKnowledge(peak).G_amplitude) && ~isempty(pk.priorKnowledge(peak2).G_amplitude)  && pk.priorKnowledge(peak).G_amplitude==pk.priorKnowledge(peak2).G_amplitude % peaks are in amplitude group
                    A2P{pkInfo.Np.all} = [A2P{pkInfo.Np.all} peak2];
                    P2A(peak2) = pkInfo.Np.all;
                    pkInfo.Names.all{pkInfo.Np.all} = pkInfo.Names.all{pkInfo.Np.all}(1:end-1);
                end
            end
        end
    end
end

% Determine which peaks are multiplets
pkInfo.multiplets.all = zeros(pkInfo.Np.all,1,'int8');
pkInfo.multiplets.pk = ones(pkInfo.Np.pk,1,'int8');
for peak = 1:pkInfo.Np.pk
    if iscell(pk.initialValues(peak).peakName)
        pkInfo.multiplets.all(P2A(peak)) = pkInfo.multiplets.all(P2A(peak)) + length(pk.initialValues(peak).peakName);
        pkInfo.multiplets.pk(peak) = length(pk.initialValues(peak).peakName);
    else
        pkInfo.multiplets.all(P2A(peak)) = pkInfo.multiplets.all(P2A(peak)) + 1;
    end
end

Params = ["chemShift", "linewidth", "amplitude", "phase", "sigma"];
ParamsFit = ["chemShift","linewidth","sigma"]; % parameters used for initial fitting

% Scale amplitude
maxA = max(abs(data.inputFid));

% Grab reference amplitude
refDx = find([pk.priorKnowledge.refPeak],1,'first');
refA = pk.initialValues(refDx).amplitude;

% Scale IV amplitudes & set bounds
for peak = 1:pkInfo.Np.pk
    pk.initialValues(peak).amplitude = pk.initialValues(peak).amplitude/refA*maxA;
    pk.bounds(peak).amplitude = [pk.bounds(peak).amplitude(1),maxA];

    % If lower bound is 0, force to be higher
    if pk.bounds(peak).amplitude(1) == 0
        pk.bounds(peak).amplitude(1) = refA*maxA*1e-4;
    end
end

% Find indices of chemShift bounds
csbIndices = cell(1,pkInfo.Np.all);
for peak = 1:pkInfo.Np.all
    ppmmin = exptParams.samples;
    ppmmax = 0;
    for subpeak = 1:numel(A2P{peak})
        [~,min_tmp] = min(abs(exptParams.ppmAxis - pk.bounds(A2P{peak}(subpeak)).chemShift(1)));
        [~,max_tmp] = min(abs(exptParams.ppmAxis - pk.bounds(A2P{peak}(subpeak)).chemShift(2)));

        % Apply LW (in # points) modifier (add/subtract 1/2 LW to max/min) - 5/15/26
        min_tmp = min_tmp - round((pk.initialValues(A2P{peak}(subpeak)).linewidth+pk.initialValues(A2P{peak}(subpeak)).sigma)*exptParams.samples*exptParams.dwellTime/2);
        max_tmp = max_tmp + round((pk.initialValues(A2P{peak}(subpeak)).linewidth+pk.initialValues(A2P{peak}(subpeak)).sigma)*exptParams.samples*exptParams.dwellTime/2);
        
        if min_tmp < ppmmin
            ppmmin = min_tmp;
        end
        if max_tmp > ppmmax
            ppmmax = max_tmp;
        end
    end
    csbIndices{peak} = ppmmin:ppmmax;
end

% Calculate noise variance
if strcmp(nuc, '1H')
    noisevar = var([data.inputSpec(exptParams.ppmAxis < pk.initialValues(1).chemShift-2); ...
                    data.inputSpec(exptParams.ppmAxis > pk.initialValues(end).chemShift+1)]);
elseif strcmp(nuc, '31P')
    if strcmp(pk.initialValues(end).peakName, 'PPA')
        noisevar = var([data.inputSpec(exptParams.ppmAxis < -20); ...
                        data.inputSpec(exptParams.ppmAxis < pk.initialValues(end).chemShift-2 & exptParams.ppmAxis > 10); ...
                        data.inputSpec(exptParams.ppmAxis > pk.initialValues(end).chemShift+2)]);
    else
        noisevar = var([data.inputSpec(exptParams.ppmAxis < -20);
                        data.inputSpec(exptParams.ppmAxis > 10)]);
    end
    % noisevar = var(twix.fids(twix.t > 0.3 & twix.t < twix.t(round(twix.sz(twix.dims.t)*0.95))));
end

% Find peak groups
if isfield(pk.priorKnowledge, 'G_sigma')
    sigmaFLAG = true;
else
    sigmaFLAG = false;
end
Gpeaks = cell(1,pkInfo.Np.all);
for peak = 1:pkInfo.Np.all
    Gpeaks{peak} = unique([[pk.priorKnowledge(A2P{peak}).G_linewidth] [pk.priorKnowledge(A2P{peak}).G_amplitude]]);
    if sigmaFLAG
        Gpeaks{peak} = unique([Gpeaks{peak} [pk.priorKnowledge(A2P{peak}).G_sigma]]);
    end
end
Gpeaks(cellfun(@(x)isempty(x),Gpeaks)) = {NaN}; % set any ungrouped peaks to NaN

% Create group mask
Groups = false(pkInfo.Np.all);
for peak = 1:pkInfo.Np.all
    if ~isnan(Gpeaks{peak})
        Groups(Gpeaks{peak}(1),peak) = true;
    end
end

% Find any peaks with no groups
for peak = 1:pkInfo.Np.all
    if ~Groups(peak,peak) && all(~Groups(:,peak))
        Groups(peak,peak) = true;
    end
end

%% Find the reference peak, calculate the offset and apply it to the PK
if ~strcmp(nuc,'1H')
    dataSpec = specApodize(exptParams.dwellTime*(0:exptParams.samples-1).',specFft(data.inputFid),30);
    dataAbs = abs(dataSpec);
    
    % Create the initial value spectrum
    initialmodelFid = AMARES.makeInitialValuesModelFid(pk, exptParams);
    
    % Change to Frequency Domain
    initialmodel = abs(specFft(initialmodelFid));
    
    % Scale it to near the spectral amplitude.
    model = (initialmodel/max(abs(initialmodel)))*max(dataAbs);
    
    % Do the convolution and find the maximum point in the result.
    convolutionRes = conv(dataAbs,flipud(model),'same');
    [~,coarseIndex] = max(abs(convolutionRes));
    
    searchRange =  round((0.5/abs(exptParams.ppmAxis(1)-exptParams.ppmAxis(2)))/2);
    
    %Refine the point by searching around the location
    searchVec = [coarseIndex-searchRange coarseIndex+searchRange];
    if searchVec(1) < 1, searchVec(1) = 1; end
    if searchVec(2) > exptParams.samples, searchVec(2) = exptParams.samples; end
    
    [~,index] = max(dataAbs(searchVec(1):searchVec(2)));
    index = coarseIndex + (index - searchRange);
    if index < 1, index = 1; elseif index > exptParams.samples, index = exptParams.samples; end
    
    exptParams.offset = exptParams.ppmAxis(index);
    
    % Apply the offset. Removed the rounding!
    for peak = 1:pkInfo.Np.all
        pk.initialValues(peak).chemShift = (pk.initialValues(peak).chemShift + exptParams.offset);
        pk.bounds(peak).chemShift = (pk.bounds(peak).chemShift + exptParams.offset);
    end
end

%% Functions
    function [PK, opts, changedBs] = Initialize(opts, nuc, ppm0, wsFLAG)
        % Set frequency
        if strcmp(nuc,'1H')
            f = 123.2030;
        elseif strcmp(nuc,'31P')
            f = 49.8719;
        end

        % Set defaults
        if strcmp(nuc, '1H')
            xlims = [0 9];
        elseif strcmp(nuc, '31P')
            xlims = [-20 25];
        end

        opts_def = struct('IVs', [], 'Bs', [], 'PKs', [], 'quiet', true, 'apodization', false, 'xlims', xlims, 'firstOrder', true, 'sigmaMax', 20, 'MaxFunEvals', 8000, 'TolFun', 1e-17, 'MaxIter', 3200, 'TolX', 1e-8);

        if strcmp(nuc, '1H')
            opts_def.CombinePeaks = struct( ...
             'MCL_CH3', true, ...
             'MCL_CH2', true);
            opts_def.KeepPeaks = struct( ...
            'IMCL_CH3', true, ...
            'EMCL_CH3', true, ...
            'IMCL_CH2', true, ...
            'EMCL_CH2', true, ...
               'Lip20', false, ...
                 'AcC', true, ...
               'Lip22', true, ...
               'Lip27', false, ...
                  'X1', true, ...
                  'X2', true, ...
                  'X3', true, ...
                 'Cho', false, ...
                 'Tau', false, ...
                  'X4', true, ...
                'Cr21', true, ...
                'Cr22', true, ...
               'Lip41', true, ...
               'Water', true, ...
               'Lip49', false, ...
               'Lip52', true, ...
               'Lip53', true, ...
              'Car_C4', false, ...
              'Car_C2', false, ...
                  'X5', false);

        elseif strcmp(nuc, '31P')
            opts_def.CombinePeaks = struct( ...
                'tNAD', true, ...
                 'PDE', false, ...
                  'Pi', false, ...
                 'PME', false);
            opts_def.KeepPeaks = struct( ...
                'bATP', true, ...
                'UDPG', false, ...
                 'NAD', true, ...
                'NADH', true, ...
                'aATP', true, ...
                'gATP', true, ...
                 'PCr', true, ...
                  'MP', false, ...
                 'GPC', true, ...
                 'GPE', true, ...
             'Unknown', false, ...
                 'Pia', true, ...
                 'Pib', true, ...
                 'DPG', false, ...
                 'G1P', false, ...
                  'PC', true, ...
                  'PE', true, ...
                 'G6P', false, ...
                 'PPA', true);
        end
       
        % Construct dictionary of peaks
        if strcmp(nuc, '1H')
            pkDict = dictionary(["IMCL_CH3", "EMCL_CH3", "IMCL_CH2", "EMCL_CH2", "Lip20", "AcC", "Lip22", "Lip27", "X1", "X2", "X3",      "Cho",         "Tau", "X4", "Cr21", "Cr22", "Lip41", "Water", "Lip49", "Lip52", "Lip53", "Car_C4",   "Car_C2", "X5"], ...
                                {         1,          2,          3,          4,       5,     6,       7,       8,    9,   10,   11, [12,18,21], [13,14,16,17],   15,     19,     20,      22,      23,      24,      25,      26,       27,         28,   29});
        elseif strcmp(nuc, '31P')
            pkDict = dictionary(["bATP", "UDPG", "NAD", "NADH", "aATP", "gATP", "PCr", "MP", "GPC", "GPE", "Unknown", "Pia", "Pib",   "DPG", "G1P", "PC", "PE", "G6P", "PPA"], ...
                                {     1,  [2,5],     3,      4,      6,      7,     8,    9,    10,    11,        12,    13,    14, [15,17],    16,   18,   19,    20,    21});
        end
        pkInd = unique(keys(pkDict),'stable');
        acDict = dictionary(["IVs",           "Bs",     "PKs"], ...
                            ["initialValues", "bounds", "priorKnowledge"]); % acronym dict
        o2nDict = dictionary(["IMCL_CH3","EMCL_CH3","IMCL_CH2","EMCL_CH2", "NAD","NADH","GPC","GPE","Pia","Pib", "PC", "PE","G6P","DPG","G1P"],... dictionary for names before & after combining
                             [ "MCL_CH3", "MCL_CH3", "MCL_CH2", "MCL_CH2","tNAD","tNAD","PDE","PDE", "Pi", "Pi","PME","PME","PME","PME","PME"]);
        n2oDict = dictionary([             "MCL_CH3",              "MCL_CH2",        "tNAD",        "PDE",         "Pi",                       "PME"], ...
                            {["IMCL_CH3","EMCL_CH3"],["IMCL_CH2","EMCL_CH2"],["NAD","NADH"],["GPC","GPE"],["Pia","Pib"],["PC","PE","G6P","DPG","G1P"]});

        % Sort through inputs
        if isempty(opts)
            opts = opts_def;
            if strcmp(nuc, '1H')
                pkFile = AMARES.priorKnowledge.PK_1H_GRMD;
            elseif strcmp(nuc, '31P')
                pkFile = AMARES.priorKnowledge.PK_31P_GRMD_nophase;
            end
        else
            if isfield(opts, 'pkFile')
                pkFile = AMARES.priorKnowledge.(opts.pkFile);
                opts = rmfield(opts, 'pkFile');
            else
                if strcmp(nuc, '1H')
                    pkFile = AMARES.priorKnowledge.PK_1H_GRMD;
                    if isfield(opts.KeepPeaks, 'Car')
                        opts.KeepPeaks.Car_C4 = opts.KeepPeaks.Car;
                        opts.KeepPeaks.Car_C2 = opts.KeepPeaks.Car;
                        opts.KeepPeaks = rmfield(opts.KeepPeaks, 'Car');
                    end
                    if isfield(opts.KeepPeaks, 'Cr')
                        opts.KeepPeaks.Cr21 = opts.KeepPeaks.Cr;
                        opts.KeepPeaks.Cr22 = opts.KeepPeaks.Cr;
                        opts.KeepPeaks = rmfield(opts.KeepPeaks, 'Cr');
                    end
                elseif strcmp(nuc, '31P')
                    pkFile = AMARES.priorKnowledge.PK_31P_GRMD_nophase;
                end
            end
            f2c = fieldnames(opts_def); % fields to check
            for i = 1:numel(f2c)
                if ~isfield(opts, f2c{i})
                    opts.(f2c{i}) = opts_def.(f2c{i});
                end
            end
        end
        opts.plotInitial = false;
        changedBs = pkFile.bounds;

        if wsFLAG % Remove water phase from group & decrease amplitude if water suppressed
            pkFile.priorKnowledge(pkDict{"Water"}).G_phase = [];
            pkFile.initialValues(pkDict{"Water"}).amplitude = 0.01;
        end

        KeepPeaks = opts.KeepPeaks;
        CombinePeaks = opts.CombinePeaks;
        opts = rmfield(opts, {'KeepPeaks', 'CombinePeaks'});
        
        % Change initial values, bounds, & prior knowledge
        fields = keys(acDict);
        fprintf('\nField\t\tSubField\tPeak\t\tValue\n--------------------------------------------------------\n')
        for i = 1:numel(fields)
            if ~isempty(opts.(fields{i})) % check if IVs/Bs/PKs is empty
                for k = 1:length(opts.(fields{i}))
                    pkValue = opts.(fields{i})(k).Value;

                    % Check for case of chemShift IV change for combined peak
                    if strcmp(opts.(fields{i})(k).Field,"chemShift") && isKey(o2nDict,opts.(fields{i})(k).Peak) && CombinePeaks.(o2nDict(opts.(fields{i})(k).Peak)) % is chemShift, is possible to combine, is combined
                        peaks = n2oDict{o2nDict(opts.(fields{i})(k).Peak)};
                        pname = o2nDict(opts.(fields{i})(k).Peak);
                    else
                        peaks = string(opts.(fields{i})(k).Peak);
                        pname = opts.(fields{i})(k).Peak;
                    end

                    % In case of bounds, ensure value is a 1x2 array
                    if strcmp(fields{i},"bounds") && ~all(size(pkValue) == [1,2])
                        error('Bounds from opts input is not the correct shape.')
                    end

                    % Check for grouped peaks for changing IVs
                    if ~strcmp(opts.(fields{i})(k).Field,"chemShift") && ... only care about non-chemShift
                            isfield(pkFile.priorKnowledge, ['G_' opts.(fields{i})(k).Field]) && ~isempty(pkFile.priorKnowledge(pkDict{pname}).(['G_' opts.(fields{i})(k).Field])) % verify parameter is grouped
                        % Search for other peaks & include for later change
                        for j = 1:numel(pkInd)
                            if ~strcmp(pkInd(j), pname) && ... not current peak
                                    any(~isempty([pkFile.priorKnowledge(pkDict{pkInd(j)}).(['G_' opts.(fields{i})(k).Field])])) && ... is grouped
                                    any([pkFile.priorKnowledge(pkDict{pkInd(j)}).(['G_' opts.(fields{i})(k).Field])]==pkFile.priorKnowledge(pkDict{pname}).(['G_' opts.(fields{i})(k).Field])) % in the same group
                                peaks(end+1) = pkInd(j);
                            end
                        end

                        % Combine peak names for display
                        if numel(peaks) > 1
                            for j = 2:numel(peaks)
                                % Check if peak is being combined
                                if isKey(o2nDict, peaks(j)) && CombinePeaks.(o2nDict(peaks(j)))
                                    pname(end+1) = o2nDict(peaks(j));
                                else
                                    pname(end+1) = peaks(j);
                                end
                            end

                            % Keep only unique names & join with slash
                            pname = join(unique(pname,"stable"),"/");
                        end
                    end

                    % Apply change
                    for subpk = peaks
                        for subsubpk = 1:numel(pkDict{subpk})
                            pkFile.(acDict(fields{i}))(pkDict{subpk}(subsubpk)).(opts.(fields{i})(k).Field) = pkValue;
                        end
                    end

                    % Mark that bounds were set manually
                    if strcmp(fields{i},'Bs')
                        for subpk = peaks
                            for subsubpk = 1:numel(pkDict{subpk})
                                changedBs(pkDict{subpk}(subsubpk)).(opts.Bs(k).Field) = true;
                            end
                        end
                    end

                    % State change
                    fprintf('%-16s%-16s%-16s', acDict(fields{i}), opts.(fields{i})(k).Field, pname)
                    if isempty(pkValue)
                        fprintf('[]')
                    else
                        fprintf('%-.3g',pkValue(1))
                        if ~isscalar(pkValue)
                            for n = 2:numel(pkValue)
                                fprintf(',%.3g',pkValue(n))
                            end
                        end
                    end
                    fprintf('\n')

                end
            end
            opts = rmfield(opts, fields{i});
        end
        fprintf('--------------------------------------------------------\n')
        
        % Mark the changes in bounds
        fields = string(fieldnames(changedBs)).';
        for i = fields(2:end)
            for j = 1:numel(changedBs)
                if ~islogical(changedBs(j).(i))
                    changedBs(j).(i) = false;
                end
            end
        end

        pk_keep = true(length(pkFile.bounds),1); % everything

        % Mark peaks to remove
        pkkeeps = string(fieldnames(opts_def.KeepPeaks)).';
        for pDx = 1:numEntries(pkDict)
            if isfield(KeepPeaks,pkkeeps(pDx))
                pk_keep(pkDict{pkkeeps(pDx)}) = KeepPeaks.(pkkeeps(pDx));
            else
                pk_keep(pkDict{pkkeeps(pDx)}) = opts_def.KeepPeaks.(pkkeeps(pDx));
            end
        end

        Combs = fieldnames(CombinePeaks);
        for i = 1:numel(Combs)
            if CombinePeaks.(Combs{i})
                % Determine new & old names
                name_New = string(Combs{i});
                switch Combs{i}
                    case 'MCL_CH3'
                        name_NewMultiplet = 'MCL_CH3';
                        name_PrimaryOld = "IMCL_CH3";
                        name_SecondaryOld = "EMCL_CH3";
                        name_TertiaryOld = name_SecondaryOld;
                        name_OtherOld = string([]);
                        upCS = true;
                        upLW = false;
                    case 'MCL_CH2'
                        name_NewMultiplet = 'MCL_CH2';
                        name_PrimaryOld = "IMCL_CH2";
                        name_SecondaryOld = "EMCL_CH2";
                        name_TertiaryOld = name_SecondaryOld;
                        name_OtherOld = string([]);
                        upCS = true;
                        upLW = false;
                    case 'tNAD'
                        name_NewMultiplet = {'tNAD1', 'tNAD2'};
                        name_PrimaryOld = "NAD";
                        name_SecondaryOld = "NADH";
                        name_TertiaryOld = name_SecondaryOld;
                        name_OtherOld = string([]);
                        upCS = false;
                        upLW = false;
                    case 'PDE'
                        name_NewMultiplet = {'PDE1', 'PDE2', 'PDE3'};
                        name_PrimaryOld = "GPC";
                        name_SecondaryOld = "GPE";
                        name_TertiaryOld = name_SecondaryOld;
                        name_OtherOld = string([]);
                        upCS = true;
                        upLW = true;
                    case 'Pi'
                        name_NewMultiplet = 'Pi';
                        name_PrimaryOld = "Pia";
                        name_SecondaryOld = "Pib";
                        name_TertiaryOld = name_SecondaryOld;
                        name_OtherOld = string([]);
                        upCS = true;
                        upLW = false;
                    case 'PME'
                        name_NewMultiplet = {'PME1', 'PME2', 'PME3'};
                        name_PrimaryOld = "PC";
                        name_SecondaryOld = "PE";
                        name_TertiaryOld = "G6P";
                        name_OtherOld = ["DPG", "G1P"];
                        upCS = true;
                        upLW = true;
                end

                % Apply changes
                pkFile.bounds(pkDict{name_PrimaryOld}).peakName = name_NewMultiplet;
                pkFile.initialValues(pkDict{name_PrimaryOld}).peakName = name_NewMultiplet;
                pkFile.priorKnowledge(pkDict{name_PrimaryOld}).peakName = name_NewMultiplet;

                if upCS
                    pkFile.bounds(pkDict{name_PrimaryOld}).chemShift(2) = pkFile.bounds(pkDict{name_TertiaryOld}).chemShift(2);
                end
                if upLW
                    pkFile.initialValues(pkDict{name_PrimaryOld}).linewidth = (pkFile.initialValues(pkDict{name_SecondaryOld}).chemShift-pkFile.initialValues(pkDict{name_PrimaryOld}).chemShift)*f + pkFile.initialValues(pkDict{name_PrimaryOld}).linewidth;
                    pkFile.bounds(pkDict{name_PrimaryOld}).linewidth = [round(pkFile.initialValues(pkDict{name_PrimaryOld}).linewidth,-1), pkFile.bounds(pkDict{name_PrimaryOld}).linewidth(2)*2];
                end

                pkDict{name_New} = pkDict{name_PrimaryOld};
                pk_keep(cell2mat(pkDict([name_SecondaryOld, name_TertiaryOld, name_OtherOld]))) = false;
                pkDict([name_PrimaryOld, name_SecondaryOld, name_TertiaryOld, name_OtherOld]) = [];
            end
        end
        
        % Correct group numbers
        kppk = find(pk_keep==true);
        initial_length = length(pk_keep);
        fields2check = convertCharsToStrings(fieldnames(pkFile.priorKnowledge))';
        fields2check = fields2check(contains(fields2check,"G"));
        for pk_rm = initial_length-1:-1:1 % loop in reverse otherwise peaks are missed
            if ~pk_keep(pk_rm)
                pkFile.priorKnowledge(pk_rm).peakName = 'bad';
                for f = fields2check % Loop through prior knowledge
                    g_rm = pkFile.priorKnowledge(pk_rm).(f); % group of bad peak
                    cs = kppk(kppk > pk_rm)';
                    for pk_c = cs % Loop through peaks after bad peak
                        g_c = pkFile.priorKnowledge(pk_c).(f); % group of current peak
                        if ~isempty(g_c) && g_c > pk_rm % current peak is grouped & group is after bad peak
                            if isempty(g_rm) % bad peak is not grouped
                                    pkFile.priorKnowledge(pk_c).(f) = g_c - 1;
                            elseif g_c > g_rm % current group is after bad peak's group
                                    pkFile.priorKnowledge(pk_c).(f) = g_c - 1;
                            end
                        end
                    end
                end
            end
        end
        
        % Remove ppm shift
        if ppm0 ~= 0
            for i = 1:numel(pkFile.initialValues)
                pkFile.initialValues(i).chemShift = pkFile.initialValues(i).chemShift - ppm0;
                pkFile.bounds(i).chemShift = pkFile.bounds(i).chemShift - ppm0;
            end
        end

        % Index which peaks to keep
        pk_idx = 1:initial_length;
        pk_idx = sort(pk_idx(pk_keep));
        
        % Remove peaks
        PK.bounds = pkFile.bounds(pk_idx);
        PK.initialValues = pkFile.initialValues(pk_idx);
        PK.priorKnowledge = pkFile.priorKnowledge(pk_idx);
        PK.svnVersion = pkFile.svnVersion;
        PK.svnHeader = pkFile.svnHeader;
        changedBs = changedBs(pk_idx);
        
        % Reset any groups of 1
        nc = length(PK.bounds); % number of peaks, multiplets combined
        for f = fields2check
            for p = 1:nc
                if ~isempty(PK.priorKnowledge(p).(f))
                    if p == nc && PK.priorKnowledge(p).(f) == p % reset if last peak
                        PK.priorKnowledge(p).(f) = [];
                    else
                        % Check for shared group
                        flag = false;
                        for p2 = 1:nc
                            if p2 ~= p
                                if PK.priorKnowledge(p2).(f) == PK.priorKnowledge(p).(f)
                                    flag = true;
                                    break
                                end
                            end
                        end
        
                        % Remove if no shared
                        if ~flag
                            PK.priorKnowledge(p).(f) = [];
                        end
                    end
                end
            end
        end

        % Check for max linewidth & sigma
        fields2check = {'linewidth','sigma'};
        for i = 1:numel(fields2check)
            if isfield(opts,[fields2check{i} 'Max'])
                for p = 1:nc
                    if PK.bounds(p).(fields2check{i})(2) > opts.([fields2check{i} 'Max'])
                        PK.bounds(p).(fields2check{i})(2) = opts.([fields2check{i} 'Max']);
                        if PK.bounds(p).(fields2check{i})(1) > opts.([fields2check{i} 'Max'])
                            PK.bounds(p).(fields2check{i})(1) = 1;
                        end
                    end
                end
            end
        end

        PK = OXSA.BoundParams(PK);
    end
end