function [results, pth] = loadAMARESresults(varargin)
%LOADAMARESRESULTS Load results .txt file from AMARES quantitation (jMRUI spectroscopy software package)
% 
%   [results, pth] = loadAMARESresults(file)
% 
%   Inputs
%          file: *optional* .txt file, with or without full path
%                   - if full path is not included, file must be in current directory
%        noplot: *optional* presence prevents figure from being displayed
% 
%   Outputs
%       results: structure containing results information & header data
%           pth: source .txt file path
%
% Jacob Degitz, Texas A&M University
% Created 12/17/2024
% Last edited 12/18/2024

% Parse inputs
[pth, npFLAG] = parseInputs(varargin{:});

% Open file
fid = fopen(pth, 'r', 'ieee-be.l64');

%% Extract header information
% File name
fgetl(fid);
fgetl(fid);
fname = fgetl(fid);
fname = fname(11:end);

% Patient name
fgetl(fid);
name = fgetl(fid);
if length(name) == 17
    name = NaN;
else
    name = name(18:end);
end

% Experiment date
date = fgetl(fid);
if length(date) == 20
    date = NaN;
else
    date = date(21:end);
end

% Spectrometer
spectrometer = fgetl(fid);
if length(spectrometer) == 14
    spectrometer = NaN;
else
    spectrometer = spectrometer(15:end);
end

% Additional info
misc = fgetl(fid);
if length(misc) == 24
    misc = NaN;
else
    misc = misc(25:end);
end

% no. points, bandwidth, a0, a1, F0, B0
fgetl(fid);
labs = strsplit(string(fgetl(fid)), '\t');
dats = strsplit(string(fgetl(fid)), '\t');
for i = 1:length(labs)
    switch labs(i)
        case "Points"
            sz = double(dats(i));
        case "Samp.Int."
            dwelltime = double(dats(i))/1000;
            spectralwidth = 1/dwelltime;
        case "ZeroOrder"
            a0 = double(dats(i));
        case "BeginTime"
            a1 = double(dats(i));
        case "Tra.Freq."
            txfrq = double(dats(i));
        case "Magn.F."
            Bo = double(dats(i));
    end
end

% Determine correct nucleus
gamma = floor((txfrq/Bo)*1e-6);
switch gamma
    case 43
        nucleus = '1H';
        gamma = 42.5774780505984; % MHz/T
    case 17
        nucleus = '31P';
        gamma = 17.2514528352478;
    case 10
        nucleus = '13C';
        gamma = 10.7083987615955;
    case 40
        nucleus = '19F';
        gamma = 40.0775824603147;
    case 11
        nucleus = '23Na';
        gamma = 11.2688453499836;
    otherwise
        nucleus = 'unknown';
        gamma = NaN;
end

% Determine ID, if applicable
if ~isnan(double(string(fname(4:6))))
    ID = fname(4:6);
else
    if ~isnan(double(string(fname(4:5))))
        ID = fname(4:5);
    else
        ID = NaN;
    end
end

% Combine into structure
results = struct('ppm', [], ...
                 't', [], ...
                 'specs', [], ...
                 'estimates', [], ...
                 'sz', sz, ...
                 'spectralwidth', spectralwidth, ...
                 'dwelltime', dwelltime, ...
                 'txfrq', txfrq, ...
                 'date', date, ...
                 'dims', struct('t', 1, 'coils', 0, 'averages', 0, 'subSpecs', 0, 'extras', 0), ...
                 'Bo', Bo, ...
                 'nucleus', nucleus, ...
                 'gamma', gamma, ...
                 'fname', fname, ...
                 'ID', ID, ...
                 'name', name, ...
                 'a1', a1, ...
                 'a0', a0, ...
                 'spectrometer', spectrometer, ...
                 'flags', [], ...
                 'misc', misc);

%% Extract actual data
% Metabolite names
fgetl(fid);
fgetl(fid);
fgetl(fid);
mets = strsplit(string(fgetl(fid)), '\t').';
if mets(end) == ""
    mets = mets(1:end-1);
end

% Loop through remaining results values
col = 0;
while true
    fgetl(fid);
    lab = fgetl(fid);

    % Check if data is over
    if strcmp(lab,'Labels : ')
        break
    else
        col = col + 1;
        vals = double(strsplit(string(fgetl(fid)), '\t')).';
        fgetl(fid);
        vals_std = double(strsplit(string(fgetl(fid)), '\t')).';

        % Remove last point
        if isnan(vals(end))
            vals = vals(1:end-1);
            vals_std = vals_std(1:end-1);
        end

        % Convert damping to linewidth
        if strcmp(lab, 'Dampings (Hz)')
            lab = 'Linewidths (Hz)';
            vals = vals/pi;
            vals_std = vals_std/pi;
        end

        % Convert labels to format compatible with tables
        lab = strsplit(string(lab)," (");
        lab = char(lab(1));

        switch lab
            case 'Frequencies'
                    f = vals;
                    fSTD = vals_std;
            case 'Amplitudes'
                    a = vals;
                    aSTD = vals_std;
            case 'Linewidths'
                    l = vals;
                    lSTD = vals_std;
            case 'Phases'
                    p = vals;
                    pSTD = vals_std;
        end
    end
end

fclose(fid);

% Import AMARES estimates
fname_estimates = [pth(1:end-4) '_estimates.mrui'];
try
    [estimates_all, ~] = loadmrui(fname_estimates);
catch
    disp('Estimates file could not be found, please select manually.')
    [estimates_all, ~] = loadmrui();
end

% Import processed spectra
pth_sep = strsplit(pth, '\');
fname = pth_sep{end};
pth_mat = [strjoin(pth_sep(1:end-1),'\') '\' fname(1:end-11) '.mat'];
matdata = load(pth_mat);
results.specs = matdata.Processed.Spectrum;

%% Organize data
% Add imporant info to header
results.ppm = estimates_all(1).ppm;
results.t = estimates_all(1).t;
results.flags = estimates_all(1).flags;

% Remove unnecessary info from estimates
estimates_all = rmfield(estimates_all, {'ppm', 't', 'spectralwidth', 'gamma', 'dwelltime', 'txfrq', 'date', 'dims', 'Bo', 'nucleus', 'fname', 'flags', 'ID', 'name', 'a1', 'a0'});

% Correlate linewidth & peak locations
mloc = zeros(length(estimates_all),1);
for i = 1:length(estimates_all)
    % Extract linewidth
    estimates_all(i).lw = l(i);

    % Obtain peak location
    [~, mloc(i)] = max(real(estimates_all(i).specs));
    estimates_all(i).pk = results.ppm(mloc(i));
end

% Check if FID
mloc = mloc/length(results.ppm);
if all(mloc < 0.1)
    for i = 1:length(estimates_all)
        nstd = mean(estimates_all(i).specs(end-round(length(results.ppm)/100):end));
        estimates_all(i).specs = fftshift(fft(estimates_all(i).specs-nstd));
    end
end

% Check for metabolite names ending in a number
possibles = false(length(estimates_all),1);
for i = 1:length(possibles)
    if ~isnan(str2double(estimates_all(i).metabolite(end)))
        possibles(i) = true;
    end
end

% Combine doublet & triplets
checked = false(length(estimates_all),1);
idxs = cell(1);
i = 0;
for k = 1:length(estimates_all)
    if ~checked(k)
        % Add to main structure
        i = i + 1;
        estimates(i) = estimates_all(k); %#ok<AGROW>
        idxs{i} = k;
        checked(k) = true;

        % Check if higher multiplicity
        if k < length(estimates_all) && length(estimates_all(k).metabolite) > 2
            for m = k+1:length(estimates_all)
                % Check if same name & last char is number or char in both
                combFLAG = false;
                if strcmp(estimates_all(k).metabolite(1:end-1), estimates_all(m).metabolite(1:end-1)) && ...
                   isnan(str2double(estimates_all(k).metabolite(end))) == isnan(str2double(estimates_all(m).metabolite(end)))
                    combFLAG = true;

                    % Longer than two char & all but last two char are same & last two char are both char or number
                elseif length(estimates_all(k).metabolite) > 3 && length(estimates_all(m).metabolite) > 3
                    if strcmp(estimates_all(k).metabolite(1:end-2), estimates_all(m).metabolite(1:end-2)) && ...
                       isnan(str2double(estimates_all(k).metabolite(end))) == isnan(str2double(estimates_all(m).metabolite(end))) && ...
                       isnan(str2double(estimates_all(k).metabolite(end-1))) == isnan(str2double(estimates_all(m).metabolite(end-1)))
                        combFLAG = true;
                    end
                end

                if combFLAG
                    checked(m) = true;

                    % Combine data
                    estimates(i).fids = [estimates(i).fids, estimates_all(m).fids];
                    estimates(i).specs = [estimates(i).specs, estimates_all(m).specs];
                    estimates(i).sz = size(estimates(i).fids);
                    estimates(i).lw = [estimates(i).lw; estimates_all(m).lw];
                    estimates(i).pk = [estimates(i).pk; estimates_all(m).pk];

                    % Trim name
                    estimates(i).metabolite = estimates_all(m).metabolite(1:end-1);

                    % Add index
                    idxs{i} = [idxs{i}, m];
                end
            end
        end
    end
end

% Integrate
A = zeros(length(estimates_all),1);
for i = 1:length(estimates)
    % Determine start & stop points
    estimates(i).lb = estimates(i).pk(1) - estimates(i).lw(1)*2/(results.Bo*results.gamma);
    estimates(i).ub = estimates(i).pk(end) + estimates(i).lw(end)*2/(results.Bo*results.gamma);

    % Integrate within range
    tmp = struct('specs', sum(estimates(i).specs,2), 'ppm', results.ppm);
    estimates(i).A = op_integrate(tmp, estimates(i).lb, estimates(i).ub, 're');
    A(idxs{i}) = estimates(i).A;
end

% Add data to output structure
mloc = zeros(length(idxs),1);
for i = 1:length(idxs)
    results.estimates(i).Name = estimates(i).metabolite;
    results.estimates(i).Spectra = estimates(i).specs;
    results.estimates(i).Frequency = f(idxs{i});
    results.estimates(i).FrequencySTD = fSTD(idxs{i});
    results.estimates(i).Amplitude = a(idxs{i});
    results.estimates(i).AmplitudeSTD = aSTD(idxs{i});
    results.estimates(i).Linewidth = l(idxs{i});
    results.estimates(i).LinewidthSTD = lSTD(idxs{i});
    results.estimates(i).Phase = p(idxs{i});
    results.estimates(i).PhaseSTD = pSTD(idxs{i});
    results.estimates(i).Area = estimates(i).A;
    [~,mloc(i)] = max(real(sum(estimates(i).specs,2)));
end

% % Remove bad peaks (peaks with STD greater than mean, assuming not a subpeak)
% good = true(length(estimates),1);
% for i = 1:length(idxs)
%     if isscalar(idxs{i})
%         if aSTD(idxs{i}) > a(idxs{i}) && lSTD(idxs{i}) > l(idxs{i})
%             good(i) = false;
%         end
%     end
% end
% results.estimates = results.estimates(good);

% Order based on location
locs = zeros(length(results.estimates),1);
for i = 1:length(locs), locs(i) = min(results.estimates(i).Frequency); end
estimates = results.estimates;
for i = 1:length(locs)
    [~,idx] = min(locs); % find min
    results.estimates(i) = estimates(idx); % place min at proper place in structure
    locs(idx) = Inf; % remove min from options
end

% Obtain overall estimate
specs_est = zeros(length(results.estimates), length(results.specs));
for i = 1:length(results.estimates), specs_est(i,:) = sum(results.estimates(i).Spectra,2); end

% Measure residual
residual = results.specs - sum(specs_est,1);

% Select noise region (curently hardcoded as [-30,-20] ppm for 31P & [-3,0] ppm for 1H
switch nucleus
    case '1H'
        idx1 = -3;
        idx2 = 0;
    case '31P'
        idx1 = -30;
        idx2 = -20;
end
[~, idx1] = min(abs(idx1 - results.ppm));
[~, idx2] = min(abs(idx2 - results.ppm));

% Measure variance & fit quality number
var_res = var(residual);
var_noise = var(results.specs(idx1:idx2));
results.FQN = var_res/var_noise;

%% Plot data
if ~npFLAG
    for k = 1:2
        if k == 1
            estdata_ind = real(specs_est);
            estdata_all = real(sum(specs_est,1));
            truedata = real(results.specs);
            resdata = real(results.specs - sum(specs_est,1));
            ylab = 'Amplitude (a.u.)';
            ylims = [min(truedata) max(truedata)];
        else
            estdata_ind = angle(specs_est);
            estdata_all = angle(sum(specs_est,1));
            truedata = angle(results.specs);
            resdata = angle(results.specs - sum(specs_est,1));
            ylab = 'Phase (rad)';
            ylims = [-pi pi];
        end

        figure;
        tiledlayout(3,1,"TileSpacing","tight","Padding","tight");

        % Plot 1
        nexttile
        c = parula(length(results.estimates)+1);
        for i = 1:length(results.estimates)
            legname = results.estimates(i).Name;
            if length(legname) == 4
                if strcmp(legname(end-2:end), 'ATP')
                    switch legname(1)
                        case 'a'
                            legname = '\alpha-ATP';
                        case 'b'
                            legname = '\beta-ATP';
                        case 'g'
                            legname = '\gamma-ATP';
                    end
                end
            end
            plot(results.ppm, estdata_ind(i,:), 'DisplayName', legname, 'Linewidth', 1.5, 'Color', c(i,:));
            if i == 1, hold on; end
        end
        xlim([results.ppm(1) results.ppm(end)])
        ylim(ylims)
        ylabel(ylab)
        set(gca, 'XDir', 'reverse')
        legend

        % Plot 2
        nexttile
        plot(results.ppm, truedata, 'DisplayName', 'Actual', 'Linewidth', 1.5); hold on;
        plot(results.ppm, estdata_all, 'DisplayName', 'AMARES', 'Linewidth', 1.5);
        xlim([results.ppm(1) results.ppm(end)])
        ylim(ylims)
        ylabel(ylab)
        set(gca, 'XDir', 'reverse')
        legend

        % Plot 3
        nexttile
        plot(results.ppm, resdata, 'DisplayName', 'Residual', 'Linewidth', 1.5);
        xlim([results.ppm(1) results.ppm(end)])
        ylim(ylims)
        xlabel('Frequency (ppm)')
        ylabel(ylab)
        set(gca, 'XDir', 'reverse')
        legend
    end
end

%% Functions
    function [pth, npFLAG] = parseInputs(varargin)
        % Set defaults
        npFLAG = false;

        % Check if no inputs
        if nargin == 0
            % Extract path via UI input
            [file, pth] = uigetfile('*.txt');
            pth = [pth file];

            % Check if input is a .mrui file
            if strcmp(pth(end-4:end),'AMARES.txt'), error('Input must be a .txt file'); end
        else
            for j = 1:length(varargin)
                % Check if input is string or character
                if isstring(varargin{j}) || ischar(varargin{j})
                    in = char(varargin{j});
                else
                    error('Input must be a string or character vector')
                end
    
                % Check if no plot flag
                if strcmpi(in, 'NoPlot')
                    npFLAG = true;
                else
                    % Check if full path or just file
                    if strcmp(in(1:2),'C:')
                        pth = in;
                    else
                        % Check if input has file type
                        if any(in == '.')
                            % Check if file is correct type
                            if strcmp(in(end-4:end),'.txt')
                                pth = [pwd '\' in];
                            else
                                error('Input must be a .txt file')
                            end
                        else
                            pth = [pwd '\' in '.txt'];
                        end
                    end
        
                    % Check if file exists
                    if ~exist(pth, 'file'), error('File could not be identified'); end
                end
            end
        end
    end
end