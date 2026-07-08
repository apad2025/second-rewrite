load("C:\Users\apad2\Desktop\Fat_water_separation\DICOM_Files\DD.mat");
datas = {'HydrogenWatSup', 'Hydrogen', 'HydrogenNoise','Phosphorus'};
dataDict = dictionary(["Hydrogen", "HydrogenWatSup"], ["Unsuppressed", "Suppressed"]);

%%
for k = 1:3
    for j = 1:5 % dates
        for i = 1:5 % dogs
            if ~ismissing(DD(i,j).Path.Main) && ~isnan(DD(i,j).(datas{k}))
                if k < 3
                    load(join([DD(i,j).Path.Main "1H" dataDict(datas{k}) string(['MID' num2str(DD(i,j).(datas{k})) '_preproc_twix.mat'])],'\'),'twix');
                    fprintf('\n%sd, %i iter & %i subiter\n', twix.results.freqcorrected.aaDomain, twix.results.freqcorrected.Iter, twix.results.freqcorrected.SubIter)
                    if ~isempty(twix.results.rmbadaverages.badaverages)
                        nbad = numel(twix.results.rmbadaverages.badaverages);
                    else
                        nbad = 0;
                    end
                    fprintf('%i averages removed',nbad)
                end
                dogexplorer(DD(i,j).Name, j, 'DataType',datas{k}, 'Preprocess', 'Override', true);
            end
        end
    end
end
%%
for k = 1:2
    for j = 1:5 % dates
        for i = 1:5 % dogs
            if ~ismissing(DD(i,j).Path.Main) && ~isnan(DD(i,j).(datas{k}))
                if ~exist(join([DD(i,j).Path.Main "1H" dataDict(datas{k}) string(['MID' num2str(DD(i,j).(datas{k})) '_preproc_rwr_twix.mat'])],'\'),'file')
                    if ~strcmp(datas{k}, 'Phosphorus') && ~((i==1 && j>3 && strcmp(datas{k},'HydrogenWatSup')) || (i==3 && j>3 && strcmp(datas{k},'HydrogenWatSup')) || (i==3 && j==5 && strcmp(datas{k},'Hydrogen')))
                        CoilCombination(DD,dataDict,datas,i,j,k)
                    end
                end
            end
        end
    end
end
%%
% Reproc
function Reproc(DD,dataDict,datas,i,j,k)
    if ~strcmp(datas{k}, 'Phosphorus')
        path = join([DD(i,j).Path.Main "1H" dataDict(datas{k})],'\');
    else
        path = join([DD(i,j).Path.Main "31P"],'\');
    end
    path_load = join([path string(['MID' num2str(DD(i,j).(datas{k})) '_raw_twix.mat'])],'\');
    path_save = join([path string(['MID' num2str(DD(i,j).(datas{k})) '_preproc_twix.mat'])],'\');
    path_save_rwr = join([path string(['MID' num2str(DD(i,j).(datas{k})) '_preproc_rwr_twix.mat'])],'\');

    if ~exist(path_save_rwr,'file')
        load(path_save,'twix'); Res = twix.results; twix_rwr = twix;
        save(path_save_rwr,'twix')
        res = struct('leftshifted', struct(), 'filtered', struct(), 'zeropadded', struct(), 'rmbadaverages', struct(), 'freqcorrected', struct(), 'phasecorrected', struct(), 'freqshifted', struct(), 'addedrcvrs', struct());
    
        cd(path)
        cd ..
        twix = loaddat(DD(i,j).(datas{k}));
        spec_comb = sum(sum(twix.specs,3),2);
        if abs(real(spec_comb(end-1))-real(spec_comb(end))) > abs(real(spec_comb(1))-real(spec_comb(end)))*10 && ...
           abs(imag(spec_comb(end-1))-imag(spec_comb(end))) > abs(imag(spec_comb(1))-imag(spec_comb(end)))*10 && ...
           abs(angle(spec_comb(end-1))-angle(spec_comb(end))) > abs(angle(spec_comb(1))-angle(spec_comb(end)))*10 && ...
           abs(abs(spec_comb(end-1))-abs(spec_comb(end))) > abs(abs(spec_comb(1))-abs(spec_comb(end)))*10
            twix.specs = circshift(twix.specs,1,1);
            twix.fids = op_ifft(twix.specs,1);
        end
        cd(path)
        save(path_load, 'twix')
    
        switch twix.nucleus
            case '31P'
                centS = 0;
                ppms.maxN = min(twix.ppm)*5/6;
                ppms.minN = min(twix.ppm);
            case '1H'
                if strcmpi(twix.hdr.Dicom.tScanOptions,'WS')
                    centS = 1.47762;
                else
                    centS = 4.65;
                end
                ppms.maxN = 0;
                ppms.minN = min(twix.ppm) + (twix.spectralwidth/twix.txfrq*1e6)/40;
                if ppms.minN < 0
                    ppms.minN = min(twix.ppm) + (twix.spectralwidth/twix.txfrq*1e6)/40;
                else
                    ppms.minN = min(twix.ppm) - (twix.spectralwidth/twix.txfrq*1e6)/40;
                end
        end
        ppms.minS = centS - (twix.spectralwidth/twix.txfrq*1e6)/20;
        ppms.maxS = centS + (twix.spectralwidth/twix.txfrq*1e6)/20;
    
        % Autophase data
        for cc = 1:twix.sz(twix.dims.coils)
            twix_1ch = IsolateChannel(twix, cc);
            twix_1ch_comb = op_averaging(twix_1ch);
            [~, phShft] = op_autophase(twix_1ch_comb,ppms.minS,ppms.maxS);
            twix_1ch = op_addphase(twix_1ch, phShft);
            twix.flags.phasecorrected(cc) = true;
            res.phasecorrected.phShift(cc) = phShft;
            twix = AppendChannel(twix, twix_1ch, cc);
        end
        clear twix_1ch twix_1ch_comb
    
        % Apply frequency & phase shift
        twix = op_zeropadSpec(twix,2);
        fids_arr = twix.fids;
        t_arr = (0:twix.dwelltime:(length(fids_arr)-1)*twix.dwelltime)';
        for cc = 1:twix.sz(twix.dims.coils)
            if size(Res.freqcorrected.fsCum,2) == 1
                CC = 1;
            else
                CC = cc;
            end
            for av = 1:twix.sz(twix.dims.averages)
                fids_arr(:,cc,av) = addphase(fids_arr(:,cc,av).*exp(1i*t_arr*Res.freqcorrected.fsCum(av,CC)*2*pi),Res.freqcorrected.phsCum(av,CC));
            end
        end
        twix.fids = fids_arr;
        twix.specs = op_fft(twix.fids,twix.dims.t);
        twix = op_zerotrimSpec(twix,twix.sz(twix.dims.t)/2);
        twix.flags.freqcorrected = true;
    
        % Remove bad averages
        mask = zeros(twix.sz(twix.dims.averages),1);
        mask(Res.rmbadaverages.badaverages) = 1;
        goodavg = find(~mask);
        twix.fids = twix.fids(:,:,goodavg,:);
        twix.specs = twix.specs(:,:,goodavg,:);
        twix.sz = size(twix.specs);
        twix.averages = length(goodavg)*twix.rawSubspecs;
    
        % Autophase data again
        for cc = 1:twix.sz(twix.dims.coils)
            twix_1ch = IsolateChannel(twix, cc);
            twix_1ch_comb = op_averaging(twix_1ch);
            [~, phShft] = op_autophase(twix_1ch_comb,ppms.minS,ppms.maxS);
            strct_ph_1ch = op_addphase(twix_1ch, phShft);
            twix.flags.phasecorrected(cc) = true;
            res.phasecorrected.phShift(cc) = res.phasecorrected.phShift(cc)+phShft;
            twix_1ch = strct_ph_1ch;
            twix = AppendChannel(twix, twix_1ch, cc);
        end
    
        % Save results
        twix.results = Res;
        save(path_save, 'twix')
    end
end

% Twix to McMRSData
function Twix2McMRSData(DD,dataDict,datas,i,j,k)
    if ~strcmp(datas{k}, 'Phosphorus')
        path = join([DD(i,j).Path.Main "1H" dataDict(datas{k})],'\');
    else
        path = join([DD(i,j).Path.Main "31P"],'\');
    end
    path_load = join([path string(['MID' num2str(DD(i,j).(datas{k})) '_preproc_twix.mat'])],'\');
    path_load_raw = join([path string(['MID' num2str(DD(i,j).(datas{k})) '_raw_twix.mat'])],'\');
    path_save = join([path string(['MID' num2str(DD(i,j).(datas{k})) '_preproc.mat'])],'\');

    Mcold = load(path_save,'McMRSData'); Mcold = Mcold.McMRSData;
    twixraw = load(path_load_raw,'twix'); twixraw = twixraw.twix;
    load(path_load,'twix');
    if twixraw.dims.coils ~= 0
        twix.rawCoils = twixraw.sz(twixraw.dims.coils);
    end
    McMRSData = FIDa2McMRS(twix, 'Save', false);
    ploc_old = (Mcold.PeakLocation - Mcold.Frequency)/Mcold.Frequency*1e6 + 4.65;
    if any(abs(McMRSData.PeakLocation-ploc_old)>0.01)
        error('Diff peak locations!')
    end

    save(path_save,'McMRSData')
end

% Coil Combination
function CoilCombination(DD,dataDict,datas,i,j,k)
    path = join([DD(i,j).Path.Main "1H" dataDict(datas{k})],'\');
    path_load = join([path string(['MID' num2str(DD(i,j).(datas{k})) '_preproc.mat'])],'\');
    path_noise = join([DD(i,j).Path.Main "1H\Noise" string(['MID' num2str(DD(i,j).HydrogenNoise) '_noiseCov.mat'])],'\');
    path_save = join([path string(['MID' num2str(DD(i,j).(datas{k})) '_proc.mat'])],'\');

    if exist(path_save,'file')
        load(path_save,'McMRSData'); Mcold = McMRSData;
        a0 = McMRSData.Operations.ZeroOrderPhase;
        if McMRSData.Operations.FirstOrderPhase ~= 0
            error('first order phase detected');
        end
    else
        a0 = 0;
    end

    load(path_load,'McMRSData')
    load(path_noise,'Covariance')

    Spectrum_avg = sum(McMRSData.Spectrum,2)./McMRSData.Size.Averages; % test the averaged acquisition, to save time

    % Find indices for signal
    if mean(McMRSData.PeakLocation) > min(McMRSData.PPMAxis) && mean(McMRSData.PeakLocation) < max(McMRSData.PPMAxis)
        Sppm = McMRSData.PeakLocation;
    else
        Sppm = (McMRSData.PeakLocation - McMRSData.Frequency)/McMRSData.Frequency*1e6 + 4.65;
    end
    [~,Sdx] = min(abs(McMRSData.PPMAxis-min(Sppm)));
    if isscalar(Sppm)
        Sdx = Sdx-9:Sdx+10;
    else
        [~,Sdx(2)] = min(abs(McMRSData.PPMAxis-max(Sppm)));
        if Sdx(2) >= Sdx(1)
            Sdx = Sdx(1)-9:Sdx(2)+10;
        else
            Sdx = Sdx(2)-9:Sdx(1)+10;
        end
    end

    % Find indices for noise
    if all(McMRSData.NoiseRegionLimits > min(McMRSData.PPMAxis)) && all(McMRSData.NoiseRegionLimits < max(McMRSData.PPMAxis))
        Nppm = McMRSData.NoiseRegionLimits;
    else
        Nppm = (McMRSData.NoiseRegionLimits - McMRSData.Frequency)/McMRSData.Frequency*1e6 + 4.65;
    end
    
    [~,Ndx] = min(abs(McMRSData.PPMAxis-Nppm(1)));
    [~,Ndx(2)] = min(abs(McMRSData.PPMAxis-Nppm(2)));
    Ndx = Ndx(1):Ndx(2);

    % Calculate noise covariance from data
    if i==4 && j==1
        Covariance_all = zeros([McMRSData.Size.Channels McMRSData.Size.Channels McMRSData.Size.Averages]);
        for a = 1:McMRSData.Size.Averages
            Covariance_all(:,:,a) = cov(squeeze(McMRSData.Spectrum(:,a,Ndx)).');
        end
        Covariance = cov(squeeze(Spectrum_avg(:,:,Ndx)).');
    end

    % generate time vector
    t = reshape((0:McMRSData.Size.Points-1)/McMRSData.SpectralWidth, [1 1 McMRSData.Size.Points]);

    % Loop variables
    lbMax = 40;
    lbMin = 0; % NEVER GO LESS THAN ZERO
    lb = round((lbMax-lbMin)/2); % starting apodization value

    % Metrics for each line broadening value
    SNRs = zeros(lbMax-lbMin+1,1); % SNR
    Qs = zeros(lbMax-lbMin+1,1); % Combination quality

    % Find optimum apodization
    for lbS = [10, 5, 2, 1]
        % Shift lb if at edges
        if lb <= lbMin+lbS*1.5
            lb = lbMin+lbS*2;
        elseif lb >= lbMax-lbS*1.5
            lb = lbMax-lbS*2;
        end

        for LB = (lb-lbS*2):lbS:(lb+lbS*2)
            if SNRs(LB+1) == 0
                Data_with_LB = exp(-pi*LB*t).*proc.IFFT(Spectrum_avg,3);
                [~, Weights_WSVD, Qs(LB+1)]  = proc.WSVD_Combination(Data_with_LB,Covariance);
                Combination = (Weights_WSVD.')*squeeze(Spectrum_avg);

                % Measure SNR
                SNRs(LB+1) = max(real(Combination(Sdx)))/std(real(Combination(Ndx)));
            end
        end

        % Find best apodization value
        [~, bestSNR] = max(SNRs);
        [~, bestQ] = max(Qs);
        if bestQ ~= bestSNR
            lb = round((bestSNR+bestQ)/2)-1;
        else
            lb = bestQ-1;
        end
    end

    Res = [(0:40)', Qs, SNRs];
    idxs = Qs ~= 0;
    Res = Res(idxs,:);
    Combined = zeros([1 McMRSData.Size.Averages McMRSData.Size.Points],'like',McMRSData.Spectrum);

    McMRSData_temp = proc.Apod(McMRSData,lb);
    if i==4 && j==1
        Quality = zeros([1 McMRSData.Size.Averages]);
        for a = 1:McMRSData.Size.Averages
            [~, Weights_WSVD, Quality(a)]  = proc.WSVD_Combination(McMRSData_temp.TimeDomain(:,a,:), Covariance_all(:,:,a));
            Combined(:,a,:) = (Weights_WSVD.')*squeeze(McMRSData.Spectrum(:,a,:));
        end
        Quality = mean(Quality);
    else
        % Quality = zeros([1 McMRSData.Size.Averages]);
        % for a = 1:McMRSData.Size.Averages
        %     [~, Weights_WSVD, Quality(a)]  = proc.WSVD_Combination(McMRSData_temp.TimeDomain(:,a,:), Covariance);
        %     Combined(:,a,:) = (Weights_WSVD.')*squeeze(McMRSData.Spectrum(:,a,:));
        % end

        [~, Weights_WSVD, Quality]  = proc.WSVD_Combination(McMRSData_temp.TimeDomain, Covariance);
        for a = 1:McMRSData.Size.Averages
            Combined(:,a,:) = (Weights_WSVD.')*squeeze(McMRSData.Spectrum(:,a,:));
        end
    end

    % Apply phasing
    Combined = Combined.*exp(1j*a0*pi/180);
    Combined_avg = sum(Combined,2)./McMRSData.Size.Averages;

    % Find new peak location
    if numel(McMRSData.PeakLocation) > 1
        if mean(McMRSData.PeakLocation) > min(McMRSData.PPMAxis) && mean(McMRSData.PeakLocation) < max(McMRSData.PPMAxis)  
            plocs_ppm = unique(McMRSData.PeakLocation);
            if isscalar(plocs_ppm)
                McMRSData.PeakLocation = plocs_ppm;
            else
                plocs = zeros(size(plocs_ppm));
                for p = 1:numel(plocs)
                    [~,plocs(p)] = min(abs(McMRSData.PPMAxis - plocs_ppm(p)));
                end
                plocmin = min(plocs)-round(McMRSData.Size.Points/20);
                if plocmin < 1
                    plocmin = 1;
                end
                plocmax = min(plocs)+round(McMRSData.Size.Points/20);
                if plocmax > McMRSData.Size.Points
                    plocmax = McMRSData.Size.Points;
                end
                [~,plocs] = max(real(squeeze(Combined_avg(:,:,plocmin:plocmax))));
                McMRSData.PeakLocation = McMRSData.PPMAxis(plocs+plocmin-1);
            end
        else
            MHzAxis = linspace(-McMRSData.SpectralWidth/2,McMRSData.SpectralWidth/2,McMRSData.Size.Points)/1e6 + McMRSData.Frequency;
    
            plocs_MHz = unique(McMRSData.PeakLocation);
            if isscalar(plocs_MHz)
                McMRSData.PeakLocation = plocs_MHz;
            else
                plocs = zeros(size(plocs_MHz));
                for p = 1:numel(plocs)
                    [~,plocs(p)] = min(abs(MHzAxis - plocs_MHz(p)));
                end
                plocmin = min(plocs)-round(McMRSData.Size.Points/20);
                if plocmin < 1
                    plocmin = 1;
                end
                plocmax = min(plocs)+round(McMRSData.Size.Points/20);
                if plocmax > McMRSData.Size.Points
                    plocmax = McMRSData.Size.Points;
                end
                [~,plocs] = max(real(squeeze(Combined_avg(:,:,plocmin:plocmax))));
                McMRSData.PeakLocation = MHzAxis(plocs+plocmin-1);
            end
        end
    end

    McMRSData.SourceFiles = McMRSData.Filename;
    McMRSData.Filename = sprintf('MID%i_WSVD + Apod_Quality%0.3f_Acq1to%i',DD(i,j).(datas{k}),Quality,McMRSData.Size.Averages);
    McMRSData.Spectrum = Combined;
    McMRSData.TimeDomain = proc.IFFT(Combined,3);
    McMRSData.Channels = [];
    McMRSData.Size.Channels = 1;
    McMRSData.Operations = struct('SourceFile', McMRSData.Operations,...
        'ChannelCombination',struct('Method','WSVD + Apod','PeakShifting',struct('Occurance',false),'WSVDApodization',lb),...
        'ZeroOrderPhase',a0,...
        'FirstOrderPhase',0,...
        'Apodization',0,...
        'ZeroPadding',0,...
        'FrequencyShift',0,...
        'BaselineCorrection',0);
    if isfield(McMRSData.Operations.SourceFile, 'ChannelCombination')
        McMRSData.Operations.SourceFile = rmfield(McMRSData.Operations.SourceFile, 'ChannelCombination');
    end

    figure; plot(McMRSData.PPMAxis, real(squeeze(Combined_avg))); axis tight; set(gca,'XDir','reverse'); 
    title(sprintf('LB=%iHz, Q=%.3f, LW=%.1f, SNR=%.1f',lb, Quality, meas.LineW(McMRSData), meas.StoN(McMRSData))); xlim([0 10])

    save(path_save,'McMRSData')
end

function [out, offset, SNR, LW] = FixPeakLoc(in)

% Find current peak location
in_avg = op_averaging(in);
[in_filt,~]=op_filter(in_avg,2);
[~, peaks] = max(real(in_filt.specs),[],in_filt.dims.t);

% Find offset from center in Hz
if in.flags.isWaterSuppressed
    centS = 1.47762;
else
    centS = 4.65;
end
[~,centDx] = min(abs(in.ppm-centS));
if numel(centDx) > 1
    centDx = mean(centDx);
end
offset = centDx - peaks;
offset = offset/in.sz(in.dims.t)*in.spectralwidth;
offset = (centS-[2.19676 2.2364 2.20469])*in.txfrq/1e6;
% Shift data
out = in;
for cc = 1:in.sz(in.dims.coils)
    in_1ch = IsolateChannel(in, cc);

    in_1ch = op_freqshift(in_1ch, offset(cc));

    out = AppendChannel(out, in_1ch, cc);
end

if ~isempty(out.results.freqshifted)
    out.results.freqshifted = out.results.freqshifted + offset;
else
    out.results.freqshifted = offset;
end
out_avg = op_averaging(out);

in_comb = op_addrcvrs(in_avg);
out_comb = op_addrcvrs(out_avg);

ylims = [min([min(real(in_comb.specs),[],'all') min(real(out_comb.specs),[],'all') min(real(in_avg.specs),[],'all') min(real(out_avg.specs),[],'all')])
         max([max(real(in_comb.specs),[],'all') max(real(out_comb.specs),[],'all') max(real(in_avg.specs),[],'all') max(real(out_avg.specs),[],'all')])];

colors = distinguishable_colors(3);

ppmrng = in.spectralwidth/in.txfrq*1e6;
[SNR(1),~,~] = op_getSNR(in_comb,centS-ppmrng/20,centS+ppmrng/20,min(in.ppm)+ppmrng/40,0,true);
LW(1) = op_getLW(in_comb,centS-ppmrng/20,centS+ppmrng/20,8,true,true);
[SNR(2),~,~] = op_getSNR(out_comb,centS-ppmrng/20,centS+ppmrng/20,min(in.ppm)+ppmrng/40,0,true);
LW(2) = op_getLW(out_comb,centS-ppmrng/20,centS+ppmrng/20,8,true,true);

figure
tiledlayout(1,2)
nexttile; plot(in_comb.ppm, real(in_comb.specs), 'Color', 'k'); hold on; 
for cc = 1:in.sz(in.dims.coils)
    plot(in_avg.ppm, real(in_avg.specs(:,cc)), 'Color', colors(cc,:)); 
end
axis tight; set(gca,'XDir','reverse'); ylabel('Amplitude (A.u.)'); xlabel('Frequency (ppm)'); ylim(ylims); xlim([0 8]); hold off; title('Original'); subtitle(sprintf('SNR=%.4g, LW=%.3gHz', SNR(1), LW(1)))
nexttile; plot(out_comb.ppm, real(out_comb.specs), 'Color', 'k'); hold on; 
for cc = 1:in.sz(in.dims.coils)
    plot(out_avg.ppm, real(out_avg.specs(:,cc)), 'Color', colors(cc,:)); 
end
axis tight; set(gca,'XDir','reverse'); xlabel('Frequency (ppm)'); ylim(ylims); xlim([0 8]); hold off; title('Shifted'); subtitle(sprintf('SNR=%.4g, LW=%.3gHz', SNR(2), LW(2)))

    function out = IsolateChannel(in, cc)
        % Inputs
        %   in: main data structure
        %   cc: channel number
        %
        % Outputs
        %  out: individual channel data

        out = in;

        % Isolate data
        switch in.dims.coils
            case 1
                out.fids = squeeze(in.fids(cc,:,:));
                out.specs = squeeze(in.specs(cc,:,:));

                % Update temporary size
                out.sz = out.sz([2,3]);
            case 2
                out.fids = squeeze(in.fids(:,cc,:));
                out.specs = squeeze(in.specs(:,cc,:));

                % Update temporary size
                out.sz = out.sz([1,3]);
            case 3
                out.fids = squeeze(in.fids(:,:,cc));
                out.specs = squeeze(in.specs(:,:,cc));

                % Update temporary size
                out.sz = out.sz([1,2]);
        end

        % Update temporary dimensions
        if out.dims.averages ~= 0 && out.dims.averages > out.dims.coils
            out.dims.averages = out.dims.averages-1;
        end
        if out.dims.subSpecs ~= 0 && out.dims.subSpecs > out.dims.coils
            out.dims.subSpecs = out.dims.subSpecs-1;
        end
        if out.dims.extras ~= 0 && out.dims.extras > out.dims.coils
            out.dims.extras = out.dims.extras-1;
        end

        % Set coil dims to zero
        out.dims.coils = 0;
    end

    function out = AppendChannel(out, in, cc)
        % Inputs
        %   in: individual channel data
        %   cc: channel number
        %  out: main data structure

        % Check if first channel
        if cc == 1
            dims = out.dims;

            % Set output structure as input structure
            out = in;
            out.dims = dims;

            % Reshape to 3D array & append data
            switch out.dims.coils
                case 1
                    out.fids = zeros(1, in.sz(1), in.sz(2));
                    out.specs = out.fids;
                case 2
                    out.fids = zeros(in.sz(1), 1, in.sz(2));
                    out.specs = out.fids;
                case 3
                    out.fids = zeros(in.sz(1), in.sz(2), 1);
                    out.specs = out.fids;
            end
        end

        % Append data
        switch out.dims.coils
            case 1
                out.fids(cc,:,:) = in.fids;
                out.specs(cc,:,:) = in.specs;
            case 2
                out.fids(:,cc,:) = in.fids;
                out.specs(:,cc,:) = in.specs;
            case 3
                out.fids(:,:,cc) = in.fids;
                out.specs(:,:,cc) = in.specs;
        end

        % Update size
        out.sz = size(out.fids);
    end
end

% Isolate channel data
function out = IsolateChannel(in, cc)
%{
        Inputs
          in: main data structure
          cc: channel number

        Outputs
         out: individual channel data
%}

out = in;

% Return with no change if all coils are selected
if numel(cc) ~= in.sz(in.dims.coils)
    % Isolate data
    switch in.dims.coils
        case 1
            out.fids = squeeze(in.fids(cc,:,:));
            out.specs = squeeze(in.specs(cc,:,:));

            % Update temporary size
            out.sz = out.sz([2,3]);
        case 2
            out.fids = squeeze(in.fids(:,cc,:));
            out.specs = squeeze(in.specs(:,cc,:));

            % Update temporary size
            out.sz = out.sz([1,3]);
        case 3
            out.fids = squeeze(in.fids(:,:,cc));
            out.specs = squeeze(in.specs(:,:,cc));

            % Update temporary size
            out.sz = out.sz([1,2]);
    end

    % Update temporary dimensions
    if out.dims.averages ~= 0 && out.dims.averages > out.dims.coils
        out.dims.averages = out.dims.averages-1;
    end
    if out.dims.subSpecs ~= 0 && out.dims.subSpecs > out.dims.coils
        out.dims.subSpecs = out.dims.subSpecs-1;
    end
    if out.dims.extras ~= 0 && out.dims.extras > out.dims.coils
        out.dims.extras = out.dims.extras-1;
    end

    % Set coil dims to zero
    out.dims.coils = 0;
end
end

% Append channel data
function out = AppendChannel(out, in, cc)
%{
        Inputs
          in: individual channel data
          cc: channel number
         out: main data structure
%}

% Append data
switch out.dims.coils
    case 1
        out.fids(cc,:,:) = in.fids;
        out.specs(cc,:,:) = in.specs;
    case 2
        out.fids(:,cc,:) = in.fids;
        out.specs(:,cc,:) = in.specs;
    case 3
        out.fids(:,:,cc) = in.fids;
        out.specs(:,:,cc) = in.specs;
end

% Update size
out.sz = size(out.fids);
end
