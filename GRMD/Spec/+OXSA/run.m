function [fitResults, CRBResults, fitStatus, FQN, FQN_global, FQN_SNR, MAD, data, pk] = run(twix, opts)
% Prepare data
[data, exptParams, pk, opts, A2P, P2A, pkInfo, Params, ParamsFit, SNR, changedBs, csbIndices, noisevar, groupMask, group] = OXSA.prep(twix, opts);

% Grab noisevar from TD as well
noisevarTD = var(data.inputFid(exptParams.timeAxis > 0.3*max(exptParams.timeAxis) & exptParams.timeAxis < exptParams.timeAxis(round(exptParams.samples*0.95))));

plotAll = true;
plotPhase = false;

%% Generate initial fit
[data, fitResultsB, CRBResultsB, fitStatusB, pk] = OXSA.runInit(pk, exptParams, data);
[PeakIndices,MDxsA] = rangePeakCorrected(fitResultsB,fitStatusB,data.modelSpecs,pk);
MDxsP = cellfun(@(y)[MDxsA{y}],A2P,UniformOutput=false);
if plotAll
    if ~isfield(opts,'inputFid') %#ok<UNRCH>
        opts.inputFid = data.inputFid;
    end
    PlotResults(fitStatusB, fitResultsB, data, opts, A2P);
    Fs = findall(groot);
    Fs(2).Name = 'OXSA0';
    clear Fs
end

% Calculate base metrics based on bounds
FQN_global_base = CalFQN(data.residual,noisevar);
MAD_global_base = CalMAD(data.residual);
[FQN_base,MAD_base,CoV_base] = EvalFit(data, noisevar, fitResultsB, CRBResultsB, PeakIndices, csbIndices, MDxsP, A2P);
met_base = CalMET(FQN_base,MAD_base,FQN_global_base,MAD_global_base);

% Construct results structure
R = struct('peakName', "", 'initialValues', [], 'fitResults', [], 'FQN', [], 'MAD', [], 'met', [], 'Quality', []);
for p = 1:pkInfo.Np.all
    R(p).peakName = pkInfo.Names.all{p};
    R(p).initialValues.chemShift = mean([pk.initialValues(A2P{p}).chemShift]);
    R(p).initialValues.linewidth = mean([pk.initialValues(A2P{p}).linewidth]);
    R(p).initialValues.sigma = mean([pk.initialValues(A2P{p}).sigma]);
end

%% Fit data
m = struct('s', struct('FQN',Inf,'MAD',Inf,'CoV',Inf,'met',Inf));
m.s(2:11) = m.s(1);
m.b = m.s;
m.p = m.s;
m.p = sDeal(m.p,"Quality",{zeros([1 11])},2,1);
for p = 2:pkInfo.Np.all
    m.p(p,:) = m.p(p-1,:);
end
b_mask = false(11,pkInfo.Np.all);
s_p = ones(1,pkInfo.Np.all,'int8'); % Previous best step
s_b = ones(1,pkInfo.Np.all,'int8')*11; % Current best step
labels = strings([1 pkInfo.Np.all]);

prevUp_met = false(10,pkInfo.Np.all);
prevUp_FQN = false(10,pkInfo.Np.all);
prevUp_MAD = false(10,pkInfo.Np.all);
prevUp_CoV = false(10,pkInfo.Np.all);

b_mask(11,:) = true;

m.p = sDeal(m.p,["FQN","MAD","CoV","met"],{FQN_base,MAD_base,CoV_base,met_base},1,11);

Qmax = 16;

s = 0;
ylims = [Inf -Inf];

% Construct table for presenting results
pklab = '';
for P = 1:pkInfo.Np.all
    pklab = [pklab '%4i']; %#ok<AGROW>
end
bot = sprintf(pklab, 1:pkInfo.Np.all);
if 4*pkInfo.Np.all > 12
    top = pad(' Peak Change',floor(4*pkInfo.Np.all),'both');
elseif 4*pkInfo.Np.all == 12
    top = ' Peak Change';
else
    top = ' Peak Change';
    bot = pad(bot,12,'both');
end
fprintf('KEY\n\t-\tDecreased metric\n\t+\tIncreased metric\n\t~\tPartial metric change\n\tg\tSaving behavior changed due to group\n\tC\tSaving behavior changed due to CoV override\n\tG\tSaving behavior changed due to global metric increase\n\t--\tIncorrect metric calculations due to crazy baseline\n')
fprintf(['Iter ||              Current              ||              Global               ||' top '  ||\n'])
fprintf(['     ||    FQN    |    MAD    |    CoV    ||    FQN    |    MAD    |    CoV    ||' bot '  ||\n'])
while s < 10
% while mean(min(abs(FQN_all-1),[],1)) > 0.5 && s < 10
    s = s + 1;
    initialValues = pk.initialValues;

    % Fit data
    [fitResults, fitStatus, CRBResults, data, MDxsA, PeakIndices] = FitandPlot(data, exptParams, pk, opts, noisevarTD, s, plotAll, A2P);
    switch plotAll
        case true
            ylims = UpdateYlims(ylims);
    end

    % Present warning if fit doesn't occur
    if fitStatus.OUTPUT.iterations == 0
        warning(fitStatus.OUTPUT.message)
    end
    
    % Calculate metrics
    if s > 1
        [FQN_step, MAD_step, CoV_step] = EvalFit(data, noisevar, fitResults, CRBResults, PeakIndices, csbIndices, MDxsP, A2P, m.b(s-1).FQN, m.b(s-1).MAD, dataB);
    else
        [FQN_step, MAD_step, CoV_step] = EvalFit(data, noisevar, fitResults, CRBResults, PeakIndices, csbIndices, MDxsP, A2P);
    end
    m.s(s) = CalAll(data.inputFid-data.modelFid, noisevarTD, fitResults, CRBResults, FQN_global_base, MAD_global_base);
    m.p = sDeal(m.p,["FQN","MAD","CoV"],{FQN_step,MAD_step,CoV_step},1,s);
    met_step = CalMET([m.p(:,s).FQN],[m.p(:,s).MAD],FQN_global_base,MAD_global_base);
    m.p = sDeal(m.p,"met",{met_step},1,s);
    clear FQN_step MAD_step CoV_step met_step
    
    % Append fitting results
    for p = 1:pkInfo.Np.all
        for param = Params
            if isfield(fitResults,param)
                % Check if current parameter can have a group & current peaks are grouped
                if numel(A2P{p}) > 1
                    if isfield(pk.priorKnowledge, join(["G_",param],"")) && ~all(isempty([pk.priorKnowledge(A2P{p}).(join(["G_",param],""))])) && numel([pk.priorKnowledge(A2P{p}).(join(["G_",param],""))])==numel(A2P{p})
                        R(p).fitResults.(param)(s) = fitResults.(param)(MDxsP{p}(1));
                    else
                        R(p).fitResults.(param)(:,s) = [fitResults.(param)([MDxsP{p}])]';
                    end
                else
                    R(p).fitResults.(param)(s) = fitResults.(param)(MDxsP{p}(1));
                end
            end
        end
        R(p).FQN(s) = m.p(p,s).FQN;
        R(p).MAD(s) = m.p(p,s).MAD;
        R(p).CoV(s) = m.p(p,s).CoV;
        R(p).met(s) = m.p(p,s).met;
    end

    %% Evaluate Metrics
    % Initial check to see for any improvement whatsoever
    if s == 1
        prevUp_met(1,:) = [m.p(:,1).met] < [m.p(:,11).met]; % smaller is better
        prevUp_FQN(1,:) = abs([m.p(:,1).FQN]-1) < abs([m.p(:,11).FQN]-1); % closer to unity is better
        prevUp_MAD(1,:) = [m.p(:,1).MAD] < [m.p(:,11).MAD]; % smaller is better
        prevUp_CoV(1,:) = [m.p(:,1).CoV] < [m.p(:,11).CoV]; % smaller is better
    else
        prevUp_met(s,:) = [m.p(:,s).met] < min(reshape([m.p(:,1:s-1).met],[pkInfo.Np.all s-1]),[],2).'; % smaller is better
        prevUp_FQN(s,:) = abs([m.p(:,s).FQN]-1) < min(reshape(abs([m.p(:,1:s-1).FQN]-1),[pkInfo.Np.all s-1]),[],2).'; % closer to unity is better
        prevUp_MAD(s,:) = [m.p(:,s).MAD] < min(reshape([m.p(:,1:s-1).MAD],[pkInfo.Np.all s-1]),[],2).'; % smaller is better
        prevUp_CoV(s,:) = [m.p(:,s).CoV] < min(reshape([m.p(:,1:s-1).CoV],[pkInfo.Np.all s-1]),[],2).'; % smaller is better
    end

    % Set labels
    labels(prevUp_met(s,:)) = "-";
    labels(~prevUp_met(s,:)) = "+";

    p2check = 1:pkInfo.Np.all;
    p2check = p2check(prevUp_met(s,:));

    if s > 1
        %% Baseline Check
        for p = p2check
            % Temporarily calculate FQN using this data
            temp.modelSpecs = dataB.modelSpecs;
            temp.modelSpecs(:,[MDxsP{p}]) = data.modelSpecs(:,[MDxsP{p}]);

            % Check if incorrect calculations due to crazy baseline
            if CalFQN(dataB.inputSpec-sum(temp.modelSpecs,2),noisevar) > CalFQN(dataB.inputSpec-sum(dataB.modelSpecs,2),noisevar)*10 % if global FQN is much much worse, baseline was crazy
                prevUp_FQN(s,p) = false;
                prevUp_MAD(s,p) = false;
                prevUp_CoV(s,p) = false;
                prevUp_met(s,p) = false;
                m.p(p,s).FQN = NaN;
                m.p(p,s).CoV = NaN;
                m.p(p,s).CoV = NaN;
                m.p(p,s).met = NaN;
                labels(p) = "--";

                p2check = p2check(p2check ~= p);
            end
            clear temp
        end

        %% Global Check
        % Temporarily insert all passing peaks into best data and see if global FQN decreases
        fitStatus_tmp = fitStatusB;
        for p = p2check
            [fitStatus_tmp,~] = insertBest(initialValues, initialValuesB, fitStatus, fitStatus_tmp, Params, A2P{p}, MDxsA);
        end

        [~, fitResults_tmp, CRBResults_tmp, fitStatus_tmp] = CalResults(fitStatus_tmp, data.inputFid, exptParams);
        m_tmp = CalAll(fitStatus_tmp.residual, noisevarTD, fitResults_tmp, CRBResults_tmp, FQN_global_base, MAD_global_base);

        % Compare results
        if m_tmp.met < m.b(s-1).met
            % Temporarily mark peaks to skip other checks
            for p = p2check
                labels(p) = "-T";
            end
        else
            % Add each peak to see which ones cause an increase and if adding all others causes a decrease
            pD = 1;
            while pD < numel(p2check)
                % Add peaks that overlap
                pDs = FindOverGroup(pD, p2check, PeakIndices, group);
                fitStatus_tmp = fitStatusB;
                for pD2 = pDs
                    [fitStatus_tmp,~] = insertBest(initialValues, initialValuesB, fitStatus, fitStatus_tmp, Params, A2P{p2check(pD2)}, MDxsA);
                end
                [~, fitResults_tmp, CRBResults_tmp, fitStatus_tmp] = CalResults(fitStatus_tmp, data.inputFid, exptParams);
                m_tmp = CalAll(fitStatus_tmp.residual, noisevarTD, fitResults_tmp, CRBResults_tmp, FQN_global_base, MAD_global_base);

                % Check if better without current peak
                if m_tmp.met > m.b(s-1).met
                    prevUp_met(s,p2check(pD):p2check(pD2)) = false;
                    labels(p2check(pD):p2check(pD2)) = labels(p2check(pD):p2check(pD2)) + "G";

                    pDflag = true(size(p2check));
                    pDflag(pD:pD2) = false;
                    p2check = p2check(pDflag);
                    clear pDflag
                else
                    pD = pD2+1;
                end
            end
        end
        clear m_tmp m_tmp1 m_tmp2 p2 fitStatus_tmp

        %% Minor Change Check
        for p = p2check
            if ~prevUp_CoV(s,p) && ~strcmp(labels(p),"-T") % only important for peaks with improved metrics but worsened CRBs
                % Check to see if improvement is insignificant & causes increased CRB
                if ((abs(m.p(p,s).met-m.p(p,s_b(p)).met)/m.p(p,s_b(p)).met)/(abs(m.p(p,s).CoV-m.p(p,s_b(p)).CoV)/m.p(p,s_b(p)).CoV)<1 && m.p(p,s).CoV>1) ...
                        || m.p(p,s).CoV>100
                    prevUp_met(s,p) = false;
                    labels(p) =  labels(p) + "C";
                    p2check = p2check(p2check ~= p);

                else % Check to see if change is insignificant & causes increased CRB
                    % Determine percent change in parameters
                    tmpResults = R(p).fitResults;
                    pPerc = zeros([1 numel(Params)-1]);
                    for param = 1:numel(Params)-1
                        if strcmp(Params(param),'linewidth') % Measure net change in combined Lorentzian & Gaussian linewidth
                            pPerc(param) = abs((tmpResults.linewidth(s)-tmpResults.linewidth(s-1))/((tmpResults.linewidth(s)+tmpResults.linewidth(s-1))/2) + ...
                                (tmpResults.sigma(s)-tmpResults.sigma(s-1))/((tmpResults.sigma(s)+tmpResults.sigma(s-1))/2));
                        else
                            pPerc(param) = abs((tmpResults.(Params(param))(s)-tmpResults.(Params(param))(s-1))/((tmpResults.(Params(param))(s)+tmpResults.(Params(param))(s-1))/2));
                        end
                    end
                    if sum(pPerc) < 0.05
                        prevUp_met(s,p) = false;
                        labels(p) =  labels(p) + "C";
                        p2check = p2check(p2check ~= p);
                    end
                end
            end
        end
    end

    % Label mismatches
    for p = 1:pkInfo.Np.all
        if prevUp_met(s,p) && (~prevUp_FQN(s,p) || ~prevUp_MAD(s,p))
            labels(p) = "~" + labels(p); % partial improvement

            % Check if CoV increases or decreases
            if ~prevUp_CoV(s,p) && ~strcmp(labels(p),"-T")
                prevUp_met(s,p) = false;
                labels(p) =  labels(p) + "C";
            end
        elseif ~prevUp_met(s,p) && (prevUp_FQN(s,p) || prevUp_MAD(s,p)) && ~strcmp(labels(p),"-G")
            labels(p) = "~" + labels(p); % partial unimprovement

            % Do quick global check
            if s > 1
                % Add peaks that overlap
                pDs = FindOverGroup(p, 1:pkInfo.Np.all, PeakIndices, group);
                fitStatus_tmp = fitStatusB;
                for pD2 = pDs
                    [fitStatus_tmp,~] = insertBest(initialValues, initialValuesB, fitStatus, fitStatus_tmp, Params, A2P{pD2}, MDxsA);
                end
                [~, fitResults_tmp, CRBResults_tmp, fitStatus_tmp] = CalResults(fitStatus_tmp, data.inputFid, exptParams);
                m_tmp = CalAll(fitStatus_tmp.residual, noisevarTD, fitResults_tmp, CRBResults_tmp, FQN_global_base, MAD_global_base);

                % Check if better with current peak
                if m_tmp.met < m.b(s-1).met
                    prevUp_met(s,p) = true;
                    labels(p) =  "~+G";
                end
            end
        end
    end

    % Remove temporary flag
    labels = strrep(labels, "-T", "-");

    %% Evaluate Groups
    for p = 1:pkInfo.Np.all
        % Now, create temporary xFit if prevUp mismatch
        if sum(groupMask(p,:)) > 1 && (~all(prevUp_met(s,groupMask(p,:))) || all(~prevUp_met(s,groupMask(p,:))))
            if s == 1
                if mean([m.p(groupMask(p,:),s).met]) < mean([met_base(groupMask(p,:))]) && mean(abs([m.p(groupMask(p,:),s).FQN]-1)) < mean(abs([FQN_base(groupMask(p,:))]-1))
                    prevUp_met(s,groupMask(p,:)) = true;
                    labels(groupMask(p,:)) = labels(groupMask(p,:)) + "g";
                else
                    prevUp_met(s,groupMask(p,:)) = false;
                    labels(groupMask(p,:)) = labels(groupMask(p,:)) + "g";
                end
            else
                data_tmp = struct('inputFid',data.inputFid,'inputSpec',data.inputSpec);
                fitStatus_tmp = fitStatusB;

                for peak2 = [MDxsA{[A2P{groupMask(p,:)}]}]
                    for fs = fieldnames(fitStatus.constraintsCellArray)'
                         xFitDx = findxFitInd(fitStatus.constraintsCellArray.(fs{:}){peak2});
                         if ~isempty(xFitDx)
                             fitStatus_tmp.xFit(xFitDx) = fitStatus.xFit(xFitDx);
                         end
                    end
                end

                % Check if grouped peaks don't have improvement
                [~,~,data_tmp.modelFids] = AMARES.makeModelFidAndJacobianReIm(fitStatus_tmp.xFit,fitStatus_tmp.constraintsCellArray,exptParams.beginTime,exptParams.dwellTime,exptParams.imagingFrequency,exptParams.samples, 'complexOutput', true);
                data_tmp.modelSpecs = specFft(data_tmp.modelFids,1);
                data_tmp.residual = data_tmp.inputSpec - sum(data_tmp.modelSpecs,2);

                % If not all are improved, take average of grouped parameters
                if ~all(CalMET(CalFQN(data_tmp.residual,noisevar),CalMAD(data_tmp.residual),FQN_global_base,MAD_global_base) <= m.b(s).met)
                    for param = ["linewidth", "amplitude", "phase", "sigma"]
                        xFitDx = findxFitInd(fitStatus.constraintsCellArray.(param){p});

                        if ~strcmp(fitStatus.constraintsCellArray.(param){p}{1}, '@(a)a;') % Ignore fixed parameters
                            if isfield(fitStatus.constraintsCellArray,param) && (~isfield(pk.priorKnowledge,join(["G_" param],"")) || (~isempty(pk.priorKnowledge(p).(join(["G_" param],""))) && any(pk.priorKnowledge(p).(join(["G_" param],"")) == group{p})))
                                if strcmp(param,"phase")
                                    for subP = 1:numel(xFitDx)
                                        if abs(fitStatus.xFit(xFitDx(subP))-fitStatusB.xFit(xFitDx(subP))) > 180 % taking the mean of phase is more complicated when phase wraps
                                            if fitStatus.xFit(xFitDx(subP)) > fitStatusB.xFit(xFitDx(subP))
                                                fitStatus.xFit(xFitDx(subP)) = mean([fitStatus.xFit(xFitDx(subP))-360 fitStatusB.xFit(xFitDx(subP))]);
                                            else
                                                fitStatus.xFit(xFitDx(subP)) = mean([fitStatus.xFit(xFitDx(subP)) fitStatusB.xFit(xFitDx(subP))-360]);
                                            end
                                        end
                                    end
                                else
                                    fitStatus.xFit(xFitDx) = mean([fitStatus.xFit(xFitDx); fitStatusB.xFit(xFitDx)],1);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    clear fitStatus_tmp fitResults_tmp data_tmp PeakIndices_tmp csbIndices_tmp FQN_tmp MAD_tmp met_tmp fs peak2

    %% Evaluate overall fit quality
    Quality = (prevUp_FQN(s,:) + prevUp_MAD(s,:))*Qmax/2;
    for p = 1:pkInfo.Np.all
        if Quality(p) > 0
            if ~prevUp_met(s,p) % one metric slightly better, one slightly worse, but overall worse
                Quality(p) = 0;
            elseif all(PercDiff(R(p).fitResults.chemShift(s),fitResultsB.chemShift(p)) < 0.05) && ... very little change
                    all(PercDiff(R(p).fitResults.amplitude(s),fitResultsB.amplitude(p)) < 0.05) && ...
                    all(PercDiff(R(p).fitResults.linewidth(s),fitResultsB.linewidth(p)) < 0.05) && ...
                    all(PercDiff(R(p).fitResults.phase(s),fitResultsB.phase(p)) < 0.05) && ...
                    isfield(fitResults,'sigma') && all(PercDiff(R(p).fitResults.sigma(s),fitResultsB.sigma(p)) < 0.05)
                Quality(p) = Quality(p)/(Qmax/4);
            elseif s > 1 && ...
                    R(p).initialValues.chemShift(s) == R(p).initialValues.chemShift(s_b(p)) && ... % no change from inputs
                    R(p).initialValues.linewidth(s) == R(p).initialValues.linewidth(s_b(p)) && ...
                    R(p).initialValues.sigma(s) == R(p).initialValues.sigma(s_b(p))
                Quality(p) = Quality(p)/(Qmax/2);
            end
        else % check if results are really bad
            if m.p(p,s).FQN > FQN_base(p) % check if step is worse than input peak
                Quality(p) = Quality(p)-(Qmax/4);
                if FQN_base(p) > FQN_global_base
                    if m.p(p,s).FQN > FQN_base(p)*1.5 % if input peak is greater than input data, check if step is over 1.5x bad than input peak
                        Quality(p) = Quality(p)-(Qmax/4);
                    end
                elseif m.p(p,s).FQN > FQN_global_base % check if step is worse than entire input data (only use when input peak is less than input data)
                    Quality(p) = Quality(p)-(Qmax/4);
                end
            elseif FQN_base(p) > FQN_global_base && m.p(p,s).FQN > FQN_global_base % if input peak is greater than input data, check if step is worse than entire input data
                Quality(p) = Quality(p)-(Qmax/8);
            end
            if m.p(p,s).MAD > MAD_base(p)
                Quality(p) = Quality(p)-(Qmax/4);
                if MAD_base(p) > MAD_global_base
                    if m.p(p,s).MAD > MAD_base(p)*1.5
                        Quality(p) = Quality(p)-(Qmax/4);
                    end
                elseif m.p(p,s).MAD > MAD_global_base
                    Quality(p) = Quality(p)-(Qmax/4);
                end
            elseif MAD_base(p) > MAD_global_base && m.p(p,s).MAD > MAD_global_base
                Quality(p) = Quality(p)-(Qmax/8);
            end
        end
        R(p).Quality(s) = Quality(p);
    end
    m.p = sDeal(m.p,"Quality",{Quality},1,s);

    %% Save Best Results
    if s == 1
        fitStatusB = fitStatus;
        dataB = data;
        initialValuesB = pk.initialValues;
        s_p(:) = ones(1,pkInfo.Np.all,'int8');
        s_b(:) = s;
        b_mask(1,:) = true;

        % Trim metric arrays
        b_mask = b_mask(1:10,:);
        m.s = m.s(1:10);
        m.b = m.b(1:10);
        m.p = m.p(:,1:10);
    else
        for p = 1:pkInfo.Np.all
            if prevUp_met(s,p)
                s_p(p) = s_b(p);
                s_b(p) = s;
                b_mask(:,p) = false;
                b_mask(s,p) = true;

                [fitStatusB,initialValuesB] = insertBest(initialValues, initialValuesB, fitStatus, fitStatusB, Params, A2P{p}, MDxsA);
            end
        end
        clear p
    end

    %% Calculate new best data
    [dataB, fitResultsB, CRBResultsB, fitStatusB] = CalResults(fitStatusB, data.inputFid, exptParams);
    m.b(s) = CalAll(dataB.inputFid-dataB.modelFid, noisevarTD, fitResultsB, CRBResultsB, FQN_global_base, MAD_global_base);
    fprintf([' ||%8.2f   |%10.2e |%8.2f   ||' strrep(pklab,'i','s') '  ||\n'], m.b(s).FQN, m.b(s).MAD, m.b(s).CoV, labels)

    %% Update initial conditions
    % Determine what to do next
    if s>4 && all(~prevUp_met(s-2:s,:),'all')
        disp('Fit has stopped improving.')
        break
    elseif s < 10
        % Update fitting parameters
        pk.initialValues = initialValuesB; % set initial values as best initial values
        for p = 1:pkInfo.Np.all
            % Select which data to selectively override initial values
            Quality = R(p).Quality(s);
            if Quality == Qmax
                tmpfitResults = A2S(fitResultsB); % best fit results
            elseif Quality > 0
                tmpfitResults = A2S(fitResults); % recent fit results
            elseif Quality <= 0
                tmpfitResults = A2S(fitResults); % recent fit results
                tmpfitResults2 = A2S(fitResultsB); % for step 1, this is equal to IVs
            end

            % Override intial values with fit results
            for pD = 1:numel(A2P{p})
                tmpfitResults_pk = tmpfitResults(MDxsA{A2P{p}(pD)}); % fit results for the current peak
                for param = ParamsFit

                    % Extract result
                    switch param
                        case "chemShift"
                            if isscalar(tmpfitResults_pk)
                                tmpfitResult = tmpfitResults_pk.(param);
                            else
                                switch pkInfo.multiplets.all(p)
                                    case 1
                                        tmpfitResult = tmpfitResults_pk.(param);
                                    case 2
                                        tmpfitResult = (tmpfitResults_pk(1).(param) + tmpfitResults_pk(2).(param))/2;
                                    case 3
                                        tmpfitResult = tmpfitResults_pk(2).(param);
                                end
                            end
                        otherwise
                            if isfield(pk.initialValues, param)
                                if any([tmpfitResults(MDxsA{A2P{p}(1)}).(param)]~=[tmpfitResults(MDxsA{A2P{p}(pD)}).(param)])
                                    error('claude was right')
                                end
                                tmpfitResult = tmpfitResults(MDxsA{A2P{p}(pD)}).(param);
                            end
                    end

                    % Only partial change, so average update with best values
                    if Quality < Qmax && Quality > 0
                        pk.initialValues(A2P{p}(pD)).(param) = (initialValuesB(A2P{p}(pD)).(param) + tmpfitResult)/2;
                    elseif Quality <= 0
                        % Ensure the parameter is not fixed
                        if ~isempty(pk.bounds(A2P{p}(pD)).(param))
                            % Determine which direction data needs to be shifted
                            if mean([tmpfitResults2(MDxsA{A2P{p}(pD)}).(param)]) - tmpfitResult > 0
                                idx = 2;
                            else
                                idx = 1;
                            end

                            if Quality == 0
                                pk.initialValues(A2P{p}(pD)).(param) = (initialValuesB(A2P{p}(pD)).(param)*3/2 + (pk.bounds(A2P{p}(pD)).(param)(idx) + tmpfitResult)/4)/2;
                            elseif Quality > -Qmax
                                pk.initialValues(A2P{p}(pD)).(param) = (initialValuesB(A2P{p}(pD)).(param) + (pk.bounds(A2P{p}(pD)).(param)(idx) + tmpfitResult)/2)/2;
                            else
                                pk.initialValues(A2P{p}(pD)).(param) = (initialValuesB(A2P{p}(pD)).(param) + pk.bounds(A2P{p}(pD)).(param)(idx) + tmpfitResult)/3;
                            end
                        end
                    end
                end
            end
        end
        clear ampR tmpfitResults pD tmpfitResults_pk modvar

        if s == 1
            % Determine which peaks can be used for LW mean & std
            if strcmp(opts.nucleus, '1H')
                notsameLW = {'Lip20', 'Lip22', 'Cr3'};
                if ~twix.flags.isWaterSuppressed
                    notsameLW{end+1} = 'Water'; %#ok<AGROW>
                end
            elseif strcmp(opts.nucleus, '31P')
                notsameLW = {'GPC','GPE','PDE','MP','Unknown','PC','PE','PME','PPA'};
            end
            sameLW = 1:length(fitResults.chemShift);
            sameLWc = 1:pkInfo.Np.pk;
            for p = 1:pkInfo.Np.pk
                if any(strcmp(pkInfo.Names.pk{p}, notsameLW)) % all but these peaks should have roughly the same linewidth
                    sameLW(MDxsA{P2A(p)}) = NaN;
                    sameLWc(p) = NaN;
                end
            end
            sameLW = sameLW(~isnan(sameLW));
            sameLWc = sameLWc(~isnan(sameLWc));
        end

        % Determine mean & std linewidth
        LWmean = mean(fitResultsB.linewidth(sameLW));
        LWstd = std(fitResultsB.linewidth(sameLW));
        
        % Update linewidth bounds
        for p = 1:pkInfo.Np.pk
            if ~changedBs(p).linewidth % Only change peaks with no manual linewidth
                if any(p == sameLWc)
                    pk.bounds(p).linewidth = LWmean + [-LWstd, LWstd]*1.5;
                else
                    pk.bounds(p).linewidth = LWmean + [-LWstd, LWstd]*2.5;
                end

                % Verify no negative
                if pk.bounds(p).linewidth(1) <= 1
                    pk.bounds(p).linewidth(1) = 1;
                end
                if isfield(opts,'linewidthMin') && pk.bounds(p).linewidth(1) < opts.linewidthMin
                    pk.bounds(p).linewidth(1) = opts.linewidthMin;
                end
                if isfield(opts,'linewidthMax') && pk.bounds(p).linewidth(2) > opts.linewidthMax
                    pk.bounds(p).linewidth(2) = opts.linewidthMax;
                end
            end
        end

        % Update sigma bounds as well
        SGbounds = [1, LWstd*3];
        if SGbounds(2) > opts.sigmaMax
            SGbounds(2) = opts.sigmaMax;
        end
        for p = 1:pkInfo.Np.pk
            if ~changedBs(p).sigma
                pk.bounds(p).sigma = SGbounds;
            end
        end

        % Ensure all parameters are within bounds
        pk = OXSA.BoundParams(pk);

        % Append to R structure
        for p = 1:pkInfo.Np.all
            R(p).initialValues.chemShift(s+1) = mean([pk.initialValues(A2P{p}).chemShift]);
            R(p).initialValues.linewidth(s+1) = mean([pk.initialValues(A2P{p}).linewidth]);
            R(p).initialValues.sigma(s+1) = mean([pk.initialValues(A2P{p}).sigma]);
        end
    end
end
m.s = m.s(1:s);
m.b = m.b(1:s);
m.p = m.p(:,1:s);
prevUp_met = prevUp_met(1:s,:); %#ok<NASGU>
prevUp_FQN = prevUp_FQN(1:s,:); %#ok<NASGU>
prevUp_MAD = prevUp_MAD(1:s,:); %#ok<NASGU>
prevUp_CoV = prevUp_CoV(1:s,:); %#ok<NASGU>


clear MAD_prev FQN_prev fs LWbounds LWmean LWstd sameLW sameLWc notsameLW initialValues

% Save best results
data = dataB;
b_mask = b_mask(1:s,:);
fitStatus = fitStatusB;
fitResults = fitResultsB;
pk.initialValues = initialValuesB;
if all(fitStatus.EXITFLAG == fitStatus.EXITFLAG(1))
    fitStatus.EXITFLAG = fitStatus.EXITFLAG(1);
end

% Calculate SD
CRBResults = AMARES.estimateCRB(exptParams.imagingFrequency, exptParams.dwellTime, exptParams.beginTime, fitStatus.noise_var, fitStatus.xFit, fitStatus.constraintsCellArray);
fitResults.covariance = CRBResults.covariance;
CRBResults = rmfield(CRBResults, 'covariance');

[FQN, MAD, CoV] = EvalFit(data, noisevar, fitResults, CRBResults, PeakIndices, csbIndices, MDxsP, A2P, m.b(s-1).FQN, m.b(s-1).MAD, dataB);

clear dataB fitStatusB fitResultsB initialValuesB CRBResultsB

%% Display results
% Plot final
if s > 1 || ~plotAll
    hFig = PlotResults(fitStatus, fitResults, data, opts, A2P);
    UpdateYlims(ylims);
end

if s > 1 || ~plotAll
    Results = formOutput(fitResults, CRBResults, fitStatus, exptParams, hFig);
else
    Results = formOutput(fitResults, CRBResults, fitStatus, exptParams);
end
results = AMARES.sortFitData([],Results,exptParams,pk,exptParams.beginTime,'fitOptions',opts);

if ~isfield(exptParams,'offset')
    exptParams.offset = 0;
else
    for p = 1:pkInfo.Np.all
        R(p).initialValues.chemShift = R(p).initialValues.chemShift - exptParams.offset;
        R(p).fitResults.chemShift = R(p).fitResults.chemShift - exptParams.offset;
    end
end
switch plotPhase
    case true
        [PDxs, ~, ~, combSpecs] = AMARES.rangePeak(fitResults, fitStatus, data.modelSpecs);
        indivColors = distinguishable_colors(numel(fitStatus.pkWithLinLsq.bounds));
        figure(Theme='light')
        tiledlayout(2,1, 'TileSpacing', 'compact')
        nexttile; hold on;
        plot(exptParams.ppmAxis-exptParams.offset, angle(data.inputSpec), 'color', 'k'); 
        plot(exptParams.ppmAxis-exptParams.offset, angle(data.modelSpec), 'color', 'r'); 
        axis tight; set(gca, 'XDir', 'reverse'); xlim(opts.xlims); set(gca, 'XTickLabel', []); ylabel('Spectrum fit'); box off; hold off;
        nexttile; hold on
        plot(exptParams.ppmAxis-exptParams.offset, angle(data.inputSpec), 'color', 'k');
        for px = 1:size(combSpecs,2)
            plot(exptParams.ppmAxis(PDxs{px})-exptParams.offset, angle(combSpecs(PDxs{px},px)), 'color', indivColors(px,:)); 
        end
        set(gca, 'XDir', 'reverse'); xlim(opts.xlims); ylim([-pi pi]); ylabel('Individual Peaks'); xlabel('\delta / ppm'); box off; hold off;
end

% Correct for 1H zero ppm
if strcmp(opts.nucleus,'1H')
    fitResults.chemShift = fitResults.chemShift - exptParams.offset;
    for p = 1:pkInfo.Np.pk
        pk.bounds(p).chemShift = pk.bounds(p).chemShift - exptParams.offset;
        pk.initialValues(p).chemShift = pk.initialValues(p).chemShift - exptParams.offset;
    end
end

% Calculate FQN/SNR
FQN_global = CalFQN(data.inputFid-data.modelFid,noisevarTD);
FQN_SNR = FQN_global/SNR;
CoV = CalCOV(fitResults, CRBResults, MDxsA);
if FQN_SNR < 1
    fprintf('FQN/SNR = %.5f\n', FQN_SNR);
elseif FQN/SNR < 10
    fprintf('FQN/SNR = %.3f\n', FQN_SNR);
else
    fprintf('FQN/SNR = %.1f\n', FQN_SNR);
end
fprintf('CoV = %.1f+/-%.1f\n', mean(CoV), std(CoV))

%% Check that fitted parameters did not hit bound limits
% This is to check that the fit was not limited unduly by the prior knowledge. The flag results in the fitted parameters being displayed in red in the eventual excel/pdf record.

% Preallocate
for param = Params
    fitResults.boundFlag.(param) = false([1 sum(pkInfo.multiplets.all)]);
end

for p = 1:pkInfo.Np.all
    for param = Params
        for mult = 1:pkInfo.multiplets.all(p)
            if ~isempty(pk.bounds(p).(param)) && ... bounds are not empty
                    (isempty(fitResults.boundFlag.(param)(MDxsP{p}(mult))) || ~fitResults.boundFlag.(param)(MDxsP{p}(mult))) % current parameter has not been checked/is false so far
                switch param
                    case 'chemShift'
                        if fitResults.(param)(MDxsP{p}(mult)) > pk.bounds(p).(param)(2)*1.05 || fitResults.(param)(MDxsP{p}(mult)) == pk.bounds(p).(param)(2)
                            fitResults.boundFlag.(param)(MDxsP{p}(mult)) = true;
                        end
                    case 'amplitude'
                        if fitResults.(param)(MDxsP{p}(mult)) > max(pk.bounds(p).(param).*[0.95 1.05]) || ...
                                fitResults.(param)(MDxsP{p}(mult)) < min(pk.bounds(p).(param).*[0.95 1.05]) || ...
                                any(abs((fitResults.(param)(MDxsP{p}(mult)) - pk.bounds(p).(param))./pk.bounds(p).(param)) < 1e-5) % difference is less than 1e-3%
                            fitResults.boundFlag.(param)(MDxsP{p}(mult)) = true;
                        end
                    otherwise
                        if fitResults.(param)(MDxsP{p}(mult)) > pk.bounds(p).(param)(2)*1.05 || ...
                                fitResults.(param)(MDxsP{p}(mult)) < pk.bounds(p).(param)(1)*0.95 || ...
                                any(abs((fitResults.(param)(MDxsP{p}(mult)) - pk.bounds(p).(param))./pk.bounds(p).(param)) < 1e-5) % difference is less than 1e-3%
                            fitResults.boundFlag.(param)(MDxsP{p}(mult)) = true;
                        end
                end
            end
        end
    end
end

%% Functions

    % Grab overlapping/grouped indices
    function pDxs = FindOverGroup(pDx, p2check, peakIndices, group)
        pDxs = pDx;
        pDx2 = pDx+1;

        % Add peaks that overlap or are grouped
        while pDx2<=numel(p2check) && ... less than end
                (peakIndices{p2check(pDxs(end))}(end)>peakIndices{p2check(pDx2)}(1) || any(group{pDx}==group{pDx2})) % overlap or grouped
            pDxs = [pDxs pDx2]; %#ok<AGROW>
            pDx2 = pDx2 + 1;
        end
    end

    % Calculate results from fitStatus
    function [data, fitResults, CRBResults, fitStatus] = CalResults(fitStatus, inputFid, exptParams)
        data.inputFid = inputFid;
        data.inputSpec = specFft(inputFid);

        fitResults = AMARES.applyModelConstraints(fitStatus.xFit,fitStatus.constraintsCellArray);
        [data.modelSpec,data.modelSpecs,data.modelFid,data.modelFids] = AMARES.makeModelSpec(fitStatus,struct('firstOrder',false));
        fitStatus.residual = data.inputFid - data.modelFid;
        fitStatus.noise_var = var(real(fitStatus.residual)); % from amaresFit
        data.residual = data.inputSpec - data.modelSpec;
        CRBResults = AMARES.estimateCRB(exptParams.imagingFrequency, exptParams.dwellTime, exptParams.beginTime, fitStatus.noise_var, fitStatus.xFit, fitStatus.constraintsCellArray);
    end

    % Insert best results into fitStatus
    function [fitStatusB,initialValuesB] = insertBest(initialValues, initialValuesB, fitStatus, fitStatusB, Params, A2P, MDxsA)
        for par = Params
            for subpks = [A2P] %#ok<NBRAK2>
                if isfield(initialValues, par)
                    initialValuesB(subpks).(par) = initialValues(subpks).(par);
                end
                for subsubpks = [MDxsA{subpks}]
                    xFitdx = findxFitInd(fitStatus.constraintsCellArray.(par){subsubpks});
                    if ~isempty(xFitdx)
                        fitStatusB.xFit(xFitdx) = fitStatus.xFit(xFitdx);
                    end
                    %%%%%%% THIS SECTION MAY NOT BE NECESSARY %%%%%%%
                    if any(cellfun(@(x,y) any(x ~= y), fitStatus.constraintsCellArray.(par){subsubpks}, fitStatusB.constraintsCellArray.(par){subsubpks})) % check if any of the constraint equation has changed
                        warning('this section should be kept')
                    end
                    fitStatusB.constraintsCellArray.(par){subsubpks} = fitStatus.constraintsCellArray.(par){subsubpks};
                    %%%%%%% THIS SECTION MAY NOT BE NECESSARY %%%%%%%
                end
            end
        end

        fitStatusB.EXITFLAG(A2P) = fitStatus.EXITFLAG;
        fitStatusB.OUTPUT(A2P) = fitStatus.OUTPUT;
        fitStatusB.relativeNorm = fitStatus.relativeNorm;
        fitStatusB.resNormSq = fitStatus.resNormSq;
        fitStatusB.pkWithLinLsq.initialValues(A2P) = fitStatus.pkWithLinLsq.initialValues(A2P);
    end

    % Deal out values to a structure
    function S = sDeal(S,field,value,dim,I)
        switch dim % sub-for loops for efficiency
            case 1
                for f = 1:numel(field)
                    for i = 1:numel(value{f})
                        S(i,I).(field(f)) = value{f}(i);
                    end
                end
            case 2
                for f = 1:numel(field)
                    for i = 1:numel(value{f})
                        S(I,i).(field(f)) = value{f}(i);
                    end
                end
        end
    end

    % Determine xFit indices in constraintsCellArray
    function ind = findxFitInd(parConstraint)
        switch parConstraint{1}
            case '@(a)a;' % skip fixed parameters
                ind = [];
            case '@(x,a,b)x(a)*x(b);'
                ind = [parConstraint{[2 3]}];
            case '@(x,a,b,c)x(a)+b*x(c);'
                ind = [parConstraint{[2 4]}];
            otherwise
                ind = parConstraint{2};
        end
    end

    % AMARES.amares output
    function Results = formOutput(fitResults, CRBResults, fitStatus, exptParams, figureHandle)
        Results.Linewidths = fitResults.linewidth;
        Results.Dampings = fitResults.linewidth .* pi; % TODO: This is only valid for Lorentzian lines.
        Results.Phases = fitResults.phase;
        Results.Amplitudes = fitResults.amplitude ;
        Results.GaussianSigma = fitResults.sigma ;

        %The position of the peaks expressed in all the different ways possible.
        Results.ChemicalShifts = fitResults.chemShift;
        Results.ChemicalShiftsIncOffset = fitResults.chemShift - exptParams.offset;
        Results.FrequenciesHz = fitResults.chemShift.*exptParams.imagingFrequency;
        Results.FrequenciesHzIncOffset = (fitResults.chemShift - exptParams.offset).*exptParams.imagingFrequency;
        Results.offsetPPM = exptParams.offset; % in ppm
        Results.offsetHz = exptParams.offset*exptParams.imagingFrequency; % in Hz

        Results.Standard_deviation_of_Amplitudes = CRBResults.amplitude;
        Results.Standard_deviation_of_Phases = CRBResults.phase;
        Results.Standard_deviation_of_Dampings = CRBResults.linewidth*pi;
        Results.Standard_deviation_of_FrequenciesHz = CRBResults.chemShift.*exptParams.imagingFrequency;
        Results.Standard_deviation_of_ChemicalShifts = CRBResults.chemShift;
        Results.Standard_deviation_of_Linewidths = CRBResults.linewidth;
        Results.Standard_deviation_of_GaussianSigma = CRBResults.sigma;

        Results.resFigureHandle = figureHandle;
        Results.relativeNorm = fitStatus.relativeNorm;

        Results.fitStatus = fitStatus;
    end

    function pd = PercDiff(A,B)
        if B == 0
            pd = 0;
        else
            pd = abs(A-B)/B;
        end
    end

    function S = A2S(A)
        % Convert 1x1 structure with array values to multidimensional structure
        fnames = fieldnames(A);
        np = numel(A.(fnames{1}));

        % Create structure
        S = struct();
        for f = fnames'
            for n = 1:np
                if any(ischar(A.(f{:})))
                    S(n).(f{:}) = A.(f{:}){n};
                else
                    S(n).(f{:}) = A.(f{:})(n);
                end
            end
        end
    end

    function m = CalAll(residual, noisevar, fitResults, CRBResults, FQN_global_base, MAD_global_base)
        m.FQN = CalFQN(residual,noisevar);
        m.MAD = CalMAD(residual);
        m.CoV = CalCOV(fitResults, CRBResults);
        m.met = CalMET(m.FQN, m.MAD, FQN_global_base, MAD_global_base);
    end
    
    function FQN = CalFQN(r, noisevar)
        FQN = var(r)/noisevar;
    end

    function MAD = CalMAD(r)
        MAD = mean(abs(r));
    end

    function met = CalMET(fqn,mad,fqn_base,mad_base)
        met = sqrt((abs(fqn-1)./abs(fqn_base-1)).*(mad./mad_base));
    end

    function CoV = CalCOV(r, crb, MDxs)
        CoV = crb.amplitude./r.amplitude;
        if nargin == 2
            CoV = mean(CoV);
        else
            CoV = cellfun(@(x) mean(CoV(x)), MDxs);
        end
    end

    function [fitResults, fitStatus, CRBResults, data, MDxs, PeakIndices] = FitandPlot(data, exptParams, pk, opts, noisevarTD, s, pFLAG, a2p)
        % Do quick DC correction
        DCcor = polyfit(exptParams.timeAxis(end-round(size(exptParams.timeAxis,1)*0.25):end),(data.inputFid(end-round(size(exptParams.timeAxis,1)*0.25):end)),0);

        % Fit data
        [fitResults, fitStatus, ~, CRBResults] = AMARES.amaresFit(data.inputFid - DCcor, exptParams, pk, false, opts);

        if fitStatus.EXITFLAG ~= 3
            error(fitStatus.OUTPUT.message)
        end

        [data.modelFid,~,data.modelFids] = AMARES.makeModelFidAndJacobianReIm(fitStatus.xFit,fitStatus.constraintsCellArray,exptParams.beginTime,exptParams.dwellTime,exptParams.imagingFrequency,exptParams.samples, 'complexOutput', true);
        data.modelSpec = specFft(data.modelFid);
        data.modelSpecs = specFft(data.modelFids,1);
        data.residual = data.inputSpec - data.modelSpec;

        if pFLAG
            PlotResults(fitStatus, fitResults, data, opts, a2p);
        end

        % Grab peak ranges
        [PeakIndices,MDxs] = rangePeakCorrected(fitResults,fitStatus,data.modelSpecs,pk);

        % Calculate fit quality number (noise_var from OXSA is actually variance of the residual) - https://doi.org/10.1002/mrm.1910390607
        fprintf('%+4s ||%8.2f   |%10.2e |%8.2f  ', string(s), CalFQN(data.inputFid-data.modelFid, noisevarTD), CalMAD(data.inputFid-data.modelFid), CalCOV(fitResults, CRBResults))
    end

    % Grab peak ranges
    function [PeakIndices,MDxs,indivSpectrum] = rangePeakCorrected(fitResults,fitStatus,modelSpecs,pk)
        [PeakIndices_old, MDxs, ~, ~] = AMARES.rangePeak(fitResults, fitStatus, modelSpecs);

        % Combine any separate multiplets
        indivSpectrum = [];
        PeakIndices = {};
        pDx = 0; % peakIndices iterator
        pDx1 = 0; % pk iterator
        while pDx1 < numel(pk.priorKnowledge)
            pDx1 = pDx1 + 1;
            if ~isempty(pk.priorKnowledge(pDx1).peakName)
                pDx = pDx + 1;
                if iscell(pk.priorKnowledge(pDx1).peakName)
                    peakName1 = pk.priorKnowledge(pDx1).peakName{1}(1:end-1);
                else
                    peakName1 = pk.priorKnowledge(pDx1).peakName;
                end
                PeakIndices{pDx} = PeakIndices_old{pDx1}; %#ok<AGROW>
                indivSpectrum(:,pDx) = modelSpecs(:,pDx1); %#ok<AGROW>
                if ~iscell(pk.priorKnowledge(pDx1).peakName) && pDx1 ~= numel(pk.priorKnowledge) % peak is not multiplet & not at very end
                    for pDx2 = pDx1+1:numel(pk.priorKnowledge)
                        if iscell(pk.priorKnowledge(pDx2).peakName)
                            peakName2 = pk.priorKnowledge(pDx2).peakName{1}(1:end-1);
                        else
                            peakName2 = pk.priorKnowledge(pDx2).peakName;
                        end
                        if ~strcmp(peakName1(1),'X') && numel(peakName1)>=3 && ~strcmp(peakName1(1:3),'Lip') && strcmp(peakName1(1:end-1), peakName2(1:end-1)) && ... peaks are not unknown/fat & names are equal except for last number
                                ~isempty(pk.priorKnowledge(pDx1).G_linewidth) && ~isempty(pk.priorKnowledge(pDx2).G_linewidth) && pk.priorKnowledge(pDx1).G_linewidth==pk.priorKnowledge(pDx2).G_linewidth && ... peaks are in linewidth group
                                ~isempty(pk.priorKnowledge(pDx1).G_amplitude) && ~isempty(pk.priorKnowledge(pDx2).G_amplitude)  && pk.priorKnowledge(pDx1).G_amplitude==pk.priorKnowledge(pDx2).G_amplitude % peaks are in amplitude group
                            PeakIndices{pDx} = min([PeakIndices_old{pDx}(1) PeakIndices_old{pDx2}(1)]):max([PeakIndices_old{pDx}(end) PeakIndices_old{pDx2}(end)]);
                            indivSpectrum(:,pDx) = indivSpectrum(:,pDx) + modelSpecs(:,pDx2);
                            pk.priorKnowledge(pDx2).peakName = [];
                        end
                    end
                end
            end
        end
    end

    % % Trim evaluation indices
    % function csbIndices2 = ModEvalIndices(PeakIndices,csbIndices,P2A,A2P,refDx,fitResults,exptParams)
    %     maxIndices = numel(PeakIndices{P2A(refDx)});
    %     csbIndices2 = csbIndices;
    %     for pDx = 1:numel(PeakIndices)
    %         % ensure csbIndices cover linewidth
    %         if numel(csbIndices2{pDx}) < mean(fitResults.linewidth(A2P{pDx}))/(exptParams.dwellTime*exptParams.samples)
    %             diffInd = ceil((mean(fitResults.linewidth(A2P{pDx})) - numel(csbIndices2{pDx}))/2);
    %             csbIndices2{pDx} = csbIndices2{pDx}(1)-diffInd:csbIndices2{pDx}(end)+diffInd;
    %         end
    % 
    %         % Check if larger than max size
    %         if numel(csbIndices2{pDx}) > maxIndices
    %             diffInd = ceil((numel(csbIndices2{pDx}) - maxIndices)/2);
    %             csbIndices2{pDx} = csbIndices2{pDx}(diffInd+1:end-diffInd);
    %         end
    %     end
    % end

    % Calculate FQN, AD & CoV for each peak
    function [FQN,MAD,CoV] = EvalFit(data, noisevar, fit, CRB, PeakIndices, csbIndices, MDxsP, A2P, FQN_global, MAD_global, best_data)
        %{
        Metrics have biases when calculated using:
        - Residual of global fit are biased by bad fits of neighboring peaks (method 1 & 2 & 5)
        - Residual of single peak fit are biased by presence of neighboring peaks (method 3 & 4 & 6)
        - peakIndices are biased by fit; inaccurate linewidth fit may not cover entire peak, or may cover too much (methods 1 & 3)
        - csbIndices are biased by prior knowledge; peaks that have large ranges (AcC) will have indices much larger than the actual peak (methods 2 & 4)
        Taking mean of all methods should help in minimizing biases
        %}
        Np = numel(PeakIndices);
        FQN = zeros(6,Np,'double');
        MAD = zeros(6,Np,'double');

        % Just calculate CoV normally
        CoV = CalCOV(fit, CRB, MDxsP);

        % See if global FQN % MAD are below threshold
        if nargin>8 && FQN_global<100 && MAD_global<1
            % Iterate through peaks
            for pDx = 1:Np
                mask = true(1,size(best_data.modelSpecs,2));
                mask(A2P{pDx}) = false;

                %{ 
                Method 1a: 
                - Calculate residual using best fit for all peaks except for current peak
                - Evaluate metrics over peak indices
                %}
                residual = real(data.inputSpec - sum(best_data.modelSpecs(:,mask),2) - sum(data.modelSpecs(:,A2P{pDx}),2));
                FQN(1,pDx) = CalFQN(residual(PeakIndices{pDx}), noisevar);
                MAD(1,pDx) = CalMAD(residual(PeakIndices{pDx}));

                %{ 
                Method 2: 
                - Calculate residual using best fit for all peaks except for current peak
                - Evaluate metrics over csb indices
                %}
                FQN(2,pDx) = CalFQN(residual(csbIndices{pDx}), noisevar);
                MAD(2,pDx) = CalMAD(residual(csbIndices{pDx}));

                %{
                Method 5: 
                - Calculate residual using best fit for all peaks except for current peak
                - Evaluate metrics entire range
                %}
                FQN(5,pDx) = CalFQN(residual, noisevar);
                MAD(5,pDx) = CalMAD(residual);
            end
        else
            %{
            Method 1b: 
            - Calculate residual using current fit for all peaks
            - Evaluate metrics over peak indices
            %}
            FQN(1,:) = cellfun(@(x)CalFQN(data.residual(x),noisevar), PeakIndices);
            MAD(1,:) = cellfun(@(x)CalMAD(data.residual(x)), PeakIndices);

            %{
            Method 2: 
            - Calculate residual using current fit for all peaks
            - Evaluate metrics over csb indices
            %}
            FQN(2,:) = cellfun(@(x)CalFQN(data.residual(x),noisevar), csbIndices);
            MAD(2,:) = cellfun(@(x)CalMAD(data.residual(x)), csbIndices);

            %{
            Method 5: 
            - Calculate residual using current fit for all peaks
            - Evaluate metrics over entire range
            %}
            FQN(5,:) = CalFQN(data.residual, noisevar);
            MAD(5,:) = CalMAD(data.residual);
        end

        % Calculate results via methods 3 & 4 & 6
        for pDx = 1:Np
            %{
            Method 3: 
            - Calculate residual using only current peak
            - Evaluate metrics over peak indices
            %}
            residual = real(data.inputSpec - sum(data.modelSpecs(:,A2P{pDx}),2));
            FQN(3,pDx) = CalFQN(residual(PeakIndices{pDx}), noisevar);
            MAD(3,pDx) = CalMAD(residual(PeakIndices{pDx}));

            %{
            Method 4: 
            - Calculate residual using only current peak
            - Evaluate metrics over csb indices
            %}
            % Calculate results using chemShift bounds
            FQN(4,pDx) = CalFQN(residual(csbIndices{pDx}), noisevar);
            MAD(4,pDx) = CalMAD(residual(csbIndices{pDx}));

            %{
            Method 6: 
            - Calculate residual using only current peak
            - Evaluate metrics over entire range
            %}
            FQN(6,pDx) = CalFQN(residual, noisevar);
            MAD(6,pDx) = CalMAD(residual);
        end

        % Take mean of three different techniques
        FQN = mean(FQN,1);
        MAD = mean(MAD,1);
    end

    function hFig = PlotResults(fitStatus, fitResults, data, options, a2p)
        % Perform phase correction
        if ~isfield(options,'firstOrder') || options.firstOrder % On by default
            if ~isfield(fitStatus.exptParams,'freqAxis')
                fitStatus.exptParams.freqAxis = fitStatus.exptParams.ppmAxis*fitStatus.exptParams.imagingFrequency;
            end
            
            % Set zero-order phase relative to the reference peak
            actualRefPeak = AMARES.getActualRefPeakDx(fitStatus.pkWithLinLsq);
            if ~isempty(actualRefPeak)
                zeroOrderPhaseRad = fitResults.phase(actualRefPeak)*pi/180;
            else
                zeroOrderPhaseRad = 0;
            end
            firstOrderCorrection = exp(-1i*(zeroOrderPhaseRad + 2*pi*fitStatus.exptParams.freqAxis*fitStatus.exptParams.beginTime));
        else
            firstOrderCorrection = 1;
        end
        
        data.inputSpec = data.inputSpec.*firstOrderCorrection;
        data.modelSpec = data.modelSpec.*firstOrderCorrection;
        data.modelSpecs = data.modelSpecs.*firstOrderCorrection;
        data.residual = data.inputSpec - data.modelSpec;

        % Plot with OXSA function
        [hFig, hAx] = AMARES.amaresPlot(fitStatus, 'hFig', figure(Theme='light'), options);

        % Update figure name
        FIGS = findall(groot);
        FIGS(2).Name = ['OXSA' num2str(FIGS(2).Number)];
        FIGS(2).NumberTitle = 'off';

        % Fix individual plot by combining multiplets
        if ~isfield(options,'plotIndividual') || options.plotIndividual
            % sp = 2;
            % hAx(sp).YAxis.Visible = 'off'; hAx(sp).XAxis.Visible = 'off';
            % hAx(sp) = subplotSetBorder(length(hAx), 1, [.1 .03 0.03 0], [0 0 0.05 0.02], sp);
        
            indivColours = distinguishable_colors(numel(a2p));
            
            % Correct for 1H zero ppm
            if strcmp(options.nucleus,'1H')
                fitStatus.xFit(1:fitStatus.constraintsCellArray.chemShift{end}{2}) = fitStatus.xFit(1:fitStatus.constraintsCellArray.chemShift{end}{2}) - fitStatus.exptParams.offset;
                fitStatus.exptParams.ppmAxis = fitStatus.exptParams.ppmAxis - fitStatus.exptParams.offset;
                fitResults.chemShift = fitResults.chemShift - fitStatus.exptParams.offset;
            end

            % Combine multiplets & calculate range to plot
            [peakPlotIndices, ~, ~, indivSpectrums] = AMARES.rangePeak(fitResults, fitStatus, data.modelSpecs);

            % Determine line indices
            lDx = 3; % move past root
            nFs = 1; % number of figures
            while strcmp(FIGS(lDx).Type,'figure') && numel(FIGS(lDx).Name)>4 && strcmp(FIGS(lDx).Name(1:4),'OXSA') % move past OXSA figures
                lDx = lDx + 1;
                nFs = nFs + 1;
            end
            nP = size(data.modelSpecs,2); % changed from size(indivSpectrums,2) bc was not deleting all subpeaks in 31P data - 5/15/26
            while lDx+nP<numel(FIGS) && ~(strcmp(FIGS(lDx+nP).Type,'text') && strcmp(FIGS(lDx+nP).String,'Individual Peaks')) % move to first line on individual peaks
                lDx = lDx + 1;
            end
            lDx = lDx+nP-1:-1:lDx;
            % Plot
            for pDx = 1:nP
                if pDx <= numel(a2p)
                    if isscalar(a2p{pDx})
                        peakPlotDxs = peakPlotIndices{a2p{pDx}};
                        indivSpectrum = indivSpectrums(peakPlotDxs,a2p{pDx});
                    else
                        peakPlotDxs = min([peakPlotIndices{a2p{pDx}(1)}(1) peakPlotIndices{a2p{pDx}(end)}(1)]);
                        peakPlotDxs = peakPlotDxs:max([peakPlotIndices{a2p{pDx}(1)}(end) peakPlotIndices{a2p{pDx}(end)}(end)]);
                        indivSpectrum = real(sum(indivSpectrums(peakPlotDxs,a2p{pDx}),2));
                    end

                    set(FIGS(lDx(pDx)), {'XData','YData','color'}, {fitStatus.exptParams.ppmAxis(peakPlotDxs),real(indivSpectrum),indivColours(pDx,:)})
                else % Remove remaining plots
                    delete(FIGS(lDx(pDx)))
                end
            end
        end
    end

    % Update ylimits for all
    function ylims = UpdateYlims(ylims)
        FIGS = findall(groot);
        lDx = 3; % root & first figure
        nFs = 1;
        while strcmp(FIGS(lDx).Type,'figure') && numel(FIGS(lDx).Name)>4 && strcmp(FIGS(lDx).Name(1:4),'OXSA') % move past OXSA figures
            lDx = lDx + 1;
            nFs = nFs + 1;
        end
        
        % Move through figure axes
        AX = [];
        for i = 1:nFs
            % Move past non-axes or axes not associated with OXSA figures
            while ~strcmp(FIGS(lDx).Type,'axes') || ~(strcmp(FIGS(lDx).Parent.Type,'figure') && numel(FIGS(lDx).Parent.Name)>=4 && strcmp(FIGS(lDx).Parent.Name(1:4),'OXSA')) 
                lDx = lDx + 1;
            end

            % Ensure axes is not a phase plot
            if ~any(abs(round(FIGS(lDx).YLim)) == 3)
                % Check ylimits (ylimits for all axes in a single figure are the same)
                ylims(1) = min([ylims(1) FIGS(lDx).YLim(1)]);
                ylims(2) = max([ylims(2) FIGS(lDx).YLim(2)]);
    
                % Add axes index to list
                while strcmp(FIGS(lDx).Type,'axes')
                    AX = [AX lDx]; %#ok<AGROW>
                    lDx = lDx + 1;
                end
            end
        end
        
        % Update axes to all have the same range
        for ax = AX
            FIGS(ax).YLim = ylims;
        end
    end

    % function PrintResults(Names,fitResults,P2A,Quality)
    %     fprintf('Peak\t\tchemShift\tlinewidth\tamplitude\tphase\n')
    %     if isfield(fitResults,'sigma') && any(fitResults.sigma~=0)
    %         fprintf('\b\t\tsigma\n')
    %         sigmaFLAG = true;
    %     else
    %         sigmaFLAG = false;
    %     end
    %     for pDx = 1:numel(Names)
    %         if nargin<5 || Quality(P2A(pDx))>0
    %             fprintf('%s\t\t%0.2f\t\t%0.1f\t\t%0.2e\t%0.0f\n',Names{pDx},fitResults.chemShift(pDx),fitResults.linewidth(pDx),fitResults.amplitude(pDx),fitResults.phase(pDx))
    %             if sigmaFLAG
    %                 fprintf('\b\t\t%0.2e\n',fitResults.sigma(pDx))
    %             end
    %         end
    %     end
    % end
end