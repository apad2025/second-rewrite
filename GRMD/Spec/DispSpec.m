clear

dtype = 'Phosphorus';
dogs = ["Waylon","Sushi","Selene","Aphrodite"];
dates = [0 76 195 297];
params = ["chemShift", "linewidth", "amplitude", "phase", "sigma"];
params_all = [params, "chemShiftDelta"];

switch dtype
    case 'Phosphorus'
        renamePeaks = dictionary(["bATP",  "NAD", "NADH", "tNAD", "aATP", "gATP", "PCr", "MP", "GPC", "GPE", "PDE", "Unknown", "Pia", "Pib", "Pi",  "PC",  "PE", "PME", "G6P", "PPA", "Ref"], ...
                                 [ "ATP", "tNAD", "tNAD", "tNAD",  "ATP",  "ATP", "PCr", "MP", "PDE", "PDE", "PDE", "Unknown",  "Pi",  "Pi", "Pi", "PME", "PME", "PME", "PME", "PPA", "PPA"]);
    case 'Hydrogen'
        error('Not added yet!')
end
pkNames_all = unique(values(renamePeaks),'stable')';
pkIndices_all = dictionary(pkNames_all, 1:numel(pkNames_all));

switch dtype
    case 'Phosphorus'
        T1 = struct('bATP',struct('Mean',mean([3.5 3.9]), 'Err', 0.5*EP('addition',[1 1],[1.1 0.4],0)));
        T1.aATP = struct('Mean',mean([2.6 3.4]), 'Err', 0.5*EP('addition',[1 1],[0.9 0.3],0));
        T1.gATP = struct('Mean',mean([4.5 5.5]), 'Err', 0.5*EP('addition',[1 1],[0.3 0.4],0));
        T1.PCr = struct('Mean',mean([6.4 6.7]), 'Err', 0.5*EP('addition',[1 1],[0.2 0.4],0));
        T1.PDE = struct('Mean',8.6, 'Err',1.2);
        T1.Pi = struct('Mean',mean([5.2 6.9]), 'Err', 0.5*EP('addition',[1 1],[0.6 1.0],0));
        T1.PME = struct('Mean',8.1, 'Err',1.7);
    case 'Hydrogen'
        error('Not added yet!')
end

% Flag which peaks need to be removed later
flag = true(size(pkNames_all));
for p = 1:numel(pkNames_all)
    if ~isfield(T1, pkNames_all(p)) && ~strcmp(pkNames_all(p),"ATP")
        flag(p) = false;
    end
end

modResults_all = cell(numel(dogs),1);
pkResults_all = modResults_all;
twix_all = modResults_all;
for s = 1:numel(dogs)
    [~,twix_all{s},modResults_all{s},pkResults_all{s}] = dogexplorer(dogs(s), 'Analyze', 'DataType', dtype, 'verbose',false);
end

%% Combine & remove first day
dxDiff = cellfun(@numel, modResults_all) - min(cellfun(@numel, modResults_all));
modResults = modResults_all{1}(dxDiff(1)+1:end);
pkResults = pkResults_all{1}(:,1,dxDiff(1)+1:end);
twix = twix_all{1}(1,dxDiff(1)+1:end);
for s = 2:numel(dogs)
    modResults(s,:) = modResults_all{s}(dxDiff(s)+1:end);
    pkResults(:,s,:) = pkResults_all{s}(:,1,dxDiff(s)+1:end);
    twix(s,:) = twix_all{s}(1,dxDiff(s)+1:end);
end

% Quantify date
ds = datetime([modResults(:,:).Date], InputFormat="yyMMdd");
d0 = min(ds);
for s = 1:numel(dogs)
    for d = 1:size(modResults,2)
        modResults(s,d).Date = days(datetime(modResults(s,d).Date, InputFormat="yyMMdd") - d0);
    end
end
clear dxDiff ds d0

%% Remove peaks with unknown T1 values
% Remove from xFit, constraintsCellArray, & pk
for s = 1:numel(dogs)
    for d = 1:numel(dates)
        % Check if flagged peaks are present
        subFlag = flag;
        for p = 1:numel(pkNames_all)
            if isempty(pkResults(p,s,d).peakDx) && ~subFlag(p)
                subFlag(p) = true;
            end
        end

        % Indices to remove
        mDxFlag = true([1 numel(modResults(s,d).constraintsCellArray.chemShift)]);
        is = [pkResults(~subFlag,s,d).multDx];
        mDxFlag([is{:}]) = false;

        % Remove from modelSpecs
        modResults(s,d).modelSpecs = modResults(s,d).modelSpecs(:,mDxFlag);
        % Remove from xFit
        xFitFlag = true(size(modResults(s,d).xFit));
        for param = params
            for p = 1:numel(mDxFlag)
                if ~mDxFlag(p)
                    switch modResults(s,d).constraintsCellArray.(param){p}{1}
                        case '@(a)a;'% fixed values aren't present in xFit
                        case '@(x,a,b)x(a)*x(b);' % two parameters need to be removed from xFit
                            xFitFlag([modResults(s,d).constraintsCellArray.(param){p}{[2 3]}]) = false;
                        case '@(x,a,b,c)x(a)+b*x(c);'% two parameters need to be removed from xFit
                            xFitFlag([modResults(s,d).constraintsCellArray.(param){p}{[2 4]}]) = false;
                        otherwise % only one parameter needs to be removed from xFit
                            xFitFlag(modResults(s,d).constraintsCellArray.(param){p}{2}) = false;
                    end
                end
            end
        end
        modResults(s,d).xFit = modResults(s,d).xFit(xFitFlag);

        % Correct indices in constraintsCellArray
        P = 1:numel(mDxFlag); % indices to check
        for param = params
            for p = P(mDxFlag)
                % Grab index
                switch modResults(s,d).constraintsCellArray.(param){p}{1}
                    case '@(a)a;' % no indices need to be fixed
                        constDx = [];
                    case '@(x,a,b)x(a)*x(b);' % two indices need to be fixed
                        constDx = [2 3];
                    case '@(x,a,b,c)x(a)+b*x(c);' % two indices need to be fixed
                        constDx = [2 4];
                    otherwise % only one index needs to be fixed
                        constDx = 2;
                end

                % Apply correction
                for pI = 1:numel(constDx)
                    if sum(~xFitFlag(1:modResults(s,d).constraintsCellArray.(param){p}{constDx(pI)})) ~= 0
                        modResults(s,d).constraintsCellArray.(param){p}{constDx(pI)} = modResults(s,d).constraintsCellArray.(param){p}{constDx(pI)} - sum(~xFitFlag(1:modResults(s,d).constraintsCellArray.(param){p}{constDx(pI)}));
                    end
                end
            end

            % Remove from constraintsCellArray
            modResults(s,d).constraintsCellArray.(param) = modResults(s,d).constraintsCellArray.(param)(mDxFlag);
        end

        % Indices to remove
        pDxFlag = true(size(modResults(s,d).pk.bounds));
        pDxFlag([pkResults(~subFlag,s,d).peakDx]) = false;

        % Remove from pk
        modResults(s,d).pk.bounds = modResults(s,d).pk.bounds(pDxFlag);
        modResults(s,d).pk.initialValues = modResults(s,d).pk.initialValues(pDxFlag);
        modResults(s,d).pk.priorKnowledge = modResults(s,d).pk.priorKnowledge(pDxFlag);

        % Recalculate covariance matrix
        CRBResults = AMARES.estimateCRB(modResults(s,d).exptParams.imagingFrequency, modResults(s,d).exptParams.dwellTime, modResults(s,d).exptParams.beginTime, modResults(s,d).noise_var, modResults(s,d).xFit, modResults(s,d).constraintsCellArray);
        modResults(s,d).covariance = CRBResults.covariance;

        % Correct multDx & peakDx in pkResults
        for p = 1:numel(pkNames_all)
            % Ignore peaks that are going to be removed
            if subFlag(p)
                % Loop through peakDx
                for pI = 1:numel(pkResults(p,s,d).peakDx)
                    pkResults(p,s,d).peakDx(pI) = pkResults(p,s,d).peakDx(pI) - sum(~pDxFlag(1:pkResults(p,s,d).peakDx(pI)));
                end

                % Loop through multDx
                for mI = 1:numel(pkResults(p,s,d).multDx)
                    for mI2 = 1:numel(pkResults(p,s,d).multDx{mI})
                        pkResults(p,s,d).multDx{mI}(mI2) = pkResults(p,s,d).multDx{mI}(mI2) - sum(~mDxFlag(1:pkResults(p,s,d).multDx{mI}(mI2)));
                    end
                end
            end
        end
    end
end

% Remove from pkResults & pkIndices
pkResults = pkResults(flag,:,:);
pkIndices = remove(pkIndices_all, pkNames_all(~flag));
pkNames = pkNames_all(flag);

% Correct Dxs in pkIndices
for p = 1:numel(pkNames)
    pkIndices(pkNames(p)) = p;
end
clear subFlag pDxFlag mDxFlag xFitFlag mI mI2 param is constDx P

%% Apply T1 correction
for s = 1:numel(dogs)
    for d = 1:numel(dates)
        % Extract data from structure
        ampConstraints = modResults(s,d).constraintsCellArray.amplitude;
        exptParams = modResults(s,d).exptParams;
        for p = 1:numel(pkNames)
            if strcmp(pkNames(p),"ATP") % check for special case of ATP, only major group where T1 values for smaller peaks are known
                pk = {'b','a','g'};
                for atp = 1:3
                    % Grab index
                    atpDx = ampConstraints{pkResults(p,s,d).multDx{atp}(1)}{2};

                    if atp==2 && atpDx==ampConstraints{pkResults(p,s,d).multDx{atp+1}(1)}{2} % if gATP & aATP have equal amplitude
                        % Find average and propagated uncertainty
                        t1_mn = mean([T1.aATP.Mean T1.gATP.Mean]);
                        t1_err = 0.5*EP('addition', [1 1], [T1.aATP.Err T1.gATP.Err], 0);

                        [pkResults(p,s,d).amplitude(atp:atp+1), pkResults(p,s,d).amplitudeCRB(atp:atp+1)] = T1cor(...
                            t1_mn, t1_err, ...
                            twix(s,d).tr/1000, twix(s,d).fa, ...
                            pkResults(p,s,d).amplitude(atp), pkResults(p,s,d).amplitudeCRB(atp));

                        modResults(s,d).xFit(atpDx) = T1cor(...
                            t1_mn, t1_err, ...
                            twix(s,d).tr/1000, twix(s,d).fa, ...
                            modResults(s,d).xFit(atpDx));

                        break
                    else
                        [pkResults(p,s,d).amplitude(atp), pkResults(p,s,d).amplitudeCRB(atp)] = T1cor(...
                            T1.([pk{atp} 'ATP']).Mean, T1.([pk{atp} 'ATP']).Err, ...
                            twix(s,d).tr/1000, twix(s,d).fa, ...
                            pkResults(p,s,d).amplitude(atp), pkResults(p,s,d).amplitudeCRB(atp));

                        modResults(s,d).xFit(atpDx) = T1cor(...
                            T1.([pk{atp} 'ATP']).Mean, T1.([pk{atp} 'ATP']).Err, ...
                            twix(s,d).tr/1000, twix(s,d).fa, ...
                            modResults(s,d).xFit(atpDx));
                    end
                end
            elseif isfield(T1, pkNames(p)) % check if name is not present
                for subp = 1:numel(pkResults(p,s,d).multDx)
                    % Grab index
                    Dx = ampConstraints{pkResults(p,s,d).multDx{subp}(1)}{2};

                    [pkResults(p,s,d).amplitude(subp), pkResults(p,s,d).amplitudeCRB(subp)] = T1cor(...
                        T1.(pkNames(p)).Mean, T1.(pkNames(p)).Err, ...
                        twix(s,d).tr/1000, twix(s,d).fa, ...
                        pkResults(p,s,d).amplitude(subp), pkResults(p,s,d).amplitudeCRB(subp));

                    modResults(s,d).xFit(Dx) = T1cor(...
                        T1.(pkNames(p)).Mean, T1.(pkNames(p)).Err, ...
                        twix(s,d).tr/1000, twix(s,d).fa, ...
                        modResults(s,d).xFit(Dx));
                end
            else
                error('Missing metabolite')
            end
        end

        % Re-calculate CRBs & covariance with new point values
        CRBResults = AMARES.estimateCRB(exptParams.imagingFrequency, exptParams.dwellTime, exptParams.beginTime, modResults(s,d).noise_var, modResults(s,d).xFit, modResults(s,d).constraintsCellArray);
        modResults(s,d).covariance = CRBResults.covariance;
    end
end
clear ampConstraints exptParams atp atpDx t1_mn t1_err subp pk

%% Quantify data in ratios
switch dtype
    case 'Phosphorus'
        refPeak = "ATP";
    case 'Hydrogen'
        error('Not added yet!')
end

for s = 1:numel(dogs)
    for d = 1:numel(dates)
        % Check if reference peak is present
        if ~isempty(pkResults(pkIndices(refPeak),s,d).multDx)
            % Select peak indices
            i = 1:size(pkResults,1);
            mask = ~cellfun(@(x) isempty(x) || (isscalar(x) && strcmp(x,"")), {pkResults(:,s,d).subPeakNames});
            mask(pkIndices(refPeak)) = false;
            i = i(mask);
    
            % Determine numerator & denominator
            num = {pkResults(i,s,d).subPeakNames};
            denom = pkResults(pkIndices(refPeak),s,d).subPeakNames;
    
            % Grab total number of peaks
            mults_num = cellfun(@(x) cellfun(@numel, x), {pkResults(i,s,d).multDx}, 'UniformOutput',false);
            num = cellfun(@(x,y) arrayfun(@(a,b) repmat(a,1,b), x, y, 'UniformOutput',false), num, mults_num,'UniformOutput',false);
            for n = 1:numel(mults_num)
                for subn = 1:numel(mults_num{n})
                    if mults_num{n}(subn) > 1
                        num{n}{subn} = num{n}{subn} + string(1:mults_num{n}(subn)); % add numbering to multiplets
                    end
                end
            end
            mults_denom = cellfun(@numel, pkResults(pkIndices(refPeak),s,d).multDx);
            denom = arrayfun(@(x,y) repmat(x,1,y), denom, mults_denom, 'UniformOutput', false);
            for subn = 1:numel(mults_denom)
                if mults_denom(subn) > 1
                    denom{subn} = denom{subn} + string(1:mults_denom(subn));
                end
            end
    
            % Extract from the added bottom cell layer
            num = cellfun(@(x) horzcat(x{:}), num, 'UniformOutput', false);
            denom = horzcat(denom{:});
    
            % Add amplitude symbol
            num = cellfun(@(x) x + "_am", num, 'UniformOutput', false);
            denom = denom + "_am";
    
            % Add summation symbol, where applicable
            for n = 1:numel(num)
                if ~isscalar(num{n})
                    num{n} = join(num{n}," + ");
                end
            end
            if ~isscalar(denom)
                denom = join(denom," + ");
            end
    
            % Add parenthesis
            num = cellfun(@(x) "(" + x + ")", num, 'UniformOutput', false);
            denom = "(" + denom + ")";
    
            % Combine & convert to char
            derivedStr = cellfun(@(x) x + " / " + denom, num, 'UniformOutput', false);
            derivedStr = cellfun(@char, derivedStr, 'UniformOutput', false);
    
            % Calculate ratio
            [rVal, rCRB] = AMARES.estimateDerivedParamAndCRB(modResults(s,d).pk, ...
                                                             modResults(s,d).xFit, ...
                                                             modResults(s,d).constraintsCellArray, ...
                                                             modResults(s,d).covariance, ...
                                                             derivedStr);
        
            % Reorganize by peak
            Dx = 0;
            for p = 1:numel(pkNames)
                if mask(p)
                    Dx = Dx + 1;
                    pkResults(p,s,d).ratio = rVal(Dx);
                    pkResults(p,s,d).ratioCRB = rCRB(Dx);
                else
                    pkResults(p,s,d).ratio = NaN;
                    pkResults(p,s,d).ratioCRB = NaN;
                end
            end
        else
            for p = 1:numel(pkNames)
                pkResults(p,s,d).ratio = NaN;
                pkResults(p,s,d).ratioCRB = NaN;
            end
        end
    end
end
clear d i mask num denom n subn mults_num mults_denom derivedStr rVal rCRB Dx

%% Plot data
figure(Theme='light'); tiledlayout();
for p = 1:numel(pkNames)
    if p ~= pkIndices(refPeak) && all(~strcmp(pkNames(p), ["PPA","Unknown"])) % Skip ref peak & PPA index
        nexttile; hold on;
        for s = 1:size(modResults,1)
            errorbar([modResults(s,:).Date], [pkResults(p,s,:).ratio], [pkResults(p,s,:).ratioCRB]);
        end
        xlabel('Day'); ylabel('A.U.'); title(join([pkNames(p) "/" refPeak],"")); xlim([-5 305]); hold off;

        % Check for negative
        ylims = ylim;
        if ylims(1) < 0
            ylim([0 ylims(2)])
        end

        if any(ylim > 0.25)
            ylim(round(ylim,1))
        end
    end
end
legend(dogs)
clear p s ylims

%% Functions
function CRB = EP(type, ind, crb, cv, data)
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

function [Ac, Acstd] = T1cor(T1, T1std, TR, alpha, A, Astd)
% T1_CORRECTION Computes the T1-corrected amplitude and its propagated uncertainty.
%
% INPUTS:
%   A         - Measured peak amplitude
%   sigma_A   - Standard deviation of peak area
%   T1        - T1 relaxation time (same units as TR)
%   sigma_T1  - Standard deviation of T1
%   TR        - Repetition time (same units as T1)
%   alpha     - Flip angle in radians
%
% OUTPUTS:
%   Ac        - T1-corrected amplitude
%   sigma_Ac  - Propagated uncertainty (standard deviation) of corrected amplitude
%
% USAGE EXAMPLE:
%   [Ac, sigma_Ac] = T1_correction(100, 2, 1.5, 0.1, 3.0, pi/3)

    % --- Precompute common terms ---
    E        = exp(-TR./T1);
    cos_a    = cosd(alpha);

    % --- Corrected area ---
    Ac = A.*(1 - cos_a.*E)./(1 - E);

    % --- Partial derivative: dAc/dA ---
    dAc_dA = (1 - cos_a.*E)./(1 - E);

    % --- Partial derivative: dAc/dT1 ---
    dAc_dT1 = A.*(1 - cos_a).*E.*TR./(T1.^2.*(1 - E).^2);

    % --- Propagated uncertainty ---
    if nargin == 6
        Acstd = sqrt((dAc_dA.*Astd).^2 + (dAc_dT1.*T1std).^2);
    end
end

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
        case "chemShift"
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
    if ~isscalar(pkNames)
        pkNames = join(pkNames," + ");
    end
    
    % Add parenthesis
    derivedStr = "(" + pkNames + ")";
end