function [fitStatus] = adjust(fitStatus, peak, param, mult, plotFLAG)
peak = char(peak);

data.inputSpec = specFft(fitStatus.inputFid);
data.inputFid = fitStatus.inputFid;

% Generate initial data
data_og = data;
[data_og.modelSpec,data_og.modelSpecs,data_og.modelFid,data_og.modelFids] = AMARES.makeModelSpec(fitStatus,struct('firstOrder',false));

% Find index of peak
dx = 0;
for i = 1:numel(fitStatus.pkWithLinLsq.initialValues)
    if ~iscell(fitStatus.pkWithLinLsq.initialValues(i).peakName)
        if strcmp(fitStatus.pkWithLinLsq.initialValues(i).peakName, peak)
            dx = i;
            break
        end
    elseif length(fitStatus.pkWithLinLsq.initialValues(i).peakName{1}) > length(peak)
        if strcmp(fitStatus.pkWithLinLsq.initialValues(i).peakName{1}(1:length(peak)), peak)
            dx = i;
            break
        end
    end
end

if dx == 0
    error('Peak could not be found!')
end

% Modify parameter
mDx = OXSA.findMDX(fitStatus.pkWithLinLsq, numel(fitStatus.pkWithLinLsq.bounds));
if iscell(param)
    param = convertCharsToStrings(param);
end
for par = 1:numel(param)
    if any(strcmp(param(par), 'phase'))
        if strcmp(fitStatus.constraintsCellArray.(param(par)){mDx{dx}}{1}, '@(a)a;')
            fitStatus.constraintsCellArray.(param(par)){mDx{dx}}{2} = mult(par);
        else
            fitStatus.xFit(fitStatus.constraintsCellArray.(param(par)){mDx{dx}}{2}) = mult(par);
        end
    else
        for sub = 1:numel(mDx{dx})
            fitStatus.xFit(fitStatus.constraintsCellArray.(param(par)){mDx{dx}(sub)}{2}) = fitStatus.xFit(fitStatus.constraintsCellArray.(param(par)){mDx{dx}(sub)}{2})*mult(par);
        end
    end
end

[data.modelSpec,data.modelSpecs,data.modelFid,data.modelFids] = AMARES.makeModelSpec(fitStatus,struct('firstOrder',false));
fitStatus.residual = data.inputFid - data.modelFid;
fitStatus.noise_var = var(fitStatus.residual);
data.residual = data.inputSpec - data.modelSpec;

if nargin == 4 || plotFLAG
    firstOrderCorrection = exp(-1i*(2*pi*fitStatus.exptParams.ppmAxis*fitStatus.exptParams.imagingFrequency*fitStatus.exptParams.beginTime));

    figure(Theme='light'); plot(fitStatus.exptParams.ppmAxis, real(data.inputSpec.*firstOrderCorrection)); axis tight; hold on
    plot(fitStatus.exptParams.ppmAxis, real(data_og.modelSpec.*firstOrderCorrection), ":");
    plot(fitStatus.exptParams.ppmAxis, real(data.modelSpec.*firstOrderCorrection), "--");
    ylabel('Spectrum fit')
    xlabel('\delta / ppm')
    set(gca,'XDir','reverse')
end
legend(["True" "Fit - Original" "Fit - New"])