function [results, results_c] = OXSAResultsTable(res)
% OXSARESULTSTABLE converts results from OXSA AMARES and outputs a table

exptParams = res.Status.exptParams;
pk = res.Status.pkWithLinLsq;
fitResults = res.Parameters;
fitResultSDs = res.ParameterSDs;

nu = length(fitResults.linewidth); % number of peaks, multiplets uncombined
nc = length(pk.initialValues); % number of peaks, multiplets combined

Peak = cell(nu,1);
ChemicalShift = zeros(nu,1);
Linewidth = zeros(nu,1);
Amplitude = zeros(nu,1);
Phase = zeros(nu,1);
Sigma = zeros(nu,1);
Peak2 = cell(nc,1);
ChemicalShift2 = zeros(nc,1);
Linewidth2 = zeros(nc,1);
Amplitude2 = zeros(nc,1);
Phase2 = zeros(nc,1);
Sigma2 = zeros(nc,1);

idx = 0;
for pDx = 1:nc
    % Determine number of peaks to add
    idx0 = idx + 1;
    if iscell(pk.bounds(pDx).peakName)
        idx = idx0 + length(pk.bounds(pDx).peakName) - 1;
    else
        idx = idx0;
    end

    if iscell(pk.bounds(pDx).peakName)
        for j = 1:idx-idx0+1
            Peak{idx0+j-1} = pk.bounds(pDx).peakName{j};
        end
    else
        Peak{idx0:idx} = pk.bounds(pDx).peakName;
    end
    ChemicalShift(idx0:idx) = fitResults.chemShift(idx0:idx);
    Linewidth(idx0:idx) = fitResults.linewidth(idx0:idx);
    Amplitude(idx0:idx) = fitResults.amplitude(idx0:idx);
    Phase(idx0:idx) = fitResults.phase(idx0:idx);
    Sigma(idx0:idx) = fitResults.sigma(idx0:idx);

    if idx0 ~= idx
        Peak2{pDx} = Peak{idx0}(1:end-1);
        ChemicalShift2(pDx) = mean(ChemicalShift(idx0:idx));
        jcoup = (ChemicalShift(idx0+1)-ChemicalShift(idx0))*exptParams.imagingFrequency;
        Linewidth2(pDx) = Linewidth(idx0); % Linewidth2(pDx) = Linewidth(idx0) + jcoup*(idx-idx0);
    else
        Peak2{pDx} = Peak{idx0:idx};
        ChemicalShift2(pDx) = ChemicalShift(idx0:idx);
        Linewidth2(pDx) = Linewidth(idx0:idx);
    end
    Amplitude2(pDx) = max(Amplitude(idx0:idx));
    Phase2(pDx) = Phase(idx0);
    Sigma2(pDx) = Sigma(idx0);
end

Amplitude2 = Amplitude2./max(Amplitude2);
results = table(round(ChemicalShift,2), round(Linewidth,2), round(Phase), Amplitude, round(Sigma,1), 'RowNames', Peak, 'VariableNames', {'ChemSh', 'LW', 'Phase', 'Amp', 'Sigma'});
results_c = table(round(ChemicalShift2,2), round(Linewidth2,2), round(Phase2), Amplitude2, round(Sigma2,1), 'RowNames', Peak2, 'VariableNames', {'ChemSh', 'LW', 'Phase', 'Amp', 'Sigma'});
end