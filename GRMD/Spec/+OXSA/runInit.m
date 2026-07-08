function [data, fitResults, CRBResults, fitStatus, pk] = runInit(pk, exptParams, data)

% Get first results
[data.modelFid, data.modelFids] = AMARES.makeInitialValuesModelFid(pk, exptParams);

% Find difference
adiff = max(abs(data.inputFid))/max(abs(data.modelFid));

% Apply difference to initial values
for j = 1:numel(pk.initialValues)
    pk.initialValues(j).amplitude = pk.initialValues(j).amplitude*adiff;
end

% Generate fit with corrected amplitude
[data.modelFid, data.modelFids] = AMARES.makeInitialValuesModelFid(pk, exptParams);
data.modelSpec = specFft(data.modelFid);
data.modelSpecs = specFft(data.modelFids);

% Construct fitStatus structure
fitStatus = struct('exptParams', exptParams, 'pkWithLinLsq', pk);

% Create xFit & constraints solely based on pk initial values
[fitStatus.xFit, fitStatus.constraintsCellArray] = xFitFake(fitStatus.pkWithLinLsq);
fitResults = AMARES.applyModelConstraints(fitStatus.xFit,fitStatus.constraintsCellArray);

fitStatus.residual = data.inputFid - data.modelFid;
fitStatus.noise_var = var(fitStatus.residual);
data.residual = data.inputSpec - data.modelSpec;
CRBResults = AMARES.estimateCRB(exptParams.imagingFrequency, exptParams.dwellTime, exptParams.beginTime, fitStatus.noise_var, fitStatus.xFit, fitStatus.constraintsCellArray);

% Generate fake xFit data
function [xFit, constraintsCellArray] = xFitFake(pk)
    % set initial values, prior knowledge and lower/upper pkWithLinLsq.bounds
    [xFit, ~, ~, optimIndex] = AMARES.initializeOptimization(pk);
    
    % create cell arrays of constraints for the makeModelFid function
    constraintsCellArray = AMARES.createModelConstraints(pk, optimIndex);
end
end
