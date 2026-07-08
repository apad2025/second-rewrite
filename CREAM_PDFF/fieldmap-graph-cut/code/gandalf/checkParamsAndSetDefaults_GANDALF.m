function VARPROParamsout = checkParamsAndSetDefaults_GANDALF(imDataParams, algoParams, verbose)

    % Check for verbose flag
    if nargin < 3
        verbose = false;
    end

    validParams = 1;

    % Start by checking validity of provided data and recon parameters
    if size(imDataParams,3) > 1
      disp('ERROR: 2D recon -- please format input data as array of size SX x SY x 1 X nCoils X nTE')
      validParams = 0;
    end

    if length(algoParams.species) > 2
      disp('ERROR: Water=fat recon -- use a multi-species function to separate more than 2 chemical species')
      validParams = 0;
    end

    if length(imDataParams.TE) < 3
      disp('ERROR: 3+ point recon -- please use a different recon for acquisitions with fewer than 3 TEs')
      validParams = 0;
    end
    VARPROParamsout.validParams = validParams;
    
    dTE = diff(imDataParams.TE);
    if sum(abs(dTE - dTE(1))) < 1e-6 % If we have uniform TE spacing
        dt = imDataParams.TE(2) - imDataParams.TE(1);
        if verbose
            fprintf('Uniform TE spacing: Period = TE(2) - TE(1)\n');
        end
    else % options:dt = #1) TE2-TE1, #2) TE7-TE6, #3) (#1+#2)/2 #4 (#1*6+#2)/7
        dt = imDataParams.TE(4) - imDataParams.TE(3);
        if verbose
            fprintf('Non-uniform TE spacing detected, assuming 2 UTE echos: Period = TE(4) - TE(3)\n');
        end
    end    
    
    VARPROParamsout.species = algoParams.species;
    VARPROParamsout.useCUDA = set_option(algoParams, 'useCUDA', 1);
    VARPROParamsout.range_r2star = set_option(algoParams, 'range_r2star', [0 500]);
    VARPROParamsout.NUM_R2STARS = set_option(algoParams, 'NUM_R2STARS', 26);
    VARPROParamsout.sampling_stepsize = set_option(algoParams, 'sampling_stepsize', 2);
    VARPROParamsout.nSamplingPeriods = set_option(algoParams, 'nSamplingPeriods', 1);
    VARPROParamsout.airSignalThreshold_percent = set_option(algoParams, 'airSignalThreshold_percent', 5);
    VARPROParamsout.gyro = set_option(algoParams, 'gyro', 42.5774780505984);
    VARPROParamsout.period = abs(1/dt);
    
    if isfield(algoParams, 'range_fm')
        VARPROParamsout.range_fm = algoParams.range_fm;
        VARPROParamsout.nSamplingPeriods = abs(diff(VARPROParamsout.range_fm))/VARPROParamsout.period;
    else
        VARPROParamsout.nSamplingPeriods = set_option(algoParams, 'nSamplingPeriods', 1);
        VARPROParamsout.range_fm = [(-VARPROParamsout.nSamplingPeriods*VARPROParamsout.period/2) (2*VARPROParamsout.sampling_stepsize + VARPROParamsout.nSamplingPeriods*VARPROParamsout.period/2)];
    end
    
    disctreizationintervall = ceil(diff(VARPROParamsout.range_fm));
    Numlayers = ceil(disctreizationintervall/VARPROParamsout.sampling_stepsize);
    if verbose
        fprintf('Numlayers = %i', Numlayers)
    end
    
    VARPROParamsout.NUM_FMS = Numlayers;
    t = linspace(VARPROParamsout.range_fm(1), VARPROParamsout.range_fm(2), VARPROParamsout.NUM_FMS);
    gridspacing = t(2)-t(1);

    if verbose
        fprintf('\ngridspacing = %f', gridspacing)
        switch length(VARPROParamsout.range_fm)
            case 2
                fprintf('\n              %i %i\n', VARPROParamsout.range_fm(1), VARPROParamsout.range_fm(2))
            case 3
                fprintf('\n              %i %i %i\n', VARPROParamsout.range_fm(1), VARPROParamsout.range_fm(2), VARPROParamsout.range_fm(3))
        end
    end
    VARPROParamsout.gridspacing = gridspacing;
end