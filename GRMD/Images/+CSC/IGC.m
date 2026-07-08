% Perform magnitude & complex fitting via hernando's iterative graph cut
function [magn, comp] = IGC(images, dataParams, algoParams_mag, verboseFLAG, mnw, mixedFLAG)
    SL = size(images,3);

    % Data preallocation
    R2Smagn = cell(SL,1);
    FMmagn = R2Smagn;
    Fmagn = R2Smagn;
    Wmagn = R2Smagn;
    dataParams_tmp = R2Smagn;

    % Isolate mixed parameters
    if mixedFLAG
        mixedparams = algoParams_mag.mixed;
        algoParams_mag = rmfield(algoParams_mag, 'mixed');
        algoParams_mix = algoParams_mag;
        algoParams_mix.NUM_MAGN = mixedparams.NUM_MAGN;
        algoParams_mix.THRESHOLD = mixedparams.THRESHOLD;
        algoParams_mix.range_r2star = mixedparams.range_r2star;

        % Data preallocation
        R2Scomp = R2Smagn;
        FMcomp = R2Smagn;
        Fcomp = R2Smagn;
        Wcomp = R2Smagn;
        algoParams_mix_tmp = R2Smagn;
    end

    for sl = 1:SL
        dataParams_tmp{sl} = dataParams;
        dataParams_tmp{sl}.images = images(:,:,sl,1,:);
    end

    % Create parallel pool
    parpool;

    % Perform magnitude fitting
    parfor (sl = 1:SL, mnw)
    % for sl = 1:SL
        fprintf('\n\nInitialize magnitude fitting for slice %i...', sl); tic;
        outparams_mag = fw_i2cm1i_3pluspoint_hernando_graphcut(dataParams_tmp{sl}, algoParams_mag, verboseFLAG);
        tst = toc; fprintf('Done (%.2f sec)', tst);

        % Extract data
        R2Smagn{sl} = outparams_mag.r2starmap;
        FMmagn{sl} = outparams_mag.fieldmap;
        Fmagn{sl} = outparams_mag.species(2).amps;
        Wmagn{sl} = outparams_mag.species(1).amps;
    end

    % Perform complex fitting
    if mixedFLAG
        parfor (sl = 1:SL, mnw)
            % for sl = 1:SL
            fprintf('\n\nInitialize mixed fitting for slice %i...', sl); tic;
            outparams_comp = fw_i2xm1c_3pluspoint_hernando_mixedfit(dataParams_tmp{sl}, algoParams_mix_tmp, verboseFLAG);
            tst = toc; fprintf('Done (%.2f sec)', tst);

            % Extract data
            R2Scomp{sl} = outparams_comp.r2starmap;
            FMcomp{sl} = outparams_comp.fieldmap;
            Fcomp{sl} = outparams_comp.species(2).amps;
            Wcomp{sl} = outparams_comp.species(1).amps;
        end
    end

    % Delete parallel pool
    delete(gcp('nocreate'))

    % Revert to arrays
    magn = struct('FM', zeros(size(images,1), size(images,2), SL), ...
                   'F', zeros(size(images,1), size(images,2), SL), ...
                   'W', zeros(size(images,1), size(images,2), SL), ...
             'R2starM', zeros(size(images,1), size(images,2), SL));

    for sl = 1:SL
        magn.R2starM(:,:,sl) = R2Smagn{sl};
        magn.FM(:,:,sl) = FMmagn{sl};
        magn.F(:,:,sl) = Fmagn{sl};
        magn.W(:,:,sl) = Wmagn{sl};
    end

    if mixedFLAG
        comp = magn;
        for sl = 1:SL
            comp.R2starM(:,:,sl) = R2Scomp{sl};
            comp.FM(:,:,sl) = FMcomp{sl};
            comp.F(:,:,sl) = Fcomp{sl};
            comp.W(:,:,sl) = Wcomp{sl};
        end
    else
        comp = [];
    end
end