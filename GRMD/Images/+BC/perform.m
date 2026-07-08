% Correct for bipolar readout error
function data = perform(data, data_raw, flags)
    if flags.verbose; fprintf('\nCorrecting for bipolar phase error...'); tic; end
    switch flags.bipolarcorrection.method
        case 'MEDI'
            Mask = data.Data.Mask;
            Mask = imdilate(Mask, strel('disk', 4));
            data.Data.Image = iField_correction_new(data.Data.Image, data.VoxelSize, Mask, data.InPlanePhaseEncodingDirection);

        case 'SEPIAslow'
            [data.Data.Image, data.Data.FIT3D] = BipolarEddyCorrect_custom(data);

        case 'SEPIAfast'
            [data.Data.Image,~] = FastBipolarCorrect(data.Data.Image, data.Data.Mask);

        case 'hernando' % https://doi.org/10.1002/mrm.26228
            % Generate phase offset
            phi0 = 0; % degrees
            phi1 = 0; % degrees
            rowFLAG = false;
            if strcmp(data.InPlanePhaseEncodingDirection, 'ROW'), rowFLAG = true; end
            if rowFLAG
                x = repmat(linspace(-1, 1, data.Size(2)).', [1 data.Size(2) data.Size(3)]);
            else % COLUMN
                x = repmat(linspace(-1, 1, data.Size(1)), [data.Size(1) 1 data.Size(3)]);
            end

            % Apply phase offset
            % for ec = 2:2:D.Size(4), data.Data.Image(:,:,:,ec) = data_raw.Image(:,:,:,ec).*exp(1i*((phi0 + phi1*x*(ec/2))*(pi/180))); end
            data.Data.Image(:,:,:,2:2:end) = data_raw.Image(:,:,:,2:2:end).*exp(1i*((phi0 + phi1*x)*(pi/180)));

            % Create data & parameter structures
            dataParams = struct('FieldStrength', data.B0, 'TE', data.TE, 'PrecessionIsClockwise', 1);
            flags_tmp = flags; flags_tmp.cscorrection = struct('method', 'IGC', 'subsample', 2);
            algoParams = CSC.grabPars(flags_tmp);
            images = zeros(data.Size(1), data.Size(2), data.Size(3), 1, data.Size(4));
            images(:,:,:,1,:) = data.Data.Image;

            % Perform magnitude & mixed fitting
            if flags.zipped, mnw = 4; else, mnw = 6; end
            [magn, comp] = CSC.IGC(images, dataParams, algoParams, false, mnw, true);

            % Calculate RMSE between magnitude & complex fitting
            R2starq = magn.R2starM < 100;
            data.Data.compWq = comp.W.*R2starq; data.Data.compFq = comp.F.*R2starq;
            data.Data.magnWq = magn.W.*R2starq; data.Data.magnFq = magn.F.*R2starq;
            RMSE = sqrt(sum((abs(data.Data.compWq) - abs(data.Data.magnWq)).^2 + (abs(data.Data.compFq) - abs(data.Data.magnFq)).^2, 'all'));
    end

    % Update flags
    data.Flags.CorrectedBipolarPhase = true;
    data.BipolarCorrection = struct('Method', flags.bipolarcorrection.method);
    if strcmp(flags.bipolarcorrection.method, 'hernando'), data.BipolarCorrection.RMSE = RMSE; end

    if flags.verbose, tm = toc; fprintf('Done (%0.2f sec)', tm); end

    % Bipolar eddy correction
    function [Image, BipolarReadoutInducedPhase] = BipolarEddyCorrect_custom(data) 
        bipolarCplxME = data.Data.Image;
        bipolarCplxME(isnan(bipolarCplxME)) = 0;
    
        % successive phase difference mean made so that the last echo taken into
        % account has to be an odd number ()
        % to= dims(4)-(-mod(dims(4),2)+1);
        % Phase of mean evens minus mean odds - contains Eddy Current and
        to = data.Size(4)-(mod(data.Size(4),2));
        PhaseDiffEvensMinusOdds = angle(mean(bipolarCplxME(:,:,:,2:2:to)./bipolarCplxME(:,:,:,1:2:to-1),4));
        PhaseDiffEvensMinusOdds(isnan(PhaseDiffEvensMinusOdds)) = 0;
    
        PhaseDiffEvens  = angle(mean(bipolarCplxME(:,:,:,4:2:end)./bipolarCplxME(:,:,:,2:2:end-2),4));
        PhaseDiffOdds   = angle(mean(bipolarCplxME(:,:,:,3:2:end)./bipolarCplxME(:,:,:,1:2:end-2),4));
        PhaseDiffEvens(isnan(PhaseDiffEvens)) = 0;
        PhaseDiffOdds(isnan(PhaseDiffOdds)) = 0;
        
        PhaseDiff = angle(exp(1i*PhaseDiffEvens) + exp(1i*PhaseDiffOdds));
    
        %%% a more robust eddy current estimation has to do unwrapping
        magnitude   = double(mean(abs(bipolarCplxME(:,:,:,1:2:end-2)),4));
        fieldMap = estimateTotalField_custom(PhaseDiff, magnitude, data.Data.Mask, data.VoxelSize);
        
        % the last echo to be taken into account has to be even
        magnitude   = double(mean(abs(bipolarCplxME(:,:,:,1:2:to)),4));
        fieldMapOddEven = estimateTotalField_custom(PhaseDiffEvensMinusOdds, magnitude, data.Data.Mask, data.VoxelSize);
        
        meanEddycurrent2 = fieldMapOddEven - fieldMap;
        
        % the chosen fitting strategy takes the magnitude image as weights... and
        % nothing else fancy.. first order correction gave the best results
        meanEddycurrent2(isnan(meanEddycurrent2)) = 0;
    
        [FIT3D,~,~] = PolyFit(double(meanEddycurrent2),data.Data.Mask,1);
    
        Image = zeros(data.Size);
        for echo = 1:data.Size(4)
            Image(:,:,:,echo) = bipolarCplxME(:,:,:,echo).*exp(-1i *(-1)^echo * 0.5 * FIT3D);
        end
        
        BipolarReadoutInducedPhase = FIT3D.* data.Data.Mask;
    
        function TotalField = estimateTotalField_custom(fieldMap, magn, mask, voxelSize)
            % Input
            % --------------
            % fieldMap      : original single-/multi-echo wrapped phase image, in rad
            % mask          : signal mask
            % matrixSize    : size of the input image
            % voxelSize     : spatial resolution of each dimension of the data, in mm
            % algorParam    : structure contains fields with algorithm-specific parameter(s)
            % headerAndExtraData : structure contains extra header info/data for the algorithm
            %
            % Output
            % --------------
            % TotalField     : unwrapped total field, in rad
    
            disp('--------------------');
            disp('Total field recovery');
            disp('--------------------');
            disp('Calculating field map...');
    
            fprintf('Temporal phase unwrapping: Optimum Weights \n');
    
            TotalField = zeros(size(fieldMap));
            for sl = 1:size(fieldMap,3)
                fprintf('slice %i...\n', sl); 
                TotalField(:,:,sl) = unwrapping_gc(fieldMap(:,:,sl), magn(:,:,sl), voxelSize, 'no', 1);
            end
            TotalField = TotalField-round(mean(TotalField(mask == 1))/(2*pi))*2*pi;
    
            % Correct output
            TotalField(isnan(TotalField)) = 0;
            TotalField(isinf(TotalField)) = 0;
    
            disp('The resulting field map with the following unit: Hz' );
        end
    end
end