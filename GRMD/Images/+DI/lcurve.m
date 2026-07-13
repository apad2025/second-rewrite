% L-curve
function [lambda, Kappa, cost_data, cost_reg] = lcurve(D, method, params, plrange)

    % Set method
    flags.method = method;
    flags.verbose = false;
    
    % Check if Bouwman
    if ~strcmp(method, 'Magnitude Weighted L1')
        if isfield(params, 'alpha') && ~isfield(params, 'lambda')
            lambda = params.alpha;
        else
            lambda = params.lambda;
        end
    
        % Preallocate data
        cost_data = zeros(1, length(lambda));
        cost_reg = zeros(1, length(lambda));
    
        % Move all data to gpu so it doesn't have to be sent back and forth
        gpuFLAG = false;
        if isfield(params, 'gpuFLAG') && params.gpuFLAG
            nfields = fieldnames(params);
            for i = 1:numel(nfields)
                params.(nfields{i}) = gpuArray(params.(nfields{i}));
            end
            nfields = {'deltaTE', 'VoxelSize', 'F0', 'B0Direction'};
            for i = 1:numel(nfields)
                D.(nfields{i}) = gpuArray(D.(nfields{i}));
            end
            nfields = {'TotalField', 'LocalField', 'NoiseSTD', 'WeightedMagnitude', 'Mask'};
            for i = 1:numel(nfields)
                D.Data.(nfields{i}) = gpuArray(D.Data.(nfields{i}));
            end
            params.gpuFLAG = false;
            gpuFLAG = true;
        end

        % Iterate through lambdas
        for i = 1:length(lambda)
            % Isolate current lambda
            if isfield(params, 'alpha') && ~isfield(params, 'lambda')
                params.alpha = lambda(i);
            else
                params.lambda = lambda(i);
            end
    
            % Solve for QSM
            fprintf('\nIter %i/%i, lambda = %.3g', i-1, length(lambda)-2, lambda(i))
    
            D = dipoleinversion(D, flags, params);
    
            % Extract cost values
            cost_data(i) = D.DipoleInversion.DataConsistencyCost;
            cost_reg(i) = D.DipoleInversion.RegularizationCost;
        end
        flags = rmfield(flags, 'verbose');

        % Gather GPU data
        if gpuFLAG
            nfields = fieldnames(params);
            for i = 1:numel(nfields)
                params.(nfields{i}) = gather(params.(nfields{i}));
            end
            nfields = {'deltaTE', 'VoxelSize', 'F0', 'B0Direction'};
            for i = 1:numel(nfields)
                D.(nfields{i}) = gather(D.(nfields{i}));
            end
            nfields = {'TotalField', 'LocalField', 'NoiseSTD', 'WeightedMagnitude', 'Mask'};
            for i = 1:numel(nfields)
                D.Data.(nfields{i}) = gather(D.Data.(nfields{i}));
            end
            params.gpuFLAG = true;
            cost_data = gather(cost_data);
            cost_reg = gather(cost_reg);
        end
    
        % Calculate curvature
        Kappa = calc_curv_spline(lambda, cost_reg, cost_data, false);
    
        % Remove the first & last points
        Kappa = Kappa(2:end-1);
        lambda = lambda(2:end-1);
        cost_data = cost_data(2:end-1);
        cost_reg = cost_reg(2:end-1);
    
        % Filter curve data
        Kappaf = medfilt2(Kappa, [1 5]);
    
        % Find max and zero curvature points
        optmax = find(Kappaf == max(Kappaf),1); % Should always be on right-hand side of graph
        optzero = find(abs(Kappaf) == min(abs(Kappaf)),1); % should always be after max-curvature
        % optmax = find(Kappa == max(Kappa(1:end-round(length(Kappa)*2/3))),1); % Should always be on right-hand side of graph
        % optzero = find(abs(Kappa) == min(abs(Kappa(optmax:end))),1); % should always be after max-curvature

        % Now find max & zero curvature points
        params_max = params;
        params_zero = params;
    
        % Isolate current lambda
        if isfield(params, 'alpha') && ~isfield(params, 'lambda')
            params_max.alpha = lambda(optmax);
            params_zero.alpha = lambda(optzero);
        else
            params_max.lambda = lambda(optmax);
            params_zero.lambda = lambda(optzero);
        end
    
        % Obtain QSM
        data_max = dipoleinversion(D, flags, params_max);
        data_zero = dipoleinversion(D, flags, params_zero);
    
    else
        % Ensure B0 direction is positive
        [D, ~] = Operations.Flip(D);

        % Convert to rad from Hz (rad = Hz*2pi*sec)
        lfs = D.Data.LocalField.*(2*pi*D.deltaTE);

        % Pad data
        xpad = 32;
        ypad = 32;
        zpad = 8;
        lfs = padarray(lfs, [xpad, ypad, zpad], 0, 'both');
        msk = padarray(D.Data.Mask, [xpad, ypad, zpad], 0, 'both');
        mag = padarray(D.Data.WeightedMagnitude, [xpad, ypad, zpad], 0, 'both');

        % Find optimal L2
        disp('Finding the optimal L2-parameter:');
        [chiL2, lambda.L2] = reconstructSusceptibility(lfs, msk, D.VoxelSize, params.lambda_L2, D.B0Direction); % ppm

        % Find optimal L1
        disp('L1-parameter sweep (pairwise calculation of Lcurve-nodes)');
        lambda.L1 = L1_sweep(lfs, msk, D.VoxelSize, params.lambda_L1, lambda.L2, D.B0Direction);

        % The final magnitude weighted reconstruction
        disp('The final magnitude weighted reconstruction');
        chiL1 = L1_magnWeighted_reconstruction(lfs, mag, msk, D.VoxelSize, lambda.L1, lambda.L2);

        % Remove padding
        chiL1 = chiL1(xpad+1:end-xpad,ypad+1:end-ypad,zpad+1:end-zpad);
        chiL2 = chiL2(xpad+1:end-xpad,ypad+1:end-ypad,zpad+1:end-zpad);

        Kappa = []; cost_data = []; cost_reg = [];
    end

    % Display results
    if ~strcmp(method, 'Magnitude Weighted L1')
        % First plot the L-curve in the linear domain
        figure; tiledlayout(1,3,"TileSpacing","compact", "Padding","compact");
        nexttile;
        plot(cost_data,         cost_reg,           'LineStyle', '-',    'Marker', '.'); axis tight; hold on;
        plot(cost_data(optmax), cost_reg(optmax),   'LineStyle', 'none', 'Marker', 'o','Color',[0.301 0.745 0.933],'LineWidth',2, 'MarkerSize',10);
        plot(cost_data(optzero),cost_reg(optzero),  'LineStyle', 'none', 'Marker', 'o','Color',[0.850 0.325 0.098],'LineWidth',2, 'MarkerSize',10); hold off
        set(gcf,'Color','white')
        xlabel('Fidelity cost'); ylabel('Regularization cost');
        legend('','Max-curvature','Zero-curvature'); legend('Location','southeast')
    
        % Now plot the L-curve in the log domain
        nexttile;
        plot(log(cost_data),            log(cost_reg),          'LineStyle', '-',    'Marker', '.'); axis tight; hold on
        plot(log(cost_data(optmax)),    log(cost_reg(optmax)),  'LineStyle', 'none', 'Marker', 'o','Color',[0.301 0.745 0.933],'LineWidth',2, 'MarkerSize',10);
        plot(log(cost_data(optzero)),   log(cost_reg(optzero)), 'LineStyle', 'none', 'Marker', 'o','Color',[0.850 0.325 0.098],'LineWidth',2, 'MarkerSize',10); hold off
        set(gcf,'Color','white');
        xlabel('Log(Fidelity cost)'); ylabel('Log(Regularization cost)');
    
        % Plot the curvature as function of the regularization weight
        nexttile;
        semilogx(lambda,         Kappaf,          'LineStyle', '-',    'Marker', '.'); axis tight; hold on;
        semilogx(lambda,         Kappa,           'LineStyle', '--');
        semilogx(lambda(optmax), Kappaf(optmax),  'LineStyle', 'none', 'Marker', 'o','Color',[0.301 0.745 0.933],'LineWidth',2, 'MarkerSize',10);
        semilogx(lambda(optzero),Kappaf(optzero), 'LineStyle', 'none', 'Marker', 'o','Color',[0.850 0.325 0.098],'LineWidth',2, 'MarkerSize',10); hold off
        set(gcf,'Color','white');
        xlabel('Regularization weight'); ylabel('Curvature');
        legend('Filtered','Raw','',''); legend('Location','northeast')

        if isfield(params, 'alpha') && ~isfield(params, 'lambda')
            pltitle_max = ['QSM, \alpha = ' num2str(params_max.alpha)];
            pltitle_zero = ['QSM, \alpha = ' num2str(params_zero.alpha)];
        else
            pltitle_max = ['QSM, \lambda = ' num2str(params_max.lambda)];
            pltitle_zero = ['QSM, \lambda = ' num2str(params_zero.lambda)];
        end
        plotmygraph(real(data_max.Data.SusceptibilityMap), 'PlotTitle', [pltitle_max, ' (Max-Curvature)'], 'ColorbarTitle', 'Magnetic Susceptibility (ppm)', 'DataRange', plrange);
        plotmygraph(real(data_zero.Data.SusceptibilityMap), 'PlotTitle', [pltitle_zero, ' (Zero-Curvature)'], 'ColorbarTitle', 'Magnetic Susceptibility (ppm)', 'DataRange', plrange);
    else
        % Extract MIPs for viewing
        chiL2_tra =        squeeze(        max(chiL2(:, :, round(D.Size(3)/2)+(-2:2)),[],3));
        chiL1_tra =        squeeze(        max(chiL1(:, :, round(D.Size(3)/2)+(-2:2)),[],3));
        chiL2_cor = flipud(squeeze(permute(max(chiL2(round(D.Size(1)/2)+(-2:2), :, :),[],1),[2,3,1])));
        chiL1_cor = flipud(squeeze(permute(max(chiL1(round(D.Size(1)/2)+(-2:2), :, :),[],1),[2,3,1])));
        chiL2_sag = flipud(squeeze(permute(max(chiL2(:, round(D.Size(2)/2)+(-2:2), :),[],2),[3,1,2])));
        chiL1_sag = flipud(squeeze(permute(max(chiL1(:, round(D.Size(2)/2)+(-2:2), :),[],2),[3,1,2])));

        % Determine colorbar range
        cmin_L2 = min([min(chiL2_tra,[],'all'), min(chiL2_cor,[],'all'), min(chiL2_sag,[],'all')]);
        cmin_L1 = min([min(chiL1_tra,[],'all'), min(chiL1_cor,[],'all'), min(chiL1_sag,[],'all')]);
        cmax_L2 = max([max(chiL2_tra,[],'all'), max(chiL2_cor,[],'all'), max(chiL2_sag,[],'all')]);
        cmax_L1 = max([max(chiL1_tra,[],'all'), max(chiL1_cor,[],'all'), max(chiL1_sag,[],'all')]);

        figureFULL; tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        nexttile; imagesc(chiL2_tra); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); ylabel('L2 QSM'); 
        pbaspect([D.Size(1)*D.VoxelSize(1), D.Size(2)*D.VoxelSize(2), 1]);
        colormap gray; clim([cmin_L2, cmax_L2]); 
        
        nexttile; imagesc(chiL2_cor); axis off;
        pbaspect([D.Size(2)*D.VoxelSize(2), D.Size(3)*D.VoxelSize(3), 1])
        colormap gray; clim([cmin_L2, cmax_L2]); 

        nexttile; imagesc(chiL2_sag); axis off;
        pbaspect([D.Size(1)*D.VoxelSize(1), D.Size(3)*D.VoxelSize(3), 1])
        colormap gray; clim([cmin_L2, cmax_L2]); cb = colorbar; cb.Label.String = 'Susceptibility (ppm)';

        nexttile; imagesc(chiL1_tra); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); ylabel('L1 Mag Weighted QSM');
        pbaspect([D.Size(1)*D.VoxelSize(1), D.Size(2)*D.VoxelSize(2), 1]);
        colormap gray; clim([cmin_L1, cmax_L1]); 
        
        nexttile; imagesc(chiL1_cor); axis off;
        pbaspect([D.Size(2)*D.VoxelSize(2), D.Size(3)*D.VoxelSize(3), 1])
        colormap gray; clim([cmin_L1, cmax_L1]); 

        nexttile; imagesc(chiL1_sag); axis off;
        pbaspect([D.Size(1)*D.VoxelSize(1), D.Size(3)*D.VoxelSize(3), 1])
        colormap gray; clim([cmin_L1, cmax_L1]); cb = colorbar; cb.Label.String = 'Susceptibility (ppm)';
    end
end