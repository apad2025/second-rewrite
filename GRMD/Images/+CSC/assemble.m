% Stack slice-wise chemical shift corrections into a single dataset
function D = assemble(D, pth_slices, sname, verbose)
% Inputs:
%             D: data structure (preprocessed, i.e. the input to CSC.performSlice)
%    pth_slices: folder holding the per-slice saves
%         sname: base name of the per-slice saves
%       verbose: display flag
%
% Outputs:
%             D: data structure, matching what CSC.perform would have returned

    if nargin < 4, verbose = true; end
    if verbose; fprintf('\nAssembling slice-wise chemical shift corrections...\n'); end

    % Check that every slice was processed
    missing = [];
    for sl = 1:D.Size(3)
        if ~isfile(fullfile(pth_slices, sprintf('%s_%03i.mat', sname, sl)))
            missing(end+1) = sl; %#ok<AGROW>
        end
    end
    if ~isempty(missing)
        error('Missing %i of %i slices. Resubmit the array job with --array=%s', numel(missing), D.Size(3), strjoin(string(missing), ','))
    end

    % Combine slices
    for sl = 1:D.Size(3)
        % Loaded without a variable list so that saves made before the solver
        % statistics were recorded (Data & TE only) do not warn
        S = load(fullfile(pth_slices, sprintf('%s_%03i.mat', sname, sl)));

        % Preallocate from the shape of the first slice
        if sl == 1
            fnames = fieldnames(S.Data);
            Data = struct();
            for i = 1:numel(fnames)
                Data.(fnames{i}) = zeros([size(S.Data.(fnames{i}), [1 2]), D.Size(3), size(S.Data.(fnames{i}), 4)], 'like', S.Data.(fnames{i}));
            end

            % Solver statistics, one entry per slice (absent in saves made
            % before the counts were recorded)
            if isfield(S, 'Stats')
                snames_stats = fieldnames(S.Stats);
                Stats = struct();
                for i = 1:numel(snames_stats)
                    Stats.(snames_stats{i}) = nan(D.Size(3), 1);
                end
            else
                snames_stats = {};
            end

            % Echo times may have been trimmed during the correction
            D.TE = S.TE;
            D.deltaTE = mean(diff(D.TE));
        end

        for i = 1:numel(fnames)
            Data.(fnames{i})(:,:,sl,:) = S.Data.(fnames{i});
        end

        for i = 1:numel(snames_stats)
            Stats.(snames_stats{i})(sl) = S.Stats.(snames_stats{i});
        end
    end
    clear S

    % Combine into main structure
    for i = 1:numel(fnames)
        D.Data.(fnames{i}) = Data.(fnames{i});
    end
    clear Data

    % Record the per-slice iteration counts alongside the other correction
    % parameters, matching the field names CSC.perform uses for a whole volume
    if ~isempty(snames_stats)
        D.CSCorrection.GraphCuts.BipolarIterations = Stats.BipolarIterations;
        D.CSCorrection.GraphCuts.Iterations = Stats.Iterations;
        if isfield(Stats, 'OTIterations')
            D.CSCorrection.OptimizationTransfer.Iterations = Stats.OTIterations;
        end
    end
    clear Stats

    % Update data size
    D.Size = size(D.Data.TotalField);

    % Update flags
    D.Flags.CorrectedChemicalShift = true;
    D.Flags.UnwrappedPhase = true;
    D.CSCorrection.FatWaterSwapChecked.Automatic = false;
    D.CSCorrection.FatWaterSwapChecked.Manual = false;
end
