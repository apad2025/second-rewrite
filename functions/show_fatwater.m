function D = show_fatwater(src, opts)
%SHOW_FATWATER Display the maps in a saved fat/water result (.mat) or D struct.
%
%   show_fatwater()            % file picker for a FatWater_*.mat
%   show_fatwater(matFile)     % load and display that .mat
%   show_fatwater(D)           % display an in-memory D struct
%   D = show_fatwater(...)     % also return the loaded D
%
%   Displays whichever maps are present in D.Data:
%       Water magnitude, Fat magnitude, PDFF (|f|/(|f|+|w|)),
%       Field map (Hz), R2* map (s^-1), and -- for the vectorized Bipolar_GC
%       backend -- the bipolar Phase (PhiMap) and Amplitude (EpsMap) modulation
%       maps. Works on output from either hernando_fat_water or
%       vectorized_fat_water.
%
%   opts fields (all optional):
%       .slices : slice indices to show (default: the slices that actually
%                 contain data, so a run over e.g. 25:26 shows just those).
%       .maps   : cellstr subset of {'water','fat','pdff','field','r2star',
%                 'phi','eps'} to restrict which maps are shown.
%
%   Uses the project's plotmygraph, so each figure is a scrollable/zoomable
%   tiled montage.

    if nargin < 2, opts = struct(); end

    % ---- Resolve the source into a D struct ----
    if nargin < 1 || isempty(src)
        [f, p] = uigetfile('*.mat', 'Select a FatWater result .mat');
        if isequal(f, 0), D = []; return; end
        src = fullfile(p, f);
    end

    if ischar(src) || isstring(src)
        matFile = char(src);
        S = load(matFile);
        if isfield(S, 'D')
            D = S.D;
        else
            error('show_fatwater:noDstruct', ...
                  'No variable ''D'' found in %s.', matFile);
        end
        fprintf('Loaded %s\n', matFile);
    elseif isstruct(src)
        D = src;
    else
        error('show_fatwater:badInput', ...
              'Input must be a .mat path or a D struct.');
    end

    if ~isfield(D, 'Data') || ~isfield(D.Data, 'Water') || ~isfield(D.Data, 'Fat')
        error('show_fatwater:noMaps', ...
              'D has no fat/water maps (D.Data.Water / D.Data.Fat missing).');
    end

    % ---- Report how the file was produced ----
    if isfield(D, 'CSCorrection') && isfield(D.CSCorrection, 'Method')
        fprintf('Separation method: %s\n', D.CSCorrection.Method);
    end

    W = D.Data.Water;
    F = D.Data.Fat;
    [nx, ny, nz] = size(W);

    % ---- Choose slices: default to those that actually contain data ----
    if isfield(opts, 'slices') && ~isempty(opts.slices)
        sl = opts.slices(:)';
    else
        energy = squeeze(sum(sum(abs(W) + abs(F), 1), 2));
        sl = find(energy > 0)';
        if isempty(sl), sl = 1:nz; end   % fall back to all slices
    end
    plr = {1:nx, 1:ny, sl, 1};
    if numel(sl) == 1
        fprintf('Showing slice %d.\n', sl);
    else
        fprintf('Showing slices %s.\n', mat2str(sl));
    end

    % ---- Which maps to show ----
    if isfield(opts, 'maps') && ~isempty(opts.maps)
        want = lower(string(opts.maps));
    else
        want = ["water" "fat" "pdff" "field" "r2star" "phi" "eps"];
    end
    show = @(name) any(want == name);

    % ---- PDFF (from water/fat, masked if a mask is available) ----
    if show("pdff")
        if isfield(D.Data, 'Mask')
            pdff = compute_pdff(W, F, D.Data.Mask);
        else
            pdff = compute_pdff(W, F);
        end
    end

    % ---- Draw ----
    if show("water")
        plotmygraph(abs(W), 'PlotTitle', 'Water Map', ...
                    'ColorbarTitle', '[a.u.]', 'DataRange', plr);
    end
    if show("fat")
        plotmygraph(abs(F), 'PlotTitle', 'Fat Map', ...
                    'ColorbarTitle', '[a.u.]', 'DataRange', plr);
    end
    if show("pdff")
        plotmygraph(pdff, 'PlotTitle', 'PDFF Map', 'Colormap', 'hot', ...
                    'ColormapLimits', [0 1], 'ColorbarTitle', 'Fat fraction', ...
                    'DataRange', plr);
    end
    if show("field") && isfield(D.Data, 'TotalField')
        plotmygraph(D.Data.TotalField, 'PlotTitle', 'Field Map', ...
                    'ColorbarTitle', 'Frequency (Hz)', 'DataRange', plr);
    end
    if show("r2star") && isfield(D.Data, 'R2StarMap')
        plotmygraph(D.Data.R2StarMap, 'PlotTitle', 'R2* Map', ...
                    'ColorbarTitle', 's^{-1}', 'DataRange', plr);
    end

    % ---- Bipolar_GC-only correction maps ----
    if show("phi") && isfield(D.Data, 'PhiMap')
        plotmygraph(D.Data.PhiMap, 'PlotTitle', 'Bipolar Phase Modulation', ...
                    'ColorbarTitle', '[rad]', 'DataRange', plr);
    end
    if show("eps") && isfield(D.Data, 'EpsMap')
        plotmygraph(D.Data.EpsMap, 'PlotTitle', 'Bipolar Amplitude Modulation', ...
                    'ColorbarTitle', '[a.u.]', 'DataRange', plr);
    end

    if nargout == 0, clear D; end
end
