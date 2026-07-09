% PLOTMYGRAPH plots data based on specified parameters
%   PLOTMYGRAPH(data) plots the given data as 2D grayscale images within a 
%   tiled layout figure, when data is a 2D or 3D matrix of data points. If 
%   data is a gpuarray, it is loaded onto the CPU before plotting. If data 
%   is a 3D matrix, each image is labeled with it's index.
%
%   PLOTMYGRAPH(___, 'PlotTitle', pltitle, 'ImageTitles', imtitles, 
%       'ColorbarTitle', ctitle, 'ColormapLimits', clims, 'Colormap', cmap,
%       'IsotropicVoxel', false, 'Fullscreen', fscreenFLAG, 'DataRange', 
%       plrange) plots the given data with the following options:
%           * pltitle, character/string, is used as the main figure title. 
%             A shortened version of pltitle is used as the figure title.
%           * imtitles, either string vector or cell array of 
%             strings/characters with length equal to the number of
%             individual plots, adds a custom title to each 2D plot.
%           * ctitle, character/string, adds a colorbar for the main figure
%             using the given input as the colorbar title.
%           * clim, 1x2 numeric, forces the colormap limits
%           * cmap, character/string, uses the colormap specified by the
%             input, which must be a built-in colormap.
%           * isoFLAG, logical, forces each plot to be a square.
%           * fscreenFLAG, logical, forces the figure to initialize as a
%             fullscreen window.
%           * plrange, Nx1 cell array for a data matrix with N dimensions,
%             only plots the data indicated by the indices included in the
%             cell array.
% 
% Jacob Degitz, Texas A&M University
% Created 10/3/2023
% Last edited 03/13/2026

function fig = plotmygraph(data, varargin)

% Create cell array of colormaps
cmaps = cellstr(colormaplist);

% Parse inputs
p = inputParser;

addRequired(p, 'data', @(x) (isnumeric(x) || islogical(x)) && (ismatrix(x) || ndims(x)<=4))
addParameter(p, 'PlotTitle', ' ', @(x) (isStringScalar(x) || ischar(x)))
addParameter(p, 'ImageTitles', [], @(x) isvector(x) || iscellstr(x))
addParameter(p, 'ColorbarTitle', [], @(x) (isStringScalar(x) || ischar(x)))
addParameter(p, 'ColormapLimits', [], @(x) (isnumeric(x) && numel(x)==2))
addParameter(p, 'Colormap', 'gray', @(x) (isStringScalar(x) || ischar(x)) && any(strcmpi(x,cmaps)))
addParameter(p, 'DataRange', [], @(x) iscell(x))
addParameter(p, 'FullScreen', false, @(x) islogical(x) && isscalar(x))
addParameter(p, 'IsotropicVoxel', true, @(x) islogical(x) && isscalar(x))

parse(p, data, varargin{:})

% Further parse DataRange
data = p.Results.data;
if ~isempty(p.Results.DataRange)
    plrange = p.Results.DataRange;
    if length(plrange) > ndims(data)
        if ~isscalar(plrange{end})
            warning('Too many indices are present in DataRange. Excess indices will be trimmed.')
        end
        for i = length(plrange):-1:(ndims(data)+1)
            plrange(i) = [];
        end
    end

    if any(cellfun(@(x)any(~isnumeric(x))||any(isnan(x))||any(isinf(x)),plrange))
        error('DataRange values cannot exceed size of input data.')
    end
else
    % Set default plot range
    if isempty(p.Results.DataRange)
        plrange = cell(1,ndims(data));
        for i = 1:ndims(data)
            plrange{i} = 1:size(data,i);
        end
    
        if ndims(data)==4
            plrange{4} = 1;
        end
    else
        plrange = p.Results.DataRange;
    
        % Append extra indices
        if length(plrange) < ndims(data)
            for i = (length(plrange)+1):(ndims(data)-length(plrange))
                plrange{i} = 1;
            end
        end
    end
end

% Further parse image titles
if ~isempty(p.Results.ImageTitles)
    if length(p.Results.ImageTitles) == length(plrange{3})
        imtitles = p.Results.ImageTitles;
    else
        error('The number of ImageTitles must be equal to the number of 2D plots present.')
    end
elseif length(plrange) > 2
    imtitles = string(plrange{3});
else
    imtitles = "1";
end

%% Extract inputs
ctitle = char(p.Results.ColorbarTitle);
cmap = p.Results.Colormap;
pltitle = char(p.Results.PlotTitle);
fsFLAG = p.Results.FullScreen;

% Set empty figure title if no input plot title
if strcmp(pltitle, ' ')
    ftitle = [];
elseif strcmp(pltitle(end), ')')
    plt_split = split(pltitle(1:strfind(pltitle, '(')), ' ');
    plt_split = plt_split(1:end-1);
    if length(plt_split) < 2
        ftitle = [char(plt_split)];
    else
        ftitle = [char(plt_split(end-1)), char(plt_split(end))];
    end
else
    plt_split = split(pltitle, ' ');
    if size(plt_split, 1) > 1
        ftitle = [char(plt_split(end-1)), char(plt_split(end))];
    else
        ftitle = char(plt_split);
    end
end

if length(ftitle) > 10 % Shorten figure title, if necessary
    ftitle = ftitle(1:10);
end

% Set FOV type
if p.Results.IsotropicVoxel
    fovtype = 'image';
else
    fovtype = 'square';
end

%% Prepare data
% Gather any data from GPU
if isgpuarray(data); data = gather(data); end
if isgpuarray(plrange); plrange = gather(plrange); end

% Determine datasize & trim data
switch length(plrange)
    case 2
        % Set flag
        flag2D = true;

        % Trim data
        data = data(plrange{1}, plrange{2});
    case 3
        % Set flag
        flag2D = false;

        % Trim data
        data = data(plrange{1}, plrange{2}, plrange{3});
    case 4
        % Set flag
        flag2D = false;

        % Trim data
        data = squeeze(data(plrange{1}, plrange{2}, plrange{3}, plrange{4}));
end

% Obtain colormap bounds
clims = p.Results.ColormapLimits;
if isempty(clims)
    clims = [double(min(data(:))) double(max(data(:)))];
end

% Ensure character data type
pltitle = char(pltitle);

%% Create figure
if isempty(ftitle)
    if fsFLAG
        fig = figure('WindowState', 'fullscreen');
    else
        fig = figure;
    end
else
    if fsFLAG
        fig = figure('Name', ftitle, 'NumberTitle', 'off', 'WindowState', 'fullscreen');
    else
        fig = figure('Name', ftitle, 'NumberTitle', 'off');
    end
end

% Plot image(s)
if flag2D
    imagesc(data, clims); axis(fovtype, 'off'); title(pltitle); colormap(cmap);
    axtoolbar(gca, {'datacursor','zoomin','zoomout','restoreview'});
else
    % Create figure
    tl = tiledlayout('flow', 'TileSpacing', 'tight', 'Padding', 'tight');

    for k = 1:size(data,3)
        ax(k) = nexttile;
        imagesc(data(:,:,k), clims); title(imtitles(k)); colormap(cmap);
    end
    axis(ax, fovtype, 'off')
    axtoolbar(tl, {'datacursor','zoomin','zoomout','restoreview'});
    title(tl, pltitle);
end

% Add colorbar
if ~isempty(ctitle)
    cb = colorbar; 
    if ~flag2D
        cb.Layout.Tile = 'east'; 
    end
    cb.Label.String = ctitle;
end
end