% Keep only D.Data, and reduce every volume in it to a single slice,
% e.g. 358x358x50 -> 358x358 and 358x358x50x6 (Image) -> 358x358x1x6
clear; clc; close all

slice = 12;

path = "";

% The file may still hold the original D, or an already-extracted Data
S = load(path);
if isfield(S, 'D')
    Data = S.D.Data;
else
    Data = S.Data;
end
clear S

% Trim each field along the slice (3rd) dimension, keeping any trailing
% dimensions such as the echoes in Image
fnames = fieldnames(Data);
for i = 1:numel(fnames)
    v = Data.(fnames{i});
    if ndims(v) >= 3 && size(v, 3) >= slice
        Data.(fnames{i}) = v(:, :, slice, :);
    end
end

save(path, 'Data');
