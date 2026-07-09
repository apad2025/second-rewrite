function D = select_slices(D, sliceRange)
%SELECT_SLICES Keep only a subset of slices (z-dimension) in the data struct.
%
%   D = select_slices(D, sliceRange)
%
%   Inputs
%       D          : pipeline data struct (crops every D.Data.* volume whose
%                    3rd dimension matches the current slice count)
%       sliceRange : slice indices to keep, e.g. 1:2, [10 20 30], or a logical
%                    mask. Use [] (empty) to keep all slices.
%
%   Applied right after loading so all downstream stages (masking, bipolar
%   correction, graph-cut separation) run only on the selected slices -- handy
%   for quickly testing on a couple of slices before processing all 50.

    if nargin < 2 || isempty(sliceRange)
        return   % keep everything
    end

    nz = D.Size(3);
    if islogical(sliceRange)
        sliceRange = find(sliceRange);
    end
    sliceRange = sliceRange(:).';

    if any(sliceRange < 1 | sliceRange > nz) || any(mod(sliceRange,1) ~= 0)
        error('select_slices:range', ...
              'sliceRange must be integer indices within 1..%d.', nz);
    end

    % Crop every Data field that spans the slice dimension
    fn = fieldnames(D.Data);
    for k = 1:numel(fn)
        A = D.Data.(fn{k});
        if size(A,3) == nz
            D.Data.(fn{k}) = A(:,:,sliceRange,:,:);
        end
    end

    D.Size = size(D.Data.Image);
    fprintf('select_slices: keeping %d of %d slices (%s).\n', ...
            numel(sliceRange), nz, mat2str(sliceRange));
end
