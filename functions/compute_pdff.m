function ff = compute_pdff(w, f, m)
%COMPUTE_PDFF Proton density fat fraction from water and fat maps.
%
%   ff = compute_pdff(w, f)        % unmasked
%   ff = compute_pdff(w, f, m)     % masked by m
%
%   Inputs
%       w : water map (complex or magnitude)
%       f : fat map   (complex or magnitude)
%       m : optional binary mask
%
%   Output
%       ff : fat fraction in [0, 1]   (ff = |f| / (|f| + |w|))
%
%   Extracted from PDFFCalc in the original DogAnalysis.m pipeline.

    ff = abs(f)./(abs(f)+abs(w));

    if nargin == 3
        ff = ff.*m;
    end

    ff(isnan(ff)) = 0;
end
