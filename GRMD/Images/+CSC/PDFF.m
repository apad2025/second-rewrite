% Calculate fat fraction
function ff = PDFF(w, f, m) 
% Inputs
%   w: water
%   f: fat
%   m: mask, optional
% 
% Outputs
%  ff: fat fraction

    ff = abs(f)./(abs(f)+abs(w));

    if nargin == 3
        ff = ff.*m;
    end

    ff(isnan(ff)) = 0;
end