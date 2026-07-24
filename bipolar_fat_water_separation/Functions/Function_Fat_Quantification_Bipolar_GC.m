%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright Véronique Fortier 2020 - MIT License
%
% Fat quantification                                                      %
%                                                                 %
% Based on: Berglund J, Johansson L, Ahlström H, Kullberg J. Three-point  %
% Dixon method enables whole-body water and fat imaging of obese subjects.%
% Magn Reson Med. 2010, Jun;63(6):1659-1668.
% Changed by Jorge Campos based on: Liu, C.-Y., McKenzie, C.A., Yu, H., Brittain, J.H. and Reeder, S.B. (2007), Fat quantification with IDEAL gradient echo imaging: Correction of bias from T1 and noise. Magn. Reson. Med., 58: 354-364. https://doi.org/10.1002/mrm.21301
% ----------------------------------------------------------------------- %
% Input argument:
%   -Sf: complex fat signal
%   -Sw: complex water signal
%
% Output argument:
%   -ns is the fat fraction map (not in percentage)
%   -wf is the water fraction map (not in percentage)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ns, wf] = Function_Fat_Quantification_Bipolar_GC(Sf, Sw)

ns = zeros(size(Sf));

% Fat-dominant mask
fDom_mask = abs(Sf) > abs(Sw);

% Calculate fat fraction in fat-dominant pixels
ns(fDom_mask) = abs(Sf(fDom_mask))./abs(Sf(fDom_mask) + Sw(fDom_mask));

% Calculate fat fraction in water-dominant pixels
ns(~fDom_mask) = 1 - abs(Sw(~fDom_mask))./abs(Sf(~fDom_mask) + Sw(~fDom_mask));

ns(isnan(ns)) = 0;
wf = 1 - ns;

end