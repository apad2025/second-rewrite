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
%   -wfat: complex fat signal
%   -wwater: complex water signal
%
% Output argument:
%   -MD_ff is the fat fraction map (not in percentage)
%   -MD_wf is the water fraction map (not in percentage)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ MD_ff, MD_wf ] = Function_Fat_Quantification_Bipolar_GC( wfat, wwater)

ff_init = abs(wfat) ./ ( abs(wfat) + abs(wwater) );

high_ff_mask = zeros(size(ff_init));
high_ff_mask(ff_init>0.5) = 1;

MD_ff = zeros(size(ff_init));
MD_wf = zeros(size(ff_init));

MD_ff(high_ff_mask==1) = abs(wfat(high_ff_mask==1)) ./ abs( wfat(high_ff_mask==1) + wwater(high_ff_mask==1) );
MD_ff(high_ff_mask==0) = 1 - abs(wwater(high_ff_mask==0)) ./ abs( wfat(high_ff_mask==0) + wwater(high_ff_mask==0) );

MD_ff(isnan(MD_ff)) = 0;

MD_wf = 1 - MD_ff;


end

