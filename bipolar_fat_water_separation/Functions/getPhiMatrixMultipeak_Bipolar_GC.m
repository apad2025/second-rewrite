% Function: getPhiMatrixMultipeak
% 
% Returns the Phi matrix of complex exponentials from the VARPRO
% formulation of Dixon imaging
% 
% Arguments:
%   - deltafs: chemical shifts of the different chemical species
%   - relAmps: relative amplitudes of different fat peaks
%   - t: echo times used
% 
% Returns: 
% 
%   - Phi: the Phi matrix, Phi_{i,j} =  exp(j 2 \pi deltaf_j t_i )
%
%
% Author: Diego Hernando
% Date created: Aug 27, 2008
% Last modified: Aug 27, 2008
%
function Phi = getPhiMatrixMultipeak_Bipolar_GC( deltafs,relAmps, t )

numte = length(t);

% Portion of the matrix for odd echoes

[DF,T_odd] = meshgrid( deltafs,t(1:1:numte/2) );
[A,T2_odd] = meshgrid( relAmps,t(1:1:numte/2) );


Phi1 = exp(j*2*pi*T_odd.*DF);
Phi_odd = [Phi1(:,1) , sum(Phi1(:,2:end).*A,2) , zeros(size(Phi1,1),2) ];

% Portion of the matrix for even echoes

[DF,T_even] = meshgrid( deltafs,t(numte/2+1:1:numte) );
[A,T2_even] = meshgrid( relAmps,t(numte/2+1:1:numte) );

Phi2 = exp(j*2*pi*T_even.*DF);
Phi_even = [zeros(size(Phi1,1),2) , Phi2(:,1) , sum(Phi2(:,2:end).*A,2) ];


Phi = [Phi_odd;Phi_even];

