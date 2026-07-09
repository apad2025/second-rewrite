function revstr = UpdatePercent(perc, revstr)
%UPDATEPERCENT display progress bar in command window
% 
%   revstr = UpdatePercent(perc, revstr)
% 
%   Inputs
%          perc: current progress in percentage format
%        revstr: revised string, only used to remove previous percentage &
%                replace with new percentage. Should be initialized as 
%                empty char outside for loop.
% 
%   Outputs
%        revstr: revised string
%
% Jacob Degitz, Texas A&M University
% Created 11/17/2024
% Last edited 3/12/2026

    arguments
        perc {mustBeNumeric}
        revstr {mustBeText}
    end

    msg = sprintf('%.2f', perc);
    fprintf([revstr, msg, '%%']);
    revstr = repmat(sprintf('\b'), 1, length(msg)+1);
end