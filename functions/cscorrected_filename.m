function name = cscorrected_filename(dicomPath)
%CSCORRECTED_FILENAME  Build "<date>_<series>_CScorrected.mat" from a DICOM path.
%   Given a DICOM series folder such as
%       ...\DICOM\20240709\GRE2D_FATWATER_WAYLON_0012
%   returns
%       20240709_GRE2D_FATWATER_WAYLON_CScorrected.mat
%   i.e. the date parent folder, the series folder with its trailing
%   "_<number>" stripped, and the CScorrected suffix. Matches the naming of
%   the pre-refactor exports. If the parent folder is not an 6-8 digit date,
%   the date prefix is omitted.
dicomPath = char(dicomPath);
dicomPath = regexprep(dicomPath, '[\\/]+$', '');        % strip trailing slash
[parent, series] = fileparts(dicomPath);                % series = GRE2D_FATWATER_WAYLON_0012
[~, dateStr] = fileparts(parent);                       % dateStr = 20240709
series = regexprep(series, '_\d+$', '');                % drop trailing _0012
if ~isempty(regexp(dateStr, '^\d{6,8}$', 'once'))
    name = sprintf('%s_%s_CScorrected.mat', dateStr, series);
else
    name = sprintf('%s_CScorrected.mat', series);
end
end
