function [data, HDR, hdr] = loadima(dcm)
%LOADIMA Import DICOM *.IMA files & header data.
%
%   [data, info] = loadima(dcm)
% 
%   Inputs
%           dcm: folder or file containing DICOM files
% 
%   Outputs
%          data: image data
%          info: header data
%
% Jacob Degitz, Texas A&M University
% Created 10/17/2023
% Last edited on 12/18/2024

switch nargin
    case 0
        direct = pwd;

        % Obtain file
        dcm = uigetdir(direct);
end

% Ensure input is char
if ~ischar(dcm)
    try
        dcm = char(dcm);
    catch
        error('Input must be either a string or character array.')
    end
end
% Check if folder or file
if strcmp(char(dcm(end-3:end)), '.IMA')
    folderFLAG = false;
else
    folderFLAG = true;
end

% Import first file
if folderFLAG
    % List files within folder
    filelist = dir(dcm);

    % Skip any folders in main folder
    i=1;
    while i<=length(filelist)
        if filelist(i).isdir==1
            filelist = filelist([1:i-1 i+1:end]);
        else
            i=i+1;
        end
    end

    % Import info of first file
    hdr = dicominfo([dcm '/' filelist(1).name]);

else % one file
    data = dicomread(dcm);
    hdr = dicominfo(dcm);
end

if isfield(hdr, 'ImageType')
    imtype = strsplit(hdr.ImageType,'\');

    % Ignore original/derived and primary/secondary
    imtype = imtype(3:end);

    if strcmp(imtype{1},'RAWDATA')
        imtype = imtype{2};
    else
        imtype = imtype{1};
    end
end

%% Extract data
if folderFLAG
    % Check if complex
    complexFLAG = false;
    rmFLAG = false;
    if isfield(hdr, 'ScanningSequence')
        switch hdr.ScanningSequence
            case 'GR' % gradient echo
                % Check scan options
                if length(hdr.ScanOptions) > 3
                    if strcmp(hdr.ScanOptions, 'DIX')
                        switch hdr.ScanOptions(4)
                            case 'F'
                            case 'W'
                            case 'IN'
                            case 'OPP'
                        end
                    end
                else
                    if strcmpi(imtype, 'M') % magnitude data
                        dcm2 = [dcm(1:end-1) num2str(str2double(dcm(end))+1)];
                    elseif strcmpi(imtype, 'P') % phase data
                        dcm2 = [dcm(1:end-1) num2str(str2double(dcm(end))-1)];
                    end
            
                    % Validate folder exists & get files
                    if isfolder(dcm2)
                        filelist2 = dir(dcm2);
                        i=1;
                        while i<=length(filelist2)
                            if filelist2(i).isdir==1
                                filelist2 = filelist2([1:i-1 i+1:end]);
                            else
                                i=i+1;
                            end
                        end
    
                        % Check if folder is opposite data type
                        hdr_tmp = dicominfo([dcm2 '/' filelist2(1).name]);
    
                        if strcmpi(imtype, 'M') && strcmpi(hdr_tmp.ImageType(18), 'P')
                            filelist_p = filelist2; dcm_p = dcm2;
                            complexFLAG = true;
                        elseif strcmpi(imtype, 'P') && strcmpi(hdr_tmp.ImageType(18), 'M')
                            filelist_p = filelist; dcm_p = dcm;
                            filelist = filelist2; dcm = dcm2;
                            complexFLAG = true;
                        end
                        hdr = dicominfo([dcm '/' filelist(1).name]);
    
                        clear hdr_tmp filelist2
                    end
                    clear dcm2
                end

            case 'RM' % RF map, magnitude & phase in same folder
                complexFLAG = false;
                rmFLAG = true;
        end

    % In the event the dicom file has nothing
    elseif strcmp(imtype, 'ND')
        % There should be another folder where the actual data is contained
        dcm2 = [dcm(1:end-1) num2str(str2double(dcm(end))+1)];
        if exist(dcm2,"dir")
            dcm = dcm2;
            filelist = dir(dcm);
            i=1;
            while i<=length(filelist)
                if filelist(i).isdir==1
                    filelist = filelist([1:i-1 i+1:end]);
                else
                    i=i+1;
                end
            end
            hdr = dicominfo([dcm '/' filelist(1).name]);
            complexFLAG = false;
        else
            error('There is no data in this dicom!')
        end
    end

    % Set 3D flag
    d3FLAG = false;
    if strcmp(hdr.MRAcquisitionType,'3D')
        d3FLAG = true;
    end

    % Determine image size
    msize(1) = single(hdr.Height);
    msize(2) = single(hdr.Width);

    % Check for 3rd or 4th dimension
    if isfield(hdr, 'SpacingBetweenSlices') % 2D dataset
        msize(3) = length(filelist);
        msize(4) = 1;
    elseif isfield(hdr, 'CardiacNumberOfImages')
        if hdr.CardiacNumberOfImages > 1
            msize(3) = length(filelist);
            msize(4) = 1;
        else
            msize(3) = 1;
            msize(4) = 1;
        end
    else
        msize(3) = 1;
        msize(4) = 1;
    end

        % Special case for rf map
    if rmFLAG
        msize(4) = 2;
        msize(3) = msize(3)/2;
    end

    % if hdr(1).ImageType(end-3) Need to add in code to separate magnitude & phase data
    % Preallocate data
    data = zeros(msize(1), msize(2), msize(3), msize(4), 'uint16');
    if complexFLAG
        data_p = data;
        RSslp = zeros(msize(3), msize(4)); % rescale slopes
        RSint = zeros(msize(3), msize(4)); % rescale intercepts
    end

    % load all images in folder
    fprintf('\nLoading images... ');
    reverseStr = UpdatePercent((1/(msize(3)*msize(4)))*100, '');
    if d3FLAG
        for j = 1:msize(4)
            reverseStr = UpdatePercent((j/msize(4))*100, reverseStr);
            data(:,:,:,j) = dicomread([dcm '\' filelist(j).name]);
            hdr(j) = dicominfo([dcm '\' filelist(j).name]);
            if complexFLAG
                data_p(:,:,:,j) = dicomread([dcm_p '\' filelist_p(j).name]);
                hdr_tmp = dicominfo([dcm_p '\' filelist_p(j).name]);
                if isfield(hdr_tmp, 'RescaleSlope')
                    RSslp(:,j) = hdr_tmp.RescaleSlope;
                    RSint(:,j) = hdr_tmp.RescaleIntercept;
                else
                    RSslp(:,j) = 2;
                    RSint(:,j) = -4096;
                end
            end
        end
    else
        idx = 0;
        for i = 1:msize(3)
            for j = 1:msize(4)
                idx = idx + 1;
                reverseStr = UpdatePercent((idx/(msize(3)*msize(4)))*100, reverseStr);
                data(:,:,i,j) = dicomread([dcm '\' filelist(idx).name]);
                hdr(i,j) = dicominfo([dcm '\' filelist(idx).name]);
                if complexFLAG
                    data_p(:,:,i,j) = dicomread([dcm_p '\' filelist_p(idx).name]);
                    hdr_tmp = dicominfo([dcm_p '\' filelist_p(idx).name]);
                    if isfield(hdr_tmp, 'RescaleSlope')
                        RSslp(i,j) = hdr_tmp.RescaleSlope;
                        RSint(i,j) = hdr_tmp.RescaleIntercept;
                    else
                        RSslp(i,j) = 2;
                        RSint(i,j) = -4096;
                    end
                end
            end
        end
    end
    fprintf('\n');
end

%% Check for error in first line/column
for j = 1:msize(4)
    if d3FLAG
        for k = 1:3
            % Extract data
            switch k
                case 1
                    check = all(data(1,:,:,j) == 0);
                    data_tmp = data(1,:,:,j);
                    if complexFLAG
                        data_p_tmp = data_p(1,:,:,j);
                    end
                case 2
                    check = all(data(2:end,1,:,j) == 0);
                    data_tmp = data(:,1,:,j);
                    if complexFLAG
                        data_p_tmp = data_p(:,1,:,j);
                    end
                case 3
                    check = all(data(2:end,2:end,1,j) == 0);
                    data_tmp = data(:,:,1,j);
                    if complexFLAG
                        data_p_tmp = data_p(:,:,1,j);
                    end
            end

            % Estimate data
            if check
                data_tmp = randi([min(data(round(end*7/8):end,round(end*7/8):end,round(end*7/8):end,j),[],"all"),max(data(round(end*7/8):end,round(end*7/8):end,round(end*7/8):end,j),[],"all")], size(data_tmp), "uint16");
                % Do the same for phase data
                if complexFLAG
                    data_p_tmp = randi([min(data_p(2:end,2:end,2:end,j),[],"all"),max(data_p(2:end,2:end,2:end,j),[],"all")], size(data_p_tmp), "uint16");
                end

                % Reincorporate data
                switch k
                    case 1
                        data(1,:,:,j) = data_tmp;
                        if complexFLAG
                            data_p(1,:,:,j) = data_p_tmp;
                        end
                    case 2
                        data(:,1,:,j) = data_tmp;
                        if complexFLAG
                            data_p(:,1,:,j) = data_p_tmp;
                        end
                    case 3
                        data(:,:,1,j) = data_tmp;
                        if complexFLAG
                            data_p(:,:,1,j) = data_p_tmp;
                        end
                end
            end
        end
    else
        for i = 1:msize(3)
            for k = 1:2
                % Extract data
                switch k
                    case 1
                        check = all(data(1,:,i,j) == 0);
                        data_tmp = data(1,:,i,j);
                        if complexFLAG
                            data_p_tmp = data_p(1,:,i,j);
                        end
                    case 2
                        check = all(data(2:end,1,i,j) == 0);
                        data_tmp = data(:,1,i,j);
                        if complexFLAG
                            data_p_tmp = data_p(:,1,i,j);
                        end
                end

                % Estimate data
                if check
                    data_tmp = randi([min(data(round(end*7/8):end,round(end*7/8):end,i,j),[],"all"),max(data(round(end*7/8):end,round(end*7/8):end,i,j),[],"all")], size(data_tmp), "uint16");
                    % Do the same for phase data
                    if complexFLAG
                        data_p_tmp = randi([min(data_p(2:end,2:end,i,j),[],"all"),max(data_p(2:end,2:end,i,j),[],"all")], size(data_p_tmp), "uint16");
                    end

                    % Reincorporate data
                    switch k
                        case 1
                            data(1,:,i,j) = data_tmp;
                            if complexFLAG
                                data_p(1,:,i,j) = data_p_tmp;
                            end
                        case 2
                            data(:,1,i,j) = data_tmp;
                            if complexFLAG
                                data_p(:,1,i,j) = data_p_tmp;
                            end
                    end
                end
            end
        end
    end
end

%% Reshape data
% Convert to double
data = double(data);
if complexFLAG, data_p = double(data_p); end

% Check rescale, if applicable
if complexFLAG
    if all(RSslp==RSslp(1)) && all(RSint==RSint(1))
        RSslp = RSslp(1);
        RSint = RSint(1);
        data_p = (data_p.*RSslp + RSint)./max(data_p,[],'all').*pi; % [-pi, pi]
    else
        error('Rescale Slope/Intercept not the same across scan data')
    end
elseif rmFLAG
    % Correct phase rf maps
    data(:,:,:,2) = (data(:,:,:,2)-2048)*180/2048; % now in units of degrees
end

% Combine complex data
if complexFLAG
    % Account for magnitude points with intensity == 0
    % When combined with phase data, the phase points are set to 0 where
    % the magnitude is set to 0. To correct for this, add a very very small
    % number to the magnitude data.
    if min(data,[],'all') == 0
        data = data + 1/max(data,[],'all');
    end
    data = data.*exp(1i.*data_p);
end

data_old = data;

% Check for other dimensions
nTE = hdr(end).EchoNumbers;
nSl = size(data,3)/nTE;
if rmFLAG
    nTE = 2;
end
data = zeros(size(data,1), size(data,2), nSl, nTE);
for i = 1:nTE
    data(:,:,:,i) = data_old(:,:,((i-1)*size(data,3)+1):(i*size(data,3)));
end

% Re-order slices to proper location
sLoc = zeros(nSl,1);
for i = 1:nSl
    sLoc(i) = hdr(i,1).SliceLocation;
end
[~,sLocNew] = sort(sLoc);

data_old = data;
for i = 1:nSl
    data(:,:,sLocNew(i),:) = data_old(:,:,i,:);
end

%% Rename fields
% Info obtained from:   https://viewer.mathworks.com/?viewer=plain_code&url=https%3A%2F%2Fwww.mathworks.com%2Fmatlabcentral%2Fmlc-downloads%2Fdownloads%2Fe5a13851-4a80-11e4-9553-005056977bd0%2Fc5ce193c-31e2-4149-ac28-c6d14a9bd6f4%2Ffiles%2Fdicm_dict.m&embed=web

% Isolate field names
fnames = fieldnames(hdr(1));

% Iterate through field names
idx = 0;
fnames_private = {};
for i = 1:length(fnames)
    fname = fnames{i};

    % Check if private field
    if length(fname) > 7 && strcmp(fname(1:7), 'Private')
        idx = idx + 1;
        fnames_private{idx} = fname; %#ok<AGROW>
    end
end

% Now iterate through private fields 
HDR = hdr;
for i = 1:length(fnames_private)
    fname = fnames_private{i};

    % Match, extract & rename field names
    switch fname
        case 'Private_0019_100a'
            newname = 'NumberOfImagesInMosaic';
        case 'Private_0019_100b'
            newname = 'SliceMeasurementDuration';
        case 'Private_0019_100c'
            newname = 'B_value';
        case 'Private_0019_100d'
            newname = 'DiffusionDirectionality';
        case 'Private_0019_100e'
            newname = 'DiffusionGradientDirection';
        case 'Private_0019_100f'
            newname = 'GradientMode';
        case 'Private_0019_1011'
            newname = 'FlowCompensation';
        case 'Private_0019_1012'
            newname = 'TablePositionOrigin';
        case 'Private_0019_1013'
            newname = 'ImaAbsTablePosition';
        case 'Private_0019_1014'
            newname = 'ImaRelTablePosition';
        case 'Private_0019_1015'
            newname = 'SlicePosition_PCS';
        case 'Private_0019_1016'
            newname = 'TimeAfterStart';
        case 'Private_0019_1017'
            newname = 'SliceResolution';
        case 'Private_0019_1018'
            newname = 'RealDwellTime';
        case 'Private_0019_1027'
            newname = 'B_matrix';
        case 'Private_0019_1028'
            newname = 'BandwidthPerPixelPhaseEncode';
        case 'Private_0019_1029'
            newname = 'MosaicRefAcqTimes';
        case 'Private_0029_1008'
            newname = 'CSAImageHeaderType';
        case 'Private_0029_1009'
            newname = 'CSAImageHeaderVersion';
        case 'Private_0029_1010'
            newname = 'CSAImageHeaderInfo';
        case 'Private_0029_1018'
            newname = 'CSASeriesHeaderType';
        case 'Private_0029_1019'
            newname = 'CSASeriesHeaderVersion';
        case 'Private_0029_1020'
            newname = 'CSASeriesHeaderInfo';
        case 'Private_0029_1060'
            newname = 'SeriesWorkflowStatus';
        case 'Private_0051_100a'
            newname = 'TimeOfAcquisition';
        case 'Private_0051_100b'
            newname = 'AcquisitionMatrixText';
        case 'Private_0051_100c'
            newname = 'FieldOfViewText';
        case 'Private_0051_100d'
            newname = 'SlicePositionText';
        case 'Private_0051_100e'
            newname = 'ImageOrientationText';
        case 'Private_0051_100f'
            newname = 'CoilString';
        case 'Private_0051_1011'
            newname = 'ImaPATModeText';
        case 'Private_0051_1012'
            newname = 'TablePositionText';
        case 'Private_0051_1013'
            newname = 'PositivePCSDirections';
        case 'Private_0051_1016'
            newname = 'ImageTypeText';
        case 'Private_0051_1017'
            newname = 'SliceThicknessText';
        case 'Private_0051_1019'
            newname = 'ScanOptionsText';
    end

    if exist("newname", 'var')
        HDR = rmfield(HDR, fname);
        for j = 1:length(HDR)
            HDR(j).(newname) = hdr(j).(fname);
        end
        clear newname
    end
end

% Trim header file
fields = {'FileModDate', ...
    'FileSize', ...
    'ColorType', ...
    'FileMetaInformationGroupLength', ...
    'FileMetaInformationVersion', ...
    'SOPInstanceUID', ...
    'MediaStorageSOPClassUID', ...
    'MediaStorageSOPInstanceUID', ...
    'TransferSyntaxUID', ...
    'ImplementationClassUID', ...
    'SpecificCharacterSet', ...
    'SOPClassUID', ...
    'StudyDate', ...
    'SeriesDate', ...
    'AcquisitionDate', ...
    'StudyTime', ...
    'SeriesTime', ...
    'AcquisitionTime', ...
    'AccessionNumber', ...
    'InstitutionName', ...
    'InstitutionAddress', ...
    'ReferringPhysicianName', ...
    'StationName' ...
    'PerformingPhysicianName', ...
    'DeviceSerialNumber', ...
    'StudyInstanceUID' ...
    'SeriesInstanceUID', ...
    'FrameOfReferenceUID', ...
    'PerformedProcedureStepStartDate', ...
    'PerformedProcedureStepStartTime', ...
    'PerformedProcedureStepID', ...
    'PerformedProcedureStepDescription', ...
    'Format', ...
    'FormatVersion', ...
    'EchoTrainLength', ...
    'Modality', ...
    'ContentTime', ...
    'PatientComments', ...
    'ContentDate', ...
    'InstitutionalDepartmentName', ...
    'OperatorsName', ...
    'ReferencedImageSequence', ...
    'AngioFlag', ...
    'dBdt', ...
    'PhotometricInterpretation', ...
    'DataSetTrailingPadding', ...
    'StudyID', ...
    'BitsAllocated', ...
    'BitsStored', ...
    'HighBit', ...
    'PixelRepresentation', ...
    'SmallestImagePixelValue', ...
    'LargestImagePixelValue', ...
    'WindowCenterWidthExplanation', ...
    'Private_0019_10xx_Creator', ...
    'Private_0019_1008', ...
    'Private_0019_1009', ...
    'Private_0029_10xx_Creator', ...
    'Private_0029_11xx_Creator', ...
    'Private_0029_1160', ...
    'Private_0051_10xx_Creator', ...
    'Private_0051_1008', ...
    'Private_0051_1009'};

truefields = isfield(HDR, fields);
for f = 1:length(fields)
    if truefields(f)
        HDR = rmfield(HDR, fields{f});
    end
end
end