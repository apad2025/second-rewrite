%% Select data & analysis pipeline
clear; clc; close all

% Dog 3, date 3 dims
% dims.x = [52,300];
% dims.y = [48,296];
% dims.z = [7,9];

% Dog 4, date 3 dims
% dims.x = [42,290];
% dims.y = [42,290];
% dims.z = [7,9];

% Dog 4, date 2 dims
% dims.x = [12,260];
% dims.y = [12,260];
% dims.z = [7,9];

% Main path
% pth_main = "G:\Shared drives\McDougall Grad Students\Projects Grants\GRMD R01\In Vivo Studies\Scanner data\3T\Dog Dates";
pth_main = "C:\Users\apad2\Desktop\Fat_water_separation\DICOM_Files";

% Select dog (high PDFF are 4,3 & 2,3) (3,3 isn't finished being processed)
% 1 - Waylon
% 2 - Sushi
% 3 - Selene
% 4 - Aphrodite
% 5 - EOS
dog = 4;
date = [2, 3];

% Select flags
flags = struct('trimmed', true, ...
                'zipped', true, ...
     'bipolarcorrection', struct('method', 'SEPIAfast'), ...
               'verbose', true, ...
                  'plot', true, ...
            'unwrapping', struct('method', 'RegGrow', ...
                                 'subsample', 1, ...
                                 'corrected', true), ...
          'cscorrection', struct('method', 'vlGC', ...
                              'subsample', 1), ...
             'bfremoval', struct('D2', struct('method','V-SHARP - SEPIA'), ...
                                 'D3', struct('method','V-SHARP - SEPIA')), ...
       'dipoleinversion', struct('method', 'MEDI'));
saveFLAG = true;
cscFLAG = true;

% Force 3D background field removal to be V-SHARP (STISuite) if that is used for 2D
if strcmp(flags.bfremoval.D2.method, 'V-SHARP - STISuite')
    flags.bfremoval.D3.method = 'V-SHARP - STISuite';
end

% Select Comparison type
if length(date) > 1
    comptype = 'longitudinal';
elseif length(dog) > 1
    comptype = 'crosssection';
end

%% Load data
switch comptype
    case 'longitudinal'
        dog = ones(length(date),1)*dog(1);
    case 'crosssection'
        date = ones(length(dog),1)*date;
end
D_raw_all = cell(length(date),1);
D_csc_all = cell(length(date),1);
D_sus_all = cell(length(date),1);
pth_all = cell(length(date),1);
snames_all = cell(length(date),1);
for i = 1:length(date)
    [snames_all{i}, DD, HDR, pth_all{i}, plrange, vout, filelocs] = DogInitialize(pth_main, dog(i), date(i), flags, cscFLAG);

    % Extract individual saves
    for j = 1:length(vout)
        D = vout{j};
        if filelocs.Preprocessed == j
            D_raw_all{i} = D;
        elseif filelocs.CorrectedCS == j
            D_csc_all{i} = D;
        elseif filelocs.SusceptibilityMap == j
            D_sus_all{i} = D;
        end
    end
    clear vout filelocs D
end

%% Construct PDFF (with eroded mask)
PDFF_all = cell(1,length(D_raw_all));
for i = 1:length(PDFF_all)
    PDFF_all{i} = (abs(D_csc_all{i}.Data.Fat)./(abs(D_csc_all{i}.Data.Fat)+abs(D_csc_all{i}.Data.Water)));
    PDFF_all{i} = PDFF_all{i}.*D_csc_all{i}.Data.Mask;
    % PDFF_all{i} = PDFF_all{i}.*imerode(D_csc_all{i}.Mask, strel('disk',1));
end

% Isolate region of interest
dims.x = [12,260; 42,290];
dims.y = [12,260; 42,290];
dims.z = [7,9; 7,9];
for i = 1:length(date)
    mag_all{i} = D_raw_all{i}.Data.WeightedMagnitude(dims.x(i,1):dims.x(i,2),dims.y(i,1):dims.y(i,2),dims.z(i,1):dims.z(i,2));
    PDFF_all{i} = PDFF_all{i}(dims.x(i,1):dims.x(i,2),dims.y(i,1):dims.y(i,2),dims.z(i,1):dims.z(i,2));
    QSM_all{i} = D_sus_all{i}.Data.SusceptibilityMap(dims.x(i,1):dims.x(i,2),dims.y(i,1):dims.y(i,2),dims.z(i,1):dims.z(i,2));

    switch comptype
        case 'longitudinal'
            titlenum = num2str(date(i));
        case 'crosssection'
            titlenum = num2str(dog(i));
    end
    plotmygraph(mag_all{i}, 'PlotTitle', ['Magnitude' titlenum], 'IsotropicVoxel',true);
    plotmygraph(PDFF_all{i}, 'PlotTitle', ['PDFF' titlenum], 'IsotropicVoxel',true);
    plotmygraph(QSM_all{i}, 'PlotTitle', ['QSM' titlenum], 'IsotropicVoxel',true);
end

%% Check for previous masks
% List files within folder
for i = 1:length(date)
    filelist = dir(pth_all{i});
    
    i_old = 1; loadFLAG = false;
    while i_old <= length(filelist)
        % Ensure not a folder & mask file
        if ~filelist(i_old).isdir
            if strcmp(filelist(i_old).name(1:3), 'ROI')
                fprintf('Loading saved data: %s\n', filelist(i_old).name);
                load(char(fullfile(pth_all{i}, filelist(i_old).name)));
                ROI{i} = roi;
                loadFLAG = true;
                break
            else
                i_old = i_old+1;
            end
        else
            i_old = i_old+1;
        end
    end
    
    % Check if no roi loaded
    if ~loadFLAG
        m_fatL = zeros(size(PDFF_all{1}));
        m_fatR = zeros(size(PDFF_all{1}));
        m_mixedL = zeros(size(PDFF_all{1}));
        m_mixedR = zeros(size(PDFF_all{1}));
        m_muscleL = zeros(size(PDFF_all{1}));
        m_muscleR = zeros(size(PDFF_all{1}));
        for j = 1:(dims.z(i,2)-dims.z(i,1)+1)
            plotmygraph(mag_all{i}(:,:,j), 'PlotTitle', 'Magnitude');
            plotmygraph(PDFF_all{i}(:,:,j), 'PlotTitle', 'PDFF');
            disp('Draw left fat')
            roi_fatL = drawfreehand();
            disp('Draw left muscle w/ high fat')
            roi_mixedL = drawfreehand();
            disp('Draw left muscle w/ low fat')
            roi_muscleL = drawfreehand();
            disp('Draw right fat')
            roi_fatR = drawfreehand();
            disp('Draw right muscle w/ high fat')
            roi_mixedR = drawfreehand();
            disp('Draw right muscle w/ low fat')
            roi_muscleR = drawfreehand();

            % Create masks
            m_fatL(:,:,j) = createMask(roi_fatL);
            m_fatR(:,:,j) = createMask(roi_fatR);
            m_mixedL(:,:,j) = createMask(roi_mixedL);
            m_mixedR(:,:,j) = createMask(roi_mixedR);
            m_muscleL(:,:,j) = createMask(roi_muscleL);
            m_muscleR(:,:,j) = createMask(roi_muscleR);

            close; close
        end
    
        % Save masks
        roi = struct('dims', struct('x', dims.x(i,:), ...
                                    'y', dims.y(i,:), ...
                                    'z', dims.z(i,:)), ...
                      'roi', struct('Name', {"Muscle_HighFF_L", "Muscle_HighFF_R", "Muscle_LowFF_L", "Muscle_LowFF_R", "Fat_L", "Fat_R"}, ...
                                    'Mask', {m_mixedL, m_mixedR, m_muscleL, m_muscleR, m_fatL, m_fatR}));
        save(char(fullfile(pth_all{i}, ['ROI_' snames_all{i}.SusceptibilityMap(4:end) '.mat'])), "roi");
        ROI{i} = roi;
    end
    
    % Plot masks
    mtotal = ROI{i}.roi(1).Mask; 
    for j = 2:length(ROI{i}.roi)
        mtotal = mtotal + ROI{i}.roi(j).Mask; 
    end
    plotmygraph(mag_all{i}.*mtotal, 'PlotTitle', 'Magnitude');
    plotmygraph(PDFF_all{i}.*mtotal, 'PlotTitle', 'PDFF');
    plotmygraph(QSM_all{i}.*mtotal, 'PlotTitle', 'QSM'); 
    
    mbounds = imerode(boundarymask(sum(mtotal,3),8), strel('square',2));
    figure; imagesc(labeloverlay(sum(PDFF_all{i},3), mbounds,'Transparency',0)); axis square; axis off;
end

%% Analysis
for I = 1:length(date)
    % Select groups
    G{I} = ROI{I}.roi;
    
    for j = 1:length(G{I})
        % Isolate ROIs
        G{I}(j).FF.AllData = PDFF_all{I}.*G{I}(j).Mask; 
        G{I}(j).QSM.AllData = QSM_all{I}.*G{I}(j).Mask; 
    
        % Preallocate data
        G{I}(j).idx = 0;
    end
    
    % Iterate through points
    for i = 1:size(PDFF_all{I},1)
        for j = 1:size(PDFF_all{I},2)
            for k = 1:size(PDFF_all{I},3)-1
                % Check if included in which group
                if G{I}(1).Mask(i,j,k)
                    group = 1; gFLAG = true;
                elseif G{I}(2).Mask(i,j,k)
                    group = 2; gFLAG = true;
                elseif G{I}(3).Mask(i,j,k)
                    group = 3; gFLAG = true;
                elseif G{I}(4).Mask(i,j,k)
                    group = 4; gFLAG = true;
                elseif G{I}(5).Mask(i,j,k)
                    group = 5; gFLAG = true;
                elseif G{I}(6).Mask(i,j,k)
                    group = 6; gFLAG = true;
                else
                    gFLAG = false;
                end
    
                if gFLAG
                    % Shift iterator
                    G{I}(group).idx = G{I}(group).idx + 1;
    
                    % Add data to array
                    G{I}(group).FF.Data(G{I}(group).idx) = G{I}(group).FF.AllData(i,j,k); %#ok<*SAGROW>
                    G{I}(group).QSM.Data(G{I}(group).idx) = G{I}(group).QSM.AllData(i,j,k);
                end
            end
        end
    end
    
    % Calculate mean & standard deviation
    for i = 1:length(G{I})
        [G{I}(i).FF.STD, G{I}(i).FF.Mean] = std(G{I}(i).FF.Data);
        [G{I}(i).QSM.STD, G{I}(i).QSM.Mean] = std(G{I}(i).QSM.Data);
    end
end

%% Combine dog data
% Select test cases
fatFLAG = false;
mhfFLAG = true;
mlfFLAG = true;

clear g
for i = 1:length(G) % Iterate through groups
    if i == 1
        i1 = 1;
    else
        i1 = length(g)+1;
    end
    i2 = i1-1 + length(G{i});
    g(i1:i2) = G{i};

    % Correct names
    for j = i1:i2
        name = char(g(j).Name);
        switch comptype
            case 'longitudinal'
                g(j).Name = string([name num2str(date(i))]);
            case 'crosssection'
                g(j).Name = string([name num2str(dog(i))]);
        end
    end
end

% Reorganize
g_old = g;
clear g
i_old = 0;
i_lr = 0;
i_new = 0;
Catch = length(date);
for i = 1:length(g_old)
    i_old = i_old + 1;
    if i_old > Catch
        i_lr = i_lr + 1;
        g_tmp = g_old(i_old-length(date)+length(g_old)/length(date));
        if i_lr == 2 % for left & right
            Catch = Catch + length(date);
            i_old = i_old - length(date);
            i_lr = 0;
        end
    else
        g_tmp = g_old(i_old);
    end

    % Determine roi type
    n_tmp = char(g_tmp.Name);
    switch n_tmp(1:3)
        case 'Fat'
            if fatFLAG
                i_new = i_new + 1;
                g(i_new) = g_tmp;
            end
        case 'Mus'
            switch n_tmp(8:10)
                case 'Low'
                    if mlfFLAG
                        i_new = i_new + 1;
                        g(i_new) = g_tmp;
                    end
                case 'Hig'
                    if mhfFLAG
                        i_new = i_new + 1;
                        g(i_new) = g_tmp;
                    end
            end
    end
end

%% Do you want to test left leg, right leg, or both?

leg = 'left';
dataFF = cell(1,1);
dataQSM = cell(1,1);
labelsFF = cell(1,1);
labelsQSM = cell(1,1);
switch leg
    case 'left'
        iter = 1:2:length(g);
    case 'right'
        iter = 2:2:length(g);
    case 'both'
        iter = 1:1:length(g);
end
% If only muscle regions, split QSM data by muscle
muscFLAG = false;
if ~fatFLAG && mhfFLAG && mlfFLAG
    muscFLAG = true;
    dataFF = cell(1,2);
    labelsFF = cell(1,2);
    dataQSM = cell(1,2);
    labelsQSM = cell(1,2);
end
for J = length(iter):-1:1
    j = iter(J);
    if muscFLAG
        if J > length(iter)/2 % Swapper
            K = 1;
        else
            K = 2;
        end
        dataQSM{K} = [dataQSM{K}; g(j).QSM.Data.'];
        labelsQSM{K} = [labelsQSM{K}; repmat(g(j).Name, [length(g(j).QSM.Data) 1])];
        dataFF{K} = [dataFF{K}; g(j).FF.Data.'];
        labelsFF{K} = [labelsFF{K}; repmat(g(j).Name, [length(g(j).FF.Data) 1])];
    else
        dataQSM{1} = [dataQSM{1}; g(j).QSM.Data.'];
        labelsQSM{1} = labelsFF;
        dataFF{1} = [dataFF{1}; g(j).FF.Data.'];
        labelsFF{1} = [labelsFF; repmat(g(j).Name, [length(g(j).FF.Data) 1])];
    end
end

% Flip data (if necessary)
name1 = char(labelsFF{1}(1));
name2 = char(labelsFF{1}(end));
if str2double(name1(end)) > str2double(name2(end))
    for i = 1:length(dataQSM)
        dataQSM{i} = [dataQSM{i}(labelsQSM{i}==labelsQSM{i}(end)); dataQSM{i}(labelsQSM{i}==labelsQSM{i}(1))];
        labelsQSM{i} = [labelsQSM{i}(labelsQSM{i}==labelsQSM{i}(end)); labelsQSM{i}(labelsQSM{i}==labelsQSM{i}(1))];
        dataFF{i} = [dataFF{i}(labelsFF{i}==labelsFF{i}(end)); dataFF{i}(labelsFF{i}==labelsFF{i}(1))];
        labelsFF{i} = [labelsFF{i}(labelsFF{i}==labelsFF{i}(end)); labelsFF{i}(labelsFF{i}==labelsFF{i}(1))];
    end
end

% Convert to categorical
for i = 1:length(dataQSM)
    labelsQSM{i} = categorical(labelsQSM{i});
    labelsFF{i} = categorical(labelsFF{i});
end

%% Determine limits
limFF = [Inf,-Inf];
limQSM = limFF;
for i = 1:length(dataFF)
    if min(dataFF{i}) < limFF(1), limFF(1) = min(dataFF{i}); end
    if max(dataFF{i}) > limFF(2), limFF(2) = max(dataFF{i}); end
    if min(dataQSM{i}) < limQSM(1), limQSM(1) = min(dataQSM{i}); end
    if max(dataQSM{i}) > limQSM(2), limQSM(2) = max(dataQSM{i}); end
end
limFF = limFF + (diff(limFF)/20)*[-1,1];
limQSM = limQSM + (diff(limQSM)/20)*[-1,1];
% Show data
figure; boxchart(labelsFF{1}, dataFF{1}, MarkerStyle='.', MarkerSize=15, LineWidth=2, JitterOutliers='on'); set(gca, FontSize=20); ylabel('Fat Fraction', FontSize=24, FontWeight='bold'); ylim(limFF);
% title('PDFF', FontSize=28);
% if muscFLAG; title('PDFF1', FontSize=28); end
title('M_{lowF}', FontSize=28); xticklabels(["Initial","2 Month"]);

figure; boxchart(labelsQSM{1}, dataQSM{1}, MarkerStyle='.', MarkerSize=15, LineWidth=2, JitterOutliers='on'); set(gca, FontSize=20); ylabel('Susceptibility (ppm)', FontSize=24, FontWeight='bold'); ylim(limQSM);
% title('QSM', FontSize=28);
% if muscFLAG; title('QSM1', FontSize=28); end
title('M_{lowF}', FontSize=28); xticklabels(["Initial","2 Month"]);

if muscFLAG
    figure; boxchart(labelsFF{2}, dataFF{2}, MarkerStyle='.', MarkerSize=15, LineWidth=2, JitterOutliers='on'); set(gca, FontSize=20); ylabel('Fat Fraction', FontSize=24, FontWeight='bold'); ylim(limFF);
    % title('PDFF2', FontSize=28);
    title('M_{highF}', FontSize=28); xticklabels(["Initial","2 Month"]);

    figure; boxchart(labelsQSM{2}, dataQSM{2}, MarkerStyle='.', MarkerSize=15, LineWidth=2, JitterOutliers='on'); set(gca, FontSize=20); ylabel('Susceptibility (ppm)', FontSize=24, FontWeight='bold'); ylim(limQSM);
    % title('QSM2', FontSize=28);
    title('M_{highF}', FontSize=28); xticklabels(["Initial","2 Month"]);
end

%% Statistical analysis
[g, pFF, statsFF] = RunStats(iter, 'PDFF', g, dataFF, labelsFF, muscFLAG);
[g, pQSM, statsQSM] = RunStats(iter, 'QSM', g, dataQSM, labelsQSM, muscFLAG);

%% Combine into table
name = strings(length(g),1);
n = zeros(length(g),1);
mFF = zeros(length(g),1);
mQSM = zeros(length(g),1);
stdFF = zeros(length(g),1);
stdQSM = zeros(length(g),1);

for i = 1:length(g)
    name(i) = g(i).Name;
    n(i) = g(i).idx;
    mFF(i) = g(i).FF.Mean;
    mQSM(i) = g(i).QSM.Mean;
    stdFF(i) = g(i).FF.STD;
    stdQSM(i) = g(i).QSM.STD;
end

tFF = table(name, n, mFF, stdFF, 'VariableNames',{'Name', 'Sample Size', 'Mean', 'Standard Deviation'});
tQSM = table(name, n, mQSM, stdQSM, 'VariableNames',{'Name', 'Sample Size', 'Mean', 'Standard Deviation'});
%%
pout = sampsizepwr('p',0.5,[],1600);
%%
function [g, p, stats] = RunStats(iter, dtype, g, data, labels, muscFLAG)
    % First, ensure all groups are normally distributed
    normalityFLAG = true;
    for i = iter
        switch dtype
            case 'PDFF'
                g(i).FF.KS = kstest(g(i).FF.Data);
                % Flag if any are not normally distributed
                if g(i).FF.KS && normalityFLAG
                    normalityFLAG = false;
                    disp('Data is not normally distributed!')
                end
            case 'QSM'
                g(i).QSM.KS = kstest(g(i).QSM.Data);
                % Flag if any are not normally distributed
                if g(i).QSM.KS && normalityFLAG
                    normalityFLAG = false;
                    disp('Data is not normally distributed!')
                end
        end
    end
    
    % Only use Bartlett if data is normally distributed: https://www.itl.nist.gov/div898/handbook/eda/section3/eda357.htm
    if normalityFLAG
        vartest = 'Bartlett';
    else
        vartest = 'LeveneAbsolute';
    end
    
    % Final stats
    if muscFLAG
        for i = 1:length(data)
            % Test to see if variance is different between any datasets
            VARp = vartestn(data{i}, labels{i}, 'TestType', vartest, 'Display', 'off');
    
            % Flag if data does not have equal variance
            if VARp < 0.05
                heteroscedasticFLAG = true;
                disp('Data does not have equal variance!')
            else
                heteroscedasticFLAG = false;
                disp('Data does have equal variance!')
            end
    
            if normalityFLAG
                % Two-sample t-test
                if heteroscedasticFLAG
                    [~, p{i}, ~, stats{i}] = ttest2(data{i}(labels{i}==labels{i}(1)), data{i}(labels{i}==labels{i}(end)), 'Vartype', 'unequal'); %#ok<*AGROW>
                else
                    [~, p{i}, ~, stats{i}] = ttest2(data{i}(labels{i}==labels{i}(1)), data{i}(labels{i}==labels{i}(end)), 'Vartype', 'equal');
                end
            
            else
                % Mann-Whitney U Test
                [p{i}, ~, stats{i}] = ranksum(data{i}(labels{i}==labels{i}(1)), data{i}(labels{i}==labels{i}(end)));
    
                % Brunner-Munzel test
                p{i}(2) = brunner_munzel(data{i}(labels{i}==labels{i}(1)), data{i}(labels{i}==labels{i}(end)));

                if (p{i}(1) < 0.05 && p{i}(2) >= 0.05) || (p{i}(1) >= 0.05 && p{i}(2) < 0.05)
                    error('diff values');
                end
            end
        end
    end
end