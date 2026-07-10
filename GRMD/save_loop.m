% save the nth slice of every 358x358x50 struct
% so it becomes 358x358
clear; clc; close all

slice = 12;

% Drop all fields except Data
myslice = "C:\Users\apad2\OneDrive\Desktop\mat_backup\FWbigc_zip.mat";
S = load(myslice, 'D');
Data = S.D.Data;
save(myslice, 'Data');

% loop over Data and trim down each array
clear myslice;
myslice = load("C:\Users\apad2\OneDrive\Desktop\mat_backup\FWbigc_zip.mat");
path = "C:\Users\apad2\OneDrive\Desktop\mat_backup\FWbigc_zip.mat";
Data = myslice.Data;
fnames = fieldnames(Data);

Data.(fnames{1}) = Data.(fnames{1})(:,:,slice);
for i = 3:6
    Data.(fnames{i}) = Data.(fnames{i})(:,:,slice);
end

for i = 8:14
    Data.(fnames{i}) = Data.(fnames{i})(:,:,slice);
end

% trim species
Data.species(1).amps = Data.species(1).amps(:,:,slice);
Data.species(2).amps = Data.species(2).amps(:,:,slice);

% trim corrected_bipolar
Data.corrected_bipolar_signal = Data.corrected_bipolar_signal(:,:,slice,:,:);
save(path, "Data");