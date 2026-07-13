% save the nth slice of every 358x358x50 struct
% so it becomes 358x358
clear; clc; close all

slice = 12;

% Drop all fields except Data
path = "C:\Users\apad2\Desktop\second_rewrite\results\FWbigc_zip.mat";
Data = load(path, 'D').D.Data;
save(path, 'Data');

% loop over Data and trim down each array
fnames = fieldnames(Data);
for i = [1, 3:6, 8:14]
    Data.(fnames{i}) = Data.(fnames{i})(:,:,slice);
end
for i = 1:2
    Data.species(i).amps = Data.species(i).amps(:,:,slice);
end
Data.corrected_bipolar_signal = Data.corrected_bipolar_signal(:,:,slice,:,:);
save(path, 'Data');
