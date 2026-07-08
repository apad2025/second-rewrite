% save the 25th slice of every 358x358x50 struct
% so it becomes 358x358

Data = D.Data;
fnames = fieldnames(Data);

Data.(fnames{1}) = Data.(fnames{1})(:,:,25);

for i = 3:14
    Data.(fnames{i}) = Data.(fnames{i})(:,:,25);
end

Data.species(1).amps = Data.species(1).amps(:,:,25);
Data.species(2).amps = Data.species(2).amps(:,:,25);