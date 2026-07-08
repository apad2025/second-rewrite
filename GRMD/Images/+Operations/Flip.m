% Flip data to account for negative B0 direction
function [data, flippedFLAG] = Flip(data)
flippedFLAG = false(length(data.B0Direction),1);
for i = 1:length(data.B0Direction)
    if data.B0Direction(i) < 0
        % Flip data along that axis
        fprintf('\nTemporarily flipping data along dimension %i to correct for negative B0 direction\n', i)
        datafields = fieldnames(data.Data);
        for j = 1:numel(datafields)
            data.Data.(datafields{j}) = flip(data.Data.(datafields{j}),i);
        end
        data.B0Direction(i) = data.B0Direction(i)*-1;
        flippedFLAG(i) = true;
    else
        flippedFLAG(i) = false;
    end
end
end