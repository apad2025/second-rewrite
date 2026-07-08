% Correct for flipped data
function data = unFlip(data, flippedFLAG)
for i = 1:length(data.B0Direction)
    if flippedFLAG(i)
        % Flip data along that axis
        datafields = fieldnames(data.Data);
        for j = 1:numel(datafields)
            data.Data.(datafields{j}) = flip(data.Data.(datafields{j}),i);
        end
        data.B0Direction(i) = data.B0Direction(i)*-1;
    end
end
end