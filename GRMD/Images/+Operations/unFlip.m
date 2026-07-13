% Correct for flipped data
function D = unFlip(D, flippedFLAG)
for i = 1:length(D.B0Direction)
    if flippedFLAG(i)
        % Flip data along that axis
        datafields = fieldnames(D.Data);
        for j = 1:numel(datafields)
            D.Data.(datafields{j}) = flip(D.Data.(datafields{j}),i);
        end
        D.B0Direction(i) = D.B0Direction(i)*-1;
    end
end
end