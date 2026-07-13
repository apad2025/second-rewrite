% Flip data to account for negative B0 direction
function [D, flippedFLAG] = Flip(D)
flippedFLAG = false(length(D.B0Direction),1);
for i = 1:length(D.B0Direction)
    if D.B0Direction(i) < 0
        % Flip data along that axis
        fprintf('\nTemporarily flipping data along dimension %i to correct for negative B0 direction\n', i)
        datafields = fieldnames(D.Data);
        for j = 1:numel(datafields)
            D.Data.(datafields{j}) = flip(D.Data.(datafields{j}),i);
        end
        D.B0Direction(i) = D.B0Direction(i)*-1;
        flippedFLAG(i) = true;
    else
        flippedFLAG(i) = false;
    end
end
end