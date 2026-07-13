% Save data
function D = Save(D, FLAG, path, sname, verbose)
if FLAG
    % Check data to ensure correct size
    % fn = fieldnames(D.Data);
    % for i = 1:length(fn)
    %     s = size(D.Data.(fn{i}));
    %     if ~all(s == D.Size(1:length(s)))
    %         error('Data size does not match')
    %     end
    % end
    if verbose; fprintf('\nSaving data\n'); end
    if isfield(D.Data, 'Mask')
        D.Data.Mask = D.Data.Mask > 0;
    end
    save(char(fullfile(path, [sname, '.mat'])), "D");
end
end