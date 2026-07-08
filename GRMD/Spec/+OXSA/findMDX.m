function mDx = findMDX(pk, np)
mDx = cellfun(@(x) [0,0], cell(1, np), 'UniformOutput', false); % multiplet indices (along results structures/arrays)

for p = 1:np
    if p > 1
        mDx{p}(1) = mDx{p-1}(end) + 1;
    else
        mDx{p}(1) = 1;
    end

    % Check if multiplet
    if iscell(pk.bounds(p).peakName)
        mDx{p} = mDx{p}(1):(mDx{p}(1)+length(pk.bounds(p).peakName)-1);
    else
        mDx{p} = mDx{p}(1);
    end
end
end