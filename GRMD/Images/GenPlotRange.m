% Generate plot range
function plrng = GenPlotRange(datasize)
    % Determine z step size
    zstep = floor(datasize(3)/6);

    plrng = {1:datasize(1), 1:datasize(2), 1:zstep:datasize(3), 1};
end