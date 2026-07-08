% Fat water plotting
function Z = showSwaps(raw, ff, z, ss) 
% Inputs:
%      raw: raw magnitude data
%       ff: magnitude fat fraction map
%        z: loop iteration
%       ss: step size
%
% Outputs:
%        Z: end index

    % Ensure not too close to end
    if z+ss-1 > size(raw,3)
        Z = size(raw,3);
    else
        Z = z+ss-1;
    end

    % Plot data
    figure('Name', ['Slice' num2str(z) 'to' num2str(Z)], 'NumberTitle', 'off','WindowState','maximized'); 
    tiledlayout(2,Z-z+1, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile; imagesc(raw(:,:,z)); set(gca,'XTick',[]); set(gca,'YTick',[]); axis square; colormap gray; ylabel('Raw'); title(num2str(z));
    for zz = z+1:Z
        nexttile; imagesc(raw(:,:,zz)); axis off; axis square; colormap gray; title(num2str(zz));
    end
    nexttile; imagesc(ff(:,:,z)); set(gca,'XTick',[]); set(gca,'YTick',[]); axis square; colormap gray; ylabel('FF');
    for zz = z+1:Z
        nexttile; imagesc(ff(:,:,zz)); axis off; axis square; colormap gray;
    end
    pause(0.5)
end