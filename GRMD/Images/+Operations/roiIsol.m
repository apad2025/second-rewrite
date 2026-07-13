% Trim FOV
function D = roiIsol(D)
% Inputs:
%     data: data structure
% 
% Outputs:
%      data: data structure

    % Reset plotrange
    newsize = [D.Size(1), 1; D.Size(2), 1; 1, D.Size(3)];
    
    % Sum mask slices
    m_temp = sum(D.Data.Mask, [1, 2]);

    % Find mask limits
    [~, newsize(3,1)] = find(m_temp, 1, 'first');
    [~, newsize(3,2)] = find(m_temp, 1, 'last');

    % Iterate through slices
    for k = newsize(3,1):newsize(3,2)
        m_temp = D.Data.Mask(:,:,k);
            
        % Find mask limits
        [~, j_f] = find(m_temp, 1, 'first');
        [~, j_l] = find(m_temp, 1, 'last');
        [~, i_f] = find(m_temp.', 1, 'first');
        [~, i_l] = find(m_temp.', 1, 'last');

        % Update mask limits as necessary
        if j_f < newsize(2,1)
            newsize(2,1) = j_f;
        end
        if j_l > newsize(2,2)
            newsize(2,2) = j_l; 
        end
        if i_f < newsize(1,1)
            newsize(1,1) = i_f; 
        end
        if i_l > newsize(1,2)
            newsize(1,2) = i_l; 
        end
    end
    
    % Find difference in sides
    roidiff = abs((newsize(2,2) - newsize(2,1)) - (newsize(1,2) - newsize(1,1)));
    
    % Change to square ROI
    if roidiff ~= 0
        if (newsize(2,2) - newsize(2,1)) > (newsize(1,2) - newsize(1,1))
            % Add side difference to smaller side
            newsize(1,1) = newsize(1,1) - floor(roidiff/2);
            newsize(1,2) = newsize(1,2) + floor(roidiff/2);
    
            % pad the larger side to even number to match smaller side
            if rem(roidiff, 2) > 0
                newsize(1,2) = newsize(1,2) + 1;
            end
        else
            % Add side difference to smaller side
            newsize(2,1) = newsize(2,1) - floor(roidiff/2);
            newsize(2,2) = newsize(2,2) + floor(roidiff/2);
    
            % pad the larger side to even number to match smaller side
            if rem(roidiff, 2) > 0
                newsize(2,2) = newsize(2,2) + 1;
            end
        end
    end
    
    % Check if not an even number
    if rem(newsize(1,2) - newsize(1,1) + 1, 2) > 0
        % Subtract one index from each dimension
        newsize(1:2,2) = newsize(1:2,2) - 1;
    end

    % Determine padding amount
    % padsize = (2^nextpow2(newsize(1,2) - newsize(1,1)) - (newsize(1,2) - newsize(1,1)))/2;
    if D.Flags.Interpolated
        padsize = 20;
    else
        padsize = 10;
    end

    % Ensure pad size is not larger than trimmed pixels
    if padsize > newsize(1,1)
        padsize = newsize(1,1) - 1;
    end
    if padsize > newsize(2,1)
        padsize = newsize(2,1) - 1;
    end
    if padsize > size(m_temp,2) - newsize(1,2)
        padsize = size(m_temp,2) - newsize(1,2) - 1;
    end
    if padsize > size(m_temp,2) - newsize(2,2)
        padsize = size(m_temp,2) - newsize(2,2) - 1;
    end

    % Pad ROI so tissue is not on edge of figures
    newsize(1,1) = newsize(1,1) - padsize;
    newsize(1,2) = newsize(1,2) + padsize;
    newsize(2,1) = newsize(2,1) - padsize;
    newsize(2,2) = newsize(2,2) + padsize;

    % Trim data
    D.Data.Image = D.Data.Image(newsize(1,1):newsize(1,2), newsize(2,1):newsize(2,2), newsize(3,1):newsize(3,2), :);
    D.Data.Mask = D.Data.Mask(newsize(1,1):newsize(1,2), newsize(2,1):newsize(2,2), newsize(3,1):newsize(3,2));
    D.Data.WeightedMagnitude = D.Data.WeightedMagnitude(newsize(1,1):newsize(1,2), newsize(2,1):newsize(2,2), newsize(3,1):newsize(3,2));
    D.Size = size(D.Data.Image);

    % Update flags
    D.Flags.Trimmed.InPlane = true;
    D.TrimmedIndices = newsize;
end