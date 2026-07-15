%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright Véronique Fortier 2017 - MIT License
%
% Quality guided unwrapping                                               
%
% Inputs:                                                                 
% -Phase image (im_phase)                                                 
% -Object mask (im_mask): mask of the object used to reduce the unwrapping time.
% **The mask should not reach the edges of the field of view in the in plane 
% dimension.                                               
% -Quality threshold (qualityCutoff): This threshold has to be in the form
% of 100/(percent threshold). Default is 3.5 (equivalent to 29%)          
%
% Outputs:                                                                
% -Unwrapped result (im_unwrapped)                                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ im_unwrapped ] = qualityGuidedUnwrapping(im_phase, im_mask, im_mag, qualityCutoff )

matrixSize=size(im_mask);

im_unwrapped=zeros(matrixSize);               % starting matrix for unwrapped phase
adjoin=zeros(matrixSize);                     % starting matrix for adjoin matrix
unwrapped_binary=zeros(matrixSize);           % Binary image to mark unwrapped pixels

if nargin<4
    qualityCutoff=3.5;
end

%% Determine if 2D or 3D
if numel(matrixSize) >= 3
    volume = true;
else
    volume = false;
end

%% Calculate second difference quality map

im_phase_hor_right=zeros(matrixSize);
im_phase_hor_left=zeros(matrixSize);
im_phase_ver_top=zeros(matrixSize);
im_phase_ver_bot=zeros(matrixSize);
im_phase_norm_back=zeros(matrixSize);
im_phase_norm_front=zeros(matrixSize);

if volume
    im_phase_hor_right(2:(end-1),2:(end-1),2:(end-1))=im_phase((2:(end-1))-1,2:(end-1),2:(end-1));
    im_phase_hor_left(2:(end-1),2:(end-1),2:(end-1))=im_phase((2:(end-1))+1,2:(end-1),2:(end-1));
    
    im_phase_ver_top(2:(end-1),2:(end-1),2:(end-1))=im_phase(2:(end-1),(2:(end-1))+1,2:(end-1));
    im_phase_ver_bot(2:(end-1),2:(end-1),2:(end-1))=im_phase(2:(end-1),(2:(end-1))-1,2:(end-1));
    
    im_phase_norm_back(2:(end-1),2:(end-1),2:(end-1))=im_phase(2:(end-1),2:(end-1),(2:(end-1))-1);
    im_phase_norm_front(2:(end-1),2:(end-1),2:(end-1))=im_phase(2:(end-1),2:(end-1),(2:(end-1))+1);
else
    im_phase_hor_right(2:(end-1),2:(end-1))=im_phase((2:(end-1))-1,2:(end-1));
    im_phase_hor_left(2:(end-1),2:(end-1))=im_phase((2:(end-1))+1,2:(end-1));

    im_phase_ver_top(2:(end-1),2:(end-1))=im_phase(2:(end-1),(2:(end-1))+1);
    im_phase_ver_bot(2:(end-1),2:(end-1))=im_phase(2:(end-1),(2:(end-1))-1);
end

H=wrapToPi(im_phase_hor_right-im_phase)-wrapToPi(im_phase-im_phase_hor_left);
V=wrapToPi(im_phase_ver_bot-im_phase)-wrapToPi(im_phase-im_phase_ver_top);
if volume
    N=wrapToPi(im_phase_norm_back-im_phase)-wrapToPi(im_phase-im_phase_norm_front);
else
    N=0;
end

SD=(sqrt(H.^2+V.^2+N.^2));


%% Identify the starting seed point on the phase quality map (central region of the central slice)
im_phase_quality=SD.*im_mask;
SE=strel('square',5);

if volume
    if all(im_mask, "all")
        im_phase_qualityMask=imgaussfilt(im_phase_quality(:,:,round(matrixSize(3)/2)),5).*imerode(im_mask(:,:,round(matrixSize(3)/2)),SE);
    else
        im_phase_qualityMask=im_phase_quality(:,:,round(matrixSize(3)/2)).*imerode(im_mask(:,:,round(matrixSize(3)/2)),SE);
    end
else
    if all(im_mask, "all")
        im_phase_qualityMask=imgaussfilt(im_phase_quality,5).*imerode(im_mask,SE);
    else
        im_phase_qualityMask=im_phase_quality.*imerode(im_mask,SE);
    end
end
im_phase_qualityMask=uint16(im_phase_qualityMask*100);
im_phase_qualityMask(find(im_phase_qualityMask==0))=1000;   % Set region outside of im_mask to an arbitrary very low quality 

[a, loc]=min(im_phase_qualityMask(:));  % The minimum is the highest quality point
[xpoint,ypoint]=ind2sub(size(im_phase_qualityMask),loc);

colref=round(ypoint);
rowref=round(xpoint);
if volume
    zpoint=round(matrixSize(3)/2);
end
 
%% Define the seed point as unwrapped
if volume
    im_unwrapped(rowref,colref,zpoint)=im_phase(rowref,colref,zpoint);                        
    unwrapped_binary(rowref,colref,zpoint)=1;
else
    im_unwrapped(rowref,colref)=im_phase(rowref,colref);                        
    unwrapped_binary(rowref,colref)=1;
end

%% Add the adjoining voxels of the seed point to the adjoin matrix
if volume
    if im_mask(rowref-1, colref, zpoint)==1 adjoin(rowref-1, colref, zpoint)=1; end       
    if im_mask(rowref+1, colref, zpoint)==1 adjoin(rowref+1, colref, zpoint)=1; end
    if im_mask(rowref, colref-1, zpoint)==1 adjoin(rowref, colref-1, zpoint)=1; end
    if im_mask(rowref, colref+1, zpoint)==1 adjoin(rowref, colref+1, zpoint)=1; end
    if im_mask(rowref, colref, zpoint-1)==1 adjoin(rowref, colref, zpoint-1)=1; end
    if im_mask(rowref, colref, zpoint+1)==1 adjoin(rowref, colref, zpoint+1)=1; end
else
    if im_mask(rowref-1, colref)==1 adjoin(rowref-1, colref)=1; end       
    if im_mask(rowref+1, colref)==1 adjoin(rowref+1, colref)=1; end
    if im_mask(rowref, colref-1)==1 adjoin(rowref, colref-1)=1; end
    if im_mask(rowref, colref+1)==1 adjoin(rowref, colref+1)=1; end
end

%% Quality-guided unwrapping
if volume
    im_unwrapped=GuidedFloodFill_3D(im_phase, im_unwrapped, unwrapped_binary, im_phase_quality, adjoin, im_mask,qualityCutoff);
else
    derivative_variance=PhaseDerivativeVariance(im_phase, im_mask);
    im_unwrapped=GuidedFloodFill_r1(im_phase, im_mag, im_unwrapped, unwrapped_binary, derivative_variance, adjoin, im_mask);
end

end

