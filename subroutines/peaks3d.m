%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @name peaks3d.m
% @author Parker Mills, Ahrens Lab, Carnegie Mellon 2009
% @brief Finds all peaks in a 3D image matrix, by finding local maxima
%        first in a series of one dimensional lines (2D), then rotates
%        a volume to look for all local maxima along the third dimension (3D)
%
% @depend Requires localmaxmin, which highlights local maxima
%
% ==================== INPUT PARAMETERS ======================
% @param    image        (3D Float) Data to be searched for peaks
% @param    threshold    (Float)    Only peaks with values above this value are valid
%
% ==================== RETURNED DATA =========================
% @return   peaks    (3D Logical)  Binary image marked with locations of each peak
%
%
% ==================== ASSUMPTIONS ======================
% @assume   Peaks have a value > 0.0
% 
%
% ==================== DISPLAYED PRODUCTS ====================
% @product
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function peaks = peaks3d(image,threshold)

%% Initialization

% Sanity checks
if(~isreal(image)) %Make sure image is real-valued
    error('peaks3d: Image is not real-valued. Quitting.');
end

if(size(threshold) ~= [1 1]) %Make threshold is a single value
    error('peaks3d: threshold not a single value. Quitting.');
end

% Store sizes
[x_dim y_dim z_dim] = size(image);
pass1 = single(zeros(x_dim,y_dim,z_dim)); % 2D pass
peaks = single(zeros(x_dim,y_dim,z_dim)); % Final 3D peaks result

%% Check for peaks
% First in 2D
for j=1:z_dim
    slice_x = rot90(rot90(rot90(localmaxmin(rot90(image(:,:,j).*(image(:,:,j) > threshold)),'max'))));
    slice_y = localmaxmin(image(:,:,j).*(image(:,:,j) > threshold),'max');
    pass1(:,:,j) = slice_x .* slice_y;
end

% Then in 3D
for j=1:y_dim
    slice_z = rot90(rot90(rot90(localmaxmin(     rot90(reshape(image(:,j,:),x_dim,z_dim) .* reshape((image(:,j,:) > threshold),x_dim,z_dim))     ,'max'))));
    peaks(:,j,:)=slice_z .* reshape(pass1(:,j,:),x_dim,z_dim);
end



%%%%%%%%%%%%
%    EOF
%%%%%%%%%%%%