%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @file images.m
% @author Parker Mills, Ahrens Lab, Carnegie Mellon
% @brief Utility for quickly visualizing 2D or 3D data in 256 grayscale
%
% ==================== INPUT PARAMETERS ======================
% @param  data       (2D/3D Float)    Image to be visualized
% @param  bounds     (1D Float)       Lower- and upper-bound of colormap window
% @param  fig_size   (String)         'normal': Normal figure size
%                                     'fullscreen': Full-screen figure
%                                     
% ==================== RETURNED DATA =========================
% @return fig        (Float)      Figure handle
%
%
% ====================== ASSUMPTIONS =========================
% @assume If a 3D matrix is provided, only the central slice is shown
% @assume If data is complex, a warning is provided and the magnitude image is calculated and displayed
%
%
% ==================== SAMPLE USAGE ==========================
% fig = images(data, bounds, fig_size)
% fig = images(phase_image, [-2 2], 'fullscreen')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function fig = images(data, bounds, fig_size)

%% Check inputs

% Ensure full_screen is correct
if(exist('fig_size','var'))
	if(~strcmp(fig_size,'fullscreen') && ~strcmp(fig_size,'normal'))
		error('images: fig_size must be specified as either "fullscreen" or "normal"');
	end
end

% If is 4D data, alert user and select first volume only
if(ndims(data) == 4)
	warning('images: Supplied image is 4-D. Only showing first 3D volume.');
	data = data(:,:,:,1);
end



%% If this is a 3D data, show center slice
if(ndims(data) == 3)
	[x_dim y_dim z_dim] = size(data);
	slice_num = ceil(z_dim/2.0);
	data = data(:,:,slice_num);
end



%% Create new figure, apply bounds
if(exist('fig_size','var'))
	if(strcmp(fig_size,'fullscreen'))
		scrsz = get(0,'ScreenSize');
		fig = figure('Position',[1 1 scrsz(3)/1.05 scrsz(4)/1.05]); %[Left bottom width height]
	end
else
	fig = figure;
end



%% If data is complex, calculate its magnitude
if (~isreal(data))
	warning('images: Provided image is complex, displaying only its magnitude.');
	data = abs(data);
end



%% View image with or without image intensity map scale
if(exist('bounds','var') && ~isempty(bounds))
    imagesc(data,[bounds(1) bounds(2)]);
else
    imagesc(data);
end



%% Apply colormap, fix axes, stylize frame
colormap(gray(256)); % 256 Grays
axis equal; % Correct aspect ratio
set(gca, 'color', 'black'); %Background black
set(gcf, 'color', 'white'); %Frame white



%%%%%%%%%%%
%   EOF
%%%%%%%%%%%