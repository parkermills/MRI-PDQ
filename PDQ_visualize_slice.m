%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @name PDQ_visualize_slice.m
% @author Parker Mills, Ahrens Lab, Carnegie Mellon 2009
% @brief Visualizes a specified slice from 3D phase volume, marking
%        the location and susceptibility of each dipole present in
%        the slice.
%
% ==================== INPUT PARAMETERS ======================
% @param   dataset         (MRIdata)    Volumetric phase dataset
% @param   slice           (1D Float)    Slice to be visualized
% @param   slice_thickness (1D Float)    Maximum distance, in slices, that a dipole can be from
%                                        currently visualized slice in order for it to be displayed.
%                                        Default is half template size (i.e., half dipole size)
%
% @param   xcorr_cutoff    (1D Float)    Cross-correlation quality threshold value; [0.0 1.0]
% @param   full_screen     (1D Float)    1: Full screen view of stuff, 0: Normal window
%
% ==================== RETURNED DATA =========================
% @return  fig        (Float)       Figure handle
%
% ====================== ASSUMPTIONS =========================
% @assume
%
% ==================== DISPLAYED PRODUCTS ====================
%
% ==================== EXAMPLE USAGE =========================
% fig1 = PDQ_visualize_slice(MRIdata, slice, slice_thickness, xcorr_cutoff);
% fig1 = PDQ_visualize_slice(MRIdata, 16, 2, 0.3);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function fig1 = PDQ_visualize_slice(MRIdata, slice, slice_thickness, xcorr_cutoff, mdm_range, full_screen)




%% Determine what phase image is being used for visualization

if(exist('MRIdata.rauscher','var'))
    backdrop = MRIdata.rauscher;
else
    backdrop = MRIdata.high_pass;
end

[x_dim y_dim z_dim] = size(backdrop);

%% Sanity checks
if(slice < 1 || slice > z_dim)
	error(['visualize_PDQ_slice_slice: Selected slice number, ',num2str(slice),', is outside of volume']);
end

%% Show the desired slice
if (full_screen)
   fig1 = images(backdrop(:,:,slice),[-1 1],'fullscreen');
else
    fig1 = images(backdrop(:,:,slice),[-1 1]);
end

% If slice thickness not provided, set default: half dipole's diameter
if(~exist('slice_thickness','var'))
    slice_thickness = floor(max(size(MRIdata.PDQ.dipoles_hp(1).phase))/2);
else
   if (isempty(slice_thickness))
       slice_thickness = floor(max(size(MRIdata.PDQ.dipoles_hp_hp(1).phase))/2);
   end
end

% Go through all dipoles
for j_1 = 1:length(MRIdata.PDQ.dipoles_hp)
	
	% If dipole within "z-depth" of user-selected slice, label it!
	if(  (MRIdata.PDQ.dipoles_hp(j_1).z >= (slice - slice_thickness))...
			&& (MRIdata.PDQ.dipoles_hp(j_1).z <= (slice + slice_thickness))...
			&& (MRIdata.PDQ.dipoles_hp(j_1).xcorr > xcorr_cutoff)...
            && (MRIdata.PDQ.dipoles_hp(j_1).m_core > mdm_range(1)))

		arrow_x = [(MRIdata.PDQ.dipoles_hp(j_1).y + rand*30-15) MRIdata.PDQ.dipoles_hp(j_1).y];
		arrow_y = [(x_dim+1)-(MRIdata.PDQ.dipoles_hp(j_1).x + rand*30-15) (x_dim+1)-MRIdata.PDQ.dipoles_hp(j_1).x];
		[arrow_x,arrow_y] = dsxy2figxy(gca, arrow_x, arrow_y);

		annotation(fig1,'textarrow', arrow_x, arrow_y,...
			'TextEdgeColor','none',...
			'TextLineWidth',2,...
			'TextColor',[1 0 0],...
			'FontSize',10,...
            'String',{...
            strcat('#'     ,num2str(j_1)                                ),...
            strcat('xcorr=',num2str(MRIdata.PDQ.dipoles_hp(j_1).xcorr ,2)),...
            strcat('mdm='  ,num2str(MRIdata.PDQ.dipoles_hp(j_1).m_core,2))...
            },'LineWidth',1,...
            'Color',[1 0 0]);
	end
end

%%%%%%%%%%%
% EOF
%%%%%%%%%%%