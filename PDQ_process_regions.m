%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates statistics on regions, if specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MRIdata = PDQ_process_regions(MRIdata, xcorr_cutoff)

%% Preferences
pref_num_random_points = 9000; % Number of random points selected for calculating distance from ROI

%% Initialize

% Extract information from MRIdata file
regions_dims = size(MRIdata.regions(1).mask);
num_regions = length(MRIdata.regions);
num_dipoles = length([MRIdata.PDQ.dipoles_hp]);

% Create coordinate grid
coordinate_grid_x = 1:regions_dims(1);
coordinate_grid_y = 1:regions_dims(2);
coordinate_grid_z = 1:regions_dims(3);



%% For each region, perform analysis
for curr_region = 1:num_regions
    
    % Invert the region's mask
    warning off MATLAB:divideByZero
    inverted_mask = 1./MRIdata.regions(curr_region).mask;
    warning on MATLAB:divideByZero
    
    
  
    %% Generate random grid of points in sample volume for later normalization
    if(~isfield(MRIdata.regions, 'random_grid'))
        % Create random x,y,z coordinate points
        random_points_x = round(  rand(pref_num_random_points, 1) .* (regions_dims(1) - 1)  )  + 1;
        random_points_y = round(  rand(pref_num_random_points, 1) .* (regions_dims(2) - 1)  )  + 1;
        random_points_z = round(  rand(pref_num_random_points, 1) .* (regions_dims(3) - 1)  )  + 1;
        
        % For each random point
        for jj = 1:pref_num_random_points
            
            % If the point is within entire MRIdata's mask, but not in current region
            if(MRIdata.mask(random_points_x(jj),random_points_y(jj),random_points_z(jj)) && ~MRIdata.regions(curr_region).mask(random_points_x(jj),random_points_y(jj),random_points_z(jj)))
                
                % Generate 3D distance grid for this specific point
                x_dist = ((coordinate_grid_x - random_points_x(jj)).^2)'* ones(1, regions_dims(1));
                y_dist = rot90((   (coordinate_grid_y - random_points_y(jj)).^2)'* ones(1, regions_dims(2)) );
                z_dist = (coordinate_grid_z - random_points_z(jj)).^2;
                xy_dist = x_dist + y_dist;
                for j_4 = 1:regions_dims(3) % Z-dim
                    random_points_xyz_dist(:,:,j_4) = sqrt(xy_dist + z_dist(j_4));
                end
                
                % Find this point's minimum distance to the current region
                distance_masked = random_points_xyz_dist .* inverted_mask;
                MRIdata.regions(curr_region).random_grid(jj) = min(min(min(distance_masked)));
                
            end
        end
    end
    
    
    
    %% Go through all dipoles, counting those located within this region
    if(~isfield(MRIdata.regions, 'dipole_min_dist'))
        
        % Set counts to zero
        MRIdata.regions(curr_region).dipole_count = 0;
        dipole_not_in_region = 0;
        
        for curr_dipole = 1:num_dipoles
            
            % Give status report
            if(mod(curr_dipole,100)==0)
                disp([num2str(curr_dipole),' Dipoles analyzed. (',num2str(curr_dipole/num_dipoles * 100),'% done with Region ',num2str(curr_region),') '])
            end
            
            % If current dipole is above xcorr_cutoff
            if (MRIdata.PDQ.dipoles_hp(curr_dipole).xcorr > xcorr_cutoff)
                location = [MRIdata.PDQ.dipoles_hp(curr_dipole).x    MRIdata.PDQ.dipoles_hp(curr_dipole).y    MRIdata.PDQ.dipoles_hp(curr_dipole).z];
                
                % If within region, just add one to the dipole count for this dipole and visualize that slice, since it's interesting
                if(MRIdata.regions(curr_region).mask(location(1),location(2),location(3)))
                    MRIdata.regions(curr_region).dipole_count = MRIdata.regions(curr_region).dipole_count + 1;
                    %%%PDQ_visualize_slice(MRIdata, MRIdata.PDQ.dipoles_hp(curr_dipole).z, 0, xcorr_cutoff, 1);
                    %%%images(MRIdata.mag(:,:,MRIdata.PDQ.dipoles_hp(curr_dipole).z));
                    
                    % If not within region, increment and determin distance of this dipole to the region
                else
                    dipole_not_in_region = dipole_not_in_region + 1;
                    
                    % Generate 3D distance grid
                    x_dist = ((coordinate_grid_x - location(1)).^2)'* ones(1,regions_dims(1));
                    y_dist = rot90((   (coordinate_grid_y - location(2)).^2)'* ones(1,regions_dims(2)) );
                    z_dist = (coordinate_grid_z - location(3)).^2;
                    xy_dist = x_dist+y_dist;
                    for j_4 = 1:regions_dims(3)
                        xyz_dist(:,:,j_4) = sqrt(xy_dist + z_dist(j_4));
                    end
                    
                    % Find nearest distance to other regions
                    distance_masked = xyz_dist .* inverted_mask;
                    MRIdata.regions(curr_region).dipole_min_dist(dipole_not_in_region) = min(min(min(distance_masked)));
                    
                end
            end
        end
    end
    
    % Prepare vectors by removing zeros
    MRIdata.regions(curr_region).random_grid = vectorize_and_remove_zeros(MRIdata.regions(curr_region).random_grid);
    
    % Visualize
    random_grid_hist = hist([zeros(1,999) MRIdata.regions(1).random_grid max(regions_dims).*ones(1,999)], 30);
    distribution_grid_hist = hist([zeros(1,1) MRIdata.regions(1).dipole_min_dist max(regions_dims).*ones(1,1)], 30);
    figure; bar(distribution_grid_hist ./ random_grid_hist);
    
    random_grid_hist = hist([zeros(1,999) MRIdata.regions(1).random_grid max(regions_dims).*ones(1,999)], 20);
    distribution_grid_hist = hist([zeros(1,1) MRIdata.regions(1).dipole_min_dist max(regions_dims).*ones(1,1)], 20);
    figure; bar(distribution_grid_hist ./ random_grid_hist);
    
    random_grid_hist = hist([zeros(1,999) MRIdata.regions(1).random_grid max(regions_dims).*ones(1,999)], 10);
    distribution_grid_hist = hist([zeros(1,1) MRIdata.regions(1).dipole_min_dist max(regions_dims).*ones(1,1)], 10);
    figure; bar(distribution_grid_hist ./ random_grid_hist);
    
end % End for each region



%%%%%%%%%%%%%%%%
%     EOF
%%%%%%%%%%%%%%%%