%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @name PDQ_visualize.m
% @author Parker Mills, Ahrens Lab, Carnegie Mellon 2009
% @brief Visualizes the results of an MRIdata that has been processed by PDQ
% @trusted no
% @robust no
% @commented no
% @optimized no
% @parallelized no
%
% ==================== INPUT PARAMETERS ======================
% @param MRIdata                  (1D MRIdata)      MRI dataset(s) to be visualized by PDQ (see README.txt for datatype details)
%
% ==================== RETURNED DATA =========================
% @return
%
% ==================== ASSUMPTIONS ===========================
% @assume If a vector of MRIdata is provided, visualize_PDQ will visualize all datasets TOGETHER\
% @assume PDQ results for high-passed phase are used, not ramp-removed phase results or other phase product results
%
% ==================== EXAMPLE USAGE =========================
% visualize_PDQ(MRIdata)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PDQ_visualize(MRIdata, fixed_xcorr_cutoff, var_xcorr_cutoff, fixed_mdm_range, var_mdm_range, just_count)



%% Preferences
pref_fixed_xcorr_cutoff_default = 0.3;
pref_var_xcorr_cutoff_default = 0.6;
pref_fixed_mdm_range_default = [-inf inf];
pref_var_mdm_range_default = [-inf inf];
pref_display_gaussian_fits = 1;



%% Initializations and sanity-checks

% If xcorr_cutoff values are not specified, set them to default
if(~exist('fixed_xcorr_cutoff','var'))
    warning(['Fixed XCORR cutoff not specified, setting to default of ',num2str(pref_fixed_xcorr_cutoff_default)]);
    fixed_xcorr_cutoff = pref_fixed_xcorr_cutoff_default;
end
if(~exist('var_xcorr_cutoff','var'))
    warning(['Variable XCORR cutoff not specified, setting to default of ',num2str(pref_fixed_xcorr_cutoff_default)]);
    var_xcorr_cutoff = pref_var_xcorr_cutoff_default;
end
if(~exist('fixed_mdm_range','var'))
    warning(['Range of acceptable fixed-radius mdm values not specified. Setting default to ',num2str(pref_fixed_mdm_range_default)]);
    fixed_mdm_range = pref_fixed_mdm_range_default;
end
if(~exist('var_mdm_range','var'))
    warning(['Range of acceptable variable-radius mdm values not specified. Setting default to ',num2str(pref_var_mdm_range_default)]);
    var_mdm_range = pref_var_mdm_range_default;
end
if(~exist('just_count','var'))
    just_count = 0;
end

% Check to see if this is a single MRIdata, or multiple MRIdata
num_datasets = length(MRIdata);

% Check to see if fixed radii are consistent
fixed_radius = MRIdata(1).PDQ.detection_template.radius;
all_fixed_radii_equal = 1;
% for i=1:num_datasets
%     if(fixed_radius ~= MRIdata(i).PDQ.detection_template.radius)
%         all_fixed_radii_equal = 0;
%     end
% end

% For each MRIdata, check to ensure they've been analyzed by PDQ already
for i = 1:num_datasets
    if(~isfield(MRIdata(i),'PDQ'))
        error('visualize_PDQ: Provided MRIdata has no PDQ field. PDQ needs to be run first. Quitting.');
    end
end

% Append all PDQ vectors
dipoles = [];
dipoles_multiradii = [];
for i = 1:num_datasets
    dipoles = [dipoles MRIdata(i).PDQ.dipoles_hp];
    dipoles_multiradii = [dipoles_multiradii MRIdata(i).PDQ.dipoles_multiradii_hp];
end

% Remove dipoles from both lists that are below xcorr_cutoffs, and those that aren't in mdm range
dipoles            = vectorize_and_remove_zeros(dipoles,            'dipoles', fixed_xcorr_cutoff, fixed_mdm_range);
dipoles_multiradii = vectorize_and_remove_zeros(dipoles_multiradii, 'dipoles', var_xcorr_cutoff,   var_mdm_range);

% Store lengths, sizes, and screen size
num_dipoles = length(dipoles);
num_dipoles_multiradii = length(dipoles_multiradii);
scrsz = get(0,'ScreenSize');



%% If there's only one dataset, show spatial distribution of dipoles
if(~just_count)
    if(num_datasets == 1)
        data_dims = size(MRIdata(1).mask);
        x_dim = data_dims(1);
        y_dim = data_dims(2);
        z_dim = data_dims(3);
        
        [aa bb cc] = ndgrid(1:x_dim,1:y_dim,1:z_dim);
        x_grid = aa .* MRIdata.mask;
        y_grid = bb .* MRIdata.mask;
        z_grid = cc .* MRIdata.mask;
        
        warning off MATLAB:divideByZero
        % Volume sampling for normalization purposes
        [sample_dist_x sdx] = hist([vectorize_and_remove_zeros(x_grid(:                     , 1:ceil(y_dim/10):y_dim, 1:ceil(z_dim/10):z_dim)) ones(1,1e4) x_dim * ones(1,1e4)], 30);
        [sample_dist_y sdy] = hist([vectorize_and_remove_zeros(y_grid(1:ceil(x_dim/10):x_dim, :                     , 1:ceil(z_dim/10):z_dim)) ones(1,1e4) y_dim * ones(1,1e4)], 30);
        [sample_dist_z sdz] = hist([vectorize_and_remove_zeros(z_grid(1:ceil(x_dim/10):x_dim, 1:ceil(y_dim/10):y_dim, :                     )) ones(1,1e4) z_dim * ones(1,1e4)], 30);
        
        % Dipole sampling
        [dipole_dist_x ddx] = hist([dipoles.x 1 x_dim], 30);
        [dipole_dist_y ddy] = hist([dipoles.y 1 y_dim], 30);
        [dipole_dist_z ddz] = hist([dipoles.z 1 z_dim], 30);
        
        % Create figures
        figure('Position',[10 scrsz(4)/2 scrsz(3)/2.2 scrsz(4)/2.5]);
        bar(sdx, dipole_dist_x ./ sample_dist_x);
        title('Dipole distribution in x-dimension'); xlabel('X Coordinate'); ylabel('# of dipoles with X Coordinate');
        
        figure('Position',[10 scrsz(4)/2 scrsz(3)/2.2 scrsz(4)/2.5]);
        bar(sdy, dipole_dist_y ./ sample_dist_y);
        title('Dipole distribution in y-dimension'); xlabel('Y Coordinate'); ylabel('# of dipoles with Y Coordinate');
        
        figure('Position',[10 scrsz(4)/2 scrsz(3)/2.2 scrsz(4)/2.5]);
        bar(sdz, dipole_dist_z ./ sample_dist_z);
        title('Dipole distribution in z-dimension'); xlabel('Z Coordinate'); ylabel('# of Dipoles with Z Coordinate');
        
        warning on MATLAB:divideByZero
    end
end

if(~just_count && all_fixed_radii_equal)
    %% If all datasets share the same radius, show fixed_radius visuals
    % XCORR
    figure('Position',[10 20 scrsz(3)/2.2 scrsz(4)/2.5]);
    hist([dipoles.xcorr],  30);
    title(['XCORR value for fixed radius of ',num2str(fixed_radius),' microns']); xlabel('XCORR (unitless)'); ylabel('# of Dipoles');
    
    % MDM SPHERE
    figure('Position',[10 20 scrsz(3)/2.2 scrsz(4)/2.5]);
        export_m_sphere = [dipoles.m_sphere]';
        save export_m_sphere export_m_sphere
    hist([dipoles.m_sphere],   30);
    title(['Sphere magnetic dipole moment for fixed radius of ',num2str(fixed_radius),' microns']); xlabel('Magnetic dipole moment (A*m^2)'); ylabel('# of Dipoles');
    
    % MDM CORE
    figure('Position',[10 20 scrsz(3)/2.2 scrsz(4)/2.5]);
    export_m_core = [dipoles.m_core]';
    save export_m_core export_m_core
    hist(export_m_core,   30);
    title(['Core magnetic dipole moment for fixed radius of ',num2str(fixed_radius),' microns']); xlabel('Magnetic dipole moment (A*m^2)'); ylabel('# of Dipoles');
    
    % SUSCEPT SPHERE
    figure('Position',[10 20 scrsz(3)/2.2 scrsz(4)/2.5]);
    hist([dipoles.suscept_sphere],  30);
    title(['Sphere magnetic susceptibility for fixed radius of ',num2str(fixed_radius),' microns']); xlabel('Susceptibility (unitless)'); ylabel('# of Dipoles');
    
    % SUSCEPT CORE
    figure('Position',[10 20 scrsz(3)/2.2 scrsz(4)/2.5]);
    hist([dipoles.suscept_core],  30);
    title(['Core magnetic susceptibility for fixed radius of ',num2str(fixed_radius),' microns']); xlabel('Susceptibility (unitless)'); ylabel('# of Dipoles');
    
    % VOLUME MAGNETIZATION SPHERE
    figure('Position',[10 20 scrsz(3)/2.2 scrsz(4)/2.5]);
    hist([dipoles.m_vol_sphere],  60);
    title(['Sphere volume magnetization for fixed radius of ',num2str(fixed_radius),' microns']); xlabel('Volume Magnetization (A/m)'); ylabel('# of Dipoles');
    
    % VOLUME MAGNETIZATION CORE
    figure('Position',[10 20 scrsz(3)/2.2 scrsz(4)/2.5]);
    hist([dipoles.m_vol_core],  60);
    title(['Core volume magnetization for fixed radius of ',num2str(fixed_radius),' microns']); xlabel('Volume Magnetization (A/m)'); ylabel('# of Dipoles');
    
    % MDM vs XCORR SPHERE
    figure('Position',[10 20 scrsz(3)/2.2 scrsz(4)/2.5]);
    scatter([dipoles.xcorr], [dipoles.m_sphere],'.','blue');
    title(['Fixed radius (',num2str(fixed_radius),' microns): Sphere magnetic dipole moment VS XCORR value']); xlabel('XCORR (0.0, 1.0)'); ylabel('Magnetic dipole moment (A*m^2)');
    
    % MDM vs XCORR CORE
    figure('Position',[10 20 scrsz(3)/2.2 scrsz(4)/2.5]);
    export_xcorr = [dipoles.xcorr];
    scatter([dipoles.xcorr], [dipoles.m_core],'.','blue');
    title(['Fixed radius (',num2str(fixed_radius),' microns): Core magnetic dipole moment VS XCORR value']); xlabel('XCORR (0.0, 1.0)'); ylabel('Magnetic dipole moment (A*m^2)');
    save export_xcorr export_xcorr
    
end


%% Fixed radius Calculations
if(all_fixed_radii_equal && ~just_count)
    [fixed_m_sphere_mean       fixed_m_sphere_std]       = gaussian_fit([dipoles.m_sphere],       30, pref_display_gaussian_fits);
    [fixed_m_core_mean         fixed_m_core_std]         = gaussian_fit([dipoles.m_core],         30, pref_display_gaussian_fits);
    [fixed_suscept_sphere_mean fixed_suscept_sphere_std] = gaussian_fit([dipoles.suscept_sphere], 30, pref_display_gaussian_fits);
    [fixed_suscept_core_mean   fixed_suscept_core_std]   = gaussian_fit([dipoles.suscept_core],   30, pref_display_gaussian_fits);
    [fixed_m_vol_sphere_mean   fixed_m_vol_sphere_std]   = gaussian_fit([dipoles.m_vol_sphere],   30, pref_display_gaussian_fits);
    [fixed_m_vol_core_mean     fixed_m_vol_core_std]     = gaussian_fit([dipoles.m_vol_core],     30, pref_display_gaussian_fits);
end

if(~just_count && (~isempty([dipoles_multiradii])))
    %% Variable radius visuals
    % RADIUS
    figure('Position',[scrsz(3)/2 20 scrsz(3)/2.2 scrsz(4)/2.5]);
    hist([dipoles_multiradii.radius],       20);
    title('Sphere estimated radius'); xlabel('Radius (microns)'); ylabel('# of Dipoles');
    
    % XCORR
    figure('Position',[scrsz(3)/2 20 scrsz(3)/2.2 scrsz(4)/2.5]);
    hist([dipoles_multiradii.xcorr],  30);
    title('XCORR for variable radius'); xlabel('XCORR (unitless)'); ylabel('# of Dipoles');
    
    % MDM SPHERE
    figure('Position',[scrsz(3)/2 20 scrsz(3)/2.2 scrsz(4)/2.5]);
    hist([dipoles_multiradii.m_sphere],   40);
    title('Sphere magnetic dipole moment for variable radius'); xlabel('Magnetic dipole moment  (A*m^2)'); ylabel('# of Dipoles');
    
    % MDM CORE
    figure('Position',[scrsz(3)/2 20 scrsz(3)/2.2 scrsz(4)/2.5]);
    hist([dipoles_multiradii.m_core],   40);
    title('Core magnetic dipole moment for variable radius'); xlabel('Magnetic dipole moment (A*m^2)'); ylabel('# of Dipoles');
    
    % SUSCEPT SPHERE
    figure('Position',[scrsz(3)/2 20 scrsz(3)/2.2 scrsz(4)/2.5]);
    hist([dipoles_multiradii.suscept_sphere],  60);
    title('Sphere magnetic susceptibility for variable radius'); xlabel('Susceptibility (unitless)'); ylabel('# of Dipoles');
    
    % SUSCEPT CORE
    figure('Position',[scrsz(3)/2 20 scrsz(3)/2.2 scrsz(4)/2.5]);
    hist([dipoles_multiradii.suscept_core],  60);
    title('Core magnetic susceptibility for variable radius'); xlabel('Susceptibility (unitless)'); ylabel('# of Dipoles');
    
    % VOLUME MAGNETIZATION SPHERE
    figure('Position',[scrsz(3)/2 20 scrsz(3)/2.2 scrsz(4)/2.5]);
    hist([dipoles_multiradii.m_vol_sphere],  60);
    title(['Sphere volume magnetization for variable radius']); xlabel('Volume Magnetization (A/m)'); ylabel('# of Dipoles');
    
    % VOLUME MAGNETIZATION CORE
    figure('Position',[scrsz(3)/2 20 scrsz(3)/2.2 scrsz(4)/2.5]);
    hist([dipoles_multiradii.m_vol_core],  60);
    title(['Core volume magnetization for variable radius']); xlabel('Volume Magnetization (A/m)'); ylabel('# of Dipoles');
    
    % MDM vs RADIUS
    figure('Position',[scrsz(3)/2 20 scrsz(3)/2.2 scrsz(4)/2.5]);
    scatter([dipoles_multiradii.radius], [dipoles_multiradii.m_core], '.','blue');
    title('Variable radius: Core magnetic dipole moment VS radius'); xlabel('Radius (microns)'); ylabel('Magnetic dipole moment (A*m^2)');
    
    % Sphere MDM vs XCORR
    figure('Position',[scrsz(3)/2 20 scrsz(3)/2.2 scrsz(4)/2.5]);
    scatter([dipoles_multiradii.xcorr], [dipoles_multiradii.m_sphere],'.','blue');
    title('Variable radius: Sphere magnetic dipole moment VS XCORR value'); xlabel('XCORR value [0.0, 1.0]'); ylabel('Magnetic dipole moment  (A*m^2)');
    
    % Core MDM vs XCORR
    figure('Position',[scrsz(3)/2 20 scrsz(3)/2.2 scrsz(4)/2.5]);
    scatter([dipoles_multiradii.xcorr], [dipoles_multiradii.m_core],'.','blue');
    title('Variable radius: Core magnetic dipole moment VS XCORR value'); xlabel('XCORR value [0.0, 1.0]'); ylabel('Magnetic dipole moment  (A*m^2)');
    
    % RADIUS vs XCORR
    figure('Position',[scrsz(3)/2 20 scrsz(3)/2.2 scrsz(4)/2.5]);
    scatter([dipoles_multiradii.xcorr], [dipoles_multiradii.radius],'.','blue');
    title('Variable radius: Radius VS XCORR value'); xlabel('XCORR value (0.0, 1.0)'); ylabel('Radius (microns)');
end


%% Variable radius calculations
if(~just_count && (~isempty([dipoles_multiradii])))
    [var_radius_mean         var_radius_std]         = gaussian_fit([dipoles_multiradii.radius],         30, pref_display_gaussian_fits);
    [var_xcorr_mean          var_xcorr_std]          = gaussian_fit([dipoles_multiradii.xcorr],          30, pref_display_gaussian_fits);
    [var_mdm_mean            var_mdm_std]            = gaussian_fit([dipoles_multiradii.m_sphere],       30, pref_display_gaussian_fits);
    [var_m_core_mean         var_m_core_std]         = gaussian_fit([dipoles_multiradii.m_core],         30, pref_display_gaussian_fits);
    [var_suscept_sphere_mean var_suscept_sphere_std] = gaussian_fit([dipoles_multiradii.suscept_sphere], 30, pref_display_gaussian_fits);
    [var_suscept_core_mean   var_suscept_core_std]   = gaussian_fit([dipoles_multiradii.suscept_core],   30, pref_display_gaussian_fits);
    [var_m_vol_sphere_mean   var_m_vol_sphere_std]   = gaussian_fit([dipoles_multiradii.m_vol_sphere],   30, pref_display_gaussian_fits);
    [var_m_vol_core_mean     var_m_vol_core_std]     = gaussian_fit([dipoles_multiradii.m_vol_core],     30, pref_display_gaussian_fits);
end


%% Neighbor-related Visualization
if(~just_count)
    figure('Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/2.2 scrsz(4)/2.5]);
    scatter([dipoles.num_neighbors] + 1,[dipoles.xcorr]);
    hold on
    %plot(sum_xcorr_per_bin)
    xlabel('Number of Neighbors INCLUDING self'); ylabel('Average XCORR value');
    hold off
    
    % (Histogram)  X_neighbors VS num_dipoles_with_X_neighbors
    % (Image)      Consensus dipoles from different numbers of neighbors
    if(max([dipoles.num_neighbors]))
        figure('Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/2.2 scrsz(4)/2.5]);
        hist([dipoles.num_neighbors], max([dipoles.num_neighbors])); xlabel('X Neighbors'); ylabel('Number of Dipoles with X Neighbors within Template Dimensions');
    end
    
    if(num_datasets == 1)
        for i = 1:max([dipoles.num_neighbors])
            images(MRIdata.PDQ.neighbor_consensus(:,:,:,i)); title(['Consensus of All Dipoles with ',num2str(i),' Neighbors.']);
        end
    end
end


%% Mahalanobis Furthest-Distance Clustering
% max_clusts = 7;
% x_dim_template = template_dims(1);
% y_dim_template = template_dims(2);
% z_dim_template = template_dims(3);
%
% % Calculate mahalanobis distances
% %distance_mahal = pdist( [[dipoles_multiradii.radius]; ([dipoles_multiradii.m_sphere])].', 'mahalanobis');
% distance_mahal = pdist( [dipoles_multiradii.m_sphere ]  .', 'mahalanobis');
%
% % Create linkage based on distances
% linkage_mahal_complete = linkage(distance_mahal, 'complete');
%
% % Create clusters from linkages
% for i = 2:max_clusts
%     clusts(:,i) = cluster(linkage_mahal_complete,'maxclust',i);
% end
%
% % Store clusters into dipole data structure
% for i = 1:min(length(dipoles), length(clusts(:,2)) )
%     dipoles_multiradii(i).cluster_2 = clusts(i,2);
%     dipoles_multiradii(i).cluster_3 = clusts(i,3);
%     dipoles_multiradii(i).cluster_4 = clusts(i,4);
%     dipoles_multiradii(i).cluster_5 = clusts(i,5);
%     dipoles_multiradii(i).cluster_6 = clusts(i,6);
%     dipoles_multiradii(i).cluster_7 = clusts(i,7);
% end

% Prepare
% cluster_consensus_2 = single(zeros(x_dim_template, y_dim_template, z_dim_template, 2));
% cluster_consensus_3 = single(zeros(x_dim_template, y_dim_template, z_dim_template, 3));
% cluster_consensus_4 = single(zeros(x_dim_template, y_dim_template, z_dim_template, 4));
% cluster_consensus_5 = single(zeros(x_dim_template, y_dim_template, z_dim_template, 5));
%
% for i = 1:2
%     count = 0;
%     for j = 1:length(dipoles)
%         if (dipoles(j).cluster_2 == i)
%             cluster_consensus_2(:,:,:,i) = cluster_consensus_2(:,:,:,i) + dipoles(j).phase;
%             count = count + 1;
%         end
%     end
%     cluster_consensus_2(:,:,:,i) = cluster_consensus_2(:,:,:,i)/count;
%     images(cluster_consensus_2(:,:,:,i)); title(['2-Cluster Furthest-Distance Mahalanobis. Dipole consensus from ',num2str(count),' Dipoles in Cluster #',num2str(i)]);
% end
%
% for i = 1:3
%     count = 0;
%     for j = 1:length(dipoles)
%         if (dipoles(j).cluster_3 == i)
%             cluster_consensus_3(:,:,:,i) = cluster_consensus_3(:,:,:,i) + dipoles(j).phase;
%             count = count + 1;
%         end
%     end
%     cluster_consensus_3(:,:,:,i) = cluster_consensus_3(:,:,:,i)/count;
%     images(cluster_consensus_3(:,:,:,i)); title(['3-Cluster Furthest-Distance Mahalanobis. Dipole consensus from ',num2str(count),' Dipoles in Cluster #',num2str(i)]);
% end
% for i = 1:4
%     count = 0;
%     for j = 1:length(dipoles)
%         if (dipoles(j).cluster_4 == i)
%             cluster_consensus_4(:,:,:,i) = cluster_consensus_4(:,:,:,i) + dipoles(j).phase;
%             count = count + 1;
%         end
%     end
%     cluster_consensus_4(:,:,:,i) = cluster_consensus_4(:,:,:,i)/count;
%     images(cluster_consensus_4(:,:,:,i)); title(['4-Cluster Furthest-Distance Mahalanobis. Dipole consensus from ',num2str(count),' Dipoles in Cluster #',num2str(i)]);
% end
% for i = 1:5
%     count = 0;
%     for j = 1:length(dipoles)
%         if (dipoles(j).cluster_5 == i)
%             cluster_consensus_5(:,:,:,i) = cluster_consensus_5(:,:,:,i) + dipoles(j).phase;
%             count = count + 1;
%         end
%     end
%     cluster_consensus_5(:,:,:,i) = cluster_consensus_5(:,:,:,i)/count;
%     images(cluster_consensus_5(:,:,:,i)); title(['5-Cluster Furthest-Distance Mahalanobis. Dipole consensus from ',num2str(count),' Dipoles in Cluster #',num2str(i)]);
% end

% % Display dendrogram
% figure('Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/2.2 scrsz(4)/2.5]);
% dendrogram(linkage_mahal_complete,5);
%
% % For each set of clusters
% for i = 1:7
%     figure('Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/2.2 scrsz(4)/2.5]);
%     scatter([dipoles_multiradii.radius]' .* (clusts(:,i)==1), ([dipoles_multiradii.m_sphere])' .* (clusts(:,i)==1) ,'r.');
%     title(['Mahalanobis Farthest-Distance Clustering: ', num2str(i), ' clusters']); xlabel('Radius'); ylabel('mdm (A*m^2)');
%     hold on
%     scatter([dipoles_multiradii.radius].' .* (clusts(:,i)==2), ([dipoles_multiradii.m_sphere]).' .* (clusts(:,i)==2),'g.');
%     scatter([dipoles_multiradii.radius].' .* (clusts(:,i)==3), ([dipoles_multiradii.m_sphere]).' .* (clusts(:,i)==3),'b.');
%     scatter([dipoles_multiradii.radius].' .* (clusts(:,i)==4), ([dipoles_multiradii.m_sphere]).' .* (clusts(:,i)==4),'k.');
%     scatter([dipoles_multiradii.radius].' .* (clusts(:,i)==5), ([dipoles_multiradii.m_sphere]).' .* (clusts(:,i)==5),'c.');
%     scatter([dipoles_multiradii.radius].' .* (clusts(:,i)==6), ([dipoles_multiradii.m_sphere]).' .* (clusts(:,i)==6),'y.');
%     scatter([dipoles_multiradii.radius].' .* (clusts(:,i)==7), ([dipoles_multiradii.m_sphere]).' .* (clusts(:,i)==7),'m.');
%     hold off
% end



%% Region Visualizations
% for curr_region = 1:num_regions
%     disp(['Dipoles in Region ',num2str(curr_region), ' = ', num2str(MRIdata.regions(curr_region).dipole_count)]);
% end
%
% volume_distribution = hist(MRIdata.random_grid,50);
% dipole_distribution = hist(MRIdata.regions.dipole_min_dist,50);
% figure; plot(dipole_distribution./volume_distribution);



%% Print dipole counts and statistics
disp('############################ Counting Result #############################');
disp(['Analyzed volume: ', num2str(PDQ_calc_volume(MRIdata(1),1 )), ' cc']);
disp([num2str(length(dipoles)),            ' dipoles found with: XCORR > ', num2str(fixed_xcorr_cutoff), ' & ', num2str(fixed_mdm_range(1)), ' < mdm < ', num2str(fixed_mdm_range(2))]);
disp([num2str(length(dipoles_multiradii)), ' dipoles found with: XCORR > ', num2str(var_xcorr_cutoff),   ' & ', num2str(var_mdm_range(1)), ' < mdm < ', num2str(var_mdm_range(2))]);
disp('##########################################################################');
disp(' ');
disp(' ');
if(~just_count)
    if(all_fixed_radii_equal)
        disp(['####### Magnetic Properties, Assuming fixed radius of ',num2str(fixed_radius), ' microns ##########']);
        disp(['# XCORR: ', num2str(mean([dipoles.xcorr])), ' +/- ', num2str(std([dipoles.xcorr]))]);
        disp(['#### Assuming Background Susceptibility Value = User-Provided Value ####']);
        disp(['# Magnetic Dipole Moment (m): ', num2str(fixed_m_sphere_mean*1e12), ' +/- ', num2str(fixed_m_sphere_std*1e12)]);
        disp(['# Susceptibility: (chi): ', num2str(fixed_suscept_sphere_mean*1e6), ' +/- ', num2str(fixed_suscept_sphere_std*1e6)]);
        disp(['# Volume Magnetization (M): ', num2str(fixed_m_vol_sphere_mean), ' +/- ', num2str(fixed_m_vol_sphere_std)]);
        disp('#');
        disp(['#### Assuming Background Susceptibility = 0.0 ##']);
        disp(['# Magnetic Dipole Moment (m): ', num2str(fixed_m_core_mean*1e12), ' +/- ', num2str(fixed_m_core_std*1e12)]);
        disp(['# Susceptibility: (chi): ', num2str(fixed_suscept_core_mean*1e6), ' +/- ', num2str(fixed_suscept_core_std*1e6)]);
        disp(['# Volume Magnetization (M): ', num2str(fixed_m_vol_core_mean), ' +/- ', num2str(fixed_m_vol_core_std)]);
        disp('###########################################################################');
        disp(' ');
        disp(' ');
    end
    if((~isempty([dipoles_multiradii])))
    disp(['########## Magnetic Properties, Allowing variable radius ##########']);
    disp(['# XCORR: ', num2str(var_xcorr_mean), ' +/- ', num2str(var_xcorr_std)]);
    disp(['# Radius: ', num2str(var_radius_mean), ' +/- ', num2str(var_radius_std)]);
    disp(['#### Assuming Background Susceptibility Value = User-Provided Value ####']);
    disp(['# Magnetic Dipole Moment (m): ', num2str(var_mdm_mean*1e12), ' +/- ', num2str(var_mdm_std*1e12)]);
    disp(['# Susceptibility: (chi): ', num2str(var_suscept_sphere_mean*1e6), ' +/- ', num2str(var_suscept_sphere_std*1e6)]);
    disp(['# Volume Magnetization (M): ', num2str(var_m_vol_sphere_mean), ' +/- ', num2str(var_m_vol_sphere_std)]);
    disp('#');
    disp(['#### Assuming Background Susceptibility = 0.0 ##']);
    disp(['# Magnetic Dipole Moment (m): ', num2str(var_m_core_mean*1e12), ' +/- ', num2str(var_m_core_std*1e12)]);
    disp(['# Susceptibility: (chi): ', num2str(var_suscept_core_mean*1e6), ' +/- ', num2str(var_suscept_core_std*1e6)]);
    disp(['# Volume Magnetization (M): ', num2str(var_m_vol_core_mean), ' +/- ', num2str(var_m_vol_core_std)]);
    disp('############################################################################');
    end
    end




%%%%%%%%%%%%%
%   EOF
%%%%%%%%%%%%%