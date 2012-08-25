%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @name PDQ_register.m
% @author Parker Mills, Ahrens Lab, Carnegie Mellon 2010
% @author T. Kevin Hitchens, Pittsburgh NMR Center for Biomedical Research, Carnegie Mellon 2010
% @brief Performs rigid registration, dipole registration, and dipole movement measurement on a set of MRI datasets
% @trusted no
% @robust no
% @commented no
% @optimized no
% @parallelized no
%
% ==================== INPUT PARAMETERS ======================
% @param MRIdataSeries          (1D MRIdata)       Datasets to be registered across time
% @param max_speed              (Float)            Maximum speed moved in microns/hour
%                                                  (e.g., 40 microns/hour = 960 microns/day)
% @param max_mdm_change         (Float)            Cutoff in maximum mdm change
%                                                  (e.g., Found to be 0.4e-12 A m^2 for in vivo mouse brains at 256x190x190 resoluton)
%
% ==================== RETURNED DATA =========================
% @return registered            (1D MRIdata)       Series of co-registered datasets
% @return registration_data     (Registration)     Structure containing registration information between datasets
%
% ==================== ASSUMPTIONS ===========================
% @assume 3D volumes are already roughly registered
% @assume 3D volumes have same dimensions
% @assume All MRIdata files must contain a field entitled 'date', which is in units of hours
%
% ==================== EXAMPLE USAGE =========================
% [TBI_brains_registered reg_data] = PDQ_register(TBI_brains,50, 0.4e-12); % 50 microns/hour
% [registered reg_data] = PDQ_register([data1 data2], 35, 0.4e-12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [registered registration_data] = PDQ_register(MRIdataSeries, max_speed, max_mdm_change)

%% Checks & Initializations

% Start parallelization
%matlabpool close force
if(matlabpool('size') == 0)
    matlabpool open
end

% Hard-Coded Preferences
pref_show_registration_steps = 0;
pref_registration_basis = 'magnitude'; % 'magnitude' for registration to be on mag images, 'mask' for it to be done on image masks

% Global Vars
timepoints = length(MRIdataSeries);


%% Sanity Checks

% Check there's more than one dataset in the series
if(timepoints < 2)
    error('PDQ_register: Require more than one dataset. Quitting.');
end

% Compare different timepoints
for i=1:timepoints-1
    
    % Check datasets have same dimensions
    if(size(MRIdataSeries(i).mag) ~= size(MRIdataSeries(i+1).mag))
        error('PDQ_register: Input MRIdata datasets have different image dimensions. Quitting.');
    end
    
    % Check datasets have same resolutions
    if(MRIdataSeries(i).resolution ~= MRIdataSeries(i+1).resolution)
        error('PDQ_register: Input MRIdata datasets have different image resolutions. Quitting.');
    end
    
    % Check datasets have different time stamps
    if(MRIdataSeries(i).date == MRIdataSeries(i+1).date)
        error('PDQ_register: Input MRIdata datasets have same time stamp. Quitting.');
    end
    
    % Check datasets have dipole verity field marked
    if(~isfield(MRIdataSeries(i).PDQ.dipoles_hp,'verity'))
        error('PDQ_register: Input MRIdata datasets need dipole verity set. Run PDQ_inspect first.');
    end
end




%% Perform rigid transformations registration between dataset pairs, going from first dataset to last dataset
registered(1) = MRIdataSeries(1);
for i = 2:timepoints
    registered(i) = register(registered(i-1), MRIdataSeries(i), pref_registration_basis, pref_show_registration_steps);
end


%% Disqualify dipoles outside of anatomy
for i = 1:timepoints
    registered(i) = PDQ_inspect(registered(i),0.3,[0 3e-12],0,0);
end




%% Prepare dipoles to be registered
for timepoint = 1:timepoints-1
    
    % Initialize distance matrices
    distance_matrix_size = max(length(registered(timepoint).PDQ.dipoles_hp), length(registered(timepoint+1).PDQ.dipoles_hp));
    empty_distance_matrix = inf .* ones(distance_matrix_size, distance_matrix_size);
    theta = empty_distance_matrix;
    phi = empty_distance_matrix;
    r = empty_distance_matrix;
    mdm_change = empty_distance_matrix;
    xcorr_change = empty_distance_matrix;
    
    for i = 1:length(registered(timepoint).PDQ.dipoles_hp) % Early Dipoles
        for j = 1:length(registered(timepoint+1).PDQ.dipoles_hp) % Late Dipoles
            p1 = registered(timepoint).PDQ.dipoles_hp(i);
            p2 = registered(timepoint+1).PDQ.dipoles_hp(j);
            
            % Calculate and store all changes ("distances")
            spatial_change.x = (p2.x - p1.x) * registered(1).resolution(1);
            spatial_change.y = (p2.y - p1.y) * registered(1).resolution(2);
            spatial_change.z = (p2.z - p1.z) * registered(1).resolution(3);
            % tbi_initial_distance(i,j) = sqrt((tbi_p(1)-p1.x)^2 + (tbi_p(2)-p1.y)^2 + (tbi_p(3)-p1.z)^2);
            % tbi_final_distance(i,j) = sqrt((tbi_p(1)-p2.x)^2 + (tbi_p(2)-p2.y)^2 + (tbi_p(3)-p2.z)^2);
            % tbi_movement(i,j) = tbi_final_distance(i,j) - tbi_initial_distance(i,j);
            mdm_change(i,j) = p2.m_core - p1.m_core;
            xcorr_change(i,j) = p2.xcorr - p1.xcorr;
            [theta_p, phi_p, r_p] = cart2sph(spatial_change.x, spatial_change.y, spatial_change.z);
            theta(i,j) = theta_p;
            phi(i,j) = phi_p;
            
            % Disqualify distant dipoles
            max_distance = max_speed * (registered(timepoint+1).date - registered(timepoint).date);
            if(r_p > max_distance)
                r(i,j) = inf;
            else
                r(i,j) = r_p;
            end
            
            % Disqualify dipoles with large mdm changes
            if (abs(mdm_change(i,j)) > max_mdm_change)
                mdm_change(i,j) = inf;
            end
            
            % Disqualify invalid dipoles
            if(strcmp(p1.verity,'invalid') || strcmp(p2.verity,'invalid'))
                r(i,j) = inf;
            end
            
        end
    end
    
    weighted_dist_matrix = r/max_distance + mdm_change/max_mdm_change;
    
    % Reverse sign of distance matrix for processing by Hungarian Algorithm, since it tries to maximize, not minimize
    weighted_dist_matrix = -weighted_dist_matrix;
    
    % Feed distance matrix to Hungarian Algorithm to find matching dipoles
    matches = libmaxmatching(weighted_dist_matrix);
    
    
    
    %% Print out results for interpretation by user
    export_matrix = zeros(max(length(registered(timepoint).PDQ.dipoles_hp),length(registered(timepoint+1).PDQ.dipoles_hp)));
    for i=1:length(registered(timepoint).PDQ.dipoles_hp) % Early Dipoles
        for j=1:length(registered(timepoint+1).PDQ.dipoles_hp) % Late Dipoles
            
            % Mark dipole as not matched if the distance equals infinity
            if(matches(i,j) && (r(i,j) ~= inf))
                disp(['Dipole ',num2str(j),' matches Early Dipole ',num2str(i),'. Movement (microns): ',num2str(r(i,j)),...
                    ', Theta: ',num2str(theta(i,j)),', Phi: ',num2str(phi(i,j)),...
                    ', mdm_change: ',num2str(1e12 * abs(mdm_change(i,j))),', xcorr_change: ',num2str(abs(xcorr_change(i,j)))]);
                export_matrix(i,1:3) = [registered(timepoint).PDQ.dipoles_hp(i).x registered(timepoint).PDQ.dipoles_hp(i).y registered(timepoint).PDQ.dipoles_hp(i).z];
                export_matrix(i,4:6) = [registered(timepoint+1).PDQ.dipoles_hp(j).x registered(timepoint+1).PDQ.dipoles_hp(j).y registered(timepoint+1).PDQ.dipoles_hp(j).z];
                export_matrix(i,7) = registered(timepoint).PDQ.dipoles_hp(i).m_core;
            end
        end
    end
    csvwrite(['PDQ',num2str(timepoint),'.csv'], export_matrix);
    
    
    
    %% Store results in Registration structure
    registration_data(timepoint).matches = matches;
    registration_data(timepoint).theta = theta;
    registration_data(timepoint).phi = phi;
    registration_data(timepoint).r = r;
    registration_data(timepoint).mdm_change = mdm_change;
    registration_data(timepoint).xcorr_change = xcorr_change;
    registration_data(timepoint).dipoles = matches .* (r ~= inf);
end

curr_point = 1;
outline = abs(diff(registered(1).mask,1,1));
for i = 1:1:size(registered(1).mask,1)-1
for j = 1:5:size(registered(1).mask,2)-1
for k = 1:5:size(registered(1).mask,3)-1
if(outline(i,j,k))
points(curr_point,1) = i;
points(curr_point,2) = j;
points(curr_point,3) = k;
curr_point = curr_point + 1;
end
end
end
end
csvwrite('mask.csv',points)


% Discover which dipoles were registered across all three datasets
%[x y] = find(registration_data(1).dipoles)
% [x y] = find(registration_data(1).dipoles)
% across3=[];
% for(i=1:length(y))
% for(j=1:length(x2))
% if(y(i) == x2(j))
% across3 = [across3 y(i)]
% end
% end
% end

%%%%%%%%%%%%%
%   EOF
%%%%%%%%%%%%%