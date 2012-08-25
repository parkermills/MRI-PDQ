%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @name calclulate_radii_and_chi.m
% @author Parker Mills, Ahrens Lab, Carnegie Mellon 2009
% @brief Calculates radii and chi for all dipoles in PDQ dipole list above xcorr_cutoff
% @trusted yes
% @robust no
% @commented no
% @optimized no
% @parallelized no
% 
% ==================== INPUT PARAMETERS ======================
% @param template           (Template)      Original template used to create the list of dipoles
% @param dipoles            (1D Dipole)     List of dipoles to be analyzed
% @param xcorr_cutoff       (Float)         Radii/susceptibility will not be calculated for dipoles with XCORR below this value
% @param template_spectrum  (3D Template)   3D matrix filled with Template datatypes, for comparison with detected dipoles 
%
% ==================== RETURNED DATA =========================
% @return dipole_list       (1D Dipole_multiradii)     List of dipoles for all radii
% 
% ==================== ASSUMPTIONS ===========================
% @assume 
%
% ==================== EXAMPLE USAGE =========================
% dipole_list = calculate_radii_and_chi(template, dipoles, xcorr_cutoff, template_spectrum)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function dipole_list = PDQ_calculate_radii_and_chi(template, dipoles, xcorr_cutoff, template_spectrum, CHI_background)


%% Preferences

pref_radii_bottom_threshold = 0.02; % Percent of lowest radii in template_spectrum that are not considered


%% Initializations

% Dipole structure
dipole_list = struct('index',{},'xcorr',{}, 'radius',{}, 'suscept_sphere',{}, 'suscept_core', {}, 'm_sphere',{}, 'm_core', {}, 'm_vol_sphere', {}, 'm_vol_core',{}, 'shift',{}, 'spectrum',{},'verified',{});



%% If a dipole is above the xcorr_cutoff, and has no neighbors
parfor i = 1:length(dipoles)
    if( (dipoles(i).xcorr > xcorr_cutoff) && (~dipoles(i).num_neighbors) )
        
        %% Fit dipole
        % Go through radii
        xcorr_list = zeros(3,3,3, length(template_spectrum(2,2,2,:)));
        index_start = ceil(pref_radii_bottom_threshold * length(template_spectrum(2,2,2,:)));
        for x_shift = 1:3
            for y_shift = 1:3
                for z_shift = 1:3
                    for radii = index_start:length(template_spectrum(2,2,2,:))
                        xcorr_list(x_shift, y_shift, z_shift, radii) = normxcorr3(dipoles(i).phase, template_spectrum(x_shift, y_shift, z_shift, radii).template, 'valid');
                    end
                end
            end
        end
        
        % Find optimal shifts
        [optimal_radius_xcorr_1, optimal_index_1] = max(xcorr_list); % Find maximum x_shift
        [optimal_radius_xcorr_2, optimal_index_2] = max(optimal_radius_xcorr_1); % Find maximum y_shift
        [optimal_radius_xcorr_3, optimal_index_3] = max(optimal_radius_xcorr_2); % Find maxiumum z_shift
        [optimal_radius_xcorr, optimal_index_4] = max(optimal_radius_xcorr_3); % Find optimal radius!
        
        % Find optimal index for the optimal shift
        matrix_index_4 = optimal_index_4;
        matrix_index_3 = optimal_index_3(:, :, :, matrix_index_4);
        matrix_index_2 = optimal_index_2(:, :, matrix_index_3,matrix_index_4);
        matrix_index_1 = optimal_index_1(:, matrix_index_2,matrix_index_3,matrix_index_4);
        
        optimal_template = template_spectrum(matrix_index_1, matrix_index_2, matrix_index_3, matrix_index_4);
        
        
        
        %% Store dipole and its various properties
        
        % Index, XCORR, etc.
        dipole = struct;
        dipole.index = i;
        dipole.xcorr = optimal_radius_xcorr;
        dipole.radius = optimal_template.radius;
        
        % Magnetic properties
        radius_3 = (dipole.radius * 1e-6)^3.0;
        volume = 4/3 * pi * radius_3;
        dipole.suscept_sphere = (optimal_template.d_Chi .* minsumsquares(dipoles(i).phase, optimal_template.template))   +   CHI_background;
        dipole.suscept_core   = (optimal_template.d_Chi .* minsumsquares(dipoles(i).phase, optimal_template.template))                     ;
        dipole.m_sphere = 1e7 * (optimal_template.B0/3) * radius_3 * dipole.suscept_sphere / (1 + dipole.suscept_sphere);
        dipole.m_core   = 1e7 * (optimal_template.B0/3) * radius_3 * dipole.suscept_core   / (1 + dipole.suscept_core  );
        dipole.m_vol_sphere = dipole.m_sphere / volume;
        dipole.m_vol_core   = dipole.m_core / volume;
        
        % Dipole radius fit information
        dipole.shift = [matrix_index_1 matrix_index_2 matrix_index_3];
        dipole.spectrum = reshape(xcorr_list(matrix_index_1, matrix_index_2, matrix_index_3,:),1,length(template_spectrum(2,2,2,:)));
        
        % Verification flag
        dipole.verified = 0;
        
        % Append to vector
        dipole_list = [dipole_list dipole];
        
    end
end



%%%%%%%%%%%%%%
%     EOF
%%%%%%%%%%%%%%