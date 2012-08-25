%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @name dipole_analysis.m
% @author Parker Mills, Ahrens Lab, Carnegie Mellon 2009
% @brief Analyzes XCORR peaks found throughout a phase volume, creating a list of dipoles
% @trusted yes
% @commented yes
% @optimized probably
% @parallelized yes
% 
% ==================== INPUT PARAMETERS ======================
% @param phase            (2D/3D Float)        Phase dataset
% @param raw_peaks        (2D/3D Float)        Peaks found in phase dataset
% @param xcorr_cutoff     (Float)              Cutoff value for a peak to be considered to be a dipole
% @param template         (Template)           Template used to form the raw peaks
%
% ==================== RETURNED DATA =========================
% @return dipoles         (1D Dipole)          List of dipoles found in volume
%
% ==================== ASSUMPTIONS ===========================
% @assume 
%
% ==================== EXAMPLE USAGE =========================
% optimal_template = dipole_analysis(phase, raw_peaks, xcorr_cutoff, template)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function dipoles = PDQ_dipole_analysis(phase, raw_peaks, xcorr_cutoff, template, CHI_background)


%% Initializations

%TEMPORARY PREF
template_shift = sqrt(3)/4.0;

% Dimensions
phase_dims = size(phase);
template_dims = size(template.template);

% Dipole structure
dipoles = struct('x',{},'y',{},'z',{},'radius',{},'phase',{},'xcorr',{},'suscept_sphere',{},'suscept_core',{},'m_sphere',{},'m_core',{},'m_vol_sphere',{},'m_vol_core',{},'verified',{});

%Create shifted templates
for z_shift = 1:3
    for y_shift = 1:3
        for x_shift = 1:3
            optimal_templates(x_shift,y_shift,z_shift) = PDQ_generate_template(template.resolution, template.orientation, template.radius, template.B0, template.TE, template.d_Chi, [template_shift*(x_shift-2) template_shift*(y_shift-2)  template_shift*(z_shift-2)]);
        end
    end
end

%% Go through all peaks, adding them as dipoles to the dipole structure
%for j_1 = 1:phase_dims(1)
for j_1 = 1:phase_dims(1) %doesn't work with optimal templates matrix
    for j_2 = 1:phase_dims(2)
        for j_3 = 1:phase_dims(3)
            if( raw_peaks(j_1,j_2,j_3) > xcorr_cutoff ) % If peak magnitude (dipole) is over xcorr_cutoff threshold

                %Determine boundary coordinates of peak (dipole)
                low_x  = j_1  -  (template_dims(1) - 1)/2;
                high_x = j_1  +  (template_dims(1) - 1)/2;
                low_y  = j_2  -  (template_dims(2) - 1)/2;
                high_y = j_2  +  (template_dims(2) - 1)/2;
                low_z  = j_3  -  (template_dims(3) - 1)/2;
                high_z = j_3  +  (template_dims(3) - 1)/2;

                % Only address peaks (dipoles) that are not perched on the edge of the dataset
                if((low_x > 0)  &&  (low_y > 0)  &&  (low_z > 0)  &&  (high_x <= phase_dims(1))  &&  (high_y <= phase_dims(2))  &&  (high_z <= phase_dims(3)))
                    
                    % Find optimal shift
                    best_xcorr = -1.0;
                    for z_shift = 1:3
                        for y_shift = 1:3
                            for x_shift = 1:3
                                current_xcorr = real(normxcorr3(optimal_templates(x_shift,y_shift,z_shift).template, phase(low_x:high_x, low_y:high_y, low_z:high_z), 'valid'));
                                if(current_xcorr > best_xcorr)
                                    template_optimal = optimal_templates(x_shift,y_shift,z_shift);
                                    best_xcorr = current_xcorr;
                                end
                            end
                        end
                    end
                    
                    
                    
                    % Store information about this dipole in dipole listing
                    dipole = struct;
                    dipole.x = j_1;
                    dipole.y = j_2;
                    dipole.z = j_3;
                    dipole.radius = template_optimal.radius;
                    dipole.phase = phase(low_x:high_x, low_y:high_y, low_z:high_z); % Store dipole phase impression
                    dipole.xcorr = best_xcorr; % Store XCORR value
                    
                    % Calculate dipole's magnetic properties
                    radius_3 = (dipole.radius * 1e-6)^3.0;
                    volume = 4/3 * pi * radius_3;
                    dipole.suscept_sphere = (template_optimal.d_Chi * minsumsquares(phase(low_x:high_x, low_y:high_y, low_z:high_z), template_optimal.template))   +   CHI_background;
                    dipole.suscept_core   = (template_optimal.d_Chi * minsumsquares(phase(low_x:high_x, low_y:high_y, low_z:high_z) ,template_optimal.template))                     ;
                    dipole.m_sphere = 1e7 * (template_optimal.B0/3) * radius_3 * dipole.suscept_sphere / (1 + dipole.suscept_sphere);
                    dipole.m_core   = 1e7 * (template_optimal.B0/3) * radius_3 * dipole.suscept_core   / (1 + dipole.suscept_core  );
                    dipole.m_vol_sphere = dipole.m_sphere / volume;
                    dipole.m_vol_core   = dipole.m_core / volume;
                    
                    % Verification structure
                    dipole.verified = 0;
                    
                    % Append this dipole to dipole vector
                    dipoles = [dipoles dipole];
                end
            end
        end
    end
end



%%%%%%%%%%%%%
%    EOF
%%%%%%%%%%%%%