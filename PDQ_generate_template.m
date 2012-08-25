%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @name generate_PDQ_template.m
% @author Parker Mills, Ahrens Lab, Carnegie Mellon 2009
% @brief Creates 3D matrix of magnetic dipole phase impression
% @trust yes
% @robust yes
% @commented yes
% @optimized no (oversampling size may be too large)
% @parallelized no
%
% ==================== INPUT PARAMETERS ======================
%
% @param    resolution     (1D Float)       Pixel resolution (microns): [x_res y_res z_res]
% @param    orientation    (Float)          Direction of B0:  1=Up-down (Along X-axis (MATLAB))
%                                                             2=In-out  (Along Z-axis)
%                                                             3=Left-right (Along Y-axis (MATLAB))
%
% @param    radius         (Float)          Spheroid radius (microns)
% @param    B0             (Float)          Field strength (Tesla or MHz)
% @param    TE             (Float)          Echo time (ms)
% @param    d_Chi          (Float)          Volume magnetic susceptibility of spheroid's material, surrounded by background media (unitless)
% @param    shift          (1D Float)       [x_shift y_shift z_shift] Shift, in pixels, of dipole from center of template
%                                           Typically, shifts are <= 0.5, since a shift of more than that would mean the dipole is at a location one pixel away
%
% ==================== RETURNED DATA =========================
% @return   template           (template)       The desired template (for Template data structure see README.txt)
%
% ==================== ASSUMPTIONS ===========================
% @assume   B0 is assumed to be in:  Tesla if B0 <= pref_tesla_mhz_decision
%                                    MHz   if B0  > pref_tesla_mhz_decision
%
% @assume   Gyromagnetic ratio is for Hydrogen, not 19F or some other nucleus
%
% ==================== DISPLAYED PRODUCTS ====================
% @product
%
% ==================== SAMPLE USAGE ==========================
% template = generate_PDQ_template(resolution, orientation, radius, B0,  TE, d_Chi,      shift);
% template = generate_PDQ_template([51 51 98],           3, 300,   11.7, 11, 1.03e-4,  [0 0 0]);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function template = PDQ_generate_template(resolution, orientation, radius, B0, TE, d_Chi, shift)


%% User-set preferences
pref_gridsize = 7; % Size of dipole template, in pixels (e.g. default is 7, resulting in template sized 7x7x7)
% Usually set to N=7 or N=9, as dipole impression falls off like r^3. Should be set to
% smallest reasonable value so that the template doesn't factor in phase perturbations in other areas!

pref_invert_dipole_template = 1; % For some reason, the entire dipole template must be inverted.

pref_samples = 10000; % Samples per pixel in template. 25,000 for excellent sampling. 10,000 is decent. 2,500 is good for rough sampling. Assuming gridsize of 7x7x7, total_samples = 343 * pref_samples;

pref_print_magnetic_information = 0; % Print out template's magnetic specifications

pref_display_figures = 0; % Visualize template by showing figures for user

pref_model_inside = 'nothing'; % Model Bz inside the sphere of SPIO? Set to 'griffiths' or 'nothing'
% Signal is usually too low inside dipoles to model their inside.
% Also, this is set only if the material is homogeneous.

gamma = 2 * pi * 42576000; % Gyromagnetic ratio (Default is for Hydrogen = 42,576,000 (2*pi*Hz)/T )
gamma_2pi = 42576000; % Gyromagnetic ratio (Default is for Hydrogen = 42,576,000 Hz/T )


% Technical settings
pref_tesla_mhz_decision = 24; % B0 assumed to be in Tesla if below this value, MHz if above this value



%% Set physical constants
mu_0 = 4 * pi * 1e-7; % Permeability of free space (Tesla * meter / Amp)



%% Sanity checks
if(size(resolution)~=3) % Check resolution
    error('generate_PDQ_template: Resolution must be vector with 3 elements. Quitting.');
end

% Check a, B0, TE, d_Chi are provided
if(~exist('radius','var'))
    error('generate_PDQ_template: radius must be provided. Quitting');
end
if(~exist('B0','var'))
    error('generate_PDQ_template: B0 (MHz or T) must be provided. Quitting');
end
if(~exist('TE','var'))
    error('generate_PDQ_template: TE (milliseconds) must be provided. Quitting');
end
if(~exist('d_Chi','var'))
    error('generate_PDQ_template: d_Chi must be provided. Quitting');
end




%% Convert all units into meters, seconds, Tesla
resolution = resolution ./ 1e6; % Convert resolution from microns to meters
TE = TE ./ 1000.0; % Convert TE from milliseconds into seconds
radius = radius ./ 1e6; % Convert radius from microns to meters

% If B0 is in MHz, convert it to Tesla.
if(abs(B0) > pref_tesla_mhz_decision)
    B0 = B0 .* 1e6 ./ gamma_2pi;
end

% If we want to invert the dipole
if(pref_invert_dipole_template)
    if(B0 > 0.0)
        B0 = -B0;
    end
end



%% Perform magnetic calculations related to template and its surrounding environment
volume = (4.0 / 3.0) * pi * radius ^ 3.0; % Calculate volume of a sphere (m^3)
mdm = 1e7 * B0 * (radius^3.0 / 3) * (d_Chi /(1 + d_Chi)); % Magnetic Dipole Moment (Amp * meter^2)
M_vol  = mdm / volume; % ASSUMES CHI_BACKGROUND = 0. Bulk magnetization, magnetic dipole moment per unit volume, (A/m)



%% If preferred, print out all the magnetic calculations
if(pref_print_magnetic_information)
    disp('You have created a dipole with:');
    disp(['Volume susceptibility (dimensionless SI units): ',num2str(d_Chi)]);
    disp(['Volume: ',num2str(volume),' (m^3) or ', num2str(volume * 1e18),' (microns^3)']);
    disp(['Bulk magnetization: ',num2str(M_vol),' A/m']);
    disp(['Magnetic dipole moment: ',num2str(mdm),' A*m^2']);
    disp(' ');
end



%% Create grid representing dipole in space. X, Y, Z values are in meters
% For each voxel, generate a high-resolution grid containing 'pref_sampling_multiplier' points, then average to make template
lower_bound_x = ((-pref_gridsize + 2.0 * shift(1)) * resolution(1)) / 2.0;
lower_bound_y = ((-pref_gridsize + 2.0 * shift(2)) * resolution(2)) / 2.0;
lower_bound_z = ((-pref_gridsize + 2.0 * shift(3)) * resolution(3)) / 2.0;



%% Calculate magnetic field inside sphere (OBSOLETE FOR THIS PROJECT, MAY BE OF VALUE IN FUTURE)
% Math note:
% Bz = (2/3) * mu_0 * M     (Griffiths, SI units)

% Set Bz inside of sphere to user preference
switch pref_model_inside
    case 'griffiths'
        inside_Bz = mu_0 .* (2/3) .* M_vol;   % Bz inside of sphere if it was a homogeneous media (Tesla)
    case 'nothing'
        inside_Bz = 0.0;
end

template.template = zeros(pref_gridsize, pref_gridsize, pref_gridsize); % Allocate template
for j_1 = 1:pref_gridsize
    x_min = lower_bound_x + (j_1 - 1) * resolution(1);
    for j_2 = 1:pref_gridsize
        y_min = lower_bound_y + (j_2 - 1) * resolution(2);
        for j_3 = 1:pref_gridsize
            z_min = lower_bound_z + (j_3 - 1) * resolution(3);
            
            samples = ceil(pref_samples / (abs(j_1-3.5)*abs(j_2-3.5)*abs(j_3-3.5)));
            
            for j_4 = 1:samples
                x = x_min + rand * resolution(1);
                y = y_min + rand * resolution(2);
                z = z_min + rand * resolution(3);
                radius2 = x^2 + y^2 + z^2;
                % Does not yet factor in ANYTHING about the material (e.g., dCHI, B0, 1/3, or a^3)
                % Math note:
                % r = sqrt(x.^2 + y.^2 + z.^2)
                % cos(theta).^2 = z.^2 / r.^2
                dipole = (3 * ((z^2)/radius2) - 1) / ( sqrt(radius2)^3 );
                spin = gamma * TE * dipole * d_Chi * (radius .^ 3.0) * B0 / 3.0;
                switch orientation
                    case 1
                        if(sqrt(radius2) >= radius )%&& abs(spin) < 4*pi)
                            template.template(j_3,j_2,j_1) = template.template(j_3,j_2,j_1) + spin;
                        else
                            template.template(j_3,j_2,j_1) = template.template(j_3,j_2,j_1) + gamma * TE * inside_Bz; % Units: 1/(T*s) * s * T
                        end
                    case 2
                        if(sqrt(radius2) >= radius )%&& abs(spin) < 4*pi)
                            template.template(j_1,j_2,j_3) = template.template(j_1,j_2,j_3) + spin;
                        else
                            template.template(j_1,j_2,j_3) = template.template(j_1,j_2,j_3) + gamma * TE * inside_Bz; % Units: 1/(T*s) * s * T
                        end
                    otherwise
                        if(sqrt(radius2) >= radius )%&& abs(spin) < 4*pi)
                            template.template(j_1,j_3,j_2) = template.template(j_1,j_3,j_2) + spin;
                        else
                            template.template(j_1,j_3,j_2) = template.template(j_1,j_3,j_2) + gamma * TE * inside_Bz; % Units: 1/(T*s) * s * T
                        end
                end
            end
            
            switch orientation
                case 1
                    template.template(j_3,j_2,j_1) = template.template(j_3,j_2,j_1) ./ samples;
                case 2
                    template.template(j_1,j_2,j_3) = template.template(j_1,j_2,j_3) ./ samples;
                otherwise
                    template.template(j_1,j_3,j_2) = template.template(j_1,j_3,j_2) ./ samples;
            end
            
            % ADVANCED: Possible future feature: If susceptibility is WELL-KNOWN, then one can normalize for phase unwrapping (i.e., can't have phase unwrapped beyond 2*pi)
            %              template.template = template.template + 2 * pi;
            %              template.template = mod(template.template, 4*pi);
            %              template.template = template.template - 2 * pi;
            %template.template(4,4,4) = 0;
        end
    end
end



%% Export template and template metadata so its parameters are known
% Convert back to microns, milliseconds from meters, seconds
template.template = single(template.template);
template.resolution = resolution .* 1e6;
template.orientation = orientation;
template.radius = radius .* 1e6; % Convert meters back into microns
if(pref_invert_dipole_template) % Correct for inversion
    template.B0 = -B0;
else
    template.B0 = B0;
end
template.TE = TE .* 1000;
template.d_Chi = d_Chi;
template.model_inside = pref_model_inside;


%% Visualize high-resolution dipole model and central slice of template
if(pref_display_figures)
    % 2D Visualization
    images(template.template(:,:,4));
    
    % 3D Visualization
    %     figure;
    %     p1 = patch(isosurface(template.template,-max(max(max(template.template)))/20.0),'FaceColor','black','EdgeColor','none');
    %     p2 = patch(isosurface(template.template, max(max(max(template.template)))/20.0),'FaceColor',[0.95,0.95,0.95],'EdgeColor','none');
    %     isonormals(template.template,p1);
    %     isonormals(template.template,p2);
    %     view(2); % Set view in 2D plane so user can see dipole orientation clearly
    %     axis tight; axis equal; camlight; camlight(-80,-10); lighting phong; title('Dipole'); % Axis stuff
end


%%%%%%%%%%%
% EOF
%%%%%%%%%%%