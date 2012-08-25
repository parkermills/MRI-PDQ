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
% @param    a              (Float)          Spheroid radius (microns)
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
% template = generate_PDQ_template(resolution, orientation, a, B0, TE, d_Chi, shift);
% template = generate_PDQ_template([51 51 98],    3, 200, 11.7, 11, 1.03e-4, [0 0 0]);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function template = PDQ_generate_template_old(resolution, orientation, radius, B0, TE, d_Chi, shift)


%% User-set preferences
pref_gridsize = 7; % Size of dipole template, in pixels (e.g. default is 7, resulting in template sized 7x7x7)
                   % Usually set to N=7 or N=9, as dipole impression falls off like r^3. Should be set to
                   % smallest reasonable value so that the template doesn't factor in phase perturbations in other areas!
                   
pref_invert_dipole_template = 1; % For some reason, the entire dipole template must be inverted.

pref_sampling_multiplier = 13; % Default is >14. Dipole is oversampled by this factor before template is calculated.
                               % 20 eliminates most artifacts, 24 ideal if templates are saved for later
                               % Do NOT use an easily-divisible number. This creates artifacts. Use a prime number like 5, 7, 13, 17, 19, or 23
                               % Large values are slow, but are more accurate. Must be > 5.0 for accurate results.
                               
pref_print_magnetic_information = 0; % Print out template's magnetic specifications

pref_display_figures = 1; % Visualize template by showing figures for user

pref_model_inside = 'nothing'; % Model Bz inside the sphere of SPIO? Set to 'schenck', 'griffiths', or 'nothing'
                               % Signal is usually too low inside dipoles to model their inside.
                               % Also, this is set only if the material is homogeneous.
                               
gamma = 2.0 * pi * 42576000; % Gyromagnetic ratio (Default is for Hydrogen = 42,576,000 (2*pi*Hz)/T)



% Technical settings
pref_CHI_default = 1.0; % Just a default CHI value if none is provided by user
pref_tesla_mhz_decision = 24; % B0 assumed to be in Tesla if below this value, MHz if above this value



%% Set physical constants
mu_0 = 4e-7 * pi; % Permeability of free space (Tesla * meter / Amp)



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
    d_Chi = pref_CHI_default;
end
if(d_Chi == 0.0) % If user doesn't provide susceptibility values, issue a default
    d_Chi = pref_CHI_default;
end



%% Convert all units into meters, seconds, Tesla
resolution = resolution ./ 1e6; % Convert resolution from microns to meters
TE = TE ./ 1000.0; % Convert TE from milliseconds into seconds
radius = radius ./ 1e6; % Convert radius from microns to meters

% If B0 is in MHz, convert it to Tesla.
if(abs(B0) > pref_tesla_mhz_decision)
    B0 = B0 .* 1e6 ./ gamma;
end

% If we want to invert the dipole
if(pref_invert_dipole_template)
    if(B0 > 0.0)
        B0 = -B0;
    end
end



%% Perform magnetic calculations related to template and its surrounding environment
% Formula for mdm is: m(emu) = [(1000 * 4 * pi * a^3 * B0) / (3 * mu_0)] / [1/d_Chi + 1]
volume = (4.0 / 3.0) * pi * radius.^3.0; % Calculate volume of a sphere (m^3)
M_true     = (abs(B0) / mu_0) / ( (1/d_Chi)  + 1.0 ); % Bulk magnetization (A/m)
CHI_true     = 1.0 / ((abs(B0)/(mu_0 * M_true)) - 1.0); % CHI (unitless SI) reverse-calculated
mdm_true     = M_true * volume; % Magnetic Dipole Moment (Amp * meter^2)
mdm_true_emu     = mdm_true     * 1000.0; % MDM (emu)



%% If preferred, print out all the magnetic calculations
if(pref_print_magnetic_information)
    disp(['You have created a dipole with:']);
    disp(['True volume susceptibility (dimensionless SI units): ',num2str(d_Chi)]);
    disp(['Volume: ',num2str(volume * 1e18),' cubic microns or ',num2str(volume * 1e9),' cubic millimeters']);
    disp(['True bulk magnetization: ',num2str(M_true),' A/m']);
    disp(['True Magnetic Moment: ',num2str(mdm_true),' A*m^2 or ',num2str(mdm_true_emu),' emu']);
    disp(['If this is saturated magnetite (92 emu/g[Fe]), then you have ',num2str(mdm_true_emu .* 1e9 ./ 92.0),' ng of Fe in this sphere']);
    disp([' ']);
end



%% Create grid representing dipole in space. X, Y, Z values are in meters
% For each voxel, generate a high-resolution grid containing 'pref_sampling_multiplier' points, then average to make template
upper_bound_x = ((pref_gridsize + 2.0 * shift(1)) * resolution(1)) / 2.0;
upper_bound_y = ((pref_gridsize + 2.0 * shift(2)) * resolution(2)) / 2.0;
upper_bound_z = ((pref_gridsize + 2.0 * shift(3)) * resolution(3)) / 2.0;
lower_bound_x = ((-pref_gridsize + 2.0 * shift(1)) * resolution(1)) / 2.0;
lower_bound_y = ((-pref_gridsize + 2.0 * shift(2)) * resolution(2)) / 2.0;
lower_bound_z = ((-pref_gridsize + 2.0 * shift(3)) * resolution(3)) / 2.0;

x_supersample = resolution(1) / pref_sampling_multiplier;
y_supersample = resolution(2) / pref_sampling_multiplier;
z_supersample = resolution(3) / pref_sampling_multiplier;


switch orientation
    case 1
        [x,z,y] = meshgrid(lower_bound_y:y_supersample:upper_bound_y, lower_bound_x:x_supersample:upper_bound_x, lower_bound_z:z_supersample:upper_bound_z);
    case 2
        [x,y,z] = meshgrid(lower_bound_y:y_supersample:upper_bound_y, lower_bound_x:x_supersample:upper_bound_x, lower_bound_z:z_supersample:upper_bound_z);
    case 3
        [z,y,x] = meshgrid(lower_bound_y:y_supersample:upper_bound_y, lower_bound_x:x_supersample:upper_bound_x, lower_bound_z:z_supersample:upper_bound_z);
    otherwise
        error('generate_PDQ_template: Orientation must be value of 1, 2, or 3. Quitting.');
end



%% Form dipole-shape from spherical coordinate geometry
% Math note:
% r = sqrt(x.^2 + y.^2 + z.^2)
% cos(theta).^2 = z.^2 / r.^2

% Does not yet factor in ANYTHING about the material (e.g., dCHI, B0, 1/3, or a^3)
radius2 = x.^2 + y.^2 + z.^2; % Calculate radius squared
warning off MATLAB:divideByZero
dipole = ( 3.*((z.^2)./radius2) - 1 ) ./ ( sqrt(radius2).^3 ); % UNITS: meters
warning on MATLAB:divideByZero



%% Calculate magnetic field inside sphere (OBSOLETE)
% Math note:
% Bz = (2/3 * dChi * B_0)   (Schenck)
% Bz = mu_0 (2/3 * M)       (Griffiths)

% Set Bz inside of sphere to user preference
inside_Bz_1 = pi*gamma*TE*B0   .* (2/3) .* d_Chi;
inside_Bz_2 = pi*gamma*TE*mu_0 .* (2/3) .* M_true;
switch pref_model_inside
    case 'schenck'
        inside_Bz = B0   .* (2/3) .* d_Chi;   % Bz inside of sphere if it was a homogeneous media (Tesla)
    case 'griffiths'
        inside_Bz = mu_0 .* (2/3) .* M_true;                 % Bz inside of sphere if it was a homogeneous media (Tesla)
    case 'nothing'
        inside_Bz = 0.0;
end



%% Calculate magnetic field outside sphere
outside_Bz = dipole .* d_Chi .* (radius .^ 3.0) .* B0 ./ 3.0;


% PHASE MULTIPLIED BY 2.0 * PI
%% Calculate phase-offset template values from Bz, inside and outside of sphere
inside_phase  =  gamma .* TE .* inside_Bz; % Units: 1/(T*s) * s * T
outside_phase =  gamma .* TE .* outside_Bz; %Units: 1/(T*s) * s * T
dipole_inside  = inside_phase  .* (sqrt(x.^2 + y.^2 + z.^2) < radius);
dipole_outside = outside_phase .* (sqrt(x.^2 + y.^2 + z.^2) >=  radius);
dipole_model = dipole_outside + dipole_inside;
dipole_model(isnan(dipole_model)) = inside_phase; % Set R = 0.0 to dipole_inside



%% Shrink dipole model to desired size
template.template = ones(pref_gridsize, pref_gridsize, pref_gridsize); % Allocate template

for j_1 = 1:pref_gridsize
    for j_2 = 1:pref_gridsize
        for j_3 = 1:pref_gridsize
            
            % Map template pixels to super-sampled pixels
            x_start = 1 + (j_1 - 1) * pref_sampling_multiplier;
            x_stop = 1 + (j_1) * pref_sampling_multiplier;
            y_start = 1 + (j_2 - 1) * pref_sampling_multiplier;
            y_stop = 1 + (j_2) * pref_sampling_multiplier;
            z_start = 1 + (j_3 - 1) * pref_sampling_multiplier;
            z_stop = 1 + (j_3) * pref_sampling_multiplier;
            phase_average = sum(sum(sum(dipole_model(x_start:x_stop,y_start:y_stop,z_start:z_stop))));

            % Calculate average
            x_pixels = x_stop - x_start + 1;
            y_pixels = y_stop - y_start + 1;
            z_pixels = z_stop - z_start + 1;
            total_pixels = x_pixels * y_pixels * z_pixels;
            template.template(j_1,j_2,j_3) = phase_average / total_pixels;
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
template.inside1 = inside_Bz_1;
template.inside2 = inside_Bz_2;



%% Visualize high-resolution dipole model and central slice of template
if(pref_display_figures)
    % 2D Visualization
%     images(dipole_model);
    images(template.template(:,:,5));
    
    % 3D Visualization
%     figure;
%     p1 = patch(isosurface(dipole_model,-max(max(max(dipole_model)))/20.0),'FaceColor','black','EdgeColor','none');
%     p2 = patch(isosurface(dipole_model, max(max(max(dipole_model)))/20.0),'FaceColor',[0.95,0.95,0.95],'EdgeColor','none');
%     isonormals(dipole_model,p1);
%     isonormals(dipole_model,p2);
%     view(2); % Set view in 2D plane so user can see dipole orientation clearly
%     axis tight; axis equal; camlight; camlight(-80,-10); lighting phong; title('Dipole'); % Axis stuff
end


%%%%%%%%%%%
% EOF
%%%%%%%%%%%