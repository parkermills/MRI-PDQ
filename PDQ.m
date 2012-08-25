%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @file PDQ.m
% @brief Runs entire PDQ operation on an MRIdata object
% @author Parker Mills, Ahrens Lab, Carnegie Mellon 2009
%
% ==================== INPUT PARAMETERS ======================
% @param MRIdata                  (MRIdata)        MRI dataset to be analyzed by PDQ (see README.txt for datatype details)
% @param CHI_background           (Float)          Magnetic susceptibility of background material surrounding SPIO deposits.
%                                                  Necessary for accurate determination of dipole's susceptibility and magnetic dipole moment
%
% @param orientation              (Float)          OPTIONAL Orientation axis, if already known beforehand
% @param user_radius              (Float)          OPTIONAL Estimated radius range of sphere of SPIO. 2 element vector if upper-bound and lower-bound known
%                                                           3 element vector if optimal_radius also known: [lower_bound upper_bound] -OR- [lower_bound upper_bound optimal_radius]
%                    
% @param dual_gaussian            (Float)          OPTIONAL 
%
% ==================== RETURNED DATA =========================
% @return MRIdata                 (MRIdata)            PDQ result (see README.txt)
%
% ==================== ASSUMPTIONS ===========================
% @assume If a data item is already present in the MRIdata structure, PDQ will not re-calculate those items.
%         You must delete them in order for them to be recalculated.
%         For example, "MRIdata = rmfield(MRIdata, 'PDQ');" will erase all PDQ results so they may be calculated from scratch.
%
% ==================== EXAMPLE USAGE =========================
% MRIdata = PDQ(MRIdata, CHI_background, orientation, user_radius, dual_gaussian);
% jhu_gel3 = PDQ(jhu_gel3, -9.035e-6, 1, [65 540], 0);
% lesleybrain0 = PDQ(lesleybrain0, -9.035e-6, 1, [0 60 32], 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MRIdata = PDQ(MRIdata, CHI_background, orientation, user_radius, dual_gaussian)


%% Preferences
platform = 'unix';  % Set to either 'windows' or 'unix' or 'mac'. On unix and mac platforms, automatic phase-unwrapping will be attempted.
                    %On windows, phase maps will be exported for manual unwrapping.

xcorr_cutoff = 0.3; % Dipoles with XCORR similarity values below this threshold are ignored.
                    % Typically set to 0.3 (30% similar). Must be between [0.0 1.0] (positive similarity).
                    % Setting below 0.1 is a false-positive disaster. Setting below 0.3 causes a moderatly-dense class of false-positive results.

xcorr_cutoff_multiradii = 0.6; % Dipoles with XCORR values below this threshold are not include in the dipole radii fit.
                               % Typically set to 0.6, meaning dipoles with < 60% resemblance to the template won't be fit
                               % This is because we assume that any dipole with <60% resemblance will get an inaccurate radius fit
                               
maximum_unwrap_hours = 24;      % At most we want to wait this many hours for a phase dataset to be unwrapped by PRELUDE. This may result in
                               % some phase unwrapping errors, but it keeps the maximum time down

pref_run_phase_ramp = 0;   % Run PDQ on phase-ramp-removed data. Worth doing to seeing how results differ from high-passed
                           % phase if ramp is very homogenous. But of course, takes twice as long.


% Technical Parameters
peak_detection_threshold = 0.0; % Peaks in XCORR similariry below this value will not be considered. We only desire peaks
                                % with a positive similarity, so this is set to 0.0
pref_phase_std = 0.0;           % When generating mask, threshold out phase incoherence below this threshold (default is 0.0)
pref_default_noise_std = 1.0;   % When generating mask, this is the default sigma value for eliminated noise
pref_num_of_radii = 50;  % When it comes to fitting dipole radius, how many would you like to test?
pref_template_shift = sqrt(3)/4.0;     % What fraction of a pixel do you want to shift the templates for detecting off-center dipoles?
                                       % This value is approximately 0.433 (sqrt(3)/4.0)
                                       % Choose this number to be irrational may reduce chance of artifacts.



%% Initializations

% Reboot matlab pool, if available, to allow for parallel function execution
disp('PDQ: Initializing');
matlabpool close force
matlabpool open

% Add PDQ field to MRIdata, if not present
if(~isfield(MRIdata,'PDQ'))
    MRIdata.PDQ = [];
end

% Store dimensions
[x_dim y_dim z_dim] = size(MRIdata.mag);
MRIdata.PDQ.xcorr_cutoff = xcorr_cutoff;

% Prompt user for key variables, if not given in function call
if(~exist('orientation','var'))
    if(isfield(MRIdata,'unwrapped'))
        fig1 = images(MRIdata.unwrapped); title('MRI data raw phase image');
    else
        fig1 = images(MRIdata.phase); title('MRI data raw phase image');
    end
    orientation = input('What direction is B0?   1: Up-Down   2: In-Out   3: Left-Right  [1-3]: ');
    close(fig1);
end
if(~exist('user_radius','var'))
    user_radius = input('What is estimated radius range of SPIO spheres? [microns microns]: ');
end



%% Calculate mask from magnitude image
disp('PDQ: Calculating masks');
if(~isfield(MRIdata,'mask'))
    
    % Establish whether the dataset is a single or dual-Gaussian image
    if(~exist('dual_gaussian','var'))
        dual_gaussian = 0;
        reply_dual_gaussian = input('Is this a dual-Gaussian image (e.g., sample in liquid)? y/n [n]:', 's');
        if(strcmp(reply_dual_gaussian,'y'))
            dual_gaussian = 1;
        end
    end
    
    MRIdata.mask = generate_mask(MRIdata, pref_default_noise_std, pref_phase_std, dual_gaussian,1);
end



%% Unwrap phase images
disp('PDQ: Unwrapping phase images');
if(~isfield(MRIdata,'unwrapped'))
    MRIdata = unwrap_phase_image(MRIdata, platform, maximum_unwrap_hours);
    if(strcmp(platform,'windows'))
        error('PDQ: You are using Windows. Data has been exported using prelude_export. You must manually unwrap your images using PRELUDE. Data saved to disc as MRIdata.mat. Quitting.');
    end
end



%% If phase images were unwrapped in 2D, normalize phase through 3rd dimension
if(strcmp(MRIdata.k_space_type, '2d'))
    MRIdata.unwrapped = Normalize3Dunwrap(MRIdata.unwrapped,MRIdata.mag);
    
    % Prompt user for sanity check
    images(reshape(MRIdata.unwrapped(:,ceil(y_dim/2),:),x_dim,z_dim)); title('Normalized volume'); % Show normalized volume as sanity check
    reply = input('Stack Normalization OK? Y/N [Y]: ', 's');
    if strcmp(reply,'n')
        error('PDQ: Stack normalization failed - requires some sort of user intervention. Quitting.');
    end
end



%% Calculate high-passed, Rausher, and ramp-removed phase images
disp('PDQ: Removing low-frequency phase changes from phase images');
if(~isfield(MRIdata,'high_pass')) MRIdata.high_pass = hp(MRIdata.unwrapped); end
if(pref_run_phase_ramp)
    if(~isfield(MRIdata,'ramp')) MRIdata = phase_ramp_remove(MRIdata, 1); end
end


% Temporary Checkpoint
save PDQ_checkpoint MRIdata;



%% Generate templates for large range of radii and template shifts
disp('PDQ: Generating templates spanning range of radii based on user-estimate');
approx_template = PDQ_generate_template(MRIdata.resolution, orientation, mean(user_radius), MRIdata.B0, MRIdata.TE, 1.0e-5, [0 0 0]);
if(~isfield(MRIdata,'template_spectrum'))
    MRIdata.template_spectrum = PDQ_generate_template_spectrum(approx_template, pref_num_of_radii, user_radius, pref_template_shift);
end

% Temporary Checkpoint
save PDQ_checkpoint MRIdata;



%% Run dipole search using the 'approximate template', using it to create the 'detection template'
disp('PDQ: Performing preliminary dipole search based on user-provided estimates');
if(~isfield(MRIdata.PDQ,'detection_template'))
    if(length(user_radius) == 3)
        MRIdata.PDQ.detection_template =PDQ_generate_template(MRIdata.resolution, orientation, user_radius(3), MRIdata.B0, MRIdata.TE, 7.0e-5, [0 0 0]);
    else
        x_correlate_hp = real(normxcorr3(approx_template.template, MRIdata.high_pass, 'same')) .* MRIdata.mask; % Command ordering: (template, image, shape)
        raw_peaks_hp = logical(peaks3d(x_correlate_hp, peak_detection_threshold)) .* x_correlate_hp;
        dipoles = PDQ_dipole_analysis(MRIdata.high_pass, raw_peaks_hp, xcorr_cutoff, approx_template, CHI_background);
        disp('PDQ: Generating optimized dipole detection template from preliminary search results');
        MRIdata.PDQ.detection_template = PDQ_generate_optimal_template(dipoles, approx_template, MRIdata.template_spectrum);
    end
end



%% Calculate 3D XCORR volumes and find peaks for high-passed phase and (optionally) ramp-removed phase
% High-pass (default)
if(~isfield(MRIdata.PDQ,'dipoles_hp'))
    disp('PDQ: Performing dipole search on high-passed phase using optimized dipole detection template');
    x_correlate_hp = real(normxcorr3(MRIdata.PDQ.detection_template.template, MRIdata.high_pass, 'same')) .* MRIdata.mask; % Command ordering: (template, image, shape)
    raw_peaks_hp = logical(peaks3d(x_correlate_hp, peak_detection_threshold)) .* x_correlate_hp;
    disp('PDQ: Calculating susceptibility of high-passed dipoles, assuming fixed radius');
    MRIdata.PDQ.dipoles_hp = PDQ_dipole_analysis(MRIdata.high_pass, raw_peaks_hp, xcorr_cutoff, MRIdata.PDQ.detection_template, CHI_background);
end

% Ramp
if(pref_run_phase_ramp)
    if(~isfield(MRIdata.PDQ,'dipoles_ramp'))
        disp('PDQ: Performing dipole search on ramp-removed phase using optimized dipole detection template');
        x_correlate_ramp = real(normxcorr3(MRIdata.PDQ.detection_template.template, MRIdata.ramp, 'same')) .* MRIdata.mask; % Command ordering: (template, image, shape)
        raw_peaks_ramp = logical(peaks3d(x_correlate_ramp, peak_detection_threshold)) .* x_correlate_ramp;
        disp('PDQ: Calculating susceptibility of ramp-removed dipoles, assuming fixed radius');
        MRIdata.PDQ.dipoles_ramp = PDQ_dipole_analysis(MRIdata.ramp, raw_peaks_ramp, xcorr_cutoff, MRIdata.PDQ.detection_template);
    end
end



%% Check to see if any dipoles were found. If not, return!
if(isempty(MRIdata.PDQ.dipoles_hp))
   return
end



%% Find neighbors of dipoles (within template dimensions, which are significant pixels if they overlap)
disp('PDQ: Finding neighbors for each dipole');
if(~isfield(MRIdata.PDQ,'neighbor_consensus'))
    [MRIdata.PDQ.dipoles_hp neighbor_consensus_hp] = PDQ_calculate_neighbors(MRIdata.PDQ.dipoles_hp);
    if(pref_run_phase_ramp)
        [MRIdata.PDQ.dipoles_ramp neighbor_consensus_ramp] = PDQ_calculate_neighbors(MRIdata.PDQ.dipoles_ramp);
    end
    MRIdata.PDQ.neighbor_consensus = neighbor_consensus_hp; % Consensus result comes from high-passed phase
end



%% Calculate Radii and Chi for non-fixed radii
disp('PDQ: Calculating susceptibility of dipoles, assuming variable radii');
if(~isfield(MRIdata.PDQ,'dipoles_multiradii_hp'))
    MRIdata.PDQ.dipoles_multiradii_hp = PDQ_calculate_radii_and_chi(MRIdata.PDQ.detection_template, MRIdata.PDQ.dipoles_hp, xcorr_cutoff_multiradii, MRIdata.template_spectrum, CHI_background);
end
if(pref_run_phase_ramp)
    if(~isfield(MRIdata.PDQ,'dipoles_multiradii_ramp'))
        MRIdata.PDQ.dipoles_multiradii_ramp = PDQ_calculate_radii_and_chi(MRIdata.PDQ.detection_template, MRIdata.PDQ.dipoles_ramp, xcorr_cutoff_multiradii,MRIdata.template_spectrum, CHI_background);
    end
end



%% If ROIs (regions of interest) are specified, count dipoles within each region, and determine distance of dipoles outside of regions
disp('PDQ: Calculating statistics for user-specified ROIs (regions of interest)');
if(isfield(MRIdata,'regions'))
    MRIdata = PDQ_process_regions(MRIdata, xcorr_cutoff_multiradii);
end



%% Display minimal output from PDQ Run
disp('PDQ: ###### Finished ######');
disp(['PDQ result: ',num2str(length(MRIdata.PDQ.dipoles_hp)), ' dipoles found in high-passed phase with XCORR greater than ',num2str(xcorr_cutoff)]);
disp(['PDQ result: Optimal fixed radius was found to be ',num2str(MRIdata.PDQ.detection_template.radius),' microns']);
disp(['PDQ result: ',num2str(length(MRIdata.PDQ.dipoles_multiradii_hp)), ' dipoles found in high-passed phase with XCORR greater than ',num2str(xcorr_cutoff_multiradii)]);



%% Close worker pool
matlabpool close force



%%%%%%%%%%%%%
% EOF
%%%%%%%%%%%%%
