%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @name generate_optimal_template.m
% @author Parker Mills, Ahrens Lab, Carnegie Mellon 2009
% @brief Generates an optimal PDQ template based on dipole list and the template used to generate the dipole list
% @trusted yes
% @robust yes
% @commented yes
% @optimized yes
% @parallelized yes
% 
% ==================== INPUT PARAMETERS ======================
% @param dipoles            (1D dipole)        List of dipoles
% @param template           (template)         Template used to form dipole list
%
% ==================== RETURNED DATA =========================
% @return optimal_template   (template)         The optimal template
%
% ==================== ASSUMPTIONS ===========================
% @assume 
%
% ==================== EXAMPLE USAGE =========================
% optimal_template = generate_optimal_template(dipoles, template)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function optimal_template = PDQ_generate_optimal_template(dipoles, template, template_spectrum)


%% Preferences

pref_min_dipoles = 10; % Minimum number of dipoles allowed as candidates for the optimal template
pref_max_dipoles = 100; % Maximum number of dipoles allowed as candidates for the optimal template

pref_xcorr_increment = 0.025; % Value by which xcorr_cutoff is incremented/decremented in order to
                             % reach a number of dipoles between [pref_min_dipoles pref_max_dipoles]

pref_radii_bottom_threshold = 0.02; % Percent of lowest radii in template_spectrum that are not considered
                                    %

              
%% Initialize

% Set a starting cutoff value. Value doesn't matter since it is soon optimized
xcorr_cutoff = 0.7;

% If there are more than pref_max_dipoles dipoles at this xcorr_cutoff value, make cutoff more strict
while( sum(  [dipoles.xcorr] > xcorr_cutoff  ) > pref_max_dipoles && (xcorr_cutoff < 1.0))
   xcorr_cutoff = xcorr_cutoff + pref_xcorr_increment;
end

% If there are fewer than pref_min_dipoles at this xcorr_cutoff value, loosen things up!
while( sum(  [dipoles.xcorr] > xcorr_cutoff  ) < pref_min_dipoles )
   xcorr_cutoff = xcorr_cutoff - pref_xcorr_increment;
   if(xcorr_cutoff < 0.3)
       warning(['generate_optimal_template: Need to find more than ',num2str(pref_min_dipoles), ' dipoles. Try a different estimate. Returning same template.']);
       optimal_template = template;
       return;
   end
end



%% Find dipoles with XCORR greater than the established cutoff.
num_dipoles = length(dipoles);
great_dipole_count = sum([dipoles.xcorr] > xcorr_cutoff);
great_dipoles = zeros(great_dipole_count, 1);
great_dipole_xcorr = zeros(great_dipole_count, great_dipole_count);
great_dipole_index = 0;
for j1 = 1:num_dipoles
    if(dipoles(j1).xcorr > xcorr_cutoff)
        great_dipole_index = great_dipole_index + 1;
        great_dipoles(great_dipole_index) =  j1;
    end
end



%% Xcorr each 'great dipole' with every other 'great dipole'. Store results.
for j1 = 1:great_dipole_count
    for j2 = j1:great_dipole_count
        great_dipole_xcorr(j1,j2) = normxcorr3(dipoles(great_dipoles(j1)).phase, dipoles(great_dipoles(j2)).phase,'valid');
    end
end



%% Find top 3 pairs and combine each pair to create "top 3 template types"

% Sort pairs
[mm,im] = max(great_dipole_xcorr);
[ss,is] = sort(mm);

% First, Second, and Third place pairs
firstpair = dipoles(great_dipoles(im(length(ss)  ))).phase + dipoles(great_dipoles(is(length(ss)  ))).phase;
secondpair = dipoles(great_dipoles(im(length(ss)-1))).phase + dipoles(great_dipoles(is(length(ss)-1))).phase;
thirdpair = dipoles(great_dipoles(im(length(ss)-2))).phase + dipoles(great_dipoles(is(length(ss)-2))).phase;

% Normalize each dipole pair, since they were added together
firstpair = firstpair ./ 2.0;
secondpair = secondpair ./ 2.0;
thirdpair = thirdpair ./ 2.0;
combined = (firstpair + secondpair + thirdpair) ./ 3.0;
images(combined); title('Top three dipole pairs, combined into one template');



%% Fit the top 3 pairs using the dipole generation engine
value_spectrum = [];
index_start = ceil(pref_radii_bottom_threshold * length(template_spectrum(2,2,2,:)));
parfor radii = index_start:length(template_spectrum(2,2,2,:))
    value = normxcorr3(combined, template_spectrum(2,2,2,radii).template, 'valid');
    value_spectrum = [value_spectrum value];
end

[optimal_radius_xcorr, optimal_index] = max(value_spectrum);
optimal_radius = template_spectrum(2,2,2, optimal_index + index_start - 1).radius;



%% Create optimal template with desired radius
optimal_template = PDQ_generate_template(template.resolution, template.orientation, optimal_radius, template.B0, template.TE, 1.0e-5, [0 0 0]);


end

%%%%%%%%%%%%%%%%
%    EOF
%%%%%%%%%%%%%%%%