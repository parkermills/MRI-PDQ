%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @name gaussian_fit.m
% @author Parker Mills, Ahrens Lab, Carnegie Mellon 2009
% @brief Creates a 1, or 2 gaussian fit
%
% ==================== INPUT PARAMETERS ======================
% @param data                         (1D/2D/3D Float)     Data to ungergo Gaussian fit
% @param num_bins                     (Float)              Number of bins
% @param pref_display_gaussian_fits   (Float)              2: Display lots of info 1: Display Gaussian fits, 0: Don't display anything
%
% ==================== RETURNED DATA =========================
% @return mask     (3D Float)     Mask
%
% ====================== ASSUMPTIONS =========================
% @assume Nothing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ true_mean true_std        true_mean_1 true_mean_2 true_std_1 true_std_2    ] = gaussian_fit(data, num_bins, pref_display_gaussian_fits)

%% Preferences
pref_pixels_per_bin = 5000; % Pixels per bin, if user doesn't provide num_bins
pref_minimum_bins = 100; % Minimum number of bins, if user doesn't provide num_bins



%% Initialize
% Store dimensions, min&max values
[x_dim y_dim z_dim] = size(data);
min_data = min(min(min(data)));
max_data = max(max(max(data)));
pixels_data = x_dim * y_dim * z_dim;
if(~exist('num_bins','var'))
    num_bins = max(pixels_data / pref_pixels_per_bin, pref_minimum_bins);
end

% Parameters are [a b c d], where a = amplitude, b = mean, c = stdev, d = scalar offset
% Parameters are [a1 a2 b1 b2 c1 c2 d], where a = amplitude, b = mean, c = stdev, d = scalar offset
% Set lower bounds. All parameters are contrained to be > 0.0
lower_bound_1 = [0.0, 0.0, 0.0, 0.0];
lower_bound_2 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];

% Set upper bounds
% Amplitude must be < pixels_data, mean must be < num_bins, stdev must be < num_bins, scalar offset must be < pixels_data
upper_bound_1 = [pixels_data,               num_bins,               num_bins,             pixels_data];
upper_bound_2 = [pixels_data,  pixels_data, num_bins,  num_bins,    num_bins,  num_bins,  pixels_data];

% Set initial values
initial_values_1 = [pixels_data                 (num_bins/2.0)                     (num_bins/4.0)                         1.0];
initial_values_2 = [pixels_data   pixels_data   (num_bins/2.0)    (num_bins/2.0)   (num_bins/4.0)    (num_bins/4.0)       1.0];



%% Use MATLAB histogram tool to summarize intensity distribution
data_hist = hist( reshape(data, 1 ,x_dim * y_dim * z_dim), num_bins).';



%% Process the data_hist in the curve fitting tool

% Set fit options - this fit assumes a gaussian signal distribution
fit_options_1 = fitoptions('Method','NonlinearLeastSquares','Lower',lower_bound_1,'Upper',upper_bound_1,'Startpoint',initial_values_1);
fit_options_2 = fitoptions('Method','NonlinearLeastSquares','Lower',lower_bound_2,'Upper',upper_bound_2,'Startpoint',initial_values_2);
fit_type_1 = fittype('   a  * (1 / (c * sqrt(2 * pi)))  * exp(-((x-b)^2)/(2*(c^2)))                                                                            +  d','options', fit_options_1);
fit_type_2 = fittype('   a1 * (1 / (c1 * sqrt(2 * pi))) * exp(-((x-b1)^2)/(2*(c1^2)))     +     a2 * (1 / (c2 * sqrt(2 * pi))) * exp(-((x-b2)^2)/(2*(c2^2)))     +  d','options', fit_options_2);
x_data = [1:num_bins].';

% Run fits
[fit_summary_1, fit_quality_1, fit_results_1] = fit(x_data, data_hist, fit_type_1);
[fit_summary_2, fit_quality_2, fit_results_2] = fit(x_data, data_hist, fit_type_2);

% Display fit information
%disp(fit_summary_1);
fit_summary_values_1 = coeffvalues(fit_summary_1);
a = fit_summary_values_1(1);
b = fit_summary_values_1(2);
c = fit_summary_values_1(3);
d = fit_summary_values_1(4);

%disp(fit_summary_2);
fit_summary_values_2 = coeffvalues(fit_summary_2);
a1 = fit_summary_values_2(1);
a2 = fit_summary_values_2(2);
b1 = fit_summary_values_2(3);
b2 = fit_summary_values_2(4);
c1 = fit_summary_values_2(5);
c2 = fit_summary_values_2(6);
dd = fit_summary_values_2(7);

range = max_data - min_data;
bin_size = range / num_bins;

% One fit
true_mean = min_data + b * bin_size;
true_std = c * bin_size;

% Two fits
true_mean_1 = min_data + b1 * bin_size;
true_mean_2 = min_data + b2 * bin_size;
true_std_1 = c1 * bin_size;
true_std_2 = c2 * bin_size;
%disp(['Actual values are: Mean: ',num2str(true_mean),' Std: ',num2str(true_std)]);



%% A little bit of Gaussian Fundaments
% PDF = 1/(sigma sqrt(2pi)) exp (- (x-u)^2 / (2 sigma^2) )
% mean = u
% mode = u
% var = sigma^2



%% Plot histogram, fit, and threshold for user sanity-check

if(pref_display_gaussian_fits)
    figure;
    bar(data_hist);
    hold on;
    plot(a .* (1.0 / (c .* sqrt(2.0 * pi))) .* exp(-((x_data-b).^2)/(2.0 .* (c.^2))) + d,'r','LineWidth',2);
    hold off;
    
    % Display additional information (dual-gaussian fits and actual results of fit)
    if(pref_display_gaussian_fits == 2)
        disp(['Single fit: mean +/- std = ',num2str(true_mean),' +/- ',num2str(true_std)]);
        disp(' ');
        figure;
        bar(data_hist);
        hold on;
        plot(   a1 .* (1.0 / (c1 .* sqrt(2.0 * pi))) .* exp(-((x_data-b1).^2)/(2.0 .* (c1.^2)))     +       a2 .* (1.0 / (c2 .* sqrt(2.0 * pi))) .* exp(-((x_data-b2).^2)/(2.0 .* (c2.^2)))       + dd,'r','LineWidth',2);
        hold off;
        disp(['Dual fit: mean1 +/- std1 = ',num2str(true_mean_1),' +/- ',num2str(true_std_1)]);
        disp(['Dual fit: mean2 +/- std2 = ',num2str(true_mean_2),' +/- ',num2str(true_std_2)]);
        disp(' ');
    end

end



%%%%%%%%%%%%%%%%%%
%      EOF
%%%%%%%%%%%%%%%%%%