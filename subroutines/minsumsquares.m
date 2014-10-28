%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @name minsumsquares.m
% @author Parker Mills, Ahrens Lab, Carnegie Mellon 2008
% @brief Calculates multiple of template using Least-Squares-Fit to data
% 
% ==================== INPUT PARAMETERS ======================
% @param   A                 (2D/3D Float)    Image
% @param   T                 (2D/3D Float)    Template
% 
% ==================== RETURNED DATA =========================
% @return  multiple    (Float)          LSF-fit multiple of Template that fits input Image
% 
% ==================== ASSUMPTIONS ======================
% @assume  User defines pref_num_iter and init_value
% 
% ==================== EXAMPLE USAGE =========================
% multiple = minsumsquares(A, T);
% multiple = minsumsquares(image, template);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function multiple = minsumsquares(A, T)

%% Preferences
pref_num_iter = 50; % Number of iterations to perform.
                    % 30 iterations does good convergence, 60 performs high-precision.
pref_init_value = 1.0; % Initial value for multiple
pref_init_value_2 = 1.0; % Initial value for constant offset



%% Initialize Simplex algorithm
[params, state] = Simplex('init', [pref_init_value, pref_init_value_2]);  % Perform Simplex Initialization



%% Run Simplex algorithm
for iters = 1:pref_num_iter
    value = mean(mean(mean(     ( ( (params(1) .* T) - params(2)) - A )   .^(2.0)    ))); %
    [params, state] = Simplex(value);   % Perform Simplex Iteration
    paramsy(iters) = params(1);
    paramsy2(iters) = params(2);
end


%% Check that dimension 


%% Report final value
multiple = Simplex('centroid');
multiple = multiple(1);

%%%%%%%%%%%%%%
%    EOF
%%%%%%%%%%%%%%