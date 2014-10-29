%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @name vectorize_and_remove_zeros.m
% @author Parker Mills, Ahrens Lab, Carnegie Mellon
% @brief Turns a matrix into a vector, removing all zero values from that vector
%
% ==================== INPUT PARAMETERS ======================
% @param   A        (ND Float)   Matrix to be vectorized and zero-cleansed
%
% ==================== RETURNED DATA =========================
% @return  vector   (1D Float)   Vector with no zero values in it
%
% ==================== ASSUMPTIONS ======================
% @assume
%
% ==================== DISPLAYED PRODUCTS ====================
% @product
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out_vector = vectorize_and_remove_zeros(A, type, xcorr_cutoff, mdm_range_cutoff, chi_range_cutoff)

%% Initialize

% Get information on Matrix A
switch ndims(A)
    case 1
        return
    case 2
        [x_dim y_dim] = size(A);
        v_length = x_dim * y_dim;
    case 3
        [x_dim y_dim z_dim] = size(A);
        v_length = x_dim * y_dim * z_dim;
    case 4
        [x_dim y_dim z_dim e_dim] = size(A);
        v_length = x_dim * y_dim * z_dim * e_dim;
    otherwise
        error('vectorize_and_remove_zeros: Provided matrix must be 1, 2, 3, or 4D. Quitting.');
end

% If chi range isn't specified, just set it to be all values
if(~exist('chi_range_cutoff','var'))
    chi_range_cutoff = [-inf inf];
end

% Check dimensionality
if(ndims(A) > 4)
    error('vectorize_and_remove_zeros: Provided matrix must be 1, 2, 3, or 4D. Quitting.');
end

% Store dimension data
[x_dim y_dim z_dim e_dim] = size(A);

% Create default mdm & chi range cutoffs if they are not specified
if(~exist('xcorr_cutoff','var'))
    xcorr_cutoff = 0.3;
end
if(~exist('chi_range_cutoff','var'))
    chi_range_cutoff = [-inf inf];
end
if(~exist('mdm_range_cutoff','var'))
    mdm_range_cutoff = [-inf inf];
end



%% Reshape matrix into vector
A = A(:);
vector_length = 0;

% If a dipole vector, add to output vector only if dipole's xcorr value is above xcorr_cutoff and mdm is within mdm_range_cutoff
if(exist('type','var'))
    if(strcmp(type,'dipoles'))
        for j = 1:v_length
            if(  (A(j).xcorr > xcorr_cutoff)  &&...
                    (A(j).m_core > mdm_range_cutoff(1)) &&...
                    (A(j).m_core < mdm_range_cutoff(2))  &&...
                    (strcmp(A(j).verity,'invalid'))...
                    )
                vector_length = vector_length + 1;
                out_vector(vector_length) = A(j);
            end
        end
    end
    % Otherwise, if just a plain old vector
else
    for j = 1:length(A)
        if(A(j) ~= 0)
            vector_length = vector_length + 1;
            out_vector(vector_length) = A(j);
        end
    end
end


% If there's no out_vector (e.g., no dipoles above the xcorr_cutoff value)
if(~exist('out_vector','var'))
    out_vector = [];
end


%%%%%%%%%%%%%%
%   EOF
%%%%%%%%%%%%%%