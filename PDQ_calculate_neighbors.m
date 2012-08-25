


function [dipoles neighbor_consensus] = PDQ_calculate_neighbors(dipoles)

%% Initialize
num_dipoles = length(dipoles);
template_dims = size(dipoles(1).phase);
x_dim_template = template_dims(1);
y_dim_template = template_dims(2);
z_dim_template = template_dims(3);


%% Go through dipole list to find neighbors of each dipole
for i = 1:num_dipoles
    
    % Determine range as half-template size
    low_x  = dipoles(i).x - (x_dim_template-1)/2;
    high_x = dipoles(i).x + (x_dim_template-1)/2;
    
    % Go through list for comparison with other dipoles
    num_neighbors = 0;
    for k = 1:num_dipoles
        neigh_location_x = dipoles(k).x;
        
        % Time-reduction sanity check (checks x-dimension first)
        if((neigh_location_x > low_x) && (neigh_location_x < high_x))
            
            low_y  = dipoles(i).y - (y_dim_template-1)/2;
            high_y = dipoles(i).y + (y_dim_template-1)/2;
            low_z  = dipoles(i).z - (z_dim_template-1)/2;
            high_z = dipoles(i).z + (z_dim_template-1)/2;
            
            % If the alternate dipole is unique, and within a template-sized region of the first dipole
            if((i ~= k) && (dipoles(k).y > low_y) && (dipoles(k).y < high_y)  && (dipoles(k).z > low_z) && (dipoles(k).z < high_z)  )
                % Neighbor found! Put this dipole's index into the original dipoles index list
                num_neighbors = num_neighbors + 1;
                dipoles(i).neighbors(num_neighbors) = k;
            end
            dipoles(i).num_neighbors = num_neighbors;
        end
        
        
    end % All comparison dipoles
    
    % Put in blank neighbors field if has no neighbors
    if(~isfield(dipoles(i),'neighbors'))
       dipoles(i).neighbors = 0; 
    end
    
end % All dipoles

% Put in blank neighbors field if has no neighbors



%% Generate consensus
num_bins = max([dipoles.num_neighbors])+1;
sum_xcorr_per_bin = zeros(1,num_bins);
sum_dipoles_per_bin = zeros(1,num_bins);
neighbor_consensus = single(zeros(x_dim_template, y_dim_template, z_dim_template, num_bins));

for j_7 = 1:num_dipoles
    temp_num_neighbors = dipoles(j_7).num_neighbors+1; %Offset by one to avoid zero-indexing
    sum_dipoles_per_bin(temp_num_neighbors) = sum_dipoles_per_bin(temp_num_neighbors) + 1;
    neighbor_consensus(:,:,:,temp_num_neighbors) = neighbor_consensus(:,:,:,temp_num_neighbors) + dipoles(j_7).phase;
    sum_xcorr_per_bin(temp_num_neighbors) = sum_xcorr_per_bin(temp_num_neighbors) + dipoles(j_7).xcorr;
end
for j_8=1:length(sum_xcorr_per_bin)
    sum_xcorr_per_bin(j_8) = sum_xcorr_per_bin(j_8) / sum_dipoles_per_bin(j_8);
    neighbor_consensus(:,:,:,j_8) = neighbor_consensus(:,:,:,j_8) / sum_dipoles_per_bin(j_8);
end


%%%%%%%%%%%%%
%    EOF
%%%%%%%%%%%%%