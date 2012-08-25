%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @name PDQ_inspect.m
% @author Parker Mills, Ahrens Lab, Carnegie Mellon 2010
% @brief Allows user to manually verify dipoles in tissue volume.
%
% ==================== INPUT PARAMETERS ======================
% @param   XXXXXXXXXXXXXXX
%
% ==================== RETURNED DATA =========================
% @return  MRIdata with dipole list 'verity' field with following values:
%                  'unverified': Dipole unverified (Default after running PDQ)
%                  'probable'  : Dipole automatically determined to be probable dipole
%                  'genuine'   : Dipole manually verified as genuine
%                  'invalid'   : Dipole automatically or manually determined to be invalid
%
% ====================== ASSUMPTIONS =========================
% @assume
%
% ==================== DISPLAYED PRODUCTS ====================
%
% ==================== EXAMPLE USAGE =========================
% MRIdata = PDQ_inspect(XXXXXXXXXXXXXXXXXXXXXXX)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function MRIdata = PDQ_inspect(MRIdata, xcorr_cutoff, mdm_range_cutoff, inspect_every_dipole, full_screen)

%% Hard-coded preferences
border_validity = 2; % if dipole within this many pixels from border, it is declared invalid



%% Initializations
data_size = size(MRIdata.high_pass);
dipole_list_length = length(MRIdata.PDQ.dipoles_hp);
multiradii_dipole_list_length = length(MRIdata.PDQ.dipoles_multiradii_hp);



%% Mark all dipoles as 'unverified'
for i = 1:dipole_list_length
    MRIdata.PDQ.dipoles_hp(i).verity = 'unverified';
end



%% Go through dipole list, performing automatic disqualifications
for i = 1:dipole_list_length
    current_dipole = MRIdata.PDQ.dipoles_hp(i);
    
    % Disqualify dipoles on borders of anatomy
    within_tissue = 1;
    
    % Mark as not within tissue if close to border
    for j=current_dipole.x - border_validity:current_dipole.x + border_validity
        for k=current_dipole.y - border_validity:current_dipole.y + border_validity
            for l=current_dipole.z - border_validity:current_dipole.z + border_validity
                if(~MRIdata.mask(min(abs(j)+1,data_size(1)),min(abs(k)+1,data_size(2)),min(abs(l)+1,data_size(3))))
                    within_tissue = 0;
                end
            end
        end
    end
    
    % Mark dipole as invalid if not within tissue
    if (~within_tissue)
        MRIdata.PDQ.dipoles_hp(i).verity = 'invalid';
        
        % Also mark dipole in the high-quality dipole list
        for j=1:multiradii_dipole_list_length
            if(MRIdata.PDQ.dipoles_multiradii_hp(j).index == i)
                MRIdata.PDQ.dipoles_multiradii_hp(j).verity = 'invalid';
            end
        end
        
        % Otherwise mark dipole as probable
    else
        MRIdata.PDQ.dipoles_hp(i).verity = 'probable';
        for j=1:multiradii_dipole_list_length
            if(MRIdata.PDQ.dipoles_multiradii_hp(j).index == i)
                MRIdata.PDQ.dipoles_multiradii_hp(j).verity = 'probable';
            end
        end
    end
    
    
    %% Disqualify dipoles that neighbor another dipole (discard the one of the pair with lower XCORR value)
    % Gather XCORR values for all neighbors
    neighboring_xcorrs = zeros(current_dipole.num_neighbors+1,1);
    neighboring_xcorrs(1) = current_dipole.xcorr;
    if(current_dipole.num_neighbors)
        for j = 2:current_dipole.num_neighbors+1
            neighboring_xcorrs(j) = MRIdata.PDQ.dipoles_hp(current_dipole.neighbors(j-1)).xcorr;
        end
    end
    
    % Determine which neighbor has the highest XCORR value
    [max_xcorr max_xcorr_index] = max(neighboring_xcorrs);
    
    % If the original dipole is not the largest, label it as invalid -- the largest one will later be marked valid
    if (max_xcorr_index ~= 1)
        MRIdata.PDQ.dipoles_hp(i).verity = 'invalid';
    end
    
    
    
    %% Disqualify dipoles that have too small XCORR values
    if (current_dipole.xcorr < xcorr_cutoff)
        MRIdata.PDQ.dipoles_hp(i).verity = 'invalid';
    end
    
    
    
    %% Disqualify dipoles that have MDM values that are out-of-bounds
    if(  (current_dipole.m_core < mdm_range_cutoff(1))   ||   (current_dipole.m_core > mdm_range_cutoff(2))  )
        MRIdata.PDQ.dipoles_hp(i).verity = 'invalid';
    end
    
    
end




%% Go through dipole list and ask user to verify dipoles that haven't been automatically eliminated
if (inspect_every_dipole)
    for i = 1:dipole_list_length
        
        current_dipole = MRIdata.PDQ.dipoles_hp(i);
        
        % If  current dipole hasn't been labeled invalid, present it to user for manual inspection
        if(~strcmp(MRIdata.PDQ.dipoles_hp(i).verity, 'invalid'))
            
            % Show slice for current dipole
            inspection_figure = PDQ_visualize_slice(MRIdata, current_dipole.z, 1, xcorr_cutoff, mdm_range_cutoff, full_screen);
            
            % Prompt user to verify current dipole's validity
            valid = input(['Is dipole #',num2str(i),' valid? (y/n) [Y]: '], 's');
            if isempty(valid)
                valid = 'y';
            end
            
            % Store result
            if(strcmp(valid,'y'))
                MRIdata.PDQ.dipoles_hp(i).verity = 'genuine';
            else
                MRIdata.PDQ.dipoles_hp(i).verity = 'invalid';
            end
            
            close(inspection_figure);
            
        end
    end
    
end



%%%%%%%%%%%
%   EOF
%%%%%%%%%%%