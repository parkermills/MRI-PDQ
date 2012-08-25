%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @name PDQ_rerun.m
% @author Parker Mills, Ahrens Lab, Carnegie Mellon 2009
% @brief  Runs PDQ on an MRIdata with 3 features:
%         1) Minimizes memory usage (good for large datasets)
%         2) Deletes previous PDQ results, if present
%         3) Saves the MRIdata to disk
%
% @trusted yes
% @robust yes
% @commented yes
% @optimized yes
% @parallelized N/A
%
% ==================== INPUT PARAMETERS ======================
% @param MRIdataname              (String)        Name of MRIdata dataset to be analyzed by PDQ
% @param orientation              (Float)         Orientation axis, if already known beforehand
% @param a                        (Float)         Estimated radius of sphere of SPIO
%
% ==================== RETURNED DATA =========================
% @return 
%
% ==================== ASSUMPTIONS ===========================
%
%
% ==================== EXAMPLE USAGE =========================
% PDQ_rerun(MRIdataname, orientation, estimate);
% PDQ_rerun('jhu_gel3', 1, 300);
% PDQ_rerun('lesleybrain0', 2, 40);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PDQ_rerun(MRIdataname, CHI_background, orientation, estimate, dual_gaussian)

% Appent ".mat" to MRIdata file name and import
MRIdataname = strcat(MRIdataname, '.mat');
MRIdata = importdata(MRIdataname);

% Trim MRIdata of k_space (not needed for PDQ process) to save memory.
MRIdata.k_space = 1.0;

% Remove for redo: high-computational cost PDQ results
% MRIdata = rmfield_if_present(MRIdata,'unwrapped');
% MRIdata = rmfield_if_present(MRIdata,'mask');
% MRIdata = rmfield_if_present(MRIdata,'template_spectrum');

% Remove for redo: low-computational cost PDQ results
MRIdata = rmfield_if_present(MRIdata,'high_pass');
MRIdata = rmfield_if_present(MRIdata,'PDQ');

% (Re-)run PDQ
MRIdata = PDQ(MRIdata, CHI_background, orientation, estimate, dual_gaussian);

% Reload original MRIdata, restoring removed data items
MRIdata_temp = importdata(MRIdataname);
MRIdata.k_space = MRIdata_temp.k_space;

clear MRIdata_temp;

% Save PDQ results, uncompressed so we don't have to wait 10 years
save(MRIdataname, 'MRIdata', '-v6');

% Clean up
clear MRIdata
close all hidden

end


function MRIdata = rmfield_if_present(MRIdata,fieldname)
if(isfield(MRIdata,fieldname))
    MRIdata = rmfield(MRIdata,fieldname);
end
end


%%%%%%%%%%%%%%
%    EOF
%%%%%%%%%%%%%%
