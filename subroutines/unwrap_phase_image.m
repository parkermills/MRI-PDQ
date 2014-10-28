%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @file unwrap_phase_image.m
% @author Parker Mills, Ahrens Lab, Carnegie Mellon
% @brief Unwraps the phase image contained in an MRIdata file
%
% ==================== INPUT PARAMETERS ======================
% @param    MRIdata                 (MRIdata)    MRI data file containing phase image
% @param    platform                (String)     'unix' or 'windows' or 'mac'
% @param    maximum_unwrap_hours    (Float)      Time limit for phase unwrapping.
%                                                If unwrapping takes longer, quality will be reduced.
%
% ==================== RETURNED DATA =========================
% @return MRIdata                   (MRIdata)
%
% ==================== ASSUMPTIONS ======================
% @assume Input data is single/double floating point in range from [-Pi to Pi]
%
% ==================== EXAMPLE USAGE ========================
% MRIdata = unwrap_phase_image(MRIdata,platform,maximum_unwrap_hours)
% MRIdata = 
% See also PDQ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MRIdata = unwrap_phase_image(MRIdata, platform, maximum_unwrap_hours)

%% Preferences
% Starting value of number of regions for PRELUDE phase unwrapping
pref_prelude_regions = 2;  % Number of unwrapping regions for PRELUDE. If set too small, you will get poor results.
                           % 20 = Great if you have sample w/ distinct phase wraps, but would like really robust unwrapping.
                           %      If phase wraps are dense, this will take days, even if dataset is small.
                           %      It is better to use N=15 and suffer a few phase wraps.
                           % 15 = (Default) Sometimes required if dataset has many phase wraps or is noisy. Can take up to 24 hours.
                           % 10 = Produces few, if any, unwrapping errors in homogeneous images
                           % 4 = Minimum needed to unwrap phase wraps
                           % 2 = Sometimes works when on high-SNR samples w/ distinct phase wrap lines. If you can "unwrap by eye", try this.

                           
                           
%% Export image data for PRELUDE
save_nii(make_nii(single(MRIdata.mag)),'mag');
save_nii(make_nii(single(MRIdata.phase)),'phase');
save_nii(make_nii(single(MRIdata.mask)),'mask');



%% If on UNIX, run PRELUDE directly
if(strcmp(platform,'unix') || strcmp(platform,'mac'))
    
    % Remove any old files, if they exist
    if(exist('prelude.nii.gz','file'))
        system('rm prelude.nii.gz');
    end
    if(exist('~/prelude.nii.gz','file'))
        system('rm ~/prelude.nii.gz');
    end
    if(exist('prelude.nii','file'))
        system('rm prelude.nii');
    end
    if(exist('~/prelude.nii','file'))
        system('rm ~/prelude.nii');
    end
    
    % Wait for prelude result
    % Starting condition: we've been waiting more than we want for these images
    been_waiting_for = 60 * 60 * maximum_unwrap_hours;
    pref_prelude_regions = pref_prelude_regions + 2; % Correct for how we immediately run Prelude after subtracting a region
    while(~exist('prelude.nii.gz','file'))
        pause(2);
        been_waiting_for = been_waiting_for + 2;
        if(been_waiting_for > 60 * 60 * maximum_unwrap_hours)
            
            % Kill prelude if it's still running
            if(strcmp(platform,'unix'))
                system('pkill prelude');
            end
            if(strcmp(platform,'mac'))
                system(['for X in `ps acx | grep -i prelude | awk {' char(39) 'print $1' char(39) '}`; do kill $X; done']);
            end
                
            % If a Mac, set PATH and FSLOUTPUTTYPE environment variables
            if(strcmp(platform,'mac'))
                setenv('PATH', [getenv('PATH') ':/Applications/fsl/bin']);
                setenv('FSLOUTPUTTYPE','NIFTI_GZ');
            end
            
            pref_prelude_regions = pref_prelude_regions - 2; % Reduce quality of 3D unwrapping by 2 regions
            
            % Run prelude on data exported by save_nii function
            % 3D Case
            if(strcmp(MRIdata.k_space_type,'3d'))
                system(['prelude -n ', num2str(pref_prelude_regions),' -p phase.hdr -a mag.hdr -m mask.hdr -f -o prelude&']);
                % 2D Case
            else
                system(['prelude -n ', num2str(pref_prelude_regions),' -p phase.hdr -a mag.hdr -m mask.hdr -o prelude&']);
            end
            
            been_waiting_for = 0;
        end
    end
    
    
    % Unzip the files, load them in
    disp(['unwrap_phase_image: Phase image unwrapped using ',num2str(pref_prelude_regions),' regions.']);
    pause(60); % Wait a minute for PRELUDE to compress and write to disk
    system('cp prelude.nii.gz ~/');
    system('gunzip ~/prelude.nii.gz');
    pause(30); % Wait a bit for gunzip to complete unzipping the file
    prelude_result         = load_nii('~/prelude.nii');
    MRIdata.unwrapped = single(prelude_result.img);
    
    % Delete files
    system('rm prelude.nii.gz');
    system('rm mag.img mag.hdr mag.mat phase.img phase.hdr phase.mat mask.img mask.hdr mask.mat');
    
    
    %% If Windows was selected
else
    error('unwrap_phase_image: Windows not currently supported for 3D unwrapping');
end


%%%%%%%%%%%%%%%%%%
%     EOF
%%%%%%%%%%%%%%%%%%