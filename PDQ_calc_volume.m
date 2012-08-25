

%% Calculates volume of MRIdata from resolution and mask data fields

function volume = PDQ_calc_volume(MRIdata, suppress_output)

volume_cubic_microns = MRIdata.resolution(1) * MRIdata.resolution(2) * MRIdata.resolution(3) * sum(sum(sum(MRIdata.mask)));
volume_cc = volume_cubic_microns * 1e-12;

if(~exist('suppress_output','var'))
    disp(['calc_volume: Dataset mask has volume of ',num2str(volume_cc),' cubic centimeters.']);
end

volume = volume_cc;

%%%%%%%%%%%%%
%    EOF
%%%%%%%%%%%%%