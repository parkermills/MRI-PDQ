%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs PDQ on all datasets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Uses command 'PDQ_rerun' with arguments (dataname, B0 field (dipole) orientation, estimated dipole radius)



%% JHU Magnetocapsule Gels

% % 20% Feridex - Measured to have r=254 +/- 10um
% PDQ_rerun('jhu_gel0',      -9.035e-6,1,[220 540 254],0); % Note: jhu_gel0 has non-standard FOV
% 
% % 10% Feridex - Measured to have r=298 +/- 19um
% PDQ_rerun('jhu_gel1_high', -9.035e-6,1,[220 420 298],0);
% PDQ_rerun('jhu_gel2',      -9.035e-6,1,[220 420 298],0);
% PDQ_rerun('jhu_gel3',      -9.035e-6,1,[220 420 298],0);
% PDQ_rerun('jhu_gel4',      -9.035e-6,1,[220 420 298],0);
% 
% % 5% Feridex - Measured to have r=271 +/- 32um
%  PDQ_rerun('jhu_gel5',    -9.035e-6,1,[125 325 271],0);
%  PDQ_rerun('jhu_gel6',    -9.035e-6,1,[125 325 271],0);
%  PDQ_rerun('jhu_gel7',    -9.035e-6,1,[125 325 271],0);
%  PDQ_rerun('jhu_gel8',    -9.035e-6,1,[125 325 271],0);
%  PDQ_rerun('jhu_gel9',    -9.035e-6,1,[125 325 271],0);
% 
% % 2.5% Feridex - Measured to have r=286 +/- 11um
%  PDQ_rerun('jhu_gel10',    -9.035e-6,1,[65 280 286],0);
%  PDQ_rerun('jhu_gel11',    -9.035e-6,1,[65 280 286],0); % Note: jhu_gel11 has non-standard FOV
%  PDQ_rerun('jhu_gel12',    -9.035e-6,1,[65 280 286],0);
% 
% % 10% Feridex; Different SNR levels - Measured to have r=298 +/- 19um
% PDQ_rerun('jhu_gel1_high',    -9.035e-6,1,[220 420 298],0);
%  PDQ_rerun('jhu_gel1_medium',  -9.035e-6,1,[220 420 298],0);
%  PDQ_rerun('jhu_gel1_low',     -9.035e-6,1,[220 420 298],0);
%  PDQ_rerun('jhu_gel1_preview', -9.035e-6,1,[220 420 298],0);


%% Porcine Liver
%PDQ_rerun('porcine_liver_exvivo', -8.835e-6,1,[200 400 326],0);


%% CMU SPIO Gels
%  PDQ_rerun('cmu_gel2',-9.035e-6,1,[0 60 5],0); % 1.6um & 10um (used to measure 10um dipole moment)
%  PDQ_rerun('cmu_gel3',-9.035e-6,1,[0 60 5],0); % 10um
% PDQ_rerun('cmu_gel4',-9.035e-6,1,[0 60 5],0); % 1.6um
% PDQ_rerun('cmu_gel5',-9.035e-6,1,[0 60 5],0); % Control (no SPIO)



%% Lesley TBI Brains
% % % 24 Hours
% PDQ_rerun('lb3', -9.035e-6,2,[0 60 7.5],1);
%  PDQ_rerun('lb4', -9.035e-6,2,[0 60 7.5],1);
%  PDQ_rerun('lb5', -9.035e-6,2,[0 60 7.5],1);
%  PDQ_rerun('lb6', -9.035e-6,2,[0 60 7.5],1);
%  PDQ_rerun('lb7', -9.035e-6,2,[0 60 7.5],1);
% 
% % %48 Hours
%  PDQ_rerun('lb8', -9.035e-6,2,[0 60 7.5],1);
%  PDQ_rerun('lb9', -9.035e-6,2,[0 60 7.5],1);
%  PDQ_rerun('lb10',-9.035e-6,2,[0 60 7.5],1);
%  PDQ_rerun('lb11',-9.035e-6,2,[0 60 7.5],1);
%  PDQ_rerun('lb12',-9.035e-6,2,[0 60 7.5],1);
% 
% % 72 Hours
% PDQ_rerun('lb13',-9.035e-6,2,[0 60 7.5],1);
% PDQ_rerun('lb14',-9.035e-6,2,[0 60 7.5],1);
% PDQ_rerun('lb15',-9.035e-6,2,[0 60 7.5],1);
% PDQ_rerun('lb16',-9.035e-6,2,[0 60 7.5],1);
% 
% % 6 Days
% PDQ_rerun('lb17',-9.035e-6,2,[0 60 7.5],1);
% PDQ_rerun('lb18',-9.035e-6,2,[0 60 7.5],1);
% PDQ_rerun('lb19',-9.035e-6,2,[0 60 7.5],1);
% PDQ_rerun('lb20',-9.035e-6,2,[0 60 7.5],1);
% PDQ_rerun('lb21',-9.035e-6,2,[0 60 7.5],1);
% 
% % Naive
% PDQ_rerun('lb22',-9.035e-6,2,[0 60 7.5],1);
% PDQ_rerun('lb23',-9.035e-6,2,[0 60 7.5],1);
% PDQ_rerun('lb24',-9.035e-6,2,[0 60 7.5],1);
% PDQ_rerun('lb25',-9.035e-6,2,[0 60 7.5],1);
% PDQ_rerun('lb26',-9.035e-6,2,[0 60 7.5],1);
% PDQ_rerun('lb27',-9.035e-6,2,[0 60 7.5],1);
% PDQ_rerun('lb28',-9.035e-6,2,[0 60 7.5],1);
% PDQ_rerun('lb29',-9.035e-6,2,[0 60 7.5],1);



%% In Vivo / Ex Vivo Lesley Brains

%  PDQ_rerun('lb1_invivo_naive',      -9.035e-6,1,[0 60 7.5],1);
% % PDQ_rerun('lb1_exvivo_naive',      -9.035e-6,2,[0 60 7.5],1);
% % PDQ_rerun('lb1_exvivo_naive_10ms', -9.035e-6,2,[0 60 7.5],1);
% % 
%  PDQ_rerun('lb2_invivo_naive_scan1',  -9.035e-6,1,[0 60 7.5],1);
%  PDQ_rerun('lb2_invivo_naive_scan2',  -9.035e-6,1,[0 60 7.5],1);
% % PDQ_rerun('lb2_exvivo_naive_scan1',  -9.035e-6,2,[0 60 7.5],0);
% % PDQ_rerun('lb2_exvivo_naive_scan2',  -9.035e-6,1,[0 60 7.5],0);
% % 
%  PDQ_rerun('lb3_invivo_72h',        -9.035e-6,3,[0 60 7.5],1);
% % 
%  PDQ_rerun('lb4_invivo_96h',  -9.035e-6,1,[0 60 7.5],1);
% % 
%  PDQ_rerun('lb5_invivo_day6',        -9.035e-6,1,[0 60 7.5],1);
% PDQ_rerun('lb5_exvivo_naive',        -9.035e-6,1,[0 60 7.5],1);
% 
%  PDQ_rerun('lb6_invivo_naive',        -9.035e-6,1,[0 60 7.5],1);
% PDQ_rerun('lb6_exvivo_naive',     -9.035e-6,1,[0 60 7.5],0);
% 
%  PDQ_rerun('lb7_invivo_naive',        -9.035e-6,1,[0 60 7.5],1);
%  PDQ_rerun('lb7_invivo_96h',     -9.035e-6,1,[0 60 7.5],1);
% PDQ_rerun('lb7_exvivo_96h',    -9.035e-6,1,[0 60 7.5],0);
% 
%  PDQ_rerun('lb8_invivo_naive',   -9.035e-6,1,[0 60 7.5],1);
%  PDQ_rerun('lb8_invivo_48h',     -9.035e-6,1,[0 60 7.5],1);
%  PDQ_rerun('lb8_invivo_72h',     -9.035e-6,1,[0 60 7.5],0);
% PDQ_rerun('lb8_exvivo_72h',    -9.035e-6,1,[0 60 7.5],0);

%  PDQ_rerun('lb10_invivo_naive',    -9.035e-6,1,[0 60 7.5],1);
% PDQ_rerun('lb10_exvivo_naive',  -9.035e-6,1,[0 60 7.5],0);

%  PDQ_rerun('lb11_invivo_naive',    -9.035e-6,1,[0 60 7.5],1);
%  PDQ_rerun('lb11_invivo_48h',   -9.035e-6,1,[0 60 7.5],1);
%  PDQ_rerun('lb11_invivo_72h',    -9.035e-6,1,[0 60 7.5],1);
%  PDQ_rerun('lb11_invivo_96h',    -9.035e-6,1,[0 60 7.5],1);
% PDQ_rerun('lb11_exvivo_96h',    -9.035e-6,1,[0 60 7.5],1);

%  PDQ_rerun('lb12_invivo_naive',  -9.035e-6,1,[0 60 7.5],1);
%  PDQ_rerun('lb12_invivo_naive_d6',  -9.035e-6,1,[0 60 7.5],1);
 
%% New Coil TBI brains
%  PDQ_rerun('lb13_invivo_48h',  -9.035e-6, 1, [0 60 7.5], 1);
%  PDQ_rerun('lb13_invivo_72h',  -9.035e-6, 1, [0 60 7.5], 1);
%  PDQ_rerun('lb13_invivo_96h',  -9.035e-6, 1, [0 60 7.5], 1);
%PDQ_rerun('lb14_invivo_48h',  -9.035e-6, 1, [0 60 7.5], 1);
%PDQ_rerun('lb14_invivo_72h',  -9.035e-6, 1, [0 60 7.5], 1);

%PDQ_rerun('lb14_invivo_96h',  -9.035e-6, 1, [0 60 7.5], 1);
%PDQ_rerun('lb14_exvivo_96h',  -9.035e-6, 1, [0 60 7.5], 1);
%PDQ_rerun('lb14_exvivo_96h_highres',  -9.035e-6, 1, [0 60 7.5], 1);
%PDQ_rerun('lb14_exvivo_96h_ultrares',  -9.035e-6, 1, [0 60 7.5], 1);



%% Microglia
% PDQ_rerun('microglia', -9.035e-6,1,[0 60 7.5],1);



%% YiJen Chronic Rejection Hearts

% PDQ_rerun('p09',        -9.035e-6,2,[0 60 9.5],0);
% PDQ_rerun('p09_rotated',-9.035e-6,1,[0 60 9.5],0);
% PDQ_rerun('p10',        -9.035e-6,1,[0 60 9.5],0);
% PDQ_rerun('p10_rotated',-9.035e-6,2,[0 60 9.5],0);
% PDQ_rerun('p11',        -9.035e-6,1,[0 60 9.5],0);
% PDQ_rerun('p12',        -9.035e-6,1,[0 60 9.5],0);
% PDQ_rerun('p15',        -9.035e-6,1,[0 60 9.5],0);
% PDQ_rerun('p18',        -9.035e-6,1,[0 60 9.5],0);
% PDQ_rerun('p20',        -9.035e-6,1,[0 60 9.5],0);
% PDQ_rerun('p20_rotated',-9.035e-6,2,[0 60 9.5],0);
% PDQ_rerun('p24',        -9.035e-6,3,[0 60 9.5],0);
% PDQ_rerun('p27',        -9.035e-6,1,[0 60 9.5],0);
% PDQ_rerun('p29',        -9.035e-6,1,[0 60 9.5],0);
% PDQ_rerun('p30',        -9.035e-6,1,[0 60 9.5],0);
% PDQ_rerun('p31',        -9.035e-6,1,[0 60 9.5],0);
% PDQ_rerun('p33',        -9.035e-6,1,[0 60 9.5],0);
% PDQ_rerun('p39',        -9.035e-6,1,[0 60 9.5],0);
% PDQ_rerun('p41',        -9.035e-6,2,[0 60 9.5],0);
% PDQ_rerun('p42',        -9.035e-6,1,[0 60 9.5],0);
% PDQ_rerun('p44',        -9.035e-6,1,[0 60 9.5],0);
% PDQ_rerun('p47',        -9.035e-6,1,[0 60 9.5],0);
% PDQ_rerun('p48',        -9.035e-6,1,[0 60 9.5],0);
% PDQ_rerun('p54',        -9.035e-6,2,[0 60 9.5],0);
% PDQ_rerun('p58',        -9.035e-6,1,[0 60 9.5],0);
% PDQ_rerun('p58_rotated',-9.035e-6,2,[0 60 9.5],0);
% PDQ_rerun('p59',        -9.035e-6,2,[0 60 9.5],0);
% PDQ_rerun('p63',        -9.035e-6,2,[0 60 9.5],0);
PDQ_rerun('Q09_1201_4.Ww1.scan10',        -9.035e-6,3,[0 60 9.5],1);

% Two hearts from Danielle
%PDQ_rerun('de_102209b',   -9.035e-6,2,[0 60 9.5],0);
%PDQ_rerun('de_102209b_2', -9.035e-6,2,[0 60 9.5],0);



%% Obsolete datasets
% PDQ_rerun('heart0',-9.035e-6,1[0 60],0);
% PDQ_rerun('oldheart1',-9.035e-6,1,[0 60],0); % Allograft
% PDQ_rerun('oldheart3',-9.035e-6,1,[0 60],0); % Isograft
% PDQ_rerun('lesleybrain0',-9.035e-6,2,[0 60]); % More macrophages
% PDQ_rerun('lesleybrain1',-9.035e-6,2,[0 60]); % Fewer macrophages


%%%%%%%%%%%%%
%   EOF
%%%%%%%%%%%%%
