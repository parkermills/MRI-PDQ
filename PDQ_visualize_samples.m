
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parker Mills, Ahrens Lab 2009
% Loads data for, and visualizes, specific experiment data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CMU Gels
% % Load files and reduce memory burden by clearing out certain fields
% cmu_gel2 = importdata('cmu_gel2.mat');
% cmu_gel2.k_space = 1.0; cmu_gel2.complex = 1.0;
% cmu_gel2.mag = 1.0; cmu_gel2.phase = 1.0; cmu_gel2.phase_hamming = 1.0;
% cmu_gel2.unwrapped= 1.0; cmu_gel2.unwrapped_hamming= 1.0;
% cmu_gel2.high_pass= 1.0; cmu_gel2.rauscher= 1.0;
% 
% cmu_gel5 = importdata('cmu_gel5.mat');
% cmu_gel5.k_space = 1.0; cmu_gel5.complex = 1.0;
% cmu_gel5.mag = 1.0; cmu_gel5.phase = 1.0; cmu_gel5.phase_hamming = 1.0;
% cmu_gel5.unwrapped= 1.0; cmu_gel5.unwrapped_hamming= 1.0;
% cmu_gel5.high_pass= 1.0; cmu_gel5.rauscher= 1.0;

% PDQ_visualize([cmu_gel2 cmu_gel5]);



%% JHU Gels
% jhu_gel0 = importdata('jhu_gel0.mat'); % Different voxel size

jhu_gel1_high = importdata('jhu_gel1_high.mat');
jhu_gel2= importdata('jhu_gel2.mat');
jhu_gel3 = importdata('jhu_gel3.mat');
jhu_gel4 = importdata('jhu_gel4.mat');

jhu_gel5 = importdata('jhu_gel5.mat');
jhu_gel6 = importdata('jhu_gel6.mat');
jhu_gel7 = importdata('jhu_gel7.mat');
jhu_gel8 = importdata('jhu_gel8.mat');
jhu_gel9 = importdata('jhu_gel9.mat');

jhu_gel10 = importdata('jhu_gel10.mat');
jhu_gel11 = importdata('jhu_gel11.mat'); % Different voxel size
jhu_gel12 = importdata('jhu_gel12.mat');

% jhu_gel1_high = importdata('jhu_gel1_high.mat');
% jhu_gel1_medium = importdata('jhu_gel1_medium.mat');
% jhu_gel1_low = importdata('jhu_gel1_low.mat');
% jhu_gel1_preview = importdata('jhu_gel1_preview.mat');


% % All JHU Gels
jhu_10 = [ jhu_gel1_high jhu_gel2 jhu_gel3 jhu_gel4];
jhu_5 = [ jhu_gel5 jhu_gel6 jhu_gel7 jhu_gel8 jhu_gel9];
jhu_25 = [ jhu_gel10 jhu_gel11 jhu_gel12];
all = [jhu_10 jhu_5 jhu_25];


% JHU gels of different noise levels
%PDQ_visualize(jhu_gel1_high);
%PDQ_visualize(jhu_gel1_medium);
%PDQ_visualize(jhu_gel1_low);
%PDQ_visualize(jhu_gel1_preview);



%% Old Hearts
%oldheart1 = importdata('oldheart1.mat');
%oldheart3 = importdata('oldheart3.mat');



%% Lesley brains
% lesleybrain0 = importdata('lesleybrain0.mat');
% lesleybrain1 = importdata('lesleybrain1.mat');

lb3=PDQ_inspect(importdata('lb3.mat'),0.31,[0.0 inf],0,1);%[.19e-12 2.6e-12]
lb4=PDQ_inspect(importdata('lb4.mat'),0.31,[0.0 inf],0,1);
lb5=PDQ_inspect(importdata('lb5.mat'),0.31,[0.0 inf],0,1);
lb6=PDQ_inspect(importdata('lb6.mat'),0.31,[0.0 inf],0,1);
lb7=PDQ_inspect(importdata('lb7.mat'),0.31,[0.0 inf],0,1);

PDQ_visualize(lb3, 0.31, 0.6, [-inf inf], [-inf inf], 1);
PDQ_visualize(lb4, 0.31, 0.6, [-inf inf], [-inf inf], 1);
PDQ_visualize(lb5, 0.31, 0.6, [-inf inf], [-inf inf], 1);
PDQ_visualize(lb6, 0.31, 0.6, [-inf inf], [-inf inf], 1);
PDQ_visualize(lb7, 0.31, 0.6, [-inf inf], [-inf inf], 1);

save lb3 lb3 -v6
save lb4 lb4 -v6
save lb5 lb5 -v6
save lb6 lb6 -v6
save lb7 lb7 -v6

all = [lb3 lb4 lb5 lb6 lb7];
clear lb3 lb4 lb5 lb6 lb7
for i=1:length(all)
all(i).k_space = [];
all(i).complex = [];
all(i).mag = [];
all(i).phase_hamming = [];
all(i).phase = [];
all(i).unwrapped_hamming = [];
all(i).unwrapped = [];
all(i).template_spectrum = [];
all(i).rauscher = [];
end

lb8=PDQ_inspect(importdata('lb8.mat'),0.41,[0.0 inf],0,1);
lb9=PDQ_inspect(importdata('lb9.mat'),0.41,[0.0 inf],0,1);
lb10=PDQ_inspect(importdata('lb10.mat'),0.41,[0.0 inf],0,1);
lb11=PDQ_inspect(importdata('lb11.mat'),0.41,[0.0 inf],0,1);
lb12=PDQ_inspect(importdata('lb12.mat'),0.41,[0.0 inf],0,1);

PDQ_visualize(lb8, 0.31, 0.6, [-inf inf], [-inf inf], 1);
PDQ_visualize(lb9, 0.31, 0.6, [-inf inf], [-inf inf], 1);
PDQ_visualize(lb10, 0.31, 0.6, [-inf inf], [-inf inf], 1);
PDQ_visualize(lb11, 0.31, 0.6, [-inf inf], [-inf inf], 1);
PDQ_visualize(lb12, 0.31, 0.6, [-inf inf], [-inf inf], 1);

save lb8 lb8 -v6
save lb9 lb9 -v6
save lb10 lb10 -v6
save lb11 lb11 -v6
save lb12 lb12 -v6

all = [all lb8 lb9 lb10 lb11 lb12];
clear lb8 lb9 lb10 lb11 lb12
for i=1:length(all)
all(i).k_space = [];
all(i).complex = [];
all(i).mag = [];
all(i).phase_hamming = [];
all(i).phase = [];
all(i).unwrapped_hamming = [];
all(i).unwrapped = [];
all(i).template_spectrum = [];
all(i).rauscher = [];
end

lb13=PDQ_inspect(importdata('lb13.mat'),0.41,[0.0 inf],0,1);
lb14=PDQ_inspect(importdata('lb14.mat'),0.41,[0.0 inf],0,1);
lb15=PDQ_inspect(importdata('lb15.mat'),0.41,[0.0 inf],0,1);
lb16=PDQ_inspect(importdata('lb16.mat'),0.41,[0.0 inf],0,1);

PDQ_visualize(lb13, 0.31, 0.6, [-inf inf], [-inf inf], 1);
PDQ_visualize(lb14, 0.31, 0.6, [-inf inf], [-inf inf], 1);
PDQ_visualize(lb15, 0.31, 0.6, [-inf inf], [-inf inf], 1);
PDQ_visualize(lb16, 0.31, 0.6, [-inf inf], [-inf inf], 1);

save lb13 lb13 -v6
save lb14 lb14 -v6
save lb15 lb15 -v6
save lb16 lb16 -v6

all = [all lb13 lb14 lb15 lb16];
clear lb13 lb14 lb15 lb16
for i=1:length(all)
all(i).k_space = [];
all(i).complex = [];
all(i).mag = [];
all(i).phase_hamming = [];
all(i).phase = [];
all(i).unwrapped_hamming = [];
all(i).unwrapped = [];
all(i).template_spectrum = [];
all(i).rauscher = [];
end

lb17=PDQ_inspect(importdata('lb17.mat'),0.31,[0.0 inf],0,1);
lb18=PDQ_inspect(importdata('lb18.mat'),0.31,[0.0 inf],0,1);
lb19=PDQ_inspect(importdata('lb19.mat'),0.31,[0.0 inf],0,1);
lb20=PDQ_inspect(importdata('lb20.mat'),0.31,[0.0 inf],0,1);
lb21=PDQ_inspect(importdata('lb21.mat'),0.31,[0.0 inf],0,1);

PDQ_visualize(lb17, 0.31, 0.6, [-inf inf], [-inf inf], 1);
PDQ_visualize(lb18, 0.31, 0.6, [-inf inf], [-inf inf], 1);
PDQ_visualize(lb19, 0.31, 0.6, [-inf inf], [-inf inf], 1);
PDQ_visualize(lb20, 0.31, 0.6, [-inf inf], [-inf inf], 1);
PDQ_visualize(lb21, 0.31, 0.6, [-inf inf], [-inf inf], 1);

save lb17 lb17 -v6
save lb18 lb18 -v6
save lb19 lb19 -v6
save lb20 lb20 -v6
save lb21 lb21 -v6

all = [all lb17 lb18 lb19 lb20 lb21];
clear lb17 lb18 lb19 lb20 lb21
for i=1:length(all)
all(i).k_space = [];
all(i).complex = [];
all(i).mag = [];
all(i).phase_hamming = [];
all(i).phase = [];
all(i).unwrapped_hamming = [];
all(i).unwrapped = [];
all(i).template_spectrum = [];
all(i).rauscher = [];
end


lb22=PDQ_inspect(importdata('lb22.mat'),0.31,[0.0 inf],0,1);
lb23=PDQ_inspect(importdata('lb23.mat'),0.31,[0.0 inf],0,1);
lb24=PDQ_inspect(importdata('lb24.mat'),0.31,[0.0 inf],0,1);
lb25=PDQ_inspect(importdata('lb25.mat'),0.31,[0.0 inf],0,1);
lb26=PDQ_inspect(importdata('lb26.mat'),0.31,[0.0 inf],0,1);
lb27=PDQ_inspect(importdata('lb27.mat'),0.31,[0.0 inf],0,1);
lb28=PDQ_inspect(importdata('lb28.mat'),0.31,[0.0 inf],0,1);
lb29=PDQ_inspect(importdata('lb29.mat'),0.31,[0.0 inf],0,1);

PDQ_visualize(lb22, 0.31, 0.6, [-inf inf], [-inf inf], 1);
PDQ_visualize(lb23, 0.31, 0.6, [-inf inf], [-inf inf], 1);
PDQ_visualize(lb24, 0.31, 0.6, [-inf inf], [-inf inf], 1);
PDQ_visualize(lb25, 0.31, 0.6, [-inf inf], [-inf inf], 1);
PDQ_visualize(lb26, 0.31, 0.6, [-inf inf], [-inf inf], 1);
PDQ_visualize(lb27, 0.31, 0.6, [-inf inf], [-inf inf], 1);
PDQ_visualize(lb28, 0.31, 0.6, [-inf inf], [-inf inf], 1);
PDQ_visualize(lb29, 0.31, 0.6, [-inf inf], [-inf inf], 1);

save lb22 lb22 -v6
save lb23 lb23 -v6
save lb24 lb24 -v6
save lb25 lb25 -v6
save lb26 lb26 -v6
save lb27 lb27 -v6
save lb28 lb28 -v6
save lb29 lb29 -v6

all = [all lb22 lb23 lb24 lb25 lb26 lb27 lb28 lb29];
clear lb22 lb23 lb24 lb25 lb26 lb27 lb28 lb29
for i=1:length(all)
all(i).k_space = [];
all(i).complex = [];
all(i).mag = [];
all(i).phase_hamming = [];
all(i).phase = [];
all(i).unwrapped_hamming = [];
all(i).unwrapped = [];
all(i).template_spectrum = [];
all(i).rauscher = [];
end

lb_all = all;
save lb_all lb_all
clear



%% In Vivo lesley brains
lb1_invivo_naive=PDQ_inspect(importdata('lb1_invivo_naive.mat'),0.31,[0.0 inf],0,1);
lb2_invivo_naive_scan1=PDQ_inspect(importdata('lb2_invivo_naive_scan1.mat'),0.31,[0.0 inf],0,1);
lb2_invivo_naive_scan2=PDQ_inspect(importdata('lb2_invivo_naive_scan2.mat'),0.31,[0.0 inf],0,1);
lb3_invivo_72h=PDQ_inspect(importdata('lb3_invivo_72h.mat'),0.31,[0.0 inf],0,1);
lb4_invivo_96h=PDQ_inspect(importdata('lb4_invivo_96h.mat'),0.31,[0.0 inf],0,1);
lb5_invivo_day6=PDQ_inspect(importdata('lb5_invivo_day6.mat'),0.31,[0.0 inf],0,1);
lb6_invivo_naive=PDQ_inspect(importdata('lb6_invivo_naive.mat'),0.31,[0.0 inf],0,1);
lb7_invivo_naive=PDQ_inspect(importdata('lb7_invivo_naive.mat'),0.31,[0.0 inf],0,1);
lb7_invivo_96h=PDQ_inspect(importdata('lb7_invivo_96h.mat'),0.31,[0.0 inf],0,1);
all = [lb1_invivo_naive lb2_invivo_naive_scan1 lb2_invivo_naive_scan2 lb3_invivo_72h lb4_invivo_96h lb5_invivo_day6 lb6_invivo_naive lb7_invivo_naive lb7_invivo_96h];
clear lb1_invivo_naive lb2_invivo_naive_scan1 lb2_invivo_naive_scan2 lb3_invivo_72h lb4_invivo_96h lb5_invivo_day6 lb6_invivo_naive lb7_invivo_naive lb7_invivo_96h


for i=1:length(all)
all(i).k_space = [];
all(i).complex = [];
all(i).mag = [];
all(i).phase_hamming = [];
all(i).phase = [];
all(i).unwrapped_hamming = [];
all(i).unwrapped = [];
all(i).template_spectrum = [];
all(i).rauscher = [];
end


lb8_invivo_naive=PDQ_inspect(importdata('lb8_invivo_naive.mat'),0.31,[0.0 inf],0,1);
lb8_invivo_48h=PDQ_inspect(importdata('lb8_invivo_48h.mat'),0.31,[0.0 inf],0,1);
lb8_invivo_72h=PDQ_inspect(importdata('lb8_invivo_72h.mat'),0.31,[0.0 inf],0,1);
lb10_invivo_naive=PDQ_inspect(importdata('lb10_invivo_naive.mat'),0.31,[0.0 inf],0,1);
lb11_invivo_naive=PDQ_inspect(importdata('lb11_invivo_naive.mat'),0.31,[0.0 inf],0,1);
lb11_invivo_48h=PDQ_inspect(importdata('lb11_invivo_48h.mat'),0.31,[0.0 inf],0,1);
lb11_invivo_72h=PDQ_inspect(importdata('lb11_invivo_72h.mat'),0.31,[0.0 inf],0,1);
lb11_invivo_96h=PDQ_inspect(importdata('lb11_invivo_96h.mat'),0.31,[0.0 inf],0,1);
lb12_invivo_naive=PDQ_inspect(importdata('lb12_invivo_naive.mat'),0.31,[0.0 inf],0,1);
lb12_invivo_naive_d6=PDQ_inspect(importdata('lb12_invivo_naive_d6.mat'),0.31,[0.0 inf],0,1);

all = [all lb8_invivo_naive lb8_invivo_48h lb8_invivo_72h lb10_invivo_naive lb11_invivo_naive lb11_invivo_48h lb11_invivo_72h lb11_invivo_96h lb12_invivo_naive lb12_invivo_naive_d6];
clear lb8_invivo_naive lb8_invivo_48h lb8_invivo_72h lb10_invivo_naive lb11_invivo_naive lb11_invivo_48h lb11_invivo_72h lb11_invivo_96h lb12_invivo_naive lb12_invivo_naive_d6


for i=1:length(all)
all(i).k_space = [];
all(i).complex = [];
all(i).mag = [];
all(i).phase_hamming = [];
all(i).phase = [];
all(i).unwrapped_hamming = [];
all(i).unwrapped = [];
all(i).template_spectrum = [];
all(i).rauscher = [];
end

lb_invivo_all = all;
save lb_invivo_all lb_invivo_all
clear



%% Chronic rejection hearts
p09 = importdata('p09.mat');
p10 = importdata('p10.mat');
p11 = importdata('p11.mat');
p12 = importdata('p12.mat');
p15 = importdata('p15.mat');
p18 = importdata('p18.mat');
p20 = importdata('p20.mat');
all = [p09 p10 p11 p12 p15 p18 p20];
clear p09 p10 p11 p12 p15 p18 p20;

for i=1:length(all)
    all(i).k_space = [];
    all(i).complex = [];
    all(i).mag = [];
    all(i).phase_hamming = [];
    all(i).phase = [];
    all(i).unwrapped_hamming = [];
    all(i).unwrapped = [];
    all(i).template_spectrum = [];
    all(i).rauscher = [];
end

p24 = importdata('p24.mat');
p27 = importdata('p27.mat');
p29 = importdata('p29.mat');
p30 = importdata('p30.mat');
p31 = importdata('p31.mat');
p33 = importdata('p33.mat');
all = [all p24 p27 p29 p30 p31 p33];
clear p24 p27 p29 p30 p31 p33;
for i=1:length(all)
    all(i).k_space = [];
    all(i).complex = [];
    all(i).mag = [];
    all(i).phase_hamming = [];
    all(i).phase = [];
    all(i).unwrapped_hamming = [];
    all(i).unwrapped = [];
    all(i).template_spectrum = [];
    all(i).rauscher = [];
end


p39 = importdata('p39.mat');
p41 = importdata('p41.mat');
p42 = importdata('p42.mat');
p44 = importdata('p44.mat');
p47 = importdata('p47.mat');
p48 = importdata('p48.mat');
p54 = importdata('p54.mat');
all = [all p39 p41 p42 p44 p47 p48 p54];
clear p39 p41 p42 p44 p47 p48 p54
for i=1:length(all)
    all(i).k_space = [];
    all(i).complex = [];
    all(i).mag = [];
    all(i).phase_hamming = [];
    all(i).phase = [];
    all(i).unwrapped_hamming = [];
    all(i).unwrapped = [];
    all(i).template_spectrum = [];
    all(i).rauscher = [];
end


p58 = importdata('p58.mat');
p59 = importdata('p59.mat');
p63 = importdata('p63.mat');
p09_rotated = importdata('p09_rotated.mat');
p10_rotated = importdata('p10_rotated.mat');
p20_rotated = importdata('p20_rotated.mat');
p58_rotated = importdata('p58_rotated.mat');
all = [all p58 p59 p63 p09_rotated p10_rotated p20_rotated p58_rotated];
clear p58 p59 p63 p09_rotated p10_rotated p20_rotated p58_rotated
for i=1:length(all)
    all(i).k_space = [];
    all(i).complex = [];
    all(i).mag = [];
    all(i).phase_hamming = [];
    all(i).phase = [];
    all(i).unwrapped_hamming = [];
    all(i).unwrapped = [];
    all(i).template_spectrum = [];
    all(i).rauscher = [];
end

%%%%%%%%%%%%%%
%    EOF
%%%%%%%%%%%%%%