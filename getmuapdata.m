
%% 1. Load file and look for a decent channel to generate masks from peaks.
TMSI_SAGA_MATLAB_SDK_FOLDER = 'C:/MyRepos/NML/MWB_TMSi';
INPUT_FILE = 'C:/Data/MetaWB/MCP04_2024_09_17/TMSi/MCP04_2024_09_17_A_PROX_3.poly5';
Y_OFFSET = 50; % Set spacing between traces

addpath(TMSI_SAGA_MATLAB_SDK_FOLDER);
x = TMSiSAGA.Poly5.read(INPUT_FILE);
[b,a] = butter(3,100/2000,'high');
uni = filtfilt(b,a,x.samples(1:64,:)')';
t = 0:(1/x.sample_rate):((size(uni,2)-1)/x.sample_rate);
figure('Color','w','WindowState','maximized'); 
plot(t,uni' + (0:Y_OFFSET:(Y_OFFSET*63)));
set(gca,'ColorOrder',jet(64),'YTick',0:Y_OFFSET:(Y_OFFSET*63),'YTickLabel',1:64);
[~,f,~] = fileparts(INPUT_FILE);
title(gca,strrep(f,'_','\_'));

%% 2. Find peak locations of a good channel
% GOOD_CH = 1;
GOOD_CH = 25;
[~,locs] = findpeaks(-uni(GOOD_CH,:),'MinPeakHeight',100);
figure; 
plot(uni(GOOD_CH,:),'MarkerIndices',locs,'Marker','*','MarkerEdgeColor','r');
vec = -20:20;
mask = locs' + vec;
tmp = uni(1,:);
snips = tmp(mask);
figure; plot(mean(snips,1));
title(sprintf('Ch-%02d Templates', GOOD_CH));

% Need to get component for all channels
muaps = zeros(64,41);
for i = 1:64
    tmp = uni(i,:);
    snips = tmp(mask);
    % figure; plot(mean(snips,1));
    muaps(i,:) = mean(snips,1);
end

save(sprintf('Data/%s_muaps_mask.mat',f), 'mask');


