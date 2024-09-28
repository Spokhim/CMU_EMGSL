
%% 1. Load file and look for a decent channel to generate masks from peaks.
TMSI_SAGA_MATLAB_SDK_FOLDER = 'C:/MyRepos/NML/MWB_TMSi';
INPUT_FILES = 'C:/Data/MetaWB/MCP04_2024_09_17/TMSi/MCP04_2024_09_17_A_PROX_3.poly5';
Y_OFFSET = 50; % Set spacing between traces

addpath(TMSI_SAGA_MATLAB_SDK_FOLDER);
x = TMSiSAGA.Poly5.read(INPUT_FILE);
[b,a] = butter(3,100/2000,'high');
uni = filtfilt(b,a,x.samples(1:64,:)')';
t = 0:(1/x.sample_rate):((size(uni,2)-1)/x.sample_rate);
figure('Color','w','WindowState','maximized'); 
plot(uni' + (0:Y_OFFSET:(Y_OFFSET*63)));
set(gca,'ColorOrder',jet(64),'YTick',0:Y_OFFSET:(Y_OFFSET*63),'YTickLabel',1:64);
[~,f,~] = fileparts(INPUT_FILE);
title(gca,strrep(f,'_','\_'));

%% 2. Find peak locations of a good channel
GOOD_CH = 25;
THRESHOLD_UV = 45;
[~,locs] = findpeaks(-uni(GOOD_CH,:),'MinPeakHeight',THRESHOLD_UV);
fig = figure('Color','w','Name',sprintf('Channel-%02d Peaks', GOOD_CH)); 
ax = axes(fig,'NextPlot','add','FontName','Tahoma');
plot(ax,uni(GOOD_CH,:), ...
    'MarkerIndices',locs, ...
    'Marker','*', ...
    'MarkerEdgeColor','r');

% Get snippets of individual waveforms
vec = -20:20;
mask = locs' + vec;

snips = cell(numel(locs),1);
for i = 1:numel(locs)
    snips{i} = uni(:,mask(i,:));
end
snips = horzcat(snips{:}); 
R = (1/(size(snips,2))) * (snips * snips');
P = pinv(R);

wsnips = zeros(64,numel(vec),numel(locs));
for i = 1:numel(locs)
    wsnips(:,:,i) = P  * uni(:,mask(i,:));
end

% Need to get component for all channels
muaps = mean(wsnips,3);

save(sprintf('Data/%s_muaps_mask.mat',f), 'mask');
save(sprintf('Data/%s_muap_templates.mat',f), 'muaps');


