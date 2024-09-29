%%EXAMPLE_POLY5_TO_MUAP_PIPELINE This script illustrates how to convert TMSi Poly5 EMG recordings to MUAPs for use with the Jupyter notebook pipeline.
close all force;
clear;
clc;

%% 1. Load file and look for a decent channel to generate masks from peaks.
DATA_INPUT_ROOT = "C:/Data/MetaWB";
SUBJ = "MCP05";
YYYY = 2024;
MM = 8;
DD = 29;
BLOCK = 3;
TRIGGER_BIT = 0;
PLOT_ALL_CHANNELS = true;
PLOT_THRESHOLDED_CHANNEL = true;
MUAP_RELATIVE_PEAK_SAMPLES = -20:20;
CMAP = [cm.umap([166 25 46]./255,32); ...
        cm.umap([51 63 72]./255,32); ...
        cm.umap([240 90 40]./255,32); ...
        cm.umap([164 97 164]./255,32);...
        cm.umap([255 184 28]./255,32); ...
        cm.umap([0 178 169]./255,32);...
        cm.umap([124 194 66]./255,32); ...
        cm.umap([108 172 228]./255,32);];
CMAP = double(CMAP)./255;
REMAP = [1:128, 193:256, 129:192];
TANK = sprintf("%s_%04d_%02d_%02d", SUBJ, YYYY, MM, DD);
INPUT_FILE = sprintf("%s/%s/MotorUnits Decomposition/Decomposition Input/%s_%d_synchronized.mat",DATA_INPUT_ROOT, TANK, TANK, BLOCK);
Y_OFFSET = 50; % Set spacing between traces
load(INPUT_FILE,'uni','sync');
uni = uni(REMAP,:);
uni = apply_sd_textiles(uni);
t = 0:(1/2000):((size(uni,2)-1)/2000);

%% 2. (Optional) Plot all synchronized channels
if PLOT_ALL_CHANNELS
    fig = figure('Color','w','WindowState','maximized'); 
    ax = axes(fig,'NextPlot','add','FontName','Tahoma', ...
        'ColorOrder',CMAP, ...
        'YTick',0:(4*Y_OFFSET):(Y_OFFSET*255), ...
        'YTickLabel',1:4:64, ...
        'YLim',[-Y_OFFSET, 256*Y_OFFSET]);
    plot(ax, uni' + (0:Y_OFFSET:(Y_OFFSET*255)));
    yline(ax,(63.5:64:191.5)*Y_OFFSET,'k--');
    f = sprintf('%s_%d', TANK, BLOCK);
    title(ax, strrep(f,'_','\_'), 'FontName','Tahoma','Color','k');
end

%% 3. Find peak locations of a good channel
GOOD_CH = 69;
THRESHOLD_UV = [10,20];
[pks,locs] = findpeaks(-uni(GOOD_CH,:),'MinPeakHeight',THRESHOLD_UV(1));
i_remove = pks > THRESHOLD_UV(2);
pks(i_remove) = [];
locs(i_remove) = [];

%% 4. (Optional) Plot the peaks we just thresholded
if PLOT_THRESHOLDED_CHANNEL
    fig = figure('Color','w','Name',sprintf('Channel-%02d Peaks', GOOD_CH)); 
    ax = axes(fig,'NextPlot','add','FontName','Tahoma');
    plot(ax,uni(GOOD_CH,:), ...
        'MarkerIndices',locs, ...
        'Marker','*', ...
        'MarkerEdgeColor','r');
end

%% 5. Extract/save MUAP indexing mask and MUAP whitened waveforms.
% Get snippets of individual waveforms
mask = locs' + MUAP_RELATIVE_PEAK_SAMPLES;
mask(any((mask<1) | (mask > size(uni,2)),2),:) = [];
gesture_mask = bitand(sync,2^TRIGGER_BIT) == 0;
n_snips = size(mask,1);

snips = cell(n_snips,1);
for i = 1:n_snips
    snips{i} = uni(:,mask(i,:));
end
% signal = horzcat(snips{:}); 
signal = uni(:,gesture_mask);
data_cov = (1/(size(signal,2))) * (signal * signal');
noise = uni(:,~gesture_mask);
noise_cov = (1/(size(noise,2))) * (noise * noise');
snips = cat(3,snips{:});
muaps = mean(snips,3);

% P = pinv(data_cov);
% wsnips = cat(3,snips{:});
% for i = 1:n_snips
%     wsnips(:,:,i) = P * wsnips(:,:,i);
% end
% muaps = mean(wsnips,3);

save(sprintf('Data/%s_muaps_mask.mat',f), 'mask');
save(sprintf('Data/%s_muaps_templates.mat',f), 'muaps');
save(sprintf('Data/%s_muaps_covariance.mat',f), 'data_cov', 'noise_cov');
fprintf(1,'\n<strong>COMPLETE</strong>: MUAP Masks/Templates/Covariances saved: Data/%s_muaps_*.mat\n\n', f);


