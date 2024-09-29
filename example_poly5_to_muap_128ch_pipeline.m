%%EXAMPLE_POLY5_TO_MUAP_PIPELINE This script illustrates how to convert TMSi Poly5 EMG recordings to MUAPs for use with the Jupyter notebook pipeline.
close all force;
clear;
clc;

%% 1. Load file and look for a decent channel to generate masks from peaks.
TMSI_SAGA_MATLAB_SDK_FOLDER = 'C:/MyRepos/NML/MWB_TMSi';
DATA_INPUT_ROOT = "C:/Data/MetaWB";
SUBJ = "MCP04";
YYYY = 2024;
MM = 9;
DD = 17;
BLOCK = 1;

A_TAG = "A_PROX";
B_TAG = "B_DIST";

PLOT_ALL_CHANNELS = true;
PLOT_THRESHOLDED_CHANNEL = true;
MUAP_RELATIVE_PEAK_SAMPLES = -20:20;

TANK = sprintf("%s_%04d_%02d_%02d", SUBJ, YYYY, MM, DD);
INPUT_FILES = [...
    sprintf("%s/%s/TMSi/%s_%s_%d.poly5",DATA_INPUT_ROOT, TANK, TANK, A_TAG, BLOCK); ...
    sprintf("%s/%s/TMSi/%s_%s_%d.poly5",DATA_INPUT_ROOT, TANK, TANK, B_TAG, BLOCK) ...
    ];
Y_OFFSET = 50; % Set spacing between traces
TRIGGER_BIT = 1;

addpath(TMSI_SAGA_MATLAB_SDK_FOLDER);
x = io.load_align_saga_data_many(INPUT_FILES,...
    "ApplySpatialFilter", true, ...
    "SpatialFilterMode", "SD Columns", ...
    "TriggerBitMask", 2);
nCh = numel(x.channels)/2;
TRIGGER_CHANNEL = nCh-2;
uni = x.samples([1:64;(nCh+1):(nCh+64)],:);
uni(rms(uni,2) < 1,:) = randn(nnz(rms(uni,2)<1),size(uni,2)).*5;
t = 0:(1/x.sample_rate):((size(uni,2)-1)/x.sample_rate);

%% 2. (Optional) Plot all synchronized channels
if PLOT_ALL_CHANNELS
    fig = figure('Color','w','WindowState','maximized'); 
    ax = axes(fig,'NextPlot','add','FontName','Tahoma', ...
        'ColorOrder',[winter(64);summer(64)], ...
        'YTick',0:(4*Y_OFFSET):(Y_OFFSET*127), ...
        'YTickLabel',1:4:64, ...
        'YLim',[-Y_OFFSET, 128*Y_OFFSET]);
    plot(ax, uni' + (0:Y_OFFSET:(Y_OFFSET*127)));
    yline(ax,63.5*Y_OFFSET,'k--');
    f = sprintf('%s_%d', TANK, BLOCK);
    title(ax, strrep(f,'_','\_'), 'FontName','Tahoma','Color','k');
end

%% 3. Find peak locations of a good channel
GOOD_CH = 49;
THRESHOLD_UV = [40,90];
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
gesture_mask = bitand(x.samples(TRIGGER_CHANNEL,:),2^TRIGGER_BIT) == 0;
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


