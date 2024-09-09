x = TMSiSAGA.Poly5.read('/Users/pokhims/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/Documents/Coding/CMU_EMGSL/Data/Pok_2024_08_21_A_PROX_8.poly5');
[b,a] = butter(3,100/2000,'high');
uni = filtfilt(b,a,x.samples(1:64,:)')';
figure; plot(uni(1,:));
figure; plot(uni(1,:));
[~,locs] = findpeaks(-uni(1,:),'MinPeakHeight',100);
figure; plot(uni(1,:),'MarkerIndices',locs,'Marker','*','MarkerEdgeColor','r');
vec = -20:20;
mask = locs' + vec;
tmp = uni(1,:);
snips = tmp(mask);
figure; plot(mean(snips,1));

% Need to get component for all channels
muaps = zeros(64,41);
for i = 1:64
    tmp = uni(i,:);
    snips = tmp(mask);
    figure; plot(mean(snips,1));
    muaps(i,:) = mean(snips,1);
end



