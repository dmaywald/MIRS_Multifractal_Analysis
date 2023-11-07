close all; clear all; clc
%addpath('/Users/dixon/Documents/TAMU/DemosNew')
%addpath('/Users/dixon/Documents/TAMU/DistVar/supCodes/function/')

addpath('Matlab Codes/')

MIRS = readtable('MIRS_Data/MIRS_spectra.csv');
MIRS = table2array(MIRS);

K = 10;l = 2^K;

% wavelet spectra

filt= [-0.075765714789341  -0.029635527645954   0.497618667632458 ...
          0.803738751805216   0.297857795605542  -0.099219543576935 ...
        -0.012603967262261   0.032223100604071];   

ismean = 2;  isplot = 1;

data = MIRS(50,:);
N=length(data);  LN= log2(N); NN=floor(LN);
y=log10(1./data(1:2^NN));

fig = figure(1);

subplot(1,2,1)
[slope, levels, log2spec ] = waveletspectra(y, 1, filt, 3, 8,ismean, isplot);

subplot(1,2,2)
 q = -5:.1:6;
[a, b, c] = mfstriangle(y, 1, filt, q, 3, 9, 1);

% 

q = -5:.1:6;

W_features = zeros(size(MIRS,1), 12);

for i = 1: size(MIRS,1)
data = MIRS(i,:);
N=length(data);  LN= log2(N); NN=floor(LN);
y= log10(1./data(1:2^NN));

[slope, levels, log2spec ] = waveletspectra(y, 1, filt, 3, 9,ismean,  0);

W_features(i,1) = slope;

[a, b, c] = mfstriangle(y, 1, filt, q, 3, 9, 0);
% c = [H LS RS LT RT B pt_L pt_R, MC, K, KC];

W_features(i,2:end) = c;
end 

% writematrix( W_features,'Wavelet_Features.csv');
% 
% %%
% 
% MIRS_feature = readtable('MIRS_Features.csv');
% MIRS_feature = table2array(MIRS_feature);
% 
% 
% y1 = MIRS_feature(:,3);
% 
% 
% figure(2)
% [f1,xi1] = ksdensity(y1); 
% plot(xi1,f1, 'r-','linewidth', 2);
% 
% hold on
% [f2,xi2] = ksdensity(W_features(:,1)); 
% plot(xi2,f2,'b--','linewidth', 2);
% hold off
% 
% 
% plot(MIRS_feature(:,2),  W_features(:,2),'.')
% xlim([0,1]); ylim([0, 1])