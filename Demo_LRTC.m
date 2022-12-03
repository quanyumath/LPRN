currentFolder = pwd;
addpath(genpath(currentFolder));
clear;
close all;
clc
indimgs = [1:500];

%% read data and produce mask
i = 3;
id = indimgs(i);
pic_name = ['ColorImage/', num2str(id), '.jpg'];
I = double(imread(pic_name));
Z = I / max(I(:)); % imshow(Z)
Nway = size(Z);
N = ndims(Z);

%% 观察到得元素 \Omege
sr = 0.3;
p = round(sr*prod(Nway));
Omega = randsample(prod(Nway), p);
data = Z(Omega);
H = 0 * Z;
H(Omega) = Z(Omega);  % imshow(H)
maxIter = 200;
tol = 1e-3;

%% LPRN
tic
[X_LPRN, RSE_tubal, Iter_tubal] = LRTC_LPRN(0.6, data, Omega, maxIter, tol/2, Z);
time_LPRN = toc;
[psnr_LPRN, ssim_LPRN, fsim_LPRN] = quality(Z, X_LPRN); %imshow(X_LPRN)
