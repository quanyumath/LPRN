currentFolder = pwd;
addpath(genpath(currentFolder));
clear;
close all;
clc
load Fountain % imshow(X(:,:,40))
X = L(:, :, 21:70);
S0 = E(:, :, 21:70);
%X=L;S0=E;
mO = X;
[n1, n2, n3] = size(mO);

%% LPRN
% 生成时空张量
[T1SQ, I1, T2SQ, I2, T3SQ, I3] = ST(X);
tic
opts = [];
opts.I1 = I1;
opts.I2 = I2;
opts.I3 = I3;
opts.T1SQ = T1SQ;
opts.T2SQ = T2SQ;
opts.T3SQ = T3SQ;
lambda = 0.05 / sqrt(max(n1, n2)*n3); %ST
maxIter = 200;
tol = 1e-3;
[X_Our, E_Our, V_Our, rel_error] = TRPCA_LPRN(0.6, lambda, mO, maxIter, 1e-5, opts);
time_LPRN = toc; %imshow(V_Our(:,:,42)+0.5)
[S_LPRN, P_LPRN, R_LPRN, F_LPRN] = S_value(X, X_Our, E_Our, S0);
%imshow(S_LPRN(:,:,42))


