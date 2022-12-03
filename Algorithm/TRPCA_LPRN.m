function [Xe, mE, mV, rel_error] = TRPCA_LPRN(p, lambda, mO, maxIter, tol, opts)

dimX = size(mO);
orderX = 3;
% 加载时空张量
I1 = opts.I1;
I2 = opts.I2;
I3 = opts.I3;
T1SQ = opts.T1SQ;
T2SQ = opts.T2SQ;
T3SQ = opts.T3SQ;

% --- Parameter Setting ---
beta = 1e-04;
rho_beta = 1.2;
beta_max = 1e8; % stop increasing when this value is reached
ADJUST_PARAMETERS = 1; % do not adjust parameters in the first iteration
COUNT_ADJUST_PARAM = 0; % number of changes of lambda/beta in sequence
BETA_MAX_REACHED = 0;
INCREASE_BETA = 1;
TOLERANCE_REL_DELTA_X = 1.3;


% --- Initialization(assignment space) ---
Xe = randn(dimX); % solution tensor
Yi = zeros([orderX, dimX]);
Zi = zeros([orderX, dimX]); % variables from variable splitting, put together in a tensor
Ki = zeros([orderX, dimX]);
Wi = zeros([orderX, dimX]); % corresponding Lagrange multiplier
mS = zeros(dimX);
mE = zeros(dimX);
mV = zeros(dimX);
r = 100;
eta1 = 100;
eta2 = 100;
eta3 = 500;
rho = 0.05;
delta = 1e-2;
%eta1=1;eta2=100;eta3=400;rho=1;

%% --- Main algorithm ---
tic
for iter = 1:maxIter
    if ~ADJUST_PARAMETERS % did not adjust parameters in the last iteration
        if (rel_deltaX > 1.2 || rel_deltaX < TOLERANCE_REL_DELTA_X)
            ADJUST_PARAMETERS = 1;
        else
            COUNT_ADJUST_PARAM = 0;
        end
    else
        ADJUST_PARAMETERS = 0; % adjusted params in the last step, keep for this iteration
        COUNT_ADJUST_PARAM = COUNT_ADJUST_PARAM + 1;
    end

    %adjust parameters
    if ADJUST_PARAMETERS
        if INCREASE_BETA
            if beta < beta_max
                beta = min(rho_beta*beta, beta_max);
                r = max(r/2, eps); %delta=min(beta,0.005); %all beta values are increased identically
                eta1 = max(eta1/1.1, 1e-03);
                eta2 = max(eta2/1.1, 1e-03);
                eta3 = max(eta3/1.2, 1e-06);
            else
                BETA_MAX_REACHED = 1;
            end
        end
        if ((~INCREASE_BETA) || BETA_MAX_REACHED)
            break;
        end
    end

    nu = 1 / beta;
    % ----------------------------- update mEk -----------------------------
    OXS = mO - Xe - mV - mS / beta;
    OXS = OXS / 3 + reshape(Zi(2, :)+Zi(3, :)+Wi(2, :)/beta+Wi(3, :)/beta, dimX) / 3;
    oldmE = mE;
    mE = prox_l1((3 * beta * OXS - delta * mV)/3/beta, lambda/beta/3); %内积
    mE = mE + 0.4 * (mE - oldmE);
    % ---------------- Update mV ----------------------
    mV = (beta * (mO - Xe - mE - mS / beta) - delta * mE) / (beta + 2 * rho); %内积
    % ----------------------------- update Xk -----------------------------
    Xe0 = Xe;
    YZ_tmp = zeros(dimX);
    for iYZ = 1:orderX
        YZ_tmp = YZ_tmp + reshape(beta*Yi(iYZ, :), dimX);
    end
    YZ_tmp = YZ_tmp + reshape(beta*Zi(1, :), dimX);
    sumW_plus_sumBetaY = reshape(sum(Ki, 1), dimX) + reshape(Wi(1, :), dimX) + YZ_tmp;

    Xe = sumW_plus_sumBetaY / (beta * 5) - (mE + mV - mO + mS / beta) / 5;
    Xe = Xe + 0.4 * (Xe - Xe0);

    % ---------------------------------------------------------------------
    for i = 1:orderX
        % --- compute Yi ---
        Yi_unfold = fft(Xe-reshape(Ki(i, :)/beta, dimX), [], i);
        for u = 1:ceil((dimX(i) + 1)/2)
            if i == 1
                [U, Sigma, V] = svd(squeeze(Yi_unfold(u, :, :)), 'econ');
            elseif i == 2
                [U, Sigma, V] = svd(squeeze(Yi_unfold(:, u, :)), 'econ');
            else
                [U, Sigma, V] = svd(Yi_unfold(:, :, u), 'econ');
            end
            Sigma = diag(Sigma);
            S = Sigma - r * (1 + r) * p * ((r + 1) .* Sigma ./ (Sigma + r)).^(p - 1) * nu ./ (r + Sigma).^2;
            S(S < 0) = 0;
            S = diag(S);
            Yi_shrink = U * S * V';
            if i == 1
                Y_tmp(u, :, :) = Yi_shrink;
            elseif i == 2
                Y_tmp(:, u, :) = Yi_shrink;
            else
                Y_tmp(:, :, u) = Yi_shrink;
            end
        end
        for u = ceil((dimX(i) + 1)/2) + 1:dimX(i)
            if i == 1
                Y_tmp(u, :, :) = conj(Y_tmp(dimX(i)-u+2, :, :));
            elseif i == 2
                Y_tmp(:, u, :) = conj(Y_tmp(:, dimX(i)-u+2, :));
            else
                Y_tmp(:, :, u) = conj(Y_tmp(:, :, dimX(i)-u+2));
            end
        end
        Yi(i, :, :, :) = ifft(Y_tmp, [], i);
    end

    for i = 1:orderX
        % --- compute Zi ---
        if i == 1
            Zi_unfold = fft(beta*Xe-reshape(Wi(i, :), dimX), [], i);
        else
            Zi_unfold = fft(beta*mE-reshape(Wi(i, :), dimX), [], i);
        end
        for u = 1:ceil((dimX(i) + 1)/2)
            if i == 1
                Z_tmp(u, :, :) = squeeze(Zi_unfold(u, :, :)) / (eta3 * T3SQ{u} + beta * I3);
            elseif i == 2
                Z_tmp(:, u, :) = (eta1 * T1SQ{u} + beta * I1) \ squeeze(Zi_unfold(:, u, :));
            else
                Z_tmp(:, :, u) = squeeze(Zi_unfold(:, :, u)) / (eta2 * T2SQ{u} + beta * I2);
            end
        end
        for u = ceil((dimX(i) + 1)/2) + 1:dimX(i)
            if i == 1
                Z_tmp(u, :, :) = conj(Z_tmp(dimX(i)-u+2, :, :));
            elseif i == 2
                Z_tmp(:, u, :) = conj(Z_tmp(:, dimX(i)-u+2, :));
            else
                Z_tmp(:, :, u) = conj(Z_tmp(:, :, dimX(i)-u+2));
            end
        end
        Zi(i, :, :, :) = ifft(Z_tmp, [], i);
    end
    % ---------------- Update Lagrange multiplier mS ----------------------
    mS = mS + beta * (Xe + mE + mV - mO);
    % ---------------- Update Lagrange multiplier Ki ----------------------
    for ik = 1:orderX
        Ki(ik, :) = Ki(ik, :) - beta * (Xe(:).' - Yi(ik, :));
    end
    % ---------------- Update Lagrange multiplier Wi ----------------------
    for iw = 1:orderX
        if iw == 1
            Wi(iw, :) = Wi(iw, :) - beta * (Xe(:).' - Zi(iw, :));
        else
            Wi(iw, :) = Wi(iw, :) - beta * (mE(:).' - Zi(iw, :));
        end
    end
    % --- Convergence criterion ---
    rel_deltaX = norm(Xe(:)-Xe0(:), inf);
    stop(4) = rel_deltaX;
    stop(5) = norm(Xe(:)-Zi(1, :)', inf);
    for iy = 1:orderX stop(iy) = norm(Xe(:)-Yi(iy, :)', inf); end
    for is = 1:5 rel_error{is}(iter) = stop(is); end
    if max(stop(:)) < tol
        break
    end
    gh(iter) = rel_deltaX
end
