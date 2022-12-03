function [Xe, rel_error, iter] = LRTC_LPRN(p, B, Idx, maxIter, tol, X)

%%
%Non-convex Alternative Direction Method of Multipliers for low-Turcker rank recovery problem
%
% Input
%      X      --- the true tensor
%      Idx    --- vector containing the locations of the known entries
%      B      --- vector containing the values at the known locations
%      p      --- the value of p for Lp-norm

% Output
%      Xk     --- solution tensor returned by NCADMM
%      Y:     --- set of tensors, representing the extra Y variables used in the algorithm
%      iter:     --- number of iterations after which Xk was obtained (i <= maxiter)

% Copyright by Kun Shang

dimX = size(X);
orderX = 3;
normX = norm(X(:));

% --- Parameter Setting ---
%beta=0.01; rho_beta = 2; beta_max = 1e10;
beta = 0.001;
rho_beta = 10;
beta_max = 1e10;
ADJUST_PARAMETERS = 1; % do not adjust parameters in the first iteration
COUNT_ADJUST_PARAM = 0; % number of changes of lambda/beta in sequence
BETA_MAX_REACHED = 0;
MAX_COUNT_ADJUST_PARAM = 40;
INCREASE_BETA = 1;
TOLERANCE_REL_DELTA_X = 9e-3;

% --- Initialization(assignment space) ---
Xe = zeros(dimX);
Xe(Idx) = B; % solution tensor
Yi = zeros([orderX, dimX]); % variables from variable splitting, put together in a tensor
Ki = zeros([orderX, dimX]); % corresponding Lagrange multiplier
r = 100;
tau = 1e-4;

%% --- Main algorithm ---
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
                beta = rho_beta * beta; %all beta values are increased identically
                r = max(r/1.05, 1e-5);
                tau = max(tau/2, 1e-5);
                % fprintf('beta adjusted, beta(1) =%g,\n', beta);
            else
                BETA_MAX_REACHED = 1;
            end
        end
        if ((~INCREASE_BETA) || BETA_MAX_REACHED)
            % fprintf('\nmaximal parameter values reached, stopping\n');
            break;
        end
    end

    if COUNT_ADJUST_PARAM == MAX_COUNT_ADJUST_PARAM
        % fprintf('\nadjusted parameters %i times in a row, stopping', MAX_COUNT_ADJUST_PARAM);
        break;
    end

    nu = 1 / beta;
    % ----------------------------- update Xk -----------------------------
    Xe0 = Xe;
    Y_tmp = zeros(dimX);
    for iY = 1:orderX
        Y_tmp = Y_tmp + reshape(beta*Yi(iY, :), dimX);
    end
    sumW_plus_sumBetaY = reshape(sum(Ki, 1), dimX) + Y_tmp;
    BS = sumW_plus_sumBetaY / (orderX * beta);
    Xe = sign(BS) .* max(0, abs(BS)-tau*nu/orderX);
    Xe(Idx) = B;
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
            if i == 3
                SP = ((1 + r) * Sigma ./ (r + Sigma)).^p;
            else
                SP = 1;
            end
            %             SP=1;
            S = Sigma - r * (1 + r) * p * ((r + 1) .* SP .* Sigma ./ (Sigma + r)).^(p - 1) * nu ./ (r + Sigma).^2;
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

    % ---------------- Update Lagrange multiplier Ki ----------------------
    for ik = 1:orderX
        Ki(ik, :) = Ki(ik, :) - beta * (Xe(:).' - Yi(ik, :));
    end

    % --- Convergence criterion ---
    rel_deltaX = norm(Xe(:)-Xe0(:), inf);
    stop(4) = rel_deltaX;
    for iy = 1:orderX stop(iy) = norm(Xe(:)-Yi(iy, :)', inf); end
    for is = 1:4 rel_error{is}(iter) = stop(is); end
    if max(stop(:)) < tol
        break
    end
    error(iter) = norm(X(:)-Xe(:)) / normX;
end
