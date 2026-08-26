function [] = explainability_validation()
%% ============================================================
%  Explainability Validation (Corrected)
%  Validates CONSTRAINT SATISFACTION and INTERPRETABLE QUANTITY
%  for each case. Does NOT compare against Q* directly since
%  Q* may not be PSD in general.
%% ============================================================

%% --- Setting 1 Recovered Weights at Convergence ---
Q_S1_C1 = [1.0354 0.4304 0.2037;
            0.4304 6.2287 2.5237;
            0.2037 2.5237 3.0859];  R_S1_C1 = 2.0152;

Q_S1_C2 = [1.0628 0.5260 0.3444;
            0.5260 5.4053 2.8955;
            0.3444 2.8955 3.0321];  R_S1_C2 = 2.0000;

Q_S1_C3 = [0.0831  0.1175 -0.0013;
            0.1175  0.5045  0.0690;
           -0.0013  0.0690  0.2483]; R_S1_C3 = 0.1641;

%% --- Setting 2 Recovered Weights at Convergence ---
Q_S2_C1 = [7.5918 1.0796 3.7102;
            1.0796 1.1791 0.6158;
            3.7102 0.6158 3.1165];  R_S2_C1 = 2.3743;

Q_S2_C2 = [7.9174 1.2615 3.7544;
            1.2615 1.2301 0.6847;
            3.7544 0.6847 3.0377];  R_S2_C2 = 2.0000;

Q_S2_C3 = [0.3425 -0.0447  0.2124;
           -0.0447  0.3488  0.0439;
            0.2124  0.0439  0.1834]; R_S2_C3 = 0.1282;

%% ============================================================
fprintf('\n======= EXPLAINABILITY VALIDATION =======\n\n');

%% --- Case 1: Spectral Box ---
% Constraint: lambda_min(diag(Q,R)) = 1
% Interpretable quantity: kappa = lambda_max / lambda_min
fprintf('--- Case 1: Spectral Box ---\n');
[lmin_S1, lmax_S1, kap_S1] = spectral_metrics(Q_S1_C1, R_S1_C1);
[lmin_S2, lmax_S2, kap_S2] = spectral_metrics(Q_S2_C1, R_S2_C1);
fprintf('Setting 1: lambda_min = %.4f | lambda_max = %.4f | kappa = %.4f\n', ...
    lmin_S1, lmax_S1, kap_S1);
fprintf('Setting 2: lambda_min = %.4f | lambda_max = %.4f | kappa = %.4f\n\n', ...
    lmin_S2, lmax_S2, kap_S2);

%% --- Case 2: Fixed R ---
% Constraint: R_tau = R* (fixed)
% Interpretable quantity: ||Q||_F (complexity of state penalty)
fprintf('--- Case 2: Fixed R ---\n');
fprintf('Setting 1: R_tau = %.4f | ||Q||_F = %.4f\n', R_S1_C2, norm(Q_S1_C2,'fro'));
fprintf('Setting 2: R_tau = %.4f | ||Q||_F = %.4f\n\n', R_S2_C2, norm(Q_S2_C2,'fro'));

%% --- Case 3: Trace Normalization ---
% Constraint: tr(Q) + tr(R) = 1
% Interpretable quantity: rho = tr(Q)/tr(R)
fprintf('--- Case 3: Trace Normalization ---\n');
tr_sum_S1 = trace(Q_S1_C3) + R_S1_C3;
rho_S1    = trace(Q_S1_C3) / R_S1_C3;
tr_sum_S2 = trace(Q_S2_C3) + R_S2_C3;
rho_S2    = trace(Q_S2_C3) / R_S2_C3;
fprintf('Setting 1: tr(Q)+tr(R) = %.4f | rho = tr(Q)/tr(R) = %.4f\n', ...
    tr_sum_S1, rho_S1);
fprintf('Setting 2: tr(Q)+tr(R) = %.4f | rho = tr(Q)/tr(R) = %.4f\n\n', ...
    tr_sum_S2, rho_S2);

fprintf('Physical interpretation of rho:\n');
fprintf('Setting 1: expert prioritizes state %.1fx more than control\n', rho_S1);
fprintf('Setting 2: expert prioritizes state %.1fx more than control\n', rho_S2);
end

%% ---------------------------------------------------------------
function [lmin, lmax, kap] = spectral_metrics(Q, R)
    eigs_all = [eig(Q); R];
    eigs_pos = sort(eigs_all(eigs_all > 1e-10));
    lmin = eigs_pos(1);
    lmax = eigs_pos(end);
    kap  = lmax / lmin;
end