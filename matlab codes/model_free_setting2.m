%% ===============================================================
% FACT 2 - ALGORITHM 3
% Input-output-driven model-free CO-PI
%
% Manuscript convention:
%
%       u_k = -Kbar*z_k
%
% Expert:
%       x_k^*, u_k^*, y_k^*
%
% Learner:
%       x_k, u_k, y_k
%
% Requires:
%   Control System Toolbox
%   YALMIP
%   MOSEK
%% ===============================================================

close all;
clc;
clear;
yalmip('clear');


%% ===============================================================
% AIRCRAFT SYSTEM
%% ===============================================================

A = [ ...
    0.90649   0.08160  -0.00050;
    0.07413   0.90121  -0.00071;
    0         0         0.13266];

B = [ ...
   -0.00151;
   -0.00960;
    0.86730];

C = [1 0 1];

n = size(A,1);
m = size(B,2);
p = size(C,1);


%% ===============================================================
% EXPERT COST
%% ===============================================================

Q1 = 3;
R1 = 2;


%% ===============================================================
% STATE RECONSTRUCTION
%% ===============================================================

N = 3;

Obs = [ ...
    C*A*A;
    C*A;
    C];

My = A^3/Obs;

Un = [ ...
    B, ...
    A*B, ...
    A^2*B];

Tn = [ ...
    0, C*B, C*A*B;
    0, 0,   C*B;
    0, 0,   0];

Theta = [ ...
    Un-My*Tn, ...
    My];

dz = size(Theta,2);     % = 6

dH = dz + m;            % = 7


%% ===============================================================
% EXPERT POLICY
%
% MATLAB dlqr:
%
%       u = -K*x
%
% Since
%
%       x = Theta*z
%
% then
%
%       u = -(K*Theta)z
%
% Therefore:
%
%       Kbar* = K*Theta
%% ===============================================================

[Kstar,~,~] = dlqr( ...
    A, ...
    B, ...
    C'*Q1*C, ...
    R1);

Kbarstar = Kstar*Theta;


%% ===============================================================
% DIRECT INITIAL HISTORY GAIN
%
% This corresponds to the original simulation initialization,
% but is written DIRECTLY instead of using Kbar*.
%% ===============================================================

Kbar0 = [ ...
     0.0553400841, ...
    -1.0117702780, ...
     0.5473341028, ...
     1.1756110151, ...
    -0.7856520842, ...
     0.0836874216];


%% ===============================================================
% SETTINGS
%% ===============================================================

SampleNumber = 200;

maxIter = 100;

epsl = 1e-2;

x1_0 = [ ...
    1;
   -1;
    0.1];


% Numerical tolerances used only for the SDP implementation.
% These do not change the theoretical Algorithm 3 constraints.

% Case 3 small positive-definite lower bound
epsQ = 1e-6;
epsR = 1e-6;

% Bellman residual tolerance after row scaling
bellTol = 1e-7;

% Small PSD feasibility tolerance used to prevent roundoff-level
% negative eigenvalues from causing MOSEK numerical failure
psdTol = 1e-8;


%% ===============================================================
% BASIC INFORMATION
%% ===============================================================

fprintf('============================================================\n');

fprintf('FACT 2 - ALGORITHM 3\n');

fprintf('============================================================\n');


fprintf('rank(Obs) = %d / %d\n', ...
    rank(Obs),n);


fprintf('cond(Obs) = %.6e\n', ...
    cond(Obs));


fprintf('rank(Theta) = %d\n', ...
    rank(Theta));


fprintf('\nKstar =\n');

disp(Kstar);


fprintf('Kbarstar =\n');

disp(Kbarstar);


fprintf('Kbar0 =\n');

disp(Kbar0);


fprintf( ...
    'Initial ||Kbar0-Kbarstar|| = %.6e\n', ...
    norm(Kbar0-Kbarstar,'fro'));


%% ===============================================================
% STORAGE
%% ===============================================================

errK = cell(3,1);

valR = cell(3,1);

valQ = cell(3,1);

policyStepHist = cell(3,1);


Qfinal = cell(3,1);

Rfinal = cell(3,1);

Kfinal = cell(3,1);


Kbar_final = Kbarstar;


%% ===============================================================
% THREE CASES
%% ===============================================================

for caseID = 1:3


    fprintf('\n============================================================\n');

    fprintf('CASE %d\n',caseID);

    fprintf('============================================================\n');


    %% Initial policy

    Kbar = Kbar0;


    errorK = [];

    valueQ = [];

    valueR = [];

    valueStep = [];


    Q_v_last = [];

    R_v_last = [];


    %% ===========================================================
    % POLICY ITERATION
    %% ===========================================================

    for iter = 1:maxIter


        Kbar_old = Kbar;


        %% =======================================================
        % DATA COLLECTION
        %% =======================================================

        x = zeros(n,SampleNumber+1);

        xe = zeros(n,SampleNumber+1);


        y = zeros(p,SampleNumber);

        ye = zeros(p,SampleNumber);


        u = zeros(m,SampleNumber);

        ue = zeros(m,SampleNumber);


        z = zeros(dz,SampleNumber);

        ze = zeros(dz,SampleNumber);


        x(:,1) = x1_0;

        xe(:,1) = x1_0;


        for k = 1:SampleNumber


            %% Current outputs

            y(:,k) = C*x(:,k);

            ye(:,k) = C*xe(:,k);


            if k >= N+1


                %% Learner history

                z(:,k) = [ ...
                    u(:,k-1);
                    u(:,k-2);
                    u(:,k-3);
                    y(:,k-1);
                    y(:,k-2);
                    y(:,k-3)];


                %% Expert history

                ze(:,k) = [ ...
                    ue(:,k-1);
                    ue(:,k-2);
                    ue(:,k-3);
                    ye(:,k-1);
                    ye(:,k-2);
                    ye(:,k-3)];


                %% Probing noise

                excitation = 0.3*sum([ ...
                    sin(k), ...
                    sin(10*k), ...
                    sin(0.1*k), ...
                    sin(4.5*k), ...
                    sin(6*k), ...
                    sin(18*k), ...
                    sin(0.01*k)]);


                %% Learner input
                %
                % Paper convention:
                %
                %       u = -Kbar*z
                %

                u(:,k) = ...
                    -Kbar*z(:,k) ...
                    + excitation;


                %% Expert input

                ue(:,k) = ...
                    -Kbarstar*ze(:,k) ...
                    + excitation;


            end


            %% State propagation

            x(:,k+1) = ...
                A*x(:,k) ...
                + B*u(:,k);


            xe(:,k+1) = ...
                A*xe(:,k) ...
                + B*ue(:,k);


        end


        %% =======================================================
        % SDP VARIABLES
        %% =======================================================

        H = sdpvar( ...
            dH,dH,'symmetric');


        He = sdpvar( ...
            dH,dH,'symmetric');


        Q = sdpvar( ...
            n,n,'symmetric');


        R = sdpvar( ...
            m,m,'symmetric');


        eta = sdpvar(1);


        %% =======================================================
        % HISTORY-SPACE COST
        %% =======================================================

        Qbar = ...
            Theta'*Q*Theta;


        LambdaBar = ...
            blkdiag(Qbar,R);


        %% =======================================================
        % SCALED BELLMAN EQUATIONS
        %
        % The theoretical equations are exact.  Numerically, each
        % residual is scaled by the energy of its data row and then
        % enforced within +/- bellTol.  This avoids artificial MOSEK
        % failures when learner and expert equations become nearly
        % linearly dependent close to policy convergence.
        %% =======================================================

        bellL_scaled = [];

        bellE_scaled = [];


        for k = N+1:SampleNumber-1


            %% Current learner augmented vector

            zeta = [ ...
                z(:,k);
                u(:,k)];


            %% Current expert augmented vector

            zeta_e = [ ...
                ze(:,k);
                ue(:,k)];


            %% Next histories

            znext = ...
                z(:,k+1);


            zenext = ...
                ze(:,k+1);


            %% Learner target-policy successor

            zeta_next = [ ...
                znext;
                -Kbar*znext];


            %% Expert successor

            zeta_e_next = [ ...
                zenext;
                -Kbarstar*zenext];


            %% ---------------------------------------------------
            % Learner Bellman residual
            %% ---------------------------------------------------

            learnerResidual = ...
                zeta'*H*zeta ...
                - zeta_next'*H*zeta_next ...
                - zeta'*LambdaBar*zeta;


            %% ---------------------------------------------------
            % Expert Bellman residual
            %% ---------------------------------------------------

            expertResidual = ...
                zeta_e'*He*zeta_e ...
                - zeta_e_next'*He*zeta_e_next ...
                - zeta_e'*LambdaBar*zeta_e;


            %% ---------------------------------------------------
            % Numerical row scaling
            %% ---------------------------------------------------

            scaleL = max([ ...
                norm(zeta)^2, ...
                norm(zeta_next)^2, ...
                1e-3]);


            scaleE = max([ ...
                norm(zeta_e)^2, ...
                norm(zeta_e_next)^2, ...
                1e-3]);


            % Scaled Bellman residuals

            bellL_scaled = [ ...
                bellL_scaled;
                learnerResidual/scaleL];


            bellE_scaled = [ ...
                bellE_scaled;
                expertResidual/scaleE];


        end


        %% =======================================================
        % CASE CONSTRAINTS
        %% =======================================================

        switch caseID


            %% ===================================================
            % CASE 1
            % Spectral box
            %% ===================================================

            case 1


                objective = ...
                    eta^2;


                structuralConstraints = [ ...
                    blkdiag(Q,R) >= eye(n+m), ...
                    blkdiag(Q,R) <= eta*eye(n+m), ...
                    eta >= 1, ...
                    eta <= 100];


            %% ===================================================
            % CASE 2
            % Fixed R
            %% ===================================================

            case 2


                objective = ...
                    norm(Q,'fro')^2;


                structuralConstraints = [ ...
                    R == R1, ...
                    Q >= eye(n), ...
                    Q <= eta*eye(n), ...
                    eta >= 1, ...
                    eta <= 50];


            %% ===================================================
            % CASE 3
            % Trace normalization
            %% ===================================================

            case 3


                objective = ...
                    norm(Q,'fro')^2 ...
                    + norm(R,'fro')^2;


                structuralConstraints = [ ...
                    trace(Q)+trace(R) == 1, ...
                    Q >= epsQ*eye(n), ...
                    R >= epsR*eye(m)];


        end


        %% =======================================================
        % COMMON CONSTRAINTS
        %
        % The small -psdTol allowance is purely numerical.  It lets
        % matrices with roundoff-level negative eigenvalues (for
        % example -1e-10) be treated as PSD by the finite-precision
        % solver.
        %% =======================================================

        commonConstraints = [ ...
            H + psdTol*eye(dH) >= 0, ...
            He + psdTol*eye(dH) >= 0, ...
            H-He + psdTol*eye(dH) >= 0];


        %% =======================================================
        % MOSEK OPTIONS
        %% =======================================================

        options = sdpsettings( ...
            'solver','mosek', ...
            'verbose',0, ...
            'savesolverinput',1, ...
            'savesolveroutput',1);


        %% =======================================================
        % NUMERICALLY ROBUST SOLVE
        %
        % Instead of requiring floating-point residuals to be exactly
        % zero, enforce the scaled theoretical Bellman equalities to
        % solver accuracy: |residual| <= bellTol.
        %% =======================================================

        constraints = [ ...
            commonConstraints, ...
            bellL_scaled <= bellTol, ...
            bellL_scaled >= -bellTol, ...
            bellE_scaled <= bellTol, ...
            bellE_scaled >= -bellTol, ...
            structuralConstraints];


        sol = optimize( ...
            constraints, ...
            objective, ...
            options);


        %% =======================================================
        % FINAL SOLVER CHECK
        %% =======================================================

        if sol.problem ~= 0


            fprintf('\n');

            fprintf( ...
                'Solver failed in Case %d at iteration %d\n', ...
                caseID,iter);


            fprintf( ...
                'YALMIP code = %d\n', ...
                sol.problem);


            fprintf( ...
                'YALMIP info = %s\n', ...
                sol.info);


            fprintf('\n');


            break;


        end


        %% =======================================================
        % GET SOLUTION
        %% =======================================================

        Q_v = value(Q);

        R_v = value(R);

        H_v = value(H);

        He_v = value(He);


        % Numerical diagnostics for the accepted SDP solution
        maxBellL = max(abs(value(bellL_scaled)));
        maxBellE = max(abs(value(bellE_scaled)));
        minEigH  = min(eig((H_v+H_v')/2));
        minEigHe = min(eig((He_v+He_v')/2));
        minEigD  = min(eig(((H_v-He_v)+(H_v-He_v)')/2));


        Q_v_last = Q_v;

        R_v_last = R_v;


        %% =======================================================
        % POLICY IMPROVEMENT
        %
        % Partition:
        %
        %       H = [ Hzz  Hzu
        %             Huz  Huu ]
        %
        % Paper convention:
        %
        %       Kbar_new = Huu^{-1} Huz
        %% =======================================================

        Huz = ...
            H_v(end,1:dz);


        Huu = ...
            H_v(end,end);


        %% Safety check

        if ~isfinite(Huu) || abs(Huu) < 1e-10


            fprintf( ...
                'Huu invalid at Case %d, iteration %d\n', ...
                caseID,iter);


            break;


        end


        %% New history gain

        Kbar_new = ...
            Huz/Huu;


        %% =======================================================
        % ERRORS
        %% =======================================================

        gainError = ...
            norm( ...
            Kbar_new-Kbarstar, ...
            'fro');


        policyStep = ...
            norm( ...
            Kbar_new-Kbar_old, ...
            'fro');


        %% Store

        errorK = [ ...
            errorK;
            gainError];


        valueQ = [ ...
            valueQ;
            norm(Q_v,'fro')];


        valueR = [ ...
            valueR;
            norm(R_v,'fro')];


        valueStep = [ ...
            valueStep;
            policyStep];


        %% Display

        fprintf( ...
            ['Iter %2d | ' ...
             '||Kbar-Kbar*|| = %.6f | ' ...
             '||Knew-Kold|| = %.6f | ' ...
             'R = %.4f | ' ...
             'Bell = %.2e\n'], ...
            iter, ...
            gainError, ...
            policyStep, ...
            R_v, ...
            max(maxBellL,maxBellE));


        %% Update gain

        Kbar = ...
            Kbar_new;


        %% =======================================================
        % ALGORITHM 3 STOPPING CRITERION
        %% =======================================================

        if policyStep <= epsl


            fprintf( ...
                'Converged at iteration %d\n', ...
                iter);


            fprintf( ...
                'Final gain error = %.6f\n', ...
                gainError);

            fprintf( ...
                'Max scaled Bellman residual = %.3e\n', ...
                max(maxBellL,maxBellE));

            fprintf( ...
                'min eig(H) = %.3e, min eig(He) = %.3e, min eig(H-He) = %.3e\n', ...
                minEigH,minEigHe,minEigD);


            break;


        end


    end


    %% ===========================================================
    % STORE CASE RESULTS
    %% ===========================================================

    errK{caseID} = ...
        errorK;


    valR{caseID} = ...
        valueR;


    valQ{caseID} = ...
        valueQ;


    policyStepHist{caseID} = ...
        valueStep;


    Qfinal{caseID} = ...
        Q_v_last;


    Rfinal{caseID} = ...
        R_v_last;


    Kfinal{caseID} = ...
        Kbar;


    %% Use Case 1 for trajectory comparison

    if caseID == 1

        Kbar_final = ...
            Kbar;

    end


end


%% ===============================================================
% SUMMARY
%% ===============================================================

fprintf('\n============================================================\n');

fprintf('FACT 2 SUMMARY\n');

fprintf('============================================================\n');


for caseID = 1:3


    fprintf('\nCASE %d\n',caseID);


    if isempty(errK{caseID})


        fprintf('No successful iteration.\n');


        continue;


    end


    fprintf( ...
        'Iterations = %d\n', ...
        length(errK{caseID}));


    fprintf( ...
        'Final gain error = %.6f\n', ...
        errK{caseID}(end));


    fprintf( ...
        'Final policy change = %.6f\n', ...
        policyStepHist{caseID}(end));


    fprintf( ...
        'Final R = %.6f\n', ...
        Rfinal{caseID});


    fprintf('Final Q =\n');

    disp(Qfinal{caseID});


end


%% ===============================================================
% TRAJECTORY COMPARISON
%% ===============================================================

T = 250;


x0 = [ ...
    1;
   -1;
    0.1];


xe = zeros(n,T+1);

xi = zeros(n,T+1);


ue_t = zeros(m,T);

ui_t = zeros(m,T);


ye_t = zeros(p,T);

yi_t = zeros(p,T);


ue_h = zeros(m,T);

ui_h = zeros(m,T);


ye_h = zeros(p,T);

yi_h = zeros(p,T);


xe(:,1) = x0;

xi(:,1) = x0;


for k = 1:T


    %% Current outputs

    ye_t(:,k) = ...
        C*xe(:,k);


    yi_t(:,k) = ...
        C*xi(:,k);


    if k >= N+1


        %% Expert history

        ze_traj = [ ...
            ue_h(:,k-1);
            ue_h(:,k-2);
            ue_h(:,k-3);
            ye_h(:,k-1);
            ye_h(:,k-2);
            ye_h(:,k-3)];


        %% Learner history

        zi_traj = [ ...
            ui_h(:,k-1);
            ui_h(:,k-2);
            ui_h(:,k-3);
            yi_h(:,k-1);
            yi_h(:,k-2);
            yi_h(:,k-3)];


        %% Expert input

        ue_t(:,k) = ...
            -Kbarstar*ze_traj;


        %% Learner input

        ui_t(:,k) = ...
            -Kbar_final*zi_traj;


    end


    %% State propagation

    xe(:,k+1) = ...
        A*xe(:,k) ...
        + B*ue_t(:,k);


    xi(:,k+1) = ...
        A*xi(:,k) ...
        + B*ui_t(:,k);


    %% Save history

    ue_h(:,k) = ...
        ue_t(:,k);


    ui_h(:,k) = ...
        ui_t(:,k);


    ye_h(:,k) = ...
        ye_t(:,k);


    yi_h(:,k) = ...
        yi_t(:,k);


end


t = 0:T;

tk = 0:T-1;


%% ===============================================================
% FIGURE COLORS
%% ===============================================================

c_K = [0 0.447 0.741];

c_R = [0.85 0.325 0.098];

c_Q = [0.466 0.674 0.188];


c_e = [0.2 0.4 0.8];

c_i = [0.85 0.33 0.10];


sc = { ...
    [0 0.45 0.74], ...
    [0.85 0.33 0.10], ...
    [0.93 0.69 0.13]};


sci = { ...
    [0.49 0.18 0.56], ...
    [0 0 0], ...
    [0 0.75 0.75]};


lw = 1.5;

ms = 4;


gset = @(ax) set( ...
    ax, ...
    'GridLineStyle','--', ...
    'FontName','Times New Roman', ...
    'FontSize',13, ...
    'FontWeight','bold', ...
    'LineWidth',0.8);


%% ---------------------------------------------------------------
%  FACT 2 FIGURE
%  SAME STYLE AS FACT 1
%
%  Left column:
%     Kbar error, R, Q
%
%  Right column:
%     States, Inputs, Outputs
%
%  Expert:
%     x_k^*, u_k^*, y_k^*
%
%  Learner:
%     x_k, u_k, y_k
%% ---------------------------------------------------------------

fig = figure( ...
    'Color','w', ...
    'Position',[100 100 900 780]);

tl = tiledlayout( ...
    9, ...
    2, ...
    'TileSpacing','tight', ...
    'Padding','compact');


case_titles = { ...
    'Case 1: spectral box', ...
    'Case 2: fixed $R$', ...
    'Case 3: trace normalization'};


for row = 1:3

    base = (row-1)*6;


    %% ===========================================================
    % LEFT COLUMN
    %% ===========================================================


    %% -----------------------------------------------------------
    % Kbar convergence
    %% -----------------------------------------------------------

    nexttile(base+1)


    plot( ...
        1:length(errK{row}), ...
        errK{row}, ...
        '-o', ...
        'Color',c_K, ...
        'LineWidth',lw, ...
        'MarkerSize',ms, ...
        'MarkerFaceColor','w', ...
        'MarkerEdgeColor',c_K);


    grid on

    gset(gca)


    %% Safe x-axis
    npts = length(errK{row});

    if npts == 1

        xlim([0.5 1.5])

        set(gca,'XTick',1)

    elseif npts > 1

        xlim([1 npts])

        set(gca,'XTick',1:npts)

    end


    if npts > 0
        legend( ...
            '$\|\bar{K}_\tau-\bar{K}^*\|$', ...
            'Interpreter','latex', ...
            'Location','northeast', ...
            'FontSize',13)
    end


    title( ...
        case_titles{row}, ...
        'Interpreter','latex', ...
        'FontName','Times New Roman', ...
        'FontSize',13)



    %% -----------------------------------------------------------
    % R convergence
    %% -----------------------------------------------------------

    nexttile(base+3)


    plot( ...
        1:length(valR{row}), ...
        valR{row}, ...
        '-o', ...
        'Color',c_R, ...
        'LineWidth',lw, ...
        'MarkerSize',ms, ...
        'MarkerFaceColor','w', ...
        'MarkerEdgeColor',c_R);


    grid on

    gset(gca)


    %% Safe x-axis
    npts = length(valR{row});

    if npts == 1

        xlim([0.5 1.5])

        set(gca,'XTick',1)

    elseif npts > 1

        xlim([1 npts])

        set(gca,'XTick',1:npts)

    end


    if npts > 0
        legend( ...
            '$R_\tau$', ...
            'Interpreter','latex', ...
            'Location','southeast', ...
            'FontSize',13)
    end



    %% -----------------------------------------------------------
    % Q convergence
    %% -----------------------------------------------------------

    nexttile(base+5)


    plot( ...
        1:length(valQ{row}), ...
        valQ{row}, ...
        '-o', ...
        'Color',c_Q, ...
        'LineWidth',lw, ...
        'MarkerSize',ms, ...
        'MarkerFaceColor','w', ...
        'MarkerEdgeColor',c_Q);


    grid on

    gset(gca)


    %% Safe x-axis
    npts = length(valQ{row});

    if npts == 1

        xlim([0.5 1.5])

        set(gca,'XTick',1)

    elseif npts > 1

        xlim([1 npts])

        set(gca,'XTick',1:npts)

    end


    if npts > 0
        legend( ...
            '$\|Q_\tau\|_F$', ...
            'Interpreter','latex', ...
            'Location','southeast', ...
            'FontSize',13)
    end


    xlabel( ...
        'Iteration $\tau$', ...
        'Interpreter','latex', ...
        'FontName','Times New Roman', ...
        'FontSize',13)



    %% ===========================================================
    % RIGHT COLUMN
    %
    % Same organization as Fact 1
    %% ===========================================================

    nexttile(base+2,[3 1])


    %% -----------------------------------------------------------
    % STATES
    %% -----------------------------------------------------------

    if row == 1


        hold on


        %% Expert states
        for s = 1:3

            plot( ...
                t, ...
                xe(s,:), ...
                '-', ...
                'Color',sc{s}, ...
                'LineWidth',lw, ...
                'DisplayName', ...
                sprintf('$x_{k_%d}^*$',s));

        end


        %% Learner states
        for s = 1:3

            plot( ...
                t, ...
                xi(s,:), ...
                '--', ...
                'Color',sci{s}, ...
                'LineWidth',lw, ...
                'DisplayName', ...
                sprintf('$x_{k_%d}$',s));

        end


        ylabel( ...
            'States', ...
            'FontName','Times New Roman', ...
            'FontSize',13)


        legend( ...
            'Interpreter','latex', ...
            'Location','northeast', ...
            'FontSize',13, ...
            'NumColumns',2)



    %% -----------------------------------------------------------
    % INPUTS
    %% -----------------------------------------------------------

    elseif row == 2


        %% Expert input
        plot( ...
            tk, ...
            ue_t, ...
            '-', ...
            'Color',c_e, ...
            'LineWidth',lw, ...
            'DisplayName','$u_k^*$');


        hold on


        %% Learner input
        plot( ...
            tk, ...
            ui_t, ...
            '--', ...
            'Color',c_i, ...
            'LineWidth',lw, ...
            'DisplayName','$u_k$');


        max_u = ...
            max(abs([ue_t(:);ui_t(:)]))*1.15;


        if isfinite(max_u) && max_u > 1e-10

            ylim([-max_u max_u])

        end


        ylabel( ...
            'Inputs', ...
            'FontName','Times New Roman', ...
            'FontSize',13)


        legend( ...
            'Interpreter','latex', ...
            'Location','northeast', ...
            'FontSize',13)



    %% -----------------------------------------------------------
    % OUTPUTS
    %% -----------------------------------------------------------

    else


        %% Expert output
        plot( ...
            tk, ...
            ye_t, ...
            '-', ...
            'Color',c_e, ...
            'LineWidth',lw, ...
            'DisplayName','$y_k^*$');


        hold on


        %% Learner output
        plot( ...
            tk, ...
            yi_t, ...
            '--', ...
            'Color',[0 1 1], ...
            'LineWidth',lw, ...
            'DisplayName','$y_k$');


        max_y = ...
            max(abs([ye_t(:);yi_t(:)]))*1.15;


        if isfinite(max_y) && max_y > 1e-10

            ylim([-max_y max_y])

        end


        ylabel( ...
            'Outputs', ...
            'FontName','Times New Roman', ...
            'FontSize',13)


        xlabel( ...
            'Time step $k$', ...
            'Interpreter','latex', ...
            'FontName','Times New Roman', ...
            'FontSize',13)


        legend( ...
            'Interpreter','latex', ...
            'Location','northeast', ...
            'FontSize',13)

    end


    %% -----------------------------------------------------------
    % Common right-column formatting
    %% -----------------------------------------------------------

    grid on

    gset(gca)

    xlim([0 T])

end


%% ===============================================================
% SAVE
%% ===============================================================

print( ...
    fig, ...
    'fact2_fig', ...
    '-depsc', ...
    '-r300');


fprintf('Saved fact2_fig.eps\n');