%% Fact 1 - Power System
%  Saves separately:
%  s1c1.eps, s1c2.eps, s1c3.eps
%  s1_states.eps, s1_inputs.eps, s1_outputs.eps
close all; clc; clear;

A = [0.883878 0       0.047013;
     0.089407 0.904838 0.002321;
     0.056415 0       1.000939];
B = [0.11753; 0.00559; 0.00359];
Q_true = [1 0 0; 0 10 0; 0 0 3];
R_true = 2;
[~,~,Kstar] = dare(A, B, Q_true, R_true);
n = 3; m = 1;
C1 = [1 0 1];   

errK = cell(3,1); valR = cell(3,1); valQ = cell(3,1);
K_final = Kstar;

%% ---------------------------------------------------------------
%  Simulation
%% ---------------------------------------------------------------
for caseID = 1:3
    fprintf('\n===== Case %d =====\n', caseID);
    K = 0.5*Kstar; SampleNumber = 30; maxIter = 100; epsl = 1e-3;
    x = [0.1;0.5;0.5]; xe = [0.1;0.5;0.5];
    errorK=[]; valueQ=[]; valueR=[];

    for iter = 1:maxIter
        u=zeros(m,SampleNumber); ue=zeros(m,SampleNumber);
        xt=zeros(n,SampleNumber+1); xet=zeros(n,SampleNumber+1);
        xt(:,1)=x; xet(:,1)=xe;
        for i=1:SampleNumber
            u(:,i)  = -K    *xt(:,i)  + 0.1*sum([8.5*sin(1.5*i) 5.0*sin(0.1*i) 1.2*sin(0.5*i) 2.0*sin(2.6*i)]);
            ue(:,i) = -Kstar*xet(:,i) + 0.1*sum([8.5*sin(1.1*i) 5.0*sin(0.1*i) 1.0*sin(0.5*i) 2.0*sin(2.6*i)]);
            xt(:,i+1)=A*xt(:,i)+B*u(:,i); xet(:,i+1)=A*xet(:,i)+B*ue(:,i);
        end
        idx=1:SampleNumber-1;
        xx0=xt(:,idx); uu0=u(:,idx); xxe=xet(:,idx); uue=ue(:,idx);
        P=sdpvar(n+m,n+m,'symmetric'); Pe=sdpvar(n+m,n+m,'symmetric');
        Q=sdpvar(n,n,'symmetric'); R=sdpvar(m,m,'symmetric'); tau=sdpvar(1);
        Qb=blkdiag(Q,R);
        Zl=[xx0(:,1:end-1);uu0(:,1:end-1)]; Yl=[xx0(:,2:end);-K*xx0(:,2:end)];
        Zel=[xxe(:,1:end-1);uue(:,1:end-1)]; Ye=[xxe(:,2:end);-Kstar*xxe(:,2:end)];
        c11=Zl'*P*Zl-Yl'*P*Yl-Zl'*Qb*Zl; d11=Zel'*Pe*Zel-Ye'*Pe*Ye-Zel'*Qb*Zel;
        switch caseID
            case 1; obj=tau^2; cc=[eye(n+m)<=Qb<=tau*eye(n+m)];
            case 2; obj=norm(Q,'fro')^2; cc=[R==1,eye(n)<=Q,Q<=tau*eye(n),tau>=1];
            case 3; obj=norm(Q,'fro')^2+norm(R,'fro')^2; cc=[trace(Q)+trace(R)==1,Q>=0,R>=1e-6*eye(m)];
        end
        CONS=[P>=0,Pe>=0,Q>=0,P-Pe>=0,c11==0,d11==0,cc];
        sol=optimize(CONS,obj,sdpsettings('solver','mosek','verbose',0));
        if sol.problem~=0, fprintf('Infeasible at iter %d\n',iter); break; end
        Q_v=value(Q); R_v=value(R);
        H=double(P); K=inv(H(n+1:end,n+1:end))*H(1:n,n+1:end)';
        errorK=[errorK;norm(K-Kstar)]; valueQ=[valueQ;norm(Q_v)]; valueR=[valueR;norm(R_v)];
        fprintf('Iter %d | ||K-K*||=%.6f | R=%.4f\n',iter,norm(K-Kstar),norm(R_v));
        if norm(K-Kstar,'fro')<epsl, fprintf('Converged at iter %d\n',iter); break; end
    end
    errK{caseID}=errorK; valR{caseID}=valueR; valQ{caseID}=valueQ;
    if caseID==1, K_final=K; end
end

%% ---------------------------------------------------------------
%  Trajectory
%% ---------------------------------------------------------------
T=250; x0=[2;2;1];
xe=zeros(n,T+1); xi=zeros(n,T+1);
ue_t=zeros(m,T); ui_t=zeros(m,T);
xe(:,1)=x0; xi(:,1)=x0;
for k=1:T
    ue_t(k)=-Kstar *xe(:,k); ui_t(k)=-K_final*xi(:,k);
    xe(:,k+1)=A*xe(:,k)+B*ue_t(k); xi(:,k+1)=A*xi(:,k)+B*ui_t(k);
end
ye=C1*xe(:,1:end-1); yi=C1*xi(:,1:end-1);   % 1xT vectors
t=0:T; tk=0:T-1;

%% ---------------------------------------------------------------
%  Colors
%% ---------------------------------------------------------------
c_K=[0 0.447 0.741]; c_R=[0.85 0.325 0.098]; c_Q=[0.466 0.674 0.188];
c_e=[0.2 0.4 0.8];   c_i=[0.85 0.33 0.10];
sc ={[0 0.45 0.74],[0.85 0.33 0.10],[0.93 0.69 0.13]};
sci={[0.49 0.18 0.56],[0 0 0],[0 0.75 0.75]};
oc_e={[0 0.45 0.74],[0.85 0.33 0.10]};
oc_i={[0.49 0.18 0.56],[0.93 0.69 0.13]};
lw=1.5; ms=4;

gset = @(ax) set(ax,'GridLineStyle','--','FontName','Times New Roman',...
    'FontSize',13,'FontWeight','bold','LineWidth',0.8);

%% ---------------------------------------------------------------
%  COMPACT FIGURE: 9x2 grid
%
%  Left column:
%  K error, R, Q for Cases 1-3
%
%  Right column:
%  States, Inputs, Outputs
%
%  Modified to reduce vertical spacing between subplots.
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
    % K convergence
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

    xlim([1 length(errK{row})])

    set(gca, ...
        'XTick',1:length(errK{row}))

    legend( ...
        '$\|K_\tau-K^*\|$', ...
        'Interpreter','latex', ...
        'Location','northeast', ...
        'FontSize',13)

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

    xlim([1 length(valR{row})])

    set(gca, ...
        'XTick',1:length(valR{row}))

    legend( ...
        '$R_\tau$', ...
        'Interpreter','latex', ...
        'Location','southeast', ...
        'FontSize',13)


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

    xlim([1 length(valQ{row})])

    set(gca, ...
        'XTick',1:length(valQ{row}))

    legend( ...
        '$\|Q_\tau\|$', ...
        'Interpreter','latex', ...
        'Location','southeast', ...
        'FontSize',13)

    xlabel( ...
        'Iteration $\tau$', ...
        'Interpreter','latex', ...
        'FontName','Times New Roman', ...
        'FontSize',13)


    %% ===========================================================
    % RIGHT COLUMN
    %
    % Each right-side plot spans the same three rows as the
    % corresponding three convergence plots on the left.
    %% ===========================================================

    nexttile(base+2,[3 1])


    %% -----------------------------------------------------------
    % STATES
    %% -----------------------------------------------------------

    if row == 1

        hold on

        for s = 1:3

            plot( ...
                t, ...
                xe(s,:), ...
                '-', ...
                'Color',sc{s}, ...
                'LineWidth',lw, ...
                'DisplayName',sprintf('$x_{k_%d}^*$',s));

        end


        for s = 1:3

            plot( ...
                t, ...
                xi(s,:), ...
                '--', ...
                'Color',sci{s}, ...
                'LineWidth',lw, ...
                'DisplayName',sprintf('$x_{k_%d}$',s));

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

        plot( ...
            tk, ...
            ue_t, ...
            '-', ...
            'Color',c_e, ...
            'LineWidth',lw, ...
            'DisplayName','$u_k^*$');

        hold on

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

        plot( ...
            tk, ...
            ye, ...
            '-', ...
            'Color',c_e, ...
            'LineWidth',lw, ...
            'DisplayName','$y_k^*$');

        hold on

        plot( ...
            tk, ...
            yi, ...
            '--', ...
            'Color',c_i, ...
            'LineWidth',lw, ...
            'DisplayName','$y_k$');


        max_y = ...
            max(abs([ye(:);yi(:)]))*1.15;


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
    % Common right-side formatting
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
    'fact1_fig_compact', ...
    '-depsc', ...
    '-r300');

fprintf('Saved fact1_fig_compact.eps\n');