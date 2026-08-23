close all;
clc;
clear;
yalmip('clear');
rng(2,'twister');

A = [0.90649  0.08160  -0.00050;
     0.07413  0.90121  -0.00071;
     0        0         0.13266];

B = [-0.00151;
     -0.00960;
      0.86730];

C_original = [1 0 1];
Q_original = 3;
R_original = 2;

[Kstar,~,~] = dlqr(A,B,C_original'*Q_original*C_original,R_original);

Ky_star = norm(Kstar,'fro');
Ccmp = Kstar/Ky_star;
policyMismatch = norm(Kstar-Ky_star*Ccmp,'fro');

Obs = [Ccmp*A*A;
       Ccmp*A;
       Ccmp];

tol = 1e-2;

x0 = [1;
     -1;
      0.1];

K0 = [0.4463640849 0.0264790968 0.0648775431];

alpha_set = [1.00 1.25 1.50 2.00 3.00 5.00];

fprintf('\n');
fprintf('============================================================\n');
fprintf('COMPARISON B: HIGHER-ALPHA SENSITIVITY\n');
fprintf('============================================================\n');
fprintf('Kstar = [%.8f %.8f %.8f]\n',Kstar);
fprintf('Ky_star = %.8f\n',Ky_star);
fprintf('Ccmp = [%.8f %.8f %.8f]\n',Ccmp);
fprintf('||Kstar-Ky_star*Ccmp|| = %.3e\n',policyMismatch);
fprintf('rank(Obs) = %d / 3\n',rank(Obs));
fprintf('cond(Obs) = %.3e\n',cond(Obs));
fprintf('rho(A-B*Kstar) = %.6f\n',max(abs(eig(A-B*Kstar))));
fprintf('rho(A-B*K0) = %.6f\n',max(abs(eig(A-B*K0))));
fprintf('Initial proposed error = %.6e\n',norm(K0-Kstar,'fro'));
fprintf('Common threshold = %.1e\n',tol);

fprintf('\n');
fprintf(['NOTE: alpha > 1 is outside the admissible range ' ...
         '(0,1] of the reference paper.\n']);
fprintf('Values above 1 are tested only as a numerical diagnostic.\n');

fprintf('\nAlpha set = ');
fprintf('%.2f ',alpha_set);
fprintf('\n');

fprintf('\n');
fprintf('============================================================\n');
fprintf('PROPOSED ALGORITHM 3\n');
fprintf('============================================================\n');

[err_prop,propConverged] = proposed_algorithm3_same_policy( ...
    A,B,Ccmp,Kstar,K0,x0,tol);

iter_prop = first_hit(err_prop,tol);

nAlpha = numel(alpha_set);

err_ref_all = cell(nAlpha,1);
refConverged = false(nAlpha,1);
iter_ref = nan(nAlpha,1);
final_ref_error = nan(nAlpha,1);
final_Q = nan(nAlpha,1);

fprintf('\n');
fprintf('============================================================\n');
fprintf('HIGHER-ALPHA SWEEP: MODIFIED REFERENCE METHOD\n');
fprintf('============================================================\n');

for ia = 1:nAlpha

    alpha = alpha_set(ia);

    fprintf('\n');
    fprintf('------------------------------------------------------------\n');

    if alpha <= 1
        fprintf('REFERENCE METHOD: alpha = %.2f  [paper-admissible]\n',alpha);
    else
        fprintf('REFERENCE METHOD: alpha = %.2f  [diagnostic only]\n',alpha);
    end

    fprintf('------------------------------------------------------------\n');

    [err_ref_all{ia},refConverged(ia),final_Q(ia)] = ...
        reference_inverseQ_modified( ...
        A,B,Ccmp,Kstar,Ky_star,x0,tol,alpha);

    iter_ref(ia) = first_hit(err_ref_all{ia},tol);
    final_ref_error(ia) = err_ref_all{ia}(end);

    if refConverged(ia)
        fprintf(['alpha = %.2f reached tolerance ' ...
                 'in %d iterations.\n'], ...
                 alpha,iter_ref(ia));
    else
        fprintf(['alpha = %.2f did NOT reach tolerance ' ...
                 'within the maximum iterations.\n'], ...
                 alpha);
    end

end

fprintf('\n');
fprintf('============================================================\n');
fprintf('HIGHER-ALPHA SWEEP SUMMARY\n');
fprintf('============================================================\n');

if propConverged

    fprintf(['Proposed Algorithm 3: %d iterations, ' ...
             'final error = %.6e\n\n'], ...
             iter_prop,err_prop(end));

else

    fprintf(['Proposed Algorithm 3 did not reach tolerance. ' ...
             'Final error = %.6e\n\n'], ...
             err_prop(end));

end

fprintf(' alpha      status        iterations       final error        final Q\n');
fprintf('--------------------------------------------------------------------------\n');

for ia = 1:nAlpha

    if alpha_set(ia) <= 1
        statusText = 'valid ';
    else
        statusText = 'test  ';
    end

    if isfinite(iter_ref(ia))

        fprintf(' %.2f       %s          %5d          %.6e      %.6e\n', ...
            alpha_set(ia),statusText,iter_ref(ia), ...
            final_ref_error(ia),final_Q(ia));

    else

        fprintf(' %.2f       %s        > maxIter      %.6e      %.6e\n', ...
            alpha_set(ia),statusText, ...
            final_ref_error(ia),final_Q(ia));

    end

end

validIndex = find(isfinite(iter_ref));

if ~isempty(validIndex)

    [fastestIterations,tempIndex] = min(iter_ref(validIndex));

    fastestIndex = validIndex(tempIndex);
    fastestAlpha = alpha_set(fastestIndex);

    fprintf('\n');
    fprintf('Fastest tested alpha = %.2f\n',fastestAlpha);
    fprintf('Fastest tested reference iterations = %d\n',fastestIterations);

    if fastestAlpha > 1
        fprintf(['WARNING: this fastest alpha is outside the ' ...
                 'reference paper''s admissible range.\n']);
    end

    if isfinite(iter_prop)
        fprintf(['Fastest tested reference / proposed ' ...
                 'iteration ratio = %.2fx\n'], ...
                 fastestIterations/iter_prop);
    end

end

indexAlpha1 = find(abs(alpha_set-1) < 1e-12,1);

if ~isempty(indexAlpha1)

    fprintf('\n');
    fprintf('------------------------------------------------------------\n');
    fprintf('REFERENCE-PAPER-COMPLIANT SETTING: alpha = 1\n');
    fprintf('------------------------------------------------------------\n');

    if isfinite(iter_ref(indexAlpha1))
        fprintf('alpha = 1 converged in %d iterations.\n', ...
            iter_ref(indexAlpha1));
    else
        fprintf('alpha = 1 did not reach tolerance.\n');
    end

    fprintf('Final alpha = 1 error = %.6e\n', ...
        final_ref_error(indexAlpha1));

    fprintf('Final alpha = 1 Q = %.6e\n', ...
        final_Q(indexAlpha1));

end

fprintf('\n');
fprintf('Expert-policy mismatch = %.3e\n',policyMismatch);
fprintf('============================================================\n');

fig1 = figure( ...
    'Color','w', ...
    'Position',[80 80 1000 600]);

ax1 = axes( ...
    'Parent',fig1, ...
    'Position',[0.10 0.13 0.86 0.80]);

hold(ax1,'on');
box(ax1,'on');
grid(ax1,'on');

plot( ...
    ax1, ...
    1:numel(err_prop), ...
    err_prop, ...
    '-o', ...
    'LineWidth',2.4, ...
    'MarkerSize',7, ...
    'DisplayName','Proposed Algorithm 3');

lineStyles = {'--','-.',':','--','-.','-'};

for ia = 1:nAlpha

    errorCurve = err_ref_all{ia};

    if alpha_set(ia) <= 1

        legendName = sprintf( ...
            'Reference, $\\alpha=%.2f$', ...
            alpha_set(ia));

    else

        legendName = sprintf( ...
            '$\\alpha=%.2f$ (diagnostic)', ...
            alpha_set(ia));

    end

    plot( ...
        ax1, ...
        1:numel(errorCurve), ...
        errorCurve, ...
        lineStyles{ia}, ...
        'LineWidth',1.35, ...
        'DisplayName',legendName);

end

yline( ...
    ax1, ...
    tol, ...
    ':', ...
    '$10^{-2}$', ...
    'Interpreter','latex', ...
    'LineWidth',1.2, ...
    'HandleVisibility','off');

set( ...
    ax1, ...
    'YScale','log', ...
    'FontName','Times New Roman', ...
    'FontSize',11, ...
    'LineWidth',1);

xlabel( ...
    ax1, ...
    'Iteration', ...
    'Interpreter','latex', ...
    'FontSize',13);

ylabel( ...
    ax1, ...
    '$\|\hat K-K^*\|$', ...
    'Interpreter','latex', ...
    'FontSize',13);

title( ...
    ax1, ...
    'Numerical sensitivity to $\alpha$', ...
    'Interpreter','latex', ...
    'FontSize',14);

legend( ...
    ax1, ...
    'Interpreter','latex', ...
    'Location','northeast', ...
    'FontSize',8);

referenceLengths = cellfun(@numel,err_ref_all);

xmax = max([numel(err_prop); referenceLengths(:)]);

xlim(ax1,[0 max(20,ceil(1.02*xmax))]);

ax2 = axes( ...
    'Parent',fig1, ...
    'Position',[0.18 0.55 0.31 0.29]);

hold(ax2,'on');
box(ax2,'on');
grid(ax2,'on');

nZoom = min(15,numel(err_prop));

plot( ...
    ax2, ...
    1:nZoom, ...
    err_prop(1:nZoom), ...
    '-o', ...
    'LineWidth',1.8, ...
    'MarkerSize',4);

for ia = 1:nAlpha

    errorCurve = err_ref_all{ia};

    nr = min(15,numel(errorCurve));

    plot( ...
        ax2, ...
        1:nr, ...
        errorCurve(1:nr), ...
        lineStyles{ia}, ...
        'LineWidth',1);

end

yline( ...
    ax2, ...
    tol, ...
    ':', ...
    'LineWidth',1);

set( ...
    ax2, ...
    'YScale','log', ...
    'FontName','Times New Roman', ...
    'FontSize',8, ...
    'LineWidth',0.8);

xlim(ax2,[1 15]);

xlabel(ax2,'Iteration','FontSize',8);
ylabel(ax2,'Gain error','FontSize',8);
title(ax2,'First 15 iterations','FontSize',9);

print( ...
    fig1, ...
    'comparison_B_higher_alpha', ...
    '-depsc', ...
    '-r300');

try

    exportgraphics( ...
        fig1, ...
        'comparison_B_higher_alpha.png', ...
        'Resolution',300);

catch
end

fig2 = figure( ...
    'Color','w', ...
    'Position',[150 120 720 480]);

bar(alpha_set,iter_ref);

grid on;

set( ...
    gca, ...
    'FontName','Times New Roman', ...
    'FontSize',11, ...
    'LineWidth',1);

xlabel( ...
    '$\alpha$', ...
    'Interpreter','latex', ...
    'FontSize',14);

ylabel( ...
    'Iterations to reach $\|\hat K-K^*\|\leq10^{-2}$', ...
    'Interpreter','latex', ...
    'FontSize',13);

title( ...
    'Numerical effect of $\alpha$ on convergence', ...
    'Interpreter','latex', ...
    'FontSize',14);

xline( ...
    1, ...
    ':', ...
    '$\alpha=1$', ...
    'Interpreter','latex', ...
    'LineWidth',1.2);

print( ...
    fig2, ...
    'reference_higher_alpha_iterations', ...
    '-depsc', ...
    '-r300');

try

    exportgraphics( ...
        fig2, ...
        'reference_higher_alpha_iterations.png', ...
        'Resolution',300);

catch
end

fprintf('\n');
fprintf('Saved comparison_B_higher_alpha.eps\n');
fprintf('Saved reference_higher_alpha_iterations.eps\n');


function [errK,converged] = proposed_algorithm3_same_policy( ...
    A,B,C,Kstar,K0,x0,tol)

n = size(A,1);
m = size(B,2);
p = size(C,1);

N = 3;

Obs = [C*A*A;
       C*A;
       C];

if rank(Obs) < n
    error('Algorithm 3 observability matrix is rank deficient.');
end

My = A^3/Obs;

Un = [B A*B A^2*B];

Tn = [0 C*B C*A*B;
      0 0   C*B;
      0 0   0];

Theta = [Un-My*Tn My];

Kbar_star = -Kstar*Theta;
Kbar = -K0*Theta;

SampleNumber = 200;
maxIter = 100;

errK = [];
converged = false;

fprintf('Theta size = %d x %d\n', ...
    size(Theta,1),size(Theta,2));

fprintf('rank(Theta) = %d / %d\n', ...
    rank(Theta),n);

fprintf('Initial ||K0-K*|| = %.6e\n', ...
    norm(K0-Kstar,'fro'));

for iter = 1:maxIter

    x = zeros(n,SampleNumber+1);
    xe = zeros(n,SampleNumber+1);

    y = zeros(p,SampleNumber);
    ye = zeros(p,SampleNumber);

    u = zeros(m,SampleNumber);
    ue = zeros(m,SampleNumber);

    z = zeros(m*N+p*N,SampleNumber);
    ze = zeros(m*N+p*N,SampleNumber);

    x(:,1) = x0;
    xe(:,1) = x0;

    for k = 1:SampleNumber

        y(:,k) = C*x(:,k);
        ye(:,k) = C*xe(:,k);

        if k >= N+1

            z(:,k) = [ ...
                u(:,k-1);
                u(:,k-2);
                u(:,k-3);
                y(:,k-1);
                y(:,k-2);
                y(:,k-3)];

            ze(:,k) = [ ...
                ue(:,k-1);
                ue(:,k-2);
                ue(:,k-3);
                ye(:,k-1);
                ye(:,k-2);
                ye(:,k-3)];

            excitation = 0.3*sum([ ...
                sin(k), ...
                sin(10*k), ...
                sin(0.1*k), ...
                sin(4.5*k), ...
                sin(6*k), ...
                sin(18*k), ...
                sin(0.01*k)]);

            u(:,k) = Kbar*z(:,k) + excitation;

            ue(:,k) = Kbar_star*ze(:,k) + excitation;

        end

        x(:,k+1) = A*x(:,k) + B*u(:,k);

        xe(:,k+1) = A*xe(:,k) + B*ue(:,k);

    end

    dz = m*N+p*N;
    dH = dz+m;

    P = sdpvar(dH,dH,'symmetric');
    Pe = sdpvar(dH,dH,'symmetric');

    Q = sdpvar(n,n,'symmetric');
    R = sdpvar(m,m,'symmetric');

    tau = sdpvar(1);

    Qbar = Theta'*Q*Theta;

    Lambda = blkdiag(Qbar,R);

    X1 = [];
    Y1 = [];

    Xe1 = [];
    Ye1 = [];

    BellL_scaled = [];
    BellE_scaled = [];

    for j = N+1:SampleNumber-1

        v_now = [ ...
            z(:,j);
            u(:,j)];

        ve_now = [ ...
            ze(:,j);
            ue(:,j)];

        z_next = z(:,j+1);
        ze_next = ze(:,j+1);

        v_next = [ ...
            z_next;
            Kbar*z_next];

        ve_next = [ ...
            ze_next;
            Kbar_star*ze_next];

        XL = v_now'*P*v_now - v_next'*P*v_next;

        YL = v_now'*Lambda*v_now;

        X1 = [X1; XL];
        Y1 = [Y1; YL];

        XE = ve_now'*Pe*ve_now - ve_next'*Pe*ve_next;

        YE = ve_now'*Lambda*ve_now;

        Xe1 = [Xe1; XE];
        Ye1 = [Ye1; YE];

        scaleL = max([ ...
            norm(v_now)^2, ...
            norm(v_next)^2, ...
            1e-3]);

        scaleE = max([ ...
            norm(ve_now)^2, ...
            norm(ve_next)^2, ...
            1e-3]);

        BellL_scaled = [ ...
            BellL_scaled;
            (XL-YL)/scaleL];

        BellE_scaled = [ ...
            BellE_scaled;
            (XE-YE)/scaleE];

    end

    objective = tau^2;

    structuralConstraints = [ ...
        blkdiag(Q,R) >= eye(n+m), ...
        blkdiag(Q,R) <= tau*eye(n+m), ...
        tau >= 1, ...
        tau <= 100];

    constraintsOriginal = [ ...
        P >= 0, ...
        Pe >= 0, ...
        Q >= 0, ...
        P-Pe >= 0, ...
        X1-Y1 == 0, ...
        Xe1-Ye1 == 0, ...
        structuralConstraints];

    options = sdpsettings( ...
        'solver','mosek', ...
        'verbose',0);

    solution = optimize( ...
        constraintsOriginal, ...
        objective, ...
        options);

    if solution.problem ~= 0

        fprintf(['Original SDP solver code %d ' ...
                 'at iteration %d.\n'], ...
                 solution.problem,iter);

        fprintf('Retrying with scaled exact Bellman equations...\n');

        constraintsScaled = [ ...
            P >= 0, ...
            Pe >= 0, ...
            Q >= 0, ...
            P-Pe >= 0, ...
            BellL_scaled == 0, ...
            BellE_scaled == 0, ...
            structuralConstraints];

        solution = optimize( ...
            constraintsScaled, ...
            objective, ...
            options);

    end

    if solution.problem ~= 0

        fprintf('Solver code = %d at iteration %d\n', ...
            solution.problem,iter);

        fprintf('YALMIP information: %s\n', ...
            solution.info);

        break;

    end

    Hvalue = value(P);

    Hzu = Hvalue(1:dz,end);
    Huu = Hvalue(end,end);

    if ~isfinite(Huu) || abs(Huu) < 1e-10

        warning( ...
            'Invalid Huu at Algorithm 3 iteration %d.', ...
            iter);

        break;

    end

    Kbar_new = -(Hzu')/Huu;

    Khat = -Kbar_new*pinv(Theta);

    gainError = norm(Khat-Kstar,'fro');

    historyError = norm( ...
        Kbar_new-Kbar_star, ...
        'fro');

    policyStep = norm( ...
        Kbar_new-Kbar, ...
        'fro');

    errK(end+1,1) = gainError;

    fprintf(['Iter %3d | ' ...
             '||K-K*|| = %.6e | ' ...
             '||Kbar-Kbar*|| = %.6e | ' ...
             'step = %.6e | ' ...
             'tau = %.6f\n'], ...
             iter, ...
             gainError, ...
             historyError, ...
             policyStep, ...
             value(tau));

    Kbar = Kbar_new;

    if gainError <= tol

        fprintf(['Proposed Algorithm 3 reached ' ...
                 'tolerance in %d iterations.\n'], ...
                 iter);

        converged = true;

        break;

    end

end

if isempty(errK)
    error('Proposed Algorithm 3 produced no feasible iterate.');
end

end


function [errK,converged,Qfinal] = reference_inverseQ_modified( ...
    A,B,C,Kstar,Ky_star,x0,tol,alpha)

N = 3;

Q = 0;
Re = 2;

T = 300;

x = zeros(3,T+1);
y = zeros(1,T);
u = zeros(1,T);
omega = zeros(1,T);

x(:,1) = x0;

for k = 1:T

    y(k) = C*x(:,k);

    omega(k) = 0.01*( ...
        sin(0.31*k) ...
        + 0.70*sin(1.10*k) ...
        + 0.50*sin(2.30*k) ...
        + 0.25*sin(3.70*k) ...
        + 0.20*sin(5.20*k));

    u_clean = -Ky_star*y(k);

    u(k) = u_clean + omega(k);

    x(:,k+1) = A*x(:,k) + B*u(k);

end

fprintf('alpha = %.2f | expert mismatch = %.3e\n', ...
    alpha,norm(Kstar-Ky_star*C,'fro'));

M = T-N-1;

PhiH = zeros(M,28);

Uhist = zeros(3,M);
Yhist = zeros(3,M);
Ohist = zeros(3,M);

ycur = zeros(1,M);
ucur = zeros(1,M);
ocur = zeros(1,M);

row = 0;

for k = N+1:T-1

    row = row+1;

    U = [ ...
        u(k-1);
        u(k-2);
        u(k-3)];

    Y = [ ...
        y(k-1);
        y(k-2);
        y(k-3)];

    O = [ ...
        omega(k-1);
        omega(k-2);
        omega(k-3)];

    zk = [ ...
        U;
        Y;
        u(k)];

    uTargetNext = u(k+1)-omega(k+1);

    zk1 = [ ...
        u(k);
        u(k-1);
        u(k-2);
        y(k);
        y(k-1);
        y(k-2);
        uTargetNext];

    PhiH(row,:) = (svecq(zk)-svecq(zk1))';

    Uhist(:,row) = U;
    Yhist(:,row) = Y;
    Ohist(:,row) = O;

    ycur(row) = y(k);
    ucur(row) = u(k);
    ocur(row) = omega(k);

end

rankH = rank(PhiH);

fprintf('H rank = %d / %d\n', ...
    rankH,size(PhiH,2));

if rankH < size(PhiH,2)
    error('Reference H regression is rank deficient.');
end

columnScale = sqrt(sum(PhiH.^2,1));

columnScale(columnScale < 1e-12) = 1;

PhiScaled = PhiH ./ columnScale;

fprintf('cond(PhiH) = %.3e\n',cond(PhiH));

fprintf('cond(scaled PhiH) = %.3e\n', ...
    cond(PhiScaled));

Lx = zeros(1,M);

maxIter = 5000;

errK = [];

converged = false;

Qfinal = Q;

for iter = 1:maxIter

    beta = Q*(ycur.^2) + Re*(ucur.^2);

    hScaled = PhiScaled \ beta';

    hvec = hScaled ./ columnScale';

    Hbar = sym_from_vem(hvec,7);

    HuU = Hbar(7,1:3);
    HuY = Hbar(7,4:6);
    Huu = Hbar(7,7);

    if ~isfinite(Huu) || abs(Huu) < 1e-12

        warning( ...
            'Reference Huu invalid at iteration %d.', ...
            iter);

        break;

    end

    base = HuU*(Uhist-Ohist) + HuY*Yhist;

    numerator = ycur*(base+Lx)';

    denominator = Huu*(ycur*ycur');

    if ~isfinite(denominator) || abs(denominator) < 1e-14

        warning('Reference policy denominator invalid.');

        break;

    end

    Ky = numerator/denominator;

    Lx_new = -base + Huu*Ky*ycur;

    uClean = ucur-ocur;

    policyResidual = Ky*ycur + uClean;

    rhsQ = alpha*Huu*(policyResidual.^2) ...
         + Q*(ycur.^2);

    phiQ = ycur.^2;

    qDenominator = phiQ*phiQ';

    if ~isfinite(qDenominator) || abs(qDenominator) < 1e-14

        warning('Reference Q denominator invalid.');

        break;

    end

    Qnew = (phiQ*rhsQ')/qDenominator;

    Khat = Ky*C;

    gainError = norm(Khat-Kstar,'fro');

    errK(end+1,1) = gainError;

    if iter <= 10 || mod(iter,100) == 0

        fprintf(['alpha = %.2f | ' ...
                 'Iter %4d | ' ...
                 'Ky = %.8f | ' ...
                 'Q = %.8f | ' ...
                 'error = %.6e | ' ...
                 'Huu = %.8f\n'], ...
                 alpha, ...
                 iter, ...
                 Ky, ...
                 Qnew, ...
                 gainError, ...
                 Huu);

    end

    Q = Qnew;

    Lx = Lx_new;

    Qfinal = Q;

    if gainError <= tol

        fprintf(['alpha = %.2f reached tolerance ' ...
                 'at iteration %d.\n'], ...
                 alpha,iter);

        converged = true;

        break;

    end

    if ~isfinite(Q) || ...
       ~isfinite(Ky) || ...
       abs(Q) > 1e12 || ...
       abs(Ky) > 1e6

        warning(['Reference method became numerically unstable ' ...
                 'for alpha = %.2f at iteration %d.'], ...
                 alpha,iter);

        break;

    end

end

if isempty(errK)
    error('Reference method produced no iterations.');
end

end


function k = first_hit(err,tol)

idx = find(err <= tol,1,'first');

if isempty(idx)
    k = NaN;
else
    k = idx;
end

end


function phi = svecq(v)

v = v(:);

n = numel(v);

phi = zeros(n*(n+1)/2,1);

counter = 1;

for i = 1:n

    phi(counter) = v(i)^2;

    counter = counter+1;

    for j = i+1:n

        phi(counter) = 2*v(i)*v(j);

        counter = counter+1;

    end

end

end


function M = sym_from_vem(v,n)

v = v(:);

M = zeros(n);

counter = 1;

for i = 1:n

    M(i,i) = v(counter);

    counter = counter+1;

    for j = i+1:n

        M(i,j) = v(counter);

        M(j,i) = v(counter);

        counter = counter+1;

    end

end

end