close all;
clc;
clear;

rng(2);

A = [0.883878 0        0.047013;
     0.089407 0.904838 0.002321;
     0.056415 0        1.000939];

B = [0.11753;
     0.00559;
     0.00359];

Qe = eye(3);
Re = 5;
tol = 1e-2;

[~,~,Kstar] = dare(A,B,Qe,Re);

fprintf('Expert gain K* = [%.4f %.4f %.4f]\n',Kstar);
fprintf('Comparison threshold = %.1e\n\n',tol);

K0 = [0.4910 0.2849 0.7923];

rho0 = max(abs(eig(A-B*K0)));

fprintf('rho(A-B*K0) = %.6f\n\n',rho0);

fprintf('============================================================\n');
fprintf('Proposed Algorithm 2\n');
fprintf('============================================================\n');

err_prop = proposed_alg2(A,B,Kstar,K0,tol);

alpha_set = [0.2 0.4 0.6 0.8 1.0];

err_base = cell(numel(alpha_set),1);

for ia = 1:numel(alpha_set)

    fprintf('\n============================================================\n');
    fprintf('Off-policy inverse Q-learning: alpha = %.1f\n',alpha_set(ia));
    fprintf('============================================================\n');

    err_base{ia} = offpolicy_inverseQ( ...
        A,B,Kstar,K0,Re,alpha_set(ia),tol);

end

fig = figure( ...
    'Color','w', ...
    'Position',[100 100 820 520]);

ax1 = axes( ...
    'Parent',fig, ...
    'Position',[0.11 0.13 0.83 0.79]);

hold(ax1,'on');
box(ax1,'on');
grid(ax1,'on');

plot( ...
    ax1, ...
    1:numel(err_prop), ...
    err_prop, ...
    '-o', ...
    'LineWidth',2.3, ...
    'MarkerSize',7, ...
    'DisplayName','Proposed Algorithm 2');

lineStyles = {'--','-.',':','--','-'};

for ia = 1:numel(alpha_set)

    e = err_base{ia};

    plot( ...
        ax1, ...
        1:numel(e), ...
        e, ...
        lineStyles{ia}, ...
        'LineWidth',1.4, ...
        'DisplayName', ...
        sprintf('Off-policy Q-learning, $\\alpha=%.1f$',alpha_set(ia)));

end

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
    '$\|\hat{K}_{\tau}-K^*\|$', ...
    'Interpreter','latex', ...
    'FontSize',13);

title( ...
    ax1, ...
    'Algorithm 2 vs Off policy', ...
    'Interpreter','latex', ...
    'FontSize',14);

legend( ...
    ax1, ...
    'Interpreter','latex', ...
    'Location','northeast', ...
    'FontSize',9);

xlim(ax1,[0 1100]);

ax2 = axes( ...
    'Parent',fig, ...
    'Position',[0.20 0.55 0.30 0.29]);

hold(ax2,'on');
box(ax2,'on');
grid(ax2,'on');

plot( ...
    ax2, ...
    1:numel(err_prop), ...
    err_prop, ...
    '-o', ...
    'LineWidth',2.0, ...
    'MarkerSize',5);

for ia = 1:numel(alpha_set)

    e = err_base{ia};

    nZoom = min(10,numel(e));

    plot( ...
        ax2, ...
        1:nZoom, ...
        e(1:nZoom), ...
        lineStyles{ia}, ...
        'LineWidth',1.1);

end

set( ...
    ax2, ...
    'YScale','log', ...
    'FontName','Times New Roman', ...
    'FontSize',8, ...
    'LineWidth',0.8);

xlim(ax2,[1 10]);

xlabel( ...
    ax2, ...
    'Iteration', ...
    'FontName','Times New Roman', ...
    'FontSize',8);

ylabel( ...
    ax2, ...
    '$\|\hat{K}-K^*\|$', ...
    'Interpreter','latex', ...
    'FontSize',8);

title( ...
    ax2, ...
    'First 10 iterations', ...
    'FontName','Times New Roman', ...
    'FontSize',9);

print( ...
    fig, ...
    'comparison_A_inset', ...
    '-depsc', ...
    '-r300');

fprintf('\nSaved comparison_A_inset.eps\n');


function errK = proposed_alg2(A,B,Kstar,K0,tol)

n = size(A,1);
m = size(B,2);

K = K0;

SampleNumber = 40;
maxIter = 30;

x0 = [0.1;
      0.5;
      0.5];

errK = [];

for iter = 1:maxIter

    xt = zeros(n,SampleNumber+1);
    xet = zeros(n,SampleNumber+1);

    u = zeros(m,SampleNumber);
    ue = zeros(m,SampleNumber);

    xt(:,1) = x0;
    xet(:,1) = x0;

    for k = 1:SampleNumber

        eL = 0.1*( ...
            8.5*sin(1.5*k) + ...
            5.0*sin(0.1*k) + ...
            1.2*sin(0.5*k) + ...
            2.0*sin(2.6*k));

        eE = 0.1*( ...
            8.5*sin(1.1*k) + ...
            5.0*sin(0.1*k) + ...
            1.0*sin(0.5*k) + ...
            2.0*sin(2.6*k));

        u(:,k) = -K*xt(:,k) + eL;

        ue(:,k) = -Kstar*xet(:,k) + eE;

        xt(:,k+1) = A*xt(:,k) + B*u(:,k);

        xet(:,k+1) = A*xet(:,k) + B*ue(:,k);

    end

    H = sdpvar(n+m,n+m,'symmetric');

    He = sdpvar(n+m,n+m,'symmetric');

    Q = sdpvar(n,n,'symmetric');

    R = sdpvar(m,m,'symmetric');

    eta = sdpvar(1);

    Lambda = blkdiag(Q,R);

    bellL = [];
    bellE = [];

    for k = 1:SampleNumber-1

        xi = [ ...
            xt(:,k);
            u(:,k)];

        xi_next = [ ...
            xt(:,k+1);
            -K*xt(:,k+1)];

        xie = [ ...
            xet(:,k);
            ue(:,k)];

        xie_next = [ ...
            xet(:,k+1);
            -Kstar*xet(:,k+1)];

        scaleL = max([ ...
            norm(xi)^2, ...
            norm(xi_next)^2, ...
            1e-3]);

        scaleE = max([ ...
            norm(xie)^2, ...
            norm(xie_next)^2, ...
            1e-3]);

        bellL = [ ...
            bellL;
            (xi'*H*xi ...
            - xi_next'*H*xi_next ...
            - xi'*Lambda*xi)/scaleL];

        bellE = [ ...
            bellE;
            (xie'*He*xie ...
            - xie_next'*He*xie_next ...
            - xie'*Lambda*xie)/scaleE];

    end

    objective = eta^2;

    bellTol = 1e-7;

    constraints = [ ...
        H >= 0, ...
        He >= 0, ...
        H-He >= 0, ...
        -bellTol <= bellL <= bellTol, ...
        -bellTol <= bellE <= bellTol, ...
        eye(n+m) <= Lambda, ...
        Lambda <= eta*eye(n+m), ...
        eta >= 1, ...
        eta <= 100];

    options = sdpsettings( ...
        'solver','mosek', ...
        'verbose',0, ...
        'savesolverinput',1, ...
        'savesolveroutput',1);

    solution = optimize( ...
        constraints, ...
        objective, ...
        options);

    if solution.problem ~= 0

        fprintf( ...
            'Solver code = %d at iteration %d\n', ...
            solution.problem, ...
            iter);

        fprintf( ...
            'YALMIP information: %s\n', ...
            solution.info);

        warning('Stopping proposed method because SDP was not solved.');

        break;

    end

    Hvalue = value(H);

    Hux = Hvalue(n+1:end,1:n);

    Huu = Hvalue(n+1:end,n+1:end);

    Knew = Huu\Hux;

    gainError = norm(Knew-Kstar,'fro');

    policyStep = norm(Knew-K,'fro');

    errK(end+1,1) = gainError;

    fprintf( ...
        ['Iter %2d | ' ...
         '||K-K*|| = %.6e | ' ...
         '||Knew-Kold|| = %.6e | ' ...
         'eta = %.6f\n'], ...
         iter, ...
         gainError, ...
         policyStep, ...
         value(eta));

    K = Knew;

    if policyStep <= tol

        fprintf( ...
            'Proposed Algorithm 2 converged in %d iterations.\n', ...
            iter);

        break;

    end

end

end


function errK = offpolicy_inverseQ(A,B,Ke,K0,R,alpha,tol)

n = size(A,1);

K = K0;
Q = zeros(n);

maxIter = 1100;

nTraj = 16;
stepsPerTraj = 20;

errK = [];

for iter = 1:maxIter

    Phi = [];
    rhsH = [];

    PhiQ = [];
    rhsQ = [];

    for trajectory = 1:nTraj

        x = [ ...
            1 + 0.08*trajectory;
           -1 + 0.04*trajectory;
            0.4 - 0.015*trajectory];

        for k = 1:stepsPerTraj

            kk = (iter-1)*nTraj*stepsPerTraj ...
               + (trajectory-1)*stepsPerTraj ...
               + k;

            d = 0.20*( ...
                sin(0.31*kk) + ...
                0.6*sin(1.70*kk) + ...
                0.3*sin(2.30*kk));

            u = -Ke*x + d;

            xNext = A*x + B*u;

            Kx = K*x;
            KxNext = K*xNext;

            regressionH = [ ...
                (svecq(x)-svecq(xNext))', ...
                2*kron(xNext',KxNext') ...
                + 2*kron(x',u'), ...
                u^2-KxNext^2];

            targetH = ...
                x'*(K'*R*K+Q)*x ...
                + (u-Kx)'*R*(u+Kx) ...
                + alpha*(u+Kx-d)'*R*(u+Kx-d);

            Phi = [Phi; regressionH];

            rhsH = [rhsH; targetH];

            PhiQ = [PhiQ; svecq(x)'];

            rhsQ = [ ...
                rhsQ;
                x'*Q*x ...
                + alpha*(u+Kx-d)'*R*(u+Kx-d)];

            x = xNext;

        end

    end

    if rank(Phi) < size(Phi,2)

        warning( ...
            'Off-policy regression lost rank at iteration %d.', ...
            iter);

        break;

    end

    h = Phi\rhsH;

    Hxu = h(7:9);
    Huu = h(10);

    Knew = Hxu.'/Huu;

    qVector = PhiQ\rhsQ;

    Qnew = sym_from_vem(qVector,n);

    gainError = norm(Knew-Ke,'fro');

    policyStep = norm(Knew-K,'fro');

    errK(end+1,1) = gainError;

    K = Knew;
    Q = Qnew;

    if iter == 1 || mod(iter,25) == 0

        fprintf( ...
            ['Iter %4d | ' ...
             '||K-K*|| = %.6e | ' ...
             '||Knew-Kold|| = %.6e\n'], ...
             iter, ...
             gainError, ...
             policyStep);

    end

    if gainError <= tol

        fprintf( ...
            ['Baseline alpha=%.1f reached ' ...
             '||K-K*|| <= %.1e in %d iterations.\n'], ...
             alpha, ...
             tol, ...
             iter);

        break;

    end

end

end


function phi = svecq(v)

v = v(:);

n = numel(v);

numberTerms = n*(n+1)/2;

phi = zeros(numberTerms,1);

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