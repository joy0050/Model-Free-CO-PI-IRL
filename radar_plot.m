%% ============================================================
%  Explainability Radar — 3-axis, circular grid, case colours
%  Saves: explainability_spider.eps
%% ============================================================
close all; clc; clear;

%% ============================================================
%%  CONFIG
%% ============================================================
FONT      = 'Times New Roman';
FS_LABEL  = 13;
FS_RING   = 13;
FS_LEG    = 13;

FIG_W  = 3.4;
FIG_H  = 3.9;
AX_POS = [0.02 0.02 0.96 0.96];
R_LBL  = 1.14;

NUDGE = [
    0.00,  -0.08;
    -0.08,  0.00;
   0.08,  0.00;
];

CC = {
    [0.18 0.45 0.72];
    [0.85 0.38 0.10];
    [0.18 0.62 0.33];
};
COL2 = [0.17 0.43 0.70];
COL3 = [0.80 0.22 0.17];

%% ============================================================
%%  DATA
%% ============================================================
kap2 = 7.665;  kap3 = 9.865;
nq2  = 7.556;  nq3  = 10.283;
rho2 = 5.094;  rho3 = 6.823;

sc2 = [(12-kap2)/6,  (12-nq2)/6,  rho2/10];
sc3 = [(12-kap3)/6,  (12-nq3)/6,  rho3/10];

LABELS = {
    {'\bf Cost Balance',          'Case 1 \ \ $\kappa$'             };
    {'\bf Explanation Simplicity','Case 2 \ \ $\|\hat{Q}_\tau\|_F$' };
    {'\bf Expert Priority',       'Case 3 \ \ $\rho$'               };
};

N = 3;

%% ============================================================
%%  FIGURE
%% ============================================================
fig = figure('Color','w');
fig.Units    = 'inches';
fig.Position = [1 1 FIG_W FIG_H];
ax = axes('Position', AX_POS);
hold on; axis equal off;
xlim([-1.85 1.85]);
ylim([-1.60 1.95]);

theta = pi/2 - (0:N-1)*(2*pi/N);

%% --- 1. Circular rings (bottom layer) ---
th_c = linspace(0, 2*pi, 300);
for r = [0.25 0.50 0.75 1.00]
    if r == 1.0
        plot(r*cos(th_c), r*sin(th_c), '-', ...
            'Color',[0.55 0.55 0.55],'LineWidth',0.9)
    else
        plot(r*cos(th_c), r*sin(th_c), ':', ...
            'Color',[0.78 0.78 0.78],'LineWidth',0.6)
    end
end

%% --- 2. Coloured spokes ---
for k = 1:N
    plot([0 cos(theta(k))],[0 sin(theta(k))],'-', ...
        'Color',[0.80 0.80 0.80],'LineWidth',0.7)
    plot([0 cos(theta(k))],[0 sin(theta(k))],'-', ...
        'Color',[CC{k} 0.55],'LineWidth',2.2)
    plot(cos(theta(k)), sin(theta(k)), 'o', ...
        'MarkerSize',5,'MarkerFaceColor',CC{k}, ...
        'MarkerEdgeColor','w','LineWidth',0.8)
end

%% --- 3. Filled polygons ---
tc = [theta theta(1)];
s2 = [sc2 sc2(1)];   s3 = [sc3 sc3(1)];
x2 = s2.*cos(tc);    y2 = s2.*sin(tc);
x3 = s3.*cos(tc);    y3 = s3.*sin(tc);

fill(x3,y3,COL3,'FaceAlpha',0.14,'EdgeColor','none')
fill(x2,y2,COL2,'FaceAlpha',0.14,'EdgeColor','none')
plot(x3,y3,'-','Color',COL3,'LineWidth',2.0)
plot(x2,y2,'-','Color',COL2,'LineWidth',2.0)
plot(x2(1:N),y2(1:N),'o','MarkerSize',7.0, ...
    'MarkerFaceColor',COL2,'MarkerEdgeColor','w','LineWidth',1.1)
plot(x3(1:N),y3(1:N),'s','MarkerSize',6.5, ...
    'MarkerFaceColor',COL3,'MarkerEdgeColor','w','LineWidth',1.1)

%% --- 4. Ring scale labels AFTER polygons so they sit on top ---
for r = [0.25 0.50 0.75]
    text(r*cos(theta(1))+0.05, r*sin(theta(1)), ...
        sprintf('%.2g', r), ...
        'HorizontalAlignment','left','VerticalAlignment','middle', ...
        'FontName',FONT,'FontSize',FS_RING,'Color',[0.35 0.35 0.35])
end

%% --- 5. Axis labels ---
for k = 1:N
    cx = cos(theta(k));   cy = sin(theta(k));
    dx = NUDGE(k,1);      dy = NUDGE(k,2);

    if cx >  0.3;  ha='left';   elseif cx < -0.3;  ha='right';  else; ha='center'; end
    if cy >  0.3;  va='bottom'; elseif cy < -0.3;  va='top';    else; va='middle'; end

    text(R_LBL*cx+dx, R_LBL*cy+dy, LABELS{k}, ...
        'HorizontalAlignment',ha,'VerticalAlignment',va, ...
        'FontName',FONT,'FontSize',FS_LABEL, ...
        'Interpreter','latex','Color',CC{k}*0.72)
end

%% --- 6. Legend (vertical, right side) ---
h1 = plot(nan,nan,'o-','Color',COL2,'MarkerFaceColor',COL2, ...
    'MarkerSize',6,'LineWidth',2.0,'DisplayName','Algorithm~2');
h2 = plot(nan,nan,'s-','Color',COL3,'MarkerFaceColor',COL3, ...
    'MarkerSize',5.5,'LineWidth',2.0,'DisplayName','Algorithm~3');
lg = legend([h1 h2],'Interpreter','latex', ...
    'FontName',FONT,'FontSize',FS_LEG,'NumColumns',1,'Box','off');
%% Tweak lg.Position = [left bottom width height] (normalised figure units)
%% to reposition the legend — increase 'left' to move right, adjust 'bottom'
lg.Position = [0.76 0.42 0.50 0.14];

%% ============================================================
%%  EXPORT
%% ============================================================
exportgraphics(fig,'explainability_spider.eps', ...
    'ContentType','vector','BackgroundColor','white');
fprintf('Saved  %.1f x %.1f in\n', FIG_W, FIG_H);