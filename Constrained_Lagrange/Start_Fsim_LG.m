% Constrained Lagrange Approach

clear all;
close all;
clc;

%% Build Lagrange dynamics

dim = 2;
e = eye(dim);

M = 5*e;
UM = 1*e;
I = 10;
g = [0 -9.8]';
k = 5000;
b = 5000;
s0 = 0.3;

syms t
% theta = sym('theta',[1,3]);
CG = sym('CG',[2,3]);
LW = sym('LW',[2,3]);
RW = sym('RW',[2,3]);
LA = sym('LA',[2,3]);
RA = sym('RA',[2,3]);

q = [CG; LW; RW; LA; RA]; %;theta;

W = blkdiag(M, UM, UM, UM, UM); %,I
G = [g; g; g; g; g]; %0; 
K = diag([k,k]);

dx = [normv(LA(:,1)-LW(:,1)) - s0;
      normv(RA(:,1)-RW(:,1)) - s0];
dxp =[LA(:,2)-LW(:,2);
      RA(:,2)-RW(:,2)];

T = 0.5 * q(:,2).' * W * q(:,2);
V = - G.' * W * q(:,1) + 0.5 * dx.' * K * dx;
D = 0.5 * dxp.' * b * dxp;

% acosd( dot(u,v)/(normv(u)*normv(v)) )
C = [ normv(LA(:,1)-CG(:,1))^2 - 0.5^2;
      normv(CG(:,1)-RA(:,1))^2 - 0.5^2;
      normv(RA(:,1)-LA(:,1))^2 - 0.8^2;
      normv(CG(:,1)-LW(:,1))^2 - 0.5^2;
      normv(CG(:,1)-RW(:,1))^2 - 0.5^2] /2;
dC = jacobian(C,q(:,1))*q(:,2); % == 0
ddC = jacobian(dC,q(:,1))*q(:,2) + jacobian(dC,q(:,2))*q(:,3); % == 0

z = sym('z',size(C));

Fe = sym('Fe',[size(q,1),1]);
Qnc = Fe - gradient(D,q(:,1));

Lag = T - V - z.'*C;
Lag_dq = gradient(Lag,q(:,2));
LAGEQ = jacobian(Lag_dq,q(:,1))*q(:,2) + jacobian(Lag_dq,q(:,2))*q(:,3) - gradient(Lag,q(:,1)) - Fe;

wn=500;
Fdae = [LAGEQ; ddC+2*wn*dC+wn^2*C];

% Monitor consistency conditions
C_fun  = matlabFunction(C,'Vars',  {q(:,1),q(:,2)});
dC_fun = matlabFunction(dC,'Vars', {q(:,1),q(:,2)});


%% Prepare simulation


qt = sym2symfunt(q);
zt = sym2symfunt(z);

Fdae_t = subs(Fdae, [z(:,1);q(:,1);q(:,2);q(:,3)], [zt;qt;diff(qt);diff(qt,2)]);
[Fdae_red, q_red, Rel] = reduceDifferentialOrder(Fdae_t, [qt;zt]);

if ~isLowIndexDAE(Fdae_red, q_red)
    error('DAE Index > 1 !!!')
end

[M_ode,b_ode] = massMatrixForm(Fdae_red, q_red);
M_ode_fun = odeFunction(M_ode, q_red, Fe);
b_ode_fun = odeFunction(b_ode, q_red, Fe);
M_ode_fun_k = @(t,q_red) M_ode_fun(t, q_red, inputs_k(q_red));
b_ode_fun_k = @(t,q_red) b_ode_fun(t, q_red, inputs_k(q_red));          
                


%% Simulate

tspan = [0,2];
cgy = 0;
cgz = 1;
q0  = [cgy cgz cgy-0.5 cgz cgy+0.5 cgz cgy-0.4 cgz+0.3 cgy+0.4 cgz+0.3]';
dq0 = 0*q0;

y0 = [q0; zeros(size(z)); 0*q0];
dy0 = 0*y0;
  
consistencyConditions = [ C_fun(q0, dq0); 
                         dC_fun(q0, dq0) ];
if any(abs(consistencyConditions) > 1e-5)
    warning('Initial conditions not consistent!!!')
end


opt = odeset('Mass', M_ode_fun_k, 'InitialSlope', dy0, ...
             'MassSingular','yes', ...
             'OutputFcn','odeplot', 'OutputSel', 1:10,...
             'RelTol', 1e-4, 'AbsTol' ,1e-5,'Stats','on',...
             'Events',@events);
  
         
% BUG OCCURS WHEN BOTH EVENTS ARE CALLED AT THE SAME TIME.
% NEED TO DECIDE WHAT TO DO IN THOSE CASES.
% PROBABLY RESTART FROM THE STEP WITH SMALLEST TIME
% ERASE ALL STEPS BEFORE THE FIRST EVENT
% 
% Need to check why sometimes the masses are going below the ground
tsol = [];
ysol = [];
tstart = tspan(1);
tend = tspan(2);
while tstart < tend
    [ti, yi, te,ye,ie] = ode23t(b_ode_fun_k, [tstart tend], y0, opt);
    if ti(end) == tend
        break;
    end
    
    tvalid = ti<=min(te);
    ti = ti(tvalid);
    yi = yi(tvalid,:);
    y0 = ye(1,:);
    if ie(1)==1 % ye(:,4) ye(:,19)
        y0(19) = -0.8*y0(19);
    end
    if ie(1)==2 % ye(:,6) ye(:,21)
        y0(21) = -0.8*y0(21);
    end
    tstart = ti(end);
    tsol = [tsol;ti(2:end)];
    ysol = [ysol;yi(2:end,:)];
end


for i=1:size(ysol,1)
    isol(i,:)  = inputs_k(ysol(i,:)'); % inputs(t)
    csol(i,:)  = C_fun(ysol(i,1:10)', ysol(i,end-10:end)');
    dcsol(i,:) = dC_fun(ysol(i,1:10)', ysol(i,end-10:end)');
end

%% Analysis

% plot forces
figure('Color','white'); hold on; grid on;
plot(tsol, isol(:,4));
plot(tsol, isol(:,6));
ylabel('External forces')

% monitor consistency conditions
figure('Color','white');
subplot(2,1,1);
plot(tsol, csol); grid on;
ylabel('C(q)')
subplot(2,1,2);
plot(tsol, dcsol); grid on;
ylabel('dC(q)')

%% Plot

fps = 60;
tsol_plot = linspace(min(tsol), max(tsol), floor(fps*max(tsol)) )';
ysol_plot = interp1(tsol,ysol,tsol_plot,'pchip');

figure('Color','white'); hold on; grid on;
axis equal
xlim([-1 1])
ylim([-3 3])
%     scatter(ysol(i,[1,3,5,7,9]), ysol(i,[2,4,6,8,10]), 50, 'filled');
pw = plot( ysol_plot(1,[3,1,5]), ysol_plot(1,[4,2,6]),'LineWidth',3);
pb = plot( ysol_plot(1,[1,7,9,1]), ysol_plot(1,[2,8,10,2]),'LineWidth',3);
wl = plot( [ysol_plot(1,3) ysol_plot(1,3)], [0 ysol_plot(1,4)]  );
wr = plot( [ysol_plot(1,5) ysol_plot(1,5)], [0 ysol_plot(1,6)]  );
sl = plot( ysol_plot(1,[3,7]), ysol_plot(1,[4,8]),'LineWidth',3);
sr = plot( ysol_plot(1,[5,9]), ysol_plot(1,[6,10]),'LineWidth',3);

realtime = 0.1;
tic
for i=1:numel(tsol_plot)
    pw.XData = ysol_plot(i,[3,1,5]);
    pw.YData = ysol_plot(i,[4,2,6]);
    pb.XData = ysol_plot(i,[1,7,9,1]);
    pb.YData = ysol_plot(i,[2,8,10,2]);
    wl.XData = [ysol_plot(i,3) ysol_plot(i,3)];
    wl.YData = [0 ysol_plot(i,4)];
    wr.XData = [ysol_plot(i,5) ysol_plot(i,5)];
    wr.YData = [0 ysol_plot(i,6)];
    sl.XData = ysol_plot(i,[3,7]);
    sl.YData = ysol_plot(i,[4,8]);
    sr.XData = ysol_plot(i,[5,9]);
    sr.YData = ysol_plot(i,[6,10]);
    drawnow();
    pause( tsol_plot(i) - toc*realtime )
end





%% Help functions

function [value,isterminal,direction] = events(t,q)
    value   = [q(4);q(6)];
    direction = [-1;-1];
    isterminal = [1;1];
end

function qt = sym2symfunt(q)
    varsqt='';
    for i=1:size(q,1)
        eval(sprintf('syms %s(t)',char(q(i,1))));
        varsqt = [varsqt,' ',char(q(i,1))];             %#ok<*AGROW>
    end
    eval(['qt=[',varsqt,'];']);
    qt = qt.';
end

function Qnc = inputs_k(q_red)
    S0 = 0.1;
    k = 500;
    c = 10;
%     if q_red(4)<=S0
%         FL = -k*(q_red(4)-S0) - c*(q_red(18));
%     else
        FL = 0;
%     end
%     if q_red(6)<=S0
%         FR = k*(S0 - q_red(6)) - c*(q_red(20));
%     else
        FR = 0;
%     end
    Qnc = [0 0 0 FL 0 FR 0 0 0 0].';
end

function n = normv(v)
    n = (v.'*v)^.5;
end


%% arch

% SET MANUALLY DAE FUNCTION
% qp = sym('qp',[size(q,1),1]);
% Fdae_red = [qp-q(:,2) ;Fdae];
% q_red  = [q(:,1); q(:,2); z];
% dq_red = [qp;     q(:,3)];
% det(jacobian(F,[qp;q(:,3);z]))    % if ~=0 then is index-1 DAE
% F = matlabFunction(Fdae_red, 'Vars', {t, q_red, dq_red, Qnc});
% F15i = @(t, q_red, dq_red)  F(t, q_red, dq_red, inputs_k(q_red));

% [newEqs,constraintEqs,oldIndex] = reduceDAEToODE(Fdae_red, q_red)
% [DAEs,DAEvars] = reduceDAEIndex(Fdae_red, q_red);
% [DAEs,DAEvars] = reduceRedundancies(DAEs,DAEvars)

% F = daeFunction(Fdae_red, q_red, Qnc);
% F15i = @(t, q_red, dq_red)  F(t, q_red, dq_red, inputs_k(q_red));

% [y0, dy0] = decic(F15i, 0, q0, [], 0*q0, []);

% opt = odeset('MassSingular','yes', 'OutputFcn','odeplot', ...
%              'RelTol', 1e-3, 'AbsTol' ,1e-3,'Stats','on');
% [tsol, ysol] = ode15i(F15i, tspan, y0, dy0, opt);
