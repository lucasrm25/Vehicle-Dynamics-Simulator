% Constrained Lagrange Approach

clear all;
close all;
clc;

%% Build Lagrange dynamics

dim = 2;
e = eye(dim);

M = 100*e;
UM = 10*e;
I = 10;
g = [0 -9.8]';
k = 100;
s0 = 0.5;

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

T = 0.5 * q(:,2).' * W * q(:,2);
V = - G.' * W * q(:,1) + 0.5 * dx.' * K * dx;

% acosd( dot(u,v)/(normv(u)*normv(v)) )
C = [ normv(LA(:,1)-CG(:,1))^2 - 0.5^2;
      normv(CG(:,1)-RA(:,1))^2 - 0.5^2;
      normv(RA(:,1)-LA(:,1))^2 - 0.8^2;
      normv(CG(:,1)-LW(:,1))^2 - 0.5^2;
      normv(CG(:,1)-RW(:,1))^2 - 0.5^2];
z = sym('z',size(C));

Qnc = sym('Qnc',[size(q,1),1]);

Lag = T - V - z.'*C;
Lag_dq = gradient(Lag,q(:,2));
LAGEQ = jacobian(Lag_dq,q(:,1))*q(:,2) + jacobian(Lag_dq,q(:,2))*q(:,3) - gradient(Lag,q(:,1)) - Qnc;

dC = jacobian(C,q(:,1))*q(:,2); % == 0
ddC = jacobian(dC,q(:,1))*q(:,2) + jacobian(dC,q(:,2))*q(:,3); % == 0


Fdae = [LAGEQ; ddC];



%% Prepare simulation


varsqt='';
for i=1:size(q,1)
    eval(sprintf('syms %s(t)',char(q(i,1))));
    varsqt = [varsqt,' ',char(q(i,1))];             %#ok<*AGROW>
end
eval(['qt=[',varsqt,'];']);
qt = qt.';

varszt='';
for i=1:size(z,1)
    eval(sprintf('syms %s(t)',char(z(i,1))));
    varszt = [varszt,' ',char(z(i,1))];
end
eval(['zt=[',varszt,'];']);
zt = zt.';


Fdae_t = subs(Fdae, [z(:,1);q(:,1);q(:,2);q(:,3)], [zt;qt;diff(qt);diff(qt,2)]);
C_t    = subs(C,    [z(:,1);q(:,1);q(:,2);q(:,3)], [zt;qt;diff(qt);diff(qt,2)]);
dC_t   = subs(dC,   [z(:,1);q(:,1);q(:,2);q(:,3)], [zt;qt;diff(qt);diff(qt,2)]);

[Fdae_red, q_red, Rel] = reduceDifferentialOrder(Fdae_t, [qt;zt]);

if ~isLowIndexDAE(Fdae_red, q_red)
    error('DAE Index > 1 !!!')
end

F = daeFunction(Fdae_red, q_red, Qnc);
f = @(t, q_red, dq_red)  F(t, q_red, dq_red, inputs_k(q_red));


tspan = [0,2];
cgy = 0;
cgz = 2;
q0 = [cgy cgz cgy-0.5 cgz cgy+0.5 cgz cgy-0.4 cgz+0.3 cgy+0.4 cgz+0.3 ...
      0 0 0 0 0 ...
      0 0 0 0 0 0 0 0 0 0]'; 
[y0, yp0] = decic(f, 0, q0, [], 0*q0, [])

opt = odeset('RelTol', 1e-3, 'AbsTol' ,1e-3, 'OutputFcn','odeplot');
[tsol, ysol] = ode15i(f, tspan, y0, yp0, opt)


figure; hold on; grid on;
axis manual
xlim([-5 5])
ylim([-5 5])
%     scatter(ysol(i,[1,3,5,7,9]), ysol(i,[2,4,6,8,10]), 50, 'filled');
pw = plot( ysol(1,[3,1,5]), ysol(1,[4,2,6]),'LineWidth',3);
pb = plot( ysol(1,[1,7,9,1]), ysol(1,[2,8,10,2]),'LineWidth',3);
for i=1:numel(tsol)
    tic
    pw.XData = ysol(i,[3,1,5]);
    pw.YData = ysol(i,[4,2,6]);
    pb.XData = ysol(i,[1,7,9,1]);
    pb.YData = ysol(i,[2,8,10,2]);
    drawnow();
    pause( 1/20 - toc )
end

% f(0,q0,q0)
% inputs_k(q0)
% 
% 
% [M_ode,b_ode] = massMatrixForm(Fdae_red, q_red);
% M_ode_fun = odeFunction(M_ode, q_red, Qnc);
% b_ode_fun = odeFunction(b_ode, q_red, Qnc);
% M_ode_fun_k = @(t,q_red) M_ode_fun(t, q_red, inputs_k(q_red));
% b_ode_fun_k = @(t,q_red) b_ode_fun(t, q_red, inputs_k(q_red));          
%                 
% opt = odeset('Mass', M_ode_fun_k, 'RelTol', 1e-3, 'AbsTol' ,1e-4, 'OutputFcn','odeplot','MassSingular','yes');
% [tSol,ySol] = ode23t(b_ode_fun_k, tspan, q0, opt);


function Qnc = inputs_k(q_red)
    S0 = 1;
    if q_red(4)<=S0
        FL = 100*(S0 - q_red(4));
    else
        FL = 0;
    end
    if q_red(6)<=S0
        FR = 100*(S0 - q_red(6));
    else
        FR = 0;
    end
    Qnc = q_red(1:10).*[0 0 0 FL 0 FR 0 0 0 0].';
end

function n = normv(v)
    n = (v.'*v)^.5;
end



% [A,b] = equationsToMatrix( [LAGEQ; ddC] , [q(:,3); z] );
% 
% q_red = [q(:,1);q(:,2);z];
% A_red = blkdiag(eye(size(q,1)),A);
% b_red = [q(:,2);b];
% 
% A_fun = matlabFunction(A_red, q_red);
% b_fun = matlabFunction(b_red, q_red);
% 
% M_ode_fun_k = @(t,q_red) A_fun(q_red, inputs_k(q_red));
% b_ode_fun_k = @(t,q_red) b_fun(q_red, inputs_k(q_red));          
%                 
% 
% opt = odeset('Mass', M_ode_fun_k, 'RelTol', 1e-3, 'AbsTol' ,1e-4, 'OutputFcn','odeplot');
% [tSol,ySol] = ode23t(b_ode_fun_k, tspan, [q_k; dq_k], opt);





