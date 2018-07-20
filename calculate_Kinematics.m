% sus = Front or Rear suspension parameters

function [ sus ] = calculate_Kinematics(prefix_left, prefix_right, sus, tyr, WH_stroke, WH_steps, SR_stroke, SR_steps )


totaltime = tic;

fprintf('Initializing Kinematics...\n')

digits(4)

% LEGEND:
%   u = Rotation or Displacement vector
%   b = Rotation angle or Displacement magnitude
%   c = Rotation center
%   v = Loop position vectors
%   n = Normal vector
%   Components (x=_1, y=_2 and z=_3)

%   UA = Upper Arm
%   UR = Upright
%   LA = Lower Arm
%   CH = Chassis
%   SR = Steering Rack
%   TR = Tie Rod
%   WH = Wheel
%   TY = Tyre
%   DS = Damper/Spring
%   RK = Rocker
%   PR = Push Rod
%   AB = ARB Blade
%   AL = ARB Link
%   AR = ARB Bar


%% Define Main Coordinates

syms bSR
LAUR = sym('LAUR_%d', [3 1]);

% Main Coordinates
q = [LAUR(3); 
     bSR];


%% Front left suspension kinematics

uUA = sus.L.init.UACH1_0 - sus.L.init.UACH2_0;         
uLA = sus.L.init.LACH1_0 - sus.L.init.LACH2_0;
cFUA = (sus.L.init.UACH1_0 + sus.L.init.UACH2_0)/2;    
cFLA = (sus.L.init.LACH1_0 + sus.L.init.LACH2_0)/2;

vUA = cFUA - sus.L.init.UAUR_0;          
vUR = sus.L.init.UAUR_0 - sus.L.init.LAUR_0;       
vLA = sus.L.init.LAUR_0 - cFLA;         
vCH = cFLA - cFUA;

% New coordinate variables to be solved (8 Variables)
UAUR = sym('UAUR_%d', [3 1]);
syms bUA bLA  

% Loop vector functions - (8 Variables, 7 Equations = 1 main CO, 7 secondary COs)
f_S = vpa([ rotateVector(vLA,uLA,bLA) + (UAUR-LAUR) + rotateVector(vUA,uUA,bUA) + vCH;
             rotatePoint(sus.L.init.UAUR_0,uUA,bUA,cFUA) - UAUR ;  %vpa(rotatePoint(sus.L.init.LAUR_0,uLA,bLA,cFLA) - LAUR) ;
             norm(UAUR-LAUR) - norm(sus.L.init.UAUR_0-sus.L.init.LAUR_0)    ]);

s_S   = [LAUR(1:2);   UAUR;   bUA; bLA];
s_S_0 = [sus.L.init.LAUR_0(1:2); sus.L.init.UAUR_0; 0;    0];
 
q_S = [LAUR(3)];
          
J_S = jacobian(f_S,s_S);
B_S = -jacobian(f_S,q_S);
K_S = (J_S\B_S);
K_S = K_S*( jacobian(q_S,q) );      % Convert to main coordinates
         

%% Steering system kinematics

uSR = [0 1 0]';

% New coordinate variables to be solved (7 Variables)
SRTR = sym('SRTR_%d', [3 1]);
TRUR = sym('TRUR_%d', [3 1]);

% Loop vector functions - (13 Variables, 6 Equations = 7 main CO, 6 secondary COs)
f_SS = vpa([ sus.L.init.SRTR_0 + bSR*uSR - SRTR;
             norm(TRUR-SRTR) - norm(sus.L.init.TRUR_0-sus.L.init.SRTR_0);
             norm(LAUR-TRUR) - norm(sus.L.init.LAUR_0-sus.L.init.TRUR_0);
             norm(UAUR-TRUR) - norm(sus.L.init.UAUR_0-sus.L.init.TRUR_0)   ]);

s_SS   = [SRTR; TRUR];
s_SS_0 = [sus.L.init.SRTR_0; sus.L.init.TRUR_0];
 
q_SS = [LAUR; UAUR; bSR];
          
J_SS = jacobian(f_SS,s_SS);
B_SS = -jacobian(f_SS,q_SS);
K_SS = (J_SS\B_SS);          % K = simplify(inv(J)*B);         % ds/dt = K * dq/dt   
K_SS = K_SS*( jacobian(q_SS,s_S)*K_S + jacobian(q_SS,q) );



%% Front Wheel kinematics


% New coordinate variables to be solved (6 Variables)
uWH = sym('uWH_%d', [3 1]);     % Wheel rotation vector
cWH = sym('cWH_%d', [3 1]);


% c0 = [DNA.DIN_CO(1) DNA.track/2+sus.L.setup.camber*tyr.unloadedRadius DNA.DIN_CO(3)]';
c0 = sus.L.init.cWH_0;
v0 = rotateVector(rotateVector([0 1 0]',[1 0 0]',sus.L.setup.camber),[0 0 1]',-sus.L.setup.toe);                       

g_WH = matlabFunction(vpa([...
    createTmatrixFunction_point(LAUR,TRUR,UAUR,sus.L.init.LAUR_0,sus.L.init.TRUR_0,sus.L.init.UAUR_0,c0);
    createTmatrixFunction_vector(LAUR,TRUR,UAUR,sus.L.init.LAUR_0,sus.L.init.TRUR_0,sus.L.init.UAUR_0,v0) ]));

s_WH   = [cWH; uWH];
q_WH = [LAUR; UAUR; TRUR];

K_WH = vpa(jacobian(g_WH,q_WH));
K_WH = K_WH*( jacobian(q_WH,s_S)*K_S + jacobian(q_WH,s_SS)*K_SS + jacobian(q_WH,q) );


%% Front tyre contact point

% New coordinate variables to be solved (7 Variables)
TYRO = sym('TYRO_%d', [3 1]);

z = [0 0 1]';
d_vector = cross(uWH,cross(-z,uWH));

g_TY = matlabFunction(vpa( d_vector/norm(d_vector)*tyr.unloadedRadius + cWH ));

s_TY   = [TYRO];
q_TY = [uWH; cWH];

K_TY = vpa(jacobian(g_TY,q_TY));
K_TY = K_TY*( jacobian(q_TY,s_WH)*K_WH );
         

%% Push Rod kinematics

PRUA = sym('PRUA_%d', [3 1]);

% g_FPR = PRUA
g_PR = matlabFunction(vpa([...
    createTmatrixFunction_point(sus.L.init.UACH1_0,sus.L.init.UACH2_0,UAUR,sus.L.init.UACH1_0,sus.L.init.UACH2_0,sus.L.init.UAUR_0,sus.L.init.PRUA_0) ]));

s_PR = PRUA;
q_PR = UAUR;

K_PR = vpa(jacobian(g_PR,q_PR));
K_PR = K_PR*( jacobian(q_PR,s_S)*K_S );

%% Rocker kinematics

RKPR = sym('RKPR_%d', [3 1]);
syms bRK

uRK = cross((sus.L.init.RKPR_0-sus.L.init.RKCH_0),(sus.L.init.RKDS_0-sus.L.init.RKCH_0));

f_RK = vpa([ rotatePoint(sus.L.init.RKPR_0,uRK,bRK,sus.L.init.RKCH_0) - RKPR;
              norm(RKPR-PRUA)-norm(sus.L.init.RKPR_0-sus.L.init.PRUA_0) ]);

s_RK   = [RKPR; bRK];
s_RK_0 = [sus.L.init.RKPR_0; 0];
q_RK   = [PRUA];

J_RK = jacobian(f_RK,s_RK);
B_RK = -jacobian(f_RK,q_RK);
K_RK = (J_RK\B_RK);
K_RK = K_RK * ( jacobian(q_RK,s_PR)*K_PR );

%% Damper/Spring Kinematics

RKDS = sym('RKDS_%d', [3 1]);
syms bDS                            % Damper/Spring length

% g_FDS = [RKDS ; bDS]
g_DS = matlabFunction(vpa([ rotatePoint(sus.L.init.RKDS_0,uRK,bRK,sus.L.init.RKCH_0);
                            norm(rotatePoint(sus.L.init.RKDS_0,uRK,bRK,sus.L.init.RKCH_0) - sus.L.init.DSCH_0) ]));

s_DS   = [RKDS; bDS];
q_DS   = bRK;

K_DS = vpa(jacobian(g_DS,q_DS));
K_DS = K_DS * ( jacobian(q_DS,s_RK)*K_RK );


%% ARB Kinematics

uAR = [0 1 0]';

RKAL = sym('RKAL_%d', [3 1]);
ALAB = sym('ALAB_%d', [3 1]);
syms bAR

f_AR = vpa([ rotatePoint(sus.L.init.RKAL_0,uRK,bRK,sus.L.init.RKCH_0) - RKAL;
             norm(ALAB-RKAL) - norm(sus.L.init.ALAB_0-sus.L.init.RKAL_0);
             rotatePoint(sus.L.init.ALAB_0,uAR,bAR,sus.L.init.ALAR_0) - ALAB ]);

s_AR   = [RKAL; ALAB; bAR];
s_AR_0 = [sus.L.init.RKAL_0; sus.L.init.ALAB_0; 0];
 
q_AR = [bRK];
          
J_AR = jacobian(f_AR,s_AR);
B_AR = -jacobian(f_AR,q_AR);
K_AR = (J_AR\B_AR);          % K = simplify(inv(J)*B);         % ds/dt = K * dq/dt         
K_AR = K_AR * ( jacobian(q_AR,s_RK)*K_RK );




%% Jacobian Matrices

s_Kmatrix = [s_S ; s_SS ; s_WH ; s_TY ; s_PR ; s_RK ; s_DS ; s_AR];
K_all     = [K_S ; K_SS ; K_WH ; K_TY ; K_PR ; K_RK ; K_DS ; K_AR];

dispstat('','init')
tic
for i=1:numel(s_Kmatrix)
    dispstat(sprintf('    Simplifying K and L Matrices %d/%d: %.1fs',i,numel(s_Kmatrix),toc))
    sus.L.Kmatrix.(char(s_Kmatrix(i))) = K_all(i,:);
    sus.L.Lmatrix.(char(s_Kmatrix(i))) = jacobian(sus.L.Kmatrix.(char(s_Kmatrix(i))),q);
end


fprintf('    Total calculation time of [K] and [L] matrices: %.1f s\n',toc)


%% Calculate secundary coordinates

dispstat('','init')
dispstat('Mapping Suspension Kinematics ','keepthis')
 

[LAUR_3_val,bSR_val] = meshgrid(linspace(sus.L.init.LAUR_0(3)+WH_stroke(1),sus.L.init.LAUR_0(3)+WH_stroke(2),WH_steps), linspace(SR_stroke(1),SR_stroke(2),SR_steps));
s_all = cat(3,LAUR_3_val,bSR_val);
s_all = syms2structArray(q,s_all);


% FS
tic
for i = 1:size(s_all,1)
    for j = 1:size(s_all,2)
        toSolve = subs(f_S, fieldnames(s_all), struct2array(s_all(i,j))' );      %     toSolve = subs(f_S, q_S, s_all(i).(char(q_S)) );
        s_S_val(i,j) = vpasolve( toSolve, s_S, s_S_0);         
        s_S_0 = struct2array(s_S_val(i,j))';
        dispstat(sprintf('    Left suspension %d/%d: %.1fs',(i-1)*size(s_all,2)+j,numel(s_all),toc))
    end
end
s_all = catstruct(s_all,s_S_val);
dispstat(' ','keepprev');

% FSS
tic
for i = 1:size(s_all,1)
    for j = 1:size(s_all,2)
        toSolve = subs(f_SS, fieldnames(s_all), struct2array(s_all(i,j))'); %     toSolve = subs(f_SS, [s_FS;q], [struct2array(s_S_val(i))';struct2array(s_all(:,i))']);
        s_SS_val(i,j) = vpasolve( toSolve, s_SS, s_SS_0);
        s_SS_0 = struct2array(s_SS_val(i,j))';
        dispstat(sprintf('    Steering System %d/%d: %.1fs',(i-1)*size(s_all,2)+j,numel(s_all),toc))
    end
end
s_all = catstruct(s_all,s_SS_val);
dispstat(' ','keepprev');

% FWH
tic
for i = 1:size(s_all,1)
    for j = 1:size(s_all,2)
        s_WH_0(i,j,:) = g_WH(s_all(i,j).LAUR_1, s_all(i,j).LAUR_2, s_all(i,j).LAUR_3, s_all(i,j).TRUR_1, s_all(i,j).TRUR_2, s_all(i,j).TRUR_3, s_all(i,j).UAUR_1, s_all(i,j).UAUR_2, s_all(i,j).UAUR_3 );   
        dispstat(sprintf('    Left Wheel %d/%d: %.1fs',(i-1)*size(s_all,2)+j,numel(s_all),toc))
    end
end
s_FWH_val = syms2structArray(s_WH,s_WH_0);
s_all = catstruct(s_all,s_FWH_val);
dispstat(' ','keepprev');


% FTY
tic
for i = 1:size(s_all,1)
    for j = 1:size(s_all,2)
        s_TY_0(i,j,:) = g_TY(s_all(i,j).cWH_1,s_all(i,j).cWH_2,s_all(i,j).cWH_3,s_all(i,j).uWH_1,s_all(i,j).uWH_2,s_all(i,j).uWH_3);
        dispstat(sprintf('    Left Tyre %d/%d: %.1fs',(i-1)*size(s_all,2)+j,numel(s_all),toc))
    end
end
s_FTY_val = syms2structArray(s_TY,s_TY_0);
s_all = catstruct(s_all,s_FTY_val);
dispstat(' ','keepprev');

% FPR
tic
for i = 1:size(s_all,1)
    for j = 1:size(s_all,2)
        s_PR_0(i,j,:) = g_PR(s_all(i,j).UAUR_1,s_all(i,j).UAUR_2,s_all(i,j).UAUR_3);
        dispstat(sprintf('    Push Rod %d/%d: %.1fs',(i-1)*size(s_all,2)+j,numel(s_all),toc))
    end
end
s_FPR_val = syms2structArray(s_PR,s_PR_0);
s_all = catstruct(s_all,s_FPR_val);
dispstat(' ','keepprev');

% FRK
tic
for i = 1:size(s_all,1)
    for j = 1:size(s_all,2)
        toSolve = subs(f_RK, fieldnames(s_all), struct2array(s_all(i,j))');
        s_RK_val(i,j) = vpasolve( toSolve, s_RK, s_RK_0);
        s_RK_0 = struct2array(s_RK_val(i,j))';
        dispstat(sprintf('    Rocker %d/%d: %.1fs',(i-1)*size(s_all,2)+j,numel(s_all),toc))
    end
end
s_all = catstruct(s_all,s_RK_val);
dispstat(' ','keepprev');

% FDS
tic
for i = 1:size(s_all,1)
    for j = 1:size(s_all,2)
        s_DS_0(i,j,:) = g_DS(s_all(i,j).bRK);
        dispstat(sprintf('    Damper/Spring %d/%d: %.1fs',(i-1)*size(s_all,2)+j,numel(s_all),toc))
    end
end
s_FDS_val = syms2structArray(s_DS,s_DS_0);
s_all = catstruct(s_all,s_FDS_val);
dispstat(' ','keepprev');


% FAR
tic
for i = 1:size(s_all,1)
    for j = 1:size(s_all,2)
        toSolve = subs(f_AR, fieldnames(s_all), struct2array(s_all(i,j))');
        s_AR_val(i,j) = vpasolve( toSolve, s_AR, s_AR_0);
        s_AR_0 = struct2array(s_AR_val(i,j))';
        dispstat(sprintf('    Anti-Roll-Bar %d/%d: %.1fs',(i-1)*size(s_all,2)+j,numel(s_all),toc))
    end
end
s_all = catstruct(s_all,s_AR_val);
dispstat(' ','keepprev');



%% FRONT SUSPENSION SUMMARY

fprintf('Converting structures and making simetry...\n'); tic;

sus.L.kin = structarray2struct(s_all);

sus.R = sus.L;

variables2inverse = {'bSR','LAUR_2','UAUR_2','cWH_2','PRUA_2','RKDS_2','RKPR_2','SRTR_2','TRUR_2','TYRO_2','UAUR_2','uWH_2','ALAB_2','RKAL_2'};
for i=1:numel(variables2inverse)
    sus.R.kin.(variables2inverse{i})     = -sus.L.kin.(variables2inverse{i});
    if isfield(sus.R.Kmatrix,variables2inverse{i})
        sus.R.Kmatrix.(variables2inverse{i}) = -sus.R.Kmatrix.(variables2inverse{i});
        sus.R.Lmatrix.(variables2inverse{i}) = -sus.R.Lmatrix.(variables2inverse{i});
    end
end
init2inverse = {'UACH1_0','UACH2_0','LACH1_0','LACH2_0','DSCH_0','RKCH_0','ALAR_0'};
for i=1:numel(init2inverse)
    sus.R.init.(init2inverse{i}) = sus.L.init.(init2inverse{i}).*[1,-1,1]';
end

fprintf('    Finished %.1fs\n',toc(totaltime));


end

