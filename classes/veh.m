classdef veh < matlab.mixin.Copyable  

% Considerations:
%
% S are considered as LPV parameters, i.e. constants in relation to derivatives
% => dK/ds = 0
%
%
% To implement:
%
% dK/dt = 0   (more radical)
%
%   that is: all S_dq_XXX = 0
%
%   POSSIVEL SOLUCAO... DEIXAR MATRIZES K e L COMO FUNCOES E NAO
%   CALCULAR AGORA POIS DEMORA MUITO E DEPOIS O PROGRAMA TENTA
%   DERIVAR ESSAS MATRIZES ENORMES...
%
%   OUTRA POSSIVEL SOLUCAO... CONSIDERAR dK/dt = 0.. ou seja
%   derivadas temporais das matrizes K como nulas... na verdade
%   eh a mesma solucao a cima.. ja que matrizes K e L se
%   tornariam agora "constantes"
%
%   PODEMOS CONFERIR qual a taxa de variacao dK/dt plotando na tela
%   Ou seja... comparar o efeito de dK/dt * qp e K * qpp

    
    
%% Properties
    properties
        DNA
        sus_fl
        sus_fr
        sus_rl
        sus_rr
        
        q
        q_red
        
        dynamicFunction
    end
    properties (GetAccess=private)
    end
    properties (Constant)
        grav = [0 0 9.81]';
    end
    properties (Dependent)
    end

%% Methods    
    methods(Static)
    end
    methods
        function obj = veh(DNA, veh_sus_fl, veh_sus_fr, veh_sus_rl, veh_sus_rr)
            if isa(veh_sus_fl,'veh_sus') && isa(veh_sus_fl,'veh_sus') && isa(veh_sus_fl,'veh_sus') && isa(veh_sus_fl,'veh_sus')
                obj.sus_fl = veh_sus_fl;
                obj.sus_fr = veh_sus_fr;
                obj.sus_rl = veh_sus_rl;
                obj.sus_rr = veh_sus_rr;
            else
                error('Wrong class type for the suspensions')
            end
            obj.DNA = DNA;
        end

        function sim(obj, q0, qp0 ,e_fl_Ftyre, e_fr_Ftyre, e_rl_Ftyre, e_rr_Ftyre)
%             https://www.mathworks.com/help/symbolic/equation-solving.html
%             https://www.mathworks.com/help/symbolic/equation-solving.html
%             https://www.mathworks.com/help/symbolic/solve-differential-algebraic-equations.html

%                     q_tm1  = [0 0 pi/4   10 10 10  obj.sus_fl.init.LAUR_0(3) obj.sus_fr.init.LAUR_0(3) obj.sus_rl.init.LAUR_0(3) obj.sus_rr.init.LAUR_0(3)  0 0];
%                     qp_tm1 = [0 0 0 0 0 0 0 0 0 0 0 0];

%             e_fl_Ftyre = [0 0 tyr.stiffness*(tyr.compressionLength-v_fl_TYRO(3))*(v_fl_TYRO(3)<tyr.compressionLength)];
%             e_fr_Ftyre = [0 0 tyr.stiffness*(tyr.compressionLength-v_fr_TYRO(3))*(v_fr_TYRO(3)<tyr.compressionLength)];
%             e_rl_Ftyre = [0 0 tyr.stiffness*(tyr.compressionLength-v_rl_TYRO(3))*(v_rl_TYRO(3)<tyr.compressionLength)];
%             e_rr_Ftyre = [0 0 tyr.stiffness*(tyr.compressionLength-v_rr_TYRO(3))*(v_rr_TYRO(3)<tyr.compressionLength)];

        end
        
        function calculateDynamics(obj)
            
            % Reshape functions
            obj.sus_fl.reshape(12,[7  11]);
            obj.sus_fr.reshape(12,[8  11]);
            obj.sus_rl.reshape(12,[9  12]);
            obj.sus_rr.reshape(12,[10 12]);

            obj.sus_fl.set_prefix('fl_');
            obj.sus_fr.set_prefix('fr_');
            obj.sus_rl.set_prefix('rl_');
            obj.sus_rr.set_prefix('rr_');
            
            fprintf('Initializing Dynamics...\n')
            
            % Tyre reaction forces at the tyre contact point - Inertial CO
            e_fl_Ftyre = sym('e_fl_Ftyre_%d', [3 1]);
            e_fr_Ftyre = sym('e_fr_Ftyre_%d', [3 1]);
            e_rl_Ftyre = sym('e_rl_Ftyre_%d', [3 1]);
            e_rr_Ftyre = sym('e_rr_Ftyre_%d', [3 1]);
            
            syms e_ang_phi(t) e_ang_teta(t) e_ang_psi(t)      e_pos_x(t) e_pos_y(t) e_pos_z(t)
            syms v_fl_LAUR_3(t) v_fr_LAUR_3(t) v_rl_LAUR_3(t) v_rr_LAUR_3(t)
            syms v_fl_bSR(t) v_fr_bSR(t)
            q     = [e_ang_phi e_ang_teta e_ang_psi e_pos_x e_pos_y e_pos_z v_fl_LAUR_3 v_fr_LAUR_3 v_rl_LAUR_3 v_rr_LAUR_3 v_fl_bSR v_fr_bSR].';
            q     = q(t);       % Convert symfun to syms array
            qp    = diff(q,t);
            qpp   = diff(qp,t);


            % Jacobian Matrix - vehicle CO to inertial CO
            e_T_v =  [cos(e_ang_psi)*cos(e_ang_teta) -sin(e_ang_psi)*cos(e_ang_phi)+cos(e_ang_psi)*sin(e_ang_teta)*sin(e_ang_phi)  sin(e_ang_psi)*sin(e_ang_phi)+cos(e_ang_psi)*sin(e_ang_teta)*cos(e_ang_phi);
                      sin(e_ang_psi)*cos(e_ang_teta)  cos(e_ang_psi)*cos(e_ang_phi)+sin(e_ang_psi)*sin(e_ang_teta)*sin(e_ang_phi) -cos(e_ang_psi)*sin(e_ang_phi)+sin(e_ang_psi)*sin(e_ang_teta)*cos(e_ang_phi);
                     -sin(e_ang_teta)                 cos(e_ang_teta)*sin(e_ang_phi)                                               cos(e_ang_teta)*cos(e_ang_phi)];

            e_T_v_dq = diff2(e_T_v,q);

            % Vehicle angular velocity - vehicle CO
            v_angp   = [-sin(e_ang_teta)                0               1;
                        cos(e_ang_teta)*sin(e_ang_phi)  cos(e_ang_phi)  0;
                        cos(e_ang_teta)*cos(e_ang_phi) -sin(e_ang_phi)  0] ...
                        * [diff(e_ang_psi) diff(e_ang_teta) diff(e_ang_phi)].';
            v_angpp    = diff(v_angp);
            v_angp_dq  = diff2(v_angp,q);
            v_angp_dqp = diff2(v_angp,qp);

            %Vehicle position - Inertial CO
            e_pos      = [e_pos_x e_pos_y e_pos_z].';
            e_pos_dq   = diff2(e_pos,q);
            e_posp     = diff(e_pos);
            e_posp_dq  = diff2(e_posp,q);
            e_posp_dqp = diff2(e_posp,qp);
            e_pospp    = diff(e_posp);


            % Set Inputs q
%             obj.sus_fl.q_val = q_tm1([7 11])';
%             obj.sus_fr.q_val = q_tm1([8 11])';
%             obj.sus_rl.q_val = q_tm1([9 12])';
%             obj.sus_rr.q_val = q_tm1([10 12])';

            % Get Outputs s
            v_fl_TYRO = [obj.sus_fl.s_prefix.TYRO_1 obj.sus_fl.s_prefix.TYRO_2 obj.sus_fl.s_prefix.TYRO_3].';
            v_fr_TYRO = [obj.sus_fr.s_prefix.TYRO_1 obj.sus_fr.s_prefix.TYRO_2 obj.sus_fr.s_prefix.TYRO_3].';
            v_rl_TYRO = [obj.sus_rl.s_prefix.TYRO_1 obj.sus_rl.s_prefix.TYRO_2 obj.sus_rl.s_prefix.TYRO_3].';
            v_rr_TYRO = [obj.sus_rr.s_prefix.TYRO_1 obj.sus_rr.s_prefix.TYRO_2 obj.sus_rr.s_prefix.TYRO_3].';

            v_fl_cWH  = [obj.sus_fl.s_prefix.cWH_1 obj.sus_fl.s_prefix.cWH_2 obj.sus_fl.s_prefix.cWH_3].';
            v_fr_cWH  = [obj.sus_fr.s_prefix.cWH_1 obj.sus_fr.s_prefix.cWH_2 obj.sus_fr.s_prefix.cWH_3].';
            v_rl_cWH  = [obj.sus_rl.s_prefix.cWH_1 obj.sus_rl.s_prefix.cWH_2 obj.sus_rl.s_prefix.cWH_3].';
            v_rr_cWH  = [obj.sus_rr.s_prefix.cWH_1 obj.sus_rr.s_prefix.cWH_2 obj.sus_rr.s_prefix.cWH_3].';

            v_fl_bDS    = obj.sus_fl.s_prefix.bDS;
            v_fr_bDS    = obj.sus_fr.s_prefix.bDS;
            v_rl_bDS    = obj.sus_rl.s_prefix.bDS;
            v_rr_bDS    = obj.sus_rr.s_prefix.bDS;


            
            v_fl_cWH_dq =  [ obj.sus_fl.reshapefun(obj.sus_fl.Kmatrix.cWH_1, q, struct2array(obj.sus_fl.s_prefix).') ;
                             obj.sus_fl.reshapefun(obj.sus_fl.Kmatrix.cWH_2, q, struct2array(obj.sus_fl.s_prefix).') ;
                             obj.sus_fl.reshapefun(obj.sus_fl.Kmatrix.cWH_3, q, struct2array(obj.sus_fl.s_prefix).') ];
            v_fl_cWH_ddq = diff2(v_fl_cWH_dq,q);

            v_fr_cWH_dq =  [ obj.sus_fr.reshapefun(obj.sus_fr.Kmatrix.cWH_1, q, struct2array(obj.sus_fr.s_prefix).') ;
                             obj.sus_fr.reshapefun(obj.sus_fr.Kmatrix.cWH_2, q, struct2array(obj.sus_fr.s_prefix).') ;
                             obj.sus_fr.reshapefun(obj.sus_fr.Kmatrix.cWH_3, q, struct2array(obj.sus_fr.s_prefix).') ];
            v_fr_cWH_ddq = diff2(v_fr_cWH_dq,q);

            v_rl_cWH_dq =  [ obj.sus_rl.reshapefun(obj.sus_rl.Kmatrix.cWH_1, q, struct2array(obj.sus_rl.s_prefix).') ;
                             obj.sus_rl.reshapefun(obj.sus_rl.Kmatrix.cWH_2, q, struct2array(obj.sus_rl.s_prefix).') ;
                             obj.sus_rl.reshapefun(obj.sus_rl.Kmatrix.cWH_3, q, struct2array(obj.sus_rl.s_prefix).') ];
            v_rl_cWH_ddq = diff2(v_rl_cWH_dq,q);

            v_rr_cWH_dq =  [ obj.sus_rr.reshapefun(obj.sus_rr.Kmatrix.cWH_1, q, struct2array(obj.sus_rr.s_prefix).') ;
                             obj.sus_rr.reshapefun(obj.sus_rr.Kmatrix.cWH_2, q, struct2array(obj.sus_rr.s_prefix).') ;
                             obj.sus_rr.reshapefun(obj.sus_rr.Kmatrix.cWH_3, q, struct2array(obj.sus_rr.s_prefix).') ];
            v_rr_cWH_ddq = diff2(v_rr_cWH_dq,q);


            v_fl_bDS_dq =  obj.sus_fl.reshapefun(obj.sus_fl.Kmatrix.bDS, q, struct2array(obj.sus_fl.s_prefix).') ;
            v_fr_bDS_dq =  obj.sus_fr.reshapefun(obj.sus_fr.Kmatrix.bDS, q, struct2array(obj.sus_fr.s_prefix).') ;
            v_rl_bDS_dq =  obj.sus_rl.reshapefun(obj.sus_rl.Kmatrix.bDS, q, struct2array(obj.sus_rl.s_prefix).') ;
            v_rr_bDS_dq =  obj.sus_rr.reshapefun(obj.sus_rr.Kmatrix.bDS, q, struct2array(obj.sus_rr.s_prefix).') ;

            v_fl_TYRO_dq = [ obj.sus_fl.reshapefun(obj.sus_fl.Kmatrix.TYRO_1,q,struct2array(obj.sus_fl.s_prefix).') ;
                             obj.sus_fl.reshapefun(obj.sus_fl.Kmatrix.TYRO_2,q,struct2array(obj.sus_fl.s_prefix).') ;
                             obj.sus_fl.reshapefun(obj.sus_fl.Kmatrix.TYRO_3,q,struct2array(obj.sus_fl.s_prefix).') ];

            v_fr_TYRO_dq = [ obj.sus_fr.reshapefun(obj.sus_fr.Kmatrix.TYRO_1,q,struct2array(obj.sus_fr.s_prefix).') ;
                             obj.sus_fr.reshapefun(obj.sus_fr.Kmatrix.TYRO_2,q,struct2array(obj.sus_fr.s_prefix).') ;
                             obj.sus_fr.reshapefun(obj.sus_fr.Kmatrix.TYRO_3,q,struct2array(obj.sus_fr.s_prefix).') ];

            v_rl_TYRO_dq = [ obj.sus_rl.reshapefun(obj.sus_rl.Kmatrix.TYRO_1,q,struct2array(obj.sus_rl.s_prefix).') ;
                             obj.sus_rl.reshapefun(obj.sus_rl.Kmatrix.TYRO_2,q,struct2array(obj.sus_rl.s_prefix).') ;
                             obj.sus_rl.reshapefun(obj.sus_rl.Kmatrix.TYRO_3,q,struct2array(obj.sus_rl.s_prefix).') ];

            v_rr_TYRO_dq = [ obj.sus_rr.reshapefun(obj.sus_rr.Kmatrix.TYRO_1,q,struct2array(obj.sus_rr.s_prefix).') ;
                             obj.sus_rr.reshapefun(obj.sus_rr.Kmatrix.TYRO_2,q,struct2array(obj.sus_rr.s_prefix).') ;
                             obj.sus_rr.reshapefun(obj.sus_rr.Kmatrix.TYRO_3,q,struct2array(obj.sus_rr.s_prefix).') ];


            e_fl_cWHp = e_posp   + e_T_v * ( v_fl_cWH_dq  * qp + cross(v_angp,v_fl_cWH) );
            e_fr_cWHp = e_posp   + e_T_v * ( v_fr_cWH_dq  * qp + cross(v_angp,v_fr_cWH) );
            e_rl_cWHp = e_posp   + e_T_v * ( v_rl_cWH_dq  * qp + cross(v_angp,v_rl_cWH) );
            e_rr_cWHp = e_posp   + e_T_v * ( v_rr_cWH_dq  * qp + cross(v_angp,v_rr_cWH) );
            e_CGp     = e_posp   + e_T_v * (                     cross(v_angp,obj.DNA.CG_s) ); 

            e_angp     = e_T_v * v_angp;
            e_angp_dq  = diff2(e_angp,q);
            e_angp_dqp = diff2(e_angp,qp);            
            
            
            % Calculation o the dd/dt^2 of each S (secondary coordinates)
            %   obs: we are ignoring de partial derivative in relation to S
            %        because or we disconsider its partial derivative or
            %        the diff only sees q as (t) dependent
            v_fl_cWHpp = cellsum(cellfun(@(a,b) a*b, v_fl_cWH_ddq, sym2cell(qp), 'UniformOutput', false)) * qp  + v_fl_cWH_dq * qpp;   % or diff(K*qp, t)
            v_fr_cWHpp = cellsum(cellfun(@(a,b) a*b, v_fr_cWH_ddq, sym2cell(qp), 'UniformOutput', false)) * qp  + v_fr_cWH_dq * qpp;
            v_rl_cWHpp = cellsum(cellfun(@(a,b) a*b, v_rl_cWH_ddq, sym2cell(qp), 'UniformOutput', false)) * qp  + v_rl_cWH_dq * qpp;
            v_rr_cWHpp = cellsum(cellfun(@(a,b) a*b, v_rr_cWH_ddq, sym2cell(qp), 'UniformOutput', false)) * qp  + v_rr_cWH_dq * qpp;           
            e_T_v_p    = simplify(diff(e_T_v));      
            v_fl_cWH_dq_p = diff(v_fl_cWH_dq);
            v_fr_cWH_dq_p = diff(v_fr_cWH_dq);
            v_rl_cWH_dq_p = diff(v_rl_cWH_dq);
            v_rr_cWH_dq_p = diff(v_rr_cWH_dq);
            e_angpp = diff(e_angp);            
            v_angp_dqp_p  = cellfun(@diff, v_angp_dqp, 'UniformOutput', false);     % function of Q only. Not S dependent
            e_angp_dqp_p  = cellfun(@diff, e_angp_dqp, 'UniformOutput', false);     % first part of v_fl_cWHpp can be also defined as this equation
            

            e_fl_cWHpp = e_pospp + e_T_v * ( cross(v_angpp,v_fl_cWH) + cross(v_angp,cross(v_angp,v_fl_cWH)) + 2*cross(v_angp,v_fl_cWH_dq*qp) + v_fl_cWHpp );
            e_fr_cWHpp = e_pospp + e_T_v * ( cross(v_angpp,v_fr_cWH) + cross(v_angp,cross(v_angp,v_fr_cWH)) + 2*cross(v_angp,v_fr_cWH_dq*qp) + v_fr_cWHpp );
            e_rl_cWHpp = e_pospp + e_T_v * ( cross(v_angpp,v_rl_cWH) + cross(v_angp,cross(v_angp,v_rl_cWH)) + 2*cross(v_angp,v_rl_cWH_dq*qp) + v_rl_cWHpp );
            e_rr_cWHpp = e_pospp + e_T_v * ( cross(v_angpp,v_rr_cWH) + cross(v_angp,cross(v_angp,v_rr_cWH)) + 2*cross(v_angp,v_rr_cWH_dq*qp) + v_rr_cWHpp );

            e_CGpp     = e_pospp + e_T_v * ( cross(v_angpp,obj.DNA.CG_s) + cross(v_angp,cross(v_angp,obj.DNA.CG_s)) );
          
            
            fprintf('    Calculating differential equations...\n'); tic;
            
            if isempty(gcp('nocreate')) parpool(4); end
            parfor_progress(pwd,numel(q));
            parfor i=1:numel(q)
                
                % We could have easily done diff(Epot, q). But here we are not considering S as q dependent
                % so we must make the calculations by hand
                Epot_dq{i} = vpa(  -obj.sus_fl.unsMass * obj.grav' * (e_T_v_dq{i} * v_fl_cWH  + e_T_v * v_fl_cWH_dq(:,i) + e_pos_dq{i}) + ...
                                   -obj.sus_fr.unsMass * obj.grav' * (e_T_v_dq{i} * v_fr_cWH  + e_T_v * v_fr_cWH_dq(:,i) + e_pos_dq{i}) + ...
                                   -obj.sus_rl.unsMass * obj.grav' * (e_T_v_dq{i} * v_rl_cWH  + e_T_v * v_rl_cWH_dq(:,i) + e_pos_dq{i}) + ...
                                   -obj.sus_rr.unsMass * obj.grav' * (e_T_v_dq{i} * v_rr_cWH  + e_T_v * v_rr_cWH_dq(:,i) + e_pos_dq{i}) + ...
                                   -obj.DNA.susMass    * obj.grav' * (e_T_v_dq{i} * obj.DNA.CG_s                         + e_pos_dq{i}) + ...
                                   obj.sus_fl.setup.kspring   * (v_fl_bDS - norm(obj.sus_fl.init.RKDS_0-obj.sus_fl.init.DSCH_0)) * v_fl_bDS_dq(:,i) + ...
                                   obj.sus_fr.setup.kspring   * (v_fr_bDS - norm(obj.sus_fr.init.RKDS_0-obj.sus_fr.init.DSCH_0)) * v_fr_bDS_dq(:,i) + ...
                                   obj.sus_rl.setup.kspring   * (v_rl_bDS - norm(obj.sus_rl.init.RKDS_0-obj.sus_rl.init.DSCH_0)) * v_rl_bDS_dq(:,i) + ...
                                   obj.sus_rr.setup.kspring   * (v_rr_bDS - norm(obj.sus_rr.init.RKDS_0-obj.sus_rr.init.DSCH_0)) * v_rr_bDS_dq(:,i) );

                e_fl_cWHp_dq{i} = e_posp_dq{i} + e_T_v_dq{i} * ( v_fl_cWH_dq     * qp + cross(v_angp,v_fl_cWH) ) ...
                                               + e_T_v       * ( v_fl_cWH_ddq{i} * qp + cross(v_angp_dq{i},v_fl_cWH) + cross(v_angp,v_fl_cWH_dq(:,i))   );
                e_fr_cWHp_dq{i} = e_posp_dq{i} + e_T_v_dq{i} * ( v_fr_cWH_dq     * qp + cross(v_angp,v_fr_cWH) ) ...
                                               + e_T_v       * ( v_fr_cWH_ddq{i} * qp + cross(v_angp_dq{i},v_fr_cWH) + cross(v_angp,v_fr_cWH_dq(:,i))   );
                e_rl_cWHp_dq{i} = e_posp_dq{i} + e_T_v_dq{i} * ( v_rl_cWH_dq     * qp + cross(v_angp,v_rl_cWH) ) ...
                                               + e_T_v       * ( v_rl_cWH_ddq{i} * qp + cross(v_angp_dq{i},v_rl_cWH) + cross(v_angp,v_rl_cWH_dq(:,i))   );
                e_rr_cWHp_dq{i} = e_posp_dq{i} + e_T_v_dq{i} * ( v_rr_cWH_dq     * qp + cross(v_angp,v_rr_cWH) ) ...
                                               + e_T_v       * ( v_rr_cWH_ddq{i} * qp + cross(v_angp_dq{i},v_rr_cWH) + cross(v_angp,v_rr_cWH_dq(:,i))   );
                e_CGp_dq{i}     = e_posp_dq{i} + e_T_v_dq{i} * (                        cross(v_angp,obj.DNA.CG_s)         ) ...
                                               + e_T_v       * (                        cross(v_angp_dq{i},obj.DNA.CG_s)   );

                Ekin_dq{i} = 0.5 * e_fl_cWHp_dq{i}.'  * obj.sus_fl.unsMass * e_fl_cWHp  + 0.5 * e_fl_cWHp.' * obj.sus_fl.unsMass * e_fl_cWHp_dq{i}  + ...
                             0.5 * e_fr_cWHp_dq{i}.'  * obj.sus_fr.unsMass * e_fr_cWHp  + 0.5 * e_fr_cWHp.' * obj.sus_fr.unsMass * e_fr_cWHp_dq{i}  + ...
                             0.5 * e_rl_cWHp_dq{i}.'  * obj.sus_rl.unsMass * e_rl_cWHp  + 0.5 * e_rl_cWHp.' * obj.sus_rl.unsMass * e_rl_cWHp_dq{i}  + ...
                             0.5 * e_rr_cWHp_dq{i}.'  * obj.sus_rr.unsMass * e_rr_cWHp  + 0.5 * e_rr_cWHp.' * obj.sus_rr.unsMass * e_rr_cWHp_dq{i}  + ...
                             0.5 * e_CGp_dq{i}.'      * obj.DNA.susMass * e_CGp      + 0.5 * e_CGp.'     * obj.DNA.susMass * e_CGp_dq{i}      + ...
                      diff2( 0.5 * v_angp.' * e_T_v.' * obj.DNA.CG_I    * e_T_v * v_angp , q(i));


                e_fl_cWHp_dqp{i} = e_posp_dqp{i} + e_T_v * ( v_fl_cWH_dq(:,i) + cross(v_angp_dqp{i},v_fl_cWH) );
                e_fr_cWHp_dqp{i} = e_posp_dqp{i} + e_T_v * ( v_fr_cWH_dq(:,i) + cross(v_angp_dqp{i},v_fr_cWH) );
                e_rl_cWHp_dqp{i} = e_posp_dqp{i} + e_T_v * ( v_rl_cWH_dq(:,i) + cross(v_angp_dqp{i},v_rl_cWH) );
                e_rr_cWHp_dqp{i} = e_posp_dqp{i} + e_T_v * ( v_rr_cWH_dq(:,i) + cross(v_angp_dqp{i},v_rr_cWH) );    
                e_CGp_dqp{i}     = e_posp_dqp{i} + e_T_v * (                    cross(v_angp_dqp{i},obj.DNA.CG_s) );

                e_fl_cWHp_dqp_p{i} =  e_T_v_p * ( v_fl_cWH_dq(:,i)   + cross(v_angp_dqp{i},v_fl_cWH) ) + ...
                                      e_T_v   * ( v_fl_cWH_dq_p(:,i) + cross(v_angp_dqp_p{i},v_fl_cWH) + cross(v_angp_dqp{i},v_fl_cWH_dq*qp) );
                e_fr_cWHp_dqp_p{i} =  e_T_v_p * ( v_fr_cWH_dq(:,i)   + cross(v_angp_dqp{i},v_fr_cWH) ) + ...
                                      e_T_v   * ( v_fr_cWH_dq_p(:,i) + cross(v_angp_dqp_p{i},v_fr_cWH) + cross(v_angp_dqp{i},v_fr_cWH_dq*qp) );
                e_rl_cWHp_dqp_p{i} =  e_T_v_p * ( v_rl_cWH_dq(:,i)   + cross(v_angp_dqp{i},v_rl_cWH) ) + ...
                                      e_T_v   * ( v_rl_cWH_dq_p(:,i) + cross(v_angp_dqp_p{i},v_rl_cWH) + cross(v_angp_dqp{i},v_rl_cWH_dq*qp) );
                e_rr_cWHp_dqp_p{i} =  e_T_v_p * ( v_rr_cWH_dq(:,i)   + cross(v_angp_dqp{i},v_rr_cWH) ) + ...
                                      e_T_v   * ( v_rr_cWH_dq_p(:,i) + cross(v_angp_dqp_p{i},v_rr_cWH) + cross(v_angp_dqp{i},v_rr_cWH_dq*qp) );
                e_CGp_dqp_p{i}     =  e_T_v_p * (                      cross(v_angp_dqp{i},obj.DNA.CG_s) ) + ...
                                      e_T_v   * (                      cross(v_angp_dqp_p{i},obj.DNA.CG_s)                                       );


                Ekin_dqp_dt{i} = 0.5 * obj.sus_fl.unsMass * ( e_fl_cWHp_dqp_p{i}.' * e_fl_cWHp + e_fl_cWHp_dqp{i}.' * e_fl_cWHpp   + ...
                                                              e_fl_cWHpp.' * e_fl_cWHp_dqp{i}  + e_fl_cWHp.' * e_fl_cWHp_dqp_p{i}) + ...
                                 0.5 * obj.sus_fr.unsMass * ( e_fr_cWHp_dqp_p{i}.' * e_fr_cWHp + e_fr_cWHp_dqp{i}.' * e_fr_cWHpp   + ...
                                                              e_fr_cWHpp.' * e_fr_cWHp_dqp{i}  + e_fr_cWHp.' * e_fr_cWHp_dqp_p{i}) + ... 
                                 0.5 * obj.sus_rl.unsMass * ( e_rl_cWHp_dqp_p{i}.' * e_rl_cWHp + e_rl_cWHp_dqp{i}.' * e_rl_cWHpp   + ...
                                                              e_rl_cWHpp.' * e_rl_cWHp_dqp{i}  + e_rl_cWHp.' * e_rl_cWHp_dqp_p{i}) + ... 
                                 0.5 * obj.sus_rr.unsMass * ( e_rr_cWHp_dqp_p{i}.' * e_rr_cWHp + e_rr_cWHp_dqp{i}.' * e_rr_cWHpp   + ...
                                                              e_rr_cWHpp.' * e_rr_cWHp_dqp{i}  + e_rr_cWHp.' * e_rr_cWHp_dqp_p{i}) + ...
                                 0.5 * obj.DNA.susMass    * ( e_CGp_dqp_p{i}.' * e_CGp + e_CGp_dqp{i}.' * e_CGpp                   + ...
                                                              e_CGpp.' * e_CGp_dqp{i}  + e_CGp.' * e_CGp_dqp_p{i}                  + ...
                                 0.5 * ( e_angp_dqp_p{i}.' * obj.DNA.CG_I * e_angp + e_angp_dqp{i}.' * obj.DNA.CG_I * e_angpp      + ...
                                         e_angpp.' * obj.DNA.CG_I * e_angp_dqp{i}  + e_angp.' * obj.DNA.CG_I * e_angp_dqp_p{i} )      ) ;

                % AINDA FALTA A FORCA DO AMORTECEDOR
                Qnc{i} =  (e_pos_dq{i}   +   e_T_v_dq{i} * v_fl_TYRO   +   e_T_v * v_fl_TYRO_dq(:,i)).' * e_fl_Ftyre + ...
                          (e_pos_dq{i}   +   e_T_v_dq{i} * v_fr_TYRO   +   e_T_v * v_fr_TYRO_dq(:,i)).' * e_fr_Ftyre + ...
                          (e_pos_dq{i}   +   e_T_v_dq{i} * v_rl_TYRO   +   e_T_v * v_rl_TYRO_dq(:,i)).' * e_rl_Ftyre + ...
                          (e_pos_dq{i}   +   e_T_v_dq{i} * v_rr_TYRO   +   e_T_v * v_rr_TYRO_dq(:,i)).' * e_rr_Ftyre ;

                % Epot_dqp_dt{i}=0  because Epot must not depend on velocities 
                Lagrange{i} = Ekin_dqp_dt{i} - Ekin_dq{i} + Epot_dq{i} - Qnc{i};

                parfor_progress(pwd);
            end
            parfor_progress(pwd,0);

            clearvars -except t q obj Lagrange e_fl_Ftyre e_fr_Ftyre e_rl_Ftyre e_rr_Ftyre


            % OBJETIVO FINAL EH CRIAR UMA FUNCAO:  equal0 = f(q,qp,s) PARA USAR COM
            % ODEFUN ... substituir qpp por qp_new e adicionar nas equacoes (qp - q_new = 0)
            % function dy = odefun(t,y,yp)
            % dy = zeros(2,1);
            % dy(1) = yp(1)-y(2);
            % dy(2) = yp(2)+1;
            %end
            
%             test{1} =  e_ang_phi;
%             test{2} = e_ang_teta;
%             test_sym = [test{:}];
%             test_sym = test_sym(t).';
            
            fprintf('    Grouping differential equations...\n'); tic;
            Lagrange_sym = [Lagrange{:}];
            Lagrange_sym = Lagrange_sym(t).';
            toc
        
%             [eqsBlocks,varsBlocks] = findDecoupledBlocks(Lagrange_sym,q)           
                     
            fprintf('    Reducing differential order...\n');
            [Lagrange_red, q_red] = reduceDifferentialOrder(Lagrange_sym, q);
            
            fprintf('    Calculating mass matrix...\n');
            [Mred,Fred] = massMatrixForm(Lagrange_red,q_red);
            
            s_fl = cell2sym(struct2cell(obj.sus_fl.s_prefix));
            s_fr = cell2sym(struct2cell(obj.sus_fr.s_prefix));
            s_rl = cell2sym(struct2cell(obj.sus_rl.s_prefix));
            s_rr = cell2sym(struct2cell(obj.sus_rr.s_prefix));
            
            
            % PREPARE DAEs WITH MATLAB FUNCTION (SEE massMatrixForm doc)...
            % CALL ode15s instead of ode15i, as now we are dealing with
            % explicit DAEs
            
            
            
%             fprintf('    Reducing differential index...\n'); tic;
%             isLowIndex = isLowIndexDAE(Lagrange_red,q_red)

%             obj.dynamicFunction = daeFunction(Lagrange_red, q_red, ...
%                                               e_fl_Ftyre,e_fr_Ftyre,e_rl_Ftyre,e_rr_Ftyre,...
%                                               s_fl,s_fr,s_rl,s_rr,...
%                                               'File', 'veh_DAEfun');           
%             tspan = [0 4*logspace(-6,6)];
%             [t,y] = ode15i(@(t,y,yp) veh.dynamicFunction(t,y,yp, q, s_fl, s_fr, s_rl, s_rr, e_fl_Ftyre, e_fr_Ftyre, e_rl_Ftyre, e_rr_Ftyre) ,tspan,y0,yp0,options);
            
        end
        
        function plotKinAnalysis(obj)                        
            fun_linspace = @(x,n) linspace(min(min(obj.sus_fl.kin.(x))),max(max(obj.sus_fl.kin.(x))),n);                        
            figure('Color','white')
            
            % SR vs toe vs heave
            subplot(2,2,1)           
            heave = fun_linspace('heave',5);
            bSR   = fun_linspace('bSR',50);            
            toe = obj.sus_fl.getKin('heave',heave,'bSR',bSR,{'toe'});            
            plot(bSR, toe*180/pi)
            grid on
            xlabel('Steering Rack Displacement [m]')
            ylabel('Toe-in angle [degree]')
            leg = legend(cellfun(@(x) sprintf('%.3f',x),num2cell(heave),'UniformOutput',false));
            title(leg,'Heave [m]')
            title('Steering Rack x Toe-in x Heave')
            
            % (spring,motion ratio) vs heave
            subplot(2,2,2) 
            heave = fun_linspace('heave',100);
            bSR = 0;
            bDS = obj.sus_fl.getKin('heave',heave,'bSR',bSR,{'bDS'}) - norm(obj.sus_fl.init.DSCH_0-obj.sus_fl.init.RKDS_0);          
            MR  = gradient(bDS,heave);
            yyaxis left
            plot(heave, bDS)
            ylabel('Spring Displacement [m]')
            yyaxis right
            plot(heave, MR)
            ylabel('Motion Ratio')
            ylim([-1.5 -0.5])
            grid on
            xlabel('Heave [m]')
            title('Heave x Spring Displacement/MotionRatio')
                        
            % heave vs camber vs toe
            subplot(2,2,3) 
            heave = fun_linspace('heave',50);
            toe   = fun_linspace('toe',10);            
            camber = obj.sus_fl.getKin('heave',heave,'toe',toe,{'camber'});            
            plot(heave, camber*180/pi)
            grid on
            xlabel('Heave [m]')
            ylabel('Camber [degree]')
            leg = legend(cellfun(@(x) sprintf('%.1f',x),num2cell(toe*180/pi),'UniformOutput',false));
            title(leg,'Toe-in angle [degree]');
            title('Heave x Camber x Toe-in')
           
            % heave vs camber vs toe
            subplot(2,2,4) 
            contour(obj.sus_fl.kin.heave,obj.sus_fl.kin.camber*180/pi,obj.sus_fl.kin.toe*180/pi,'ShowText','on')
            grid on
            xlabel('Heave [m]')
            ylabel('Camber [degree]')
            title('Heave x Camber x Toe-in')
            
            % track vs heave vs toe
        end
        
        function plot(obj, fig, init)
            % Plot all the four suspensions... add parameter to the
            % suspension plot function (CG pos and rotation)
        end
    end
    
end

