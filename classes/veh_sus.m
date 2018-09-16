classdef veh_sus < matlab.mixin.Copyable  
    
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
    
    
%% Properties

    properties
        init
        setup
        veh_tyr
        unsMass
        
        q           % main coordinates - vector of symbolic names
        q_val
        s           % secondary coordinates - vector of symbolic names
        s_prefix    % struct containing symbolic names with prefix
        kl          % kinematic loops - class
        kin         % kinematic tables                        
        
        Kmatrix
        Lmatrix
        
        Kmatrix_num
        Lmatrix_num
        
        reshapefun

    end
    properties (GetAccess=private)
        pHandles
    end
    properties (Constant)
    end
    properties (Dependent)
    end
    
%% Methods
    methods(Static)
        function [ structArray ] = syms2structArray( syms, array, dim )
            cell_array = num2cell(array);
            fieldNames = cellfun(@char, sym2cell(syms),'UniformOutput', false);
            structArray = cell2struct(cell_array,fieldNames, dim);
        end
        function [ structArray ] = syms2structArray2( syms, array )
            for i=1:numel(array)
                structArray.(char(syms(i))) = array(i);
            end
        end
    end
    
    methods
        function obj = veh_sus(init, setup, unsMass, veh_tyr, WH_stroke, WH_steps, SR_stroke, SR_steps)
            obj.init    = init;
            obj.setup   = setup;
            obj.veh_tyr = veh_tyr;
            obj.unsMass = unsMass;
            obj.calculateKinematics(WH_stroke, WH_steps, SR_stroke, SR_steps);
        end
        
        function obj = reshape(obj,new_size,new_positions)
            shift = zeros(numel(obj.q),new_size);
            shift(1,new_positions(1)) = 1;
            shift(2,new_positions(2)) = 1;
            obj.reshapefun = @(fun,q,s) vpa(fun(q(new_positions),s))*shift;
        end
        
        function Vq = getKin(obj, XFieldName, Xq, YFieldName, Yq, fieldNames)   
            X = obj.kin.(char(XFieldName));
            Y = obj.kin.(char(YFieldName));
            for i=1:numel(fieldNames)
                V  = obj.kin.(char(fieldNames(i)));
                if numel(Xq) > 1 || numel(Yq) > 1
                    Vq = griddata(X,Y,V,Xq,Yq','cubic');
                    return
                end
                if min(size(V)) == 1
                    Vq(i,1) = interp1(X,V,Xq,'pchip',0);
                else
%                     Vq(i,1) = interp2(X,Y,V,Xq,Yq,'cubic',0);
                    Vq(i,1) = griddata(X,Y,V,Xq,Yq','cubic');
                end
            end
        end

        function Vq = getKinq(obj, q, fieldNames)
            Vq = getKin(obj, obj.q(1), q(1), obj.q(2), q(2), fieldNames);
        end
        
        function obj = set_prefix(obj, prefix)
            if ~isempty(obj.s)
                for i=1:numel(obj.s)
                    obj.s_prefix.(char(obj.s(i))) = sym(strcat(prefix,char(obj.s(i))));
                end
            end
            
        end
        
        function compare_K_analytic_numeric(obj,fieldname)
            multiWaitbar( 'CloseAll' );
            
            s  = obj.kin.(fieldname);
            q1 = obj.kin.(char(obj.q(1)));
            q2 = obj.kin.(char(obj.q(2)));
            % REFINE
            [q1_ref q2_ref] = meshgrid(linspace(min(min(q1)),max(max(q1)),100) , linspace(min(min(q2)),max(max(q2)),100));
            s_num_ref = interp2(q1,q2,s,q1_ref,q2_ref,'cubic');
            % REFINE
            
            % NUMERIC            
            [dsx dsy] = gradient(s_num_ref);
            [dqx dqy] = gradient(q1_ref);
            K1_num_ref = dsx./dqx;
            [dqx dqy] = gradient(q2_ref);
            K2_num_ref = dsy./dqy;
            % NUMERIC

            % ANALYTIC
            for i=1:size(q1,1)
                for j=1:size(q1,2)
                    q = [q1(i,j) q2(i,j)]';
                    s = obj.getKinq(q,obj.s);                    
                    K_anl(i,j,1:2) = obj.Kmatrix.(fieldname)(q,s);
                    multiWaitbar('Comparing K matrices Analytical/Numeric', 'Value', (j+(i-1)*size(q1,2)) / numel(q1), 'Color', 'g' );
                end
            end
            multiWaitbar( 'CloseAll' );
            K1_anl_ref = interp2(q1,q2,K_anl(:,:,1),q1_ref,q2_ref,'cubic');
            K2_anl_ref = interp2(q1,q2,K_anl(:,:,2),q1_ref,q2_ref,'cubic');
            % ANALYTIC
            
            figure
            subplot(1,2,1)
            surf(q1_ref,q2_ref,K1_num_ref,'FaceColor','r') % surf(q1,q2,K1_num);
            hold on
            surf(q1_ref,q2_ref,K1_anl_ref,'FaceColor','b')
            xlabel(char(obj.q(1)))
            ylabel(char(obj.q(2)))
            zlabel(sprintf('d%s/d%s',fieldname,char(obj.q(1))))
            legend({'Numeric','Analytic'})
            subplot(1,2,2)            
            surf(q1_ref,q2_ref,K2_num_ref,'FaceColor','r') % surf(q1,q2,K1_num);
            hold on
            surf(q1_ref,q2_ref,K2_anl_ref,'FaceColor','b')
            xlabel(char(obj.q(1)))
            ylabel(char(obj.q(2)))
            zlabel(sprintf('d%s/d%s',fieldname,char(obj.q(2))))
%             legend({'Numeric','Analytic'})
        end
        
        function calculateKinematics(obj,WH_stroke, WH_steps, SR_stroke, SR_steps)

            multiWaitbar( 'CloseAll' );
            digits(4)
            
            bSR  = sym('bSR');
            LAUR = sym('LAUR_%d', [3 1]);

            % Main Coordinates
            obj.q = [LAUR(3); 
                     bSR];

            %% Front left suspension kinematics

            uUA = obj.init.UACH1_0 - obj.init.UACH2_0;         
            uLA = obj.init.LACH1_0 - obj.init.LACH2_0;
            cFUA = (obj.init.UACH1_0 + obj.init.UACH2_0)/2;    
            cFLA = (obj.init.LACH1_0 + obj.init.LACH2_0)/2;

            vUA = cFUA - obj.init.UAUR_0;          
            vUR = obj.init.UAUR_0 - obj.init.LAUR_0;       
            vLA = obj.init.LAUR_0 - cFLA;         
            vCH = cFLA - cFUA;

            % New coordinate variables to be solved (8 Variables)
            UAUR = sym('UAUR_%d', [3 1]);
            bUA  = sym('bUA');
            bLA  = sym('bLA');

            % Loop vector functions - (8 Variables, 7 Equations = 1 main CO, 7 secondary COs)
            f_S = [ rotateVector(vLA,uLA,bLA) + (UAUR-LAUR) + rotateVector(vUA,uUA,bUA) + vCH;
                    rotatePoint(obj.init.UAUR_0,uUA,bUA,cFUA) - UAUR ;
                    norm(UAUR-LAUR) - norm(obj.init.UAUR_0-obj.init.LAUR_0)    ];
            f_S = simplify(f_S);            
            s_S   = [LAUR(1:2);   UAUR;   bUA; bLA];
            s_S_0 = [obj.init.LAUR_0(1:2); obj.init.UAUR_0; 0;    0];            
            q_S = [LAUR(3)];
            kinLoop_S = kinLoop(f_S, q_S, s_S, s_S_0);            

            %% Steering system kinematics

            uSR = [0 1 0]';

            % New coordinate variables to be solved (7 Variables)
            SRTR = sym('SRTR_%d', [3 1]);
            TRUR = sym('TRUR_%d', [3 1]);

            % Loop vector functions - (13 Variables, 6 Equations = 7 main CO, 6 secondary COs)
            f_SS = [ obj.init.SRTR_0 + bSR*uSR - SRTR;
                     norm(TRUR-SRTR) - norm(obj.init.TRUR_0-obj.init.SRTR_0);
                     norm(LAUR-TRUR) - norm(obj.init.LAUR_0-obj.init.TRUR_0);
                     norm(UAUR-TRUR) - norm(obj.init.UAUR_0-obj.init.TRUR_0)   ];
            f_SS = simplify(f_SS);
            s_SS   = [SRTR; TRUR];
            s_SS_0 = [obj.init.SRTR_0; obj.init.TRUR_0];
            q_SS = [LAUR; UAUR; bSR];
            kinLoop_SS = kinLoop(f_SS, q_SS, s_SS, s_SS_0 );           

            %% Front Wheel kinematics

            % New coordinate variables to be solved (6 Variables)
            uWH = sym('uWH_%d', [3 1]);     % Wheel rotation vector
            cWH = sym('cWH_%d', [3 1]);


            % c0 = [DNA.DIN_CO(1) DNA.track/2+obj.setup.camber*obj.tyre.unloadedRadius DNA.DIN_CO(3)]';
            c0 = obj.init.cWH_0;
            v0 = rotateVector(rotateVector([0 1 0]',[1 0 0]',-obj.setup.camber),[0 0 1]',-obj.setup.toe);                       

            f_WH =[ createTmatrixFunction_point(LAUR,TRUR,UAUR,obj.init.LAUR_0,obj.init.TRUR_0,obj.init.UAUR_0,c0) ;
                    createTmatrixFunction_vector(LAUR,TRUR,UAUR,obj.init.LAUR_0,obj.init.TRUR_0,obj.init.UAUR_0,v0) ];
            s_WH   = [cWH; uWH];
            q_WH = [LAUR; UAUR; TRUR];
            kinLoop_WH = kinLoop(f_WH, q_WH, s_WH, 0, 'fun');           

            %% Front tyre contact point

            % New coordinate variables to be solved (7 Variables)
            TYRO = sym('TYRO_%d', [3 1]);

            z = [0 0 1]';
            d_vector = cross(uWH,cross(-z,uWH));

            f_TY = d_vector/norm(d_vector)*obj.veh_tyr.unloadedRadius + cWH ;
            s_TY = [TYRO];
            q_TY = [cWH; uWH];            
            kinLoop_TY = kinLoop(f_TY, q_TY, s_TY, 0, 'fun');            

            %% Push Rod kinematics

            PRUA = sym('PRUA_%d', [3 1]);

            f_PR = createTmatrixFunction_point(obj.init.UACH1_0,obj.init.UACH2_0,UAUR,obj.init.UACH1_0,obj.init.UACH2_0,obj.init.UAUR_0,obj.init.PRUA_0) ;
            s_PR = PRUA;
            q_PR = UAUR;            
            kinLoop_PR = kinLoop(f_PR, q_PR, s_PR, 0, 'fun');                        

            %% Rocker kinematics

            RKPR = sym('RKPR_%d', [3 1]);
            syms bRK

            uRK = cross((obj.init.RKPR_0-obj.init.RKCH_0),(obj.init.RKDS_0-obj.init.RKCH_0));

            f_RK = [ rotatePoint(obj.init.RKPR_0,uRK,bRK,obj.init.RKCH_0) - RKPR;
                     norm(RKPR-PRUA)-norm(obj.init.RKPR_0-obj.init.PRUA_0) ];
            f_RK = simplify(f_RK);
            s_RK   = [RKPR; bRK];
            s_RK_0 = [obj.init.RKPR_0; 0];
            q_RK   = [PRUA];            
            kinLoop_RK = kinLoop(f_RK, q_RK, s_RK, s_RK_0);           

            %% Damper/Spring Kinematics

            RKDS = sym('RKDS_%d', [3 1]);
            syms bDS                            % Damper/Spring length

            f_DS = [ rotatePoint(obj.init.RKDS_0,uRK,bRK,obj.init.RKCH_0);
                     norm(rotatePoint(obj.init.RKDS_0,uRK,bRK,obj.init.RKCH_0) - obj.init.DSCH_0) ];
            s_DS   = [RKDS; bDS];
            q_DS   = bRK;           
            kinLoop_DS = kinLoop(f_DS, q_DS, s_DS, 0, 'fun');

            %% ARB Kinematics

            uAR = [0 1 0]';

            RKAL = sym('RKAL_%d', [3 1]);
            ALAB = sym('ALAB_%d', [3 1]);
            bAR = sym('bAR'); % syms bAR

            f_AR = [ rotatePoint(obj.init.RKAL_0,uRK,bRK,obj.init.RKCH_0) - RKAL;
                     norm(ALAB-RKAL) - norm(obj.init.ALAB_0-obj.init.RKAL_0);
                     rotatePoint(obj.init.ALAB_0,uAR,bAR,obj.init.ALAR_0) - ALAB ];
            s_AR   = [RKAL; ALAB; bAR];
            s_AR_0 = [obj.init.RKAL_0; obj.init.ALAB_0; 0];
            q_AR = [bRK];           
            kinLoop_AR = kinLoop(f_AR, q_AR, s_AR, s_AR_0);   
            
            %% Toe and Camber

            syms toe camber heave           
            obj.init.TYRO_0 = kinLoop_TY.solve(kinLoop_WH.solve([obj.init.LAUR_0; obj.init.UAUR_0; obj.init.TRUR_0]));            
            f_TC = [ angle2vectors([1 1 0]'.*uWH,[0 1 0]',[0 0 1]')  ;
                     angle2vectors([0 1 1]'.*uWH,[0 1 0]',[1 0 0]')  ;
                     TYRO(3) - obj.init.TYRO_0(3) ];
            s_TC = [toe; camber; heave];
            q_TC = [uWH; TYRO(3)];           
            kinLoop_TC = kinLoop(f_TC, q_TC, s_TC, 0, 'fun');
            
            %% Jacobian Matrices
            
            obj.s  = [s_S ; s_SS ; s_WH ; s_TY ; s_PR ; s_RK ; s_DS ; s_AR; s_TC];
            obj.kl = [kinLoop_S ; kinLoop_SS ; kinLoop_WH ; kinLoop_TY ; kinLoop_PR ; kinLoop_RK ; kinLoop_DS ; kinLoop_AR; kinLoop_TC];
            
            for i=1:numel(obj.kl)
                jacob_K_Kmain = jacobian(obj.kl(i).q,obj.q);        % K matrix definition
                for j=1:i-1
                    jacob_K_Kmain = jacob_K_Kmain + jacobian(obj.kl(i).q,obj.kl(j).s) * obj.kl(j).Kmain;    % 
                end
                obj.kl(i).Kmain = obj.kl(i).jacobian * jacob_K_Kmain;
                for j=1:numel(obj.kl(i).s)
                    obj.Kmatrix.(char(obj.kl(i).s(j))) = matlabFunction( obj.kl(i).Kmain(j,:) ,'Vars',{obj.q,obj.s});
                    obj.Lmatrix.(char(obj.kl(i).s(j))) = matlabFunction( jacobian(obj.kl(i).Kmain(j,:),obj.q) ,'Vars',{obj.q,obj.s});
                    multiWaitbar('Simplifying K and L Matrices', 'Increment', 1/numel(obj.s), 'Color', 'g' );                               
                end
            end

            
            %% Calculate secundary coordinates
                        
            [LAUR_3_val,bSR_val] = meshgrid(linspace(obj.init.LAUR_0(3)+WH_stroke(1),obj.init.LAUR_0(3)+WH_stroke(2),WH_steps),...
                                            linspace(SR_stroke(1),SR_stroke(2),SR_steps));
            names = [obj.q; obj.s];
            val   = cat(3,LAUR_3_val,bSR_val, zeros(size(LAUR_3_val,1),size(LAUR_3_val,2),numel(obj.s)) );
            for k = 1:numel(obj.kl)
                [~, idxq] = ismember(obj.kl(k).q, names);                
                [~, idxs] = ismember(obj.kl(k).s, names);
                for i=1:size(val,1)
                    multiWaitbar('Mapping Suspension Kinematics', 'Value', (i+(k-1)*size(val,1)) / (numel(obj.kl)*size(val,1)), 'Color', 'g' );
                    for j=1:size(val,2)
                        q_val = squeeze(val(i,j,idxq(idxq>0)));
                        val(i,j,idxs) = obj.kl(k).solve(q_val);
                    end
                end
            end
            for i=1:numel(names)
                obj.kin.(char(names(i))) = val(:,:,i);
            end
%             multiWaitbar( 'CloseAll' );
        end
            
        function new = mirror(obj)
            new = copy(obj);
            variables2inverse = {'bSR','LAUR_2','UAUR_2','cWH_2','PRUA_2','RKDS_2','RKPR_2','SRTR_2','TRUR_2','TYRO_2','uWH_2','ALAB_2','RKAL_2'};                                    
            for i=1:numel(variables2inverse)
                new.kin.(variables2inverse{i})  = -obj.kin.(variables2inverse{i});
                if isfield(new.Kmatrix,variables2inverse{i})
                    aux = str2func(func2str(new.Kmatrix.(variables2inverse{i})));
                    new.Kmatrix.(variables2inverse{i}) = @(q,s) -aux(q,s);
                    aux = str2func(func2str(new.Lmatrix.(variables2inverse{i})));
                    new.Lmatrix.(variables2inverse{i}) = @(q,s) -aux(q,s);
                end
            end
            init2inverse = {'UACH1_0','UACH2_0','LACH1_0','LACH2_0','DSCH_0','RKCH_0','ALAR_0'};
            for i=1:numel(init2inverse)
                new.init.(init2inverse{i}) = new.init.(init2inverse{i}).*[1,-1,1]';
            end
        end
        
        function vec = plotSpring(~, p1, p2, radius, nCoils )
            offset = 0.033;
            springSize = norm(p2-p1)-offset*2;
            teta = 0:(pi/30):2*pi*nCoils;
            x = radius*sin(teta);
            y = radius*cos(teta);
            z = teta/(2*pi*nCoils)*springSize + offset;

            x = [0 0 x 0 0];
            y = [0 0 y 0 0];
            z = [0 offset z offset+springSize 2*offset+springSize];

            % Rotate
            direction = p2-p1;
            uInitial  = [0 0 1]';
            rotation_vector = cross(uInitial,direction);
            rotation_angle  = acos(dot(uInitial,direction)/(norm(uInitial)*norm(direction)));

            newPoint = rotatePoint([x(:),y(:),z(:)]', rotation_vector, rotation_angle, [0 0 0]');
            x(:) = newPoint(1,:); y(:) = newPoint(2,:); z(:) = newPoint(3,:);

            % Translate
            x = x + p1(1);
            y = y + p1(2);
            z = z + p1(3);

            vec = [x;y;z];
        end
        
        function plot(obj, fig, init)
            plt = @(x) obj.getKinq(obj.q_val, {[x '_1'];[x '_2'];[x '_3']} );

            LAUR_val = plt('LAUR');
            UAUR_val = plt('UAUR');
            SRTR_val = plt('SRTR');
            TRUR_val = plt('TRUR');
            uWH_val  = plt('uWH');
            cWH_val  = plt('cWH');
            TYRO_val = plt('TYRO');
            RKPR_val = plt('RKPR');
            PRUA_val = plt('PRUA');
            RKDS_val = plt('RKDS');
            RKAL_val = plt('RKAL');
            ALAB_val = plt('ALAB');

            figure(fig)

            if init
                hold on
                % Suspension
                vec = [obj.init.LACH1_0 LAUR_val obj.init.LACH2_0]; pLA = plot3( vec(1,:),vec(2,:),vec(3,:), 'k', 'LineWidth', 5, 'MarkerSize',20,'Marker','.','MarkerEdgeColor','k');
                vec = [obj.init.UACH1_0 UAUR_val obj.init.UACH2_0]; pUA = plot3( vec(1,:),vec(2,:),vec(3,:), 'k', 'LineWidth', 5, 'MarkerSize',20,'Marker','.');
                % Steering
                vec = [SRTR_val.*[1; 0; 1] SRTR_val TRUR_val]; pSR = plot3( vec(1,:),vec(2,:),vec(3,:), 'b', 'LineWidth', 3, 'MarkerSize',20,'Marker','.','MarkerEdgeColor','k','MarkerFaceColor','k');
                % Upright
                vec = [LAUR_val UAUR_val TRUR_val];  pUR = fill3(vec(1,:),vec(2,:),vec(3,:),'r');
                % Wheel
                vec = obj.veh_tyr.plot(uWH_val, cWH_val);   pWH = surf(vec(:,:,1), vec(:,:,2), vec(:,:,3), 'FaceAlpha', 0.3, 'EdgeAlpha', 0.5, 'FaceColor', 'k', 'FaceLighting','flat');
                % Tyre-Road Contact
                pTY = scatter3(TYRO_val(1),TYRO_val(2),TYRO_val(3),50,'filled','MarkerFaceColor','r');
                % Push-Rod
                vec = [RKPR_val PRUA_val]; pPR = plot3( vec(1,:),vec(2,:),vec(3,:), 'g', 'LineWidth', 5, 'MarkerSize',15,'Marker','.');
                % Rocker
                vec = [RKDS_val RKPR_val obj.init.RKCH_0];  pRK = fill3(vec(1,:),vec(2,:),vec(3,:),'g');
                % Spring/Damper
                vec = obj.plotSpring(RKDS_val, obj.init.DSCH_0, 0.02, 6); pSD = plot3(vec(1,:),vec(2,:),vec(3,:), 'Color',[1 .5 0],'LineWidth',5);
                % ARB
                vec = [obj.init.ALAR_0.*[1; 0; 1] obj.init.ALAR_0 ALAB_val RKAL_val]; pAR = plot3( vec(1,:),vec(2,:),vec(3,:), 'c', 'LineWidth', 5, 'MarkerSize',20,'Marker','.','MarkerEdgeColor','k','MarkerFaceColor','k');
                hold off

                obj.pHandles = struct('pLA',pLA,'pUA',pUA,'pSR',pSR,'pUR',pUR,'pWH',pWH,'pTY',pTY,'pPR',pPR,'pRK',pRK,'pSD',pSD,'pAR',pAR);
            else
                % Suspension
                vec = [obj.init.LACH1_0 LAUR_val obj.init.LACH2_0]; obj.pHandles.pLA.XData = vec(1,:); obj.pHandles.pLA.YData = vec(2,:); obj.pHandles.pLA.ZData = vec(3,:);
                vec = [obj.init.UACH1_0 UAUR_val obj.init.UACH2_0]; obj.pHandles.pUA.XData = vec(1,:); obj.pHandles.pUA.YData = vec(2,:); obj.pHandles.pUA.ZData = vec(3,:);
                % Steering
                vec = [SRTR_val.*[1; 0; 1] SRTR_val TRUR_val]; obj.pHandles.pSR.XData = vec(1,:); obj.pHandles.pSR.YData = vec(2,:); obj.pHandles.pSR.ZData = vec(3,:);
                % Upright
                vec = [LAUR_val UAUR_val TRUR_val]; obj.pHandles.pUR.Vertices = vec';
                % Wheel
                vec = obj.veh_tyr.plot(uWH_val, cWH_val); obj.pHandles.pWH.XData=vec(:,:,1); obj.pHandles.pWH.YData=vec(:,:,2); obj.pHandles.pWH.ZData=vec(:,:,3); 
                % Tyre-Road Contact
                obj.pHandles.pTY.XData = TYRO_val(1); obj.pHandles.pTY.YData = TYRO_val(2); obj.pHandles.pTY.ZData = TYRO_val(3);     
                % Push Rod
                vec = [RKPR_val PRUA_val]; obj.pHandles.pPR.XData = vec(1,:); obj.pHandles.pPR.YData = vec(2,:); obj.pHandles.pPR.ZData = vec(3,:);        
                % Rocker
                vec = [RKDS_val RKPR_val obj.init.RKCH_0];  obj.pHandles.pRK.Vertices = vec';       
                % Spring/Damper
                vec = obj.plotSpring(RKDS_val, obj.init.DSCH_0, 0.02, 6); obj.pHandles.pSD.XData = vec(1,:); obj.pHandles.pSD.YData = vec(2,:); obj.pHandles.pSD.ZData = vec(3,:);
                % ARB
                vec = [obj.init.ALAR_0.*[1; 0; 1] obj.init.ALAR_0 ALAB_val RKAL_val]; obj.pHandles.pAR.XData = vec(1,:); obj.pHandles.pAR.YData = vec(2,:); obj.pHandles.pAR.ZData = vec(3,:);

            end
        end
    end
    
end

