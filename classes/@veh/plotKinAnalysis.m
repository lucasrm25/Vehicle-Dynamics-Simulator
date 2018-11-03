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