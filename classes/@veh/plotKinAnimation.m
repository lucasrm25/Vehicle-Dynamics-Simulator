function plotKinAnimation(obj)
    % Plot all the four suspensions... add parameter to the
    % suspension plot function (CG pos and rotation)
    
    movie = false;

    fig = figure('Color','w'); hold on; grid on;
    xlabel('x axis')
    ylabel('y axis')
    zlabel('z axis')

    axis equal
    xlim([-0.2 2])
    ylim([-0.9 0.9])
    zlim([-0.1 0.6])

    % view(3)
    az = 210;
    el = 45;
    view(az, el);

    bSR_slider = 0;
    FLAUR_3_slider = 0;
    RLAUR_3_slider = 0;
    hscrollbar_bSR    = uicontrol('style','slider','units','normalized','position',[0.05 0 0.90 .05],'callback', @(src,evt) assignin('base','bSR_slider',src.Value) );
    hscrollbar_FLAUR  = uicontrol('style','slider','units','normalized','position',[0 .05 .05 0.95],'callback',@(src,evt) assignin('base','FLAUR_3_slider',src.Value) );
    hscrollbar_RLAUR  = uicontrol('style','slider','units','normalized','position',[0.95 .05 .05 0.95],'callback',@(src,evt) assignin('base','RLAUR_3_slider',src.Value) );


    WH_stroke = [-2 2]*0.0254;
    SR_stroke = [-0.02 0.02];

    loops = 40;
    F(loops) = struct('cdata',[],'colormap',[]);
    init=true;
    while true
    for i = 1:loops
        FXq = obj.sus_fl.init.LAUR_0(3)+WH_stroke(1) + (WH_stroke(2)-WH_stroke(1))*(sin(i/loops*2*pi)/2+0.5); %*FLAUR_3_slider;
        FYq = -(SR_stroke(1) + (SR_stroke(2)-SR_stroke(1))*bSR_slider);
        RXq = obj.sus_rl.init.LAUR_0(3)+WH_stroke(1) + (WH_stroke(2)-WH_stroke(1))*(sin(i/loops*2*pi+pi/2)/2+0.5); %*RLAUR_3_slider;
        RYq = 0;

        q_tm1 = [0 0 0   0 0 0  FXq FXq RXq RXq   FYq RYq ];

        obj.sus_fl.q_val = q_tm1([7 11])';
        obj.sus_fr.q_val = q_tm1([8 11])';
        obj.sus_rl.q_val = q_tm1([9 12])';
        obj.sus_rr.q_val = q_tm1([10 12])';

        obj.sus_fl.plot(fig,init)
        obj.sus_fr.plot(fig,init)
        obj.sus_rl.plot(fig,init)
        obj.sus_rr.plot(fig,init)

    %     axis vis3d
        drawnow();
        if movie
            F(i) = getframe(gcf);
        end
        init = false;
    end
    end

    if movie
        fig = figure;
        movie(fig,F,2);
    end
end
