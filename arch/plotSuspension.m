function suspensionHandles = plotSuspension( fig, suspensionHandlesIN, init, tyr, Fsus, XFieldname ,Xq, YFieldname,Yq )

    LAUR_val = interpStructArray(Fsus.kin,XFieldname,Xq,YFieldname,Yq,{'LAUR_1','LAUR_2','LAUR_3'});
    UAUR_val = interpStructArray(Fsus.kin,XFieldname,Xq,YFieldname,Yq,{'UAUR_1','UAUR_2','UAUR_3'});
    SRTR_val = interpStructArray(Fsus.kin,XFieldname,Xq,YFieldname,Yq,{'SRTR_1','SRTR_2','SRTR_3'});
    TRUR_val = interpStructArray(Fsus.kin,XFieldname,Xq,YFieldname,Yq,{'TRUR_1','TRUR_2','TRUR_3'});
    uWH_val  = interpStructArray(Fsus.kin,XFieldname,Xq,YFieldname,Yq,{'uWH_1','uWH_2','uWH_3'});
    cWH_val  = interpStructArray(Fsus.kin,XFieldname,Xq,YFieldname,Yq,{'cWH_1','cWH_2','cWH_3'});
    TYRO_val = interpStructArray(Fsus.kin,XFieldname,Xq,YFieldname,Yq,{'TYRO_1','TYRO_2','TYRO_3'});
    RKPR_val = interpStructArray(Fsus.kin,XFieldname,Xq,YFieldname,Yq,{'RKPR_1','RKPR_2','RKPR_3'});
    PRUA_val = interpStructArray(Fsus.kin,XFieldname,Xq,YFieldname,Yq,{'PRUA_1','PRUA_2','PRUA_3'});
    RKDS_val = interpStructArray(Fsus.kin,XFieldname,Xq,YFieldname,Yq,{'RKDS_1','RKDS_2','RKDS_3'});
    
    RKAL_val = interpStructArray(Fsus.kin,XFieldname,Xq,YFieldname,Yq,{'RKAL_1','RKAL_2','RKAL_3'});
    ALAB_val = interpStructArray(Fsus.kin,XFieldname,Xq,YFieldname,Yq,{'ALAB_1','ALAB_2','ALAB_3'});
    
    figure(fig)
  
    if init
        hold on
        % Suspension
        vec = [Fsus.init.LACH1_0 LAUR_val Fsus.init.LACH2_0]; pLA = plot3( vec(1,:),vec(2,:),vec(3,:), 'k', 'LineWidth', 5, 'MarkerSize',20,'Marker','.','MarkerEdgeColor','k');
        vec = [Fsus.init.UACH1_0 UAUR_val Fsus.init.UACH2_0]; pUA = plot3( vec(1,:),vec(2,:),vec(3,:), 'k', 'LineWidth', 5, 'MarkerSize',20,'Marker','.');
        % Steering
        vec = [SRTR_val.*[1; 0; 1] SRTR_val TRUR_val]; pSR = plot3( vec(1,:),vec(2,:),vec(3,:), 'b', 'LineWidth', 3, 'MarkerSize',20,'Marker','.','MarkerEdgeColor','k','MarkerFaceColor','k');
        % Upright
        vec = [LAUR_val UAUR_val TRUR_val];  pUR = fill3(vec(1,:),vec(2,:),vec(3,:),'r');
        % Wheel
        vec = plotTyre(uWH_val, cWH_val, tyr.unloadedRadius, tyr.aspectRatio, tyr.width);   pWH = surf(vec(:,:,1), vec(:,:,2), vec(:,:,3), 'FaceAlpha', 0.3, 'EdgeAlpha', 0.5, 'FaceColor', 'k', 'FaceLighting','flat');
        % Tyre-Road Contact
        pTY = scatter3(TYRO_val(1),TYRO_val(2),TYRO_val(3),50,'filled','MarkerFaceColor','r');
        % Push-Rod
        vec = [RKPR_val PRUA_val]; pPR = plot3( vec(1,:),vec(2,:),vec(3,:), 'g', 'LineWidth', 5, 'MarkerSize',15,'Marker','.');
        % Rocker
        vec = [RKDS_val RKPR_val Fsus.init.RKCH_0];  pRK = fill3(vec(1,:),vec(2,:),vec(3,:),'g');
        % Spring/Damper
        vec = plotSpring(RKDS_val, Fsus.init.DSCH_0, 0.02, 6); pSD = plot3(vec(1,:),vec(2,:),vec(3,:), 'Color',[1 .5 0],'LineWidth',5);
        % ARB
        vec = [Fsus.init.ALAR_0.*[1; 0; 1] Fsus.init.ALAR_0 ALAB_val RKAL_val]; pAR = plot3( vec(1,:),vec(2,:),vec(3,:), 'c', 'LineWidth', 5, 'MarkerSize',20,'Marker','.','MarkerEdgeColor','k','MarkerFaceColor','k');
        hold off
        
        suspensionHandles = struct('pLA',pLA,'pUA',pUA,'pSR',pSR,'pUR',pUR,'pWH',pWH,'pTY',pTY,'pPR',pPR,'pRK',pRK,'pSD',pSD,'pAR',pAR);
    else
        % Suspension
        vec = [Fsus.init.LACH1_0 LAUR_val Fsus.init.LACH2_0]; suspensionHandlesIN.pLA.XData = vec(1,:); suspensionHandlesIN.pLA.YData = vec(2,:); suspensionHandlesIN.pLA.ZData = vec(3,:);
        vec = [Fsus.init.UACH1_0 UAUR_val Fsus.init.UACH2_0]; suspensionHandlesIN.pUA.XData = vec(1,:); suspensionHandlesIN.pUA.YData = vec(2,:); suspensionHandlesIN.pUA.ZData = vec(3,:);
        % Steering
        vec = [SRTR_val.*[1; 0; 1] SRTR_val TRUR_val]; suspensionHandlesIN.pSR.XData = vec(1,:); suspensionHandlesIN.pSR.YData = vec(2,:); suspensionHandlesIN.pSR.ZData = vec(3,:);
        % Upright
        vec = [LAUR_val UAUR_val TRUR_val]; suspensionHandlesIN.pUR.Vertices = vec';
        % Wheel
        vec = plotTyre(uWH_val, cWH_val, tyr.unloadedRadius, tyr.aspectRatio, tyr.width); suspensionHandlesIN.pWH.XData=vec(:,:,1); suspensionHandlesIN.pWH.YData=vec(:,:,2); suspensionHandlesIN.pWH.ZData=vec(:,:,3); 
        % Tyre-Road Contact
        suspensionHandlesIN.pTY.XData = TYRO_val(1); suspensionHandlesIN.pTY.YData = TYRO_val(2); suspensionHandlesIN.pTY.ZData = TYRO_val(3);     
        % Push Rod
        vec = [RKPR_val PRUA_val]; suspensionHandlesIN.pPR.XData = vec(1,:); suspensionHandlesIN.pPR.YData = vec(2,:); suspensionHandlesIN.pPR.ZData = vec(3,:);        
        % Rocker
        vec = [RKDS_val RKPR_val Fsus.init.RKCH_0];  suspensionHandlesIN.pRK.Vertices = vec';       
        % Spring/Damper
        vec = plotSpring(RKDS_val, Fsus.init.DSCH_0, 0.02, 6); suspensionHandlesIN.pSD.XData = vec(1,:); suspensionHandlesIN.pSD.YData = vec(2,:); suspensionHandlesIN.pSD.ZData = vec(3,:);
        % ARB
        vec = [Fsus.init.ALAR_0.*[1; 0; 1] Fsus.init.ALAR_0 ALAB_val RKAL_val]; suspensionHandlesIN.pAR.XData = vec(1,:); suspensionHandlesIN.pAR.YData = vec(2,:); suspensionHandlesIN.pAR.ZData = vec(3,:);
        
        suspensionHandles = suspensionHandlesIN;       
    end
end