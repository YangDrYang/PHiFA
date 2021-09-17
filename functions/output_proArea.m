function [epoch,laser_projArea,srp_projArea,alb_projArea,ir_projArea,drag_projArea] = output_proArea(is, sim_out)

epoch(1,length(sim_out)-1) = 0;
laser_projArea(1,length(sim_out)-1) = 0;
srp_projArea(1,length(sim_out)-1) = 0;
alb_projArea(1,length(sim_out)-1) = 0;
ir_projArea(1,length(sim_out)-1) = 0;
drag_projArea(1,length(sim_out)-1) = 0;

for i = 2:length(sim_out)
    if i~=1 && sim_out(i).t==0
        epoch = epoch(:,1:i-2);
        laser_projArea = laser_projArea(:,1:i-2);
        drag_projArea = drag_projArea(:,1:i-2);
        srp_projArea = srp_projArea(:,1:i-2);
        alb_projArea = alb_projArea(:,1:i-2);
        ir_projArea = ir_projArea(:,1:i-2);
        
        break;
    end
    epoch(i-1) = sim_out(i).t(end);
    if is.a_lase
        laser_projArea(1,i-1) = sim_out(i).laser.projarea;
    else
        laser_projArea(1,i-1) = 0;
    end
    if is.a_drag
        drag_projArea(1,i-1) = sim_out(i).drag.projarea;
    else
        drag_projArea(1,i-1) = 0;
    end
    if is.a_srad
        srp_projArea(1,i-1) = sim_out(i).srad.projarea;
    else
        srp_projArea(1,i-1) = 0;
    end
    
    if is.a_erad
        alb_projArea(1,i-1) = sim_out(i).albedo.projarea;
        ir_projArea(1,i-1) = sim_out(i).infrared.projarea;
    else
        alb_projArea(1,i-1) = 0;
        ir_projArea(1,i-1) = 0;
    end

end