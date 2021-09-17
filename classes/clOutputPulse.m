classdef clOutputPulse < handle
    properties
        %log variables
        mean_intensity
        max_intensity
        photon_pressure
        fom
        pointing_error
        max_intensity_offset
        fried                   % FRIED's parameter
        rytov                   % log-amplitude variance
        number_screens
        zenith_angle
        station_height
        target_distance
        mean_cm
        mean_area
        mean_aoi
        momentum
        mean_momentum
        amomentum
        mean_angmomentum
        number_hits
        
        time
        step
        lastxv  % laserstation position
        xv      %target object position
        wwdot
        deltat
        
        %miscellenous
        sc = 1;
        nSteps
        bUseLogfile logical = false
        %logfile
        file string
        fileID
    end
    
    methods
        function obj = clOutputPulse()
            
        end
        
        function init(obj, varargin)
            if nargin > 1
                nsteps = varargin{1};
                obj.bUseLogfile = false;
                
                obj.mean_intensity(1:nsteps+1) = 0;
                obj.max_intensity(1:nsteps+1) = 0;
                obj.photon_pressure(1:nsteps+1) = 0;
                obj.fom(1:nsteps+1) = 0;
                obj.pointing_error(1:nsteps+1) = 0;
                obj.max_intensity_offset(1:nsteps+1) = 0;
                obj.fried(1:nsteps+1) = 0;
                obj.rytov(1:nsteps+1) = 0;
                obj.number_screens(1:nsteps+1) = 0;
                obj.zenith_angle(1:nsteps+1) = 0;
                obj.station_height(1:nsteps+1) = 0;
                obj.target_distance(1:nsteps+1) = 0;
                obj.mean_cm(1:nsteps+1) = 0;
                obj.mean_area(1:nsteps+1) = 0;
                obj.mean_aoi(1:nsteps+1) = 0;
                obj.momentum(1:nsteps+1) = 0;
                obj.mean_momentum(1:nsteps+1) = 0;
                obj.amomentum(1:nsteps+1) = 0;
                obj.mean_angmomentum(1:nsteps+1) = 0;
                obj.number_hits(1:nsteps+1) = 0;

                obj.time(1:nsteps+1) = 0;
                obj.step(1:nsteps+1) = 0;
                obj.xv(1:nsteps+1,1:6) = 0;
                obj.lastxv(1:nsteps+1,1:6) = 0;
                obj.wwdot(1:nsteps+1,1:6) = 0;
                obj.deltat(1:nsteps+1) = 0;
            else
                obj.file = sprintf('logfiles\\%s_%s',mfilename(),...
                    datestr(now,'yy-mm-dd_HH-MM-SS'));
                obj.fileID = fopen(obj.file, 'w');
                fprintf(obj.fileID, ''); %print header TODO
            end
        end
        
        function stepFinished(obj, time, step, xv, lastxv, wwdot, deltat)
                obj.time(obj.sc) = time;
                obj.step(obj.sc) = step;
                obj.xv(obj.sc,1:6) = xv;
                obj.lastxv(obj.sc,1:6) = lastxv;
                obj.wwdot(obj.sc,1:6) = wwdot;
                obj.deltat(obj.sc) = deltat;
                
                if ~obj.bUseLogfile
                    obj.sc=obj.sc+1;
                end
        end
        
        function ioSetNextPulse(obj, datastruct)
            obj.pointing_error(obj.sc) = datastruct.poe;
            obj.station_height(obj.sc) = datastruct.h0;
            obj.zenith_angle(obj.sc) = datastruct.za;
            obj.target_distance(obj.sc) = datastruct.tdist;
            obj.max_intensity_offset(obj.sc) = datastruct.turberror;
            obj.max_intensity(obj.sc) = datastruct.maxint;
            obj.rytov(obj.sc) = datastruct.rytov;
            obj.fried(obj.sc) = datastruct.fried;
            obj.number_screens(obj.sc) = datastruct.nscr;
        end
        
        function ioGetImpulse(obj, mom, amom, datastruct)
            obj.momentum(obj.sc,1:3) = mom;
            obj.amomentum(obj.sc,1:3) = amom;
            obj.mean_momentum(obj.sc) = datastruct.mimp;
            obj.mean_angmomentum(obj.sc) = datastruct.mamom;
            obj.mean_area(obj.sc) = datastruct.marea;
            obj.mean_cm(obj.sc) = datastruct.mcm;
            obj.mean_aoi(obj.sc) = datastruct.maoi;
            obj.mean_intensity(obj.sc) = datastruct.mint;
            obj.number_hits(obj.sc) = datastruct.nhits;
            obj.fom(obj.sc) = datastruct.fom;
            obj.photon_pressure(obj.sc) = 0; %TODO
        end
        
        function fh = plotMovement(obj)
            if ~obj.bUseLogfile
                fh = figure;
                subplot(2,2,1);
                plot(obj.time, obj.xv(:,1),obj.time, obj.xv(:,2),obj.time, obj.xv(:,3));
                title('Position');
                legend({'x','y','z'});
                xlabel('time [s]');
                ylabel('distance [m]');
                grid;

                subplot(2,2,2);
                plot(obj.time, obj.xv(:,4),obj.time, obj.xv(:,5),obj.time, obj.xv(:,6));
                title('Velocity');
                legend({'x','y','z'});
                xlabel('time [s]');
                ylabel('velocity [m/s]');
                grid;

                % subplot(2,3,3);
                % plot(obj.time, lmom(:,1),obj.time, lmom(:,2),obj.time, lmom(:,3));
                % title('Linear Momentum');
                % legend({'x','y','z'});
                % xlabel('obj.time [s]');
                % ylabel('Momentum [Nm]');

                subplot(2,2,3);
                plot(obj.time, obj.wwdot(:,3),obj.time, obj.wwdot(:,2),obj.time, obj.wwdot(:,1));
                title('Euler Angles');
                xlabel('time [s]');
                ylabel('angle [rad]');
                grid;

                subplot(2,2,4);
                plot(obj.time, obj.wwdot(:,4),obj.time, obj.wwdot(:,5),obj.time, obj.wwdot(:,6));
                title('Euler Angular Rate');
                legend({'x','y','z'});
                xlabel('time [s]');
                ylabel('anglular rate [rad/s]');
                grid;
                % subplot(2,3,6);
                % plot(obj.time, hmom(:,1),obj.time, hmom(:,2),obj.time, hmom(:,3));
                % title('Angular Momentum');
                % legend({'x','y','z'});
                % xlabel('obj.time [s]');
                % ylabel('Momentum [Nm]');

                % saveas(gcf, 'figures\test_getimpulse.png');
                savefig('figures\test_getimpulseplusbeamprop_xv.fig');
            else
                fprintf('Data written to logfile. No access at this point\n');
            end
        end
    end
end

