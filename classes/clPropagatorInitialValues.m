classdef clPropagatorInitialValues
   properties
       % starting epoch
       year = 2002
       mon = 4
       day = 24
       hour = 21
       min = 55
       sec = 28
       
       % propagation time
       Mjd_UTC % staring time  [MJD]
       mjd 
       proptime = 250 % duration [s]
       step % step length
       
       dof = 6
            
       a_grav = 1 % harmonic terms of gravity
       g_grav = 0
       DSPOSE = 1 % using DSPOSE - check if it can be deleted as it seems its not used anymore
       n = 40
       m = 40
       
       a_lase = 0 % laser engagements
       g_lase = 0
       a_drag = 0 % drag
       g_drag = 0
       draggradient = 0 %consider drag gradient
       a_srad = 0 % solar radiation pressure
       g_srad = 0
       a_erad = 0 % reflected and emmited earth radiation pressure
       g_erad = 0
       g_mag  = 0 % magnetic torque
       a_sun = 0 % gravity of sun
       a_moon = 0 % gravity of moon
       a_planets = 0 % gravity of other planets
       
       a_solidEarthTides = 0 % solid earth tides
       a_oceanTides = 0 % ocean tides
       a_relativity = 0 % relativity effects
       
       albedo_model %albedo_model
       
       b_debugging = 0
       b_logging = 1
       n_log_only_every_xth_step = 1
   end
    
end