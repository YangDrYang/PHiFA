function solarInt = getSolarIntensity(au, mjd)
% function to get solar intensity at distance to sun and given time. atm
% only const implemented which is sufficient for the model, but it can be
% easily expanded using function in SC Toolbox

p = clPropagator.instance();
if nargin == 0
    solarInt = p.const.I_Sol;
elseif nargin == 1
    solarInt = SolarFlx(au);
elseif nargin == 2
%     [aP, f, fHat, fHat400] = SolarFluxPrediction( jD, timing ) % uses SCT
%     function
    solarInt = SolarFlx(au); % TODO see comment above
end

end