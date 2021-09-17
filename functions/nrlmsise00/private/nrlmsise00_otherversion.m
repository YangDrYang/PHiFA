%{
---------------------------------------------------------------------------------
C     NRLMSISE-00
C     -----------
C        Neutral Atmosphere Empirical Model from the surface to lower
C        exosphere
C
C        NEW FEATURES:
C          *Extensive satellite drag database used in model generation
C          *Revised O2 (and O) in lower thermosphere
C          *Additional nonlinear solar activity term
C          *"ANOMALOUS OXYGEN" NUMBER DENSITY, OUTPUT D(9)
C           At high altitudes (> 500 km), hot atomic oxygen or ionized
C           oxygen can become appreciable for some ranges of subroutine
C           inputs, thereby affecting drag on satellites and debris. We
C           group these species under the term "anomalous oxygen," since
C           their individual variations are not presently separable with
C           the drag data used to define this model component.
C
C        SUBROUTINES FOR SPECIAL OUTPUTS:
C
C        HIGH ALTITUDE DRAG: EFFECTIVE TOTAL MASS DENSITY
C        (SUBROUTINE GTD7D, OUTPUT D(6))
C           For atmospheric drag calculations at altitudes above 500 km,
C           call SUBROUTINE GTD7D to compute the "effective total mass
C           density" by including contributions from "anomalous oxygen."
C           See "NOTES ON OUTPUT VARIABLES" below on D(6).
C
C        PRESSURE GRID (SUBROUTINE GHP7)
C          See subroutine GHP7 to specify outputs at a pressure level
C          rather than at an altitude.
C
C     INPUT VARIABLES:
C        IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
C              (Year ignored in current model)
C        SEC - UT(SEC)
C        ALT - ALTITUDE(km)
C        GLAT - GEODETIC LATITUDE(degree)
C        GLONG - GEODETIC LONGITUDE(degree)
C        STL - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
C        F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
C        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
C        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
C           - ARRAY CONTAINING:
C             (1) DAILY AP
C             (2) 3 HR AP INDEX FOR CURRENT TIME
C             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
C             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
C             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
C             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
C                    TO CURRENT TIME
C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
C                    TO CURRENT TIME
C        MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
C                 CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
C                 MASS 17 IS Anomalous O ONLY.)
%   FLAGS  :a numerical array of 25 values for setting particular
%          variations in calculation the output.  Setting a value to 0.0
%          removes that value's effect on the output.  Setting a value to
%          1.0 applies the main and the cross term effects of that value
%          on the output.  Setting a value to 2.0 applies only the cross
%          term effect of that value on the output.  Additionally setting
%          FLAGS(9) = -1 uses the entire matrix APH rather than just
%          APH(:,1). The variations contained in FLAGS are ordered as
%          follows:
%           FLAGS(1)  :F10.7 effect on mean
%           FLAGS(2)  :Time independent
%           FLAGS(3)  :Symmetrical annual
%           FLAGS(4)  :Symmetrical semi-annual
%           FLAGS(5)  :Asymmetrical annual
%           FLAGS(6)  :Asymmetrical semi-annual
%           FLAGS(7)  :Diurnal
%           FLAGS(8)  :Semi-diurnal
%           FLAGS(9)  :Daily AP
%           FLAGS(10) :All UT, longitudinal effects
%           FLAGS(11) :Longitudinal
%           FLAGS(12) :UT and mixed UT, longitudinal
%           FLAGS(13) :Mixed AP, UT, longitudinal
%           FLAGS(14) :Ter-diurnal
%           FLAGS(15) :Departures from diffusive equilibrium
%           FLAGS(16) :All exospheric temperature variations
%           FLAGS(17) :All variations from 120,000 meter temperature (TLB)
%           FLAGS(18) :All lower thermosphere (TN1) temperature variations
%           FLAGS(19) :All 120,000 meter gradient (S) variations
%           FLAGS(20) :All upper stratosphere (TN2) temperature variations
%           FLAGS(21) :All variations from 120,000 meter values (ZLB)
%           FLAGS(22) :All lower mesosphere temperature (TN3) variations
%           FLAGS(23) :Turbopause scale height variations
%          The default values are 1.0 for all FLAGS.
C
C     NOTES ON INPUT VARIABLES:
C        UT, Local Time, and Longitude are used independently in the
C        model and are not of equal importance for every situation.
C        For the most physically realistic calculation these three
C        variables should be consistent (STL=SEC/3600+GLONG/15).
C        The Equation of Time departures from the above formula
C        for apparent local time can be included if available but
C        are of minor importance.
c
C        F107 and F107A values used to generate the model correspond
C        to the 10.7 cm radio flux at the actual distance of the Earth
C        from the Sun rather than the radio flux at 1 AU. The following
C        site provides both classes of values:
C        ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/
C
C        F107, F107A, and AP effects are neither large nor well
C        established below 80 km and these parameters should be set to
C        150., 150., and 4. respectively.
C
C     OUTPUT VARIABLES:
C        D(1) - HE NUMBER DENSITY(m-3)
C        D(2) - O NUMBER DENSITY(m-3)
C        D(3) - N2 NUMBER DENSITY(m-3)
C        D(4) - O2 NUMBER DENSITY(m-3)
C        D(5) - AR NUMBER DENSITY(m-3)
C        D(6) - TOTAL MASS DENSITY(kg/m3)
C        D(7) - H NUMBER DENSITY(m-3)
C        D(8) - N NUMBER DENSITY(m-3)
C        D(9) - Anomalous oxygen NUMBER DENSITY(m-3)
C        T(1) - EXOSPHERIC TEMPERATURE
C        T(2) - TEMPERATURE AT ALT
C
C        O, H, and N are set to zero below 72.5 km
C
C        T(1), Exospheric temperature, is set to global average for
C        altitudes below 120 km. The 120 km gradient is left at global
C        average value for altitudes below 72 km.
C
C        D(6), TOTAL MASS DENSITY, is NOT the same for subroutines GTD7
C        and GTD7D
C
C          SUBROUTINE GTD7 -- D(6) is the sum of the mass densities of the
C          species labeled by indices 1-5 and 7-8 in output variable D.
C          This includes He, O, N2, O2, Ar, H, and N but does NOT include
C          anomalous oxygen (species index 9).
C
C          SUBROUTINE GTD7D -- D(6) is the "effective total mass density
C          for drag" and is the sum of the mass densities of all species
C          in this model, INCLUDING anomalous oxygen.
C
%   Matlab code Copy right: Changyong HE email:changyong.he@rmit.edu.au
---------------------------------------------------------------------------------
%}
function [d, t] = nrlmsise00(iyear,idoy,ut,alt,xlat,xlong,...
  f107a, f107, ap,flag, Oplus_flag)
%% Check for latest version and all required files.

if nargin < 9 || nargin>10
  error('nrlmsise00:the number of inputs is invalid!\n');
elseif nargin == 9
  Oplus_flag = false;
end
if strcmpi('Oxygen',Oplus_flag)
  Oplus_flag = true;
end

% mass = 48*ones(size(idoy),'int32');%default value
%% Call the nrlmsise00 function

% lst = sec(1)/3600 + lon(1)/15;
lst = ut/3600 + xlong/15;
iyd = mod(iyear,100)*1000 + fix(idoy);
[d,t] = nrlmsise00_mex(int32(iyd),ut,alt,xlat,xlong,lst,f107a,f107,ap,(flag));
if(Oplus_flag)
  d(:,6) = 1.66e-27*(16.0*d(:,9)) + d(:,6);
end
% d = double(d);
% t = double(t);


end
