    module nrlmsise00_constants

    public
    ! Kind helpers
    integer , parameter :: dp = kind(1.0D0)
    integer , parameter :: sp = kind(1.0)
    integer , parameter :: r8 = kind(1.0D0)
    integer , parameter :: r4 = kind(1.0)
    integer , parameter :: shortint = 2
    integer , parameter :: normint = 4
    integer , parameter :: longint = 8

    ! Numbers 1-10
    real(dp) , parameter :: d_0 = 0.0_r8
    real(dp) , parameter :: d_1 = 1.0_r8
    real(dp) , parameter :: d_2 = 2.0_r8
    real(dp) , parameter :: d_3 = 3.0_r8
    real(dp) , parameter :: d_4 = 4.0_r8
    real(dp) , parameter :: d_5 = 5.0_r8
    real(dp) , parameter :: d_6 = 6.0_r8
    real(dp) , parameter :: d_7 = 7.0_r8
    real(dp) , parameter :: d_8 = 8.0_r8
    real(dp) , parameter :: d_9 = 9.0_r8
    real(dp) , parameter :: d_10 = 10.0_r8
    real(dp) , parameter :: d_16 = 16.0_r8
    real(dp) , parameter :: d_32 = 32.0_r8

    ! Simple Fractions
    real(dp) , parameter :: d_half = d_1/d_2

    real(dp) , parameter :: d_1q2 = d_1/d_2
    real(dp) , parameter :: d_3q2 = d_3/d_2
    real(dp) , parameter :: d_5q2 = d_5/d_2
    real(dp) , parameter :: d_7q2 = d_7/d_2
    real(dp) , parameter :: d_9q2 = d_9/d_2
    real(dp) , parameter :: d_1q3 = d_1/d_3
    real(dp) , parameter :: d_2q3 = d_2/d_3
    real(dp) , parameter :: d_4q3 = d_4/d_3
    real(dp) , parameter :: d_5q3 = d_5/d_3
    real(dp) , parameter :: d_1q4 = d_1/d_4
    real(dp) , parameter :: d_3q4 = d_3/d_4
    real(dp) , parameter :: d_5q4 = d_5/d_4
    real(dp) , parameter :: d_7q4 = d_7/d_4
    real(dp) , parameter :: d_1q5 = d_1/d_5
    real(dp) , parameter :: d_2q5 = d_2/d_5
    real(dp) , parameter :: d_3q5 = d_3/d_5
    real(dp) , parameter :: d_4q5 = d_4/d_5
    real(dp) , parameter :: d_6q5 = d_6/d_5
    real(dp) , parameter :: d_7q5 = d_7/d_5

    ! 100 and 1000 and their inverse
    real(dp) , parameter :: d_100 = 100.0_r8
    real(dp) , parameter :: d_1000 = 1000.0_r8
    real(dp) , parameter :: d_r100 = 0.01_r8
    real(dp) , parameter :: d_r1000 = 0.001_r8

    ! Power of 10
    real(dp) , parameter :: d_10E1  = 10.0_r8
    real(dp) , parameter :: d_10E2  = 100.0_r8
    real(dp) , parameter :: d_10E3  = 1000.0_r8
    real(dp) , parameter :: d_10E4  = 10000.0_r8
    real(dp) , parameter :: d_10E5  = 100000.0_r8
    real(dp) , parameter :: d_10E6  = 1000000.0_r8
    real(dp) , parameter :: d_10E7  = 10000000.0_r8
    real(dp) , parameter :: d_10E8  = 100000000.0_r8
    real(dp) , parameter :: d_10E9  = 1000000000.0_r8
    real(dp) , parameter :: d_10E10 = 10000000000.0_r8
    real(dp) , parameter :: d_10E11 = 100000000000.0_r8
    real(dp) , parameter :: d_10E12 = 1000000000000.0_r8
    real(dp) , parameter :: d_10E13 = 10000000000000.0_r8
    real(dp) , parameter :: d_10EM1  = 0.1_r8
    real(dp) , parameter :: d_10EM2  = 0.01_r8
    real(dp) , parameter :: d_10EM3  = 0.001_r8
    real(dp) , parameter :: d_10EM4  = 0.0001_r8
    real(dp) , parameter :: d_10EM5  = 0.00001_r8
    real(dp) , parameter :: d_10EM6  = 0.000001_r8
    real(dp) , parameter :: d_10EM7  = 0.0000001_r8
    real(dp) , parameter :: d_10EM8  = 0.00000001_r8
    real(dp) , parameter :: d_10EM9  = 0.000000001_r8
    real(dp) , parameter :: d_10EM10 = 0.0000000001_r8
    real(dp) , parameter :: d_10EM11 = 0.00000000001_r8
    real(dp) , parameter :: d_10EM12 = 0.000000000001_r8
    real(dp) , parameter :: d_10EM13 = 0.0000000000001_r8

    ! Time related constants
    real(dp) , parameter :: d_seconds_per_hour = 3600.0_r8
    real(sp) , parameter :: s_seconds_per_hour = 3600.0
    integer  , parameter :: i_seconds_per_hour = 3600

    real(dp) , parameter :: d_seconds_per_minute = 60.0_r8
    real(sp) , parameter :: s_seconds_per_minute = 60.0
    integer  , parameter :: i_seconds_per_minute = 60

    real(dp) , parameter :: d_minutes_per_hour = 60.0_r8
    real(sp) , parameter :: s_minutes_per_hour = 60.0
    integer  , parameter :: i_minutes_per_hour = 60

    real(dp) , parameter :: d_hours_per_day = 24.0_r8
    real(sp) , parameter :: s_hours_per_day = 24.0
    integer  , parameter :: i_hours_per_day = 24

    real(dp) , parameter :: d_minutes_per_day = 1440.0_r8
    real(sp) , parameter :: s_minutes_per_day = 1440.0
    integer  , parameter :: i_minutes_per_day = 1440

    real(dp) , parameter :: d_seconds_per_day = 86400.0_r8
    real(sp) , parameter :: s_seconds_per_day = 86400.0
    integer  , parameter :: i_seconds_per_day = 86400

    real(dp) , parameter :: d_months_per_year = 12.0_r8
    real(sp) , parameter :: s_months_per_year = 12.0
    integer  , parameter :: i_months_per_year = 12

    ! Degrees helpers
    real(dp) , parameter :: d_000 = 0.0D0
    real(dp) , parameter :: d_045 = 45.0D0
    real(dp) , parameter :: d_090 = 90.0D0
    real(dp) , parameter :: d_180 = 180.0D0
    real(dp) , parameter :: d_360 = 360.0D0

    ! Trigonometric
    real(dp) , parameter :: mathpi  = 3.1415926535897932384626433832795029_r8
    real(dp) , parameter :: invpi   = d_1/mathpi
    real(dp) , parameter :: twopi   = d_2*mathpi
    real(dp) , parameter :: halfpi  = d_1q2*mathpi
    real(dp) , parameter :: deg2rad = mathpi/d_180
    real(dp) , parameter :: rad2deg = d_180/mathpi

    integer , parameter :: maxpath = 256
    ! Basic phisical constants of matter. The 10^-23 and 10^23 factors
    ! have been removed to ease calculations.

    ! Boltzman Constant - CODATA 2010
    real(8) , parameter :: boltzk = 1.380648813_r8
    ! Avogadro Constant - International Avogadro cohordination 2011
    real(8) , parameter :: navgdr = 6.0221407818_r8

    ! Gas constant in J/K/mol
    real(dp) , parameter :: rgasmol = boltzk*navgdr

    ! Molecular weight of water
    real(dp) , parameter :: watmolwgt = 18.01528_r8 ! g/mol
    ! Mean dry air molecular weight
    real(dp) , parameter :: airmolwgt = 28.96443_r8 ! g/mol
    ! Ratio of mean molecular weight of water to that of dry air
    real(dp) , parameter :: wgtfrac = watmolwgt/airmolwgt

    ! Gas constant for dry air in J/K/kg
    real(dp) , parameter :: rgas = (rgasmol/airmolwgt)*d_1000
    ! Gas constant for water in J/K/kg
    real(dp) , parameter :: rwat = (rgasmol/watmolwgt)*d_1000

    ! 0 C in Kelvin
    real(dp) , parameter :: tzero = 273.16_r8

    ! Standard Gravity (m/sec**2) 3rd CGPM
    real(dp) , parameter :: egrav = 9.80665_r8

    ! Earth radius in meters
    real(8) , parameter :: earthrad = 6.371229D+06
    real(8) , parameter :: erkm = earthrad/d_1000
    ! Angular velocity of rotation of Earth
    real(8) , parameter :: eomeg = 7.2921159D-05
    real(8) , parameter :: eomeg2 = d_2*eomeg

    ! Hydrostatic coefficient
    real(8) , parameter :: gmr = egrav*airmolwgt/rgasmol

    ! Specific heat at constant pressure for dry air J/kg/K
    real(dp) , parameter :: cpd = 1005.46_r8
    ! Specific heat at constant pressure for moist air J/kg/K
    real(dp) , parameter :: cpv = 1869.46_r8
    ! Specific heat of water at 15 Celsius J/kg/K
    real(dp) , parameter :: cpw = 4186.95_r8
    ! Specific heat of water at 0 Celsius J/kg/K
    real(dp) , parameter :: cpw0 = 4218.0_r8

    ! Derived
    real(dp) , parameter :: rgovrw = rgas/rwat
    real(dp) , parameter :: rwovrg = rwat/rgas
    real(dp) , parameter :: rgovcp = rgas/cpd
    real(dp) , parameter :: rgovg  = rgas/egrav
    real(dp) , parameter :: govrg  = egrav/rgas

    real(dp) , parameter :: regrav = d_1/egrav
    real(dp) , parameter :: rrgas = d_1/rgas
    real(dp) , parameter :: rcpd = d_1/cpd

    end module nrlmsise00_constants