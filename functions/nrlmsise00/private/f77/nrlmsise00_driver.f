      PROGRAM MAIN
      USE m_nrlmsise00
C      TEST DRIVER FOR GTD7 (ATMOSPHERIC MODEL)
      REAL(8) ::  D(9,16),T(2,16),SW(23),APH(7)
      INTEGER(4) ::  IDAY(15)
      REAL(8) ::  UT(15),ALT(15),XLAT(15),XLONG(15),XLST(15),
     & F107A(15),F107(15),AP(15)    
C      COMMON/GTS3C/DL(16)
C      COMMON/DATIME/ISDATE(3),ISTIME(2),NAME(2)
      DATA IDAY/172,81,13*172/
      DATA UT/29000.,29000.,75000.,12*29000./
      DATA ALT/400.,400.,1000.,100.,6*400.,0,10.,30.,50.,70./
      DATA XLAT/4*60.,0.,10*60./
      DATA XLONG/5*-70.,0.,9*-70./
      DATA XLST/6*16.,4.,8*16./
      DATA F107A/7*150.,70.,7*150./
      DATA F107/8*150.,180.,6*150./
      DATA AP/9*4.,40.,5*4./
      DATA APH/7*100./,SW/8*1.,-1.,14*1./
C
      XLST = UT/3600.0 + XLONG/15.0  
      CALL METERS(.TRUE.) 
      CALL TSELEC(SW)
      WRITE(6,'("-----------------------------------------------")')
      DO I=1,10
         CALL GTD7(IDAY(I),UT(I),ALT(I),XLAT(I),XLONG(I),XLST(I),
     &             F107A(I),F107(I),APH,48,D(1,I),T(1,I))
         WRITE(6,100) I, (D(J,I),J=1,9)
         WRITE(6,101) T(1,I),T(2,I)
      ENDDO
      
  100 FORMAT(1X,I2,6(1X,1PE16.10)/4X,3(1X,1PE16.10))
  101 FORMAT(4X,2F8.3/)
  200 FORMAT(//' DAY  ',5I12)
  201 FORMAT(' UT   ',5F12.0)
  202 FORMAT(' ALT  ',5F12.0)
  203 FORMAT(' LAT  ',5F12.0)
  204 FORMAT(' LONG ',5F12.0)
  205 FORMAT(' LST  ',5F12.0)
  206 FORMAT(' F107A',5F12.0)
  207 FORMAT(' F107 ',5F12.0)
  208 FORMAT(' AP   ',5F12.0)
  210 FORMAT(/' TINF ',5F12.2)
  211 FORMAT(' TG   ',5F12.2)
  212 FORMAT(' HE   ',1P5E12.3)
  213 FORMAT(' O    ',1P5E12.3)
  214 FORMAT(' N2   ',1P5E12.3)
  215 FORMAT(' O2   ',1P5E12.3)
  216 FORMAT(' AR   ',1P5E12.3)
  217 FORMAT(' H    ',1P5E12.3)
  219 FORMAT(' N    ',1P5E12.3)
  220 FORMAT(' ANM O',1P5E12.3)
  218 FORMAT(' RHO  ',1P5E12.3)
  300 FORMAT(1X,2A4,2X,3A4,2X,2A4)
      STOP
      END PROGRAM MAIN