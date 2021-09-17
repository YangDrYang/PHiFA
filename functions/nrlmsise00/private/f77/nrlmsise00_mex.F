#include "fintrf.h"

C ------------------------------------------------------------------------
C     nrlmsise00.f - MATLAB gateway to nrlmsise00
C     which is used to calculate the Jacchia-Bowman 2008 Model Atmosphere
C     04 Feb 2016
C     gateway written by Changyong He <changyong.he@rmit.edu.au>
C ------------------------------------------------------------------------
      subroutine mexfunction(nlhs,plhs,nrhs,prhs)
C     Declarations
      implicit none
C     Required declarations for MATLAB
      mwPointer   ::  plhs(*),prhs(*)
      integer     ::  nlhs,nrhs
      
      mwPointer   ::  mxGetPr
      mwPointer   ::  mxGetM, mxGetN
      mwPointer   ::  mxCreateNumericMatrix
      integer*4   ::  mxClassIDFromClassName
      integer*4   ::  mxIsInt32, mxIsDouble
C     input variables
      integer*4   ::  m, n
      integer*4   ::  classid
      integer     ::  err_mesg
      real*8, dimension(23) :: FLAG
      integer*4, dimension(:,:), allocatable ::  IYD
      real*8, dimension(:,:), allocatable ::  SEC, ALT, GLAT,GLONG,
     +    F107A,F107,AP
C     output variables
      REAL*8, dimension(:,:), allocatable :: D, T
C     Local variables
C      integer*4,parameter :: MASS = 48      
C Associate classid with single precision (FORTRAN real)
      classid = mxClassIDFromClassName('double')
C      classid_int = mxClassIDFromClassName('int32')
C Check the inputs
      if (nrhs /= 9) then
          call mexErrMsgIdAndTxt ('MATLAB:nrlmsise00:invalidInputs',
     +    '9 inputs required.')
      end if
C     Validate inputs
C     Check to see both inputs are numeric.
      if (mxIsInt32(prhs(1)) .ne. 1) then
          call mexErrMsgIdAndTxt ('MATLAB:nrlmsise00:NonInteger',
     +    'Input #1 is not a int8 variable or a array.')
      end if 
      if (mxIsDouble(prhs(2)) .ne. 1) then
          call mexErrMsgIdAndTxt ('MATLAB:nrlmsise00:NonDouble',
     +    'Input #2 is not a double variable or a array.')
      end if
      if (mxIsDouble(prhs(3)) .ne. 1) then
          call mexErrMsgIdAndTxt ('MATLAB:nrlmsise00:NonDouble',
     +    'Input #3 is not a double variable or a array.')
      end if
      if (mxIsDouble(prhs(4)) .ne. 1) then
          call mexErrMsgIdAndTxt ('MATLAB:nrlmsise00:NonDouble',
     +    'Input #4 is not a double variable or a array.')
      end if

      if (mxIsDouble(prhs(5)) .ne. 1) then
          call mexErrMsgIdAndTxt ('MATLAB:nrlmsise00:NonDouble',
     +    'Input #5 is not a double array.')
      end if
      if (mxIsDouble(prhs(6)) .ne. 1) then
          call mexErrMsgIdAndTxt ('MATLAB:nrlmsise00:NonDouble',
     +    'Input #6 is not a double array.')
      end if
      if (mxIsDouble(prhs(7)) .ne. 1) then
          call mexErrMsgIdAndTxt ('MATLAB:nrlmsise00:NonDouble',
     +    'Input #7 is not a double array.')
      end if
      if (mxIsDouble(prhs(8)) .ne. 1) then
          call mexErrMsgIdAndTxt ('MATLAB:nrlmsise00:NonDouble',
     +    'Input #8 is not a double array.')
      end if
      if (mxIsDouble(prhs(9)) .ne. 1) then
          call mexErrMsgIdAndTxt ('MATLAB:nrlmsise00:NonDouble',
     +    'Input #9 is not a double array.')
      end if
   

C     Check the size of inputs
      m = mxGetM(prhs(1));
      n = mxGetN(prhs(1));
      if (n .NE. 1) then !iday
          call mexErrMsgIdAndTxt ('MATLAB:nrlmsise00:Invalid',
     +    'Input # 1 is not a variable or a vector.')
      end if
      
      if (mxGetM(prhs(2))/=m.or.mxGetN(prhs(2))/=1.or.!ut
     +mxGetM(prhs(3))/=m.or.mxGetN(prhs(3))/=1.or.    !alt
     +mxGetM(prhs(4))/=m.or.mxGetN(prhs(4))/=1.or.    !xlat
     +mxGetM(prhs(5))/=m.or.mxGetN(prhs(5))/=1.or.    !xlong
     +mxGetM(prhs(6))/=m.or.mxGetN(prhs(6))/=1.or.    !f107a
     +mxGetM(prhs(7))/=m.or.mxGetN(prhs(7))/=1.or.    !f107
     +mxGetM(prhs(8))/=m.or.mxGetN(prhs(8))/=7.or.    !ap
     +mxGetM(prhs(9))/=1.or.mxGetN(prhs(9))/=23) then !flag
          call mexErrMsgIdAndTxt ('MATLAB:nrlmsise00:InvalidSize',
     +    'Check the size of inputs.')
      end if
C     Get variables passed from MATLAB, properly typed
C     allocate the memory to inputs
C     (IYD,SEC,ALT
      allocate(IYD(m,1), stat = err_mesg)  
      if(0 .ne. err_mesg) then  
          call mexErrMsgIdAndTxt(
     +    'MATLAB:nrlmsise00:fail to allocate momery to Input #1',
     +    err_mesg)
      end if
      call mxCopyPtrToInteger4(mxGetPr(prhs(1)),IYD,m*1)

      allocate(SEC(m,1), stat = err_mesg)  
      if(0 .ne. err_mesg) then  
          call mexErrMsgIdAndTxt(
     +    'MATLAB:nrlmsise00:fail to allocate momery to Input #2',
     +    err_mesg)
      end if
      call mxCopyPtrToReal8(mxGetPr(prhs(2)),SEC,m*1)

      allocate(ALT(m,1), stat = err_mesg)  
      if(0 .ne. err_mesg) then  
          call mexErrMsgIdAndTxt(
     +    'MATLAB:nrlmsise00:fail to allocate momery to Input #3',
     +    err_mesg)
      end if
      call mxCopyPtrToReal8(mxGetPr(prhs(3)),ALT,m*1)
C     GLAT,GLONG
      allocate(GLAT(m,1), stat = err_mesg)  
      if(0 .ne. err_mesg) then  
          call mexErrMsgIdAndTxt(
     +    'MATLAB:nrlmsise00:fail to allocate momery to Input #4',
     +    err_mesg)
      end if
      call mxCopyPtrToReal8(mxGetPr(prhs(4)),GLAT,m*1)
      
      allocate(GLONG(m,1), stat = err_mesg)  
      if(0 .ne. err_mesg) then  
          call mexErrMsgIdAndTxt(
     +    'MATLAB:nrlmsise00:fail to allocate momery to Input #5',
     +    err_mesg)
      end if
      call mxCopyPtrToReal8(mxGetPr(prhs(5)),GLONG,m*1)

      
C     F107A,F107,AP,MASS,
      allocate(F107A(m,1), stat = err_mesg)  
      if(0 .ne. err_mesg) then  
          call mexErrMsgIdAndTxt(
     +    'MATLAB:nrlmsise00:fail to allocate momery to Input #6',
     +    err_mesg)
      end if
      call mxCopyPtrToReal8(mxGetPr(prhs(6)),F107A,m*1)

      allocate(F107(m,1), stat = err_mesg)  
      if(0 .ne. err_mesg) then  
          call mexErrMsgIdAndTxt(
     +    'MATLAB:nrlmsise00:fail to allocate momery to Input #7',
     +    err_mesg)
      end if
      call mxCopyPtrToReal8(mxGetPr(prhs(7)),F107,m*1)
      
      allocate(AP(m,7), stat = err_mesg)  
      if(0 .ne. err_mesg) then  
          call mexErrMsgIdAndTxt(
     +    'MATLAB:nrlmsise00:fail to allocate momery to Input #8',
     +    err_mesg)
      end if
      call mxCopyPtrToReal8(mxGetPr(prhs(8)),AP,m*7)
      
      call mxCopyPtrToReal8(mxGetPr(prhs(9)),FLAG,23)
C     Check outputs and Create MATLAB return matrices 
      if (nlhs == 2) then
          plhs(1) = mxCreateNumericMatrix(m,9,classid,0)
          plhs(2) = mxCreateNumericMatrix(m,2,classid,0)
      elseif (nlhs == 1) then
          plhs(1) = mxCreateNumericMatrix(m,9,classid,0)
      else 
          call mexErrMsgIdAndTxt('MATLAB:nrlmsise00:invalidOutputs',
     +    'Two outputs required.')
      end if
C     allocate the memory to outputs (m x 9)
      allocate(D(m,9), stat = err_mesg)  
      if(0 .ne. err_mesg) then  
          call mexErrMsgIdAndTxt(
     +    'MATLAB:nrlmsise00:fail to allocate momery to output #1',
     +    err_mesg)
      end if

      allocate(T(m,2), stat = err_mesg)  
      if(0 .ne. err_mesg) then  
          call mexErrMsgIdAndTxt(
     +    'MATLAB:nrlmsise00:fail to allocate momery to output #2',
     +    err_mesg)
      end if

c allocate the memory to the output parameters
c Do the deed: call MSIS
      call nrlmsise00_array( m, IYD,SEC,ALT,GLAT,
     +                       GLONG,F107A,F107,AP,FLAG,
     +                       D, T)

c Put the returned MSIS data into the MATLAB return matrices we created       
C     D, T
      if (nlhs == 2) then
          call mxCopyReal8ToPtr(D,mxGetPr(plhs(1)),m*9)
          call mxCopyReal8ToPtr(T,mxGetPr(plhs(2)),m*2)
      elseif (nlhs == 1) then
          call mxCopyReal8ToPtr(D,mxGetPr(plhs(1)),m*9)
      end if 

c Done
      deallocate(IYD,SEC)
      deallocate(ALT,GLAT,GLONG)
      deallocate(F107A,F107,AP)

      deallocate(D,T)
C Done
c      return
      end
      
C***********************************************************************
      subroutine nrlmsise00_array( m, IYD,SEC,ALT,GLAT,GLONG,F107A,
     +    F107,AP,FLAG,D,T)
C
      USE m_nrlmsise00
      IMPLICIT NONE
C
      INTEGER*4 :: I, m, IYD(m,1)
      INTEGER*4,PARAMETER :: MASS = 48
      REAL*8 :: FLAG(23)
      REAL*8 :: SEC(m,1),ALT(m,1),GLAT(m,1),GLONG(m,1)
      REAL*8 :: STL,F107A(m,1),F107(m,1),AP(m,7)
      REAL*8 :: D(m,9),T(m,2)
c Return in kg/m^-3 units, not kg/cm^-3
      call METERS(.true.)
C FLAGS  :a numerical array of 25 values for setting particular
C         variations in calculation the output.
      call TSELEC(FLAG)
      DO I=1,m,1
          STL = DMOD(SEC(I,1)/3600.0D0 + GLONG(I,1)/15.0D0,24.0D0)
C     COMPUTE DENSITY KG/M3 RHO
          CALL GTD7(IYD(I,1),SEC(I,1),ALT(I,1),GLAT(I,1),GLONG(I,1),
     +        STL,F107A(I,1),F107(I,1),AP(I,1:7),MASS,
     +        D(I,1:9),T(I,1:2))
c          write(6,*) D(I,6)
      END DO
c
      RETURN
      END
