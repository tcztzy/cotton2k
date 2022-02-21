! File VersionID:
!   $Id: sptabulated.f90 299 2016-08-04 17:30:36Z kroes006 $
! ----------------------------------------------------------------------
      subroutine EvalTabulatedFunction(inverse,n,ind1,ind2,ind3,node,   &
     &                                  sptab,ientrytab,xe,ye,dyedxe)

      implicit none
      include 'arrays.fi'

      integer k, klo, khi, n, inverse, maxtry, ntry, ind1, ind2, ind3
      real(8) xe, ye, dyedxe
      real(8) x1, x2, f1, f2, d1, d2, h, delta, del1, del2, c2, c2t2,   &
     &        c3, c3t3, xx
      real(8) xlow, xhigh, xtry, ytry
      real(8) sptab(5,macp,matab)
      integer ientrytab(macp,0:matabentries), node    


!     'ye' is only used as 'inverse' = 1. The option that requires this is not yet operationel
!     (combination of sophys-tables (swsophy=1) and swdiscrvert=1).
! Therefore, 'ye' initial gets a dummy value for Forecheck
!      if (inverse == 1) ye = -99.d0
!      if (inverse == 0) xe = -99.d0
       ye = 0.d0

!      sptab(1,node,i):    Theta
!      sptab(2,node,i):    h
!      sptab(3,node,i):    dTheta / dh
!      sptab(4,node,i):    K
!      sptab(5,node,i):    dK / dh

!     bi-section search method for finding the array entry

      if(inverse .eq. 0)then
         if(xe.ge.-1.0d-9) then
            stop 'Invalid call of Function EvalTabulatedFunction'
         end if
         k   = int(1000*(log10(-xe)+1.d0))+1
         if(k.lt.1) k=1
         klo = ientrytab(node,k) 
         khi = klo + 1

!     bounds of interval
         x1 = sptab(ind1,node,klo)
         x2 = sptab(ind1,node,khi)

         if( (xe-x1) .lt. 0.0d0 )then
            klo = klo -1
            khi = khi -1
            x1 = sptab(ind1,node,klo)
            x2 = sptab(ind1,node,khi)
         else if( (xe-x1) .ge. (x2-x1) )then
            klo = klo +1
            khi = khi +1
            x1 = sptab(ind1,node,klo)
            x2 = sptab(ind1,node,khi)
         end if
      
         f1 = sptab(ind2,node,klo)
         f2 = sptab(ind2,node,khi)
         d1 = sptab(ind3,node,klo)
         d2 = sptab(ind3,node,khi)
         
!     coefficients to be used
         h     = x2 - x1
         delta = (f2 - f1)/h
         del1  = (d1 - delta)/h
         del2  = (d2 - delta)/h
         c2    = -(del1+del1 + del2)
         c2t2  = c2 + c2
         c3    = (del1 + del2)/h
         c3t3  = c3+c3+c3
!     distance to reference of interval
         xx = xe - x1
!     compute
         ye = f1 + xx*(d1 + xx*(c2 + xx*c3))

      else if(inverse.eq.1)then
!     bi-section search method for finding the array entry
         klo=1 
         khi=n
2000     if (khi-klo.gt.1) then
            k=(khi+klo)/2
            if(sptab(ind2,node,k).gt.ye)then
               khi=k
            else
               klo=k
            endif
            goto 2000
         endif 
         x1 = sptab(ind1,node,klo)
         x2 = sptab(ind1,node,khi)
         f1 = sptab(ind2,node,klo)
         f2 = sptab(ind2,node,khi)
         d1 = sptab(ind3,node,klo)
         d2 = sptab(ind3,node,khi)
         h     = x2 - x1
         delta = (f2 - f1)/h
         del1  = (d1 - delta)/h
         del2  = (d2 - delta)/h
         c2    = -(del1+del1 + del2)
         c2t2  = c2 + c2
         c3    = (del1 + del2)/h
         c3t3  = c3+c3+c3

!     bi-section search method for finding xx-value
         xlow   = x1
         xhigh  = x2
         maxtry = 20
         ntry   = 1
         xtry = 0.5d0*(xlow+xhigh)
         xx = xtry - x1
         ytry = f1 + xx*(d1 + xx*(c2 + xx*c3))
         do while (dabs(ye-ytry).gt.1.0d-05 .and.ntry.lt.maxtry)
            ntry = ntry + 1
            if(ytry.lt.ye)then
               xlow = xtry
            else
               xhigh = xtry
            end if
            xtry = 0.5d0*(xlow+xhigh)
            xx = xtry - x1
            ytry = f1 + xx*(d1 + xx*(c2 + xx*c3))
         end do
         xe = xx + x1
      end if

      dyedxe = d1 + xx*(c2t2 + xx*c3t3)

      return
      end

      subroutine PreProcTabulatedFunction(flag,n,x,y,dydx)
      implicit none
      include 'arrays.fi'

      integer  i,n, flag
      real(8)  y(matab), dydx(matab)
      integer  ic(2), INCFD, NWK, IERR
      real(8)  vc(2), switch, x(matab)
      real(8)  f(1,matab),d(1,matab),wk(2*matab)

!---- initialize the dydx computations

!---- second derivative at x(1) is assumed to be zero
      ic(1) = 2    
      vc(1) = 0.0d0

!---- first derivative at x(n) is assumed to be 1.0d-7 (saturation)
      if(flag.eq.1)then
         ic(2) = 1
         vc(2) = 1.0d-07
      else if(flag.eq.2)then
         ic(2) = 2
         vc(2) = 0.0
      end if

      SWITCH = 0.0d0
      INCFD = 1
      NWK = 2*(n-1)
      do i=1,n
         f(1,i) = y(i)
      end do

!---- perform the dydx computations

      call DPCHIC (IC, VC, SWITCH,n,x,f,d,INCFD,WK, NWK, IERR)

!---- store results in dydx array

      do i=1,n
         dydx(i) = d(1,i)
      end do

      return
      end

      SUBROUTINE DPCHIC (IC, VC, SWITCH, N, X, F, D, INCFD, WK, NWK,    &
     &   IERR)
!C***BEGIN PROLOGUE  DPCHIC
!C***PURPOSE  Set derivatives needed to determine a piecewise monotone
!C            piecewise cubic Hermite interpolant to given data.
!C            User control is available over boundary conditions and/or
!C            treatment of points where monotonicity switches direction.
!C***LIBRARY   SLATEC (PCHIP)
!C***CATEGORY  E1A
!C***TYPE      DOUBLE PRECISION (PCHIC-S, DPCHIC-D)
!C***KEYWORDS  CUBIC HERMITE INTERPOLATION, MONOTONE INTERPOLATION,
!C             PCHIP, PIECEWISE CUBIC INTERPOLATION,
!C             SHAPE-PRESERVING INTERPOLATION
!C***AUTHOR  Fritsch, F. N., (LLNL)
!C             Lawrence Livermore National Laboratory
!C             P.O. Box 808  (L-316)
!C             Livermore, CA  94550
!C             FTS 532-4275, (510) 422-4275
!C***DESCRIPTION
!C
!C         DPCHIC:  Piecewise Cubic Hermite Interpolation Coefficients.
!C
!C     Sets derivatives needed to determine a piecewise monotone piece-
!C     wise cubic interpolant to the data given in X and F satisfying the
!C     boundary conditions specified by IC and VC.
!C
!C     The treatment of points where monotonicity switches direction is
!C     controlled by argument SWITCH.
!C
!C     To facilitate two-dimensional applications, includes an increment
!C     between successive values of the F- and D-arrays.
!C
!C     The resulting piecewise cubic Hermite function may be evaluated
!C     by DPCHFE or DPCHFD.
!C
!C ----------------------------------------------------------------------
!C
!C  Calling sequence:
!C
!C        PARAMETER  (INCFD = ...)
!C        INTEGER  IC(2), N, NWK, IERR
!C        DOUBLE PRECISION  VC(2), SWITCH, X(N), F(INCFD,N), D(INCFD,N),
!C                          WK(NWK)
!C
!C        CALL DPCHIC (IC, VC, SWITCH, N, X, F, D, INCFD, WK, NWK, IERR)
!C
!C   Parameters:
!C
!C     IC -- (input) integer array of length 2 specifying desired
!C           boundary conditions:
!C           IC(1) = IBEG, desired condition at beginning of data.
!C           IC(2) = IEND, desired condition at end of data.
!C
!C           IBEG = 0  for the default boundary condition (the same as
!C                     used by DPCHIM).
!C           If IBEG.NE.0, then its sign indicates whether the boundary
!C                     derivative is to be adjusted, if necessary, to be
!C                     compatible with monotonicity:
!C              IBEG.GT.0  if no adjustment is to be performed.
!C              IBEG.LT.0  if the derivative is to be adjusted for
!C                     monotonicity.
!C
!C           Allowable values for the magnitude of IBEG are:
!C           IBEG = 1  if first derivative at X(1) is given in VC(1).
!C           IBEG = 2  if second derivative at X(1) is given in VC(1).
!C           IBEG = 3  to use the 3-point difference formula for D(1).
!C                     (Reverts to the default b.c. if N.LT.3 .)
!C           IBEG = 4  to use the 4-point difference formula for D(1).
!C                     (Reverts to the default b.c. if N.LT.4 .)
!C           IBEG = 5  to set D(1) so that the second derivative is con-
!C              tinuous at X(2). (Reverts to the default b.c. if N.LT.4.)
!C              This option is somewhat analogous to the "not a knot"
!C              boundary condition provided by DPCHSP.
!C
!C          NOTES (IBEG):
!C           1. An error return is taken if ABS(IBEG).GT.5 .
!C           2. Only in case  IBEG.LE.0  is it guaranteed that the
!C              interpolant will be monotonic in the first interval.
!C              If the returned value of D(1) lies between zero and
!C              3*SLOPE(1), the interpolant will be monotonic.  This
!C              is **NOT** checked if IBEG.GT.0 .
!C           3. If IBEG.LT.0 and D(1) had to be changed to achieve mono-
!C              tonicity, a warning error is returned.
!C
!C           IEND may take on the same values as IBEG, but applied to
!C           derivative at X(N).  In case IEND = 1 or 2, the value is
!C           given in VC(2).
!C
!C          NOTES (IEND):
!C           1. An error return is taken if ABS(IEND).GT.5 .
!C           2. Only in case  IEND.LE.0  is it guaranteed that the
!C              interpolant will be monotonic in the last interval.
!C              If the returned value of D(1+(N-1)*INCFD) lies between
!C              zero and 3*SLOPE(N-1), the interpolant will be monotonic.
!C              This is **NOT** checked if IEND.GT.0 .
!C           3. If IEND.LT.0 and D(1+(N-1)*INCFD) had to be changed to
!C              achieve monotonicity, a warning error is returned.
!C
!C     VC -- (input) real(8)array of length 2 specifying desired boundary
!C           values, as indicated above.
!C           VC(1) need be set only if IC(1) = 1 or 2 .
!C           VC(2) need be set only if IC(2) = 1 or 2 .
!C
!C     SWITCH -- (input) indicates desired treatment of points where
!C           direction of monotonicity switches:
!C           Set SWITCH to zero if interpolant is required to be mono-
!C           tonic in each interval, regardless of monotonicity of data.
!C             NOTES:
!C              1. This will cause D to be set to zero at all switch
!C                 points, thus forcing extrema there.
!C              2. The result of using this option with the default boun-
!C                 dary conditions will be identical to using DPCHIM, but
!C                 will generally cost more compute time.
!C                 This option is provided only to facilitate comparison
!C                 of different switch and/or boundary conditions.
!C           Set SWITCH nonzero to use a formula based on the 3-point
!C              difference formula in the vicinity of switch points.
!C           If SWITCH is positive, the interpolant on each interval
!C              containing an extremum is controlled to not deviate from
!C              the data by more than SWITCH*DFLOC, where DFLOC is the
!C              maximum of the change of F on this interval and its two
!C              immediate neighbors.
!C           If SWITCH is negative, no such control is to be imposed.
!C
!C     N -- (input) number of data points.  (Error return if N.LT.2 .)
!C
!C     X -- (input) real(8)array of independent variable values.  The
!C           elements of X must be strictly increasing:
!C                X(I-1) .LT. X(I),  I = 2(1)N.
!C           (Error return if not.)
!C
!C     F -- (input) real(8)array of dependent variable values to be
!C           interpolated.  F(1+(I-1)*INCFD) is value corresponding to
!C           X(I).
!C
!C     D -- (output) real(8)array of derivative values at the data
!C           points.  These values will determine a monotone cubic
!C           Hermite function on each subinterval on which the data
!C           are monotonic, except possibly adjacent to switches in
!C           monotonicity. The value corresponding to X(I) is stored in
!C                D(1+(I-1)*INCFD),  I=1(1)N.
!C           No other entries in D are changed.
!C
!C     INCFD -- (input) increment between successive values in F and D.
!C           This argument is provided primarily for 2-D applications.
!C           (Error return if  INCFD.LT.1 .)
!C
!C     WK -- (scratch) real(8)array of working storage.  The user may
!C           wish to know that the returned values are:
!C              WK(I)     = H(I)     = X(I+1) - X(I) ;
!C              WK(N-1+I) = SLOPE(I) = (F(1,I+1) - F(1,I)) / H(I)
!C           for  I = 1(1)N-1.
!C
!C     NWK -- (input) length of work array.
!C           (Error return if  NWK.LT.2*(N-1) .)
!C
!C     IERR -- (output) error flag.
!C           Normal return:
!C              IERR = 0  (no errors).
!C           Warning errors:
!C              IERR = 1  if IBEG.LT.0 and D(1) had to be adjusted for
!C                        monotonicity.
!C              IERR = 2  if IEND.LT.0 and D(1+(N-1)*INCFD) had to be
!C                        adjusted for monotonicity.
!C              IERR = 3  if both of the above are true.
!C           "Recoverable" errors:
!C              IERR = -1  if N.LT.2 .
!C              IERR = -2  if INCFD.LT.1 .
!C              IERR = -3  if the X-array is not strictly increasing.
!C              IERR = -4  if ABS(IBEG).GT.5 .
!C              IERR = -5  if ABS(IEND).GT.5 .
!C              IERR = -6  if both of the above are true.
!C              IERR = -7  if NWK.LT.2*(N-1) .
!C             (The D-array has not been changed in any of these cases.)
!C               NOTE:  The above errors are checked in the order listed,
!C                   and following arguments have **NOT** been validated.
!C
!C***REFERENCES  1. F. N. Fritsch, Piecewise Cubic Hermite Interpolation
!C                 Package, Report UCRL-87285, Lawrence Livermore Natio-
!C                 nal Laboratory, July 1982.  [Poster presented at the
!C                 SIAM 30th Anniversary Meeting, 19-23 July 1982.]
!C               2. F. N. Fritsch and J. Butland, A method for construc-
!C                 ting local monotone piecewise cubic interpolants, SIAM
!C                 Journal on Scientific and Statistical Computing 5, 2
!C                 (June 1984), pp. 300-304.
!C               3. F. N. Fritsch and R. E. Carlson, Monotone piecewise
!C                 cubic interpolation, SIAM Journal on Numerical Ana-
!C                 lysis 17, 2 (April 1980), pp. 238-246.
!!
!C***ROUTINES CALLED  DPCHCE, DPCHCI, DPCHCS, XERMSG
!!
!C***REVISION HISTORY  (YYMMDD)
!C   820218  DATE WRITTEN
!C   820804  Converted to SLATEC library version.
!C   870707  Corrected XERROR calls for d.p. name(s).
!C   870813  Updated Reference 2.
!C   890206  Corrected XERROR calls.
!C   890411  Added SAVE statements (Vers. 3.2).
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890703  Corrected category record.  (WRB)
!C   890831  Modified array declarations.  (WRB)
!C   891006  Cosmetic changes to prologue.  (WRB)
!C   891006  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!C   920429  Revised format and order of references.  (WRB,FNF)
!C***END PROLOGUE  DPCHIC
!C  Programming notes:
!C
!C     To produce a single precision version, simply:
!C        a. Change DPCHIC to PCHIC wherever it occurs,
!C        b. Change DPCHCE to PCHCE wherever it occurs,
!C        c. Change DPCHCI to PCHCI wherever it occurs,
!C        d. Change DPCHCS to PCHCS wherever it occurs,
!C        e. Change the double precision declarations to real, and
!C        f. Change the constant  ZERO  to single precision.
!C
!C  DECLARE ARGUMENTS.
!C
      INTEGER  IC(2), N, INCFD, NWK, IERR
      DOUBLE PRECISION  VC(2), SWITCH, X(*), F(INCFD,*), D(INCFD,*),    &
     & WK(NWK)
!C
!C  DECLARE LOCAL VARIABLES.
!C
      INTEGER  I, IBEG, IEND, NLESS1
      DOUBLE PRECISION  vsmall
      DATA  vsmall/1.d-15/
!C
!C  VALIDITY-CHECK ARGUMENTS.
!C
!C***FIRST EXECUTABLE STATEMENT  DPCHIC
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
!C
      IBEG = IC(1)
      IEND = IC(2)
      IERR = 0
      IF (ABS(IBEG) .GT. 5)  IERR = IERR - 1
      IF (ABS(IEND) .GT. 5)  IERR = IERR - 2
      IF (IERR .LT. 0)  GO TO 5004
!C
!C  FUNCTION DEFINITION IS OK -- GO ON.
!C
      NLESS1 = N - 1
      IF ( NWK .LT. 2*NLESS1 )  GO TO 5007
!C
!C  SET UP H AND SLOPE ARRAYS.
!C
      DO 20  I = 1, NLESS1
         WK(I) = X(I+1) - X(I)
         WK(NLESS1+I) = (F(1,I+1) - F(1,I)) / WK(I)
   20 CONTINUE
!C
!C  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
!C
      IF (NLESS1 .GT. 1)  GO TO 1000
      D(1,1) = WK(2)
      D(1,N) = WK(2)
      GO TO 3000
!C
!C  NORMAL CASE  (N .GE. 3) .
!C
 1000 CONTINUE
!C
!C  SET INTERIOR DERIVATIVES AND DEFAULT END CONDITIONS.
!C
!C     --------------------------------------
      CALL DPCHCI (N, WK(1:N-1), WK(N:2*N-2), D, INCFD)  
!C     --------------------------------------
!C
!C  SET DERIVATIVES AT POINTS WHERE MONOTONICITY SWITCHES DIRECTION.
!C
      IF (dabs(SWITCH).lt.vsmall)  GO TO 3000

!C     ----------------------------------------------------
      CALL DPCHCS (SWITCH, N, WK(1:N-1), WK(N:2*N-2), D, INCFD, IERR)
!C     ----------------------------------------------------
      IF (IERR .NE. 0)  GO TO 5008
!C
!C  SET END CONDITIONS.
!C
 3000 CONTINUE
      IF ( (IBEG.EQ.0) .AND. (IEND.EQ.0) )  GO TO 5000
!C     -------------------------------------------------------
      CALL DPCHCE (IC, VC, N, X, WK(1:N-1), WK(N:2*N-2), D, INCFD, IERR)
!C     -------------------------------------------------------
      IF (IERR .LT. 0)  GO TO 5009
!C
!C  NORMAL RETURN.
!C
 5000 CONTINUE
      RETURN
!C
!C  ERROR RETURNS.
!C
 5001 CONTINUE
!C     N.LT.2 RETURN.
      IERR = -1
!      CALL XERMSG ('SLATEC', 'DPCHIC',
!     +   'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
      RETURN
!C
 5002 CONTINUE
!C     INCFD.LT.1 RETURN.
      IERR = -2
!      CALL XERMSG ('SLATEC', 'DPCHIC', 'INCREMENT LESS THAN ONE', IERR,
!     +   1)
      RETURN
!C
 5003 CONTINUE
!C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
!      CALL XERMSG ('SLATEC', 'DPCHIC',
!     +   'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)
      RETURN
!C
 5004 CONTINUE
!C     IC OUT OF RANGE RETURN.
      IERR = IERR - 3
!      CALL XERMSG ('SLATEC', 'DPCHIC', 'IC OUT OF RANGE', IERR, 1)
      RETURN
!C
 5007 CONTINUE
!C     NWK .LT. 2*(N-1)  RETURN.
      IERR = -7
!      CALL XERMSG ('SLATEC', 'DPCHIC', 'WORK ARRAY TOO SMALL', IERR, 1)
      RETURN
!C
 5008 CONTINUE
!C     ERROR RETURN FROM DPCHCS.
      IERR = -8
!      CALL XERMSG ('SLATEC', 'DPCHIC', 'ERROR RETURN FROM DPCHCS',
!     +   IERR, 1)
      RETURN
!C
 5009 CONTINUE
!C     ERROR RETURN FROM DPCHCE.
!C   *** THIS CASE SHOULD NEVER OCCUR ***
      IERR = -9
!      CALL XERMSG ('SLATEC', 'DPCHIC', 'ERROR RETURN FROM DPCHCE',
!     +   IERR, 1)
      RETURN
!C------------- LAST LINE OF DPCHIC FOLLOWS -----------------------------
      END

      SUBROUTINE DPCHCE (IC, VC, N, X, H, SLOPE, D, INCFD, IERR)
!C***BEGIN PROLOGUE  DPCHCE
!C***SUBSIDIARY
!C***PURPOSE  Set boundary conditions for DPCHIC
!C***LIBRARY   SLATEC (PCHIP)
!C***TYPE      DOUBLE PRECISION (PCHCE-S, DPCHCE-D)
!C***AUTHOR  Fritsch, F. N., (LLNL)
!C***DESCRIPTION
!C
!C          DPCHCE:  DPCHIC End Derivative Setter.
!C
!C    Called by DPCHIC to set end derivatives as requested by the user.
!C    It must be called after interior derivative values have been set.
!C                      -----
!C
!C    To facilitate two-dimensional applications, includes an increment
!C    between successive values of the D-array.
!C
!C ----------------------------------------------------------------------
!C
!C  Calling sequence:
!C
!C        PARAMETER  (INCFD = ...)
!C        INTEGER  IC(2), N, IERR
!C        DOUBLE PRECISION  VC(2), X(N), H(N), SLOPE(N), D(INCFD,N)
!C
!C        CALL  DPCHCE (IC, VC, N, X, H, SLOPE, D, INCFD, IERR)
!C
!C   Parameters:
!C
!C     IC -- (input) integer array of length 2 specifying desired
!C           boundary conditions:
!C           IC(1) = IBEG, desired condition at beginning of data.
!C           IC(2) = IEND, desired condition at end of data.
!C           ( see prologue to DPCHIC for details. )
!C
!C     VC -- (input) real(8)array of length 2 specifying desired boundary
!C           values.  VC(1) need be set only if IC(1) = 2 or 3 .
!C                    VC(2) need be set only if IC(2) = 2 or 3 .
!C
!C     N -- (input) number of data points.  (assumes N.GE.2)
!C
!C     X -- (input) real(8)array of independent variable values.  (the
!C           elements of X are assumed to be strictly increasing.)
!C
!C     H -- (input) real(8)array of interval lengths.
!C     SLOPE -- (input) real(8)array of data slopes.
!C           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
!C                  H(I) =  X(I+1)-X(I),
!C              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
!C
!C     D -- (input) real(8)array of derivative values at the data points.
!C           The value corresponding to X(I) must be stored in
!C                D(1+(I-1)*INCFD),  I=1(1)N.
!C          (output) the value of D at X(1) and/or X(N) is changed, if
!C           necessary, to produce the requested boundary conditions.
!C           no other entries in D are changed.
!C
!C     INCFD -- (input) increment between successive values in D.
!C           This argument is provided primarily for 2-D applications.
!C
!C     IERR -- (output) error flag.
!C           Normal return:
!C              IERR = 0  (no errors).
!C           Warning errors:
!C              IERR = 1  if IBEG.LT.0 and D(1) had to be adjusted for
!C                        monotonicity.
!C              IERR = 2  if IEND.LT.0 and D(1+(N-1)*INCFD) had to be
!C                        adjusted for monotonicity.
!C              IERR = 3  if both of the above are true.
!C
!C    -------
!C    WARNING:  This routine does no validity-checking of arguments.
!C    -------
!C
!C  Fortran intrinsics used:  ABS.
!C
!C***SEE ALSO  DPCHIC
!C***ROUTINES CALLED  DPCHDF, DPCHST, XERMSG
!C***REVISION HISTORY  (YYMMDD)
!C   820218  DATE WRITTEN
!C   820805  Converted to SLATEC library version.
!C   870707  Corrected XERROR calls for d.p. name(s).
!C   890206  Corrected XERROR calls.
!C   890411  Added SAVE statements (Vers. 3.2).
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890831  Modified array declarations.  (WRB)
!C   890831  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!C   900328  Added TYPE section.  (WRB)
!C   910408  Updated AUTHOR section in prologue.  (WRB)
!C   930503  Improved purpose.  (FNF)
!C***END PROLOGUE  DPCHCE
!C
!C  Programming notes:
!C     1. The function DPCHST(ARG1,ARG2)  is assumed to return zero if
!C        either argument is zero, +1 if they are of the same sign, and
!C        -1 if they are of opposite sign.
!C     2. One could reduce the number of arguments and amount of local
!C        storage, at the expense of reduced code clarity, by passing in
!C        the array WK (rather than splitting it into H and SLOPE) and
!C        increasing its length enough to incorporate STEMP and XTEMP.
!C     3. The two monotonicity checks only use the sufficient conditions.
!C        Thus, it is possible (but unlikely) for a boundary condition to
!C        be changed, even though the original interpolant was monotonic.
!C        (At least the result is a continuous function of the data.)
!C**End
!C
!C  DECLARE ARGUMENTS.
!C
      INTEGER  IC(2), N, INCFD, IERR
      DOUBLE PRECISION  VC(2), X(*), H(*), SLOPE(*), D(INCFD,*)
!C
!C  DECLARE LOCAL VARIABLES.
!C
      INTEGER  IBEG, IEND, IERF, INDEX, J, K
      DOUBLE PRECISION  HALF, STEMP(3), THREE, TWO, XTEMP(4), ZERO
      SAVE ZERO, HALF, TWO, THREE
      DOUBLE PRECISION  DPCHDF, DPCHST, vsmall
!C
!C  INITIALIZE.
!C
      DATA  ZERO /0.D0/,  HALF/.5D0/,  TWO/2.D0/, THREE/3.D0/
      DATA  vsmall /1.D-15/
!C
!C***FIRST EXECUTABLE STATEMENT  DPCHCE
      IBEG = IC(1)
      IEND = IC(2)
      IERR = 0
!C
!C  SET TO DEFAULT BOUNDARY CONDITIONS IF N IS TOO SMALL.
!C
      IF ( ABS(IBEG).GT.N )  IBEG = 0
      IF ( ABS(IEND).GT.N )  IEND = 0
!C
!C  TREAT BEGINNING BOUNDARY CONDITION.
!C
      IF (IBEG .EQ. 0)  GO TO 2000
      K = ABS(IBEG)
      IF (K .EQ. 1)  THEN
!C        BOUNDARY VALUE PROVIDED.
         D(1,1) = VC(1)
      ELSE IF (K .EQ. 2)  THEN
!C        BOUNDARY SECOND DERIVATIVE PROVIDED.
         D(1,1) = HALF*( (THREE*SLOPE(1) - D(1,2)) - HALF*VC(1)*H(1) )
      ELSE IF (K .LT. 5)  THEN
!C        USE K-POINT DERIVATIVE FORMULA.
!C        PICK UP FIRST K POINTS, IN REVERSE ORDER.
         DO 10  J = 1, K
            INDEX = K-J+1
!C           INDEX RUNS FROM K DOWN TO 1.
            XTEMP(J) = X(INDEX)
            IF (J .LT. K)  STEMP(J) = SLOPE(INDEX-1)
   10    CONTINUE
!C                 -----------------------------
         D(1,1) = DPCHDF (K, XTEMP, STEMP, IERF)
!C                 -----------------------------
         IF (IERF .NE. 0)  GO TO 5001
      ELSE
!C        USE 'NOT A KNOT' CONDITION.
         D(1,1) = ( THREE*(H(1)*SLOPE(2) + H(2)*SLOPE(1))               &
     &             - TWO*(H(1)+H(2))*D(1,2) - H(1)*D(1,3) ) / H(2)
      ENDIF
!C
      IF (IBEG .GT. 0)  GO TO 2000
!C
!C  CHECK D(1,1) FOR COMPATIBILITY WITH MONOTONICITY.
!C
      IF (dabs(SLOPE(1)).lt.vsmall)  THEN
         IF (dabs(D(1,1)).gt.ZERO)  THEN
            D(1,1) = ZERO
            IERR = IERR + 1
         ENDIF
      ELSE IF ( DPCHST(D(1,1),SLOPE(1)) .LT. ZERO)  THEN
         D(1,1) = ZERO
         IERR = IERR + 1
      ELSE IF ( ABS(D(1,1)) .GT. THREE*ABS(SLOPE(1)) )  THEN
         D(1,1) = THREE*SLOPE(1)
         IERR = IERR + 1
      ENDIF
!C
!C  TREAT END BOUNDARY CONDITION.
!C
 2000 CONTINUE
      IF (IEND .EQ. 0)  GO TO 5000
      K = ABS(IEND)
      IF (K .EQ. 1)  THEN
!C        BOUNDARY VALUE PROVIDED.
         D(1,N) = VC(2)
      ELSE IF (K .EQ. 2)  THEN
!C        BOUNDARY SECOND DERIVATIVE PROVIDED.
         D(1,N) = HALF*( (THREE*SLOPE(N-1) - D(1,N-1)) +                &
     &                                           HALF*VC(2)*H(N-1) )
      ELSE IF (K .LT. 5)  THEN
!C        USE K-POINT DERIVATIVE FORMULA.
!C        PICK UP LAST K POINTS.
         DO 2010  J = 1, K
            INDEX = N-K+J
!C           INDEX RUNS FROM N+1-K UP TO N.
            XTEMP(J) = X(INDEX)
            IF (J .LT. K)  STEMP(J) = SLOPE(INDEX)
 2010    CONTINUE
!C                 -----------------------------
         D(1,N) = DPCHDF (K, XTEMP, STEMP, IERF)
!C                 -----------------------------
         IF (IERF .NE. 0)  GO TO 5001
      ELSE
!C        USE 'NOT A KNOT' CONDITION.
         D(1,N) = ( THREE*(H(N-1)*SLOPE(N-2) + H(N-2)*SLOPE(N-1))       &
     &             - TWO*(H(N-1)+H(N-2))*D(1,N-1) - H(N-1)*D(1,N-2) )   &
     &                                                         / H(N-2)
      ENDIF
!C
      IF (IEND .GT. 0)  GO TO 5000
!C
!C  CHECK D(1,N) FOR COMPATIBILITY WITH MONOTONICITY.
!C

      IF (dabs(SLOPE(N-1)).lt.vsmall)  THEN
         IF (dabs(D(1,N)).gt.ZERO)  THEN
            D(1,N) = ZERO
            IERR = IERR + 2
         ENDIF
      ELSE IF ( DPCHST(D(1,N),SLOPE(N-1)) .LT. ZERO)  THEN
         D(1,N) = ZERO
         IERR = IERR + 2
      ELSE IF ( ABS(D(1,N)) .GT. THREE*ABS(SLOPE(N-1)) )  THEN
         D(1,N) = THREE*SLOPE(N-1)
         IERR = IERR + 2
      ENDIF
!C
!C  NORMAL RETURN.
!C
 5000 CONTINUE
      RETURN
!C
!C  ERROR RETURN.
!C
 5001 CONTINUE
!C     ERROR RETURN FROM DPCHDF.
!C   *** THIS CASE SHOULD NEVER OCCUR ***
      IERR = -1
!      CALL XERMSG ('SLATEC', 'DPCHCE', 'ERROR RETURN FROM DPCHDF',
!     +   IERR, 1)
      RETURN
!C------------- LAST LINE OF DPCHCE FOLLOWS -----------------------------
      END

      SUBROUTINE DPCHCS (SWITCH, N, H, SLOPE, D, INCFD, IERR)
!C***BEGIN PROLOGUE  DPCHCS
!C***SUBSIDIARY
!C***PURPOSE  Adjusts derivative values for DPCHIC
!C***LIBRARY   SLATEC (PCHIP)
!C***TYPE      DOUBLE PRECISION (PCHCS-S, DPCHCS-D)
!C***AUTHOR  Fritsch, F. N., (LLNL)
!C***DESCRIPTION
!C
!C         DPCHCS:  DPCHIC Monotonicity Switch Derivative Setter.
!C
!C     Called by  DPCHIC  to adjust the values of D in the vicinity of a
!C     switch in direction of monotonicity, to produce a more "visually
!C     pleasing" curve than that given by  DPCHIM .
!C
!C ----------------------------------------------------------------------
!C
!C  Calling sequence:
!C
!C        PARAMETER  (INCFD = ...)
!C        INTEGER  N, IERR
!C        DOUBLE PRECISION  SWITCH, H(N), SLOPE(N), D(INCFD,N)
!C
!C        CALL  DPCHCS (SWITCH, N, H, SLOPE, D, INCFD, IERR)
!C
!C   Parameters:
!C
!C     SWITCH -- (input) indicates the amount of control desired over
!C           local excursions from data.
!C
!C     N -- (input) number of data points.  (assumes N.GT.2 .)
!C
!C     H -- (input) real(8)array of interval lengths.
!C     SLOPE -- (input) real(8)array of data slopes.
!C           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
!C                  H(I) =  X(I+1)-X(I),
!C              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
!C
!C     D -- (input) real(8)array of derivative values at the data points,
!C           as determined by DPCHCI.
!C          (output) derivatives in the vicinity of switches in direction
!C           of monotonicity may be adjusted to produce a more "visually
!C           pleasing" curve.
!C           The value corresponding to X(I) is stored in
!C                D(1+(I-1)*INCFD),  I=1(1)N.
!C           No other entries in D are changed.
!C
!C     INCFD -- (input) increment between successive values in D.
!C           This argument is provided primarily for 2-D applications.
!C
!C     IERR -- (output) error flag.  should be zero.
!C           If negative, trouble in DPCHSW.  (should never happen.)
!C
!C    -------
!C    WARNING:  This routine does no validity-checking of arguments.
!C    -------
!C
!C  Fortran intrinsics used:  ABS, MAX, MIN.
!C
!C***SEE ALSO  DPCHIC
!C***ROUTINES CALLED  DPCHST, DPCHSW
!C***REVISION HISTORY  (YYMMDD)
!C   820218  DATE WRITTEN
!C   820617  Redesigned to (1) fix  problem with lack of continuity
!C           approaching a flat-topped peak (2) be cleaner and
!C           easier to verify.
!C           Eliminated subroutines PCHSA and PCHSX in the process.
!C   820622  1. Limited fact to not exceed one, so computed D is a
!C             convex combination of DPCHCI value and DPCHSD value.
!C           2. Changed fudge from 1 to 4 (based on experiments).
!C   820623  Moved PCHSD to an inline function (eliminating MSWTYP).
!C   820805  Converted to SLATEC library version.
!C   870707  Corrected conversion to double precision.
!C   870813  Minor cosmetic changes.
!C   890411  Added SAVE statements (Vers. 3.2).
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890831  Modified array declarations.  (WRB)
!C   891006  Modified spacing in computation of DFLOC.  (WRB)
!C   891006  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900328  Added TYPE section.  (WRB)
!C   910408  Updated AUTHOR section in prologue.  (WRB)
!C   930503  Improved purpose.  (FNF)
!C***END PROLOGUE  DPCHCS
!C
!C  Programming notes:
!C     1. The function  DPCHST(ARG1,ARG2)  is assumed to return zero if
!C        either argument is zero, +1 if they are of the same sign, and
!C        -1 if they are of opposite sign.
!C**End
!C
!C  DECLARE ARGUMENTS.
!C
      INTEGER  N, INCFD, IERR
      DOUBLE PRECISION  SWITCH, H(*), SLOPE(*), D(INCFD,*)
!C
!C  DECLARE LOCAL VARIABLES.
!C
      INTEGER  I, INDX, K, NLESS1
      DOUBLE PRECISION  DEL(3), DEXT, DFLOC, DFMX, FACT, FUDGE, ONE,    &
     &      SLMAX, WTAVE(2), ZERO
      SAVE ZERO, ONE, FUDGE
      DOUBLE PRECISION  DPCHST
!C
!C  DEFINE INLINE FUNCTION FOR WEIGHTED AVERAGE OF SLOPES.
!C
      DOUBLE PRECISION  DPCHSD, S1, S2, H1, H2
      DPCHSD(S1,S2,H1,H2) = (H2/(H1+H2))*S1 + (H1/(H1+H2))*S2
!C
!C  INITIALIZE.
!C
      DATA  ZERO /0.D0/,  ONE/1.D0/
      DATA  FUDGE /4.D0/
!C***FIRST EXECUTABLE STATEMENT  DPCHCS
      IERR = 0
      NLESS1 = N - 1
!C
!C  LOOP OVER SEGMENTS.
!C
      DO 900  I = 2, NLESS1
! obsolete arithmetic if         IF ( DPCHST(SLOPE(I-1),SLOPE(I)) )  100, 300, 900
         if (DPCHST(SLOPE(I-1),SLOPE(I)) < 0.d0) then
            goto 100
         else if (DPCHST(SLOPE(I-1),SLOPE(I)) > 0.d0) then
            goto 900
         else
            goto 300
         end if
!C             --------------------------
!C
  100    CONTINUE
!C
!C....... SLOPE SWITCHES MONOTONICITY AT I-TH POINT .....................
!C
!C           DO NOT CHANGE D IF 'UP-DOWN-UP'.
            IF (I .GT. 2)  THEN
               IF ( DPCHST(SLOPE(I-2),SLOPE(I)) .GT. ZERO)  GO TO 900
!C                   --------------------------
            ENDIF
            IF (I .LT. NLESS1)  THEN
               IF ( DPCHST(SLOPE(I+1),SLOPE(I-1)) .GT. ZERO)  GO TO 900
!C                   ----------------------------
            ENDIF
!C
!C   ....... COMPUTE PROVISIONAL VALUE FOR D(1,I).
!C
            DEXT = DPCHSD (SLOPE(I-1), SLOPE(I), H(I-1), H(I))
!C
!C   ....... DETERMINE WHICH INTERVAL CONTAINS THE EXTREMUM.
!C
! obsolete arithmetic if            IF ( DPCHST(DEXT, SLOPE(I-1)) )  200, 900, 250
            if (DPCHST(DEXT, SLOPE(I-1)) < 0.d0) then
               goto 200
            else if (DPCHST(DEXT, SLOPE(I-1)) > 0.d0) then
               goto 250
            else
               goto 900
            end if
!C                -----------------------
!C
  200       CONTINUE
!C              DEXT AND SLOPE(I-1) HAVE OPPOSITE SIGNS --
!C                        EXTREMUM IS IN (X(I-1),X(I)).
               K = I-1
!C              SET UP TO COMPUTE NEW VALUES FOR D(1,I-1) AND D(1,I).
               WTAVE(2) = DEXT
               IF (K .GT. 1)                                            &
     &            WTAVE(1) = DPCHSD (SLOPE(K-1), SLOPE(K), H(K-1), H(K))
               GO TO 400
!C
  250       CONTINUE
!C              DEXT AND SLOPE(I) HAVE OPPOSITE SIGNS --
!C                        EXTREMUM IS IN (X(I),X(I+1)).
               K = I
!C              SET UP TO COMPUTE NEW VALUES FOR D(1,I) AND D(1,I+1).
               WTAVE(1) = DEXT
               IF (K .LT. NLESS1)                                       &
     &            WTAVE(2) = DPCHSD (SLOPE(K), SLOPE(K+1), H(K), H(K+1))
               GO TO 400
!C
  300    CONTINUE
!C
!C....... AT LEAST ONE OF SLOPE(I-1) AND SLOPE(I) IS ZERO --
!C                     CHECK FOR FLAT-TOPPED PEAK .......................
!C
            IF (I .EQ. NLESS1)  GO TO 900
            IF ( DPCHST(SLOPE(I-1), SLOPE(I+1)) .GE. ZERO)  GO TO 900
!C                -----------------------------
!C
!C           WE HAVE FLAT-TOPPED PEAK ON (X(I),X(I+1)).
            K = I
!C           SET UP TO COMPUTE NEW VALUES FOR D(1,I) AND D(1,I+1).
            WTAVE(1) = DPCHSD (SLOPE(K-1), SLOPE(K), H(K-1), H(K))
            WTAVE(2) = DPCHSD (SLOPE(K), SLOPE(K+1), H(K), H(K+1))
!C
  400    CONTINUE
!C
!C....... AT THIS POINT WE HAVE DETERMINED THAT THERE WILL BE AN EXTREMUM
!C        ON (X(K),X(K+1)), WHERE K=I OR I-1, AND HAVE SET ARRAY WTAVE--
!C           WTAVE(1) IS A WEIGHTED AVERAGE OF SLOPE(K-1) AND SLOPE(K),
!C                    IF K.GT.1
!C           WTAVE(2) IS A WEIGHTED AVERAGE OF SLOPE(K) AND SLOPE(K+1),
!C                    IF K.LT.N-1
!C
         SLMAX = ABS(SLOPE(K))
         IF (K .GT. 1)    SLMAX = MAX( SLMAX, ABS(SLOPE(K-1)) )
         IF (K.LT.NLESS1) SLMAX = MAX( SLMAX, ABS(SLOPE(K+1)) )
!C
         IF (K .GT. 1)  DEL(1) = SLOPE(K-1) / SLMAX
         DEL(2) = SLOPE(K) / SLMAX
         IF (K.LT.NLESS1)  DEL(3) = SLOPE(K+1) / SLMAX
!C
         IF ((K.GT.1) .AND. (K.LT.NLESS1))  THEN
!C           NORMAL CASE -- EXTREMUM IS NOT IN A BOUNDARY INTERVAL.
            FACT = FUDGE* ABS(DEL(3)*(DEL(1)-DEL(2))*(WTAVE(2)/SLMAX))
            D(1,K) = D(1,K) + MIN(FACT,ONE)*(WTAVE(1) - D(1,K))
            FACT = FUDGE* ABS(DEL(1)*(DEL(3)-DEL(2))*(WTAVE(1)/SLMAX))
            D(1,K+1) = D(1,K+1) + MIN(FACT,ONE)*(WTAVE(2) - D(1,K+1))
         ELSE
!C           SPECIAL CASE K=1 (WHICH CAN OCCUR ONLY IF I=2) OR
!C                        K=NLESS1 (WHICH CAN OCCUR ONLY IF I=NLESS1).
            FACT = FUDGE* ABS(DEL(2))
            D(1,I) = MIN(FACT,ONE) * WTAVE(I-K+1)
!C              NOTE THAT I-K+1 = 1 IF K=I  (=NLESS1),
!C                        I-K+1 = 2 IF K=I-1(=1).
         ENDIF
!C
!C
!C....... ADJUST IF NECESSARY TO LIMIT EXCURSIONS FROM DATA.
!C
         IF (SWITCH .LE. ZERO)  GO TO 900
!C
         DFLOC = H(K)*ABS(SLOPE(K))
         IF (K .GT. 1)    DFLOC = MAX( DFLOC, H(K-1)*ABS(SLOPE(K-1)) )
         IF (K.LT.NLESS1) DFLOC = MAX( DFLOC, H(K+1)*ABS(SLOPE(K+1)) )
         DFMX = SWITCH*DFLOC
         INDX = I-K+1
!C        INDX = 1 IF K=I, 2 IF K=I-1.
!C        ---------------------------------------------------------------
         CALL DPCHSW(DFMX, INDX, D(1,K), D(1,K+1), H(K), SLOPE(K), IERR)
!C        ---------------------------------------------------------------
         IF (IERR .NE. 0)  RETURN
!C
!C....... END OF SEGMENT LOOP.
!C
  900 CONTINUE
!C
      RETURN
!C------------- LAST LINE OF DPCHCS FOLLOWS -----------------------------
      END

      SUBROUTINE DPCHCI (N, H, SLOPE, D, INCFD)
!C***BEGIN PROLOGUE  DPCHCI
!C***SUBSIDIARY
!C***PURPOSE  Set interior derivatives for DPCHIC
!C***LIBRARY   SLATEC (PCHIP)
!C***TYPE      DOUBLE PRECISION (PCHCI-S, DPCHCI-D)
!C***AUTHOR  Fritsch, F. N., (LLNL)
!C***DESCRIPTION
!C
!C          DPCHCI:  DPCHIC Initial Derivative Setter.
!C
!C    Called by DPCHIC to set derivatives needed to determine a monotone
!C    piecewise cubic Hermite interpolant to the data.
!C
!C    Default boundary conditions are provided which are compatible
!C    with monotonicity.  If the data are only piecewise monotonic, the
!C    interpolant will have an extremum at each point where monotonicity
!C    switches direction.
!C
!C    To facilitate two-dimensional applications, includes an increment
!C    between successive values of the D-array.
!C
!C    The resulting piecewise cubic Hermite function should be identical
!C    (within roundoff error) to that produced by DPCHIM.
!C
!C ----------------------------------------------------------------------
!C
!C  Calling sequence:
!C
!C        PARAMETER  (INCFD = ...)
!C        INTEGER  N
!C        DOUBLE PRECISION  H(N), SLOPE(N), D(INCFD,N)
!C
!C        CALL  DPCHCI (N, H, SLOPE, D, INCFD)
!C
!C   Parameters:
!C
!C     N -- (input) number of data points.
!C           If N=2, simply does linear interpolation.
!C
!C     H -- (input) real(8)array of interval lengths.
!C     SLOPE -- (input) real(8)array of data slopes.
!C           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
!C                  H(I) =  X(I+1)-X(I),
!C              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
!C
!C     D -- (output) real(8)array of derivative values at data points.
!C           If the data are monotonic, these values will determine a
!C           a monotone cubic Hermite function.
!C           The value corresponding to X(I) is stored in
!C                D(1+(I-1)*INCFD),  I=1(1)N.
!C           No other entries in D are changed.
!C
!C     INCFD -- (input) increment between successive values in D.
!C           This argument is provided primarily for 2-D applications.
!C
!C    -------
!C    WARNING:  This routine does no validity-checking of arguments.
!C    -------
!C
!C  Fortran intrinsics used:  ABS, MAX, MIN.
!C
!C***SEE ALSO  DPCHIC
!C***ROUTINES CALLED  DPCHST
!C***REVISION HISTORY  (YYMMDD)
!C   820218  DATE WRITTEN
!C   820601  Modified end conditions to be continuous functions of
!C           data when monotonicity switches in next interval.
!C   820602  1. Modified formulas so end conditions are less prone
!C             to over/underflow problems.
!C           2. Minor modification to HSUM calculation.
!C   820805  Converted to SLATEC library version.
!C   870813  Minor cosmetic changes.
!C   890411  Added SAVE statements (Vers. 3.2).
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890831  Modified array declarations.  (WRB)
!C   890831  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900328  Added TYPE section.  (WRB)
!C   910408  Updated AUTHOR section in prologue.  (WRB)
!C   930503  Improved purpose.  (FNF)
!C***END PROLOGUE  DPCHCI
!C
!C  Programming notes:
!C     1. The function  DPCHST(ARG1,ARG2)  is assumed to return zero if
!C        either argument is zero, +1 if they are of the same sign, and
!C        -1 if they are of opposite sign.
!C**End
!C
!C  DECLARE ARGUMENTS.
!C
      INTEGER  N, INCFD
      DOUBLE PRECISION  H(*), SLOPE(*), D(INCFD,*)
!C
!C  DECLARE LOCAL VARIABLES.
!C
      INTEGER  I, NLESS1
      DOUBLE PRECISION  DEL1, DEL2, DMAX, DMIN, DRAT1, DRAT2, HSUM,     &
     &      HSUMT3, THREE, W1, W2, ZERO
      SAVE ZERO, THREE
      DOUBLE PRECISION  DPCHST
!C
!C  INITIALIZE.
!C
      DATA  ZERO /0.D0/, THREE/3.D0/
!C***FIRST EXECUTABLE STATEMENT  DPCHCI
      NLESS1 = N - 1
      DEL1 = SLOPE(1)
!C
!C  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
!C
      IF (NLESS1 .GT. 1)  GO TO 10
      D(1,1) = DEL1
      D(1,N) = DEL1
      GO TO 5000
!C
!C  NORMAL CASE  (N .GE. 3).
!C
   10 CONTINUE
      DEL2 = SLOPE(2)
!C
!C  SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!C     SHAPE-PRESERVING.
!C
      HSUM = H(1) + H(2)
      W1 = (H(1) + HSUM)/HSUM
      W2 = -H(1)/HSUM
      D(1,1) = W1*DEL1 + W2*DEL2
      IF ( DPCHST(D(1,1),DEL1) .LE. ZERO)  THEN
         D(1,1) = ZERO
      ELSE IF ( DPCHST(DEL1,DEL2) .LT. ZERO)  THEN
!C        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
         DMAX = THREE*DEL1
         IF (ABS(D(1,1)) .GT. ABS(DMAX))  D(1,1) = DMAX
      ENDIF
!C
!C  LOOP THROUGH INTERIOR POINTS.
!C
      DO 50  I = 2, NLESS1
         IF (I .EQ. 2)  GO TO 40
!C
         HSUM = H(I-1) + H(I)
         DEL1 = DEL2
         DEL2 = SLOPE(I)
   40    CONTINUE
!C
!C        SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.
!C
         D(1,I) = ZERO
         IF ( DPCHST(DEL1,DEL2) .LE. ZERO)  GO TO 50
!C
!C        USE BRODLIE MODIFICATION OF BUTLAND FORMULA.
!C
         HSUMT3 = HSUM+HSUM+HSUM
         W1 = (HSUM + H(I-1))/HSUMT3
         W2 = (HSUM + H(I)  )/HSUMT3
         DMAX = MAX( ABS(DEL1), ABS(DEL2) )
         DMIN = MIN( ABS(DEL1), ABS(DEL2) )
         DRAT1 = DEL1/DMAX
         DRAT2 = DEL2/DMAX
         D(1,I) = DMIN/(W1*DRAT1 + W2*DRAT2)
!C
   50 CONTINUE
!C
!C  SET D(N) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!C     SHAPE-PRESERVING.
!C
      W1 = -H(N-1)/HSUM
      W2 = (H(N-1) + HSUM)/HSUM
      D(1,N) = W1*DEL1 + W2*DEL2
      IF ( DPCHST(D(1,N),DEL2) .LE. ZERO)  THEN
         D(1,N) = ZERO
      ELSE IF ( DPCHST(DEL1,DEL2) .LT. ZERO)  THEN
!C        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
         DMAX = THREE*DEL2
         IF (ABS(D(1,N)) .GT. ABS(DMAX))  D(1,N) = DMAX
      ENDIF
!C
!C  NORMAL RETURN.
!C
 5000 CONTINUE
      RETURN
!C------------- LAST LINE OF DPCHCI FOLLOWS -----------------------------
      END

      DOUBLE PRECISION FUNCTION DPCHDF (K, X, S, IERR)
!C***BEGIN PROLOGUE  DPCHDF
!C***SUBSIDIARY
!C***PURPOSE  Computes divided differences for DPCHCE and DPCHSP
!C***LIBRARY   SLATEC (PCHIP)
!C***TYPE      DOUBLE PRECISION (PCHDF-S, DPCHDF-D)
!C***AUTHOR  Fritsch, F. N., (LLNL)
!C***DESCRIPTION
!C
!C          DPCHDF:   DPCHIP Finite Difference Formula
!C
!C     Uses a divided difference formulation to compute a K-point approx-
!C     imation to the derivative at X(K) based on the data in X and S.
!C
!C     Called by  DPCHCE  and  DPCHSP  to compute 3- and 4-point boundary
!C     derivative approximations.
!C
!C ----------------------------------------------------------------------
!C
!C     On input:
!C        K      is the order of the desired derivative approximation.
!C               K must be at least 3 (error return if not).
!C        X      contains the K values of the independent variable.
!C               X need not be ordered, but the values **MUST** be
!C               distinct.  (Not checked here.)
!C        S      contains the associated slope values:
!C                  S(I) = (F(I+1)-F(I))/(X(I+1)-X(I)), I=1(1)K-1.
!C               (Note that S need only be of length K-1.)
!C
!C     On return:
!C        S      will be destroyed.
!C        IERR   will be set to -1 if K.LT.2 .
!C        DPCHDF  will be set to the desired derivative approximation if
!C               IERR=0 or to zero if IERR=-1.
!C
!C ----------------------------------------------------------------------
!C
!C***SEE ALSO  DPCHCE, DPCHSP
!C***REFERENCES  Carl de Boor, A Practical Guide to Splines, Springer-
!C                 Verlag, New York, 1978, pp. 10-16.
!C***ROUTINES CALLED  XERMSG
!C***REVISION HISTORY  (YYMMDD)
!C   820503  DATE WRITTEN
!C   820805  Converted to SLATEC library version.
!C   870707  Corrected XERROR calls for d.p. name(s).
!C   870813  Minor cosmetic changes.
!C   890206  Corrected XERROR calls.
!C   890411  Added SAVE statements (Vers. 3.2).
!C   890411  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!C   900328  Added TYPE section.  (WRB)
!C   910408  Updated AUTHOR and DATE WRITTEN sections in prologue.  (WRB)
!C   920429  Revised format and order of references.  (WRB,FNF)
!C   930503  Improved purpose.  (FNF)
!C***END PROLOGUE  DPCHDF
!C
!C**End
!C
!C  DECLARE ARGUMENTS.
!C
      INTEGER  K, IERR
      DOUBLE PRECISION  X(K), S(K)
!C
!C  DECLARE LOCAL VARIABLES.
!C
      INTEGER  I, J
      DOUBLE PRECISION  VALUE, ZERO
      SAVE ZERO
      DATA  ZERO /0.D0/
!C
!C  CHECK FOR LEGAL VALUE OF K.
!C
!C***FIRST EXECUTABLE STATEMENT  DPCHDF
      IF (K .LT. 3)  GO TO 5001
!C
!C  COMPUTE COEFFICIENTS OF INTERPOLATING POLYNOMIAL.
!C
      DO 10  J = 2, K-1
         DO 9  I = 1, K-J
            S(I) = (S(I+1)-S(I))/(X(I+J)-X(I))
    9    CONTINUE
   10 CONTINUE
!C
!C  EVALUATE DERIVATIVE AT X(K).
!C
      VALUE = S(1)
      DO 20  I = 2, K-1
         VALUE = S(I) + VALUE*(X(K)-X(I))
   20 CONTINUE
!C
!C  NORMAL RETURN.
!C
      IERR = 0
      DPCHDF = VALUE
      RETURN
!C
!C  ERROR RETURN.
!C
 5001 CONTINUE
!C     K.LT.3 RETURN.
      IERR = -1
!      CALL XERMSG ('SLATEC', 'DPCHDF', 'K LESS THAN THREE', IERR, 1)
      DPCHDF = ZERO
      RETURN
!C------------- LAST LINE OF DPCHDF FOLLOWS -----------------------------
      END


      DOUBLE PRECISION FUNCTION DPCHST (ARG1, ARG2)
!C***BEGIN PROLOGUE  DPCHST
!C***SUBSIDIARY
!C***PURPOSE  DPCHIP Sign-Testing Routine
!C***LIBRARY   SLATEC (PCHIP)
!C***TYPE      DOUBLE PRECISION (PCHST-S, DPCHST-D)
!C***AUTHOR  Fritsch, F. N., (LLNL)
!C***DESCRIPTION
!C
!C         DPCHST:  DPCHIP Sign-Testing Routine.
!C
!C
!C     Returns:
!C        -1. if ARG1 and ARG2 are of opposite sign.
!C         0. if either argument is zero.
!C        +1. if ARG1 and ARG2 are of the same sign.
!C
!C     The object is to do this without multiplying ARG1*ARG2, to avoid
!C     possible over/underflow problems.
!C
!C  Fortran intrinsics used:  SIGN.
!C
!C***SEE ALSO  DPCHCE, DPCHCI, DPCHCS, DPCHIM
!C***ROUTINES CALLED  (NONE)
!C***REVISION HISTORY  (YYMMDD)
!C   811103  DATE WRITTEN
!C   820805  Converted to SLATEC library version.
!C   870813  Minor cosmetic changes.
!C   890411  Added SAVE statements (Vers. 3.2).
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890531  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900328  Added TYPE section.  (WRB)
!C   910408  Updated AUTHOR and DATE WRITTEN sections in prologue.  (WRB)
!C   930503  Improved purpose.  (FNF)
!C***END PROLOGUE  DPCHST
!C
!C**End
!C
!C  DECLARE ARGUMENTS.
!C
      DOUBLE PRECISION  ARG1, ARG2
!C
!C  DECLARE LOCAL VARIABLES.
!C
      DOUBLE PRECISION  ONE, ZERO, vsmall
      SAVE ZERO, ONE
      DATA  ZERO /0.D0/,  vsmall/1.D-15/,  ONE/1.D0/
      
!C
!C  PERFORM THE TEST.
!C
!C***FIRST EXECUTABLE STATEMENT  DPCHST
      DPCHST = SIGN(ONE,ARG1) * SIGN(ONE,ARG2)
      IF ((dabs(ARG1).lt.vsmall) .OR. (dabs(ARG2).lt.vsmall))           &
     &     DPCHST = ZERO
      
!C
      RETURN
!C------------- LAST LINE OF DPCHST FOLLOWS -----------------------------
      END


      SUBROUTINE DPCHSW (DFMAX, IEXTRM, D1, D2, H, SLOPE, IERR)
!C***BEGIN PROLOGUE  DPCHSW
!C***SUBSIDIARY
!C***PURPOSE  Limits excursion from data for DPCHCS
!C***LIBRARY   SLATEC (PCHIP)
!C***TYPE      DOUBLE PRECISION (PCHSW-S, DPCHSW-D)
!C***AUTHOR  Fritsch, F. N., (LLNL)
!C***DESCRIPTION
!C
!C         DPCHSW:  DPCHCS Switch Excursion Limiter.
!C
!C     Called by  DPCHCS  to adjust D1 and D2 if necessary to insure that
!C     the extremum on this interval is not further than DFMAX from the
!C     extreme data value.
!C
!C ----------------------------------------------------------------------
!C
!C  Calling sequence:
!C
!C        INTEGER  IEXTRM, IERR
!C        DOUBLE PRECISION  DFMAX, D1, D2, H, SLOPE
!C
!C        CALL  DPCHSW (DFMAX, IEXTRM, D1, D2, H, SLOPE, IERR)
!C
!C   Parameters:
!C
!C     DFMAX -- (input) maximum allowed difference between F(IEXTRM) and
!C           the cubic determined by derivative values D1,D2.  (assumes
!C           DFMAX.GT.0.)
!C
!C     IEXTRM -- (input) index of the extreme data value.  (assumes
!C           IEXTRM = 1 or 2 .  Any value .NE.1 is treated as 2.)
!C
!C     D1,D2 -- (input) derivative values at the ends of the interval.
!C           (Assumes D1*D2 .LE. 0.)
!C          (output) may be modified if necessary to meet the restriction
!C           imposed by DFMAX.
!C
!C     H -- (input) interval length.  (Assumes  H.GT.0.)
!C
!C     SLOPE -- (input) data slope on the interval.
!C
!C     IERR -- (output) error flag.  should be zero.
!C           If IERR=-1, assumption on D1 and D2 is not satisfied.
!C           If IERR=-2, quadratic equation locating extremum has
!C                       negative discriminant (should never occur).
!C
!C    -------
!C    WARNING:  This routine does no validity-checking of arguments.
!C    -------
!C
!C  Fortran intrinsics used:  ABS, SIGN, SQRT.
!C
!C  ***SEE ALSO  DPCHCS
!C***ROUTINES CALLED  D1MACH, XERMSG
!C***REVISION HISTORY  (YYMMDD)
!C   820218  DATE WRITTEN
!C   820805  Converted to SLATEC library version.
!C   870707  Corrected XERROR calls for d.p. name(s).
!C   870707  Replaced DATA statement for SMALL with a use of D1MACH.
!C   870813  Minor cosmetic changes.
!C   890206  Corrected XERROR calls.
!C   890411  1. Added SAVE statements (Vers. 3.2).
!C           2. Added DOUBLE PRECISION declaration for D1MACH.
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890531  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!C   900328  Added TYPE section.  (WRB)
!C   910408  Updated AUTHOR and DATE WRITTEN sections in prologue.  (WRB)
!C   920526  Eliminated possible divide by zero problem.  (FNF)
!C   930503  Improved purpose.  (FNF)
!C***END PROLOGUE  DPCHSW
!C
!C**End
!C
!C  DECLARE ARGUMENTS.
!C
      INTEGER  IEXTRM, IERR
      DOUBLE PRECISION  DFMAX, D1, D2, H, SLOPE
!C
!C  DECLARE LOCAL VARIABLES.
!C
      DOUBLE PRECISION  CP, FACT, HPHI, LAMBDA, NU, ONE, PHI, RADCAL,   &
     &                  RHO, SIGMA, SMALL, THAT, THIRD, THREE, TWO, ZERO
      SAVE ZERO, ONE, TWO, THREE, FACT
      SAVE THIRD
      DOUBLE PRECISION  vsmall
!C
      DATA  ZERO /0.D0/,  ONE /1.D0/,  TWO /2.D0/, THREE /3.D0/,        &
     &      FACT /100.D0/, vsmall/1.d-15/
!C        THIRD SHOULD BE SLIGHTLY LESS THAN 1/3.
      DATA  THIRD /0.33333D0/
!C
!C  NOTATION AND GENERAL REMARKS.
!C
!C     RHO IS THE RATIO OF THE DATA SLOPE TO THE DERIVATIVE BEING TESTED.
!C     LAMBDA IS THE RATIO OF D2 TO D1.
!C     THAT = T-HAT(RHO) IS THE NORMALIZED LOCATION OF THE EXTREMUM.
!C     PHI IS THE NORMALIZED VALUE OF P(X)-F1 AT X = XHAT = X-HAT(RHO),
!C           WHERE  THAT = (XHAT - X1)/H .
!C        THAT IS, P(XHAT)-F1 = D*H*PHI,  WHERE D=D1 OR D2.
!C     SIMILARLY,  P(XHAT)-F2 = D*H*(PHI-RHO) .
!C
!C      SMALL SHOULD BE A FEW ORDERS OF MAGNITUDE GREATER THAN MACHEPS.
!C***FIRST EXECUTABLE STATEMENT  DPCHSW
!      SMALL = FACT*D1MACH(4)
      SMALL = FACT*1.d-16
!C
!C  DO MAIN CALCULATION.
!C
      IF (abs(D1) .lt. vsmall)  THEN  
!C
!C        SPECIAL CASE -- D1.EQ.ZERO .
!C
!C          IF D2 IS ALSO ZERO, THIS ROUTINE SHOULD NOT HAVE BEEN CALLED.
         IF (abs(D2) .lt. vsmall)  GO TO 5001
!C
         RHO = SLOPE/D2
!C          EXTREMUM IS OUTSIDE INTERVAL WHEN RHO .GE. 1/3 .
         IF (RHO .GE. THIRD)  GO TO 5000
         THAT = (TWO*(THREE*RHO-ONE)) / (THREE*(TWO*RHO-ONE))
         PHI = THAT**2 * ((THREE*RHO-ONE)/THREE)
!C
!C          CONVERT TO DISTANCE FROM F2 IF IEXTRM.NE.1 .
         IF (IEXTRM .NE. 1)  PHI = PHI - RHO
!C
!C          TEST FOR EXCEEDING LIMIT, AND ADJUST ACCORDINGLY.
         HPHI = H * ABS(PHI)
         IF (HPHI*ABS(D2) .GT. DFMAX)  THEN
!C           AT THIS POINT, HPHI.GT.0, SO DIVIDE IS OK.
            D2 = SIGN (DFMAX/HPHI, D2)
         ENDIF
      ELSE
!C
         RHO = SLOPE/D1
         LAMBDA = -D2/D1
         IF (dabs(D2).lt.vsmall)  THEN
!C
!C           SPECIAL CASE -- D2.EQ.ZERO .
!C
!C             EXTREMUM IS OUTSIDE INTERVAL WHEN RHO .GE. 1/3 .
            IF (RHO .GE. THIRD)  GO TO 5000
            CP = TWO - THREE*RHO
            NU = ONE - TWO*RHO
            THAT = ONE / (THREE*NU)
         ELSE
            IF (LAMBDA .LE. ZERO)  GO TO 5001
!C
!C           NORMAL CASE -- D1 AND D2 BOTH NONZERO, OPPOSITE SIGNS.
!C
            NU = ONE - LAMBDA - TWO*RHO
            SIGMA = ONE - RHO
            CP = NU + SIGMA
            IF (ABS(NU) .GT. SMALL)  THEN
               RADCAL = (NU - (TWO*RHO+ONE))*NU + SIGMA**2
               IF (RADCAL .LT. ZERO)  GO TO 5002
               THAT = (CP - SQRT(RADCAL)) / (THREE*NU)
            ELSE
               THAT = ONE/(TWO*SIGMA)
            ENDIF
         ENDIF
         PHI = THAT*((NU*THAT - CP)*THAT + ONE)
!C
!C          CONVERT TO DISTANCE FROM F2 IF IEXTRM.NE.1 .
         IF (IEXTRM .NE. 1)  PHI = PHI - RHO
!C
!C          TEST FOR EXCEEDING LIMIT, AND ADJUST ACCORDINGLY.
         HPHI = H * ABS(PHI)
         IF (HPHI*ABS(D1) .GT. DFMAX)  THEN
!C           AT THIS POINT, HPHI.GT.0, SO DIVIDE IS OK.
            D1 = SIGN (DFMAX/HPHI, D1)
            D2 = -LAMBDA*D1
         ENDIF
      ENDIF
!C
!C  NORMAL RETURN.
!C
 5000 CONTINUE
      IERR = 0
      RETURN
!C
!C  ERROR RETURNS.
!C
 5001 CONTINUE
!C     D1 AND D2 BOTH ZERO, OR BOTH NONZERO AND SAME SIGN.
      IERR = -1
!      CALL XERMSG ('SLATEC', 'DPCHSW', 'D1 AND/OR D2 INVALID', IERR, 1)
      RETURN
!C
 5002 CONTINUE
!C     NEGATIVE VALUE OF RADICAL (SHOULD NEVER OCCUR).
      IERR = -2
!      CALL XERMSG ('SLATEC', 'DPCHSW', 'NEGATIVE RADICAL', IERR, 1)
      RETURN
!C------------- LAST LINE OF DPCHSW FOLLOWS -----------------------------
      END