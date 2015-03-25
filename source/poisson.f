*DECK CHKPR4
      SUBROUTINE CHKPR4 (IORDER, A, B, M, MBDCND, C, D, N, NBDCND, COFX,
     +   IDMN, IERROR)
C***BEGIN PROLOGUE  CHKPR4
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SEPX4
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (CHKPR4-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This program checks the input parameters for errors.
C
C***SEE ALSO  SEPX4
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  CHKPR4
      EXTERNAL COFX
C***FIRST EXECUTABLE STATEMENT  CHKPR4
      IERROR = 1
      IF (A.GE.B .OR. C.GE.D) RETURN
C
C     CHECK BOUNDARY SWITCHES
C
      IERROR = 2
      IF (MBDCND.LT.0 .OR. MBDCND.GT.4) RETURN
      IERROR = 3
      IF (NBDCND.LT.0 .OR. NBDCND.GT.4) RETURN
C
C     CHECK FIRST DIMENSION IN CALLING ROUTINE
C
      IERROR = 5
      IF (IDMN .LT. 7) RETURN
C
C     CHECK M
C
      IERROR = 6
      IF (M.GT.(IDMN-1) .OR. M.LT.6) RETURN
C
C     CHECK N
C
      IERROR = 7
      IF (N .LT. 5) RETURN
C
C     CHECK IORDER
C
      IERROR = 8
      IF (IORDER.NE.2 .AND. IORDER.NE.4) RETURN
C
C     CHECK THAT EQUATION IS ELLIPTIC
C
      DLX = (B-A)/M
      DO  30 I=2,M
         XI = A+(I-1)*DLX
         CALL COFX (XI,AI,BI,CI)
      IF (AI.GT.0.0) GO TO 10
      IERROR=10
      RETURN
   10 CONTINUE
   30 CONTINUE
C
C     NO ERROR FOUND
C
      IERROR = 0
      RETURN
      END
*DECK CHKSN4
      SUBROUTINE CHKSN4 (MBDCND, NBDCND, ALPHA, BETA, COFX, SINGLR)
C***BEGIN PROLOGUE  CHKSN4
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SEPX4
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (CHKSN4-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This subroutine checks if the PDE SEPX4
C     must solve is a singular operator.
C
C***SEE ALSO  SEPX4
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    SPL4
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  CHKSN4
C
      COMMON /SPL4/   KSWX       ,KSWY       ,K          ,L          ,
     1                AIT        ,BIT        ,CIT        ,DIT        ,
     2                MIT        ,NIT        ,IS         ,MS         ,
     3                JS         ,NS         ,DLX        ,DLY        ,
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4
      LOGICAL         SINGLR
      EXTERNAL COFX
C***FIRST EXECUTABLE STATEMENT  CHKSN4
      SINGLR = .FALSE.
C
C     CHECK IF THE BOUNDARY CONDITIONS ARE
C     ENTIRELY PERIODIC AND/OR MIXED
C
      IF ((MBDCND.NE.0 .AND. MBDCND.NE.3) .OR.
     1    (NBDCND.NE.0 .AND. NBDCND.NE.3)) RETURN
C
C     CHECK THAT MIXED CONDITIONS ARE PURE NEUMAN
C
      IF (MBDCND .NE. 3) GO TO  10
      IF (ALPHA.NE.0.0 .OR. BETA.NE.0.0) RETURN
   10 CONTINUE
C
C     CHECK THAT NON-DERIVATIVE COEFFICIENT FUNCTIONS
C     ARE ZERO
C
      DO  30 I=IS,MS
         XI = AIT+(I-1)*DLX
         CALL COFX (XI,AI,BI,CI)
         IF (CI .NE. 0.0) RETURN
   30 CONTINUE
C
C     THE OPERATOR MUST BE SINGULAR IF THIS POINT IS REACHED
C
      SINGLR = .TRUE.
      RETURN
      END
*DECK COSGEN
      SUBROUTINE COSGEN (N, IJUMP, FNUM, FDEN, A)
C***BEGIN PROLOGUE  COSGEN
C***SUBSIDIARY
C***PURPOSE  Subsidiary to GENBUN
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (COSGEN-S, CMPCSG-C)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This subroutine computes required cosine values in ascending
C     order.  When IJUMP .GT. 1 the routine computes values
C
C        2*COS(J*PI/L) , J=1,2,...,L and J .NE. 0(MOD N/IJUMP+1)
C
C     where L = IJUMP*(N/IJUMP+1).
C
C
C     when IJUMP = 1 it computes
C
C            2*COS((J-FNUM)*PI/(N+FDEN)) ,  J=1, 2, ... ,N
C
C     where
C        FNUM = 0.5, FDEN = 0.0, for regular reduction values.
C        FNUM = 0.0, FDEN = 1.0, for B-R and C-R when ISTAG = 1
C        FNUM = 0.0, FDEN = 0.5, for B-R and C-R when ISTAG = 2
C        FNUM = 0.5, FDEN = 0.5, for B-R and C-R when ISTAG = 2
C                                in POISN2 only.
C
C***SEE ALSO  GENBUN
C***ROUTINES CALLED  PIMACH
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  COSGEN
      DIMENSION       A(*)
C
C
C***FIRST EXECUTABLE STATEMENT  COSGEN
      PI = PIMACH(DUM)
      IF (N .EQ. 0) GO TO 105
      IF (IJUMP .EQ. 1) GO TO 103
      K3 = N/IJUMP+1
      K4 = K3-1
      PIBYN = PI/(N+IJUMP)
      DO 102 K=1,IJUMP
         K1 = (K-1)*K3
         K5 = (K-1)*K4
         DO 101 I=1,K4
            X = K1+I
            K2 = K5+I
            A(K2) = -2.*COS(X*PIBYN)
  101    CONTINUE
  102 CONTINUE
      GO TO 105
  103 CONTINUE
      NP1 = N+1
      Y = PI/(N+FDEN)
      DO 104 I=1,N
         X = NP1-I-FNUM
         A(I) = 2.*COS(X*Y)
  104 CONTINUE
  105 CONTINUE
      RETURN
      END
*DECK DEFE4
      SUBROUTINE DEFE4 (COFX, IDMN, USOL, GRHS)
C***BEGIN PROLOGUE  DEFE4
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SEPX4
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (DEFE4-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This subroutine first approximates the truncation error given by
C     TRUN1(X,Y)=DLX**2*TX+DLY**2*TY where
C     TX=AFUN(X)*UXXXX/12.0+BFUN(X)*UXXX/6.0 on the interior and
C     at the boundaries if periodic (here UXXX,UXXXX are the third
C     and fourth partial derivatives of U with respect to X).
C     TX is of the form AFUN(X)/3.0*(UXXXX/4.0+UXXX/DLX)
C     at X=A or X=B if the boundary condition there is mixed.
C     TX=0.0 along specified boundaries.  TY has symmetric form
C     in Y with X,AFUN(X),BFUN(X) replaced by Y,DFUN(Y),EFUN(Y).
C     The second order solution in USOL is used to approximate
C     (via second order finite differencing) the truncation error
C     and the result is added to the right hand side in GRHS
C     and then transferred to USOL to be used as a new right
C     hand side when calling BLKTRI for a fourth order solution.
C
C***SEE ALSO  SEPX4
C***ROUTINES CALLED  DX4, DY4
C***COMMON BLOCKS    SPL4
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  DEFE4
C
      COMMON /SPL4/   KSWX       ,KSWY       ,K          ,L          ,
     1                AIT        ,BIT        ,CIT        ,DIT        ,
     2                MIT        ,NIT        ,IS         ,MS         ,
     3                JS         ,NS         ,DLX        ,DLY        ,
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4
      DIMENSION       GRHS(IDMN,*)           ,USOL(IDMN,*)
      EXTERNAL COFX
C***FIRST EXECUTABLE STATEMENT  DEFE4
         DO  30 I=IS,MS
            XI = AIT+(I-1)*DLX
            CALL COFX (XI,AI,BI,CI)
         DO 30 J=JS,NS
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT (XI,YJ)
C
            CALL DX4(USOL,IDMN,I,J,UXXX,UXXXX)
            CALL DY4(USOL,IDMN,I,J,UYYY,UYYYY)
            TX = AI*UXXXX/12.0+BI*UXXX/6.0
             TY=UYYYY/12.0
C
C     RESET FORM OF TRUNCATION IF AT BOUNDARY WHICH IS NON-PERIODIC
C
            IF (KSWX.EQ.1 .OR. (I.GT.1 .AND. I.LT.K)) GO TO  10
            TX = AI/3.0*(UXXXX/4.0+UXXX/DLX)
   10       IF (KSWY.EQ.1 .OR. (J.GT.1 .AND. J.LT.L)) GO TO  20
            TY = (UYYYY/4.0+UYYY/DLY)/3.0
   20 GRHS(I,J)=GRHS(I,J)+DLY**2*(DLX**2*TX+DLY**2*TY)
   30    CONTINUE
C
C     RESET THE RIGHT HAND SIDE IN USOL
C
      DO  60 I=IS,MS
         DO  50 J=JS,NS
            USOL(I,J) = GRHS(I,J)
   50    CONTINUE
   60 CONTINUE
      RETURN
      END
*DECK DX4
      SUBROUTINE DX4 (U, IDMN, I, J, UXXX, UXXXX)
C***BEGIN PROLOGUE  DX4
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SEPX4
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (DX4-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This program computes second order finite difference
C     approximations to the third and fourth X
C     partial derivatives of U at the (I,J) mesh point.
C
C***SEE ALSO  SEPX4
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    SPL4
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  DX4
C
      COMMON /SPL4/   KSWX       ,KSWY       ,K          ,L          ,
     1                AIT        ,BIT        ,CIT        ,DIT        ,
     2                MIT        ,NIT        ,IS         ,MS         ,
     3                JS         ,NS         ,DLX        ,DLY        ,
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4
      DIMENSION       U(IDMN,*)
C***FIRST EXECUTABLE STATEMENT  DX4
      IF (I.GT.2 .AND. I.LT.(K-1)) GO TO  50
      IF (I .EQ. 1) GO TO  10
      IF (I .EQ. 2) GO TO  30
      IF (I .EQ. K-1) GO TO  60
      IF (I .EQ. K) GO TO  80
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=A
C
   10 IF (KSWX .EQ. 1) GO TO  20
      UXXX = (-5.0*U(1,J)+18.0*U(2,J)-24.0*U(3,J)+14.0*U(4,J)-
     1                                               3.0*U(5,J))/(TDLX3)
      UXXXX = (3.0*U(1,J)-14.0*U(2,J)+26.0*U(3,J)-24.0*U(4,J)+
     1                                      11.0*U(5,J)-2.0*U(6,J))/DLX4
      RETURN
C
C     PERIODIC AT X=A
C
   20 UXXX = (-U(K-2,J)+2.0*U(K-1,J)-2.0*U(2,J)+U(3,J))/(TDLX3)
      UXXXX = (U(K-2,J)-4.0*U(K-1,J)+6.0*U(1,J)-4.0*U(2,J)+U(3,J))/DLX4
      RETURN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=A+DLX
C
   30 IF (KSWX .EQ. 1) GO TO  40
      UXXX = (-3.0*U(1,J)+10.0*U(2,J)-12.0*U(3,J)+6.0*U(4,J)-U(5,J))/
     1       TDLX3
      UXXXX = (2.0*U(1,J)-9.0*U(2,J)+16.0*U(3,J)-14.0*U(4,J)+6.0*U(5,J)-
     1                                                      U(6,J))/DLX4
      RETURN
C
C     PERIODIC AT X=A+DLX
C
   40 UXXX = (-U(K-1,J)+2.0*U(1,J)-2.0*U(3,J)+U(4,J))/(TDLX3)
      UXXXX = (U(K-1,J)-4.0*U(1,J)+6.0*U(2,J)-4.0*U(3,J)+U(4,J))/DLX4
      RETURN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS ON THE INTERIOR
C
   50 CONTINUE
      UXXX = (-U(I-2,J)+2.0*U(I-1,J)-2.0*U(I+1,J)+U(I+2,J))/TDLX3
      UXXXX = (U(I-2,J)-4.0*U(I-1,J)+6.0*U(I,J)-4.0*U(I+1,J)+U(I+2,J))/
     1        DLX4
      RETURN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=B-DLX
C
   60 IF (KSWX .EQ. 1) GO TO  70
      UXXX = (U(K-4,J)-6.0*U(K-3,J)+12.0*U(K-2,J)-10.0*U(K-1,J)+
     1                                                 3.0*U(K,J))/TDLX3
      UXXXX = (-U(K-5,J)+6.0*U(K-4,J)-14.0*U(K-3,J)+16.0*U(K-2,J)-
     1                                     9.0*U(K-1,J)+2.0*U(K,J))/DLX4
      RETURN
C
C     PERIODIC AT X=B-DLX
C
   70 UXXX = (-U(K-3,J)+2.0*U(K-2,J)-2.0*U(1,J)+U(2,J))/TDLX3
      UXXXX = (U(K-3,J)-4.0*U(K-2,J)+6.0*U(K-1,J)-4.0*U(1,J)+U(2,J))/
     1        DLX4
      RETURN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=B
C
   80 UXXX = -(3.0*U(K-4,J)-14.0*U(K-3,J)+24.0*U(K-2,J)-18.0*U(K-1,J)+
     1                                                 5.0*U(K,J))/TDLX3
      UXXXX = (-2.0*U(K-5,J)+11.0*U(K-4,J)-24.0*U(K-3,J)+26.0*U(K-2,J)-
     1                                    14.0*U(K-1,J)+3.0*U(K,J))/DLX4
      RETURN
      END
*DECK DY4
      SUBROUTINE DY4 (U, IDMN, I, J, UYYY, UYYYY)
C***BEGIN PROLOGUE  DY4
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SEPX4
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (DY4-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This program computes second order finite difference
C     approximations to the third and fourth Y
C     partial derivatives of U at the (I,J) mesh point.
C
C***SEE ALSO  SEPX4
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    SPL4
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  DY4
C
      COMMON /SPL4/   KSWX       ,KSWY       ,K          ,L          ,
     1                AIT        ,BIT        ,CIT        ,DIT        ,
     2                MIT        ,NIT        ,IS         ,MS         ,
     3                JS         ,NS         ,DLX        ,DLY        ,
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4
      DIMENSION       U(IDMN,*)
C***FIRST EXECUTABLE STATEMENT  DY4
      IF (J.GT.2 .AND. J.LT.(L-1)) GO TO  50
      IF (J .EQ. 1) GO TO  10
      IF (J .EQ. 2) GO TO  30
      IF (J .EQ. L-1) GO TO  60
      IF (J .EQ. L) GO TO  80
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=C
C
   10 IF (KSWY .EQ. 1) GO TO  20
      UYYY = (-5.0*U(I,1)+18.0*U(I,2)-24.0*U(I,3)+14.0*U(I,4)-
     1                                                 3.0*U(I,5))/TDLY3
      UYYYY = (3.0*U(I,1)-14.0*U(I,2)+26.0*U(I,3)-24.0*U(I,4)+
     1                                      11.0*U(I,5)-2.0*U(I,6))/DLY4
      RETURN
C
C     PERIODIC AT X=A
C
   20 UYYY = (-U(I,L-2)+2.0*U(I,L-1)-2.0*U(I,2)+U(I,3))/TDLY3
      UYYYY = (U(I,L-2)-4.0*U(I,L-1)+6.0*U(I,1)-4.0*U(I,2)+U(I,3))/DLY4
      RETURN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=C+DLY
C
   30 IF (KSWY .EQ. 1) GO TO  40
      UYYY = (-3.0*U(I,1)+10.0*U(I,2)-12.0*U(I,3)+6.0*U(I,4)-U(I,5))/
     1       TDLY3
      UYYYY = (2.0*U(I,1)-9.0*U(I,2)+16.0*U(I,3)-14.0*U(I,4)+6.0*U(I,5)-
     1                                                      U(I,6))/DLY4
      RETURN
C
C     PERIODIC AT Y=C+DLY
C
   40 UYYY = (-U(I,L-1)+2.0*U(I,1)-2.0*U(I,3)+U(I,4))/TDLY3
      UYYYY = (U(I,L-1)-4.0*U(I,1)+6.0*U(I,2)-4.0*U(I,3)+U(I,4))/DLY4
      RETURN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS ON THE INTERIOR
C
   50 CONTINUE
      UYYY = (-U(I,J-2)+2.0*U(I,J-1)-2.0*U(I,J+1)+U(I,J+2))/TDLY3
      UYYYY = (U(I,J-2)-4.0*U(I,J-1)+6.0*U(I,J)-4.0*U(I,J+1)+U(I,J+2))/
     1        DLY4
      RETURN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=D-DLY
C
   60 IF (KSWY .EQ. 1) GO TO  70
      UYYY = (U(I,L-4)-6.0*U(I,L-3)+12.0*U(I,L-2)-10.0*U(I,L-1)+
     1                                                 3.0*U(I,L))/TDLY3
      UYYYY = (-U(I,L-5)+6.0*U(I,L-4)-14.0*U(I,L-3)+16.0*U(I,L-2)-
     1                                     9.0*U(I,L-1)+2.0*U(I,L))/DLY4
      RETURN
C
C     PERIODIC AT Y=D-DLY
C
   70 CONTINUE
      UYYY = (-U(I,L-3)+2.0*U(I,L-2)-2.0*U(I,1)+U(I,2))/TDLY3
      UYYYY = (U(I,L-3)-4.0*U(I,L-2)+6.0*U(I,L-1)-4.0*U(I,1)+U(I,2))/
     1        DLY4
      RETURN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=D
C
   80 UYYY = -(3.0*U(I,L-4)-14.0*U(I,L-3)+24.0*U(I,L-2)-18.0*U(I,L-1)+
     1                                                 5.0*U(I,L))/TDLY3
      UYYYY = (-2.0*U(I,L-5)+11.0*U(I,L-4)-24.0*U(I,L-3)+26.0*U(I,L-2)-
     1                                    14.0*U(I,L-1)+3.0*U(I,L))/DLY4
      RETURN
      END
*DECK GENBUN
      SUBROUTINE GENBUN (NPEROD, N, MPEROD, M, A, B, C, IDIMY, Y,
     +   IERROR, W)
C***BEGIN PROLOGUE  GENBUN
C***PURPOSE  Solve by a cyclic reduction algorithm the linear system
C            of equations that results from a finite difference
C            approximation to certain 2-d elliptic PDE's on a centered
C            grid .
C***LIBRARY   SLATEC (FISHPACK)
C***CATEGORY  I2B4B
C***TYPE      SINGLE PRECISION (GENBUN-S, CMGNBN-C)
C***KEYWORDS  ELLIPTIC, FISHPACK, PDE, TRIDIAGONAL
C***AUTHOR  Adams, J., (NCAR)
C           Swarztrauber, P. N., (NCAR)
C           Sweet, R., (NCAR)
C***DESCRIPTION
C
C     Subroutine GENBUN solves the linear system of equations
C
C          A(I)*X(I-1,J) + B(I)*X(I,J) + C(I)*X(I+1,J)
C
C          + X(I,J-1) - 2.*X(I,J) + X(I,J+1) = Y(I,J)
C
C               for I = 1,2,...,M  and  J = 1,2,...,N.
C
C     The indices I+1 and I-1 are evaluated modulo M, i.e.,
C     X(0,J) = X(M,J) and X(M+1,J) = X(1,J), and X(I,0) may be equal to
C     0, X(I,2), or X(I,N) and X(I,N+1) may be equal to 0, X(I,N-1), or
C     X(I,1) depending on an input parameter.
C
C
C     * * * * * * * *    Parameter Description     * * * * * * * * * *
C
C             * * * * * *   On Input    * * * * * *
C
C     NPEROD
C       Indicates the values that X(I,0) and X(I,N+1) are assumed to
C       have.
C
C       = 0  If X(I,0) = X(I,N) and X(I,N+1) = X(I,1).
C       = 1  If X(I,0) = X(I,N+1) = 0  .
C       = 2  If X(I,0) = 0 and X(I,N+1) = X(I,N-1).
C       = 3  If X(I,0) = X(I,2) and X(I,N+1) = X(I,N-1).
C       = 4  If X(I,0) = X(I,2) and X(I,N+1) = 0.
C
C     N
C       The number of unknowns in the J-direction.  N must be greater
C       than 2.
C
C     MPEROD
C       = 0 if A(1) and C(M) are not zero.
C       = 1 if A(1) = C(M) = 0.
C
C     M
C       The number of unknowns in the I-direction.  M must be greater
C       than 2.
C
C     A,B,C
C       One-dimensional arrays of length M that specify the
C       coefficients in the linear equations given above.  If MPEROD = 0
C       the array elements must not depend upon the index I, but must be
C       constant.  Specifically, the subroutine checks the following
C       condition
C
C             A(I) = C(1)
C             C(I) = C(1)
C             B(I) = B(1)
C
C       for I=1,2,...,M.
C
C     IDIMY
C       The row (or first) dimension of the two-dimensional array Y as
C       it appears in the program calling GENBUN.  This parameter is
C       used to specify the variable dimension of Y.  IDIMY must be at
C       least M.
C
C     Y
C       A two-dimensional array that specifies the values of the right
C       side of the linear system of equations given above.  Y must be
C       dimensioned at least M*N.
C
C     W
C       A one-dimensional array that must be provided by the user for
C       work space.  W may require up to 4*N + (10 + INT(log2(N)))*M
C       locations.  The actual number of locations used is computed by
C       GENBUN and is returned in location W(1).
C
C
C             * * * * * *   On Output     * * * * * *
C
C     Y
C       Contains the solution X.
C
C     IERROR
C       An error flag that indicates invalid input parameters.  Except
C       for number zero, a solution is not attempted.
C
C       = 0  No error.
C       = 1  M .LE. 2
C       = 2  N .LE. 2
C       = 3  IDIMY .LT. M
C       = 4  NPEROD .LT. 0 or NPEROD .GT. 4
C       = 5  MPEROD .LT. 0 or MPEROD .GT. 1
C       = 6  A(I) .NE. C(1) or C(I) .NE. C(1) or B(I) .NE. B(1) for
C            some I=1,2,...,M.
C       = 7  A(1) .NE. 0 or C(M) .NE. 0 and MPEROD = 1
C
C     W
C       W(1) contains the required length of W.
C
C *Long Description:
C
C     * * * * * * *   Program Specifications    * * * * * * * * * * * *
C
C     Dimension of   A(M),B(M),C(M),Y(IDIMY,N),W(see parameter list)
C     Arguments
C
C     Latest         June 1, 1976
C     Revision
C
C     Subprograms    GENBUN,POISD2,POISN2,POISP2,COSGEN,MERGE,TRIX,TRI3,
C     Required       PIMACH
C
C     Special        NONE
C     Conditions
C
C     Common         NONE
C     Blocks
C
C     I/O            NONE
C
C     Precision      Single
C
C     Specialist     Roland Sweet
C
C     Language       FORTRAN
C
C     History        Standardized April 1, 1973
C                    Revised August 20,1973
C                    Revised January 1, 1976
C
C     Algorithm      The linear system is solved by a cyclic reduction
C                    algorithm described in the reference.
C
C     Space          4944(decimal) = 11520(octal) locations on the NCAR
C     Required       Control Data 7600.
C
C     Timing and        The execution time T on the NCAR Control Data
C     Accuracy       7600 for subroutine GENBUN is roughly proportional
C                    to M*N*log2(N), but also depends on the input
C                    parameter NPEROD.  Some typical values are listed
C                    in the table below.  More comprehensive timing
C                    charts may be found in the reference.
C                       To measure the accuracy of the algorithm a
C                    uniform random number generator was used to create
C                    a solution array X for the system given in the
C                    'PURPOSE' with
C
C                       A(I) = C(I) = -0.5*B(I) = 1,       I=1,2,...,M
C
C                    and, when MPEROD = 1
C
C                       A(1) = C(M) = 0
C                       A(M) = C(1) = 2.
C
C                    The solution X was substituted into the given sys-
C                    tem and, using double precision, a right side Y was
C                    computed.  Using this array Y subroutine GENBUN was
C                    called to produce an approximate solution Z.  Then
C                    the relative error, defined as
C
C                       E = MAX(ABS(Z(I,J)-X(I,J)))/MAX(ABS(X(I,J)))
C
C                    where the two maxima are taken over all I=1,2,...,M
C                    and J=1,2,...,N, was computed.  The value of E is
C                    given in the table below for some typical values of
C                    M and N.
C
C
C                       M (=N)    MPEROD    NPEROD    T(MSECS)    E
C                       ------    ------    ------    --------  ------
C
C                         31        0         0          36     6.E-14
C                         31        1         1          21     4.E-13
C                         31        1         3          41     3.E-13
C                         32        0         0          29     9.E-14
C                         32        1         1          32     3.E-13
C                         32        1         3          48     1.E-13
C                         33        0         0          36     9.E-14
C                         33        1         1          30     4.E-13
C                         33        1         3          34     1.E-13
C                         63        0         0         150     1.E-13
C                         63        1         1          91     1.E-12
C                         63        1         3         173     2.E-13
C                         64        0         0         122     1.E-13
C                         64        1         1         128     1.E-12
C                         64        1         3         199     6.E-13
C                         65        0         0         143     2.E-13
C                         65        1         1         120     1.E-12
C                         65        1         3         138     4.E-13
C
C     Portability    American National Standards Institute Fortran.
C                    The machine dependent constant PI is defined in
C                    function PIMACH.
C
C     Required       COS
C     Resident
C     Routines
C
C     Reference      Sweet, R., 'A Cyclic Reduction Algorithm For
C                    Solving Block Tridiagonal Systems Of Arbitrary
C                    Dimensions,' SIAM J. on Numer. Anal.,
C                    14(Sept., 1977), PP. 706-720.
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C***REFERENCES  R. Sweet, A cyclic reduction algorithm for solving
C                 block tridiagonal systems of arbitrary dimensions,
C                 SIAM Journal on Numerical Analysis 14, (September
C                 1977), pp. 706-720.
C***ROUTINES CALLED  POISD2, POISN2, POISP2
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  GENBUN
C
C
      DIMENSION       Y(IDIMY,*)
      DIMENSION       W(*)       ,B(*)       ,A(*)       ,C(*)
C***FIRST EXECUTABLE STATEMENT  GENBUN
      IERROR = 0
      IF (M .LE. 2) IERROR = 1
      IF (N .LE. 2) IERROR = 2
      IF (IDIMY .LT. M) IERROR = 3
      IF (NPEROD.LT.0 .OR. NPEROD.GT.4) IERROR = 4
      IF (MPEROD.LT.0 .OR. MPEROD.GT.1) IERROR = 5
      IF (MPEROD .EQ. 1) GO TO 102
      DO 101 I=2,M
         IF (A(I) .NE. C(1)) GO TO 103
         IF (C(I) .NE. C(1)) GO TO 103
         IF (B(I) .NE. B(1)) GO TO 103
  101 CONTINUE
      GO TO 104
  102 IF (A(1).NE.0. .OR. C(M).NE.0.) IERROR = 7
      GO TO 104
  103 IERROR = 6
  104 IF (IERROR .NE. 0) RETURN
      MP1 = M+1
      IWBA = MP1
      IWBB = IWBA+M
      IWBC = IWBB+M
      IWB2 = IWBC+M
      IWB3 = IWB2+M
      IWW1 = IWB3+M
      IWW2 = IWW1+M
      IWW3 = IWW2+M
      IWD = IWW3+M
      IWTCOS = IWD+M
      IWP = IWTCOS+4*N
      DO 106 I=1,M
         K = IWBA+I-1
         W(K) = -A(I)
         K = IWBC+I-1
         W(K) = -C(I)
         K = IWBB+I-1
         W(K) = 2.-B(I)
         DO 105 J=1,N
            Y(I,J) = -Y(I,J)
  105    CONTINUE
  106 CONTINUE
      MP = MPEROD+1
      NP = NPEROD+1
      GO TO (114,107),MP
  107 GO TO (108,109,110,111,123),NP
  108 CALL POISP2 (M,N,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWB2),
     1             W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS),
     2             W(IWP))
      GO TO 112
  109 CALL POISD2 (M,N,1,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWW1),
     1             W(IWD),W(IWTCOS),W(IWP))
      GO TO 112
  110 CALL POISN2 (M,N,1,2,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWB2),
     1             W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS),
     2             W(IWP))
      GO TO 112
  111 CALL POISN2 (M,N,1,1,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWB2),
     1             W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS),
     2             W(IWP))
  112 IPSTOR = W(IWW1)
      IREV = 2
      IF (NPEROD .EQ. 4) GO TO 124
  113 GO TO (127,133),MP
  114 CONTINUE
C
C     REORDER UNKNOWNS WHEN MP =0
C
      MH = (M+1)/2
      MHM1 = MH-1
      MODD = 1
      IF (MH*2 .EQ. M) MODD = 2
      DO 119 J=1,N
         DO 115 I=1,MHM1
            MHPI = MH+I
            MHMI = MH-I
            W(I) = Y(MHMI,J)-Y(MHPI,J)
            W(MHPI) = Y(MHMI,J)+Y(MHPI,J)
  115    CONTINUE
         W(MH) = 2.*Y(MH,J)
         GO TO (117,116),MODD
  116    W(M) = 2.*Y(M,J)
  117    CONTINUE
         DO 118 I=1,M
            Y(I,J) = W(I)
  118    CONTINUE
  119 CONTINUE
      K = IWBC+MHM1-1
      I = IWBA+MHM1
      W(K) = 0.
      W(I) = 0.
      W(K+1) = 2.*W(K+1)
      GO TO (120,121),MODD
  120 CONTINUE
      K = IWBB+MHM1-1
      W(K) = W(K)-W(I-1)
      W(IWBC-1) = W(IWBC-1)+W(IWBB-1)
      GO TO 122
  121 W(IWBB-1) = W(K+1)
  122 CONTINUE
      GO TO 107
C
C     REVERSE COLUMNS WHEN NPEROD = 4.
C
  123 IREV = 1
      NBY2 = N/2
  124 DO 126 J=1,NBY2
         MSKIP = N+1-J
         DO 125 I=1,M
            A1 = Y(I,J)
            Y(I,J) = Y(I,MSKIP)
            Y(I,MSKIP) = A1
  125    CONTINUE
  126 CONTINUE
      GO TO (110,113),IREV
  127 CONTINUE
      DO 132 J=1,N
         DO 128 I=1,MHM1
            MHMI = MH-I
            MHPI = MH+I
            W(MHMI) = .5*(Y(MHPI,J)+Y(I,J))
            W(MHPI) = .5*(Y(MHPI,J)-Y(I,J))
  128    CONTINUE
         W(MH) = .5*Y(MH,J)
         GO TO (130,129),MODD
  129    W(M) = .5*Y(M,J)
  130    CONTINUE
         DO 131 I=1,M
            Y(I,J) = W(I)
  131    CONTINUE
  132 CONTINUE
  133 CONTINUE
C
C     RETURN STORAGE REQUIREMENTS FOR W ARRAY.
C
      W(1) = IPSTOR+IWP-1
      RETURN
      END
*DECK MERGE
      SUBROUTINE MERGE(TCOS,I1,M1,I2,M2,I3)
C***BEGIN PROLOGUE  MERGE
C***REFER TO  GENBUN
C
C     This subroutine merges two ascending strings of numbers in the
C     array TCOS.  The first string is of length M1 and starts at
C     TCOS(I1+1).  The second string is of length M2 and starts at
C     TCOS(I2+1).  The merged string goes into TCOS(I3+1).
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  MERGE
      DIMENSION       TCOS(1)
C
C
C***FIRST EXECUTABLE STATEMENT  MERGE
      J1 = 1
      J2 = 1
      J = I3
      IF (M1 .EQ. 0) GO TO 107
      IF (M2 .EQ. 0) GO TO 104
  101 J = J+1
      L = J1+I1
      X = TCOS(L)
      L = J2+I2
      Y = TCOS(L)
      IF (X-Y) 102,102,103
  102 TCOS(J) = X
      J1 = J1+1
      IF (J1 .GT. M1) GO TO 106
      GO TO 101
  103 TCOS(J) = Y
      J2 = J2+1
      IF (J2 .LE. M2) GO TO 101
      IF (J1 .GT. M1) GO TO 109
  104 K = J-J1+1
      DO 105 J=J1,M1
         M = K+J
         L = J+I1
         TCOS(M) = TCOS(L)
  105 CONTINUE
      GO TO 109
  106 CONTINUE
      IF (J2 .GT. M2) GO TO 109
  107 K = J-J2+1
      DO 108 J=J2,M2
         M = K+J
         L = J+I2
         TCOS(M) = TCOS(L)
  108 CONTINUE
  109 CONTINUE
      RETURN
      END
*DECK MINSO4
      SUBROUTINE MINSO4 (USOL, IDMN, ZN, ZM, PERTB)
C***BEGIN PROLOGUE  MINSO4
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SEPX4
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (MINSO4-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This subroutine orthogonalizes the array USOL with respect to
C     the constant array in a weighted least squares norm.
C
C     Entry at MINSO4 occurs when the final solution is
C     to be minimized with respect to the weighted
C     least squares norm.
C
C***SEE ALSO  SEPX4
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    SPL4
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  MINSO4
C
      COMMON /SPL4/   KSWX       ,KSWY       ,K          ,L          ,
     1                AIT        ,BIT        ,CIT        ,DIT        ,
     2                MIT        ,NIT        ,IS         ,MS         ,
     3                JS         ,NS         ,DLX        ,DLY        ,
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4
      DIMENSION       USOL(IDMN,*)           ,ZN(*)      ,ZM(*)
C***FIRST EXECUTABLE STATEMENT  MINSO4
      ISTR = 1
      IFNL = K
      JSTR = 1
      JFNL = L
C
C     COMPUTE WEIGHTED INNER PRODUCTS
C
      UTE = 0.0
      ETE = 0.0
      DO  20 I=IS,MS
         II = I-IS+1
         DO  10 J=JS,NS
            JJ = J-JS+1
            ETE = ETE+ZM(II)*ZN(JJ)
            UTE = UTE+USOL(I,J)*ZM(II)*ZN(JJ)
   10    CONTINUE
   20 CONTINUE
C
C     SET PERTURBATION PARAMETER
C
      PERTRB = UTE/ETE
C
C     SUBTRACT OFF CONSTANT PERTRB
C
      DO  40 I=ISTR,IFNL
         DO  30 J=JSTR,JFNL
            USOL(I,J) = USOL(I,J)-PERTRB
   30    CONTINUE
   40 CONTINUE
      RETURN
      END
*DECK ORTHO4
      SUBROUTINE ORTHO4 (USOL, IDMN, ZN, ZM, PERTRB)
C***BEGIN PROLOGUE  ORTHO4
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SEPX4
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (ORTHO4-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This subroutine orthogonalizes the array USOL with respect to
C     the constant array in a weighted least squares norm.
C
C***SEE ALSO  SEPX4
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    SPL4
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  ORTHO4
C
      COMMON /SPL4/   KSWX       ,KSWY       ,K          ,L          ,
     1                AIT        ,BIT        ,CIT        ,DIT        ,
     2                MIT        ,NIT        ,IS         ,MS         ,
     3                JS         ,NS         ,DLX        ,DLY        ,
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4
      DIMENSION       USOL(IDMN,*)           ,ZN(*)      ,ZM(*)
C***FIRST EXECUTABLE STATEMENT  ORTHO4
      ISTR = IS
      IFNL = MS
      JSTR = JS
      JFNL = NS
C
C     COMPUTE WEIGHTED INNER PRODUCTS
C
      UTE = 0.0
      ETE = 0.0
      DO  20 I=IS,MS
         II = I-IS+1
         DO  10 J=JS,NS
            JJ = J-JS+1
            ETE = ETE+ZM(II)*ZN(JJ)
            UTE = UTE+USOL(I,J)*ZM(II)*ZN(JJ)
   10    CONTINUE
   20 CONTINUE
C
C     SET PERTURBATION PARAMETER
C
      PERTRB = UTE/ETE
C
C     SUBTRACT OFF CONSTANT PERTRB
C
      DO  40 I=ISTR,IFNL
         DO  30 J=JSTR,JFNL
            USOL(I,J) = USOL(I,J)-PERTRB
   30    CONTINUE
   40 CONTINUE
      RETURN
      END
*DECK PIMACH
      FUNCTION PIMACH (DUM)
C***BEGIN PROLOGUE  PIMACH
C***SUBSIDIARY
C***PURPOSE  Subsidiary to HSTCSP, HSTSSP and HWSCSP
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (PIMACH-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This subprogram supplies the value of the constant PI correct to
C     machine precision where
C
C     PI=3.1415926535897932384626433832795028841971693993751058209749446
C
C***SEE ALSO  HSTCSP, HSTSSP, HWSCSP
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  PIMACH
C
C***FIRST EXECUTABLE STATEMENT  PIMACH
      PIMACH = 3.14159265358979
      RETURN
      END
*DECK POISD2
      SUBROUTINE POISD2 (MR, NR, ISTAG, BA, BB, BC, Q, IDIMQ, B, W, D,
     +   TCOS, P)
C***BEGIN PROLOGUE  POISD2
C***SUBSIDIARY
C***PURPOSE  Subsidiary to GENBUN
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (POISD2-S, CMPOSD-C)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Subroutine to solve Poisson's equation for Dirichlet boundary
C     conditions.
C
C     ISTAG = 1 if the last diagonal block is the matrix A.
C     ISTAG = 2 if the last diagonal block is the matrix A+I.
C
C***SEE ALSO  GENBUN
C***ROUTINES CALLED  COSGEN, S1MERG, TRIX
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   920130  Modified to use merge routine S1MERG rather than deleted
C           routine MERGE.  (WRB)
C***END PROLOGUE  POISD2
C
      DIMENSION       Q(IDIMQ,*) ,BA(*)      ,BB(*)      ,BC(*)      ,
     1                TCOS(*)    ,B(*)       ,D(*)       ,W(*)       ,
     2                P(*)
C***FIRST EXECUTABLE STATEMENT  POISD2
      M = MR
      N = NR
      JSH = 0
      FI = 1./ISTAG
      IP = -M
      IPSTOR = 0
      GO TO (101,102),ISTAG
  101 KR = 0
      IRREG = 1
      IF (N .GT. 1) GO TO 106
      TCOS(1) = 0.
      GO TO 103
  102 KR = 1
      JSTSAV = 1
      IRREG = 2
      IF (N .GT. 1) GO TO 106
      TCOS(1) = -1.
  103 DO 104 I=1,M
         B(I) = Q(I,1)
  104 CONTINUE
      CALL TRIX (1,0,M,BA,BB,BC,B,TCOS,D,W)
      DO 105 I=1,M
         Q(I,1) = B(I)
  105 CONTINUE
      GO TO 183
  106 LR = 0
      DO 107 I=1,M
         P(I) = 0.
  107 CONTINUE
      NUN = N
      JST = 1
      JSP = N
C
C     IRREG = 1 WHEN NO IRREGULARITIES HAVE OCCURRED, OTHERWISE IT IS 2.
C
  108 L = 2*JST
      NODD = 2-2*((NUN+1)/2)+NUN
C
C     NODD = 1 WHEN NUN IS ODD, OTHERWISE IT IS 2.
C
      GO TO (110,109),NODD
  109 JSP = JSP-L
      GO TO 111
  110 JSP = JSP-JST
      IF (IRREG .NE. 1) JSP = JSP-L
  111 CONTINUE
C
C     REGULAR REDUCTION
C
      CALL COSGEN (JST,1,0.5,0.0,TCOS)
      IF (L .GT. JSP) GO TO 118
      DO 117 J=L,JSP,L
         JM1 = J-JSH
         JP1 = J+JSH
         JM2 = J-JST
         JP2 = J+JST
         JM3 = JM2-JSH
         JP3 = JP2+JSH
         IF (JST .NE. 1) GO TO 113
         DO 112 I=1,M
            B(I) = 2.*Q(I,J)
            Q(I,J) = Q(I,JM2)+Q(I,JP2)
  112    CONTINUE
         GO TO 115
  113    DO 114 I=1,M
            T = Q(I,J)-Q(I,JM1)-Q(I,JP1)+Q(I,JM2)+Q(I,JP2)
            B(I) = T+Q(I,J)-Q(I,JM3)-Q(I,JP3)
            Q(I,J) = T
  114    CONTINUE
  115    CONTINUE
         CALL TRIX (JST,0,M,BA,BB,BC,B,TCOS,D,W)
         DO 116 I=1,M
            Q(I,J) = Q(I,J)+B(I)
  116    CONTINUE
  117 CONTINUE
C
C     REDUCTION FOR LAST UNKNOWN
C
  118 GO TO (119,136),NODD
  119 GO TO (152,120),IRREG
C
C     ODD NUMBER OF UNKNOWNS
C
  120 JSP = JSP+L
      J = JSP
      JM1 = J-JSH
      JP1 = J+JSH
      JM2 = J-JST
      JP2 = J+JST
      JM3 = JM2-JSH
      GO TO (123,121),ISTAG
  121 CONTINUE
      IF (JST .NE. 1) GO TO 123
      DO 122 I=1,M
         B(I) = Q(I,J)
         Q(I,J) = 0.
  122 CONTINUE
      GO TO 130
  123 GO TO (124,126),NODDPR
  124 DO 125 I=1,M
         IP1 = IP+I
         B(I) = .5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))+P(IP1)+Q(I,J)
  125 CONTINUE
      GO TO 128
  126 DO 127 I=1,M
         B(I) = .5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))+Q(I,JP2)-Q(I,JP1)+Q(I,J)
  127 CONTINUE
  128 DO 129 I=1,M
         Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
  129 CONTINUE
  130 CALL TRIX (JST,0,M,BA,BB,BC,B,TCOS,D,W)
      IP = IP+M
      IPSTOR = MAX(IPSTOR,IP+M)
      DO 131 I=1,M
         IP1 = IP+I
         P(IP1) = Q(I,J)+B(I)
         B(I) = Q(I,JP2)+P(IP1)
  131 CONTINUE
      IF (LR .NE. 0) GO TO 133
      DO 132 I=1,JST
         KRPI = KR+I
         TCOS(KRPI) = TCOS(I)
  132 CONTINUE
      GO TO 134
  133 CONTINUE
      CALL COSGEN (LR,JSTSAV,0.,FI,TCOS(JST+1))
      CALL S1MERG (TCOS,0,JST,JST,LR,KR)
  134 CONTINUE
      CALL COSGEN (KR,JSTSAV,0.0,FI,TCOS)
      CALL TRIX (KR,KR,M,BA,BB,BC,B,TCOS,D,W)
      DO 135 I=1,M
         IP1 = IP+I
         Q(I,J) = Q(I,JM2)+B(I)+P(IP1)
  135 CONTINUE
      LR = KR
      KR = KR+L
      GO TO 152
C
C     EVEN NUMBER OF UNKNOWNS
C
  136 JSP = JSP+L
      J = JSP
      JM1 = J-JSH
      JP1 = J+JSH
      JM2 = J-JST
      JP2 = J+JST
      JM3 = JM2-JSH
      GO TO (137,138),IRREG
  137 CONTINUE
      JSTSAV = JST
      IDEG = JST
      KR = L
      GO TO 139
  138 CALL COSGEN (KR,JSTSAV,0.0,FI,TCOS)
      CALL COSGEN (LR,JSTSAV,0.0,FI,TCOS(KR+1))
      IDEG = KR
      KR = KR+JST
  139 IF (JST .NE. 1) GO TO 141
      IRREG = 2
      DO 140 I=1,M
         B(I) = Q(I,J)
         Q(I,J) = Q(I,JM2)
  140 CONTINUE
      GO TO 150
  141 DO 142 I=1,M
         B(I) = Q(I,J)+.5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
  142 CONTINUE
      GO TO (143,145),IRREG
  143 DO 144 I=1,M
         Q(I,J) = Q(I,JM2)+.5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
  144 CONTINUE
      IRREG = 2
      GO TO 150
  145 CONTINUE
      GO TO (146,148),NODDPR
  146 DO 147 I=1,M
         IP1 = IP+I
         Q(I,J) = Q(I,JM2)+P(IP1)
  147 CONTINUE
      IP = IP-M
      GO TO 150
  148 DO 149 I=1,M
         Q(I,J) = Q(I,JM2)+Q(I,J)-Q(I,JM1)
  149 CONTINUE
  150 CALL TRIX (IDEG,LR,M,BA,BB,BC,B,TCOS,D,W)
      DO 151 I=1,M
         Q(I,J) = Q(I,J)+B(I)
  151 CONTINUE
  152 NUN = NUN/2
      NODDPR = NODD
      JSH = JST
      JST = 2*JST
      IF (NUN .GE. 2) GO TO 108
C
C     START SOLUTION.
C
      J = JSP
      DO 153 I=1,M
         B(I) = Q(I,J)
  153 CONTINUE
      GO TO (154,155),IRREG
  154 CONTINUE
      CALL COSGEN (JST,1,0.5,0.0,TCOS)
      IDEG = JST
      GO TO 156
  155 KR = LR+JST
      CALL COSGEN (KR,JSTSAV,0.0,FI,TCOS)
      CALL COSGEN (LR,JSTSAV,0.0,FI,TCOS(KR+1))
      IDEG = KR
  156 CONTINUE
      CALL TRIX (IDEG,LR,M,BA,BB,BC,B,TCOS,D,W)
      JM1 = J-JSH
      JP1 = J+JSH
      GO TO (157,159),IRREG
  157 DO 158 I=1,M
         Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))+B(I)
  158 CONTINUE
      GO TO 164
  159 GO TO (160,162),NODDPR
  160 DO 161 I=1,M
         IP1 = IP+I
         Q(I,J) = P(IP1)+B(I)
  161 CONTINUE
      IP = IP-M
      GO TO 164
  162 DO 163 I=1,M
         Q(I,J) = Q(I,J)-Q(I,JM1)+B(I)
  163 CONTINUE
  164 CONTINUE
C
C     START BACK SUBSTITUTION.
C
      JST = JST/2
      JSH = JST/2
      NUN = 2*NUN
      IF (NUN .GT. N) GO TO 183
      DO 182 J=JST,N,L
         JM1 = J-JSH
         JP1 = J+JSH
         JM2 = J-JST
         JP2 = J+JST
         IF (J .GT. JST) GO TO 166
         DO 165 I=1,M
            B(I) = Q(I,J)+Q(I,JP2)
  165    CONTINUE
         GO TO 170
  166    IF (JP2 .LE. N) GO TO 168
         DO 167 I=1,M
            B(I) = Q(I,J)+Q(I,JM2)
  167    CONTINUE
         IF (JST .LT. JSTSAV) IRREG = 1
         GO TO (170,171),IRREG
  168    DO 169 I=1,M
            B(I) = Q(I,J)+Q(I,JM2)+Q(I,JP2)
  169    CONTINUE
  170    CONTINUE
         CALL COSGEN (JST,1,0.5,0.0,TCOS)
         IDEG = JST
         JDEG = 0
         GO TO 172
  171    IF (J+L .GT. N) LR = LR-JST
         KR = JST+LR
         CALL COSGEN (KR,JSTSAV,0.0,FI,TCOS)
         CALL COSGEN (LR,JSTSAV,0.0,FI,TCOS(KR+1))
         IDEG = KR
         JDEG = LR
  172    CONTINUE
         CALL TRIX (IDEG,JDEG,M,BA,BB,BC,B,TCOS,D,W)
         IF (JST .GT. 1) GO TO 174
         DO 173 I=1,M
            Q(I,J) = B(I)
  173    CONTINUE
         GO TO 182
  174    IF (JP2 .GT. N) GO TO 177
  175    DO 176 I=1,M
            Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))+B(I)
  176    CONTINUE
         GO TO 182
  177    GO TO (175,178),IRREG
  178    IF (J+JSH .GT. N) GO TO 180
         DO 179 I=1,M
            IP1 = IP+I
            Q(I,J) = B(I)+P(IP1)
  179    CONTINUE
         IP = IP-M
         GO TO 182
  180    DO 181 I=1,M
            Q(I,J) = B(I)+Q(I,J)-Q(I,JM1)
  181    CONTINUE
  182 CONTINUE
      L = L/2
      GO TO 164
  183 CONTINUE
C
C     RETURN STORAGE REQUIREMENTS FOR P VECTORS.
C
      W(1) = IPSTOR
      RETURN
      END
*DECK POISN2
      SUBROUTINE POISN2 (M, N, ISTAG, MIXBND, A, BB, C, Q, IDIMQ, B, B2,
     +   B3, W, W2, W3, D, TCOS, P)
C***BEGIN PROLOGUE  POISN2
C***SUBSIDIARY
C***PURPOSE  Subsidiary to GENBUN
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (POISN2-S, CMPOSN-C)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Subroutine to solve Poisson's equation with Neumann boundary
C     conditions.
C
C     ISTAG = 1 if the last diagonal block is A.
C     ISTAG = 2 if the last diagonal block is A-I.
C     MIXBND = 1 if have Neumann boundary conditions at both boundaries.
C     MIXBND = 2 if have Neumann boundary conditions at bottom and
C     Dirichlet condition at top.  (for this case, must have ISTAG = 1.)
C
C***SEE ALSO  GENBUN
C***ROUTINES CALLED  COSGEN, S1MERG, TRI3, TRIX
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   920130  Modified to use merge routine S1MERG rather than deleted
C           routine MERGE.  (WRB)
C***END PROLOGUE  POISN2
C
      DIMENSION       A(*)       ,BB(*)      ,C(*)       ,Q(IDIMQ,*) ,
     1                B(*)       ,B2(*)      ,B3(*)      ,W(*)       ,
     2                W2(*)      ,W3(*)      ,D(*)       ,TCOS(*)    ,
     3                K(4)       ,P(*)
      EQUIVALENCE     (K(1),K1)  ,(K(2),K2)  ,(K(3),K3)  ,(K(4),K4)
C***FIRST EXECUTABLE STATEMENT  POISN2
      FISTAG = 3-ISTAG
      FNUM = 1./ISTAG
      FDEN = 0.5*(ISTAG-1)
      MR = M
      IP = -MR
      IPSTOR = 0
      I2R = 1
      JR = 2
      NR = N
      NLAST = N
      KR = 1
      LR = 0
      GO TO (101,103),ISTAG
  101 CONTINUE
      DO 102 I=1,MR
         Q(I,N) = .5*Q(I,N)
  102 CONTINUE
      GO TO (103,104),MIXBND
  103 IF (N .LE. 3) GO TO 155
  104 CONTINUE
      JR = 2*I2R
      NROD = 1
      IF ((NR/2)*2 .EQ. NR) NROD = 0
      GO TO (105,106),MIXBND
  105 JSTART = 1
      GO TO 107
  106 JSTART = JR
      NROD = 1-NROD
  107 CONTINUE
      JSTOP = NLAST-JR
      IF (NROD .EQ. 0) JSTOP = JSTOP-I2R
      CALL COSGEN (I2R,1,0.5,0.0,TCOS)
      I2RBY2 = I2R/2
      IF (JSTOP .GE. JSTART) GO TO 108
      J = JR
      GO TO 116
  108 CONTINUE
C
C     REGULAR REDUCTION.
C
      DO 115 J=JSTART,JSTOP,JR
         JP1 = J+I2RBY2
         JP2 = J+I2R
         JP3 = JP2+I2RBY2
         JM1 = J-I2RBY2
         JM2 = J-I2R
         JM3 = JM2-I2RBY2
         IF (J .NE. 1) GO TO 109
         JM1 = JP1
         JM2 = JP2
         JM3 = JP3
  109    CONTINUE
         IF (I2R .NE. 1) GO TO 111
         IF (J .EQ. 1) JM2 = JP2
         DO 110 I=1,MR
            B(I) = 2.*Q(I,J)
            Q(I,J) = Q(I,JM2)+Q(I,JP2)
  110    CONTINUE
         GO TO 113
  111    CONTINUE
         DO 112 I=1,MR
            FI = Q(I,J)
            Q(I,J) = Q(I,J)-Q(I,JM1)-Q(I,JP1)+Q(I,JM2)+Q(I,JP2)
            B(I) = FI+Q(I,J)-Q(I,JM3)-Q(I,JP3)
  112    CONTINUE
  113    CONTINUE
         CALL TRIX (I2R,0,MR,A,BB,C,B,TCOS,D,W)
         DO 114 I=1,MR
            Q(I,J) = Q(I,J)+B(I)
  114    CONTINUE
C
C     END OF REDUCTION FOR REGULAR UNKNOWNS.
C
  115 CONTINUE
C
C     BEGIN SPECIAL REDUCTION FOR LAST UNKNOWN.
C
      J = JSTOP+JR
  116 NLAST = J
      JM1 = J-I2RBY2
      JM2 = J-I2R
      JM3 = JM2-I2RBY2
      IF (NROD .EQ. 0) GO TO 128
C
C     ODD NUMBER OF UNKNOWNS
C
      IF (I2R .NE. 1) GO TO 118
      DO 117 I=1,MR
         B(I) = FISTAG*Q(I,J)
         Q(I,J) = Q(I,JM2)
  117 CONTINUE
      GO TO 126
  118 DO 119 I=1,MR
         B(I) = Q(I,J)+.5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
  119 CONTINUE
      IF (NRODPR .NE. 0) GO TO 121
      DO 120 I=1,MR
         II = IP+I
         Q(I,J) = Q(I,JM2)+P(II)
  120 CONTINUE
      IP = IP-MR
      GO TO 123
  121 CONTINUE
      DO 122 I=1,MR
         Q(I,J) = Q(I,J)-Q(I,JM1)+Q(I,JM2)
  122 CONTINUE
  123 IF (LR .EQ. 0) GO TO 124
      CALL COSGEN (LR,1,0.5,FDEN,TCOS(KR+1))
      GO TO 126
  124 CONTINUE
      DO 125 I=1,MR
         B(I) = FISTAG*B(I)
  125 CONTINUE
  126 CONTINUE
      CALL COSGEN (KR,1,0.5,FDEN,TCOS)
      CALL TRIX (KR,LR,MR,A,BB,C,B,TCOS,D,W)
      DO 127 I=1,MR
         Q(I,J) = Q(I,J)+B(I)
  127 CONTINUE
      KR = KR+I2R
      GO TO 151
  128 CONTINUE
C
C     EVEN NUMBER OF UNKNOWNS
C
      JP1 = J+I2RBY2
      JP2 = J+I2R
      IF (I2R .NE. 1) GO TO 135
      DO 129 I=1,MR
         B(I) = Q(I,J)
  129 CONTINUE
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      IP = 0
      IPSTOR = MR
      GO TO (133,130),ISTAG
  130 DO 131 I=1,MR
         P(I) = B(I)
         B(I) = B(I)+Q(I,N)
  131 CONTINUE
      TCOS(1) = 1.
      TCOS(2) = 0.
      CALL TRIX (1,1,MR,A,BB,C,B,TCOS,D,W)
      DO 132 I=1,MR
         Q(I,J) = Q(I,JM2)+P(I)+B(I)
  132 CONTINUE
      GO TO 150
  133 CONTINUE
      DO 134 I=1,MR
         P(I) = B(I)
         Q(I,J) = Q(I,JM2)+2.*Q(I,JP2)+3.*B(I)
  134 CONTINUE
      GO TO 150
  135 CONTINUE
      DO 136 I=1,MR
         B(I) = Q(I,J)+.5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
  136 CONTINUE
      IF (NRODPR .NE. 0) GO TO 138
      DO 137 I=1,MR
         II = IP+I
         B(I) = B(I)+P(II)
  137 CONTINUE
      GO TO 140
  138 CONTINUE
      DO 139 I=1,MR
         B(I) = B(I)+Q(I,JP2)-Q(I,JP1)
  139 CONTINUE
  140 CONTINUE
      CALL TRIX (I2R,0,MR,A,BB,C,B,TCOS,D,W)
      IP = IP+MR
      IPSTOR = MAX(IPSTOR,IP+MR)
      DO 141 I=1,MR
         II = IP+I
         P(II) = B(I)+.5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
         B(I) = P(II)+Q(I,JP2)
  141 CONTINUE
      IF (LR .EQ. 0) GO TO 142
      CALL COSGEN (LR,1,0.5,FDEN,TCOS(I2R+1))
      CALL S1MERG (TCOS,0,I2R,I2R,LR,KR)
      GO TO 144
  142 DO 143 I=1,I2R
         II = KR+I
         TCOS(II) = TCOS(I)
  143 CONTINUE
  144 CALL COSGEN (KR,1,0.5,FDEN,TCOS)
      IF (LR .NE. 0) GO TO 145
      GO TO (146,145),ISTAG
  145 CONTINUE
      CALL TRIX (KR,KR,MR,A,BB,C,B,TCOS,D,W)
      GO TO 148
  146 CONTINUE
      DO 147 I=1,MR
         B(I) = FISTAG*B(I)
  147 CONTINUE
  148 CONTINUE
      DO 149 I=1,MR
         II = IP+I
         Q(I,J) = Q(I,JM2)+P(II)+B(I)
  149 CONTINUE
  150 CONTINUE
      LR = KR
      KR = KR+JR
  151 CONTINUE
      GO TO (152,153),MIXBND
  152 NR = (NLAST-1)/JR+1
      IF (NR .LE. 3) GO TO 155
      GO TO 154
  153 NR = NLAST/JR
      IF (NR .LE. 1) GO TO 192
  154 I2R = JR
      NRODPR = NROD
      GO TO 104
  155 CONTINUE
C
C      BEGIN SOLUTION
C
      J = 1+JR
      JM1 = J-I2R
      JP1 = J+I2R
      JM2 = NLAST-I2R
      IF (NR .EQ. 2) GO TO 184
      IF (LR .NE. 0) GO TO 170
      IF (N .NE. 3) GO TO 161
C
C     CASE N = 3.
C
      GO TO (156,168),ISTAG
  156 CONTINUE
      DO 157 I=1,MR
         B(I) = Q(I,2)
  157 CONTINUE
      TCOS(1) = 0.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 158 I=1,MR
         Q(I,2) = B(I)
         B(I) = 4.*B(I)+Q(I,1)+2.*Q(I,3)
  158 CONTINUE
      TCOS(1) = -2.
      TCOS(2) = 2.
      I1 = 2
      I2 = 0
      CALL TRIX (I1,I2,MR,A,BB,C,B,TCOS,D,W)
      DO 159 I=1,MR
         Q(I,2) = Q(I,2)+B(I)
         B(I) = Q(I,1)+2.*Q(I,2)
  159 CONTINUE
      TCOS(1) = 0.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 160 I=1,MR
         Q(I,1) = B(I)
  160 CONTINUE
      JR = 1
      I2R = 0
      GO TO 194
C
C     CASE N = 2**P+1
C
  161 CONTINUE
      GO TO (162,170),ISTAG
  162 CONTINUE
      DO 163 I=1,MR
         B(I) = Q(I,J)+.5*Q(I,1)-Q(I,JM1)+Q(I,NLAST)-Q(I,JM2)
  163 CONTINUE
      CALL COSGEN (JR,1,0.5,0.0,TCOS)
      CALL TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
      DO 164 I=1,MR
         Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))+B(I)
         B(I) = Q(I,1)+2.*Q(I,NLAST)+4.*Q(I,J)
  164 CONTINUE
      JR2 = 2*JR
      CALL COSGEN (JR,1,0.0,0.0,TCOS)
      DO 165 I=1,JR
         I1 = JR+I
         I2 = JR+1-I
         TCOS(I1) = -TCOS(I2)
  165 CONTINUE
      CALL TRIX (JR2,0,MR,A,BB,C,B,TCOS,D,W)
      DO 166 I=1,MR
         Q(I,J) = Q(I,J)+B(I)
         B(I) = Q(I,1)+2.*Q(I,J)
  166 CONTINUE
      CALL COSGEN (JR,1,0.5,0.0,TCOS)
      CALL TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
      DO 167 I=1,MR
         Q(I,1) = .5*Q(I,1)-Q(I,JM1)+B(I)
  167 CONTINUE
      GO TO 194
C
C     CASE OF GENERAL N WITH NR = 3 .
C
  168 DO 169 I=1,MR
         B(I) = Q(I,2)
         Q(I,2) = 0.
         B2(I) = Q(I,3)
         B3(I) = Q(I,1)
  169 CONTINUE
      JR = 1
      I2R = 0
      J = 2
      GO TO 177
  170 CONTINUE
      DO 171 I=1,MR
         B(I) = .5*Q(I,1)-Q(I,JM1)+Q(I,J)
  171 CONTINUE
      IF (NROD .NE. 0) GO TO 173
      DO 172 I=1,MR
         II = IP+I
         B(I) = B(I)+P(II)
  172 CONTINUE
      GO TO 175
  173 DO 174 I=1,MR
         B(I) = B(I)+Q(I,NLAST)-Q(I,JM2)
  174 CONTINUE
  175 CONTINUE
      DO 176 I=1,MR
         T = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
         Q(I,J) = T
         B2(I) = Q(I,NLAST)+T
         B3(I) = Q(I,1)+2.*T
  176 CONTINUE
  177 CONTINUE
      K1 = KR+2*JR-1
      K2 = KR+JR
      TCOS(K1+1) = -2.
      K4 = K1+3-ISTAG
      CALL COSGEN (K2+ISTAG-2,1,0.0,FNUM,TCOS(K4))
      K4 = K1+K2+1
      CALL COSGEN (JR-1,1,0.0,1.0,TCOS(K4))
      CALL S1MERG (TCOS,K1,K2,K1+K2,JR-1,0)
      K3 = K1+K2+LR
      CALL COSGEN (JR,1,0.5,0.0,TCOS(K3+1))
      K4 = K3+JR+1
      CALL COSGEN (KR,1,0.5,FDEN,TCOS(K4))
      CALL S1MERG (TCOS,K3,JR,K3+JR,KR,K1)
      IF (LR .EQ. 0) GO TO 178
      CALL COSGEN (LR,1,0.5,FDEN,TCOS(K4))
      CALL S1MERG (TCOS,K3,JR,K3+JR,LR,K3-LR)
      CALL COSGEN (KR,1,0.5,FDEN,TCOS(K4))
  178 K3 = KR
      K4 = KR
      CALL TRI3 (MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
      DO 179 I=1,MR
         B(I) = B(I)+B2(I)+B3(I)
  179 CONTINUE
      TCOS(1) = 2.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 180 I=1,MR
         Q(I,J) = Q(I,J)+B(I)
         B(I) = Q(I,1)+2.*Q(I,J)
  180 CONTINUE
      CALL COSGEN (JR,1,0.5,0.0,TCOS)
      CALL TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
      IF (JR .NE. 1) GO TO 182
      DO 181 I=1,MR
         Q(I,1) = B(I)
  181 CONTINUE
      GO TO 194
  182 CONTINUE
      DO 183 I=1,MR
         Q(I,1) = .5*Q(I,1)-Q(I,JM1)+B(I)
  183 CONTINUE
      GO TO 194
  184 CONTINUE
      IF (N .NE. 2) GO TO 188
C
C     CASE  N = 2
C
      DO 185 I=1,MR
         B(I) = Q(I,1)
  185 CONTINUE
      TCOS(1) = 0.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 186 I=1,MR
         Q(I,1) = B(I)
         B(I) = 2.*(Q(I,2)+B(I))*FISTAG
  186 CONTINUE
      TCOS(1) = -FISTAG
      TCOS(2) = 2.
      CALL TRIX (2,0,MR,A,BB,C,B,TCOS,D,W)
      DO 187 I=1,MR
         Q(I,1) = Q(I,1)+B(I)
  187 CONTINUE
      JR = 1
      I2R = 0
      GO TO 194
  188 CONTINUE
C
C     CASE OF GENERAL N AND NR = 2 .
C
      DO 189 I=1,MR
         II = IP+I
         B3(I) = 0.
         B(I) = Q(I,1)+2.*P(II)
         Q(I,1) = .5*Q(I,1)-Q(I,JM1)
         B2(I) = 2.*(Q(I,1)+Q(I,NLAST))
  189 CONTINUE
      K1 = KR+JR-1
      TCOS(K1+1) = -2.
      K4 = K1+3-ISTAG
      CALL COSGEN (KR+ISTAG-2,1,0.0,FNUM,TCOS(K4))
      K4 = K1+KR+1
      CALL COSGEN (JR-1,1,0.0,1.0,TCOS(K4))
      CALL S1MERG (TCOS,K1,KR,K1+KR,JR-1,0)
      CALL COSGEN (KR,1,0.5,FDEN,TCOS(K1+1))
      K2 = KR
      K4 = K1+K2+1
      CALL COSGEN (LR,1,0.5,FDEN,TCOS(K4))
      K3 = LR
      K4 = 0
      CALL TRI3 (MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
      DO 190 I=1,MR
         B(I) = B(I)+B2(I)
  190 CONTINUE
      TCOS(1) = 2.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 191 I=1,MR
         Q(I,1) = Q(I,1)+B(I)
  191 CONTINUE
      GO TO 194
  192 DO 193 I=1,MR
         B(I) = Q(I,NLAST)
  193 CONTINUE
      GO TO 196
  194 CONTINUE
C
C     START BACK SUBSTITUTION.
C
      J = NLAST-JR
      DO 195 I=1,MR
         B(I) = Q(I,NLAST)+Q(I,J)
  195 CONTINUE
  196 JM2 = NLAST-I2R
      IF (JR .NE. 1) GO TO 198
      DO 197 I=1,MR
         Q(I,NLAST) = 0.
  197 CONTINUE
      GO TO 202
  198 CONTINUE
      IF (NROD .NE. 0) GO TO 200
      DO 199 I=1,MR
         II = IP+I
         Q(I,NLAST) = P(II)
  199 CONTINUE
      IP = IP-MR
      GO TO 202
  200 DO 201 I=1,MR
         Q(I,NLAST) = Q(I,NLAST)-Q(I,JM2)
  201 CONTINUE
  202 CONTINUE
      CALL COSGEN (KR,1,0.5,FDEN,TCOS)
      CALL COSGEN (LR,1,0.5,FDEN,TCOS(KR+1))
      IF (LR .NE. 0) GO TO 204
      DO 203 I=1,MR
         B(I) = FISTAG*B(I)
  203 CONTINUE
  204 CONTINUE
      CALL TRIX (KR,LR,MR,A,BB,C,B,TCOS,D,W)
      DO 205 I=1,MR
         Q(I,NLAST) = Q(I,NLAST)+B(I)
  205 CONTINUE
      NLASTP = NLAST
  206 CONTINUE
      JSTEP = JR
      JR = I2R
      I2R = I2R/2
      IF (JR .EQ. 0) GO TO 222
      GO TO (207,208),MIXBND
  207 JSTART = 1+JR
      GO TO 209
  208 JSTART = JR
  209 CONTINUE
      KR = KR-JR
      IF (NLAST+JR .GT. N) GO TO 210
      KR = KR-JR
      NLAST = NLAST+JR
      JSTOP = NLAST-JSTEP
      GO TO 211
  210 CONTINUE
      JSTOP = NLAST-JR
  211 CONTINUE
      LR = KR-JR
      CALL COSGEN (JR,1,0.5,0.0,TCOS)
      DO 221 J=JSTART,JSTOP,JSTEP
         JM2 = J-JR
         JP2 = J+JR
         IF (J .NE. JR) GO TO 213
         DO 212 I=1,MR
            B(I) = Q(I,J)+Q(I,JP2)
  212    CONTINUE
         GO TO 215
  213    CONTINUE
         DO 214 I=1,MR
            B(I) = Q(I,J)+Q(I,JM2)+Q(I,JP2)
  214    CONTINUE
  215    CONTINUE
         IF (JR .NE. 1) GO TO 217
         DO 216 I=1,MR
            Q(I,J) = 0.
  216    CONTINUE
         GO TO 219
  217    CONTINUE
         JM1 = J-I2R
         JP1 = J+I2R
         DO 218 I=1,MR
            Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
  218    CONTINUE
  219    CONTINUE
         CALL TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
         DO 220 I=1,MR
            Q(I,J) = Q(I,J)+B(I)
  220    CONTINUE
  221 CONTINUE
      NROD = 1
      IF (NLAST+I2R .LE. N) NROD = 0
      IF (NLASTP .NE. NLAST) GO TO 194
      GO TO 206
  222 CONTINUE
C
C     RETURN STORAGE REQUIREMENTS FOR P VECTORS.
C
      W(1) = IPSTOR
      RETURN
      END
*DECK POISP2
      SUBROUTINE POISP2 (M, N, A, BB, C, Q, IDIMQ, B, B2, B3, W, W2, W3,
     +   D, TCOS, P)
C***BEGIN PROLOGUE  POISP2
C***SUBSIDIARY
C***PURPOSE  Subsidiary to GENBUN
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (POISP2-S, CMPOSP-C)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Subroutine to solve Poisson equation with periodic boundary
C     conditions.
C
C***SEE ALSO  GENBUN
C***ROUTINES CALLED  POISD2, POISN2
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  POISP2
C
      DIMENSION       A(*)       ,BB(*)      ,C(*)       ,Q(IDIMQ,*) ,
     1                B(*)       ,B2(*)      ,B3(*)      ,W(*)       ,
     2                W2(*)      ,W3(*)      ,D(*)       ,TCOS(*)    ,
     3                P(*)
C***FIRST EXECUTABLE STATEMENT  POISP2
      MR = M
      NR = (N+1)/2
      NRM1 = NR-1
      IF (2*NR .NE. N) GO TO 107
C
C     EVEN NUMBER OF UNKNOWNS
C
      DO 102 J=1,NRM1
         NRMJ = NR-J
         NRPJ = NR+J
         DO 101 I=1,MR
            S = Q(I,NRMJ)-Q(I,NRPJ)
            T = Q(I,NRMJ)+Q(I,NRPJ)
            Q(I,NRMJ) = S
            Q(I,NRPJ) = T
  101    CONTINUE
  102 CONTINUE
      DO 103 I=1,MR
         Q(I,NR) = 2.*Q(I,NR)
         Q(I,N) = 2.*Q(I,N)
  103 CONTINUE
      CALL POISD2 (MR,NRM1,1,A,BB,C,Q,IDIMQ,B,W,D,TCOS,P)
      IPSTOR = W(1)
      CALL POISN2 (MR,NR+1,1,1,A,BB,C,Q(1,NR),IDIMQ,B,B2,B3,W,W2,W3,D,
     1             TCOS,P)
      IPSTOR = MAX(IPSTOR,INT(W(1)))
      DO 105 J=1,NRM1
         NRMJ = NR-J
         NRPJ = NR+J
         DO 104 I=1,MR
            S = .5*(Q(I,NRPJ)+Q(I,NRMJ))
            T = .5*(Q(I,NRPJ)-Q(I,NRMJ))
            Q(I,NRMJ) = S
            Q(I,NRPJ) = T
  104    CONTINUE
  105 CONTINUE
      DO 106 I=1,MR
         Q(I,NR) = .5*Q(I,NR)
         Q(I,N) = .5*Q(I,N)
  106 CONTINUE
      GO TO 118
  107 CONTINUE
C
C     ODD  NUMBER OF UNKNOWNS
C
      DO 109 J=1,NRM1
         NRPJ = N+1-J
         DO 108 I=1,MR
            S = Q(I,J)-Q(I,NRPJ)
            T = Q(I,J)+Q(I,NRPJ)
            Q(I,J) = S
            Q(I,NRPJ) = T
  108    CONTINUE
  109 CONTINUE
      DO 110 I=1,MR
         Q(I,NR) = 2.*Q(I,NR)
  110 CONTINUE
      LH = NRM1/2
      DO 112 J=1,LH
         NRMJ = NR-J
         DO 111 I=1,MR
            S = Q(I,J)
            Q(I,J) = Q(I,NRMJ)
            Q(I,NRMJ) = S
  111    CONTINUE
  112 CONTINUE
      CALL POISD2 (MR,NRM1,2,A,BB,C,Q,IDIMQ,B,W,D,TCOS,P)
      IPSTOR = W(1)
      CALL POISN2 (MR,NR,2,1,A,BB,C,Q(1,NR),IDIMQ,B,B2,B3,W,W2,W3,D,
     1             TCOS,P)
      IPSTOR = MAX(IPSTOR,INT(W(1)))
      DO 114 J=1,NRM1
         NRPJ = NR+J
         DO 113 I=1,MR
            S = .5*(Q(I,NRPJ)+Q(I,J))
            T = .5*(Q(I,NRPJ)-Q(I,J))
            Q(I,NRPJ) = T
            Q(I,J) = S
  113    CONTINUE
  114 CONTINUE
      DO 115 I=1,MR
         Q(I,NR) = .5*Q(I,NR)
  115 CONTINUE
      DO 117 J=1,LH
         NRMJ = NR-J
         DO 116 I=1,MR
            S = Q(I,J)
            Q(I,J) = Q(I,NRMJ)
            Q(I,NRMJ) = S
  116    CONTINUE
  117 CONTINUE
  118 CONTINUE
C
C     RETURN STORAGE REQUIREMENTS FOR P VECTORS.
C
      W(1) = IPSTOR
      RETURN
      END
*DECK S1MERG
      SUBROUTINE S1MERG (TCOS, I1, M1, I2, M2, I3)
C***BEGIN PROLOGUE  S1MERG
C***SUBSIDIARY
C***PURPOSE  Merge two strings of ascending real numbers.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (S1MERG-S, D1MERG-D, C1MERG-C, I1MERG-I)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   This subroutine merges two ascending strings of numbers in the
C   array TCOS.  The first string is of length M1 and starts at
C   TCOS(I1+1).  The second string is of length M2 and starts at
C   TCOS(I2+1).  The merged string goes into TCOS(I3+1).
C
C***SEE ALSO  GENBUN
C***ROUTINES CALLED  SCOPY
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   901120  Modified to use IF-THEN-ELSE.  Previous spaghetti code did
C           not compile correctly with optimization on the IBM RS6000.
C           (RWC)
C   920130  Code name changed from MERGE to S1MERG.  (WRB)
C***END PROLOGUE  S1MERG
      INTEGER I1, I2, I3, M1, M2
      REAL TCOS(*)
C
      INTEGER J1, J2, J3
C
C***FIRST EXECUTABLE STATEMENT  S1MERG
      IF (M1.EQ.0 .AND. M2.EQ.0) RETURN
C
      IF (M1.EQ.0 .AND. M2.NE.0) THEN
         CALL SCOPY (M2, TCOS(I2+1), 1, TCOS(I3+1), 1)
         RETURN
      ENDIF
C
      IF (M1.NE.0 .AND. M2.EQ.0) THEN
         CALL SCOPY (M1, TCOS(I1+1), 1, TCOS(I3+1), 1)
         RETURN
      ENDIF
C
      J1 = 1
      J2 = 1
      J3 = 1
C
   10 IF (TCOS(I1+J1) .LE. TCOS(I2+J2)) THEN
         TCOS(I3+J3) = TCOS(I1+J1)
         J1 = J1+1
         IF (J1 .GT. M1) THEN
            CALL SCOPY (M2-J2+1, TCOS(I2+J2), 1, TCOS(I3+J3+1), 1)
            RETURN
         ENDIF
      ELSE
         TCOS(I3+J3) = TCOS(I2+J2)
         J2 = J2+1
         IF (J2 .GT. M2) THEN
            CALL SCOPY (M1-J1+1, TCOS(I1+J1), 1, TCOS(I3+J3+1), 1)
            RETURN
         ENDIF
      ENDIF
      J3 = J3+1
      GO TO 10
      END
*DECK SCOPY
      SUBROUTINE SCOPY(N,SX,INCX,SY,INCY)
C***BEGIN PROLOGUE  SCOPY
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D1A5
C***KEYWORDS  BLAS,COPY,LINEAR ALGEBRA,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  Copy s.p. vector y = x
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       SX  single precision vector with N elements
C     INCX  storage spacing between elements of SX
C       SY  single precision vector with N elements
C     INCY  storage spacing between elements of SY
C
C     --Output--
C       SY  copy of vector SX (unchanged if N .LE. 0)
C
C     Copy single precision SX to single precision SY.
C     For I = 0 to N-1, copy  SX(LX+I*INCX) to SY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is
C     defined in a similar way using INCY.
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SCOPY
C
      REAL SX(1),SY(1)
C***FIRST EXECUTABLE STATEMENT  SCOPY
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C        CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7.
C
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SY(I) = SX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        SY(I) = SX(I)
        SY(I + 1) = SX(I + 1)
        SY(I + 2) = SX(I + 2)
        SY(I + 3) = SX(I + 3)
        SY(I + 4) = SX(I + 4)
        SY(I + 5) = SX(I + 5)
        SY(I + 6) = SX(I + 6)
   50 CONTINUE
      RETURN
C
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          SY(I) = SX(I)
   70     CONTINUE
      RETURN
      END
*DECK SEPX4
      SUBROUTINE SEPX4 (IORDER, A, B, M, MBDCND, BDA, ALPHA, BDB, BETA,
     +   C, D, N, NBDCND, BDC, BDD, COFX, GRHS, USOL, IDMN, W, PERTRB,
     +   IERROR)
C***BEGIN PROLOGUE  SEPX4
C***PURPOSE  Solve for either the second or fourth order finite
C            difference approximation to the solution of a separable
C            elliptic partial differential equation on a rectangle.
C            Any combination of periodic or mixed boundary conditions is
C            allowed.
C***LIBRARY   SLATEC (FISHPACK)
C***CATEGORY  I2B1A2
C***TYPE      SINGLE PRECISION (SEPX4-S)
C***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, SEPARABLE
C***AUTHOR  Adams, J., (NCAR)
C           Swarztrauber, P. N., (NCAR)
C           Sweet, R., (NCAR)
C***DESCRIPTION
C
C Purpose                SEPX4 solves for either the second-order
C                        finite difference approximation or a
C                        fourth-order approximation  to the
C                        solution of a separable elliptic equation
C                             AF(X)*UXX+BF(X)*UX+CF(X)*U+UYY = G(X,Y)
C
C                        on a rectangle (X greater than or equal to A
C                        and less than or equal to B; Y greater than
C                        or equal to C and less than or equal to D).
C                        Any combination of periodic or mixed boundary
C                        conditions is allowed.
C                        If boundary conditions in the X direction
C                        are periodic (see MBDCND=0 below) then the
C                        coefficients must satisfy
C                        AF(X)=C1,BF(X)=0,CF(X)=C2 for all X.
C                        Here C1,C2 are constants, C1.GT.0.
C
C                        The possible boundary conditions are
C                        in the X-direction:
C                         (0) Periodic, U(X+B-A,Y)=U(X,Y) for all Y,X
C                         (1) U(A,Y), U(B,Y) are specified for all Y
C                         (2) U(A,Y), dU(B,Y)/dX+BETA*U(B,Y) are
C                             specified for all Y
C                         (3) dU(A,Y)/dX+ALPHA*U(A,Y),dU(B,Y)/dX+
C                             BETA*U(B,Y) are specified for all Y
C                         (4) dU(A,Y)/dX+ALPHA*U(A,Y),U(B,Y) are
C                             specified for all Y
C
C                        In the Y-direction:
C                         (0) Periodic, U(X,Y+D-C)=U(X,Y) for all X,Y
C                         (1) U(X,C),U(X,D) are specified for all X
C                         (2) U(X,C),dU(X,D)/dY are specified for all X
C                         (3) dU(X,C)/DY,dU(X,D)/dY are specified for
C                            all X
C                        (4) dU(X,C)/DY,U(X,D) are specified for all X
C
C Usage                  Call SEPX4(IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,
C                                  BETA,C,D,N,NBDCND,BDC,BDD,COFX,
C                                  GRHS,USOL,IDMN,W,PERTRB,IERROR)
C
C Arguments
C
C                        IORDER
C                          = 2 If a second-order approximation is sought
C                          = 4 If a fourth-order approximation is sought
C
C                        A,B
C                          The range of the X-independent variable;
C                          i.e., X is greater than or equal to A and
C                          less than or equal to B.  A must be less than
C                          B.
C
C                        M
C                          The number of panels into which the interval
C                          [A,B] is subdivided.  Hence, there will be
C                          M+1 grid points in the X-direction given by
C                          XI=A+(I-1)*DLX for I=1,2,...,M+1 where
C                          DLX=(B-A)/M is the panel width.  M must be
C                          less than IDMN and greater than 5.
C
C                        MBDCND
C                          Indicates the type of boundary condition at
C                          X=A and X=B
C                          = 0 If the solution is periodic in X; i.e.,
C                              U(X+B-A,Y)=U(X,Y) for all Y,X
C                          = 1 If the solution is specified at X=A and
C                              X=B; i.e., U(A,Y) and U(B,Y) are
C                              specified for all Y
C                          = 2 If the solution is specified at X=A and
C                              the boundary condition is mixed at X=B;
C                              i.e., U(A,Y) and dU(B,Y)/dX+BETA*U(B,Y)
C                              are specified for all Y
C                          = 3 If the boundary conditions at X=A and X=B
C                              are mixed; i.e., dU(A,Y)/dX+ALPHA*U(A,Y)
C                              and dU(B,Y)/dX+BETA*U(B,Y) are specified
C                              for all Y
C                          = 4 If the boundary condition at X=A is mixed
C                              and the solution is specified at X=B;
C                              i.e., dU(A,Y)/dX+ALPHA*U(A,Y) and U(B,Y)
C                              are specified for all Y
C
C                        BDA
C                          A one-dimensional array of length N+1 that
C                          specifies the values of dU(A,Y)/dX+
C                          ALPHA*U(A,Y) at X=A, when MBDCND=3 or 4.
C                               BDA(J) = dU(A,YJ)/dX+ALPHA*U(A,YJ);
C                               J=1,2,...,N+1
C                          When MBDCND has any other value, BDA is a
C                          dummy parameter.
C
C On Input               ALPHA
C                          The scalar multiplying the solution in case
C                          of a mixed boundary condition AT X=A (see
C                          argument BDA).  If MBDCND = 3,4 then ALPHA is
C                          a dummy parameter.
C
C                        BDB
C                          A one-dimensional array of length N+1 that
C                          specifies the values of dU(B,Y)/dX+
C                          BETA*U(B,Y) at X=B.  when MBDCND=2 or 3
C                               BDB(J) = dU(B,YJ)/dX+BETA*U(B,YJ);
C                               J=1,2,...,N+1
C                          When MBDCND has any other value, BDB is a
C                          dummy parameter.
C
C                        BETA
C                          The scalar multiplying the solution in case
C                          of a mixed boundary condition at X=B (see
C                          argument BDB).  If MBDCND=2,3 then BETA is a
C                          dummy parameter.
C
C                        C,D
C                          The range of the Y-independent variable;
C                          i.e., Y is greater than or equal to C and
C                          less than or equal to D.  C must be less than
C                          D.
C
C                        N
C                          The number of panels into which the interval
C                          [C,D] is subdivided.  Hence, there will be
C                          N+1 grid points in the Y-direction given by
C                          YJ=C+(J-1)*DLY for J=1,2,...,N+1 where
C                          DLY=(D-C)/N is the panel width.  In addition,
C                          N must be greater than 4.
C
C                        NBDCND
C                          Indicates the types of boundary conditions at
C                          Y=C and Y=D
C                          = 0 If the solution is periodic in Y; i.e.,
C                              U(X,Y+D-C)=U(X,Y) for all X,Y
C                          = 1 If the solution is specified at Y=C and
C                              Y = D, i.e., U(X,C) and U(X,D) are
C                              specified for all X
C                          = 2 If the solution is specified at Y=C and
C                              the boundary condition is mixed at Y=D;
C                              i.e., dU(X,C)/dY and U(X,D)
C                              are specified for all X
C                          = 3 If the boundary conditions are mixed at
C                              Y= C and Y=D i.e., dU(X,D)/DY
C                              and dU(X,D)/dY are specified
C                              for all X
C                          = 4 If the boundary condition is mixed at Y=C
C                              and the solution is specified at Y=D;
C                              i.e. dU(X,C)/dY+GAMA*U(X,C) and U(X,D)
C                              are specified for all X
C
C                        BDC
C                          A one-dimensional array of length M+1 that
C                          specifies the value dU(X,C)/DY
C                          at Y=C.  When NBDCND=3 or 4
C                            BDC(I) = dU(XI,C)/DY
C                             I=1,2,...,M+1.
C                          When NBDCND has any other value, BDC is a
C                          dummy parameter.
C
C
C                        BDD
C                          A one-dimensional array of length M+1 that
C                          specifies the value of dU(X,D)/DY
C                          at Y=D.  When NBDCND=2 or 3
C                            BDD(I)=dU(XI,D)/DY
C                             I=1,2,...,M+1.
C                          When NBDCND has any other value, BDD is a
C                          dummy parameter.
C
C
C                        COFX
C                          A user-supplied subprogram with
C                          parameters X, AFUN, BFUN, CFUN which
C                          returns the values of the X-dependent
C                          coefficients AF(X), BF(X), CF(X) in
C                          the elliptic equation at X.
C                          If boundary conditions in the X direction
C                          are periodic then the coefficients
C                          must satisfy AF(X)=C1,BF(X)=0,CF(X)=C2 for
C                          all X.  Here C1.GT.0 and C2 are constants.
C
C                          Note that COFX must be declared external
C                          in the calling routine.
C
C                        GRHS
C                          A two-dimensional array that specifies the
C                          values of the right-hand side of the elliptic
C                          equation; i.e., GRHS(I,J)=G(XI,YI), for
C                          I=2,...,M; J=2,...,N.  At the boundaries,
C                          GRHS is defined by
C
C                          MBDCND   GRHS(1,J)   GRHS(M+1,J)
C                          ------   ---------   -----------
C                            0      G(A,YJ)     G(B,YJ)
C                            1         *           *
C                            2         *        G(B,YJ)  J=1,2,...,N+1
C                            3      G(A,YJ)     G(B,YJ)
C                            4      G(A,YJ)        *
C
C                          NBDCND   GRHS(I,1)   GRHS(I,N+1)
C                          ------   ---------   -----------
C                            0      G(XI,C)     G(XI,D)
C                            1         *           *
C                            2         *        G(XI,D)  I=1,2,...,M+1
C                            3      G(XI,C)     G(XI,D)
C                            4      G(XI,C)        *
C
C                          where * means these quantities are not used.
C                          GRHS should be dimensioned IDMN by at least
C                          N+1 in the calling routine.
C
C                        USOL
C                          A two-dimensional array that specifies the
C                          values of the solution along the boundaries.
C                          At the boundaries, USOL is defined by
C
C                          MBDCND   USOL(1,J)   USOL(M+1,J)
C                          ------   ---------   -----------
C                            0         *           *
C                            1      U(A,YJ)     U(B,YJ)
C                            2      U(A,YJ)        *     J=1,2,...,N+1
C                            3         *           *
C                            4         *        U(B,YJ)
C
C                          NBDCND   USOL(I,1)   USOL(I,N+1)
C                          ------   ---------   -----------
C                            0         *           *
C                            1      U(XI,C)     U(XI,D)
C                            2      U(XI,C)        *     I=1,2,...,M+1
C                            3         *           *
C                            4         *        U(XI,D)
C
C                          where * means the quantities are not used in
C                          the solution.
C
C                          If IORDER=2, the user may equivalence GRHS
C                          and USOL to save space.  Note that in this
C                          case the tables specifying the boundaries of
C                          the GRHS and USOL arrays determine the
C                          boundaries uniquely except at the corners.
C                          If the tables call for both G(X,Y) and
C                          U(X,Y) at a corner then the solution must be
C                          chosen.  For example, if MBDCND=2 and
C                          NBDCND=4, then U(A,C), U(A,D), U(B,D) must be
C                          chosen at the corners in addition to G(B,C).
C
C                          If IORDER=4, then the two arrays, USOL and
C                          GRHS, must be distinct.
C
C                          USOL should be dimensioned IDMN by at least
C                          N+1 in the calling routine.
C
C                        IDMN
C                          The row (or first) dimension of the arrays
C                          GRHS and USOL as it appears in the program
C                          calling SEPX4.  This parameter is used to
C                          specify the variable dimension of GRHS and
C                          USOL.  IDMN must be at least 7 and greater
C                          than or equal to M+1.
C
C                        W
C                          A one-dimensional array that must be provided
C                          by the user for work space.
C                          10*N+(16+INT(log2(N)))*(M+1)+23 will suffice
C                          as a length for W.  The actual length of
C                          W in the calling routine must be set in W(1)
C                          (see IERROR=11).
C
C On Output              USOL
C                          Contains the approximate solution to the
C                          elliptic equation.  USOL(I,J) is the
C                          approximation to U(XI,YJ) for I=1,2...,M+1
C                          and J=1,2,...,N+1.  The approximation has
C                          error O(DLX**2+DLY**2) if called with
C                          IORDER=2 and O(DLX**4+DLY**4) if called with
C                          IORDER=4.
C
C                        W
C                          W(1) contains the exact minimal length (in
C                          floating point) required for the work space
C                          (see IERROR=11).
C
C                        PERTRB
C                          If a combination of periodic or derivative
C                          boundary conditions (i.e., ALPHA=BETA=0 if
C                          MBDCND=3) is specified and if CF(X)=0 for all
C                          X, then a solution to the discretized matrix
C                          equation may not exist (reflecting the non-
C                          uniqueness of solutions to the PDE).  PERTRB
C                          is a constant calculated and subtracted from
C                          the right hand side of the matrix equation
C                          insuring the existence of a solution.
C                          SEPX4 computes this solution which is a
C                          weighted minimal least squares solution to
C                          the original problem.  If singularity is
C                          not detected PERTRB=0.0 is returned by
C                          SEPX4.
C
C                        IERROR
C                          An error flag that indicates invalid input
C                          parameters or failure to find a solution
C                          = 0  No error
C                          = 1  If A greater than B or C greater than D
C                          = 2  If MBDCND less than 0 or MBDCND greater
C                               than 4
C                          = 3  If NBDCND less than 0 or NBDCND greater
C                               than 4
C                          = 4  If attempt to find a solution fails.
C                               (the linear system generated is not
C                               diagonally dominant.)
C                          = 5  If IDMN is too small (see discussion of
C                               IDMN)
C                          = 6  If M is too small or too large (see
C                               discussion of M)
C                          = 7  If N is too small (see discussion of N)
C                          = 8  If IORDER is not 2 or 4
C                          = 10 If AFUN is less than or equal to zero
C                               for some interior mesh point XI
C                          = 11 If the work space length input in W(1)
C                               is less than the exact minimal work
C                               space length required output in W(1).
C                          = 12 If MBDCND=0 and AF(X)=CF(X)=constant
C                               or BF(X)=0 for all X is not true.
C
C *Long Description:
C
C Dimension of           BDA(N+1), BDB(N+1), BDC(M+1), BDD(M+1),
C Arguments              USOL(IDMN,N+1), GRHS(IDMN,N+1),
C                        W (see argument list)
C
C Latest Revision        October 1980
C
C Special Conditions     NONE
C
C Common Blocks          SPL4
C
C I/O                    NONE
C
C Precision              Single
C
C Required Library       NONE
C Files
C
C Specialist             John C. Adams, NCAR, Boulder, Colorado  80307
C
C Language               FORTRAN
C
C
C Entry Points           SEPX4,SPELI4,CHKPR4,CHKSN4,ORTHO4,MINSO4,TRIS4,
C                        DEFE4,DX4,DY4
C
C History                SEPX4 was developed by modifying the ULIB
C                        routine SEPELI during October 1978.
C                        It should be used instead of SEPELI whenever
C                        possible.  The increase in speed is at least
C                        a factor of three.
C
C Algorithm              SEPX4 automatically discretizes the separable
C                        elliptic equation which is then solved by a
C                        generalized cyclic reduction algorithm in the
C                        subroutine POIS.  The fourth order solution
C                        is obtained using the technique of
C                        deferred corrections referenced below.
C
C
C References             Keller, H.B., 'Numerical Methods for Two-point
C                          Boundary-value Problems', Blaisdel (1968),
C                          Waltham, Mass.
C
C                        Swarztrauber, P., and R. Sweet (1975):
C                          'Efficient FORTRAN Subprograms For The
C                          Solution of Elliptic Partial Differential
C                          Equations'.  NCAR Technical Note
C                          NCAR-TN/IA-109, pp. 135-137.
C
C***REFERENCES  H. B. Keller, Numerical Methods for Two-point
C                 Boundary-value Problems, Blaisdel, Waltham, Mass.,
C                 1968.
C               P. N. Swarztrauber and R. Sweet, Efficient Fortran
C                 subprograms for the solution of elliptic equations,
C                 NCAR TN/IA-109, July 1975, 138 pp.
C***ROUTINES CALLED  CHKPR4, SPELI4
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920122  Minor corrections and modifications to prologue.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SEPX4
C
      DIMENSION       GRHS(IDMN,*)           ,USOL(IDMN,*)
      DIMENSION       BDA(*)     ,BDB(*)     ,BDC(*)     ,BDD(*)     ,
     1                W(*)
      EXTERNAL COFX
C***FIRST EXECUTABLE STATEMENT  SEPX4
      CALL CHKPR4(IORDER,A,B,M,MBDCND,C,D,N,NBDCND,COFX,IDMN,IERROR)
      IF (IERROR .NE. 0) RETURN
C
C     COMPUTE MINIMUM WORK SPACE AND CHECK WORK SPACE LENGTH INPUT
C
      L = N+1
      IF (NBDCND .EQ. 0) L = N
      K = M+1
      L = N+1
C     ESTIMATE LOG BASE 2 OF N
      LOG2N=INT(LOG(REAL(N+1))/LOG(2.0)+0.5)
      LENGTH=4*(N+1)+(10+LOG2N)*(M+1)
      IERROR = 11
      LINPUT = INT(W(1)+0.5)
      LOUTPT = LENGTH+6*(K+L)+1
      W(1) = LOUTPT
      IF (LOUTPT .GT. LINPUT) RETURN
      IERROR = 0
C
C     SET WORK SPACE INDICES
C
      I1 = LENGTH+2
      I2 = I1+L
      I3 = I2+L
      I4 = I3+L
      I5 = I4+L
      I6 = I5+L
      I7 = I6+L
      I8 = I7+K
      I9 = I8+K
      I10 = I9+K
      I11 = I10+K
      I12 = I11+K
      I13 = 2
      CALL SPELI4(IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,
     1NBDCND,BDC,BDD,COFX,W(I1),W(I2),W(I3),
     2             W(I4),W(I5),W(I6),W(I7),W(I8),W(I9),W(I10),W(I11),
     3             W(I12),GRHS,USOL,IDMN,W(I13),PERTRB,IERROR)
      RETURN
      END
*DECK SPELI4
      SUBROUTINE SPELI4 (IORDER, A, B, M, MBDCND, BDA, ALPHA, BDB, BETA,
     +   C, D, N, NBDCND, BDC, BDD, COFX, AN, BN, CN, DN, UN, ZN, AM,
     +   BM, CM, DM, UM, ZM, GRHS, USOL, IDMN, W, PERTRB, IERROR)
C***BEGIN PROLOGUE  SPELI4
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SEPX4
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (SPELI4-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     SPELI4 sets up vectors and arrays for input to BLKTRI
C     and computes a second order solution in USOL.  A return jump to
C     SEPX4 occurs if IORDER=2.  If IORDER=4 a fourth order
C     solution is generated in USOL.
C
C***SEE ALSO  SEPX4
C***ROUTINES CALLED  CHKSN4, DEFE4, GENBUN, MINSO4, ORTHO4, TRIS4
C***COMMON BLOCKS    SPL4
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891009  Removed unreferenced variable.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  SPELI4
C
      DIMENSION       BDA(*)     ,BDB(*)     ,BDC(*)     ,BDD(*)     ,
     1                W(*)
      DIMENSION       GRHS(IDMN,*)           ,USOL(IDMN,*)
      DIMENSION       AN(*)      ,BN(*)      ,CN(*)      ,DN(*)      ,
     1                UN(*)      ,ZN(*)
      DIMENSION       AM(*)      ,BM(*)      ,CM(*)      ,DM(*)      ,
     1                UM(*)      ,ZM(*)
      COMMON /SPL4/   KSWX       ,KSWY       ,K          ,L          ,
     1                AIT        ,BIT        ,CIT        ,DIT        ,
     2                MIT        ,NIT        ,IS         ,MS         ,
     3                JS         ,NS         ,DLX        ,DLY        ,
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4
      LOGICAL         SINGLR
      EXTERNAL COFX
C***FIRST EXECUTABLE STATEMENT  SPELI4
      KSWX = MBDCND+1
      KSWY = NBDCND+1
      K = M+1
      L = N+1
      AIT = A
      BIT = B
      CIT = C
      DIT = D
      DLY=(DIT-CIT)/N
C
C     SET RIGHT HAND SIDE VALUES FROM GRHS IN USOL ON THE INTERIOR
C     AND NON-SPECIFIED BOUNDARIES.
C
      DO  20 I=2,M
         DO  10 J=2,N
      USOL(I,J)=DLY**2*GRHS(I,J)
   10    CONTINUE
   20 CONTINUE
      IF (KSWX.EQ.2 .OR. KSWX.EQ.3) GO TO  40
      DO  30 J=2,N
      USOL(1,J)=DLY**2*GRHS(1,J)
   30 CONTINUE
   40 CONTINUE
      IF (KSWX.EQ.2 .OR. KSWX.EQ.5) GO TO  60
      DO  50 J=2,N
      USOL(K,J)=DLY**2*GRHS(K,J)
   50 CONTINUE
   60 CONTINUE
      IF (KSWY.EQ.2 .OR. KSWY.EQ.3) GO TO  80
      DO  70 I=2,M
      USOL(I,1)=DLY**2*GRHS(I,1)
   70 CONTINUE
   80 CONTINUE
      IF (KSWY.EQ.2 .OR. KSWY.EQ.5) GO TO 100
      DO  90 I=2,M
      USOL(I,L)=DLY**2*GRHS(I,L)
   90 CONTINUE
  100 CONTINUE
      IF (KSWX.NE.2 .AND. KSWX.NE.3 .AND. KSWY.NE.2 .AND. KSWY.NE.3)
     1USOL(1,1)=DLY**2*GRHS(1,1)
      IF (KSWX.NE.2 .AND. KSWX.NE.5 .AND. KSWY.NE.2 .AND. KSWY.NE.3)
     1USOL(K,1)=DLY**2*GRHS(K,1)
      IF (KSWX.NE.2 .AND. KSWX.NE.3 .AND. KSWY.NE.2 .AND. KSWY.NE.5)
     1USOL(1,L)=DLY**2*GRHS(1,L)
      IF (KSWX.NE.2 .AND. KSWX.NE.5 .AND. KSWY.NE.2 .AND. KSWY.NE.5)
     1USOL(K,L)=DLY**2*GRHS(K,L)
C
C     SET SWITCHES FOR PERIODIC OR NON-PERIODIC BOUNDARIES
C
      MP=1
      IF(KSWX.EQ.1) MP=0
      NP=NBDCND
C
C     SET DLX,DLY AND SIZE OF BLOCK TRI-DIAGONAL SYSTEM GENERATED
C     IN NINT,MINT
C
      DLX = (BIT-AIT)/M
      MIT = K-1
      IF (KSWX .EQ. 2) MIT = K-2
      IF (KSWX .EQ. 4) MIT = K
      DLY = (DIT-CIT)/N
      NIT = L-1
      IF (KSWY .EQ. 2) NIT = L-2
      IF (KSWY .EQ. 4) NIT = L
      TDLX3 = 2.0*DLX**3
      DLX4 = DLX**4
      TDLY3 = 2.0*DLY**3
      DLY4 = DLY**4
C
C     SET SUBSCRIPT LIMITS FOR PORTION OF ARRAY TO INPUT TO BLKTRI
C
      IS = 1
      JS = 1
      IF (KSWX.EQ.2 .OR. KSWX.EQ.3) IS = 2
      IF (KSWY.EQ.2 .OR. KSWY.EQ.3) JS = 2
      NS = NIT+JS-1
      MS = MIT+IS-1
C
C     SET X - DIRECTION
C
      DO 110 I=1,MIT
         XI = AIT+(IS+I-2)*DLX
         CALL COFX (XI,AI,BI,CI)
         AXI = (AI/DLX-0.5*BI)/DLX
         BXI = -2.*AI/DLX**2+CI
         CXI = (AI/DLX+0.5*BI)/DLX
      AM(I)=DLY**2*AXI
      BM(I)=DLY**2*BXI
      CM(I)=DLY**2*CXI
  110 CONTINUE
C
C     SET Y DIRECTION
C
      DO 120 J=1,NIT
      DYJ=1.0
      EYJ=-2.0
      FYJ=1.0
         AN(J) = DYJ
         BN(J) = EYJ
         CN(J) = FYJ
  120 CONTINUE
C
C     ADJUST EDGES IN X DIRECTION UNLESS PERIODIC
C
      AX1 = AM(1)
      CXM = CM(MIT)
      GO TO (170,130,150,160,140),KSWX
C
C     DIRICHLET-DIRICHLET IN X DIRECTION
C
  130 AM(1) = 0.0
      CM(MIT) = 0.0
      GO TO 170
C
C     MIXED-DIRICHLET IN X DIRECTION
C
  140 AM(1) = 0.0
      BM(1) = BM(1)+2.*ALPHA*DLX*AX1
      CM(1) = CM(1)+AX1
      CM(MIT) = 0.0
      GO TO 170
C
C     DIRICHLET-MIXED IN X DIRECTION
C
  150 AM(1) = 0.0
      AM(MIT) = AM(MIT)+CXM
      BM(MIT) = BM(MIT)-2.*BETA*DLX*CXM
      CM(MIT) = 0.0
      GO TO 170
C
C     MIXED - MIXED IN X DIRECTION
C
  160 CONTINUE
      AM(1) = 0.0
      BM(1) = BM(1)+2.*DLX*ALPHA*AX1
      CM(1) = CM(1)+AX1
      AM(MIT) = AM(MIT)+CXM
      BM(MIT) = BM(MIT)-2.*DLX*BETA*CXM
      CM(MIT) = 0.0
  170 CONTINUE
C
C     ADJUST IN Y DIRECTION UNLESS PERIODIC
C
      DY1 = AN(1)
      FYN = CN(NIT)
      GAMA=0.0
      XNU=0.0
      GO TO (220,180,200,210,190),KSWY
C
C     DIRICHLET-DIRICHLET IN Y DIRECTION
C
  180 CONTINUE
      AN(1) = 0.0
      CN(NIT) = 0.0
      GO TO 220
C
C     MIXED-DIRICHLET IN Y DIRECTION
C
  190 CONTINUE
      AN(1) = 0.0
      BN(1) = BN(1)+2.*DLY*GAMA*DY1
      CN(1) = CN(1)+DY1
      CN(NIT) = 0.0
      GO TO 220
C
C     DIRICHLET-MIXED IN Y DIRECTION
C
  200 AN(1) = 0.0
      AN(NIT) = AN(NIT)+FYN
      BN(NIT) = BN(NIT)-2.*DLY*XNU*FYN
      CN(NIT) = 0.0
      GO TO 220
C
C     MIXED - MIXED DIRECTION IN Y DIRECTION
C
  210 CONTINUE
      AN(1) = 0.0
      BN(1) = BN(1)+2.*DLY*GAMA*DY1
      CN(1) = CN(1)+DY1
      AN(NIT) = AN(NIT)+FYN
      BN(NIT) = BN(NIT)-2.0*DLY*XNU*FYN
      CN(NIT) = 0.0
  220 IF (KSWX .EQ. 1) GO TO 270
C
C     ADJUST USOL ALONG X EDGE
C
      DO 260 J=JS,NS
         IF (KSWX.NE.2 .AND. KSWX.NE.3) GO TO 230
         USOL(IS,J) = USOL(IS,J)-AX1*USOL(1,J)
         GO TO 240
  230    USOL(IS,J) = USOL(IS,J)+2.0*DLX*AX1*BDA(J)
  240    IF (KSWX.NE.2 .AND. KSWX.NE.5) GO TO 250
         USOL(MS,J) = USOL(MS,J)-CXM*USOL(K,J)
         GO TO 260
  250    USOL(MS,J) = USOL(MS,J)-2.0*DLX*CXM*BDB(J)
  260 CONTINUE
  270 IF (KSWY .EQ. 1) GO TO 320
C
C     ADJUST USOL ALONG Y EDGE
C
      DO 310 I=IS,MS
         IF (KSWY.NE.2 .AND. KSWY.NE.3) GO TO 280
         USOL(I,JS) = USOL(I,JS)-DY1*USOL(I,1)
         GO TO 290
  280    USOL(I,JS) = USOL(I,JS)+2.0*DLY*DY1*BDC(I)
  290    IF (KSWY.NE.2 .AND. KSWY.NE.5) GO TO 300
         USOL(I,NS) = USOL(I,NS)-FYN*USOL(I,L)
         GO TO 310
  300    USOL(I,NS) = USOL(I,NS)-2.0*DLY*FYN*BDD(I)
  310 CONTINUE
  320 CONTINUE
C
C     SAVE ADJUSTED EDGES IN GRHS IF IORDER=4
C
      IF (IORDER .NE. 4) GO TO 350
      DO 330 J=JS,NS
         GRHS(IS,J) = USOL(IS,J)
         GRHS(MS,J) = USOL(MS,J)
  330 CONTINUE
      DO 340 I=IS,MS
         GRHS(I,JS) = USOL(I,JS)
         GRHS(I,NS) = USOL(I,NS)
  340 CONTINUE
  350 CONTINUE
      IORD = IORDER
      PERTRB = 0.0
C
C     CHECK IF OPERATOR IS SINGULAR
C
      CALL CHKSN4(MBDCND,NBDCND,ALPHA,BETA,COFX,SINGLR)
C
C     COMPUTE NON-ZERO EIGENVECTOR IN NULL SPACE OF TRANSPOSE
C     IF SINGULAR
C
      IF (SINGLR) CALL TRIS4 (MIT,AM,BM,CM,DM,UM,ZM)
      IF (SINGLR) CALL TRIS4 (NIT,AN,BN,CN,DN,UN,ZN)
C
C     ADJUST RIGHT HAND SIDE IF NECESSARY
C
  360 CONTINUE
      IF (SINGLR) CALL ORTHO4 (USOL,IDMN,ZN,ZM,PERTRB)
C
C     COMPUTE SOLUTION
C
C     SAVE ADJUSTED RIGHT HAND SIDE IN GRHS
      DO 444 J=JS,NS
      DO 444 I=IS,MS
      GRHS(I,J)=USOL(I,J)
  444 CONTINUE

      CALL GENBUN(NP,NIT,MP,MIT,AM,BM,CM,IDMN,USOL(IS,JS),IEROR,W)
C     CHECK IF ERROR DETECTED IN POIS
C     THIS CAN ONLY CORRESPOND TO IERROR=12
      IF(IEROR.EQ.0) GO TO 224
C     SET ERROR FLAG IF IMPROPER COEFFICIENTS INPUT TO POIS
      IERROR=12
      RETURN
  224 CONTINUE
      IF (IERROR .NE. 0) RETURN
C
C     SET PERIODIC BOUNDARIES IF NECESSARY
C
      IF (KSWX .NE. 1) GO TO 380
      DO 370 J=1,L
         USOL(K,J) = USOL(1,J)
  370 CONTINUE
  380 IF (KSWY .NE. 1) GO TO 400
      DO 390 I=1,K
         USOL(I,L) = USOL(I,1)
  390 CONTINUE
  400 CONTINUE
C
C     MINIMIZE SOLUTION WITH RESPECT TO WEIGHTED LEAST SQUARES
C     NORM IF OPERATOR IS SINGULAR
C

      IF (SINGLR) CALL MINSO4 (USOL,IDMN,ZN,ZM,PRTRB)
C
C     RETURN IF DEFERRED CORRECTIONS AND A FOURTH ORDER SOLUTION ARE
C     NOT FLAGGED
C
      IF (IORD .EQ. 2) RETURN
      IORD = 2
C
C     COMPUTE NEW RIGHT HAND SIDE FOR FOURTH ORDER SOLUTION
C
      CALL DEFE4(COFX,IDMN,USOL,GRHS)
      GO TO 360

      END
*DECK TRI3
      SUBROUTINE TRI3 (M, A, B, C, K, Y1, Y2, Y3, TCOS, D, W1, W2, W3)
C***BEGIN PROLOGUE  TRI3
C***SUBSIDIARY
C***PURPOSE  Subsidiary to GENBUN
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (TRI3-S, CMPTR3-C)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Subroutine to solve three linear systems whose common coefficient
C     matrix is a rational function in the matrix given by
C
C                  TRIDIAGONAL (...,A(I),B(I),C(I),...)
C
C***SEE ALSO  GENBUN
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890206  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  TRI3
      DIMENSION       A(*)       ,B(*)       ,C(*)       ,K(4)       ,
     1                TCOS(*)    ,Y1(*)      ,Y2(*)      ,Y3(*)      ,
     2                D(*)       ,W1(*)      ,W2(*)      ,W3(*)
      INTEGER K1P1, K2P1, K3P1, K4P1
C
C***FIRST EXECUTABLE STATEMENT  TRI3
      MM1 = M-1
      K1 = K(1)
      K2 = K(2)
      K3 = K(3)
      K4 = K(4)
      K1P1 = K1+1
      K2P1 = K2+1
      K3P1 = K3+1
      K4P1 = K4+1
      K2K3K4 = K2+K3+K4
      IF (K2K3K4 .EQ. 0) GO TO 101
      L1 = (K1+1)/(K2+1)
      L2 = (K1+1)/(K3+1)
      L3 = (K1+1)/(K4+1)
      LINT1 = 1
      LINT2 = 1
      LINT3 = 1
      KINT1 = K1
      KINT2 = KINT1+K2
      KINT3 = KINT2+K3
  101 CONTINUE
      DO 115 N=1,K1
         X = TCOS(N)
         IF (K2K3K4 .EQ. 0) GO TO 107
         IF (N .NE. L1) GO TO 103
         DO 102 I=1,M
            W1(I) = Y1(I)
  102    CONTINUE
  103    IF (N .NE. L2) GO TO 105
         DO 104 I=1,M
            W2(I) = Y2(I)
  104    CONTINUE
  105    IF (N .NE. L3) GO TO 107
         DO 106 I=1,M
            W3(I) = Y3(I)
  106    CONTINUE
  107    CONTINUE
         Z = 1./(B(1)-X)
         D(1) = C(1)*Z
         Y1(1) = Y1(1)*Z
         Y2(1) = Y2(1)*Z
         Y3(1) = Y3(1)*Z
         DO 108 I=2,M
            Z = 1./(B(I)-X-A(I)*D(I-1))
            D(I) = C(I)*Z
            Y1(I) = (Y1(I)-A(I)*Y1(I-1))*Z
            Y2(I) = (Y2(I)-A(I)*Y2(I-1))*Z
            Y3(I) = (Y3(I)-A(I)*Y3(I-1))*Z
  108    CONTINUE
         DO 109 IP=1,MM1
            I = M-IP
            Y1(I) = Y1(I)-D(I)*Y1(I+1)
            Y2(I) = Y2(I)-D(I)*Y2(I+1)
            Y3(I) = Y3(I)-D(I)*Y3(I+1)
  109    CONTINUE
         IF (K2K3K4 .EQ. 0) GO TO 115
         IF (N .NE. L1) GO TO 111
         I = LINT1+KINT1
         XX = X-TCOS(I)
         DO 110 I=1,M
            Y1(I) = XX*Y1(I)+W1(I)
  110    CONTINUE
         LINT1 = LINT1+1
         L1 = (LINT1*K1P1)/K2P1
  111    IF (N .NE. L2) GO TO 113
         I = LINT2+KINT2
         XX = X-TCOS(I)
         DO 112 I=1,M
            Y2(I) = XX*Y2(I)+W2(I)
  112    CONTINUE
         LINT2 = LINT2+1
         L2 = (LINT2*K1P1)/K3P1
  113    IF (N .NE. L3) GO TO 115
         I = LINT3+KINT3
         XX = X-TCOS(I)
         DO 114 I=1,M
            Y3(I) = XX*Y3(I)+W3(I)
  114    CONTINUE
         LINT3 = LINT3+1
         L3 = (LINT3*K1P1)/K4P1
  115 CONTINUE
      RETURN
      END
*DECK TRIS4
      SUBROUTINE TRIS4 (N, A, B, C, D, U, Z)
C***BEGIN PROLOGUE  TRIS4
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SEPX4
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (TRIS4-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This subroutine solves for a non-zero eigenvector corresponding
C     to the zero eigenvalue of the transpose of the rank
C     deficient ONE matrix with subdiagonal A, diagonal B, and
C     superdiagonal C , with A(1) in the (1,N) position, with
C     C(N) in the (N,1) position, AND all other elements zero.
C
C***SEE ALSO  SEPX4
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  TRIS4
C
      DIMENSION       A(*)       ,B(*)       ,C(*)       ,D(*)       ,
     1                U(*)       ,Z(*)
C***FIRST EXECUTABLE STATEMENT  TRIS4
      BN = B(N)
      D(1) = A(2)/B(1)
      V = A(1)
      U(1) = C(N)/B(1)
      NM2 = N-2
      DO  10 J=2,NM2
         DEN = B(J)-C(J-1)*D(J-1)
         D(J) = A(J+1)/DEN
         U(J) = -C(J-1)*U(J-1)/DEN
         BN = BN-V*U(J-1)
         V = -V*D(J-1)
   10 CONTINUE
      DEN = B(N-1)-C(N-2)*D(N-2)
      D(N-1) = (A(N)-C(N-2)*U(N-2))/DEN
      AN = C(N-1)-V*D(N-2)
      BN = BN-V*U(N-2)
      DEN = BN-AN*D(N-1)
C
C     SET LAST COMPONENT EQUAL TO ONE
C
      Z(N) = 1.0
      Z(N-1) = -D(N-1)
      NM1 = N-1
      DO  20 J=2,NM1
         K = N-J
         Z(K) = -D(K)*Z(K+1)-U(K)*Z(N)
   20 CONTINUE
      RETURN
      END
*DECK TRIX
      SUBROUTINE TRIX (IDEGBR, IDEGCR, M, A, B, C, Y, TCOS, D, W)
C***BEGIN PROLOGUE  TRIX
C***SUBSIDIARY
C***PURPOSE  Subsidiary to GENBUN
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (TRIX-S, CMPTRX-C)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Subroutine to solve a system of linear equations where the
C     coefficient matrix is a rational function in the matrix given by
C     TRIDIAGONAL  ( . . . , A(I), B(I), C(I), . . . ).
C
C***SEE ALSO  GENBUN
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  TRIX
C
      DIMENSION       A(*)       ,B(*)       ,C(*)       ,Y(*)       ,
     1                TCOS(*)    ,D(*)       ,W(*)
      INTEGER KB, KC
C***FIRST EXECUTABLE STATEMENT  TRIX
      MM1 = M-1
      KB = IDEGBR+1
      KC = IDEGCR+1
      L = (IDEGBR+1)/(IDEGCR+1)
      LINT = 1
      DO 108 K=1,IDEGBR
         X = TCOS(K)
         IF (K .NE. L) GO TO 102
         I = IDEGBR+LINT
         XX = X-TCOS(I)
         DO 101 I=1,M
            W(I) = Y(I)
            Y(I) = XX*Y(I)
  101    CONTINUE
  102    CONTINUE
         Z = 1./(B(1)-X)
         D(1) = C(1)*Z
         Y(1) = Y(1)*Z
         DO 103 I=2,MM1
            Z = 1./(B(I)-X-A(I)*D(I-1))
            D(I) = C(I)*Z
            Y(I) = (Y(I)-A(I)*Y(I-1))*Z
  103    CONTINUE
         Z = B(M)-X-A(M)*D(MM1)
         IF (Z .NE. 0.) GO TO 104
         Y(M) = 0.
         GO TO 105
  104    Y(M) = (Y(M)-A(M)*Y(MM1))/Z
  105    CONTINUE
         DO 106 IP=1,MM1
            I = M-IP
            Y(I) = Y(I)-D(I)*Y(I+1)
  106    CONTINUE
         IF (K .NE. L) GO TO 108
         DO 107 I=1,M
            Y(I) = Y(I)+W(I)
  107    CONTINUE
         LINT = LINT+1
         L = (LINT*KB)/KC
  108 CONTINUE
      RETURN
      END
