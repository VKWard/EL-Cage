         DOUBLE PRECISION :: POSA,POSB,POSCB1,POSCB2, POSDUMP
         DOUBLE PRECISION :: POSCB3,POSCB4
         DOUBLE PRECISION :: VEC, VECB, OLDVEC, OLDVECB
         DOUBLE PRECISION :: R2ADD, R2ADD_A, R2ADD_B
	 DOUBLE PRECISION :: R2ADDALL, R2ADDALL_A, R2ADDALL_B
         DOUBLE PRECISION :: POSDIST, POSDIST2, POSDISTOLD, POSDIST2B
         INTEGER :: ANV, ANVB, PAIRMIN1, PAIRMIN2

C*******************************************************************
C COMMON groups for  POSITIONS data
C*******************************************************************

        COMMON /POS/ POSA(3*N), POSB(3*N), VEC(3), VECB(3), 
     &               OLDVEC(3), OLDVECB(3), POSDUMP(3*3*N) 
        COMMON /ANV/  ANV(3), ANVB(3)
        COMMON /POSDISTS/ POSDIST(N,N), POSDIST2(N,N), POSDISTOLD(N,N), 
     &         POSDIST2B(N,N) 

        COMMON /POSCB/ POSCB1(3*N), POSCB2(3*N), POSCB3(3*N), POSCB4(3*N)
C*******************************************************************
C POSCB1/2 Positions before/after the previous cagebreak
C*******************************************************************
        COMMON /R2/ R2ADD, R2ADD_A, R2ADD_B, R2ADDALL, R2ADDALL_A, R2ADDALL_B
        COMMON /PAIRMIN/ PAIRMIN1, PAIRMIN2

