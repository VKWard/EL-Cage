C*******************************************************************
C  IO2(step_number) reads in and deals with input data
C*********************************************************************
C
         SUBROUTINE IO2(I3, NRECORD)
         IMPLICIT NONE
         DOUBLE PRECISION DUMMY 
         INTEGER I,J1,J3,J2,J4,I3,I7,I8, K, NRECORD, NDUMMY  
         INCLUDE 'system.h'     ! specifies box size and atom numbers
         INCLUDE 'pos.h'
         INCLUDE 'counter.h'
         LOGICAL DEBUG
         

         DEBUG=.FALSE. 
         
         IF (DEBUG) WRITE(11, *) 'Io2 starts'
         CALL FLUSH(11)
         
C Place results from previous step in old
         IF (I3.GT.1) THEN 
         NEIGHOLD=NEIGH
         NNEIGHOLD=NNEIGHMAT
         EXTRAOLD=NEIGH
         NEXTRAOLD=NNEIGHMAT
         POSB=POSA
         POSDISTOLD=POSDIST

         ELSE
         IF (DEBUG) WRITE(11, *) 'Old vec'
         CALL FLUSH(11)
         OLDVEC=0.0D0
         ENDIF

C Zero counters which need reseting at each step
         
         NATOMSCB=0
         NATOMSCB_A=0
         NATOMSCB_B=0
         NATOMSCBPROD_A=0
         NATOMSCBPROD_B=0
         IF (DEBUG) WRITE(11, *) 'Natomscb=0'
         CALL FLUSH(11)
C Falsify logic tests which need reseting at each step 

         CBREVERSE=.FALSE.  
         CBREVERSE2=.FALSE.

C Cycle through all the atoms, reading the positions data       
         
         IF (DEBUG) WRITE(11, *) 'NFORMAT=', NFORMAT
         CALL FLUSH(11)

C        I7=NRECORD-I3+1
         IF (NFORMAT.eq.1) THEN
         READ(12,REC=I3) (POSA(K), K=1,3*N)
         ELSE IF (NFORMAT.eq.2) THEN
         DO I=1,N        
         READ(12,*) (POSA(K), K=3*I-2,3*I) 
         ENDDO
         ELSE IF (NFORMAT.eq.3) THEN !mb2098: start reading position (3)
         ! ts.data file format

CC Vanessas version:
C         IF (MOD(I3,2).NE.0) THEN ! Read the first minimum of the pair only
C            READ(22,*) DUMMY, DUMMY, NDUMMY, PAIRMIN1, PAIRMIN2 
C            ! This is specific to a three-site rigid body?
C            READ(12,REC=PAIRMIN1) (POSA(K), K=1,3*N)
C         ELSE ! Read the second minimum of the pair this time through
C            READ(12,REC=PAIRMIN2) (POSA(K), K=1,3*N)
C         ENDIF
CC my version:
         IF (MOD(I3,2).NE.0) THEN ! Read the first minimum of the pair only
            READ(22,*) DUMMY, DUMMY, NDUMMY, PAIRMIN1, PAIRMIN2 
            ! This is specific to a three-site rigid body?
            READ(12,REC=PAIRMIN1) (POSA(K), K=1,3*N)
         ELSE ! Read the second minimum of the pair this time through
            READ(12,REC=PAIRMIN2) (POSA(K), K=1,3*N)
              !mb2098: check for box unwrapping:
             DO J1=1,N
                J3=3*(J1-1)
                VEC(1)=POSA(J3+1)-POSB(J3+1)
                VEC(2)=POSA(J3+2)-POSB(J3+2)
                VEC(3)=POSA(J3+3)-POSB(J3+3)
                  ANV(1)=NINT(VEC(1)/BOXLX)
                  ANV(2)=NINT(VEC(2)/BOXLY)
                  ANV(3)=NINT(VEC(3)/BOXLZ)
                POSA(J3+1)=POSA(J3+1)-BOXLX*ANV(1)
                POSA(J3+2)=POSA(J3+2)-BOXLY*ANV(2)
                POSA(J3+3)=POSA(J3+3)-BOXLZ*ANV(3)
             ENDDO
             !mb2098: end check for unwrapping
         ENDIF

         ELSE ! mb2098: finish reading positions (3)
         WRITE (6,*) 'Error: input format undefined or unsupported'
         ENDIF
           
C Calculate vector distances between atoms
         IF (DEBUG) WRITE(11, *) 'Calculate vectors between atoms'
         CALL FLUSH(11)
         DO J1=1,N

            POSDIST(J1,J1)=0.0D0
            POSDIST2(J1,J1)=0.0D0
            J3=3*(J1-1)
            DO J2=J1+1, N
               J4=3*(J2-1)
               VEC(1)=POSA(J3+1)-POSA(J4+1)
               VEC(2)=POSA(J3+2)-POSA(J4+2)
               VEC(3)=POSA(J3+3)-POSA(J4+3)
                  ANV(1)=NINT(VEC(1)/BOXLX)
                  ANV(2)=NINT(VEC(2)/BOXLY)
                  ANV(3)=NINT(VEC(3)/BOXLZ)
            IF (I3.GT.1) THEN
               OLDVEC(1)=POSB(J3+1)-POSB(J4+1)
               OLDVEC(2)=POSB(J3+2)-POSB(J4+2)
               OLDVEC(3)=POSB(J3+3)-POSB(J4+3)
                  ANVB(1)=NINT(OLDVEC(1)/BOXLX)
                  ANVB(2)=NINT(OLDVEC(2)/BOXLY)
                  ANVB(3)=NINT(OLDVEC(3)/BOXLZ)
C These vectors are relative to the current image.
C Needed for cage gain in the current step. 
               OLDVECB(1)=OLDVEC(1)-BOXLX*ANV(1)
               OLDVECB(2)=OLDVEC(2)-BOXLY*ANV(2)
               OLDVECB(3)=OLDVEC(3)-BOXLZ*ANV(3)
C These vectors are relative to the previous image.
C Needed for cage loss in the current step. 
               OLDVEC(1)=OLDVEC(1)-BOXLX*ANVB(1)
               OLDVEC(2)=OLDVEC(2)-BOXLY*ANVB(2)
               OLDVEC(3)=OLDVEC(3)-BOXLZ*ANVB(3)
               VECB(1)=VEC(1)-BOXLX*ANVB(1)
               VECB(2)=VEC(2)-BOXLY*ANVB(2)
               VECB(3)=VEC(3)-BOXLZ*ANVB(3)
             ENDIF
               VEC(1)=VEC(1)-BOXLX*ANV(1)
               VEC(2)=VEC(2)-BOXLY*ANV(2)
               VEC(3)=VEC(3)-BOXLZ*ANV(3)

              POSDIST(J1,J2)= VEC(1)**2 +
     &           VEC(2)**2 +
     &           VEC(3)**2
              POSDIST(J2,J1)=POSDIST(J1,J2) 
            IF (I3.GT.1) THEN
C For cage loss - from previous to current step   
              POSDIST2(J1,J2)= SQRT(VECB(1)**2 + VECB(2)**2 +
     &                 VECB(3)**2)-
     &                 SQRT(OLDVEC(1)**2+OLDVEC(2)**2 +
     &                 OLDVEC(3)**2)
              POSDIST2(J2,J1)=POSDIST2(J1,J2)
C For cage gain - from previous to current step
              POSDIST2B(J1,J2)= SQRT(OLDVECB(1)**2+OLDVECB(2)**2 + 
     &                 OLDVECB(3)**2)-
     &                 SQRT(VEC(1)**2 + VEC(2)**2 +
     &                 VEC(3)**2)
              POSDIST2B(J2,J1)=POSDIST2B(J1,J2)
           ENDIF
            ENDDO
         ENDDO
         IF (DEBUG) WRITE(11, *) 'Finish:Calculate vectors between atoms'
         CALL FLUSH(11)
         RETURN
         END
