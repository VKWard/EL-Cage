C*******************************************************************
C  CAGEBREAK(STEP_NUMBER,ATOM_NUMBER) checks for cage breaking
C*********************************************************************
C
         SUBROUTINE CAGEBREAK(I2,I)
         IMPLICIT NONE
         INCLUDE 'system.h'
         INCLUDE 'pos.h'
         INCLUDE 'timing.h'
         DOUBLE PRECISION  POSDIST3, DISTV, TIMES, MINDUR
         DOUBLE PRECISION  R2TOTAL, R2STEP1, R2STEP2
         DOUBLE PRECISION  XDIFF, YDIFF, ZDIFF, CANGLE 
         DOUBLE PRECISION  XDIFF2, YDIFF2, ZDIFF2 
         INTEGER :: CB(N), CBIN(N), CBOUT(N)
         INTEGER :: NDIVCB(N,10000),NDIVCBC(N)
         INTEGER NCAGEB,I, I2, J2, J3, J7, J8,J9, NCBOUT, NCBIN
         INTEGER NTEMP, NTEMP2, NTEMP3, NTEMP5, NDIV, NDIVCH
         LOGICAL REVING, REVING2, CAGE
         LOGICAL DEBUG
         integer :: NCUT
         INCLUDE 'counter.h'    
         INCLUDE 'criteria.h'
         INCLUDE 'switch.h' 
    

         DEBUG=.FALSE.
         NCAGEB=1         ! Number of neighbours lost/gained + 1 (each atom in each step)
         NCBIN=1          ! Number of neighbours gained + 1 (each atom in each step)
         NCBOUT=1         ! Number of neighbours lost + 1 (each atom in each step)
         CB(:)=0          ! Number list of neighbours lost 
         REVING=.FALSE.   ! If reversal type 1 takes place (each atom in each step)
         REVING2=.FALSE.  ! If reversal type 2 takes place (each atom in each step)
         CAGE=.FALSE.
        
        if(I.le.NA) then  ! change NCUT for atomtype A/B. 
          NCUT=NCUTA
        else
          NCUT=NCUTB
        endif
 
        IF (.NOT.ABSDIST) THEN
C We can check neighbours lost and new neighbours.  
          IF (CHANGEOUT) THEN

          IF (DEBUG) WRITE(11, *) 'Check all previous neighbours, atom=',I
            DO J2=1, NNEIGHOLD(I)
              IF (DEBUG) WRITE(11, *) 'Check neighbours, neighbour=', J2 
              DO J3=1,NEXTRAMAT(I)  
C Check with extended list and if the neighbour is still present, set cage to true.
        !      IF (DEBUG) WRITE(11, *) 'IF statement'
                IF (NEIGHOLD(I, J2).EQ.EXTRA(I, J3)) CAGE=.TRUE.   
        !      IF (DEBUG) WRITE(11, *) 'end of IF statement'
              ENDDO
              IF (.NOT.CAGE) THEN
                CBOUT(NCBOUT)=NEIGHOLD(I, J2)           
                NCBOUT=NCBOUT+1
              ELSE
              ENDIF
              CAGE=.FALSE.
            ENDDO
          ENDIF
          IF (CHANGEIN) THEN
            IF (DEBUG) WRITE(11, *) 'Check all new neighbours, atom=',I
            DO J2=1, NNEIGHMAT(I)
              IF (DEBUG) WRITE(11, *) 'Check neighbours, neighbour=', J2
              DO J3=1,NEXTRAOLD(I)
C Check with extended list and if the neighbour is still present, set cage to true.
!              IF (DEBUG) WRITE(11, *) 'IF statement'
                IF (NEIGH(I, J2).EQ.EXTRAOLD(I, J3)) CAGE=.TRUE.
!              IF (DEBUG) WRITE(11, *) 'end of IF statement'
              ENDDO
              IF (.NOT.CAGE) THEN
                CBIN(NCBIN)=NEIGH(I, J2)
                NCBIN=NCBIN+1
              ELSE
              ENDIF
              CAGE=.FALSE.
            ENDDO
          ENDIF   
        
          IF (CHANGEOR) THEN  
            IF (NCBOUT.GE.NCBIN) THEN
              NCAGEB=NCBOUT
              CB(:)=-CBOUT(:)
            ELSE
              NCAGEB=NCBIN
              CB(:)=CBIN(:) 
            ENDIF   
          ELSE
            NCAGEB=NCBOUT+NCBIN-1
            CB(:)=-CBOUT(:)
            DO J9=1,NCBIN-1
              CB(NCBOUT+J9-1)=CBIN(J9) 
            ENDDO

          ENDIF

       ELSE  ! i.e. IF(ABSDIST)
C Absolute distance checking rather than neighbour analysis
         DISTV=DISTVAA
           POSDIST3=(POSB(3*I-2)-POSA(3*I-2))**2 +
     &           (POSB(3*I-1)-POSA(3*I-1))**2 +
     &           (POSB(3*I)-POSA(3*I))**2
         IF (SQRT(POSDIST3).GE.DISTV) NCAGEB=NCUT+1
C        write(*, *) 'WARNING: absolute distance checking has not been
C    &                implemented for separate NCUT-paramters!'
       ENDIF  !mb2098: end of if statement to check for cagebreaks

         IF (NCAGEB.GT.NCUT) THEN
C Check for reverse move of either type for this atom.
C Total for this and previous cagebreak
100        R2TOTAL= (POSCB1(3*I-2)-POSA(3*I-2))**2
     &             +(POSCB1(3*I-1)-POSA(3*I-1))**2
     &             +(POSCB1(3*I)-POSA(3*I))**2
C Current cagebreak
           R2STEP2= (POSB(3*I-2)-POSA(3*I-2))**2
     &             +(POSB(3*I-1)-POSA(3*I-1))**2
     &             +(POSB(3*I)-POSA(3*I))**2
C Previous cagebreak
           R2STEP1= (POSCB1(3*I-2)-POSCB2(3*I-2))**2
     &             +(POSCB1(3*I-1)-POSCB2(3*I-1))**2
     &             +(POSCB1(3*I)-POSCB2(3*I))**2
           IF(R2STEP1.gt.0) THEN 
C          write(*,*) 'identified cb for atom',I,'NCUT=',NCUT
C          write(*,*) "Old old pos", POSCB1(3*I-2:3*I)
C          write(*,*) "Old new pos", POSCB2(3*I-2:3*I)
           write(19,*) I, I2, POSB(3*I-2:3*I)
           write(19,*) I, I2, POSA(3*I-2:3*I)
C          write(*,*) "R2Total, R2step2, R2step1", I, R2TOTAL, R2STEP2, R2STEP1
           ELSE
C          write(*,*) 'identified cb for atom',I,'NCUT=',NCUT
           write(19,*) I, I2, POSB(3*I-2:3*I)
           write(19,*) I, I2, POSA(3*I-2:3*I)
C          write(*,*) 'First CB for this atom, assumed productive.'
           ENDIF
C Criterion for identical steps
C If the step is identical, the reverse must have taken place as a 
C non-cage-breaking step.
            IF(ABS(R2TOTAL-R2STEP1).LT.SAMECUT.OR.ABS(R2TOTAL-R2STEP2).LT.SAMECUT) THEN
C           write(*,*) "NCB reversal detected"
            WRITE(15, 77) I2, I, R2TOTAL, R2STEP1, R2STEP2
            CBREVERSE2=.TRUE.
            NREV2ALL=NREV2ALL+1
            NCOUNTREV(I)=NCOUNTREV(I)+2
C           IF (NCOUNTREV(I).EQ.2) THEN
C            IF(I.le.NA) THEN
C               R2ADD_A=R2ADD_A-R2STEP1
C                NPROD_A=NPROD_A-1
C            ELSE
C               R2ADD_B=R2ADD_B-R2STEP1
C                NPROD_B=NPROD_B-1
C            ENDIF    
C           ENDIF        
            REVING2=.TRUE.
            ELSE 
                IF(R2TOTAL.LT.REVCUT) THEN 
C Criterion for reverse steps
C            write(*,*) "Direct reversal detected"
             WRITE(14, 77) I2, I,R2TOTAL, R2STEP1, R2STEP2
             CBREVERSE=.TRUE.
             NREV1ALL=NREV1ALL+1
             NCOUNTREV(I)=NCOUNTREV(I)+1
C           IF (NCOUNTREV(I).EQ.1) THEN
C            IF(I.le.NA) THEN
C               R2ADD_A=R2ADD_A-R2STEP1
C                NPROD_A=NPROD_A-1
C            ELSE
C               R2ADD_B=R2ADD_B-R2STEP1
C                NPROD_B=NPROD_B-1
C            ENDIF    
C           ENDIF        
             REVING=.TRUE.           
             ELSE
C            write(*,*) "No reversal detected"
C Irreversible steps
             REVING=.FALSE.
             REVING2=.FALSE.
             IF (MOD(NCOUNTREV(I),2).EQ.0) THEN
             R2ADD=R2STEP2+R2ADD
C             write(18, 99) I2, I, R2STEP1
             IF(I.le.NA) THEN
C               write(*,*) "Adding ", R2STEP2, "to R2ADD_A for atom ", I, "at step ", I2
                R2ADD_A=R2ADD_A+R2STEP2
                 NPROD_A=NPROD_A+1
C                PRINT*, 'R2ADD_A= ', R2ADD_A 
             ELSE
C               write(*,*) "Adding ", R2STEP2, "to R2ADD_B for atom ",I, "at step ", I2
                R2ADD_B=R2ADD_B+R2STEP2
                 NPROD_B=NPROD_B+1
C                PRINT*, 'R2ADD_B= ', R2ADD_B 
             ENDIF
             ENDIF
             NCOUNTREV(I)=0
             ENDIF
            ENDIF
            R2ADDALL=R2ADDALL+R2STEP2
            WRITE(20,*) SQRT(R2STEP2), I
           
            IF(I.le.NA) THEN
C               write(*,*) "Adding ", R2STEP2, "to R2ADDALL_A for atom ",I, "at step ", I2
                R2ADDALL_A=R2ADDALL_A+R2STEP2
                write(18, 99) I2, I, R2STEP2, 0.0, POSB(3*I-2)-POSA(3*I-2), 
     &                             POSB(3*I-1)-POSA(3*I-1), POSB(3*I)-POSA(3*I)
            ELSE
C               write(*,*) "Adding ", R2STEP2, "to R2ADDALL_B for atom",I, "at step ", I2
                R2ADDALL_B=R2ADDALL_B+R2STEP2
                write(18, 99) I2, I, 0.0, R2STEP2, POSB(3*I-2)-POSA(3*I-2), 
     &                              POSB(3*I-1)-POSA(3*I-1), POSB(3*I)-POSA(3*I)
            ENDIF
99           format(I8, I8, 5F15.10)
77           FORMAT(I7, I6,3F15.10)

C Record positions before and after the cage-breaking event
            POSCB3(3*I-2)=POSCB1(3*I-2)
            POSCB3(3*I-1)=POSCB1(3*I-1)
            POSCB3(3*I)=POSCB1(3*I)
            POSCB4(3*I-2)=POSCB2(3*I-2)
            POSCB4(3*I-1)=POSCB2(3*I-1)
            POSCB4(3*I)=POSCB2(3*I)
            POSCB1(3*I-2)=POSB(3*I-2)
            POSCB1(3*I-1)=POSB(3*I-1)
            POSCB1(3*I)=POSB(3*I)
            POSCB2(3*I-2)=POSA(3*I-2)
            POSCB2(3*I-1)=POSA(3*I-1)
            POSCB2(3*I)=POSA(3*I)

         NATOMSCB=NATOMSCB+1
         if(I.le.NA) then
             NATOMSCB_A=NATOMSCB_A+1
             if((.not.REVING).and.(.not.REVING))NATOMSCBPROD_A=NATOMSCBPROD_A+1
C            PRINT*, 'cagebreak', ' NATOMSCBPROD_A= ', NATOMSCBPROD_A
         else 
             NATOMSCB_B=NATOMSCB_B+1
             if((.not.REVING).and.(.not.REVING))NATOMSCBPROD_B=NATOMSCBPROD_B+1
C            PRINT*, 'cagebreak', ' NATOMSCBPROD_B= ',   NATOMSCBPROD_B
         end if
         WRITE (11, 42) I2, I, NCOUNTREV(I),NCAGEB-1,
     &                   CB(1:NCAGEB-1) 
         IF (NFORMAT.EQ.3) WRITE(6,*) PAIRMIN1, PAIRMIN2
C        WRITE (11, 42) I2, I, TIME(I),NCOUNTREV(I),NCAGEB-1, 
    
42       FORMAT(I8,I8,I5,20I8) 


         IF(TIMEDIV) THEN

C NDIVCB(I,NDIVCBC(I)) should be the time interval of the last cagebreak for this atom
         NTEMP3=NDIVCBC(I)
         NTEMP=NDIV-NDIVCB(I, NTEMP3)
CIF3 Was the previous cagebreaking event in a different time interval 
         IF (REVING.AND.NTEMP.EQ.0) REM=REM+1
         IF (REVING.AND.NTEMP.GT.0) THEN 
C Cycle over previous cagebreaks, need to be at least one div apart in the sequence
         DO J8=1, NTEMP3
         NTEMP2=NDIV-NDIVCB(I, J8)
CIF4 Check if we are looking at an even or odd correlation.
         NTEMP5=NTEMP3-J8 
         IF (MOD(NTEMP5,2).EQ.0) THEN
         NREVS(NTEMP2)=NREVS(NTEMP2)-1
         ELSE
         NREVS(NTEMP2)=NREVS(NTEMP2)+1
CFI4 End of correlation additions 
         ENDIF
         ENDDO
C End of cycle, update length of sequence to include current cagebreak
         NDIVCBC(I)=NDIVCBC(I)+1
         ELSE
         IF (REVING2.AND.NTEMP.NE.2) REM=REM+1
C         IF (REVING2.AND.NTEMP.EQ.2) THEN
C         NREVS(NTEMP)=NREVS(NTEMP)+1
C         NREVS(NTEMP-1)=NREVS(NTEMP-1)-1
C         ELSE
         NDIVCBC(I)=1
CFI3 Previous cage-break was not in a different time interval. Therefore end of sequence.
C    Start count at 1 again for this cagebreak
C         ENDIF 
         ENDIF
C     Not sure about rev2 at the moment.    
C         IF (REVING2.AND.NTEMP.EQ.2) THEN
C         NREVS(NTEMP)=NREVS(NTEMP)-1
C         ENDIF

C  Update list of cagebreaks with current break
         NTEMP2=NDIVCBC(I)  
         NDIVCB(I,NTEMP2)=NDIV
         
         ELSE
C Cagebreak reset count
         TIMES(I)=MINDUR
CFI  End of timediv if.              
         ENDIF


         DO J7=1, NCOUNTREV(I)
         IF(I.LE.NA) THEN
                 REVHIST(J7)=REVHIST(J7)+1
         ELSE        
                 REVHISTB(J7)=REVHISTB(J7)+1
         ENDIF
         ENDDO 
         IF (NCOUNTREV(I).GT.NREVHIST) NREVHIST=NCOUNTREV(I)
         IF (NCOUNTREV(I).GT.10000) WRITE(6,*) 'Error, too many revs'
         ELSE
C No break in this step
         ENDIF
        
         RETURN
         END 
