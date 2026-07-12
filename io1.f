C*******************************************************************
C  IO1(run_number) performs some IO and initial setup
C*********************************************************************
C
         SUBROUTINE IO1(I1)
         IMPLICIT NONE
         INCLUDE 'system.h' 
         INCLUDE 'pos.h'
         CHARACTER*20 change, cbreak, cbrev1, cbrev2, energy, coords, r2
         character*20 jumps, cbpos, rjumps
         INTEGER I1
         INCLUDE 'counter.h'
         INCLUDE 'switch.h'
C
C        UNIT numbers are defined here.
C        RUN numbers and step numbers for each run     10 !opened in main program
C        RESULTS files if required (RESULTS=.TRUE.)    20,30 !opened in main program
C        Input coordinates                             12
C        Input energies                                16
C        List of cagebreaks                            13
C        List of reversed cagebreaks of type 1         14
C        List of reversed cagebreaks of type 2         15
C        List of neighbour changes                     11
C        List of r2 per step                           17
         
         IF (I1.LT.10) THEN
         WRITE(change,80)'nnchange.',I1
C        WRITE(energy,83)'quench.energy.',I1
         WRITE(coords,81)'pathcoords.',I1
         WRITE(cbreak,82)'break.',I1
         WRITE(cbrev1,81)'break.rev1.',I1
         WRITE(cbrev2,81)'break.rev2.',I1
         WRITE(r2,80)'rsquared.',I1
         write(jumps,81)'jumpwidths.',I1
         write(cbpos,81)'cbposition.',I1
         write(rjumps,82)'rjump.',I1
82       FORMAT(A6,I1)
80       FORMAT(A9,I1)
81       FORMAT(A11,I1)
83       FORMAT(A14,I1)
         ELSE
         IF (I1.LT.100) THEN
         WRITE(change,70)'nnchange.',I1
C        WRITE(energy,73)'quench.energy.',I1
         WRITE(coords,71)'pathcoords.',I1
         WRITE(cbreak,72)'break.',I1
         WRITE(cbrev1,71)'break.rev1.',I1
         WRITE(cbrev2,71)'break.rev2.',I1
         WRITE(r2,70)'rsquared.',I1
72       FORMAT(A6,I2)
70       FORMAT(A9,I2)
71       FORMAT(A11,I2)
73       FORMAT(A14,I2)
         ELSE
         WRITE(change,90)'nnchange.',I1
C        WRITE(energy,93)'quench.energy.',I1
         WRITE(coords,91)'pathcoords.',I1
         WRITE(cbreak,92)'break.',I1
         WRITE(cbrev1,91)'break.rev1.',I1
         WRITE(cbrev2,91)'break.rev2.',I1
         WRITE(r2,90)'rsquared.',I1
92       FORMAT(A6,I4)
90       FORMAT(A9,I4)
91       FORMAT(A11,I4)
93       FORMAT(A14,I4)
         ENDIF  
         ENDIF  


         OPEN(UNIT=11,FILE=change,STATUS='NEW')
      IF (NFORMAT.eq.3) OPEN(UNIT=22,FILE='ts.data',STATUS='OLD')
      IF (NFORMAT.eq.3) OPEN(UNIT=12,FILE='points.min',ACCESS='DIRECT',RECL=8*N*3)
      IF (NFORMAT.eq.2) OPEN(UNIT=12,FILE=coords,STATUS='OLD')
      IF (NFORMAT.eq.1) OPEN(UNIT=12,FILE=coords,ACCESS='DIRECT',RECL=8*N*3)
         OPEN(UNIT=13,FILE=cbreak,STATUS='NEW')
         OPEN(UNIT=14,FILE=cbrev1,STATUS='NEW')
         OPEN(UNIT=15,FILE=cbrev2,STATUS='NEW')
C        OPEN(UNIT=16,FILE=energy,STATUS='OLD')
         OPEN(UNIT=17,FILE=r2,STATUS='NEW')
         open(unit=18,file=jumps, status='new')
         open(unit=19,file=cbpos, status='new')
         open(unit=20,file=rjumps, status='new')
        
C
C     Zero counters
C 
         NEIGH(:,:)=0
         EXTRA(:,:)=0
         AVNNEIGH=0.0D0 
         AVNNEIGHA=0.0D0
         AVNNEIGHB=0.0D0
         NCB=0
         NREV1=0
         NREV2=0
         NALL=0
         NALL_A=0
         NALL_B=0
         NPROD_A=0
         NPROD_B=0
         NREV1ALL=0
         NREV2ALL=0
         NREVBOTH=0
         NATOMSCB=0
         NATOMSCB_A=0
         NATOMSCB_B=0
         NATOMSCBPROD_A=0
         NATOMSCBPROD_B=0
         NCOUNTREV(:)=0
         REVHIST(:)=0
         REVHISTB(:)=0
         NREVHIST=0
         NREVS(:)=0
         REM=0 
         R2ADD=0.0D0     ! Total r2 for irreversible atomic moves only.
         R2ADD_A=0.0D0   ! ,, only A atoms
         R2ADD_B=0.0D0   ! ,, only B atoms
         R2ADDALL=0.0D0  ! Total r2 for all atomic moves.
         R2ADDALL_A=0.0D0
         R2ADDALL_B=0.0D0
C
C     Falsify logical tests 
C
         CBREVERSE=.FALSE.
         CBREVERSE2=.FALSE.

         RETURN
         END
