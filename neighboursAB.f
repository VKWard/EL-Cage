C**************************************************************
C Program to find the nearest neighbours of each atom and check 
C for changes in consecutive steps.
C**************************************************************

         PROGRAM FINDNEIGHDISTAB
         IMPLICIT NONE
         INCLUDE 'system.h'     ! specifies box size and atom numbers        
         INCLUDE 'pos.h'    ! common group for positions data
         DOUBLE PRECISION CORRL
         DOUBLE PRECISION TEMP, TEMPI, TEMP2 
         
         INTEGER I1, I, J, K
         INTEGER J2, J3, J1, J4, J7
         INTEGER J8, I4
         INTEGER NLINES, NREC, NRUN, J9, NDIV
         INTEGER I2, I3, NUMTS, NTS, I7
         INTEGER, ALLOCATABLE :: TS(:)
         LOGICAL DEBUG
   
         INCLUDE 'criteria.h'   ! specifies cutoffs etc..
         INCLUDE 'counter.h' ! counters for events such as cagebreaks
         INCLUDE 'switch.h' ! options, atom type and program options 
         DEBUG=.FALSE. 
 
         OPEN(UNIT=10,FILE='SETUP',STATUS='OLD')
         IF (RESULTS) THEN
         OPEN(UNIT=20,FILE='PEVENT',STATUS='NEW')
         OPEN(UNIT=30,FILE='PREV',STATUS='NEW')
         ENDIF
         READ(10,*)NLINES, NFORMAT
C Loop over all runs of different energies
         DO J9=1, NLINES
         NUMTS=0
         REMOVE=.FALSE. 
         NTS=0        

         IF (RESULTS) THEN
         READ(10,*) NRUN, NREC, TEMP, NUMTS
         IF (NUMTS.GT.0) THEN
         REMOVE=.TRUE.
         NTS=1
         ALLOCATE(TS(NUMTS))
         READ(10,*) (TS(I4), I4=1,NUMTS)
         ENDIF
         ELSE 
         READ(10,*) NRUN, NREC
         IF(DEBUG) WRITE(*,*) "NREC:", NREC 
         ENDIF 
         
         NDIV=1

C Opens input and output files, zeros counters and falsifies logic tests        
         
         CALL IO1(NRUN)   

C Cycle through all the records of atomic positions
         IF (NFORMAT.EQ.3) NREC=NREC*2
         DO I2=1,NREC
         IF (REMOVE) THEN
         IF (TS(NTS).LE.I2) THEN
         IF (NTS.GE.NUMTS) REMOVE=.FALSE. 
         NTS=NTS+1
         CYCLE 
         ENDIF
         ENDIF
         
         IF (DEBUG) WRITE(11, *) 'Loop over records, I2:', I2
         CALL FLUSH(11)
         IF (TIMED) CALL TIMERS(.FALSE.) 

C Change index so that we can use a different step if we need to exclude TS in a series
         I3=I2
         IF (TSCOORDS) I3=2*I2-1

C Reads in coords and calculates vectors between atoms
         IF (DEBUG) WRITE(11, *) 'Call IO2'
         CALL FLUSH(11)
         CALL IO2(I3, NREC)
         IF (DEBUG) WRITE(11, *) 'Back from IO2'
         CALL FLUSH(11)
C Either use a nearest neighbour analysis or an absolute distance criterion
C to decide if a cage-breaking event has taken place.

      IF (.NOT.ABSDIST) THEN
C***************************************************************        
C Check the neighbours of each A atom         
C***************************************************************
      IF (DEBUG) WRITE(11, *) 'Call checkneigh'
      CALL FLUSH(11)
      CALL CHECKNEIGH(DISTAA, DISTVAA, DISTAB, DISTVAB, 1, NA, I2) 

C***************************************************************        
C Check the neighbours of each B atom         
C***************************************************************
      IF (NB.GT.0) THEN
      CALL CHECKNEIGH(DISTAB, DISTVAB, DISTBB, DISTVBB, NA+1, N, I2)
      ENDIF 
         
C        IF (DEBUG) WRITE(11, *) 'Number of neighbours=', NNEIGH
C        IF (DEBUG) WRITE(11, *) 'Checking nearest neighbours'
        CALL FLUSH(11)

C Compare nearest neighbours in this step to the previous step
  
       ENDIF
       IF (NFORMAT.EQ.3.AND.MOD(I2,2).NE.0) CYCLE
       IF (I2.GT.1) THEN

C Cycle over all atoms checking for cagebreaks and reversals
        IF (DEBUG) WRITE(11, *) 'Cycle for cagebreak, record number=', I2
        DO I=NLOOPSTART,NLOOP
        CALL CAGEBREAK(I2,I)
        ENDDO 
       
       ENDIF

C Summarise for this step
        IF (NATOMSCB.GT.0) THEN
         WRITE(13, *) I2, NATOMSCB
         NCB=NCB+1 
         NALL=NALL+NATOMSCB
         if(NATOMSCB_A.GT.0) then
                 NALL_A=NALL_A+NATOMSCB_A
C                PRINT*, 'neighboursAB', ' NPROD_A= ', NPROD_A
         endif
         if(NATOMSCB_B.GT.0) then 
                NALL_B=NALL_B+NATOMSCB_B
C                PRINT*, 'neighboursAB', ' NPROD_B= ', NPROD_B
         endif
        ENDIF  
        IF (CBREVERSE) THEN 
           NREV1=NREV1+1
        IF (CBREVERSE2) NREVBOTH=NREVBOTH+1 
        ELSE 
        IF (CBREVERSE2) NREV2=NREV2+1
        ENDIF
C print R2 for this step:
        WRITE(17,111) I2,R2ADD,R2ADD_A,R2ADD_B,R2ADDALL,R2ADDALL_A,R2ADDALL_B        
111          FORMAT(I8, 6F19.10) 
        ENDDO   ! sn402: end loop over I2

         AVNNEIGH=AVNNEIGH/DBLE(NREC*N)
         WRITE(6, *) 'Average number of neighbours=', AVNNEIGH
         AVNNEIGHA=AVNNEIGHA/DBLE(NREC*NA)
         WRITE(6, *) 'Average for an A atom=', AVNNEIGHA
         AVNNEIGHB=AVNNEIGHB/DBLE(NREC*NB)
         WRITE(6, *) 'Average for an B atom=', AVNNEIGHB
         WRITE(6, *) 'Number of cage-breaking steps=', NCB
         WRITE(6, *) 'Number of cage-breaking atoms=', NALL
         write(6, *) 'Number of cage-breaks Li=', NALL_A
         write(6, *) 'Number of cage-breaks Cl=', NALL_B
         write(6, *) 'Number of productive CBs Li=', NPROD_A
         write(6, *) 'Numver of productive CBs Cl=', NPROD_B
         WRITE(6, *) 'Number of reverse steps,type1=', NREV1
         WRITE(6, *) 'Number of reverse steps,type2=', NREV2
         WRITE(6, *) 'Number of reverse atoms,type1=', NREV1ALL
         WRITE(6, *) 'Number of reverse atoms,type2=', NREV2ALL
         WRITE(6, *) 'Number of steps type1 + type2=', NREVBOTH
         WRITE(6, *) 'Total r2:  irreversible steps=', R2ADD
         WRITE(6, *) 'Total r2: irrev, only A atoms=', R2ADD_A
         WRITE(6, *) 'Total r2: irrev, only B atoms=', R2ADD_B
         WRITE(6, *) 'Total r2:           all steps=', R2ADDALL
         WRITE(6, *) 'Total r2:   all, only A atoms=', R2ADDALL_A
         WRITE(6, *) 'Total r2:   all, only B atoms=', R2ADDALL_B
         CORRL=0D0
C        Find correlation effect of reversals
         DO J7=1, NREVHIST
         CORRL=((-1)**J7)*REVHIST(J7)+CORRL
         ENDDO 
         WRITE(6,*)  'Correlation by reversal numbers=', CORRL
         CORRL=0D0
         DO J7=1, NREVHIST
         CORRL=((-1)**J7)*REVHISTB(J7)+CORRL
         ENDDO 
         WRITE(6,*)  'Correlation of B by reversal numbers=', CORRL
         WRITE(6, *) 'REVERSALS - A only'
         WRITE(6, *) 'length: number'
         DO J7=1, NREVHIST
         WRITE(6, *) J7, REVHIST(J7)
         ENDDO   
 
         IF (RESULTS) THEN
         TEMPI=1/TEMP
         TEMP=LOG(DBLE(NCB)/DBLE((NREC-1-NUMTS)))
         TEMP2=LOG(DBLE(NALL)/DBLE((NREC-1-NUMTS)*(NLOOP-NLOOPSTART+1)))
         WRITE(20,*) TEMPI, TEMP, TEMP2 
         TEMP=DBLE(-NREV1ALL-NREV2ALL)/DBLE(NALL)
         TEMP2=CORRL/DBLE(NALL)
         WRITE(30,*) TEMPI, TEMP, TEMP2 
         ENDIF
         
         IF (TIMED) THEN
         WRITE(6,*) NREVS(1:100)
         WRITE(6,*) REM 
         ENDIF
      
         CLOSE(11)
         CLOSE(12)
         CLOSE(13)
         CLOSE(14)
         CLOSE(15)
         CLOSE(16)
         CLOSE(17)
         close(18)

         ENDDO
C End of loop over all runs of different temperatures
         CLOSE(10)
         END 
