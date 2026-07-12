C*******************************************************************
C  Timers is used when TIMED is true, read in times for each min.
C*********************************************************************
C
         SUBROUTINE TIMERS(RESET)
         IMPLICIT NONE
         INCLUDE 'system.h'
         INCLUDE 'timing.h'
         DOUBLE PRECISION MINDUR, TIMES
         INTEGER NDIVCH, NDIV, I
         LOGICAL RESET
         
         INCLUDE 'switch.h' 
   

       
         READ(16, *)MINDUR
C Update time with residence time in min before current transition.
C We can work out the time before a cagebreak for each atom using this measure.
         TIMES(:)=TIMES(:)+MINDUR

C Timediv allows us to consider time intervals of minima.
C i.e. coarse grain the results to certain time intervals

           IF(TIMEDIV) THEN

CIF2 Have we reached a new time interval
             DO I=1,N
             IF (TIMES(I).GT.DIV) THEN
C Update to new interval and time in new interval
C MORE THAN ONE DIVISION?
             NDIVCH=NINT(TIMES(I)/DIV)
             IF((TIMES(I)-DIV*NDIVCH).GT.0) THEN           
             TIMES(I)=TIMES(I)-DIV*NDIVCH
           
             IF(I.EQ.1) NDIV=NDIV+NDIVCH
             ELSE
             TIMES(I)=TIMES(I)-DIV*(NDIVCH-1)
             IF(I.EQ.1) NDIV=NDIV+(NDIVCH-1)
             ENDIF
             ENDIF
             ENDDO
           ENDIF

           RETURN
           END
