C*******************************************************************
C  CHECKNEIGH(DIST, DISTV, ATOMSTART, ATOMEND) checks neighbours
C*********************************************************************
C
         SUBROUTINE CHECKNEIGH(DIST, DISTV, DISTB, DISTVB, NBEGIN, NEND, K)
         IMPLICIT NONE
         DOUBLE PRECISION DIST, DISTV, DISTB, DISTVB
         DOUBLE PRECISION DISTC, DISTVC 
         INTEGER NBEGIN, NEND, I, J, K
         LOGICAL DEBUG
         INCLUDE 'system.h'
         INCLUDE 'counter.h'
         INCLUDE 'pos.h'
         INCLUDE 'switch.h' 
          
         DEBUG=.FALSE.
          
          DO I=NBEGIN,NEND

!          IF(I.EQ.60) THEN
!            DEBUG=.TRUE.
!          ELSE
!            DEBUG=.FALSE.
!          ENDIF

          DISTVC=DISTV
          DISTC=DIST
C Default is an AA pair, change these values if AB
           
C Clear the list of neighbours and extended neighbour list for this atom 
C Reset counters for neighbour lists to 1 
            NEIGH(I,:)=0
            EXTRA(I,:)=0
            NEXTRA=1
            NNEIGH=1
            NEXTRALAST=NEXTRAOLD(I)+1
C Check for nearest neighbours  
            DO J=1,N
            IF (J.GT.NA) DISTC=DISTB
C If the distance between test atom and neighbour is less than DISTAA 
C or DISTAB, as required, classify as a nearest neighbour
              
             IF (SQRT(POSDIST(I,J)).LE.DISTC) THEN
C Add neighbour to list
               NEIGH(I,NNEIGH)=J
!               IF(I.EQ.60) WRITE(11,*) J, "is a neighbour of 60. Distance= ", POSDIST(I,J)
                CALL FLUSH(11)
               NNEIGH=NNEIGH+1
C EXTRA contains an extended list of neighbours to compare with the 
C previous step
               EXTRA(I, NEXTRA)=J
               NEXTRA=NEXTRA+1
              ELSE 
C Add to EXTRA neighbours within 0.2 of the cutoff to allow for 
C vibrations if RELDIST is false, else check that the atom and its
C neighbour have not separated by more than half the pair equilibrium
C distance.
               IF (K.GT.1) THEN 
               IF (J.GT.NA) DISTVC=DISTVB
                 IF (.NOT.RELDIST) THEN
                   IF (SQRT(POSDIST(I,J)).LE.DISTC+DISTVC) THEN
                   EXTRA(I, NEXTRA)=J
                   NEXTRA=NEXTRA+1  
                   ENDIF
                 ELSE
!                   IF(I.EQ.60) THEN 
!                     WRITE(11,*) "Checking for neighbour change, 60 ", J
!                     WRITE(11,*) SQRT(POSDIST(I,J)), DISTC*2
!                     WRITE(11,*) POSDIST2(I,J), DISTVC
!                   ENDIF
                   IF (SQRT(POSDIST(I,J)).LE.DISTC*2.AND.POSDIST2(I,J).LT.DISTVC) THEN
!                   IF(I.EQ.60) WRITE(11,*) "Adding ", J, "to the EXTRA list for 60"
                   EXTRA(I, NEXTRA)=J
                   NEXTRA=NEXTRA+1  
                   ENDIF
                 ENDIF
               ENDIF 
             ENDIF
C Update Extraold to compare previous step with this step. 
             IF (SQRT(POSDISTOLD(I,J)).GT.DISTC) THEN
             IF (K.GT.1) THEN 
             IF (J.GT.NA) DISTVC=DISTVB
               IF (.NOT.RELDIST) THEN
                 IF (SQRT(POSDIST(I,J)).LE.DISTC+DISTVC) THEN
                 EXTRAOLD(I, NEXTRALAST)=J
                 NEXTRALAST=NEXTRALAST+1  
                 ENDIF
               ELSE
                 IF(SQRT(POSDISTOLD(I,J)).LE.DISTC*2.AND.POSDIST2B(I,J).LT.DISTVC) THEN
                 EXTRAOLD(I, NEXTRALAST)=J
                 NEXTRALAST=NEXTRALAST+1  
                 ENDIF
               ENDIF
             ENDIF 
           ENDIF
              
           ENDDO
C        IF (I.EQ.58) WRITE (11, *) I, NEIGH(I,1:NNEIGH)
         IF (DEBUG) WRITE (11, *) I, NEIGH(I,1:NNEIGH)
C Create an array for the step detailing the neighbour number for each atom  
         NNEIGHMAT(I)=NNEIGH-1
         IF (DEBUG) WRITE(6, *) 'Number of neighbours=', NNEIGH-2
C Add to the average for numbers of neighbours
C -2 to subtract self.
         IF (NBEGIN.eq.1) THEN
         AVNNEIGHA=AVNNEIGHA+NNEIGH-2
         ELSE
         AVNNEIGHB=AVNNEIGHB+NNEIGH-2
         ENDIF
         AVNNEIGH=AVNNEIGH+NNEIGH-2
         NEXTRAMAT(I)=NEXTRA-1 
         NEXTRAOLD(I)=NEXTRALAST-1 
         
         ENDDO
         RETURN
         END
