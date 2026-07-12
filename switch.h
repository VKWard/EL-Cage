CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Logical switches
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         
         INTEGER NLOOPSTART, NLOOP         
         LOGICAL ABSDIST, RELDIST, AATOMS, BATOMS
         LOGICAL TSCOORDS, TIMED, TIMEDIV
         LOGICAL CHANGEIN, CHANGEOUT, CHANGEOR 
         LOGICAL RESULTS, REMOVE 

           
         ABSDIST=.FALSE.      ! For neighbour loss use absolute distance
         RELDIST=.TRUE.       ! For neighbour loss use relative distance
         AATOMS=.TRUE.        ! Include A atoms
         BATOMS=.TRUE.        ! Include B atoms 
         TSCOORDS=.FALSE.     ! Read every other record (min/ts series)
         TIMED=.FALSE.        ! Read in times for each minimum
         TIMEDIV=.FALSE.      ! Consider time intervals of minima
         CHANGEIN=.TRUE.      ! Consider loss of neighbours 
         CHANGEOUT=.TRUE.     ! Consider gain of neighbours
         CHANGEOR=.FALSE.      ! Require loss/gain of ncb not a sum.
         
         RESULTS=.FALSE.       ! Read in temperature and output 1/t vs result files
        ! Removal of erroneous ts from a minmin traj also allowed with this option

C Loop extrema depend on the choice of atoms in the calculation         
         NLOOPSTART=1
         NLOOP=N
         IF (AATOMS.AND.BATOMS) THEN
         NLOOP=N 
         ELSE
         IF (AATOMS) NLOOP=NA
         IF (BATOMS) NLOOPSTART=NA+1
         ENDIF
