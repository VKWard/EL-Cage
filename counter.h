C********************************************************************
C  Common group of counters
C********************************************************************

      DOUBLE PRECISION AVNNEIGH, AVNNEIGHA, AVNNEIGHB  
      INTEGER NCB, NREV1, NREV2, NALL, NREV1ALL, NREV2ALL, NREVBOTH
      INTEGER NALL_A, NALL_B, NPROD_A, NPROD_B
      INTEGER REM
      INTEGER*4 NEXTRA, NNEIGH, NATOMSCB, NATOMSCB_A, NATOMSCB_B, NREVHIST, NEXTRALAST
      INTEGER NATOMSCBPROD_A, NATOMSCBPROD_B
      LOGICAL CBREVERSE, CBREVERSE2

      COMMON /COUNTER/ AVNNEIGH, AVNNEIGHA, AVNNEIGHB
C********************************************************************
C AVVNEIGH/A/B  Average number of neighbours for an atom (all,A or B)
C********************************************************************

      COMMON /COUNTER2/ NCB, NREV1, NREV2, NALL, NREV1ALL, NREV2ALL, 
     &            REM, NREVBOTH, NEXTRA, NNEIGH, NATOMSCB, NREVHIST,
     &            NATOMSCB_A, NATOMSCB_B, NALL_A, NALL_B,
     &		  NPROD_A, NPROD_B

C********************************************************************
C NCB         Number of cagebreaks in each run 
C NREV1/2     Number of reversible cagebreaks of type 1/2 in each run
C NALL        Number of atomic cagebreaks in each run (all, A, B atoms)
C NPROD_A/B   Numver of productive atomic cagebreaks (A,B atoms)
C NREV1/2ALL  Number of atomic reversals of type 1/2 in total, each run
C NREVBOTH    Number of reversals with both type 1&2 occuring, each run
C REM         
C Set to zero at start of run in io1.f
C NATOMSCB    Number of atoms cagebreaking in each step.
C NREVHIST    Maximum number of reversals
C Set to zero for each step in io2.f
C********************************************************************
      INTEGER*4 :: NCOUNTREV, REVHIST, REVHISTB, NREVS
      COMMON /COUNTER3/ NCOUNTREV(N), REVHIST(10000),
     & REVHISTB(10000),NREVS(1000) 

C********************************************************************
C NCOUNTREV  Number of atomic reversals in series in total.
C REVHIST
C REVHISTB
C NREVS     
C Set to zero at start of run in io1.f and after any irreversible step.
C Array is declared and allocated in main program.
C********************************************************************
                     
      COMMON /LOGTESTS/ CBREVERSE, CBREVERSE2

C********************************************************************
C CBREVERSE/2 True if any reversals of type 1/2 take place in a step
C These are falsified for each step in io2.f
C********************************************************************  
      INTEGER*4 :: NEIGH, EXTRA, NEXTRAMAT, NNEIGHMAT      
      COMMON /NEIGHBOUR/ NEIGH(1:N,1:500), EXTRA(1:N,1:500), NEXTRAMAT(N), NNEIGHMAT(N)

C********************************************************************
C NEIGH      Array giving list of neighbours for each atom. 
C EXTRA      Array giving an extended list of neighbours for each atom.
C NNEIGHMAT  Array showing the number of records in neigh for each atom.
C NEXTRAMAT  Array showing the number of records in extra for each atom.
C********************************************************************
      INTEGER*4 :: NEIGHOLD, NNEIGHOLD
      COMMON /NEIGHOLD/ NEIGHOLD(1:N,1:500), NNEIGHOLD(N)

C********************************************************************
C NEIGHOLD   Array giving list of neighbours for each atom (last step). 
C NNEIGHOLD  Array showing the number of records in neigh for each atom
C            (last step).
C********************************************************************
      
      INTEGER*4 :: EXTRAOLD, NEXTRAOLD
      COMMON /EXTRAOLD/ EXTRAOLD(1:N,1:500), NEXTRAOLD(N)

C********************************************************************
C EXTRAOLD   Array giving extended list of neighbours for each atom 
C            (last step). 
C NEXTRAOLD  Array showing the number of records in extra for each atom
C            (last step).
C********************************************************************
 
      INTEGER*4  NFORMAT 
      COMMON /FFMAT/ NFORMAT
C*******************************************************************
C  options are in SETUP: 1:xmol,2:xyz
C*******************************************************************
