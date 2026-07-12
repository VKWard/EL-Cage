***
C
         PROGRAM REMOVEBOX
         IMPLICIT NONE
         CHARACTER(len=20) positions,msddiff, new_positions 
         DOUBLE PRECISION,ALLOCATABLE ::  POS(:,:), NEWPOS(:,:)  
         DOUBLE PRECISION BOXLX, BOXLY, BOXLZ, ddump
         DOUBLE PRECISION x1, x2, y1, y2, z1, z2 
         INTEGER N,J,J1,J2, I,I1,K2,K3,K4,K6, K7, NTIME
         INTEGER num_args, ix, ndump, nread, nrecs
         INTEGER NLi, NCl, K
         CHARACTER(len=20), dimension(:), allocatable :: args
         CHARACTER(len=100) line, testline


         OPEN(UNIT=4,FILE='boxlengths', STATUS='REPLACE') 
         OPEN(UNIT=44,FILE='nrecs', STATUS='REPLACE') 
         OPEN(UNIT=34,FILE='ionnumbers', STATUS='REPLACE') 
         NTIME=2500
         N=1720
         BOXLX=3.9498129957D0+0.0455129957
         BOXLY=3.9397447249D0+0.0354447249
         BOXLZ=12.734592886D0+0.1355928862

         WRITE(positions,'(A11)') 'oldfile.xyz'
         WRITE(new_positions,'(A11)')'newfile.xyz'
         OPEN(UNIT=1,FILE=positions,STATUS='OLD')
         OPEN(UNIT=2,FILE=new_positions,STATUS='UNKNOWN')
         READ(1,fmt='(a)') testline
         N=-8
         DO
         READ(1,fmt='(a)') line
         IF (line.eq.testline) EXIT
         N=N+1 
         ENDDO
         CLOSE(1)
         OPEN(UNIT=1,FILE=positions,STATUS='OLD')
         ALLOCATE(POS(N*3,NTIME))
         ALLOCATE(NEWPOS(N*3,NTIME))
         NEWPOS(:,:)=0
         
         nrecs=0
         DO J1=1,1
         READ(1,*, END=22)
         nrecs=nrecs+1
         READ(1,*)
         READ(1,*)
         READ(1,*)
         READ(1,*)
         READ(1,*) x1, x2
         READ(1,*) y1, y2
         READ(1,*) z1, z2
         READ(1,*)
         BOXLX=x2-x1
         BOXLY=y2-y1
         BOXLZ=z2-z1
         IF (J1.EQ.1) WRITE(4,*) '     DOUBLE PRECISION :: BOXLX='
     &                                  ,BOXLX  
     &                                ,',BOXLY=',BOXLY,',BOXLZ=', BOXLZ
         nread=0
         NLi=0 
         NCl=0
         DO K2=1,N
         nread=nread+1
         READ(1,*)J2,ndump,ddump,(POS(I,J1), I=3*nread-2,3*nread)
         POS(3*nread-2,J1)=POS(3*nread-2,J1)-x1
         POS(3*nread-1,J1)=POS(3*nread-1,J1)-y1
         POS(3*nread,J1)=POS(3*nread,J1)-z1
         IF(ndump.eq.3) THEN
                nread=nread-1
                CYCLE
         ENDIF    
         IF(ndump.eq.1) THEN
                 NLi=NLi+1
         ENDIF        
         IF(ndump.eq.2) THEN
                 NCl=NCl+1
         ENDIF        
         ENDDO
         ENDDO
         NEWPOS(:,1)=POS(:,1)
         WRITE(34,*) '     INTEGER,PARAMETER :: N=',nread,', NA=',NLi,
     &                             ', NB=',NCl 

         DO J1=2,NTIME
                READ(1,*, END=22)
                nrecs=nrecs+1
                DO K=1,8
                        READ(1,*)
                ENDDO
         DO K2=1,N
         READ(1,*)J2,ndump,ddump,(POS(I,J1), I=3*K2-2,3*K2)
         POS(3*K2-2,J1)=POS(3*K2-2,J1)-x1
         POS(3*K2-1,J1)=POS(3*K2-1,J1)-y1
         POS(3*K2,J1)=POS(3*K2,J1)-z1
         IF(ndump.eq.3) THEN
                CYCLE
         ENDIF    
         ENDDO
         ENDDO
         
22       DO K3=2,nrecs
         DO K4=1,nread
         J=3*K4-2 
         NEWPOS(J,K3)=POS(J,K3)+
     &     ANINT((NEWPOS(J,K3-1)-POS(J,K3))/BOXLX)*BOXLX
         J=3*K4-1 
         NEWPOS(J,K3)=POS(J,K3)+
     &     ANINT((NEWPOS(J,K3-1)-POS(J,K3))/BOXLY)*BOXLY
         J=3*K4 
         NEWPOS(J,K3)=POS(J,K3)+
     &     ANINT((NEWPOS(J,K3-1)-POS(J,K3))/BOXLZ)*BOXLZ
         ENDDO
         ENDDO
 
         WRITE(44,*) 1, nrecs 
         DO K6=1,nrecs
         DO K7=1,nread
         WRITE(2,*)(NEWPOS(I,K6), I=3*K7-2,3*K7)
         ENDDO 
         ENDDO
         END
