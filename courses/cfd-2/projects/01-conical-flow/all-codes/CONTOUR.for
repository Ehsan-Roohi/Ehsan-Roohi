!****************************************************************************
!
!  PROGRAM: Conical Flow
!
!
!****************************************************************************
!============================================================================
!                                MAIN PROGRAM
!============================================================================
	PROGRAM CONICALFLOW 

      !--DECLARE THE DIMENSION SIZE

      implicit DoublePrecision  (A-H,O-Z)

      integer,parameter :: N=240,N2=2618	
  
	DoublePrecision :: TETA(N2),TETAF(N),RHO(N),RHOO(N),RHOF(N),P(N)	
     &	,PO(N),PF(N),U(N),T(N2),F1(11,N),F2(11,N),F3(11,N),DT(N2) 

      JJ=11
	Open(21,FILE='TEMP.DAT')
	OPEN (31,FILE='TETA.DAT') 
	OPEN (33,FILE='DT.DAT') 

	OPEN (32,FILE='FP.DAT')
	WRITE (32,*) 'variables= "TETA" "TIME" "P" ' 
      WRITE (32,*) 'ZONE i=',238,' J=',JJ,'F=POINT'

	

	DO J=1,JJ
	
	K1=(N-2)*(J-1)+1
	K2=(N-2)*(J)

	M=1

	DO I=K1,K2

      READ (21,*) T(I)
      READ (31,*) TETA(I)
	READ (33,*) DT(I)
	F1(J,M)=TETA(I)
	F2(J,M)=DT(I)
	F3(J,M)=T(I)
C      F2(J,M)=TETA(I)      
      
	WRITE (32,*)  F1(J,M),F2(J,M),F3(J,M)
	M=M+1
      
	END DO
	WRITE (*,*) J,M !K1,K2

      END DO
	
	End program 

