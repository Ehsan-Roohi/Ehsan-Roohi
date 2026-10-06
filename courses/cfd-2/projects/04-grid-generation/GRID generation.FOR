!****************************************************************************
!
!  PROGRAM: Grid Generation
!
!
!****************************************************************************
!============================================================================
!                                MAIN PROGRAM
!============================================================================
 	PROGRAM GRID

	!--DECLARE THE DIMENSION SIZE

	ImplIcIt DoublePrecIsIon  (A-H,O-Z)

	INTEGER,PARAMETER :: NX=51
	INTEGER,PARAMETER :: NY=51
 
	INTEGER ITYPE,I,J,C,ILE,COUNT

	REAL  X(NX,NY),Y(NX,NY),DXDSI(NX,NY),DXDETHA(NX,NY),
     &      DYDSI(NX,NY),DYDETHA(NX,NY),X1(NX,NY),Y1(NX,NY),
     & 	  ALPHA(NX,NY),BETHA(NX,NY),GAMMA(NX,NY) 

	REAL  DELTAI,DELTAPI,R,D,EY,ERR,OMEGA
	REAL  B,eth,a,u,f,T
	
	DoublePrecIsIon:: PI,EPS

	real:: xsi(NX,NY) , xeth(NX,NY) , ysi(NX,NY)
	real:: yeth(NX,NY) , jac(NX, NY) 
	real:: xss(NX) , xtt(NX) , yss(NX) , ytt(NX) 
	real:: xst(NX) , yst(NX) 
	real::  w , r1 , r2 , signx , signy
	real:: xsip(NX,NY) , xethp(NX,NY) , ysip(NX,NY)
	real:: yethp(NX,NY)

		 
      ! Body of the Code

	CALL INPUT	 
	CALL BC
	CALL INITIALIZE
	CALL ELIPTIC
	CALL OUTPUT	
       
	CONTAINS
!****************************************************************
	SubroutIne Input

	C=1
	PI=3.14159
	EPS=1.d-6
	T=0.12D0

	End SubroutIne
!****************************************************************
!      Define airfoil and circle boundary conditions
!
!****************************************************************
	SubroutIne BC

      ILE=(NX+1)/2

	 !AIRFOIL (J=1 REFERS TO AIRFOIL BOUNDARY CONDITION)

       DELTAI = (2.0/(NX - 1))
	 
	 !LEADING EDGE POSITION OF THE AIRFOIL 

	 X(ILE,1) = 0.0
	 Y(ILE,1) = 0.0

	 !TRAILING EDGE POSITION OF THE AIRFOIL 

	 X(NX,1) = 1.0
	 Y(NX,1) = 0.0
       
	 !TRAILING EDGE POSITION OF THE AIRFOIL  

       X(1,1) = 1.0
	 Y(1,1) = 0.0	 

	 !AIRFOIL THICKNESS (NACA 0012)

	 DO I=2,ILE-1

	 X(I,1)=1.0-(I-1)*DELTAI

	 !AIRFOIL FORMULA, LOWER EDGE

	 Y(I,1)=-(T/0.2D0*(0.2969D0*sqrt(X(I,1))-0.126D0*X(I,1)   
     & -0.3516D0*(X(I,1)**2)+0.2843D0*(X(I,1)**3)-0.1015D0*
     & (X(I,1)**4)))

	 END DO
								    
	 DO	I=ILE+1,NX-1

	 X(I,1)=1+(I-NX)*DELTAI

	!AIRFOIL FORMULA, UPPER EDGE

	 Y(I,1)=T/0.2D0*(0.2969D0*sqrt(X(I,1))-0.126D0*X(I,1)   
     & -0.3516D0*(X(I,1)**2)+0.2843D0*(X(I,1)**3)-0.1015D0*
     & (X(I,1)**4))

	 END DO	 
	 
c     ************** FAR FIELD BOUNDARY CONDITION **********

c      Note: J=NY represents the FAR FIELD BOUNDARY CONDITION

	 D=2*C
	 R=0.5*C+D
	 deltapi=PI*DELTAI

	 do I=1,ILE-1
	 X(I,NY)=0.5*C+R*COS((I-1)*deltapi)
	 Y(I,NY)=-R*SIN((I-1)*deltapi)
	 end do

	 do	I=ILE,NX
	 X(I,NY)=0.5*C-R*COS((I-ILE)*deltapi)
	 Y(I,NY)=R*SIN((I-ILE)*deltapi)
	 end do 

	 return
	 end subroutIne BC
!****************************************************************
	subroutine INITIALIZE

c	 In this subroutine the Initialization of the flow field is done by the linear Interpolation In
c	 the direction of the etha.
	
	Do I=1,NX
	Do J=2,NY-1

	EY=(J-1.0)/(NY-1.0)
	X(I,J)=(1-EY)*X(I,1)+EY*X(I,NY)
	Y(I,J)=(1-EY)*Y(I,1)+EY*Y(I,NY)

	End Do 
      End Do 

	End subroutine Initialize
!**************************************************************** 
	subroutine Eliptic
c
c	elliptic Grid Generation 

	
	ERR=1.0
	X1=X
	Y1=Y

	omega=0.9
	
		
	COUNT=1
	
     
	DO WHILE (ERR>EPS)

	   DO J=2,NY-1

	   DXDSI(1,J)=(X(2,J)-X(NX-1,J))/2.0		
	   DXDSI(NX,J)=DXDSI(1,J)   ! ASSIGNE FIRST POINT PROPERTIES TO LAST POINT
 
         DYDSI(1,J)=(Y(2,J)-Y(NX-1,J))/2.0
	   DYDSI(NX,J)=DYDSI(1,J)

	   DXDETHA(1,J)=(X(1,J+1)-X(1,J-1))/2.0
	   DXDETHA(NX,J)=DXDETHA(1,J)

	   DYDETHA(1,J)=(Y(1,J+1)-Y(1,J-1))/2.0
	   DYDETHA(NX,J)=DYDETHA(1,J)


	   ALPHA(1,J)=DXDETHA(1,J)**2+DYDETHA(1,J)**2   !Eq. 9-65
	   ALPHA(NX,J)=ALPHA(1,J)

	   BETHA(1,J)=DXDSI(1,J)*DXDETHA(1,J)+DYDSI(1,J)*DYDETHA(1,J) !Eq. 9-66
	   BETHA(NX,J)=BETHA(1,J)

	   GAMMA(1,J)=DXDSI(1,J)**2+DYDSI(1,J)**2 !Eq. 9-67
	   GAMMA(NX,J)=GAMMA(1,J)
		
	   END DO
		 
	   DO I=2,NX-1
	   DO J=2,NY-1

	   DXDSI(I,J)=(X(I+1,J)-X(I-1,J))/2.0
	   DYDSI(I,J)=(Y(I+1,J)-Y(I-1,J))/2.0

	   DXDETHA(I,J)=(X(I,J+1)-X(I,J-1))/2.0
	   DYDETHA(I,J)=(Y(I,J+1)-Y(I,J-1))/2.0

	   ALPHA(I,J)=DXDETHA(I,J)**2+DYDETHA(I,J)**2

	   BETHA(I,J)=DXDSI(I,J)*DXDETHA(I,J)+DYDSI(I,J)*DYDETHA(I,J)

	   GAMMA(I,J)=DXDSI(I,J)**2+DYDSI(I,J)**2

	   END DO 
	   END DO 

	   DO J=2,NY-1

	   X1(1,J)=(ALPHA(1,J)*(X1(2,J)+X1(NX-1,J))+
     &   GAMMA(1,J)*(X1(1,J+1)+X1(1,J-1))-0.5*BETHA(1,J)*
     &   (X1(2,J+1)-X1(NX-1,J+1)-X1(2,J-1)+ 
     &   X1(NX-1,J-1)))/(2*(ALPHA(1,J)+GAMMA(1,J)))
	 
 	   Y1(1,J)=(ALPHA(1,J)*(Y1(2,J)+Y1(NX-1,J))+
     &   GAMMA(1,J)*(Y1(1,J+1)+Y1(1,J-1))-0.5*BETHA(1,J)*
     &   (Y1(2,J+1)-Y1(NX-1,J+1)-Y1(2,J-1)+
     &   Y1(NX-1,J-1)))/(2*(ALPHA(1,J)+GAMMA(1,J)))
			
	   X1(NX,J)=X1(1,J)
	   Y1(NX,J)=Y1(1,J)			   	
	   END DO
		 
	   DO I=2,NX-1
	   DO J=2,NY-1

	   X1(I,J)=(ALPHA(I,J)*(X1(I+1,J)+X1(I-1,J))+
     &   GAMMA(I,J)*(X1(I,J+1)+X1(I,J-1))-0.5*BETHA(I,J)*
     &   (X1(I+1,J+1)-X1(I-1,J+1)-X1(I+1,J-1)+X1(I-1,J-1)))
     &   /(2*(ALPHA(I,J)+GAMMA(I,J)))

	   Y1(I,J)=(ALPHA(I,J)*(Y1(I+1,J)+Y1(I-1,J))+
     &   GAMMA(I,J)*(Y1(I,J+1)+Y1(I,J-1))-0.5*BETHA(I,J)*
     &   (Y1(I+1,J+1)-Y1(I-1,J+1) - Y1(I+1,J-1)+Y1(I-1,J-1)))
     &   /(2*(ALPHA(I,J)+GAMMA(I,J)))
			
         END DO 
	   END DO 

	   X1=OMEGA*X1+(1-OMEGA)*X
	   Y1=OMEGA*Y1+(1-OMEGA)*Y
	   
	   ER1=MAXVAL(ABS(X1-X))
	   ER2=MAXVAL(ABS(Y1-Y))
	   Err=MAX(ER1,ER2)		 

	   X=X1
	   Y=Y1
	   COUNT=COUNT+1

	END DO	
      
	END Subroutine
!****************************************************************
      Subroutine Output

c      Printing the result

	 Open(1,file='GRID1.plt')
	 Write(1,*) 'Variable = "X","Y"'
	 Write(1,5) NX,NY
5      FORMAT('  ZONE ','   i=',I5,'   j=',I5)
	 
	 DO i=1,NX
	 DO j=1,NY
	 Write(1,*) X(i,j),Y(i,j)
	 END DO 
	 END DO 

	End Subroutine
!****************************************************************
	End Program GRID
