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

	INTEGER,PARAMETER :: NX=25
	INTEGER,PARAMETER :: NY=25
 
	INTEGER ITYPE,I,J,C,ILE,COUNT

	REAL  X(NX,NY),Y(NX,NY),DXDSI(NX,NY),DXDETHA(NX,NY),
     &      DYDSI(NX,NY),DYDETHA(NX,NY),X1(NX,NY),Y1(NX,NY),
     & 	  ALPHA(NX,NY),BETHA(NX,NY),GAMMA(NX,NY) 

	REAL  DELTAI,DELTAPI,R,D,EY,ERROR,OMEGA
	REAL  B,eth,a,u,f,T
	
	DoublePrecIsIon:: PI,EPS

	real:: xsi(NX,NY) , xeth(NX,NY) , ysi(NX,NY)
	real:: yeth(NX,NY) , jac(NX, NY) 
	real:: xss(NX) , xtt(NX) , yss(NX) , ytt(NX) 
	real:: xst(NX) , yst(NX) , p1(NX) , q1(NX)
	real::  w , r1 , r2 , signx , signy
	real:: xsip(NX,NY) , xethp(NX,NY) , ysip(NX,NY)
	real:: yethp(NX,NY)
	real:: sa,sb,seth,wp,wq,aa,bb

		 
      ! Body of the Code

	CALL INPUT	 
	CALL BC
	CALL INITIALIZE
	CALL ELIPTIC
	CALL ORTHOGONAL 
	CALL OUTPUT	
       
	CONTAINS
!****************************************************************
	SubroutIne Input

	C=1
	PI=3.14159
	EPS=1.d-6
	T=0.12D0
	      
	sa=0.004
	sb=0.004
	seth=0.005
	wp=0.00001
	wq=0.00001
	aa=1
	bb=2

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

	 Y(I,1)=-(T/0.2D0*(0.2969D0*sqrt(X(I,1))-0.126D0*X(I,1)   
     & -0.3516D0*(X(I,1)**2)+0.2843D0*(X(I,1)**3)-0.1015D0*
     & (X(I,1)**4)))

	 END DO
								    
	 DO	I=ILE+1,NX-1

	 X(I,1)=1+(I-NX)*DELTAI

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

	
	ERROR=1.0
	X1=X
	Y1=Y
	
	WRITE(*,*) 'ENTER THE RELAXATION FACTOR>' 
	READ(*,*)   OMEGA
	
	COUNT=1
	
	OPEN(2,file='ERROR.plt')
      
	DO WHILE (ERROR>EPS)

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
	   Error=MAX(MAXVAL(ABS(X1-X)),MAXVAL(ABS(Y1-Y)))
		 
	   WRITE(2,*) ERROR

	   X=X1
	   Y=Y1
	   COUNT=COUNT+1

	END DO	
      
	END Subroutine 
!****************************************************************
	Subroutine orthogonal

C	Applying orthogonality on the airfoil and circle using control functions. 

	error = 1.0
	p1 = 0
	q1 = 0
	x1 = x
	y1 = y

C	print*, 'Enter the relaxation factor>' 
C	read*, w
	
	count = 1
	open(2,file='error.plt')
      
	do while(error > EPS)
			
		xss(1) = x1(2,1) - 2*x1(1,1) + x1(NX-1,1)
 		yss(1) = y1(2,1) - 2*y1(1,1) + y1(NX-1,1)
		xsi(1,1) = (x1(2,1) - x1(NX-1,1))/2.0
		ysi(1,1) = (y1(2,1) - y1(NX-1,1))/2.0
	 	xeth(1,1) = -seth*ysi(1,1)/((xsi(1,1)**2 + ysi(1,1)**2)**0.5)
		yeth(1,1) = seth*xsi(1,1)/((xsi(1,1)**2 + ysi(1,1)**2)**0.5)
		xtt(1) = (-7.0*x1(1,1) + 8.0*x1(1,2) - x1(1,3))/2.0
     &				- (3.0*xeth(1,1))
		ytt(1) = (-7.0*y1(1,1) + 8.0*y1(1,2) - y1(1,3))/2.0
     &				- (3.0*yeth(1,1))
		xss(NX) = xss(1)
 			yss(NX) = yss(1)
		xsi(NX,1) = xsi(1,1)
		ysi(NX,1) = ysi(1,1) 
		xeth(NX,1) = xeth(1,1)
		yeth(NX,1) = yeth(1,1)
		xtt(NX) = xtt(1) 
		ytt(NX) = ytt(1)
		
		do i = 2,(NX-1)
			xss(i) = x1(i+1,1) - 2*x1(i,1) + x1(i-1,1)
			yss(i) = y1(i+1,1) - 2*y1(i,1) + y1(i-1,1)
			xsi(i,1) = (x1(i+1,1) - x1(i-1,1))/2.0
			ysi(i,1) = (y1(i+1,1) - y1(i-1,1))/2.0
			xeth(i,1) = -seth*ysi(i,1)/((xsi(i,1)**2 
     &					+ ysi(i,1)**2)**0.5)
			yeth(i,1) = seth*xsi(i,1)/((xsi(i,1)**2 
     &					+ ysi(i,1)**2)**0.5)
			xtt(i) = (-7.0*x1(i,1) + 8.0*x1(i,2) - x1(i,3))/2.0
     &					- 3.0*xeth(i,1)
			ytt(i) = (-7.0*y1(i,1) + 8.0*y1(i,2) - y1(i,3))/2.0
     &					- 3.0*yeth(i,1)
		end do 
		xst(1)	= (xeth(2,1) - xeth(NX-1,1))/2.0
		yst(1)	= (yeth(2,1) - yeth(NX-1,1))/2.0
		jac(1,1) = 1.0/(xsi(1,1)*yeth(1,1) - xeth(1,1)*ysi(1,1))
		xst(NX) = xst(1)
		yst(NX) = yst(1)
		jac(NX,1) = jac(1,1)
		do i = 2,NX-1
			xst(i)	= (xeth(i+1,1) - xeth(i-1,1))/2.0
			yst(i)	= (yeth(i+1,1) - yeth(i-1,1))/2.0
			jac(i,1) = 1.0/(xsi(i,1)*yeth(i,1) - xeth(i,1)*ysi(i,1))
		end do
		
		do i = 1,NX
			alpha(i,1) = xeth(i,1)**2 + yeth(i,1)**2
			betha(i,1) = xsi(i,1)*xeth(i,1) + ysi(i,1)*yeth(i,1)
			gamma(i,1) = xsi(i,1)**2 + ysi(i,1)**2
		end do  

		do i = 1,NX
			r1 = -1.0*(jac(i,1)**2)*(alpha(i,1)*xss(i)
     &				- 2.0*betha(i,1)*xst(i) + gamma(i,1)*xtt(i))
			r2 = -1.0*(jac(i,1)**2)*(alpha(i,1)*yss(i)
     &				- 2.0*betha(i,1)*yst(i) + gamma(i,1)*ytt(i))
			p1(i) = wp*(jac(i,1)*(yeth(i,1)*r1 - xeth(i,1)*r2))
     &		        + (1-wp)*p1(i)
			q1(i) = wq*(jac(i,1)*(-1.0*ysi(i,1)*r1 + xsi(i,1)*r2))
     &				+ (1-wq)*q1(i)
		end do

		do j = 2,(NY-1)
			xsi(1,j) = (x1(2,j) - x1(NX-1,j))/2.0		
			xsi(NX,j) = xsi(1,j)
      		ysi(1,j) = (y1(2,j) - y1(NX-1,j))/2.0
			ysi(NX,j) = ysi(1,j)
			xeth(1,j) = (x1(1,j+1) - x1(1,j-1))/2.0
c			xeth(1,j) = xeth(1,1)
			xeth(NX,j) = xeth(1,j)
			yeth(1,j) = (y1(1,j+1) - y1(1,j-1))/2.0
c			yeth(1,j) = yeth(1,1)
			yeth(NX,j) = yeth(1,j)
			jac(1,j) = 1.0/(xsi(1,j)*yeth(1,j) - xeth(1,j)*ysi(1,j))
			jac(NX,j) = jac(1,j)
			alpha(1,j) = xeth(1,j)**2 + yeth(1,j)**2
			alpha(NX,j) = alpha(1,j)
			betha(1,j) = xsi(1,j)*xeth(1,j) + ysi(1,j)*yeth(1,j)
			betha(NX,j) = betha(1,j)
			gamma(1,j) = xsi(1,j)**2 + ysi(1,j)**2
			gamma(NX,j) = gamma(1,j)
		end do 
		do i = 2,(NX-1)
			do j = 2,(NY-1)
				xsi(i,j) = (x1(i+1,j) - x1(i-1,j))/2.0
				ysi(i,j) = (y1(i+1,j) - y1(i-1,j))/2.0
				xeth(i,j) = (x1(i,j+1) - x1(i,j-1))/2.0
				yeth(i,j) = (y1(i,j+1) - y1(i,j-1))/2.0
				jac(i,j) = 1.0/(xsi(i,j)*yeth(i,j)
     &        			 - xeth(i,j)*ysi(i,j))
				alpha(i,j) = xeth(i,j)**2 + yeth(i,j)**2
				betha(i,j) = xsi(i,j)*xeth(i,j) + ysi(i,j)*yeth(i,j)
				gamma(i,j) = xsi(i,j)**2 + ysi(i,j)**2
			end do 
		end do 
		do j = 2,NY-1
		    if (p1(1)>0) then 
	 			xsip(1,j) = x1(2,j) - x1(1,j) 
				ysip(1,j) = y1(2,j) - y1(1,j)
	  		else 
		 		xsip(1,j) = x1(1,j) - x1(NX-1,j)
				ysip(1,j) = y1(1,j) - y1(NX-1,j)
	   		end if 
			if (q1(1)>0) then 
		 		xethp(1,j) = x1(1,j+1) - x1(1,j)  
				yethp(1,j) = y1(1,j+1) - y1(1,j)
	  		else 
		 		xethp(1,j) = x1(1,j) - x1(1,j-1)  
				yethp(1,j) = y1(1,j) - y1(1,j-1)
	   		end if 
			xsip(NX,j) = xsip(1,j)
			ysip(NX,j) = ysip(1,j)
			xethp(NX,j) = xethp(1,j)
			yethp(NX,j) = yethp(1,j)
		end do 
		do i = 2 , (NX-1)
			do j = 2,(NY-1)
				if (p1(i)>0) then 
	 				xsip(i,j) = x1(i+1,j) - x1(i,j) 
					ysip(i,j) = y1(i+1,j) - y1(i,j)
	  			else 
	 				xsip(i,j) = x1(i,j) - x1(i-1,j)
					ysip(i,j) = y1(i,j) - y1(i-1,j)
	   			end if 
				if (q1(i)>0) then 
	 				xethp(i,j) = x1(i,j+1) - x1(i,j)  
					yethp(i,j) = y1(i,j+1) - y1(i,j)
	  			else 
	 				xethp(i,j) = x1(i,j) - x1(i,j-1)  
					yethp(i,j) = y1(i,j) - y1(i,j-1)
	   			end if 
	 		end do
		end do 
		do j = 2,(NY-1)
			  x1(1,j) = (alpha(1,j)*(x1(2,j)+x1(NX-1,j)) +
     &			 gamma(1,j)*(x1(1,j+1) + x1(1,j-1)) - 0.5*betha(1,j)*
     &		    (x1(2,j+1) - x1(NX-1,j+1) - x1(2,j-1) + 
     &            x1(NX-1,j-1)) + (jac(1,j)**(-2.0))*
     &		(p1(1)*(exp(aa*(-j+1.0)))*xsip(1,j) + q1(1)*
     &				(exp(bb*(-j+1.0)))*xethp(1,j))
     &				)/(2*(alpha(1,j) + gamma(1,j)))
	 
 	   		  y1(1,j) =	(alpha(1,j)*(y1(2,j)+y1(NX-1,j)) +
     &			  gamma(1,j)*(y1(1,j+1) + y1(1,j-1)) - 0.5*betha(1,j)*
     &		     (y1(2,j+1) - y1(NX-1,j+1) - y1(2,j-1)  +
     &				 y1(NX-1,j-1)) + (jac(1,j)**(-2.0))*
     &		(p1(1)*(exp(aa*(-j+1.0)))*ysip(1,j) + q1(1)*
     &				(exp(bb*(-j+1.0)))*yethp(1,j))
     &				)/(2*(alpha(1,j) + gamma(1,j)))
			
			   x1(NX,j) = x1(1,j)
			   y1(NX,j) = y1(1,j)	
		end do 
		do i = 2,(NX-1)
			do j = 2,(NY-1)
				x1(i,j) = (alpha(i,j)*(x1(i+1,j)+x1(i-1,j)) +
     &			 gamma(i,j)*(x1(i,j+1) + x1(i,j-1)) - 0.5*betha(i,j)*
     &		    (x1(i+1,j+1) - x1(i-1,j+1) - x1(i+1,j-1) + 
     &			 x1(i-1,j-1)) + (jac(i,j)**(-2.0))*
     &		(p1(i)*(exp(aa*(-j+1.0)))*xsip(i,j) + q1(i)*
     &				(exp(bb*(-j+1.0)))*xethp(i,j))
     &             )/(2*(alpha(i,j) + gamma(i,j)))

				y1(i,j) = (alpha(i,j)*(y1(i+1,j)+y1(i-1,j)) +
     &			 gamma(i,j)*(y1(i,j+1) + y1(i,j-1)) - 0.5*betha(i,j)*
     &		    (y1(i+1,j+1) - y1(i-1,j+1) - y1(i+1,j-1) + 
     &			y1(i-1,j-1)) + (jac(i,j)**(-2.0))*
     &		(p1(i)*(exp(aa*(-j+1.0)))*ysip(i,j) + q1(i)*
     &				(exp(bb*(-j+1.0)))*yethp(i,j))
     &             )/(2*(alpha(i,j) + gamma(i,j)))
			end do 
		end do
		x1 = OMEGA * x1 + (1-OMEGA)*x
		y1 = OMEGA * y1 + (1-OMEGA)*y
		error = max(maxval(abs(x1 - x)),maxval(abs(y1 - y))) 
		write(*,*) error
		x = x1
		y = y1
		count = count + 1
	end do 
	    
	end subroutine orthogonal 
!****************************************************************
      Subroutine Output

c      Printing the result

	 Open(1,file='GRID.plt')
	 Write(1,*) 'Variable = "X","Y"'
	 Write(1,5) NX,NY
5      FORMAT('  ZONE ','   i=',I5,'   j=',I5)
	 
	 DO i=1,NX
	 DO j=1,NY
	 Write(1,*) X(I,J),Y(I,J)
	 END DO 
	 END DO 

	End Subroutine
!****************************************************************
	End Program GRID
