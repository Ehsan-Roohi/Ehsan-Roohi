!****************************************************************************
!
!  PROGRAM: AIRFOIL
!
!
!****************************************************************************

	program AIRFOIL	
	   
	
	PARAMETER (NN=1548)
	PARAMETER (NE=2421)
	PARAMETER (NEBC=5000)
	PARAMETER (MM=4)   !R-K INTEGRATION STEPS

	implicit DoublePrecision  (A-H,O-Z)

	! Variables

	Integer     :: NNode,NBound,Neib(NE,4),NW(NEBC),NM1(NE)

	DoublePrecision ::AK(4,4),AR(NEBC),A(NEBC,NEBC),R,EPS,CINP
     &,RHOOZ(NEBC),FUOZ(NEBC),FVOZ(NEBC),FEOZ(NEBC),CPP(NEBC),SO(NEBC)
     &,ALF

	DoublePrecision :: Xnode(NEBC),Ynode(NEBC),H(NEBC),C(NEBC),U(NEBC)
     &,V(NEBC),RHO(NEBC),RHOO(NEBC),P(NEBC),PO(NEBC),ULO(NEBC,4)
     &,VLO(NEBC,4),QLO(NEBC,4),RLO(NEBC,4),E11,E2,E3,E4,TY(NN),CC(NEBC)

	DoublePrecision::UO(NEBC),VO(NEBC),FU(NEBC),FV(NEBC),FUO(NEBC)
     &,FVO(NEBC) 
	
	DoublePrecision :: E(NEBC),EF(NEBC),EO(NEBC),HO(NEBC),FE(NEBC)	
     &,FEO(NEBC),TOLD(NEBC),TT(NE),PP(NE),RR(NE),UU(NE),VV(NE)
	     
	
	DoublePrecision::RESC1(NEBC),RESE1(NEBC),RESMR1(NEBC),RESMT1(NEBC)
	DoublePrecision::RESC2(NEBC),RESE2(NEBC),RESMR2(NEBC),RESMT2(NEBC)
      DoublePrecision ::RESC(NEBC),RESE(NEBC),RESMR(NEBC),RESMT(NEBC)
		
	DoublePrecision :: DUTDT1(NEBC),DUTDT2(NEBC),DURDT1(NEBC),
     &UR1(NEBC),UR2(NEBC),UT1(NEBC),UT2(NEBC),FMTV(NEBC),FMRV(NEBC)
     &,TIME,DF(NEBC),DFP(NEBC),DFM(NEBC),DRP(NEBC),DRM(NEBC),DEP(NEBC)
     &,DEM(NEBC),DURDT2(NEBC),FEV(NEBC)	

      DoublePrecision :: FINVT(NEBC),FINVR(NEBC),FINVE(NEBC),SCI(NEBC)
     &,SMTI(NEBC),SEI(NEBC),SEV(NEBC),TR(NEBC),TL(NEBC),XK(NEBC)

      DoublePrecision :: XM(NEBC),URHAT(NEBC),CO(NEBC),DT,T(NEBC)

      DoublePrecision :: UI,TI,URHAT1,UTHAT1,URHAT2,EADD	
     &,UTHAT2,GAMA,UHAT1,UHAT2,RE,CP,PR,RHOI,PIN,PRA,CV
	
	DoublePrecision :: TFHAT1,TFHAT2,VFHAT1,VFHAT2,UFHAT1,UFHAT2,	 
     &RESET(NEBC),RESMRT(NEBC),RESMTT(NEBC),SUM1(NEBC),SUM2(NEBC)
     &,DTDT1(NEBC),DTDT2(NEBC),XMIN,SUM4(NEBC),HLO(NEBC,4),SUM3(NEBC)

      DoublePrecision :: DNADX(4),DNADY(4),RHOLO(NEBC,4),PLO(NEBC,4)

	DoublePrecision ::AREA(NEBC),DS(NEBC,4),ALFA(10)

      DoublePrecision ::DX(NEBC,4),DY(NEBC,4),VI

      DoublePrecision ::X_CELL(NEBC),Y_CELL(NEBC)

	Integer :: NodeNum(NEBC,4),BOUNDNum(NEBC),Bound(NEBC,4),MSUM
     
     &	,NODE_VICINITY(NEBC,50)

	Integer :: JJ,JJJ,I,J,N1,N2,N3,N4,K,N11,N22,N33,METHOD,BCTYPE,ZERO


	GOTO 122
	CALL SYSTEM ("del *.obj ")
      CALL SYSTEM ("del *.plg ")
	CALL SYSTEM ("del *.opt ")
	CALL SYSTEM ("del *.dsp ")
	CALL SYSTEM ("del *.dsw ")
	CALL SYSTEM ("del *.txt ")

	! Body of AIRFOIL
	
122		CALL INPUT
		CALL MESHInput
	    CALL Output1
		CALL SetNeighbours
		CALL BCNODES	   
		CALL Initialize			    
10        CALL BC		
          CALL FLUX
		CALL UPDATE		
		CALL CALC
	    CALL BC
	    CALL CONVERGENCE
		
		DO  I=1,10
C		IF (MC.EQ.2000*I) CALL OUTPUT
C		IF (MC.EQ.4000*(I-1)) CALL POST_PROCESSING
	    END DO

		IF (RESIDUAL.GT.EPS) GOTO 10		
20	    CALL Output
		CALL POST_PROCESSING

	contains
!****************************************************************************
	subroutine Input
      
	    GAMA=1.4D0

	    R=287.D0

	    ZERO=0

	    ALFA(2)=0.0485D0 !0.0033

	    ALFA(3)=0.1120D0 !0.2069

	    ALFA(4)=0.1994D0 !0.4265

	    ALFA(5)=0.3285D0 !1.

	    ALFA(6)=0.5477D0

          ALFA(7)=1.  

	    CP=1004.5D0

	    CV=717.5D0

	    ALF=10.*3.14159/180

	    UI=1.*cos(ALF)

	    VI=1.*sin(ALF)
		
		RHOI=1.D0	

	    TIME=0

		XMIN=0.8D0

	    CINP=SQRT(UI**2+VI**2)/XMIN 
		
		PIN=RHOI*CINP**2/GAMA

	    HIN=0.5D0+1./((GAMA-1)*XMIN**2)

	    EIN=HIN-CINP**2/GAMA
		
C		TI=PIN/(R*RHOI)
		
		TI=GAMA*XMIN**2*PIN/RHOI	      	    		       	      

c	    DT=1.D-5

	    EPS=1.D-4
		
		RE=500.	

	END SUBROUTINE
!****************************************************************************
	Subroutine MESHINPUT
	
		Integer :: I,J,JJ,II,III		

		Open(1,file='Node.DAT')
		Open(102,FILE='RES.DAT')
	    
		DO I=1,NN         
	    READ(1,*)  J,XNode(I),YNode(I)
	    END DO

	    Open(2,file='connectivity1.DAT')		
		
	    DO I=1,550	   !RECTANGULAR ELEMENTS
	    READ(2,*)  JJ,II,III,NodeNum(I,1),NodeNum(I,2),NodeNum(I,3)
     &	,NodeNum(I,4)

	    END DO

	    CLOSE (2)

C	    GOTO 101

	    Open(3,file='connectivity2.DAT')		
		
		DO I=551,NE	!TRIANGULAR ELEMENTS
		READ(3,*)  J,II,III,NodeNum(I,1),NodeNum(I,2),NodeNum(I,3)
          
		NodeNum(I,4)=0

		END DO

	    CLOSE (3)


101	    DO I=1,NE

		N1=NodeNum(I,1)
		N2=NodeNum(I,2)
		N3=NodeNum(I,3)
          N4=NodeNum(I,4)

	    IF(N4.EQ.0)THEN

	    AREA(I)=ABS(XNode(N2)*YNode(N3)-XNode(N3)*YNode(N2)+
     &	XNode(N3)*YNode(N1)-XNode(N1)*YNode(N3)+
     &    XNode(N1)*YNode(N2)-XNode(N2)*YNode(N1))/2.

	    ELSE IF(N4.NE.0)THEN

          AREA(I)=1/2.*ABS(
     &    (XNODE(N2)-XNODE(N4))*(YNODE(N3)-YNODE(N1))-
     &    (XNODE(N3)-XNODE(N1))*(YNODE(N2)-YNODE(N4)))

	    END IF
	    
		DX(I,1)=-(XNode(N1)-XNode(N2))
          DY(I,1)=-(YNode(N1)-YNode(N2))
	    DS(I,1)=SQRT(DX(I,1)**2+DY(I,1)**2)
		   
	    DX(I,2)=-(XNode(N2)-XNode(N3))
          DY(I,2)=-(YNode(N2)-YNode(N3))
		DS(I,2)=SQRT(DX(I,2)**2+DY(I,2)**2)            

	    IF(N4.EQ.0)THEN
		DX(I,3)=-(XNode(N3)-XNode(N1))	    
	    DY(I,3)=-(YNode(N3)-YNode(N1))
	    DS(I,3)=SQRT(DX(I,3)**2+DY(I,3)**2)

	    ELSE IF(N4.NE.0)THEN
          DX(I,3)=-(XNODE(N3)-XNODE(N4))
          DY(I,3)=-(YNODE(N3)-YNODE(N4))
          DS(I,3)=SQRT((DX(I,3))**2+(DY(I,3))**2)

      	DX(I,4)=-(XNODE(N4)-XNODE(N1))
          DY(I,4)=-(YNODE(N4)-YNODE(N1))
          DS(I,4)=SQRT((DX(I,4))**2+(DY(I,4))**2)

          END IF

	    IF(N4.EQ.0)THEN
          X_CELL(I)=(XNODE(N1)+XNODE(N2)+XNODE(N3))/3.0
          Y_CELL(I)=(YNODE(N1)+YNODE(N2)+YNODE(N3))/3.0
          
		ELSE IF(N4.NE.0)THEN
          X_CELL(I)=(XNODE(N1)+XNODE(N2)+XNODE(N3)+XNODE(N4))/4.0
          Y_CELL(I)=(YNODE(N1)+YNODE(N2)+YNODE(N3)+XNODE(N4))/4.0
          END IF

	    END DO

	   !FINDING THE CELLS AROUND A NODE

	    MSUM=0

          DO I=1,NN
          NUM=0

          DO J=1,NE

          IF((NodeNum(J,1).EQ.I).OR.(NodeNum(J,2).EQ.I).OR.
     &       (NodeNum(J,3).EQ.I).OR.(NodeNum(J,4).EQ.I)) THEN

          NUM=NUM+1

          END IF
          END DO

          IF(NUM>MSUM)THEN
          MSUM=NUM
          END IF
          END DO

         
		DO I=1,NN
          K=0
          DO J=1,NE
          IF((NodeNum(J,1).EQ.I).OR.(NodeNum(J,2).EQ.I).OR.
     &	   (NodeNum(J,3).EQ.I).OR.(NodeNum(J,4).EQ.I)) THEN
          K=K+1
          NODE_VICINITY(I,K)=J
          END IF
          END DO
          END DO
	    
	End subroutine
!****************************************************************************
	Subroutine SetNeighbours

		    INTEGER :: I,J,N1,N2,N3,N4,M1,M2,M3,M4

		    Neib=0
		    Do I=1,NE

			N1=NodeNum(I,1)
			N2=NodeNum(I,2)
			N3=NodeNum(I,3)	
			N4=NodeNum(I,4)

			Do J=1,NE

		    M1=NodeNum(J,1)
			M2=NodeNum(J,2)
			M3=NodeNum(J,3)	
			M4=NodeNum(J,4)

	 !FOR TRINGULAR CELLS
		IF((J.NE.I).AND.(N4.EQ.0))THEN
		IF((N2.EQ.M1).OR.(N2.EQ.M2).OR.
     &	(N2.EQ.M3).OR.(N2.EQ.M4))THEN
		IF((N1.EQ.M1).OR.(N1.EQ.M2).OR.
     &	(N1.EQ.M3).OR.(N1.EQ.M4))THEN
		NEIB(I,1)=J
		ENDIF
		ENDIF

		IF((N2.EQ.M1).OR.(N2.EQ.M2).OR.
     &	(N2.EQ.M3).OR.(N2.EQ.M4))THEN
		IF((N3.EQ.M1).OR.(N3.EQ.M2).OR.
     &	(N3.EQ.M3).OR.(N3.EQ.M4))THEN
		NEIB(I,2)=J
		ENDIF
		ENDIF

		IF((N1.EQ.M1).OR.(N1.EQ.M2).OR.
     &	(N1.EQ.M3).OR.(N1.EQ.M4))THEN
		IF((N3.EQ.M1).OR.(N3.EQ.M2).OR.
     &	(N3.EQ.M3).OR.(N3.EQ.M4))THEN
		NEIB(I,3)=J
		ENDIF
		ENDIF

		NEIB(I,4)=-2

	! FOR QUAD CELLS
		ELSEIF((J.NE.I).AND.(N4.NE.0))THEN	        
		

            IF((N2.EQ.M1).OR.(N2.EQ.M2).OR.
     &		  (N2.EQ.M3).OR.(N2.EQ.M4))THEN
            IF((N1.EQ.M1).OR.(N1.EQ.M2).OR.
     &		  (N1.EQ.M3).OR.(N1.EQ.M4))THEN
            NEIB(I,1)=J
            ENDIF
            ENDIF

            IF((N2.EQ.M1).OR.(N2.EQ.M2).
     &		  OR.(N2.EQ.M3).OR.(N2.EQ.M4))THEN
            IF((N3.EQ.M1).OR.(N3.EQ.M2).
     &		  OR.(N3.EQ.M3).OR.(N3.EQ.M4))THEN
            NEIB(I,2)=J
            ENDIF
            ENDIF

            IF((N4.EQ.M1).OR.(N4.EQ.M2).
     & 		  OR.(N4.EQ.M3).OR.(N4.EQ.M4))THEN
            IF((N3.EQ.M1).OR.(N3.EQ.M2).
     &		  OR.(N3.EQ.M3).OR.(N3.EQ.M4))THEN
            NEIB(I,3)=J
            ENDIF
            ENDIF

            IF((N4.EQ.M1).OR.(N4.EQ.M2).
     &		  OR.(N4.EQ.M3).OR.(N4.EQ.M4))THEN
            IF((N1.EQ.M1).OR.(N1.EQ.M2).
     &		  OR.(N1.EQ.M3).OR.(N1.EQ.M4))THEN
            NEIB(I,4)=J
            ENDIF
            ENDIF

            ENDIF           

144	      CONTINUE
            END DO
	      END DO


		DO I=1,NE
		IF (NEIB(I,4).EQ.-2) NEIB(I,4)=-2 		
		END DO

		
		DO I=1,NE
		DO J=1,3

		IF (NEIB(I,J).EQ.0) THEN

		DO K=1,101
C     IF (I.EQ.WALL(K,1)) NEIB(I,J)=-1;
		END DO

		END IF
      
		END DO
		END DO


		OPEN(23,FILE='NEIB.DAT')
		OPEN(22,FILE='NEIGHBOUR.DAT')
		DO I=1,NE
		WRITE(22,*)I,NEIB(I,1),NEIB(I,2),NEIB(I,3),NEIB(I,4)
C	WRITE(23,*)I,NEIB(I,1),NEIB(I,2),NEIB(I,3),NEIB(I,4)
		ENDDO
		CLOSE(22)

	end subroutine 
!****************************************************************************
	Subroutine BCNODES

	    ! RECOGNIZE B.C CELLES

          Open(9,file='BC.DAT')
		Open(30,file='BC1.DAT')  !FAR FEILD 
		Open(31,file='BC2.DAT')  !UPPER WALL  
	    Open(32,file='BC3.DAT')  !LOWER WALL
	    OPEN(15,file='BE.DAT')

	
		!FAR FIELD

		I=1

		DO WHILE(I<500)

  		READ(30,*,END=101) BOUNDNum(I)

		I=I+1

		END DO

          NW=0 

101		NBOUND=I-1 	

		! ASSIGNE B.C. NODES TO THEIR ELEMENTS    
	    
	    DO J=1,NBOUND  
	    DO I=1,NE
	    
		IF ((BOUNDNum(J).EQ.NodeNum(I,1).OR.
     &	    BOUNDNum(J).EQ.NodeNum(I,2).OR.
     &        BOUNDNum(J).EQ.NodeNum(I,3))
     &	   .AND.(BOUNDNum(J+1).EQ.NodeNum(I,1).OR. 
     &	    BOUNDNum(J+1).EQ.NodeNum(I,2).OR.
     &	    BOUNDNum(J+1).EQ.NodeNum(I,3))) THEN

	        K=I
	        BOUND(K,1)=BOUNDNum(J)
	        BOUND(K,2)=BOUNDNum(J+1)
			
			WRITE(15,*) K,BOUND(K,1),BOUND(K,2)	        
			              
			END IF   

	     END DO
	     END DO 

          !RECOGNIZE UPPER WALL BC
	 
	   	I=1

		DO while(I<500)

  		READ(31,*,END=102) BOUNDNum(I)

		I=I+1

		END DO

102		NBOUND=I-1 	

			!ASSIGNE B.C. NODES TO THEIR ELEMENTS    
	    
	    DO J=1,NBOUND  
	    DO I=1,NE
	    
			IF ((BOUNDNum(J).EQ.NodeNum(I,1).OR.
     &	    BOUNDNum(J).EQ.NodeNum(I,2).OR.
     &        BOUNDNum(J).EQ.NodeNum(I,3).OR.
     &        BOUNDNum(J).EQ.NodeNum(I,4))
     &	   .AND.(BOUNDNum(J+1).EQ.NodeNum(I,1).OR. 
     &	    BOUNDNum(J+1).EQ.NodeNum(I,2).OR.
     &	    BOUNDNum(J+1).EQ.NodeNum(I,3).OR.
     &	    BOUNDNum(J+1).EQ.NodeNum(I,4))) THEN

	        K=I
	        BOUND(K,1)=BOUNDNum(J)
	        BOUND(K,2)=BOUNDNum(J+1)

	        NW(K)=1
	        NM1(K)=1
			
			WRITE(15,*) K,BOUND(K,1),BOUND(K,2)	        
			              
			END IF   

	     END DO
	     END DO 	                             	         

          !RECOGNIZE LOWER WALL BC
	 
	   	I=1

		DO while(I<500)

  		READ(32,*,END=103) BOUNDNum(I)

		I=I+1

		END DO

103		NBOUND=I-1 	

		!ASSIGNE B.C. NODES TO THEIR ELEMENTS    
	    
	    DO J=1,NBOUND  
	    DO I=1,NE
	    
			IF ((BOUNDNum(J).EQ.NodeNum(I,1).OR.
     &	    BOUNDNum(J).EQ.NodeNum(I,2).OR.
     &        BOUNDNum(J).EQ.NodeNum(I,3).OR.
     &        BOUNDNum(J).EQ.NodeNum(I,4))
     &	   .AND.(BOUNDNum(J+1).EQ.NodeNum(I,1).OR. 
     &	    BOUNDNum(J+1).EQ.NodeNum(I,2).OR.
     &	    BOUNDNum(J+1).EQ.NodeNum(I,3).OR.
     &	    BOUNDNum(J+1).EQ.NodeNum(I,4))) THEN

	        K=I
	        BOUND(K,1)=BOUNDNum(J)
	        BOUND(K,2)=BOUNDNum(J+1)

	        NW(K)=1
	        NM1(K)=-1
			
			WRITE(15,*) K,BOUND(K,1),BOUND(K,2)	        
			              
			END IF   

	     END DO
	     END DO 	                		   	   

	End subroutine
!****************************************************************************
	Subroutine Initialize

	   	U(:)=1.*COS(ALF) 

	    V(:)=1.*SIN(ALF)  

	    RHO(:)=1.

	    C(:)=1./XMIN
		  
	    P(:)=RHO(:)*CINP**2/GAMA	   

	    H(:)=0.5D0+1./((GAMA-1)*XMIN**2)

	    E(:)=H(:)-C(:)**2/GAMA	    

	    XM(:)=SQRT(U(:)**2+V(:)**2)/C(:)

	    T(:)=TI

		UO=U	    
	    VO=V	   
	    PO=P	    
	    HO=H	    
	    RHOO=RHO  
	    CO=C	    
	    EO=E
		TOLD=T
			
	End subroutine 
!****************************************************************************
	Subroutine BC     	
		 		 
101		 UO=U	    
	     VO=V	   
	     PO=P	    
	     HO=H	    
	     RHOO=RHO  
	     CO=C	    
	     EO=E   
	        		 	
	End subroutine 			
!****************************************************************************
      Subroutine FLUX       

	EADD=EPS
	
c	DO J=1,6

	DO I=1,NE

	FUO(I)=RHOO(I)*UO(I)
	FVO(I)=RHOO(I)*VO(I)
	FEO(I)=RHOO(I)*EO(I)

	END DO

c	IF (J.EQ.1) THEN

c      DO I=1,NE

c	RHOOZ(I)=RHOO(I)
c	FVOZ(I)=FVO(I)
c	FUOZ(I)=FUO(I)
c	FEOZ(I)=FEO(I)	

c	END DO
	
c	END IF

	DO I=1,NE

	QLO(I,1)=(UO(I)*DY(I,1)-VO(I)*DX(I,1))/DS(I,1)
      QLO(I,2)=(UO(I)*DY(I,2)-VO(I)*DX(I,2))/DS(I,2)
      QLO(I,3)=(UO(I)*DY(I,3)-VO(I)*DX(I,3))/DS(I,3)
	IF (NODENUM(I,4).NE.0) QLO(I,4)=(UO(I)*DY(I,4)-VO(I)*DX(I,4))
     &                                  /DS(I,4)

	RLO(I,1)=(UO(I)*DX(I,1)+VO(I)*DY(I,1))/DS(I,1)
      RLO(I,2)=(UO(I)*DX(I,2)+VO(I)*DY(I,2))/DS(I,2)
      RLO(I,3)=(UO(I)*DX(I,3)+VO(I)*DY(I,3))/DS(I,3)
      IF (NODENUM(I,4).NE.0) RLO(I,4)=(UO(I)*DX(I,4)+VO(I)*DY(I,4))
     &                                  /DS(I,4)
	END DO	

	!!! CONTINUITY  
		
	DO I=1,NE

	NE1=NEIB(I,1)
	NE2=NEIB(I,2)
	NE3=NEIB(I,3)
	NE4=NEIB(I,4)

c      NF=4

c	IF (NEIB(I,4).EQ.-2) NF=3

	IF (NE1.EQ.ZERO)  THEN 

	L=1

	NE1=NE1+NEBC

	DX(NE1,L)=DX(I,L)

	DY(NE1,L)=DY(I,L)

	DS(NE1,L)=DS(I,L)  
	
	IF (NW(I).NE.1) THEN    

      QLO(NE1,L)=(UI*DY(I,L)-VI*DX(I,L))/DS(I,L) 

	RLO(NE1,L)=(UI*DX(I,L)+VI*DY(I,L))/DS(I,L) 

	UO(NE1)=UI

	VO(NE1)=VI
      
	RHOO(NE1)=RHOI  
	
	HO(NE1)=0.5D0*SQRT(UI**2+VI**2)+(CINP**2)/(GAMA-1)
	
	EO(NE1)=HO(NE1)-CINP**2/GAMA  
	
	PO(NE1)=RHOI*CINP**2/GAMA 

	TOLD(NE1)=TI
	
	END IF    
	
	IF (NW(I).EQ.1) THEN
	
	QLO(NE1,L)=-QLO(I,L)

	RLO(NE1,L)=-RLO(I,L)

	UO(NE1)=DS(NE1,L)*(QLO(NE1,L)*DY(NE1,L)+RLO(NE1,L)*DX(NE1,L))/

     &	(DY(NE1,L)**2+DX(NE1,L)**2)

	VO(NE1)=DS(NE1,L)*(-QLO(NE1,L)*DX(NE1,L)+RLO(NE1,L)*DY(NE1,L))

     &	/(DY(NE1,L)**2+DX(NE1,L)**2)	
	
	RHOO(NE1)=RHOO(I)  
	
	HO(NE1)=HO(I)
	
	EO(NE1)=EO(I) 
	
	PO(NE1)=PO(I)

	TOLD(NE1)=TOLD(I)
		
	
	END IF  
	
	END IF  

	IF (NE2.EQ.ZERO) THEN

      L=2 
	
	NE2=NE2+NEBC

	DX(NE2,L)=DX(I,L)

	DY(NE2,L)=DY(I,L)

	DS(NE2,L)=DS(I,L)

	IF (NW(I).NE.1) THEN 

      QLO(NE2,L)=(UI*DY(I,L)-VI*DX(I,L))/DS(I,L) 

      RLO(NE2,L)=(UI*DX(I,L)+VI*DY(I,L))/DS(I,L)

	UO(NE2)=UI		

	VO(NE2)=VI
     
	RHOO(NE2)=RHOI

	HO(NE2)=0.5D0*SQRT(UI**2+VI**2)+(CINP**2)/(GAMA-1)

	EO(NE2)=HO(NE2)-CINP**2/GAMA  

	PO(NE2)=RHOI*CINP**2/GAMA

	TOLD(NE2)=TI
	
	END IF
	
	IF (NW(I).EQ.1) THEN
      
	RLO(NE2,L)=-RLO(I,L)

	QLO(NE2,L)=-QLO(I,L)

	UO(NE2)=DS(NE2,L)*((QLO(NE2,L)*DY(NE2,L))+(RLO(NE2,L)*DX(NE2,L)))

     &	/(DY(NE2,L)**2+DX(NE2,L)**2)

	VO(NE2)=DS(NE2,L)*((-QLO(NE2,L)*DX(NE2,L))+(RLO(NE2,L)*DY(NE2,L)))

     &	/(DY(NE2,L)**2+DX(NE2,L)**2)
	
	
	RHOO(NE2)=RHOO(I)  
	
	HO(NE2)=HO(I)
	
	EO(NE2)=EO(I)  
	
	PO(NE2)=PO(I)

	TOLD(NE2)=TOLD(I)	

	END IF
      
	END IF

	IF (NE3.EQ.ZERO)     THEN

      L=3
	
	NE3=NE3+NEBC

	DX(NE3,L)=DX(I,L)

	DY(NE3,L)=DY(I,L)

	DS(NE3,L)=DS(I,L)
	
	IF (NW(I).NE.1) THEN 
	
	QLO(NE3,L)=(UI*DY(I,L)-VI*DX(I,L))/DS(I,L) 	

	RLO(NE3,L)=(UI*DX(I,L)+VI*DY(I,L))/DS(I,L)

	UO(NE3)=UI

	VO(NE3)=VI

	RHOO(NE3)=RHOI
	
	HO(NE3)=0.5D0*SQRT(UI**2+VI**2)+(CINP**2)/(GAMA-1)

	EO(NE3)=HO(NE3)-CINP**2/GAMA  
	
	PO(NE3)=RHOI*CINP**2/GAMA

	TOLD(NE3)=TI

	END IF	

	IF (NW(I).EQ.1) THEN	
	
	RLO(NE3,L)=-RLO(I,L)

      QLO(NE3,L)=-QLO(I,L)	
	
	UO(NE3)=DS(NE3,L)*((QLO(NE3,L)*DY(NE3,L))+(RLO(NE3,L)*DX(NE3,L)))

     &	/(DY(NE3,L)**2+DX(NE3,L)**2)

	VO(NE3)=DS(NE3,L)*((-QLO(NE3,L)*DX(NE3,L))+(RLO(NE3,L)*DY(NE3,L)))

     &	/(DY(NE3,L)**2+DX(NE3,L)**2)
     
      RHOO(NE3)=RHOO(I)  
	
	HO(NE3)=HO(I)
	
	EO(NE3)=EO(I)  
	
	PO(NE3)=PO(I)

	TOLD(NE3)=TOLD(I)
     
	END IF

	END IF

	IF (NE4.NE.-2) THEN

	IF (NE4.EQ.ZERO) THEN

      L=4 
	
	NE4=NE4+NEBC

	DX(NE4,L)=DX(I,L)

	DY(NE4,L)=DY(I,L)

	DS(NE4,L)=DS(I,L)
	
	IF (NW(I).NE.1) THEN  
	
	QLO(NE4,L)=(UI*DY(I,L)-VI*DX(I,L))/DS(I,L) 	

	RLO(NE4,L)=(UI*DX(I,L)+VI*DY(I,L))/DS(I,L)

	UO(NE4)=UI

	VO(NE4)=VI

	RHOO(NE4)=RHOI
	
	HO(NE4)=0.5D0*SQRT(UI**2+VI**2)+(CINP**2)/(GAMA-1)

	EO(NE4)=HO(NE4)-CINP**2/GAMA  
	
	PO(NE4)=RHOI*CINP**2/GAMA

	TOLD(NE4)=TI

	END IF	

	IF (NW(I).EQ.1) THEN	
	
	RLO(NE4,L)=-RLO(I,L)

      QLO(NE4,L)=-QLO(I,L)	
	
	UO(NE4)=DS(NE4,L)*((QLO(NE4,L)*DY(NE4,L))+(RLO(NE4,L)*DX(NE4,L)))

     &	/(DY(NE4,L)**2+DX(NE4,L)**2)

	VO(NE4)=DS(NE4,L)*((-QLO(NE4,L)*DX(NE4,L))+(RLO(NE4,L)*DY(NE4,L)))

     &	/(DY(NE4,L)**2+DX(NE4,L)**2)
     
      RHOO(NE4)=RHOO(I)  
	
	HO(NE4)=HO(I)
	
	EO(NE4)=EO(I)  
	
	PO(NE4)=PO(I)
	
	TOLD(NE4)=TOLD(I)
     
	END IF

	END IF

	END IF

	! ROE'S AVERAGING

      RHOFH1=SQRT(RHOO(I)*RHOO(NE1))
	RHOS1=(SQRT(RHOO(I))+SQRT(RHOO(NE1)))
	UFHAT1=(UO(I)*SQRT(RHOO(I))+UO(NE1)*SQRT(RHOO(NE1)))/RHOS1
	VFHAT1=(VO(I)*SQRT(RHOO(I))+VO(NE1)*SQRT(RHOO(NE1)))/RHOS1
	HFHAT1=(HO(I)*SQRT(RHOO(I))+HO(NE1)*SQRT(RHOO(NE1)))/RHOS1
	CHAT1=SQRT((GAMA-1)*(HFHAT1-0.5D0*(UFHAT1**2+VFHAT1**2)))
	QFHAT1=(UFHAT1*DY(I,1)-VFHAT1*DX(I,1))/DS(I,1)
	RFHAT1=(UFHAT1*DX(I,1)+VFHAT1*DY(I,1))/DS(I,1)
	

	RHOFH2=SQRT(RHOO(I)*RHOO(NE2))
	RHOS2=(SQRT(RHOO(I))+SQRT(RHOO(NE2)))
	UFHAT2=(UO(I)*SQRT(RHOO(I))+UO(NE2)*SQRT(RHOO(NE2)))/RHOS2
	VFHAT2=(VO(I)*SQRT(RHOO(I))+VO(NE2)*SQRT(RHOO(NE2)))/RHOS2
	HFHAT2=(HO(I)*SQRT(RHOO(I))+HO(NE2)*SQRT(RHOO(NE2)))/RHOS2
	CHAT2=SQRT((GAMA-1)*(HFHAT2-0.5D0*(UFHAT2**2+VFHAT2**2)))
	QFHAT2=(UFHAT2*DY(I,2)-VFHAT2*DX(I,2))/DS(I,2)
	RFHAT2=(UFHAT2*DX(I,2)+VFHAT2*DY(I,2))/DS(I,2)


	RHOFH3=SQRT(RHOO(I)*RHOO(NE3))
	RHOS3=(SQRT(RHOO(I))+SQRT(RHOO(NE3)))
      UFHAT3=(UO(I)*SQRT(RHOO(I))+UO(NE3)*SQRT(RHOO(NE3)))/RHOS3
	VFHAT3=(VO(I)*SQRT(RHOO(I))+VO(NE3)*SQRT(RHOO(NE3)))/RHOS3 
	HFHAT3=(HO(I)*SQRT(RHOO(I))+HO(NE3)*SQRT(RHOO(NE3)))/RHOS3
      CHAT3=SQRT((GAMA-1)*(HFHAT3-0.5D0*(UFHAT3**2+VFHAT3**2))) 	
	QFHAT3=(UFHAT3*DY(I,3)-VFHAT3*DX(I,3))/DS(I,3)
	RFHAT3=(UFHAT3*DX(I,3)+VFHAT3*DY(I,3))/DS(I,3)
      
	IF (NEIB(I,4).NE.-2) THEN

	RHOFH4=SQRT(RHOO(I)*RHOO(NE4))
	RHOS4=(SQRT(RHOO(I))+SQRT(RHOO(NE4)))
      UFHAT4=(UO(I)*SQRT(RHOO(I))+UO(NE4)*SQRT(RHOO(NE4)))/RHOS4
	VFHAT4=(VO(I)*SQRT(RHOO(I))+VO(NE4)*SQRT(RHOO(NE4)))/RHOS4 
	HFHAT4=(HO(I)*SQRT(RHOO(I))+HO(NE4)*SQRT(RHOO(NE4)))/RHOS4
      CHAT4=SQRT((GAMA-1)*(HFHAT4-0.5D0*(UFHAT4**2+VFHAT4**2))) 	
	QFHAT4=(UFHAT4*DY(I,4)-VFHAT4*DX(I,4))/DS(I,4)
	RFHAT4=(UFHAT4*DX(I,4)+VFHAT4*DY(I,4))/DS(I,4)

	END IF
      
      !1st SIDE	

	MF=1

	DELX=DX(I,MF)/DS(I,MF);DELY=DY(I,MF)/DS(I,MF)

	IF (NE1.GE.NEBC) THEN  ! FACE ADJACENT TO THE BOUNDARY

	L1=1;M1=L1;K1=1
	
	GOTO 12

	END IF

	DO K=1,4

	IF (NEIB(NE1,K).EQ.I) L1=K

	IF (NEIB(I,K).EQ.NE1) M1=K

	END DO

	K1=-1

12	DP1=PO(NE1)-PO(I)

	DQ1=K1*QLO(NE1,L1)-QLO(I,M1)

	DRHO1=RHOO(NE1)-RHOO(I)

	DR1=K1*RLO(NE1,L1)-RLO(I,M1)

	
	V1=(DP1-RHOFH1*CHAT1*DQ1)/(2*CHAT1**2)

	COF1=ABS(QFHAT1-CHAT1)*V1


	V2=RHOFH1*DR1/CHAT1

	COF2=ABS(QFHAT1)*V2


	V3=DRHO1-DP1/CHAT1**2

	COF3=ABS(QFHAT1)*V3


	V4=(DP1+RHOFH1*CHAT1*DQ1)/(2*CHAT1**2)

	COF4=ABS(QFHAT1+CHAT1)*V4

	CALL VISCOUS(I,1,TAU_XX,TAU_XY,TAU_YY,Q_X,Q_Y)

C	IF (I.EQ.51) WRITE(*,*)	TAU_XX,TAU_XY,TAU_YY,Q_X,Q_Y
	
	FVT1=(-TAU_XX*DELY+TAU_XY*DELX)*DS(I,1)		

	FVV1=(-TAU_XY*DELY+TAU_YY*DELX)*DS(I,1)	

	FVE1=((-UFHAT1*TAU_XX-VFHAT1*TAU_XY+Q_X)*DELY-
     &	 (-UFHAT1*TAU_XY-VFHAT1*TAU_YY+Q_Y)*DELX)*DS(I,1)	

      
	!2nd SIDE

	MF=2

	DELX=DX(I,MF)/DS(I,MF);DELY=DY(I,MF)/DS(I,MF)

	IF (NE2.GE.NEBC) THEN

	L2=2;M2=L2;K2=1

	GOTO 14

	END IF

	DO K=1,4

	IF (NEIB(NE2,K).EQ.I) L2=K

	IF (NEIB(I,K).EQ.NE2) M2=K

	END DO

	K2=-1

14    DP2=PO(NE2)-PO(I)
	
	DQ2=K2*QLO(NE2,L2)-QLO(I,M2)

	DRHO2=RHOO(NE2)-RHOO(I)

	DR2=K2*RLO(NE2,L2)-RLO(I,M2)


	V5=(DP2-RHOFH2*CHAT2*DQ2)/(2*CHAT2**2)

	COF5=ABS(QFHAT2-CHAT2)*V5

      V6=RHOFH2*DR2/CHAT2  

	COF6=ABS(QFHAT2)*V6

	V7=DRHO2-DP2/CHAT2**2
	
	COF7=ABS(QFHAT2)*V7

	V8=(DP2+RHOFH2*CHAT2*DQ2)/(2*CHAT2**2)

	COF8=ABS(QFHAT2+CHAT2)*V8

	CALL VISCOUS(I,2,TAU_XX,TAU_XY,TAU_YY,Q_X,Q_Y)

	FVT2=(-TAU_XX*DELY+TAU_XY*DELX)*DS(I,2)		

	FVV2=(-TAU_XY*DELY+TAU_YY*DELX)*DS(I,2)	

	FVE2=((-UFHAT2*TAU_XX-VFHAT2*TAU_XY+Q_X)*DELY-
     &	 (-UFHAT2*TAU_XY-VFHAT2*TAU_YY+Q_Y)*DELX)*DS(I,2)	
 
      !3rd SIDE

	MF=3

	DELX=DX(I,MF)/DS(I,MF);DELY=DY(I,MF)/DS(I,MF)

	IF (NE3.GE.NEBC) THEN

	L3=3;M3=L3;K3=1
	
	GOTO 16

	END IF

	DO K=1,4

	IF (NEIB(NE3,K).EQ.I) L3=K

	IF (NEIB(I,K).EQ.NE3) M3=K

	END DO

	K3=-1

16    DP3=PO(NE3)-PO(I)
	
	DQ3=K3*QLO(NE3,L3)-QLO(I,M3)

	DRHO3=RHOO(NE3)-RHOO(I)

	DR3=K3*RLO(NE3,L3)-RLO(I,M3)

	V9=(DP3-RHOFH3*CHAT3*DQ3)/(2*CHAT3**2)

	COF9=ABS(QFHAT3-CHAT3)*V9

      V10=RHOFH3*DR3/CHAT3  

	COF10=ABS(QFHAT3)*V10

	V11=DRHO3-DP3/CHAT3**2
	
	COF11=ABS(QFHAT3)*V11

	V12=(DP3+RHOFH3*CHAT3*DQ3)/(2*CHAT3**2)

	COF12=ABS(QFHAT3+CHAT3)*V12

	CALL VISCOUS(I,3,TAU_XX,TAU_XY,TAU_YY,Q_X,Q_Y)

	FVT3=(-TAU_XX*DELY+TAU_XY*DELX)*DS(I,3)		

	FVV3=(-TAU_XY*DELY+TAU_YY*DELX)*DS(I,3)	

	FVE3=((-UFHAT3*TAU_XX-VFHAT3*TAU_XY+Q_X)*DELY-
     &	 (-UFHAT3*TAU_XY-VFHAT3*TAU_YY+Q_Y)*DELX)*DS(I,3)	

	!4th SIDE

	MF=4

	DELX=DX(I,MF)/DS(I,MF);DELY=DY(I,MF)/DS(I,MF)

	IF (NODENUM(I,4).NE.0) THEN

	IF (NE4.GE.NEBC) THEN

	L4=4;M4=L4;K4=1
	
	GOTO 18

	END IF

	DO K=1,4

	IF (NEIB(NE4,K).EQ.I) L4=K

	IF (NEIB(I,K).EQ.NE4) M4=K

	END DO

	K4=-1

18    DP4=PO(NE4)-PO(I)
	
	DQ4=K4*QLO(NE4,L4)-QLO(I,M4)

	DRHO4=RHOO(NE4)-RHOO(I)

	DR4=K4*RLO(NE4,L4)-RLO(I,M4)


	V13=(DP4-RHOFH4*CHAT4*DQ4)/(2*CHAT4**2)

	COF13=ABS(QFHAT4-CHAT4)*V13

      V14=RHOFH4*DR4/CHAT4  

	COF14=ABS(QFHAT4)*V14

	V15=DRHO4-DP4/CHAT4**2
	
	COF15=ABS(QFHAT4)*V15

	V16=(DP4+RHOFH4*CHAT4*DQ4)/(2*CHAT4**2)

	COF16=ABS(QFHAT4+CHAT4)*V16      
	
	CALL VISCOUS(I,4,TAU_XX,TAU_XY,TAU_YY,Q_X,Q_Y)

	FVT4=(-TAU_XX*DELY+TAU_XY*DELX)*DS(I,4)		

	FVV4=(-TAU_XY*DELY+TAU_YY*DELX)*DS(I,4)	

	FVE4=((-UFHAT4*TAU_XX-VFHAT4*TAU_XY+Q_X)*DELY-
     &	 (-UFHAT4*TAU_XY-VFHAT4*TAU_YY+Q_Y)*DELX)*DS(I,4)	

108   END IF 

      F11=COF1

	F12=0
     
	F13=COF3

	F14=COF4

      F21=COF5

	F22=0

	F23=COF7

	F24=COF8

	F31=COF9

	F32=0

	F33=COF11

	F34=COF12

	IF (NEIB(I,4).NE.-2) THEN

	FC1=COF13

	FC2=0

	FC3=COF15

	FC4=COF16

	END IF
 
      IW=0

	IF (NEIB(I,4).NE.-2) IW=1

	DEME=(ABS(QFHAT1)+CHAT1)*DS(I,1)+(ABS(QFHAT2)+CHAT2)*DS(I,2)
     &   +(ABS(QFHAT3)+CHAT3)*DS(I,3)+IW*(ABS(QFHAT4)+CHAT4)*DS(I,4)
	
	DT=0.2*AREA(I)/DEME
		
       ! ROE

	PHI1=1./2.*(RHOO(I)*QLO(I,1)+(K1)*RHOO(NE1)*QLO(NE1,L1))

     &	-1./2.*(F11+F12+F13+F14)

	PHI2=1./2.*(RHOO(I)*QLO(I,2)+(K2)*RHOO(NE2)*QLO(NE2,L2))

     &	-1./2.*(F21+F22+F23+F24)

	PHI3=1./2.*(RHOO(I)*QLO(I,3)+(K3)*RHOO(NE3)*QLO(NE3,L3))

     &	-1./2.*(F31+F32+F33+F34) 
	
	IF (NEIB(I,4).NE.-2) THEN
	
	PHI4=1./2.*(RHOO(I)*QLO(I,4)+(K4)*RHOO(NE4)*QLO(NE4,L4))

     &	-1./2.*(FC1+FC2+FC3+FC4) 

	END IF 

	DF(I)=PHI1*DS(I,1)+PHI2*DS(I,2)+PHI3*DS(I,3)+IW*PHI4*DS(I,4) 
	             
	RESC1(I)=-DF(I)/AREA(I)

	  !! UPDATE CONTINUITY

c	RHOO(I)=RHOOZ(I)+ALFA(J+1)*DT*RESC1(I)
      
	RHO(I)=RHOO(I)+DT*RESC1(I)
    !------------------------------------------------------------------------
	  !!!!! U MOMENTUN

	 ! 1st SIDE 
	
	F41=COF1*(UFHAT1-CHAT1*DY(I,1)/DS(I,1))

	F42=COF2*CHAT1*DX(I,1)/DS(I,1)

	F43=COF3*UFHAT1

	F44=COF4*(UFHAT1+CHAT1*DY(I,1)/DS(I,1))

      !2nd SIDE

      F51=COF5*(UFHAT2-CHAT2*DY(I,2)/DS(I,2))

	F52=COF6*CHAT2*DX(I,2)/DS(I,2)

	F53=COF7*UFHAT2

	F54=COF8*(UFHAT2+CHAT2*DY(I,2)/DS(I,2))

	! 3rd SIDE

      F61=COF9*(UFHAT3-CHAT3*DY(I,3)/DS(I,3))

	F62=COF10*CHAT3*DX(I,3)/DS(I,3)

	F63=COF11*UFHAT3

	F64=COF12*(UFHAT3+CHAT3*DY(I,3)/DS(I,3))

	!4TH SIDE

	IF (NEIB(I,4).NE.-2) THEN

	FU1=COF13*(UFHAT4-CHAT4*DY(I,4)/DS(I,4))

	FU2=COF14*CHAT3*DX(I,4)/DS(I,4)

	FU3=COF15*UFHAT4

	FU4=COF16*(UFHAT4+CHAT4*DY(I,4)/DS(I,4))

	END IF

    	
	SX1=1./2.*(RHOO(I)*QLO(I,1)*UO(I)+PO(I)*DY(I,1)/DS(I,1)+(K1)*
	
     &(RHOO(NE1)*QLO(NE1,L1)*UO(NE1)+PO(NE1)*DY(NE1,L1)/DS(NE1,L1)))

     &-1./2.*(F41+F42+F43+F44)

	SX2=1./2.*(RHOO(I)*QLO(I,2)*UO(I)+PO(I)*DY(I,2)/DS(I,2)+(K2)*
	
     &(RHOO(NE2)*QLO(NE2,L2)*UO(NE2)+PO(NE2)*DY(NE2,L2)/DS(NE2,L2)))

     &-1./2.*(F51+F52+F53+F54) 

	SX3=1./2.*(RHOO(I)*QLO(I,3)*UO(I)+PO(I)*DY(I,3)/DS(I,3)+(K3)*
	
     &(RHOO(NE3)*QLO(NE3,L3)*UO(NE3)+PO(NE3)*DY(NE3,L3)/DS(NE3,L3)))

     &-1./2.*(F61+F62+F63+F64)

	IF (NEIB(I,4).NE.-2) THEN

	SX4=1./2.*(RHOO(I)*QLO(I,4)*UO(I)+PO(I)*DY(I,4)/DS(I,4)+(K4)*
	
     &(RHOO(NE4)*QLO(NE4,L4)*UO(NE4)+PO(NE4)*DY(NE4,L4)/DS(NE4,L4)))

     &-1./2.*(FU1+FU2+FU3+FU4)  
     
      END IF 

	FINVT(I)=SX1*DS(I,1)+SX2*DS(I,2)+SX3*DS(I,3)+IW*SX4*DS(I,4) !INVISCID FLUX

	FVT=FVT1+FVT2+FVT3+IW*FVT4 ! VISCOUS FLUX

	! CALCULATE FLUX RESIDUAL

      RESMT1(I)=-(FINVT(I)+FVT)/AREA(I)

	!! UPDATE T-MOM

	!INVISCID
      
	FU(I)=FUO(I)+DT*RESMT1(I)

	!VISCOUS

c	FUO(I)=FUOZ(I)+ALFA(J+1)*DT*RESMT1(I)

	U(I)=FU(I)/RHO(I)
!------------------------------------------------------------------------
	!!!!!!!! V MOMENTUN

	F71=COF1*(VFHAT1+CHAT1*DX(I,1)/DS(I,1))

	F72=COF2*(CHAT1*DY(I,1)/DS(I,1))

	F73=COF3*VFHAT1

	F74=COF4*(VFHAT1-CHAT1*DX(I,1)/DS(I,1))

	! I-0.5

      F81=COF5*(VFHAT2+CHAT2*DX(I,2)/DS(I,2))

	F82=COF6*(CHAT2*DY(I,2)/DS(I,2))

	F83=COF7*VFHAT2

	F84=COF8*(VFHAT2-CHAT2*DX(I,2)/DS(I,2))

      
	F91=COF9*(VFHAT3+CHAT3*DX(I,3)/DS(I,3))

	F92=COF10*(CHAT3*DY(I,3)/DS(I,3))

	F93=COF11*VFHAT3

	F94=COF12*(VFHAT3-CHAT3*DX(I,3)/DS(I,3))

	IF (NEIB(I,4).NE.-2) THEN

	FR1=COF13*(VFHAT4+CHAT4*DX(I,4)/DS(I,4))

	FR2=COF14*(CHAT4*DY(I,4)/DS(I,4))

	FR3=COF15*VFHAT4

	FR4=COF16*(VFHAT4-CHAT4*DX(I,4)/DS(I,4))

	END IF

	DR1=1./2.*(RHOO(I)*QLO(I,1)*VO(I)-PO(I)*DX(I,1)/DS(I,1)+(K1)*
	
     &(RHOO(NE1)*QLO(NE1,L1)*VO(NE1)-PO(NE1)*DX(NE1,L1)/DS(NE1,L1)))
	
     &-1./2.*(F71+F72+F73+F74)

	DR2=1./2.*(RHOO(I)*QLO(I,2)*VO(I)-PO(I)*DX(I,2)/DS(I,2)+(K2)*
	
     &(RHOO(NE2)*QLO(NE2,L2)*VO(NE2)-PO(NE2)*DX(NE2,L2)/DS(NE2,L2)))
     
     &-1./2.*(F81+F82+F83+F84)

	DR3=1./2.*(RHOO(I)*QLO(I,3)*VO(I)-PO(I)*DX(I,3)/DS(I,3)+(K3)*
	
     &(RHOO(NE3)*QLO(NE3,L3)*VO(NE3)-PO(NE3)*DX(NE3,L3)/DS(NE3,L3)))
     
     &-1./2.*(F91+F92+F93+F94)

	IF (NEIB(I,4).NE.-2) THEN
     		
      DR4=1./2.*(RHOO(I)*QLO(I,4)*VO(I)-PO(I)*DX(I,4)/DS(I,4)+(K4)*
	
     &(RHOO(NE4)*QLO(NE4,L4)*VO(NE4)-PO(NE4)*DX(NE4,L4)/DS(NE4,L4)))
     
     &-1./2.*(FR1+FR2+FR3+FR4)

	END IF 
     		
	FINVR(I)=DR1*DS(I,1)+DR2*DS(I,2)+DR3*DS(I,3)+IW*DR4*DS(I,4) !INVISCID

	FVV=FVV1+FVV2+FVV3+IW*FVV4 !VISCOUS
    
	! CALCULATE FLUX RESIDUAL

	RESMR1(I)=-(FINVR(I)+FVV)/AREA(I)

	!! UPDATE R-MOM
      
	FV(I)=FVO(I)+DT*RESMR1(I)
        
c	FVO(I)=FVOZ(I)+ALFA(J+1)*DT*RESMR1(I)

	V(I)=FV(I)/RHO(I)
!------------------------------------------------------------------------
	!!!! ENERGY EQUATION
			
	! I+0.5 

      F101=COF1*(HFHAT1-QFHAT1*CHAT1)

	F102=COF2*RFHAT1*CHAT1

	F103=COF3*1./2.*(UFHAT1**2+VFHAT1**2) 

	F104=COF4*(HFHAT1+QFHAT1*CHAT1)

	! I-0.5

      F111=COF5*(HFHAT2-QFHAT2*CHAT2)

	F112=COF6*RFHAT2*CHAT2

	F113=COF7*1./2.*(UFHAT2**2+VFHAT2**2)

	F114=COF8*(HFHAT2+QFHAT2*CHAT2)


	F121=COF9*(HFHAT3-QFHAT3*CHAT3)

	F122=COF10*RFHAT3*CHAT3

	F123=COF11*1./2.*(UFHAT3**2+VFHAT3**2)

	F124=COF12*(HFHAT3+QFHAT3*CHAT3)
	
	IF (NEIB(I,4).NE.-2) THEN
	
	FE1=COF13*(HFHAT4-QFHAT4*CHAT4)

	FE2=COF14*RFHAT4*CHAT4

	FE3=COF15*1./2.*(UFHAT4**2+VFHAT4**2)

	FE4=COF16*(HFHAT4+QFHAT4*CHAT4)

	END IF			

	DEP1=1./2.*(RHOO(I)*QLO(I,1)*HO(I)+(K1)*(RHOO(NE1)*QLO(NE1,L1)
	
     &*HO(NE1)))-1./2.*(F101+F102+F103+F104)

	DEP2=1./2.*(RHOO(I)*QLO(I,2)*HO(I)+(K2)*(RHOO(NE2)*QLO(NE2,L2)
	
     &*HO(NE2)))-1./2.*(F111+F112+F113+F114)

	DEP3=1./2.*(RHOO(I)*QLO(I,3)*HO(I)+(K3)*(RHOO(NE3)*QLO(NE3,L3)

     &*HO(NE3)))-1./2.*(F121+F122+F123+F124)

	IF (NEIB(I,4).NE.-2) THEN
       
	DEP4=1./2.*(RHOO(I)*QLO(I,4)*HO(I)+(K4)*(RHOO(NE4)*QLO(NE4,L4)

     &*HO(NE4)))-1./2.*(FE1+FE2+FE3+FE4)
	 
      END IF
      		
	FINVE(I)=DEP1*DS(I,1)+DEP2*DS(I,2)+DEP3*DS(I,3)+IW*DEP4*DS(I,4) !INVISCID
	
	FVE=FVE1+FVE2+FVE3+IW*FVE4  !VISCOUS  
      
	!! UPDATE ENERGY 

	RESE1(I)=-(FINVE(I)+FVE)/AREA(I)     

	! UPDATE ENERGY

	FE(I)=FEO(I)+DT*RESE1(I)

	E(I)=FE(I)/RHO(I)  

C	IF (I.EQ.52) WRITE(*,*)	FVT,FVV,FVE
	
	END DO ! ELEMENTS

	CONTINUE

c	END DO ! XMULTI STAGE
     	 
      DO I=1,NE

c      WRITE (15,*) I,RHO(I),U(I)

c      RHO(I)=RHOO(I)
	
c	U(I)=FUO(I)/RHOO(I)
      
c	V(I)=FVO(I)/RHOO(I)

c	E(I)=FEO(I)/RHOO(I)

	END DO

	END SUBROUTINE
!****************************************************************************
!     CALCULATE PARAMETERS
!
!****************************************************************************      
	SUBROUTINE CELLNODE(NUM_NODE,UC,VC,TC)

	INTEGER::NUM_NODE
	DOUBLEPRECISION::UC,VC,TC
	DOUBLEPRECISION::U0,V0

	UC=0.0
	VC=0.0
	TC=0.0
	K=0

	DO M=1,MSUM

	NV=NODE_VICINITY(NUM_NODE,M)

	IF(NV.NE.0)THEN

	RHOD=RHOO(NV)

	U0=UO(NV)

	V0=VO(NV)

	UC=UC+U0

	VC=VC+V0

	PD=(GAMA-1)*(FEO(NV)-0.5*RHOD*(U0**2+V0**2))

	TC=TC+GAMA*XMIN**2*PD/RHOD

	K=K+1

	ELSE
	GOTO 26
	ENDIF

	ENDDO

26    CONTINUE
	UC=UC/K
	VC=VC/K
	TC=TC/K
	
	END SUBROUTINE 
!****************************************************************************
!     CALCULATE PARAMETERS
!
!****************************************************************************      
	Subroutine VISCOUS(ELEM,FACE,TAU_XX,TAU_XY,TAU_YY,Q_X,Q_Y)
	
	INTEGER:: ELEM,FACE	
      
	! CALCULATE THE SUM OF THE RESIDUALS

	! FACE 1

	I=ELEM
	N1=NODENUM(I,1)
	N2=NODENUM(I,2)
	N3=NODENUM(I,3)
	N4=NODENUM(I,4)

	NE1=NEIB(I,1)
	NE2=NEIB(I,2)
	NE3=NEIB(I,3)
	NE4=NEIB(I,4)

	IF (FACE.EQ.1) THEN

	CALL CELLNODE(N1,U1,V1,T1)	
	CALL CELLNODE(N2,U3,V3,T3)

	IF(NEIB(I,1)>0) THEN
	DX1=X_CELL(NE1)-XNODE(N1)
	DY1=Y_CELL(NE1)-YNODE(N1)

	DX2=XNODE(N2)-X_CELL(NE1)
	DY2=YNODE(N2)-Y_CELL(NE1)

	DX3=X_CELL(I)-XNODE(N2)
	DY3=Y_CELL(I)-YNODE(N2)

	DX4=XNODE(N1)-X_CELL(I)
	DY4=YNODE(N1)-Y_CELL(I)
	ENDIF

	! FACE 2

	ELSEIF (FACE.EQ.2) THEN

	CALL CELLNODE(N2,U1,V1,T1)
	CALL CELLNODE(N3,U3,V3,T3)
	IF(NE2>0)THEN
	DX1=X_CELL(NE2)-XNODE(N2)
	DY1=Y_CELL(NE2)-YNODE(N2)

	DX2=XNODE(N3)-X_CELL(NE2)
	DY2=YNODE(N3)-Y_CELL(NE2)

	DX3=X_CELL(I)-XNODE(N3)
	DY3=Y_CELL(I)-YNODE(N3)

	DX4=XNODE(N2)-X_CELL(I)
	DY4=YNODE(N2)-Y_CELL(I)
	ENDIF

	!FACE 3
	
	ELSEIF (FACE.EQ.3) THEN

	IF(N4.EQ.0)THEN
	CALL CELLNODE(N3,U1,V1,T1)
	CALL CELLNODE(N1,U3,V3,T3)
	IF(NE3>0)THEN
	DX1=X_CELL(NE3)-XNODE(N3)
	DY1=Y_CELL(NE3)-YNODE(N3)

	DX2=XNODE(N1)-X_CELL(NE3)
	DY2=YNODE(N1)-Y_CELL(NE3)

	DX3=X_CELL(I)-XNODE(N1)
	DY3=Y_CELL(I)-YNODE(N1)

	DX4=XNODE(N3)-X_CELL(I)
	DY4=YNODE(N3)-Y_CELL(I)
	ENDIF
	ELSEIF(N4.NE.0)THEN
	CALL CELLNODE(N3,U1,V1,T1)
	CALL CELLNODE(N4,U3,V3,T3)
	IF(NE3>0)THEN
	DX1=X_CELL(NE3)-XNODE(N3)
	DY1=Y_CELL(NE3)-YNODE(N3)

	DX2=XNODE(N4)-X_CELL(NE3)
	DY2=YNODE(N4)-Y_CELL(NE3)

	DX3=X_CELL(I)-XNODE(N4)
	DY3=Y_CELL(I)-YNODE(N4)

	DX4=XNODE(N3)-X_CELL(I)
	DY4=YNODE(N3)-Y_CELL(I)
	ENDIF
	ENDIF

	!FACE 4
	
	ELSEIF (FACE.EQ.4) THEN

	IF(N4.NE.0)THEN
	CALL CELLNODE(N4,U1,V1,T1)
	CALL CELLNODE(N1,U3,V3,T3)
	IF(NE4>0)THEN
	DX1=X_CELL(NE4)-XNODE(N4)
	DY1=Y_CELL(NE4)-YNODE(N4)

	DX2=XNODE(N1)-X_CELL(NE4)
	DY2=YNODE(N1)-Y_CELL(NE4)

	DX3=X_CELL(I)-XNODE(N1)
	DY3=Y_CELL(I)-YNODE(N1)

	DX4=XNODE(N4)-X_CELL(I)
	DY4=YNODE(N4)-Y_CELL(I)
	ENDIF
	ENDIF
	ENDIF

	IF(NEIB(I,FACE)>0)THEN

	IF ((N4.EQ.0).AND.(NODENUM(NEIB(I,FACE),4).EQ.0))THEN
	AREACV=(AREA(I)+AREA(NEIB(I,FACE)))/3.0

	ELSE IF((N4.EQ.0).AND.(NODENUM(NEIB(I,FACE),4)/=0))THEN
	AREACV=(AREA(I))/3.0+(AREA(NEIB(I,FACE)))/4.0

	ELSE IF((N4/=0).AND.(NODENUM(NEIB(I,FACE),4).EQ.0))THEN
	AREACV=(AREA(I))/4.0+(AREA(NEIB(I,FACE)))/3.0

	ELSE IF((N4/=0).AND.(NODENUM(NEIB(I,FACE),4)/=0))THEN
	AREACV=(AREA(I)+AREA(NEIB(I,FACE)))/4.0
	ENDIF
	
	J=NEIB(I,FACE)

	
	DU_DX=1.0/AREACV*
     &((U1+UO(J))/2.0*DY1+(UO(J)+U3)/2.0*DY2+
     &(U3+UO(I))/2.0*DY3+(UO(I)+U1)/2.0*DY4)

	DU_DY=-1.0/AREACV*
     &((U1+UO(J))/2.0*DX1+(UO(J)+U3)/2.0*DX2+
     &(U3+UO(I))/2.0*DX3+(UO(I)+U1)/2.0*DX4)

	DV_DX= 1.0/AREACV*
     &((V1+VO(J))/2.0*DY1+(VO(J)+V3)/2.0*DY2+
     &(V3+VO(I))/2.0*DY3+(VO(I)+V1)/2.0*DY4)

	DV_DY=-1.0/AREACV*
     &((V1+VO(J))/2.0*DX1+(VO(J)+V3)/2.0*DX2+
     &(V3+VO(I))/2.0*DX3+(VO(I)+V1)/2.0*DX4)

	DT_DX=1.0/AREACV*
     &((T1+TOLD(J))/2.0*DY1+(TOLD(J)+T3)/2.0*DY2+
     &(T3+TOLD(I))/2.0*DY3+(TOLD(I)+T1)/2.0*DY4)

	DT_DY=-1.0/AREACV*
     &((T1+TOLD(J))/2.0*DX1+(TOLD(J)+T3)/2.0*DX2+
     & (T3+TOLD(I))/2.0*DX3+(TOLD(I)+T1)/2.0*DX4)

	ELSE IF (NW(I).EQ.1) THEN !NEIB(I,FACE).EQ.-1

	AREACV=(AREA(I))/4.0
	DU_DX=-0.5/AREACV*UO(I)*DY(I,FACE)
	DU_DY=0.5/AREACV*UO(I)*DX(I,FACE)
	DV_DX=-0.5/AREACV*VO(I)*DY(I,FACE)
	DV_DY=0.5/AREACV*VO(I)*DX(I,FACE)

	DT_DX=0.0
	DT_DY=0.0
	ENDIF

	PR=0.7
	T_H=0.5*(TOLD(J)+TOLD(I))
	XMU=(1+110.4/300)*(T_H**(1.5))/(T_H+110.4/300)
	

	TAU_XX=(2.0/3.0)*(XMU/RE)*(2*DU_DX-DV_DY)
	TAU_YY=(2.0/3.0)*(XMU/RE)*(2*DV_DY-DU_DX)
	TAU_XY=(XMU/RE)*(DU_DY+DV_DX)
	QX=-(XMU/((GAMA-1)*RE*PR*XMIN**2))*DT_DX
	QY=-(XMU/((GAMA-1)*RE*PR*XMIN**2))*DT_DY

c	IF (I.EQ.51) WRITE(*,*) J !U1,U3,T1,T3

	END SUBROUTINE	
!****************************************************************************
!     CALCULATE PARAMETERS
!
!****************************************************************************      
	subroutine UPDATE	
      
	! CALCULATE THE SUM OF THE RESIDUALS

	DO I=1,NE
	      
	SUM1(I)=1. 
	SUM2(I)=1. 
	SUM3(I)=1. 
	SUM4(I)=1. 
	  
	END DO

	! PUT THE VALUES IN THE LAST TIME STEP SOLUTION

	DO I=1,NE
      
	VO(I)=V(I) 	      

	UO(I)=U(I)       

      RHOO(I)=RHO(I)

	FUO(I)=FU(I)

	FVO(I)=FV(I)	

	FEO(I)=FE(I)
	
	EO(I)=E(I)

C	TOLD(I)=T(I)
	
	END DO

	END SUBROUTINE	
!****************************************************************************
!     CALCULATE PARAMETERS
!
!****************************************************************************
      subroutine CALC
      	
	DO I=1,NE
	
	ES=E(I)-0.5D0*(U(I)**2+V(I)**2)	

	GM=GAMA*(GAMA-1)

      C(I)=SQRT(ES*GM)

	CO(I)=C(I)

      P(I)=RHO(I)*C(I)**2/GAMA

	PO(I)=P(I)

	XM(I)=SQRT(U(I)**2+V(I)**2)/C(I)
	
	DC=(GAMA-1)/(2*GAMA)*(U(I)**2+V(I)**2)	
	
	H(I)=GAMA*(E(I)-DC)

	HO(I)=H(I)      
      
      T(I)=GAMA*XMIN**2*P(I)/RHO(I) !ES/CV*(UI**2/TI)

	TOLD(I)=T(I)

	END DO		         

	END SUBROUTINE     
!****************************************************************************
!     CHECK CONVERGENCE
!
!****************************************************************************
      subroutine CONVERGENCE 

	DO I=1,NE
	
	RESC(I)=ABS(RESC1(I)/SUM1(I))

	RESMTT(I)=ABS(RESMT1(I)/SUM2(I))

      RESMRT(I)=ABS(RESMR1(I)/SUM3(I))

      RESET(I)=ABS(RESE1(I)/SUM4(I))

	END DO
	
	RESCMAX=RESC(1)
	
	RESEMAX=RESET(1)

	RESMRMAX=RESMRT(1)

	RESMTMAX=RESMTT(1)
	
	DO I=1,NE
       
      IF (RESC(I).GT.RESCMAX)   RESCMAX=RESC(I)

	IF (RESET(I).GT.RESEMAX)   RESEMAX=RESET(I)

	IF (RESMRT(I).GT.RESMRMAX) RESMRMAX=RESMRT(I)

	IF (RESMTT(I).GT.RESMTMAX) RESMTMAX=RESMTT(I)

	END DO

	E11=0.
	E2=0.
	E3=0.
	E4=0.

	DO I=1,NE
      
	E11=E11+ABS(RESC(I))**2
	E2=E2+ABS(RESMTT(I))**2
	E3=E3+ABS(RESMRT(I))**2
	E4=E4+ABS(RESET(I))**2

	END DO

      E11=SQRT(E11)/(NE)
	E2=SQRT(E2)/(NE)
	E3=SQRT(E3)/(NE)
	E4=SQRT(E4)/(NE)

c	RESIDUAL=MAX(RESCMAX,RESEMAX,RESMRMAX,RESMTMAX)	

      RESIDUAL=MAX(E11,E2,E3,E4)	
	
	WRITE (*,*)  RESIDUAL 

	WRITE (102,*)  MC,RESIDUAL
	
	MC=MC+1 
	
C	TIME=TIME+DT

	END SUBROUTINE
!****************************************************************************
	Subroutine Output1
		integer :: i,j
		
		open(45,file='GRID HYBRID.plt')

          TY=0

	   write(45,*),'VARIABLES="X" "Y" "TY"'
	   write(45,*),'ZONE F=FEPOINT,ET=QUADRILATERAL,N=', NN, ',E=', NE
   
		do I=1,NN
		write(45,*) Xnode(I),Ynode(I),TY(I)
		end do


	DO I=1,NE

      IF(NodeNum(I,4).EQ.0)THEN

      WRITE(45,'(3I10)')NodeNum(I,1),NodeNum(I,2),NodeNum(I,3)
     &	,NodeNum(I,3)
      ELSE

      WRITE(45,'(4I10)')NodeNum(I,1),NodeNum(I,2),NodeNum(I,3)
     &	,NodeNum(I,4)
      
	END IF
      END DO 

	end Subroutine 
!****************************************************************************
	Subroutine Output
		integer :: i,j
		
		open(5,file='output.plt')
	    OPEN(12,file='output.DAT')

	write(5,*),'VARIABLES="X" "Y" "U" "V" "P" "MACH" "CPP" "S" "RHO"'
	write(5,*),'ZONE F=FEPOINT,ET=QUADRILATERAL,N=', NN, ',E=', NE

	    
		UU=0
		VV=0

		PP=0
		CC=0

		
		DO I=1,NE

		K1=NodeNum(I,1)
		K2=NodeNum(I,2)
		K3=NodeNum(I,3)
	    K4=NodeNum(I,4)

		NE1=NEIB(I,1)
		NE2=NEIB(I,2)
		NE3=NEIB(I,3)	
	    NE4=NEIB(I,4)

		UU(K1)=U(I)*AREA(I)	
		VV(K1)=V(I)*AREA(I)
		PP(K1)=P(I)*AREA(I)
		CC(K1)=C(I)*AREA(I)
		RR(K1)=RHO(I)*AREA(I)

		UU(K2)=U(I)*AREA(I)	
		VV(K2)=V(I)*AREA(I)
		PP(K2)=P(I)*AREA(I)
		CC(K2)=C(I)*AREA(I)
		RR(K2)=RHO(I)*AREA(I)


		UU(K3)=U(I)*AREA(I)	
		VV(K3)=V(I)*AREA(I)
		PP(K3)=P(I)*AREA(I)
		CC(K3)=C(I)*AREA(I)
		RR(K3)=RHO(I)*AREA(I) 

	    IF (K4.NE.0) THEN
	    UU(K4)=U(I)*AREA(I)	
		VV(K4)=V(I)*AREA(I)
		PP(K4)=P(I)*AREA(I)
		CC(K4)=C(I)*AREA(I)
		RR(K4)=RHO(I)*AREA(I) 

	    END IF

		DOM1=AREA(I)
		DOM2=AREA(I)
		DOM3=AREA(I)
	    DOM4=AREA(I)


		DO J=1,NE

		IF (I.NE.J) THEN

		IF (NodeNum(J,1).EQ.K1.OR.NodeNum(J,2).EQ.K1     
     &	.OR.NodeNum(J,3).EQ.K1.OR.NodeNum(J,4).EQ.K1) THEN

		UU(K1)=UU(K1)+U(J)*AREA(J)
		VV(K1)=VV(K1)+V(J)*AREA(J)
		PP(K1)=PP(K1)+P(J)*AREA(J)
		CC(K1)=CC(K1)+C(J)*AREA(J)
		RR(K1)=RR(K1)+RHO(J)*AREA(J)
		
		DOM1=DOM1+AREA(J)

		END IF

		IF (NodeNum(J,1).EQ.K2.OR.NodeNum(J,2).EQ.K2     
     &	.OR.NodeNum(J,3).EQ.K2.OR.NodeNum(J,4).EQ.K2) THEN

		UU(K2)=UU(K2)+U(J)*AREA(J)
		VV(K2)=VV(K2)+V(J)*AREA(J)
		PP(K2)=PP(K2)+P(J)*AREA(J)
		CC(K2)=CC(K2)+C(J)*AREA(J)
		RR(K2)=RR(K2)+RHO(J)*AREA(J)
		
		DOM2=DOM2+AREA(J)

		END IF

		IF (NodeNum(J,1).EQ.K3.OR.NodeNum(J,2).EQ.K3     
     &	.OR.NodeNum(J,3).EQ.K3.OR.NodeNum(J,4).EQ.K3) THEN

		UU(K3)=UU(K3)+U(J)*AREA(J)
		VV(K3)=VV(K3)+V(J)*AREA(J)
		PP(K3)=PP(K3)+P(J)*AREA(J)
		CC(K3)=CC(K3)+C(J)*AREA(J)
		RR(K3)=RR(K3)+RHO(J)*AREA(J)

		DOM3=DOM3+AREA(J)

		END IF

	    IF (K4.NE.0) THEN

		IF (NodeNum(J,1).EQ.K4.OR.NodeNum(J,2).EQ.K4     
     &	.OR.NodeNum(J,3).EQ.K4.OR.NodeNum(J,4).EQ.K4) THEN

		UU(K4)=UU(K4)+U(J)*AREA(J)
		VV(K4)=VV(K4)+V(J)*AREA(J)
		PP(K4)=PP(K4)+P(J)*AREA(J)
		CC(K4)=CC(K4)+C(J)*AREA(J)
		RR(K4)=RR(K4)+RHO(J)*AREA(J)

		DOM4=DOM4+AREA(J)

		END IF

		END IF

		END IF	   
	  			  
		END DO
		
		UU(K1)=UU(K1)/DOM1
		VV(K1)=VV(K1)/DOM1
		PP(K1)=PP(K1)/DOM1
		CC(K1)=CC(K1)/DOM1
		RR(K1)=RR(K1)/DOM1

		UU(K2)=UU(K2)/DOM2
		VV(K2)=VV(K2)/DOM2
		PP(K2)=PP(K2)/DOM2
		CC(K2)=CC(K2)/DOM2
		RR(K2)=RR(K2)/DOM2

		UU(K3)=UU(K3)/DOM3
		VV(K3)=VV(K3)/DOM3
		PP(K3)=PP(K3)/DOM3
		CC(K3)=CC(K3)/DOM3
		RR(K3)=RR(K3)/DOM3

		IF (K4.NE.0) THEN
		UU(K4)=UU(K4)/DOM4
		VV(K4)=VV(K4)/DOM4
		PP(K4)=PP(K4)/DOM4
		CC(K4)=CC(K4)/DOM4
		RR(K4)=RR(K4)/DOM4
		END IF
		
		END DO

		
		DO I=1,NN

          CPP(I)=(PP(I)-CINP**2/GAMA)/(0.5D0*RHOI*1.**2) 
		SO(I)=(PP(I)/RR(I)**GAMA)-1
		  
		write(5,*) Xnode(I),Ynode(I),UU(I),VV(I),PP(I),
     &	SQRT(UU(I)**2+VV(I)**2)/CC(I),CPP(I),SO(I),RR(I)
	
		end do

		DO I=1,NE
	    IF(NodeNum(I,4).EQ.0)THEN

           WRITE(5,*)NodeNum(I,1),NodeNum(I,2),NodeNum(I,3)
     &	,NodeNum(I,3)
           
		 ELSE

           WRITE(5,*)NodeNum(I,1),NodeNum(I,2),NodeNum(I,3)
     &	,NodeNum(I,4)
      
      	END IF
		end do

	    Do I=1,NE
		write(12,*) SO(I),CPP(I) !,P(I)
		end do

	 !CALCULATE FORCES

	   FX=0
	   FY=0
	   
	   DO I=1,NE
         
	   IF (NW(I).EQ.1) THEN	   

	   DO K=1,3
	  
	   IF (NEIB(I,K).EQ.ZERO) THEN

	   FX=FX+NM1(I)*P(I)*DY(I,K)

	   FY=FY+NM1(I)*P(I)*DX(I,K)
	   
c	   WRITE(15,*) I,P(I)*DY(I,K),P(I)*DX(I,K)       
   
	   END IF
	   END DO	   
	   END IF
	   END DO
	   
c	   WRITE(15,*) -FX*2,FY*2	

	end Subroutine Output
!****************************************************************************
	 !POST_PROCESSING
!****************************************************************************
	SUBROUTINE POST_PROCESSING

	DOUBLEPRECISION ::DENSITY(NN),U(NN),V(NN),PRESSURE(NN),MACH(NN)
     &	,ENTROPY(NN),CP(NN),PT(NN),DENSITY1(NN)

	DOUBLE PRECISION::QUANTITY1(NN),QUANTITY2(NN),QUANTITY3(NN)
     &	,QUANTITY4(NN),SUM_AREA

	INTEGER::CELL,I,N1,N2,WALL(100,2)

	QUANTITY=0
	DO I=1,NN
	SUM_AREA=0

	DO CELL=1,NE
	IF((NODENUM(CELL,1)==I).OR.(NODENUM(CELL,2)==I).OR.
     &(NODENUM(CELL,3)==I).OR.(NODENUM(CELL,4)==I))THEN

	QUANTITY1(I)=QUANTITY1(I)+RHO(CELL)*AREA(CELL)
	QUANTITY2(I)=QUANTITY2(I)+UO(CELL)*AREA(CELL)
	QUANTITY3(I)=QUANTITY3(I)+VO(CELL)*AREA(CELL)
	QUANTITY4(I)=QUANTITY4(I)+E(CELL)*AREA(CELL)

	SUM_AREA=SUM_AREA+AREA(CELL)

	ENDIF
	ENDDO

	QUANTITY1(I)=QUANTITY1(I)/SUM_AREA
	QUANTITY2(I)=QUANTITY2(I)/SUM_AREA
	QUANTITY3(I)=QUANTITY3(I)/SUM_AREA
	QUANTITY4(I)=QUANTITY4(I)/SUM_AREA

	ENDDO


	DO I=1,NN

	DENSITY1(I)=QUANTITY1(I)
	U(I)=QUANTITY2(I)
	V(I)=QUANTITY3(I)
	
	PRESSURE(I)=(GAMA-1)*(QUANTITY4(I)*QUANTITY1(I)
     &		   -0.5*(QUANTITY2(I)**2+QUANTITY3(I)**2))
	
	MACH(I)=DSQRT(U(I)**2+V(I)**2)/DSQRT(GAMA*PRESSURE(I)/DENSITY1(I))
	

	CP(I)=2*(PRESSURE(I)-PIN)
	ENTROPY(I)=(DENSITY1(I)**(-GAMA))*PRESSURE(I)/PIN-1
	PT(I)=PRESSURE(I)+0.5*DENSITY1(I)*(U(I)**2+V(I)**2)
	ENDDO	


	OPEN(4,FILE='AIRFOILX.PLT')
	OPEN(5,FILE='AIRFOILY.PLT')
	WRITE(4,'(A)')'VARIABLES="X","DENSITY","U","V","PRESSURE","MACH"
     &	,"ENTROPY","-CP","PT"'
	WRITE(5,'(A)')'VARIABLES="Y","DENSITY","U","V","PRESSURE","MACH"
     &	,"ENTROPY","-CP","PT"'

	DO I=1,50
	WALL(I,1)=50+I
	END DO
	
	DO I=51,100
	WALL(I,1)=300-50+I
	END DO

	WALL(1:100,2)=2
	
	DO I=1,100
	IF (WALL(I,2)==1) THEN
	N1=NODENUM(WALL(I,1),1);N2=NODENUM(WALL(I,1),2)
	ELSEIF (WALL(I,2)==2) THEN
	N1=NODENUM(WALL(I,1),2);N2=NODENUM(WALL(I,1),3)
	ELSEIF (WALL(I,2)==3) THEN
	N1=NODENUM(WALL(I,1),3);N2=NODENUM(WALL(I,1),1)
	ENDIF

	WRITE(4,'(9ES20.4)') XNODE(N1),DENSITY1(N1),U(N1),V(N1)
     &,PRESSURE(N1),MACH(N1),ENTROPY(N1),-CP(N1),PT(N1)
	WRITE(4,'(9ES20.4)') XNODE(N2),DENSITY1(N2),U(N2),V(N2)
     &,PRESSURE(N2),MACH(N2),ENTROPY(N2),-CP(N2),PT(N2)
	WRITE(5,'(9ES20.4)') YNODE(N1),DENSITY1(N1),U(N1),V(N1)
     &,PRESSURE(N1),MACH(N1),ENTROPY(N1),-CP(N1),PT(N1)
	WRITE(5,'(9ES20.4)') YNODE(N2),DENSITY1(N2),U(N2),V(N2)
     &,PRESSURE(N2),MACH(N2),ENTROPY(N2),-CP(N2),PT(N2)

	ENDDO

	CLOSE(4)
	CLOSE(5)

	CL=0;CD=0
	FX=0;FY=0;
	DO I=1,100
	NUM=WALL(I,1);FACE=WALL(I,2)
	PF=(GAMA-1)*(EO(NUM)*RHOO(NUM)-0.5*(UO(NUM)**2+VO(NUM)**2))
	FX=FX+PF*DY(NUM,FACE)
	FY=FY-PF*DX(NUM,FACE)
	ENDDO
	CL=2*(FY*DCOS(ALF)-FX*DSIN(ALF))
	CD=2*(FY*DSIN(ALF)+FX*DCOS(ALF))


	WRITE(15,*) CL,CD


	END SUBROUTINE POST_PROCESSING
	End program 

