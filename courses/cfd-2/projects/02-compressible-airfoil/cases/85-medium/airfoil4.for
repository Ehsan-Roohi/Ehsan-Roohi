!****************************************************************************
!
!  PROGRAM: AIRFOIL
!
!
!****************************************************************************

	program AIRFOIL
	
	    USE MSIMSL
		USE MSFLIB
	
	PARAMETER (NN=1298)
	PARAMETER (NE=2526) 
	PARAMETER (NEBC=3000)
	PARAMETER (MM=6)   !R-K INTEGRATION STEPS

	implicit DoublePrecision  (A-H,O-Z)

	! Variables

	Integer     :: NNode,NBound,Neib(NE,3),FACE(NE,3),NW(NEBC),NM1(NE)

	DoublePrecision ::T(NEBC),AK(3,3),AR(NEBC),A(NEBC,NEBC),R,EPS,CINP,
     &RHOOZ(NEBC),FUOZ(NEBC),FVOZ(NEBC),FEOZ(NEBC),CPP(NEBC),SO(NEBC)
     &,ANG(NEBC)

	DoublePrecision :: Xnode(NEBC),Ynode(NEBC),H(NEBC),C(NEBC),U(NEBC)
     &,V(NEBC),RHO(NEBC),RHOO(NEBC),P(NEBC),PO(NEBC),ULO(NEBC,3)
     &,VLO(NEBC,3),QLO(NEBC,3),RLO(NEBC,3),E11,E2,E3,E4,TY(NN),CC(NEBC)

	DoublePrecision::UO(NEBC),VO(NEBC),FU(NEBC),FV(NEBC),FUO(NEBC)
     &,FVO(NEBC) 
	
	DoublePrecision :: E(NEBC),EF(NEBC),EO(NEBC),HO(NEBC),FE(NEBC)	
     &,FEO(NEBC),TOLD(NEBC),TT(NE),PP(NE),RR(NE),UU(NE),VV(NE)
	     
	
	DoublePrecision ::RESC1(NEBC),RESE1(NEBC),RESMR1(NEBC),RESMT1(NEBC)
	DoublePrecision ::RESC2(NEBC),RESE2(NEBC),RESMR2(NEBC),RESMT2(NEBC)
      DoublePrecision ::RESC(NEBC),RESE(NEBC),RESMR(NEBC),RESMT(NEBC)
		
	DoublePrecision :: DUTDT1(NEBC),DUTDT2(NEBC),DURDT1(NEBC),
     &UR1(NEBC),UR2(NEBC),UT1(NEBC),UT2(NEBC),FMTV(NEBC),FMRV(NEBC)
     &,TIME,DF(NEBC),DFP(NEBC),DFM(NEBC),DRP(NEBC),DRM(NEBC),DEP(NEBC)
     &,DEM(NEBC),DURDT2(NEBC),FEV(NEBC)	

      DoublePrecision :: FINVT(NEBC),FINVR(NEBC),FINVE(NEBC),SCI(NEBC)
     &,SMTI(NEBC),SEI(NEBC),SEV(NEBC),TR(NEBC),TL(NEBC),XK(NEBC)

      DoublePrecision :: XM(NEBC),URHAT(NEBC),CO(NEBC),DT

      DoublePrecision :: UI,TI,URHAT1,UTHAT1,URHAT2,EADD	
     &,UTHAT2,GAMA,UHAT1,UHAT2,RE,CP,PR,RHOI,PIN,PRA,CV

	
	DoublePrecision :: TFHAT1,TFHAT2,VFHAT1,VFHAT2,UFHAT1,UFHAT2,	 
     &RESET(NEBC),RESMRT(NEBC),RESMTT(NEBC),SUM1(NEBC),SUM2(NEBC)
     &,DTDT1(NEBC),DTDT2(NEBC),XMIN,SUM4(NEBC),HLO(NEBC,3),SUM3(NEBC)


      DoublePrecision :: DNADX(3),DNADY(3),RHOLO(NEBC,3),PLO(NEBC,3)

	DoublePrecision ::AREV(NEBC,NEBC),AREA(NEBC),DS(NEBC,3),ALFA(10)

      DoublePrecision ::H1(NEBC),C1(NEBC),DX(NEBC,3),DY(NEBC,3),VI

	Integer         :: NodeNum(NEBC,3),BOUNDNum(NEBC),Bound(NEBC,2)

	Integer :: JJ,JJJ,I,J,N1,N2,N3,K,N11,N22,N33,METHOD,BCTYPE,ZERO


	CALL SYSTEM ("del *.obj ")
      CALL SYSTEM ("del *.plg ")
	CALL SYSTEM ("del *.opt ")
	CALL SYSTEM ("del *.dsp ")
	CALL SYSTEM ("del *.dsw ")
	CALL SYSTEM ("del *.txt ")

	! Body of AIRFOIL
	
		CALL INPUT
		CALL MESHInput
		CALL SetNeighbours
		CALL BCNODES	   
		CALL Initialize	
		CALL Output1    
10        CALL OLD		
          CALL FLUX
		CALL UPDATE		
		CALL CALC
	    CALL OLD
	    CALL CONVERGENCE 	
		If (RESIDUAL.GT.EPS) GOTO 10		
20	    CALL Output

	contains
!****************************************************************************
	subroutine Input
      
	    GAMA=1.4D0

	    R=287.D0

	    ZERO=0

	    ALFA(2)=.0485D0 !0.0033

	    ALFA(3)=0.1120D0 !0.2069

	    ALFA(4)=0.1994D0 !0.4265

	    ALFA(5)=0.3285D0 !1.

	    ALFA(6)=0.5477D0

          ALFA(7)=1.  

	    CP=1004.5D0

	    CV=717.5D0

	    UI=1.*COS(1*3.14159/180)

	    VI=1.*SIN(1*3.14159/180)

	    VEL=SQRT(UI**2+VI**2)
		
		RHOI=1.D0	

	    TIME=0

		XMIN=0.85D0 !1.2D0

	    CINP=VEL/XMIN 
		
		PIN=RHOI*CINP**2/GAMA

	    HIN=0.5D0+1./((GAMA-1)*XMIN**2)

	    EIN=HIN-CINP**2/GAMA
		
		TI=PIN/(R*RHOI)	      	    		       	      

	    DT=1.D-2

	    EPS=1.D-7		

	END SUBROUTINE
!****************************************************************************
	Subroutine MESHINPUT
	
		Integer :: I,J,II,III		

		Open(1,file='Node.DAT')
	    
		DO I=1,NN         
	    READ(1,*)  J,XNode(I),YNode(I)
	    END DO

	    Open(2,file='connectivity.DAT')		
		
		DO I=1,NE	
		READ(2,*)  J,II,III,NodeNum(I,1),NodeNum(I,2),NodeNum(I,3)
		END DO

	    CLOSE (2)

	    DO I=1,NE

		N1=NodeNum(I,1)
		N2=NodeNum(I,2)
		N3=NodeNum(I,3)

	    AREA(I)=(XNode(N2)*YNode(N3)-XNode(N3)*YNode(N2)+
     &	XNode(N3)*YNode(N1)-XNode(N1)*YNode(N3)+
     &    XNode(N1)*YNode(N2)-XNode(N2)*YNode(N1))/2.
	    
		DX(I,1)=(XNode(N1)-XNode(N2))
	    DX(I,2)=(XNode(N2)-XNode(N3))
	    DX(I,3)=(XNode(N3)-XNode(N1))

	    DY(I,1)=(YNode(N1)-YNode(N2))
	    DY(I,2)=(YNode(N2)-YNode(N3))
	    DY(I,3)=(YNode(N3)-YNode(N1))
		
		! SIDE LENGTH OF THE ELEMENT
  
	    DS(I,1)=SQRT(DX(I,1)**2+DY(I,1)**2)
	    DS(I,2)=SQRT(DX(I,2)**2+DY(I,2)**2)
	    DS(I,3)=SQRT(DX(I,3)**2+DY(I,3)**2)

	    END DO
	End subroutine
!****************************************************************************
	Subroutine SetNeighbours

		    INTEGER :: I,J,N1,N2,N3,N4

		    Neib=0

		    Do i=1,NE

			N1=NodeNum(i,1)
			N2=NodeNum(i,2)
			N3=NodeNum(i,3)			

			Do J=1,NE

			IF (I.NE.J) THEN

	        IF ((N1==NodeNum(j,1)).and.(N2==NodeNum(j,2))) THEN
			NEIB(I,1)=J
	        FACE(J,1)=I
	        END IF

	        IF ((N2==NodeNum(j,1)).and.(N1==NodeNum(j,2))) THEN
			NEIB(I,1)=J
	        FACE(J,1)=I
	        END IF

	        IF ((N1==NodeNum(j,2)).and.(N2==NodeNum(j,3))) THEN
			NEIB(I,1)=J
	        FACE(J,2)=I
	        END IF

			IF ((N2==NodeNum(j,2)).and.(N1==NodeNum(j,3))) THEN 
			NEIB(I,1)=J 
			FACE(J,2)=I
	        END IF 


	        IF ((N1==NodeNum(j,3)).and.(N2==NodeNum(j,1))) THEN
			NEIB(I,1)=J
	        FACE(J,3)=I
	        END IF                


	        IF ((N2==NodeNum(j,3)).and.(N1==NodeNum(j,1))) THEN 
			NEIB(I,1)=J
	        FACE(J,3)=I
	        END IF 


	       	                   
		    IF ((N2==NodeNum(j,1)).and.(N3==NodeNum(j,2))) THEN 
			NEIB(I,2)=J
              FACE(J,1)=I
			END IF   

	        IF ((N3==NodeNum(j,1)).and.(N2==NodeNum(j,2))) THEN
			NEIB(I,2)=J
	        FACE(J,1)=I
			END IF  

	        IF ((N2==NodeNum(j,2)).and.(N3==NodeNum(j,3))) THEN 
			NEIB(I,2)=J
	        FACE(J,2)=I
			END IF  

	        IF ((N3==NodeNum(j,2)).and.(N2==NodeNum(j,3))) THEN 
			NEIB(I,2)=J
	        FACE(J,2)=I
			END IF 


	        IF ((N2==NodeNum(j,3)).and.(N3==NodeNum(j,1))) THEN 
			NEIB(I,2)=J
	        FACE(J,3)=I
			END IF 

	        IF ((N3==NodeNum(j,3)).and.(N2==NodeNum(j,1))) THEN 
			NEIB(I,2)=J
	        FACE(J,3)=I
			END IF 

	        
			IF ((N3==NodeNum(j,1)).and.(N1==NodeNum(j,2))) THEN
			NEIB(I,3)=J
	        FACE(J,1)=I
			END IF 

	        IF ((N1==NodeNum(j,1)).and.(N3==NodeNum(j,2))) THEN 
			NEIB(I,3)=J
	        FACE(J,1)=I
			END IF 

	        IF ((N3==NodeNum(j,2)).and.(N1==NodeNum(j,3))) THEN 
			NEIB(I,3)=J
	        FACE(J,2)=I
			END IF 

	        IF ((N1==NodeNum(j,2)).and.(N3==NodeNum(j,3))) THEN 
			NEIB(I,3)=J
	        FACE(J,2)=I
			END IF 

	        IF ((N3==NodeNum(j,3)).and.(N1==NodeNum(j,1))) THEN 
			NEIB(I,3)=J
	        FACE(J,3)=I
			END IF 

	        IF ((N1==NodeNum(j,3)).and.(N3==NodeNum(j,1))) THEN 
			NEIB(I,3)=J 
	        FACE(J,3)=I
			END IF 
	       
	        END IF

	        END DO
	        END DO		

	end subroutine 
!****************************************************************************
	Subroutine BCNODES

	    ! RECOGNIZE B.C CELLES

		Open(2,FILE='RES.DAT')

		Open(3,file='BC.DAT')

	    OPEN(15,file='BE.DAT')
	    
	    
		I=1

		DO while(I<500)

  		READ(3,*,END=101) BOUNDNum(I)

		I=I+1

		END DO

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
			
C			WRITE(15,*) K,BOUND(K,1),BOUND(K,2)	        
			              
			END IF   

	     END DO
	     END DO 

       ! RECOGNIZE WALL BC
	 
	   NW=0
	   
	   DO I=1,19,2

         NW(I)=1 
	   END DO

	   DO I=20,38,2

	   NW(I)=1  
	   END DO 	                              	         
		   	   
	End subroutine
!****************************************************************************
	Subroutine Initialize

	    DO I=1,NE

	    UI=1.*COS(1*3.14159/180)

	    VI=1.*SIN(1*3.14159/180)  

	    RHO(I)=1.

	    C(I)=1./XMIN
		  
	    P(I)=RHO(I)*CINP**2/GAMA	   

	    H(I)=0.5D0+1./((GAMA-1)*XMIN**2)

	    E(I)=H(I)-C(I)**2/GAMA	    

	    XM(I)=U(I)/C(I)

	    END DO
		
		UO=U	    
	    VO=V	   
	    PO=P	    
	    HO=H	    
	    RHOO=RHO  
	    CO=C	    
	    EO=E
			
	End subroutine 
!****************************************************************************
	Subroutine OLD 
    	
		 		 
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
	
	DO J=1,6

	DO I=1,NE

	FUO(I)=RHOO(I)*UO(I)

	FVO(I)=RHOO(I)*VO(I)

	FEO(I)=RHOO(I)*EO(I)

	END DO

	IF (J.EQ.1) THEN

      DO I=1,NE

	RHOOZ(I)=RHOO(I)

	FVOZ(I)=FVO(I)

	FUOZ(I)=FUO(I)

	FEOZ(I)=FEO(I)	

	END DO
	
	END IF

	DO I=1,NE

	QLO(I,1)=(UO(I)*DY(I,1)-VO(I)*DX(I,1))/DS(I,1)
      QLO(I,2)=(UO(I)*DY(I,2)-VO(I)*DX(I,2))/DS(I,2)
      QLO(I,3)=(UO(I)*DY(I,3)-VO(I)*DX(I,3))/DS(I,3)

	RLO(I,1)=(UO(I)*DX(I,1)+VO(I)*DY(I,1))/DS(I,1)
      RLO(I,2)=(UO(I)*DX(I,2)+VO(I)*DY(I,2))/DS(I,2)
      RLO(I,3)=(UO(I)*DX(I,3)+VO(I)*DY(I,3))/DS(I,3)

	END DO	

	!!! CONTINUITY  
		
	DO I=1,NE

	NE1=NEIB(I,1)
	NE2=NEIB(I,2)
	NE3=NEIB(I,3)

      N1=NodeNum(i,1)
	N2=NodeNum(i,2)
	N3=NodeNum(i,3)

	IF (NE1.EQ.ZERO)       THEN 
	
	DO K=1,3
	
	IF (FACE(I,K).EQ.ZERO) L=K

	END DO
	
	NE1=NE1+NEBC

	DX(NE1,L)=DX(I,L)

	DY(NE1,L)=DY(I,L)

	DS(NE1,L)=DS(I,L) 

	ANG(I)=ATAN(YNODE(N1)/(XNODE(N1)-0.5D0))
	
	XKAPA1=0.5D0*0.33D0*SQRT(1-XMIN**2)

      XKAPA2=2*3.14159*10*(1-XMIN**2*(SIN(ANG(I)-2*3.14159/180.))**2)

	XKAPA=XKAPA1/XKAPA2

	UO(NE1)=UI+XKAPA*SIN(ANG(I))

	VO(NE1)=VI-XKAPA*COS(ANG(I))     

      QLO(NE1,L)=(UO(NE1)*DY(I,L)-VO(NE1)*DX(I,L))/DS(I,L) 

	RLO(NE1,L)=(UO(NE1)*DX(I,L)+VO(NE1)*DY(I,L))/DS(I,L) 
	    
	RHOO(NE1)=RHOI  
	
	HO(NE1)=0.5D0*SQRT(UI**2+VI**2)+(CINP**2)/(GAMA-1)
	
	EO(NE1)=HO(NE1)-CINP**2/GAMA  
	
	PO(NE1)=RHOI*CINP**2/GAMA 	     

	
	IF (NW(I)==1) THEN
	
	QLO(NE1,L)=-QLO(I,L)

	RLO(NE1,L)=RLO(I,L)

	UO(NE1)=DS(NE1,L)*(QLO(NE1,L)*DY(NE1,L)+RLO(NE1,L)*DX(NE1,L))/

     &	(DY(NE1,L)**2+DX(NE1,L)**2)

	VO(NE1)=DS(NE1,L)*(-QLO(NE1,L)*DX(NE1,L)+RLO(NE1,L)*DY(NE1,L))

     &	/(DY(NE1,L)**2+DX(NE1,L)**2)	
	
	RHOO(NE1)=RHOO(I)  
	
	HO(NE1)=HO(I)
	
	EO(NE1)=EO(I) 
	
	PO(NE1)=PO(I)
		
	
	END IF  
	
	END IF  

	IF (NE2.EQ.ZERO)      THEN

	DO K=1,3
	
	IF (FACE(I,K).EQ.ZERO) L=K

	END DO
	
	NE2=NE2+NEBC

	DX(NE2,L)=DX(I,L)

	DY(NE2,L)=DY(I,L)

	DS(NE2,L)=DS(I,L)

	ANG(I)=ATAN(YNODE(N1)/(XNODE(N1)-0.5D0))
	
	XKAPA1=0.5D0*0.33D0*SQRT(1-XMIN**2)

      XKAPA2=2*3.14159*10*(1-XMIN**2*(SIN(ANG(I)-2*3.14159/180.))**2)

	XKAPA=XKAPA1/XKAPA2

	UO(NE2)=UI+XKAPA*SIN(ANG(I))

	VO(NE2)=VI-XKAPA*COS(ANG(I))     

	QLO(NE2,L)=(UO(NE2)*DY(I,L)-VO(NE2)*DX(I,L))/DS(I,L) 

      RLO(NE2,L)=(UO(NE2)*DX(I,L)+VO(NE2)*DY(I,L))/DS(I,L)
	     
	RHOO(NE2)=RHOI

	HO(NE2)=0.5D0*SQRT(UI**2+VI**2)+(CINP**2)/(GAMA-1)

	EO(NE2)=HO(NE2)-CINP**2/GAMA  

	PO(NE2)=RHOI*CINP**2/GAMA	

	
	IF (NW(I)==1) THEN
      
	RLO(NE2,L)=RLO(I,L)

	QLO(NE2,L)=-QLO(I,L)

	UO(NE2)=DS(NE2,L)*((QLO(NE2,L)*DY(NE2,L))+(RLO(NE2,L)*DX(NE2,L)))

     &	/(DY(NE2,L)**2+DX(NE2,L)**2)

	VO(NE2)=DS(NE2,L)*((-QLO(NE2,L)*DX(NE2,L))+(RLO(NE2,L)*DY(NE2,L)))

     &	/(DY(NE2,L)**2+DX(NE2,L)**2)
	
	
	RHOO(NE2)=RHOO(I)  
	
	HO(NE2)=HO(I)
	
	EO(NE2)=EO(I)  
	
	PO(NE2)=PO(I)
	

	END IF
      
	END IF

	IF (NE3.EQ.ZERO)     THEN

	DO K=1,3
	
	IF (FACE(I,K).EQ.ZERO) L=K

	END DO
	
	NE3=NE3+NEBC

	DX(NE3,L)=DX(I,L)

	DY(NE3,L)=DY(I,L)

	DS(NE3,L)=DS(I,L) 

      ANG(I)=ATAN(YNODE(N1)/(XNODE(N1)-0.5D0))
	
	XKAPA1=0.5D0*0.33D0*SQRT(1-XMIN**2)

      XKAPA2=2*3.14159*10*(1-XMIN**2*(SIN(ANG(I)-2*3.14159/180.))**2)

	XKAPA=XKAPA1/XKAPA2

	UO(NE3)=UI+XKAPA*SIN(ANG(I))

	VO(NE3)=VI-XKAPA*COS(ANG(I))      
	
	QLO(NE3,L)=(UO(NE3)*DY(I,L)-VO(NE3)*DX(I,L))/DS(I,L) 	

	RLO(NE3,L)=(UO(NE3)*DX(I,L)+VO(NE3)*DY(I,L))/DS(I,L)

	RHOO(NE3)=RHOI
	
	HO(NE3)=0.5D0*SQRT(UI**2+VI**2)+(CINP**2)/(GAMA-1)

	EO(NE3)=HO(NE3)-CINP**2/GAMA  
	
	PO(NE3)=RHOI*CINP**2/GAMA	

	IF (NW(I)==1) THEN	
	
	RLO(NE3,L)=RLO(I,L)

      QLO(NE3,L)=-QLO(I,L)	
	
	UO(NE3)=DS(NE3,L)*((QLO(NE3,L)*DY(NE3,L))+(RLO(NE3,L)*DX(NE3,L)))

     &	/(DY(NE3,L)**2+DX(NE3,L)**2)

	VO(NE3)=DS(NE3,L)*((-QLO(NE3,L)*DX(NE3,L))+(RLO(NE3,L)*DY(NE3,L)))

     &	/(DY(NE3,L)**2+DX(NE3,L)**2)
     
      RHOO(NE3)=RHOO(I)  
	
	HO(NE3)=HO(I)
	
	EO(NE3)=EO(I)  
	
	PO(NE3)=PO(I)
	
     
	END IF

	END IF

	! ROE'S AVERAGING

      RHOFH1=SQRT(RHOO(I)*RHOO(NE1))

	RHOFH2=SQRT(RHOO(I)*RHOO(NE2))

	RHOFH3=SQRT(RHOO(I)*RHOO(NE3))


	RHOS1=(SQRT(RHOO(I))+SQRT(RHOO(NE1)))

	RHOS2=(SQRT(RHOO(I))+SQRT(RHOO(NE2)))

	RHOS3=(SQRT(RHOO(I))+SQRT(RHOO(NE3)))



	UFHAT1=(UO(I)*SQRT(RHOO(I))+UO(NE1)*SQRT(RHOO(NE1)))/RHOS1

	UFHAT2=(UO(I)*SQRT(RHOO(I))+UO(NE2)*SQRT(RHOO(NE2)))/RHOS2

      UFHAT3=(UO(I)*SQRT(RHOO(I))+UO(NE3)*SQRT(RHOO(NE3)))/RHOS3

	
	VFHAT1=(VO(I)*SQRT(RHOO(I))+VO(NE1)*SQRT(RHOO(NE1)))/RHOS1

	VFHAT2=(VO(I)*SQRT(RHOO(I))+VO(NE2)*SQRT(RHOO(NE2)))/RHOS2

	VFHAT3=(VO(I)*SQRT(RHOO(I))+VO(NE3)*SQRT(RHOO(NE3)))/RHOS3



	HFHAT1=(HO(I)*SQRT(RHOO(I))+HO(NE1)*SQRT(RHOO(NE1)))/RHOS1

      HFHAT2=(HO(I)*SQRT(RHOO(I))+HO(NE2)*SQRT(RHOO(NE2)))/RHOS2

	HFHAT3=(HO(I)*SQRT(RHOO(I))+HO(NE3)*SQRT(RHOO(NE3)))/RHOS3


	CHAT1=SQRT((GAMA-1)*(HFHAT1-0.5D0*(UFHAT1**2+VFHAT1**2)))

	CHAT2=SQRT((GAMA-1)*(HFHAT2-0.5D0*(UFHAT2**2+VFHAT2**2)))

      CHAT3=SQRT((GAMA-1)*(HFHAT3-0.5D0*(UFHAT3**2+VFHAT3**2)))     


	QFHAT1=(UFHAT1*DY(I,1)-VFHAT1*DX(I,1))/DS(I,1)


	QFHAT2=(UFHAT2*DY(I,2)-VFHAT2*DX(I,2))/DS(I,2)

	
	QFHAT3=(UFHAT3*DY(I,3)-VFHAT3*DX(I,3))/DS(I,3)

	

      RFHAT1=(UFHAT1*DX(I,1)+VFHAT1*DY(I,1))/DS(I,1)


	RFHAT2=(UFHAT2*DX(I,2)+VFHAT2*DY(I,2))/DS(I,2)


	RFHAT3=(UFHAT3*DX(I,3)+VFHAT3*DY(I,3))/DS(I,3)
      
      ! 1st SIDE

	DO K=1,3
	
	IF (FACE(I,K).EQ.ZERO) L1=K	

	END DO

	IF (NE1.GE.NEBC) THEN 

	M1=L1
	
	GOTO 12

	END IF

	DO K=1,3

	IF (FACE(NE1,K).EQ.I) L1=K

	IF (FACE(I,K).EQ.NE1) M1=K

	END DO

12	K1=-1

	IF (NE1.EQ.NEBC) K1=1

	DP1=PO(NE1)-PO(I)

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

      
	! 2nd SIDE

	DO K=1,3
	
	IF (FACE(I,K).EQ.ZERO) L2=K	

	END DO

	IF (NE2.GE.NEBC) THEN
	
	M2=L2

	GOTO 14

	END IF

	DO K=1,3

	IF (FACE(NE2,K).EQ.I) L2=K

	IF (FACE(I,K).EQ.NE2) M2=K

	END DO


14	K2=-1

	IF (NE2.EQ.NEBC) K2=1

      DP2=PO(NE2)-PO(I)
	
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

 
      !3rd SIDE

	DO K=1,3
	
	IF (FACE(I,K).EQ.ZERO) L3=K	

	END DO

	IF (NE3.GE.NEBC) THEN

	M3=L3
	
	GOTO 16

	END IF

	DO K=1,3

	IF (FACE(NE3,K).EQ.I) L3=K

	IF (FACE(I,K).EQ.NE3) M3=K

	END DO

16	K3=-1

      IF (NE3.EQ.NEBC) K3=1

      DP3=PO(NE3)-PO(I)
	
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

	DEME=(ABS(QFHAT1)+CHAT1)*DS(I,1)+(ABS(QFHAT2)+CHAT2)*DS(I,2)
     &	+(ABS(QFHAT3)+CHAT3)*DS(I,3)

	DT=.5*AREA(I)/DEME
		
       ! ROE

	PHI1=1./2.*(RHOO(I)*QLO(I,1)+(-K1)*RHOO(NE1)*QLO(NE1,L1))

     &	-1./2.*(F11+F12+F13+F14)

	PHI2=1./2.*(RHOO(I)*QLO(I,2)+(-K2)*RHOO(NE2)*QLO(NE2,L2))

     &	-1./2.*(F21+F22+F23+F24)

	PHI3=1./2.*(RHOO(I)*QLO(I,3)+(-K3)*RHOO(NE3)*QLO(NE3,L3))

     &	-1./2.*(F31+F32+F33+F34) 

	DF(I)=PHI1*DS(I,1)+PHI2*DS(I,2)+PHI3*DS(I,3) 
	             
	RESC1(I)=-DF(I)/AREA(I)

	!! UPDATE CONTINUITY

	RHOO(I)=RHOOZ(I)+ALFA(J+1)*DT*RESC1(I)
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

    	
	SX1=1./2.*(RHOO(I)*QLO(I,1)*UO(I)+PO(I)*DY(I,1)/DS(I,1)+(-K1)*
	
     &(RHOO(NE1)*QLO(NE1,L1)*UO(NE1)+PO(NE1)*DY(NE1,L1)/DS(NE1,L1)))

     &-1./2.*(F41+F42+F43+F44)

	SX2=1./2.*(RHOO(I)*QLO(I,2)*UO(I)+PO(I)*DY(I,2)/DS(I,2)+(-K2)*
	
     &(RHOO(NE2)*QLO(NE2,L2)*UO(NE2)+PO(NE2)*DY(NE2,L2)/DS(NE2,L2)))

     &-1./2.*(F51+F52+F53+F54) 

	SX3=1./2.*(RHOO(I)*QLO(I,3)*UO(I)+PO(I)*DY(I,3)/DS(I,3)+(-K3)*
	
     &(RHOO(NE3)*QLO(NE3,L3)*UO(NE3)+PO(NE3)*DY(NE3,L3)/DS(NE3,L3)))

     &-1./2.*(F61+F62+F63+F64)

	FINVT(I)=SX1*DS(I,1)+SX2*DS(I,2)+SX3*DS(I,3)

	! CALCULATE FLUX RESIDUAL

      RESMT1(I)=-FINVT(I)/AREA(I)

	!! UPDATE T-MOM

	FUO(I)=FUOZ(I)+ALFA(J+1)*DT*RESMT1(I)

	UO(I)=FUO(I)/RHOO(I)
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


	DR1=1./2.*(RHOO(I)*QLO(I,1)*VO(I)-PO(I)*DX(I,1)/DS(I,1)+(-K1)*
	
     &(RHOO(NE1)*QLO(NE1,L1)*VO(NE1)-PO(NE1)*DX(NE1,L1)/DS(NE1,L1)))
	
     &-1./2.*(F71+F72+F73+F74)

	DR2=1./2.*(RHOO(I)*QLO(I,2)*VO(I)-PO(I)*DX(I,2)/DS(I,2)+(-K2)*
	
     &(RHOO(NE2)*QLO(NE2,L2)*VO(NE2)-PO(NE2)*DX(NE2,L2)/DS(NE2,L2)))
     
     &-1./2.*(F81+F82+F83+F84)

	DR3=1./2.*(RHOO(I)*QLO(I,3)*VO(I)-PO(I)*DX(I,3)/DS(I,3)+(-K3)*
	
     &(RHOO(NE3)*QLO(NE3,L3)*VO(NE3)-PO(NE3)*DX(NE3,L3)/DS(NE3,L3)))
     
     &-1./2.*(F91+F92+F93+F94)
     		
     		
	FINVR(I)=DR1*DS(I,1)+DR2*DS(I,2)+DR3*DS(I,3)
    
	! CALCULATE FLUX RESIDUAL

	RESMR1(I)=-FINVR(I)/AREA(I)

	!! UPDATE R-MOM

	FVO(I)=FVOZ(I)+ALFA(J+1)*DT*RESMR1(I)

	VO(I)=FVO(I)/RHOO(I)
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
		

	DEP1=1./2.*(RHOO(I)*QLO(I,1)*HO(I)+(-K1)*(RHOO(NE1)*QLO(NE1,L1)
	
     &*HO(NE1)))-1./2.*(F101+F102+F103+F104)

	DEP2=1./2.*(RHOO(I)*QLO(I,2)*HO(I)+(-K2)*(RHOO(NE2)*QLO(NE2,L2)
	
     &*HO(NE2)))-1./2.*(F111+F112+F113+F114)

	DEP3=1./2.*(RHOO(I)*QLO(I,3)*HO(I)+(-K3)*(RHOO(NE3)*QLO(NE3,L3)

     &*HO(NE3)))-1./2.*(F121+F122+F123+F124)
     		
	FINVE(I)=DEP1*DS(I,1)+DEP2*DS(I,2)+DEP3*DS(I,3)
      
	!! UPDATE ENERGY 

	RESE1(I)=-FINVE(I)/AREA(I)     

	! UPDATE ENERGY

      FEO(I)=FEOZ(I)+ALFA(J+1)*DT*RESE1(I)

	EO(I)=FEO(I)/RHOO(I)  

	END DO

	END DO
     	 
      DO I=1,NE

      RHO(I)=RHOO(I)
	
	U(I)=FUO(I)/RHOO(I)
      
	V(I)=FVO(I)/RHOO(I)

	E(I)=FEO(I)/RHOO(I)

	END DO


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

	XM(I)=U(I)/C(I)
	
	DC=(GAMA-1)/(2*GAMA)*(U(I)**2+V(I)**2)	
	
	H(I)=GAMA*(E(I)-DC)

	HO(I)=H(I)      
      
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

	RESIDUAL=MAX(RESCMAX,RESEMAX,RESMRMAX,RESMTMAX)	

C      RESIDUAL=MAX(E11,E2,E3,E4)	
	
	WRITE (*,*)  RESIDUAL 

	WRITE (2,*)  MC,RESIDUAL  
	
	MC=MC+1 
	
	TIME=TIME+DT

	END SUBROUTINE
!****************************************************************************
	Subroutine Output1
		integer :: i,j
		
		open(45,file='output1.plt')

          TY=0

		write(45,*),'VARIABLES="X" "Y" "TY"'
		write(45,*),'ZONE F=FEPOINT,ET=Triangle,N=', NN, ',E=', NE
   
		do I=1,NN
		write(45,*) Xnode(I),Ynode(I),TY(I)
		end do

		do I=1,NE
		write(45,*) NodeNum(I,1),NodeNum(I,2),NodeNum(I,3)
		end do

	end Subroutine 
!****************************************************************************
	Subroutine Output
		integer :: i,j
		
		open(5,file='output.plt')
	    OPEN(12,file='output.DAT')

		write(5,*),'VARIABLES="X" "Y" "U" "V" "P" "MACH" "CPP" "S"'
		write(5,*),'ZONE F=FEPOINT,ET=Triangle,N=', NN, ',E=', NE

	    
      UU=0
	VV=0

	PP=0
	CC=0

	
	DO I=1,NE

      K1=NodeNum(I,1)
      K2=NodeNum(I,2)
      K3=NodeNum(I,3)

      NE1=NEIB(I,1)
      NE2=NEIB(I,2)
      NE3=NEIB(I,3)	

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


	DOM1=AREA(I)
	DOM2=AREA(I)
	DOM3=AREA(I)


	DO J=1,NE

	IF (I.NE.J) THEN

	IF (NodeNum(J,1).EQ.K1.OR.NodeNum(J,2).EQ.K1     
     &	.OR.NodeNum(J,3).EQ.K1) THEN

	UU(K1)=UU(K1)+U(J)*AREA(J)
      VV(K1)=VV(K1)+V(J)*AREA(J)
	PP(K1)=PP(K1)+P(J)*AREA(J)
	CC(K1)=CC(K1)+C(J)*AREA(J)
	RR(K1)=RR(K1)+RHO(J)*AREA(J)
	
	DOM1=DOM1+AREA(J)

	END IF

      IF (NodeNum(J,1).EQ.K2.OR.NodeNum(J,2).EQ.K2     
     &	.OR.NodeNum(J,3).EQ.K2) THEN

      UU(K2)=UU(K2)+U(J)*AREA(J)
      VV(K2)=VV(K2)+V(J)*AREA(J)
	PP(K2)=PP(K2)+P(J)*AREA(J)
	CC(K2)=CC(K2)+C(J)*AREA(J)
	RR(K2)=RR(K2)+RHO(J)*AREA(J)
	
	DOM2=DOM2+AREA(J)

	END IF

	IF (NodeNum(J,1).EQ.K3.OR.NodeNum(J,2).EQ.K3     
     &	.OR.NodeNum(J,3).EQ.K3) THEN

	UU(K3)=UU(K3)+U(J)*AREA(J)
      VV(K3)=VV(K3)+V(J)*AREA(J)
	PP(K3)=PP(K3)+P(J)*AREA(J)
	CC(K3)=CC(K3)+C(J)*AREA(J)
	RR(K3)=RR(K3)+RHO(J)*AREA(J)

	DOM3=DOM3+AREA(J)

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
	
	END DO

		
		DO I=1,NN

          CPP(I)=(PP(I)-CINP**2/GAMA)/(0.5D0*RHOI*UI**2) 
		SO(I)=(PP(I)/RR(I)**GAMA)-1
		  
		write(5,*) Xnode(I),Ynode(I),UU(I),VV(I),PP(I),
     &		SQRT(UU(I)**2+VV(I)**2)/CC(I),CPP(I),SO(I)
		end do

		do I=1,NE
		write(5,*) NodeNum(I,1),NodeNum(I,2),NodeNum(I,3)
		end do

	    Do I=1,NE
		write(12,*) SO(I),CPP(I) !,P(I)
		end do

	 !CALCULATE FORCES

	   FX=0
	   FY=0
	   NM1=-1

	
	   DO I=1,19,2

         IF (NW(I).EQ.1) NM1(I)=1
	    
	   END DO

	   
	   DO I=1,NE
         
	   IF (NW(I).EQ.1) THEN	   

	   DO K=1,3
	  
	   IF (NEIB(I,K).EQ.ZERO) THEN

	   FX=FX+NM1(I)*P(I)*DY(I,K)

	   FY=FY+NM1(I)*P(I)*DX(I,K)
	   
	   WRITE(15,*) I,P(I)*DY(I,K),P(I)*DX(I,K)       
   
	   END IF 

	   END DO
	   
	   END IF	   
	   
	   END DO
	   
	   WRITE(15,*) -FX*2,FY*2	

	end Subroutine Output
!****************************************************************************
	End program 

