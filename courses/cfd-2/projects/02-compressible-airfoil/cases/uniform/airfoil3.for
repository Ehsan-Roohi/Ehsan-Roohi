!****************************************************************************
!
!  PROGRAM: AIRFOIL
!
!
!****************************************************************************

	program AIRFOIL
	
	    USE MSIMSL
		USE MSFLIB
	
	PARAMETER (NN=258)
	PARAMETER (NE=464) 
	PARAMETER (NEBC=1100)
	PARAMETER (MM=4)   !R-K INTEGRATION STEPS

	implicit DoublePrecision  (A-H,O-Z)

	! Variables

	Integer         :: NNode,NBound,Neib(NE,3),FACE(NE,3)

	DoublePrecision ::T(NEBC),AK(3,3),AR(NEBC),A(NEBC,NEBC),R,EPS,CINP

	DoublePrecision :: Xnode(NEBC),Ynode(NEBC),H(NEBC),C(NEBC),U(NEBC)
     &,V(NEBC),RHO(NEBC),RHOO(NEBC),P(NEBC),PO(NEBC),ULO(NEBC,3)
     &,VLO(NEBC,3),QLO(NEBC,3),RLO(NEBC,3),E11,E2,E3,E4

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

	DoublePrecision ::AREV(NEBC,NEBC),AREA(NEBC),DS(NEBC,3),ALFA(MM+1)

      DoublePrecision ::H1(NEBC),C1(NEBC),DX(NEBC,3),DY(NEBC,3)

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
10        CALL BC		
          CALL FLUX
		CALL UPDATE		
		CALL CALC
	    CALL BC
	    CALL CONVERGENCE 	
	    IF (MC.EQ.50000) GOTO 20
		If (RESIDUAL.GT.EPS) GOTO 10		
20	    CALL Output

	contains
!****************************************************************************
	subroutine Input
      
	    GAMA=1.4D0

	    R=287.D0

	    ZERO=0

	    ALFA(2)=0.0033

	    ALFA(3)=0.2069

	    ALFA(4)=0.4265

	    ALFA(5)=1.

	    CP=1004.5D0

	    CV=717.5D0

	    UI=1.
		
		RHOI=1.D0	

	    TIME=0

		XMIN=1.2D0

	    CINP=UI/XMIN 
		
		PIN=RHOI*CINP**2/GAMA

	    HIN=0.5D0+1./((GAMA-1)*XMIN**2)

	    EIN=HIN-CINP**2/GAMA
		
		TI=PIN/(R*RHOI)	      	    		       	      

	    DT=4.D-3

	    EPS=1.D-5 		

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

	        WRITE(15,*) K,BOUND(K,1),BOUND(K,2)
			              
			END IF   

	    END DO
	    END DO 

	End subroutine
!****************************************************************************
	Subroutine Initialize

	    DO I=1,NE

	   	U(I)=1. 

	    V(I)=0.  

	    RHO(I)=1.
		  
	    P(I)=RHO(I)*CINP**2/GAMA	   

	    H(I)=0.5D0+1./((GAMA-1)*XMIN**2)

	    E(I)=H(I)-C(I)**2/GAMA

	    C(I)=1./XMIN

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
	Subroutine BC              
	    	
		 ! APPLY B.C VALUES
		 
		 DO I=1,NE
		 
	     IF (BOUND(I,1).NE.ZERO.AND.BOUND(I,2).NE.ZERO) THEN

           XM(I)=XMIN
		 
		 U(I)=1.   

	     V(I)=0.  

	     C(I)=U(I)/XMIN	     

	     RHO(I)=RHOI
		  
	     P(I)=RHO(I)*CINP**2/GAMA	   

	     H(I)=0.5D0+1./((GAMA-1)*XMIN**2)

	     E(I)=H(I)-C(I)**2/GAMA
		 
		 END IF 

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
      Subroutine FLUX       

	EADD=EPS
	
	DO I=1,NE

	FUO(I)=RHOO(I)*UO(I)

	FVO(I)=RHOO(I)*VO(I)

	FEO(I)=RHOO(I)*EO(I)

	END DO

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

	IF (NE1.EQ.ZERO)       THEN 
	
	DO K=1,3
	
	IF (FACE(I,K).EQ.ZERO) L=K

	END DO
	
	NE1=NE1+1000
      
	RHOO(NE1)=RHOO(I) 
	
	HO(NE1)=HO(I) 

      QLO(NE1,L)=-QLO(I,1) 

      RLO(NE1,L)=-RLO(I,1)

	UO(NE1)=UO(I)

	VO(NE1)=VO(I)
	
	PO(NE1)=PO(I)

      DX(NE1,L)=DX(I,1)

	DY(NE1,L)=DY(I,1)

	DS(NE1,L)=DS(I,1)
	
	END IF  

	IF (NE2.EQ.ZERO)      THEN

	DO K=1,3
	
	IF (FACE(I,K).EQ.ZERO) L=K

	END DO
	
	NE2=NE2+1000
      
	RHOO(NE2)=RHOO(I)

	HO(NE2)=HO(I) 

      QLO(NE2,L)=-QLO(I,2) 

      RLO(NE2,L)=-RLO(I,2)

	UO(NE2)=UO(I)

	VO(NE2)=VO(I)

	PO(NE2)=PO(I)

	DX(NE2,L)=DX(I,2)

	DY(NE2,L)=DY(I,2)

	DS(NE2,L)=DS(I,2)

	END IF

	IF (NE3.EQ.ZERO)     THEN

	DO K=1,3
	
	IF (FACE(I,K).EQ.ZERO) L=K

	END DO
	
	NE3=NE3+1000	

	RHOO(NE3)=RHOO(I)
	
	HO(NE3)=HO(I) 	

	QLO(NE3,L)=-QLO(I,3)	

	RLO(NE3,L)=-RLO(I,3)
	
	UO(NE3)=UO(I)

	VO(NE3)=VO(I)    

      PO(NE3)=PO(I)

	DX(NE3,L)=DX(I,3)

	DY(NE3,L)=DY(I,3)

	DS(NE3,L)=DS(I,3)

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



C	PFHAT1=(PO(I)*SQRT(RHOO(I))+PO(NE1)*SQRT(RHOO(NE1)))/RHOS1

C	PFHAT2=(PO(I)*SQRT(RHOO(I))+PO(NE2)*SQRT(RHOO(NE2)))/RHOS2

C	PFHAT3=(PO(I)*SQRT(RHOO(I))+PO(NE3)*SQRT(RHOO(NE3)))/RHOS3


C	CHAT1=(PFHAT1*GAMA/RHOFH1)**0.5D0

C	CHAT2=(PFHAT2*GAMA/RHOFH2)**0.5D0

C     CHAT3=(PFHAT3*GAMA/RHOFH3)**0.5D0

	
      

	
C	QFHAT1=(QLO(I,1)*SQRT(RHOO(I))+QLO(NE1,1)*SQRT(RHOO(NE1)))/RHOS1
	
	QFHAT1=(UFHAT1*DY(I,1)-VFHAT1*DX(I,1))/DS(I,1)

C	QFHAT2=(QLO(I,2)*SQRT(RHOO(I))+QLO(NE2,2)*SQRT(RHOO(NE2)))/RHOS2
	
	QFHAT2=(UFHAT2*DY(I,2)-VFHAT2*DX(I,2))/DS(I,2)

C      QFHAT3=(QLO(I,3)*SQRT(RHOO(I))+QLO(NE3,3)*SQRT(RHOO(NE3)))/RHOS3
	
	QFHAT3=(UFHAT3*DY(I,3)-VFHAT3*DX(I,3))/DS(I,3)

	
C	RFHAT1=(RLO(I,1)*SQRT(RHOO(I))+RLO(NE1,1)*SQRT(RHOO(NE1)))/RHOS1

      RFHAT1=(UFHAT1*DX(I,1)+VFHAT1*DY(I,1))/DS(I,1)

C	RFHAT2=(RLO(I,2)*SQRT(RHOO(I))+RLO(NE2,2)*SQRT(RHOO(NE2)))/RHOS2

	RFHAT2=(UFHAT2*DX(I,2)+VFHAT2*DY(I,2))/DS(I,2)

C      RFHAT3=(RLO(I,3)*SQRT(RHOO(I))+RLO(NE3,3)*SQRT(RHOO(NE3)))/RHOS3

	RFHAT3=(UFHAT3*DX(I,3)+VFHAT3*DY(I,3))/DS(I,3)
      
      ! 1st SIDE

	DO K=1,3
	
	IF (FACE(I,K).EQ.ZERO) L1=K	

	END DO

	IF (NE1.GE.1000) GOTO 12

	DO K=1,3

	IF (FACE(NE1,K).EQ.I) L1=K

	END DO

	
12	DP1=PO(NE1)-PO(I)

	DQ1=-QLO(NE1,L1)-QLO(I,1)

	DRHO1=RHOO(NE1)-RHOO(I)

	DR1=-RLO(NE1,L1)-RLO(I,1)

	
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

	IF (NE2.GE.1000) GOTO 14

	DO K=1,3

	IF (FACE(NE2,K).EQ.I) L2=K

	END DO


14	DP2=PO(NE2)-PO(I)
	
	DQ2=-QLO(NE2,L2)-QLO(I,2)

	DRHO2=RHOO(NE2)-RHOO(I)

	DR2=-RLO(NE2,L2)-RLO(I,2)



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

	IF (NE3.GE.1000) GOTO 16

	DO K=1,3

	IF (FACE(NE3,K).EQ.I) L3=K

	END DO

16	DP3=PO(NE3)-PO(I)
	
	DQ3=-QLO(NE3,L3)-QLO(I,3)

	DRHO3=RHOO(NE3)-RHOO(I)

	DR3=-RLO(NE3,L3)-RLO(I,3)


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
		
       ! ROE

	PHI1=1./2.*(RHOO(I)*QLO(I,1)+RHOO(NE1)*QLO(NE1,L1))

     &	-1./2.*(F11+F12+F13+F14)

	PHI2=1./2.*(RHOO(I)*QLO(I,2)+RHOO(NE2)*QLO(NE2,L2))

     &	-1./2.*(F21+F22+F23+F24)

	PHI3=1./2.*(RHOO(I)*QLO(I,3)+RHOO(NE3)*QLO(NE3,L3))

     &	-1./2.*(F31+F32+F33+F34) 

	DF(I)=PHI1*DS(I,1)+PHI2*DS(I,2)+PHI3*DS(I,3) 
	             
	RESC1(I)=-DF(I)/AREA(I)

	!! UPDATE CONTINUITY

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

    	
	SX1=1./2.*(RHOO(I)*QLO(I,1)*UO(I)+PO(I)*DY(I,1)/DS(I,1)
	
     &  +RHOO(NE1)*QLO(NE1,L1)*UO(NE1)+PO(NE1)*DY(NE1,L1)/DS(NE1,L1))

     &  -1./2.*(F41+F42+F43+F44)

	SX2=1./2.*(RHOO(I)*QLO(I,2)*UO(I)+PO(I)*DY(I,2)/DS(I,2)
	
     &   +RHOO(NE2)*QLO(NE2,L2)*UO(NE2)+PO(NE2)*DY(NE2,L2)/DS(NE2,L2))

     &   -1./2.*(F51+F52+F53+F54) 

	SX3=1./2.*(RHOO(I)*QLO(I,3)*UO(I)+PO(I)*DY(I,3)/DS(I,3)
	
     &   +RHOO(NE3)*QLO(NE3,L3)*UO(NE3)+PO(NE3)*DY(NE3,L3)/DS(NE3,L3))

     &    -1./2.*(F61+F62+F63+F64) 
	
	FINVT(I)=SX1*DS(I,1)+SX2*DS(I,2)+SX3*DS(I,3)

	! CALCULATE FLUX RESIDUAL

      RESMT1(I)=-FINVT(I)/AREA(I)

	!! UPDATE T-MOM

	FU(I)=FUO(I)+DT*RESMT1(I)

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


	DR1=1./2.*(RHOO(I)*QLO(I,1)*VO(I)-PO(I)*DX(I,1)/DS(I,1)
	
     &   +RHOO(NE1)*QLO(NE1,L1)*VO(NE1)-PO(NE1)*DX(NE1,L1)/DS(NE1,L1))
	
     &	-1./2.*(F71+F72+F73+F74)

	DR2=1./2.*(RHOO(I)*QLO(I,2)*VO(I)-PO(I)*DX(I,2)/DS(I,2)
	
     &   +RHOO(NE2)*QLO(NE2,L2)*VO(NE2)-PO(NE2)*DX(NE2,L2)/DS(NE2,L2))
     
     &    -1./2.*(F81+F82+F83+F84)

	DR3=1./2.*(RHOO(I)*QLO(I,3)*VO(I)-PO(I)*DX(I,3)/DS(I,3)
	
     &   +RHOO(NE3)*QLO(NE3,L3)*VO(NE3)-PO(NE3)*DX(NE3,L3)/DS(NE3,L3))
     
     &    -1./2.*(F91+F92+F93+F94)
     		
	FINVR(I)=DR1*DS(I,1)+DR2*DS(I,2)+DR3*DS(I,3)
    
	! CALCULATE FLUX RESIDUAL

	RESMR1(I)=-FINVR(I)/AREA(I)

	!! UPDATE R-MOM

	FV(I)=FVO(I)+DT*RESMR1(I)

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
		

	DEP1=1./2.*(RHOO(I)*QLO(I,1)*HO(I)+RHOO(NE1)*QLO(NE1,L1)*HO(NE1))
	
     &	-1./2.*(F101+F102+F103+F104)

	DEP2=1./2.*(RHOO(I)*QLO(I,2)*HO(I)+RHOO(NE2)*QLO(NE2,L2)*HO(NE2))
	
     &	-1./2.*(F111+F112+F113+F114)

	DEP3=1./2.*(RHOO(I)*QLO(I,3)*HO(I)+RHOO(NE3)*QLO(NE3,L3)*HO(NE3))

     &	-1./2.*(F121+F122+F123+F124)
     		
	FINVE(I)=DEP1*DS(I,1)+DEP2*DS(I,2)+DEP3*DS(I,3)
      
      CONTINUE

	!! UPDATE ENERGY 

	RESE1(I)=-FINVE(I)/AREA(I)     

	! UPDATE ENERGY

      FE(I)=FEO(I)+DT*RESE1(I)

	E(I)=FE(I)/RHO(I)  

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
	Subroutine Output
		integer :: i,j
		
		open(5,file='output.plt')
	    OPEN(12,file='output.DAT')

		write(5,*),'VARIABLES="X" "Y" "U" "V" "RHO" "P"'
		write(5,*),'ZONE F=FEPOINT,ET=Triangle,N=', NN, ',E=', NE

		DO I=1,NE

	    K1=NodeNum(I,1)
	    K2=NodeNum(I,2)
	    K3=NodeNum(I,3)

	    NE1=NEIB(I,1)
	    NE2=NEIB(I,2)
	    NE3=NEIB(I,3)

		IF (NE1.EQ.ZERO)       THEN 
	
	NE1=NE1+1000
      
	RHO(NE1)=RHO(I) 
	
	U(NE1)=U(I)

	V(NE1)=V(I)

	H(NE1)=H(I) 

      P(NE1)=P(I)
	     
	
	END IF  

	IF (NE2.EQ.ZERO)      THEN
	
	NE2=NE2+1000
      
	RHO(NE2)=RHO(I)

	U(NE2)=U(I)

	V(NE2)=V(I)

	H(NE2)=H(I)      

	P(NE2)=P(I)
	

	END IF

	IF (NE3.EQ.ZERO)     THEN
	
	NE3=NE3+1000	

	RHO(NE3)=RHO(I)

	U(NE3)=U(I)

	V(NE3)=V(I)	

      P(NE3)=P(I)	

	END IF
	   
	  UU(K1)=((U(I)+U(NE1))/2.)-((U(I)+U(NE2))/2.)+((U(I)+U(NE3))/2.)

	  VV(K1)=((V(I)+V(NE1))/2.)-((V(I)+V(NE2))/2.)+((V(I)+V(NE3))/2.)

	  RR(K1)=((RHO(I)+RHO(NE1))/2.)-((RHO(I)+RHO(NE2))/2.)+
     &		((RHO(I)+RHO(NE3))/2.)

	  PP(K1)=((P(I)+P(NE1))/2.)-((P(I)+P(NE2))/2.)+((P(I)+P(NE3))/2.)

		    
	  UU(K2)=((U(I)+U(NE1))/2.)+((U(I)+U(NE2))/2.)-((U(I)+U(NE3))/2.)

	  VV(K2)=((V(I)+V(NE1))/2.)+((V(I)+V(NE2))/2.)-((V(I)+V(NE3))/2.)

	  RR(K2)=((RHO(I)+RHO(NE1))/2.)+((RHO(I)+RHO(NE2))/2.)-
     &		((RHO(I)+RHO(NE3))/2.)

        PP(K2)=((P(I)+P(NE1))/2.)+((P(I)+P(NE2))/2.)-((P(I)+P(NE3))/2.)


    
	 UU(K3)=-((U(I)+U(NE1))/2.)+((U(I)+U(NE2))/2.)+((U(I)+U(NE3))/2.)

	 VV(K3)=-((V(I)+V(NE1))/2.)+((V(I)+V(NE2))/2.)+((V(I)+V(NE3))/2.)

	 RR(K3)=-((RHO(I)+RHO(NE1))/2.)+((RHO(I)+RHO(NE2))/2.)+
     &		((RHO(I)+RHO(NE3))/2.)

	 P(K3)=-((P(I)+P(NE1))/2.)+((P(I)+P(NE2))/2.)+((P(I)+P(NE3))/2.)
		  

		END DO

		
		do I=1,NN
		write(5,*) Xnode(I),Ynode(I),UU(I),VV(I),RR(I),PP(I)
		end do

		do I=1,NE
		write(5,*) NodeNum(I,1),NodeNum(I,2),NodeNum(I,3)
		end do

	    Do I=1,NE
		write(12,*) U(I),V(I),RHO(I) !,P(I)
		end do

	end Subroutine Output
!****************************************************************************
	End program 

