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
	PARAMETER (NEBC=2000)
	PARAMETER (MM=4)   !R-K INTEGRATION STEPS

	implicit DoublePrecision  (A-H,O-Z)

	! Variables

	Integer         :: NNode,NBound,Neib(NE,3)

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
		CALL BC
		CALL Initialize	    
		CALL SetNeighbours
C	    CALL Output
10		CALL FLUX
		CALL UPDATE		
		CALL CALC
	    CALL CONVERGENCE 	
	    IF (MC.EQ.1000) GOTO 20
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

	    PRA=72.D-2      ! PRANDTL

	    UI=1.
		
		RHOI=1.D0	

	    TIME=0

		XMIN=1.2D0

	    CINP=1./XMIN
		
		PIN=CINP**2/GAMA

	    HIN=0.5D0+1./((GAMA-1)*XMIN**2)

	    EIN=HIN-CINP**2/GAMA
		
		TI=PIN/(R*RHOI)	      	    		       	      

	    DT=7.D-6

	    EPS=1.D-6 		

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

	    DO K=1,NE

		N1=NodeNum(K,1)
		N2=NodeNum(K,2)
		N3=NodeNum(K,3)

	    AREA(K)=(XNode(N2)*YNode(N3)-XNode(N3)*YNode(N2)+
     &	XNode(N3)*YNode(N1)-XNode(N1)*YNode(N3)+
     &    XNode(N1)*YNode(N2)-XNode(N2)*YNode(N1))/2.
	    
		DX(K,1)=(XNode(N1)-XNode(N2))
	    DX(K,2)=(XNode(N2)-XNode(N3))
	    DX(K,3)=(XNode(N1)-XNode(N3))

	    DY(K,1)=(YNode(N1)-YNode(N2))
	    DY(K,2)=(YNode(N2)-YNode(N3))
	    DY(K,3)=(YNode(N1)-YNode(N3))
		
		! SIDE LENGTH OF THE ELEMENT
  
	    DS(K,1)=SQRT(DX(K,1)**2+DY(K,1)**2)
	    DS(K,2)=SQRT(DX(K,2)**2+DY(K,2)**2)
	    DS(K,3)=SQRT(DX(K,3)**2+DY(K,3)**2)

	    END DO

	End subroutine
!****************************************************************************
	Subroutine BC
	                 
	    Open(2,FILE='RES.DAT')

		Open(3,file='BC.DAT')
	    
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
			              
			END IF   

	     END DO
	     END DO 
		 
		 ! APPLY B.C VALUES
		 
		 DO I=1,NE
		 
	     IF (BOUND(K,1).NE.ZERO.AND.BOUND(K,2).NE.ZERO) THEN

           XM(I)=XMIN
		 
		 U(I)=1.   

	     V(I)=0.  

	     RHO(I)=RHOI
		  
	     P(I)=CINP**2/GAMA	   

	     H(I)=0.5D0+1./((GAMA-1)*XMIN**2)

	     T(I)=(H(I)-1./2.*(U(I)**2+V(I)**2))/CP  

	     E(I)=0.5D0+1./((GAMA-1)*XMIN**2)

	     C(I)=1./XMIN	     
		 
		 END IF 

		 END DO   
	        		 	
	End subroutine 			
!****************************************************************************
	Subroutine Initialize

	    DO I=1,NE

	   	U(I)=1. 

	    V(I)=0.  

	    RHO(I)=1.
		  
	    P(I)=CINP**2/GAMA	   

	    H(I)=0.5D0+1./((GAMA-1)*XMIN**2)

	    T(I)=(H(I)-1./2.*(U(I)**2+V(I)**2))/CP  

	    E(I)=0.5D0+1./((GAMA-1)*XMIN**2)

	    C(I)=1./XMIN

	    XM(I)=U(I)/C(I)

	    END DO
		
		UO=U	    
	    VO=V	   
	    PO=P	    
	    HO=H	    
	    RHOO(1:NE)=RHO(1:NE)  
	    CO=C	    
	    EO=E
		TOLD=T	    
				
	End subroutine 
!****************************************************************************
	subroutine SetNeighbours

		integer :: i,j,N1,N2,N3,N4

		Neib=0

		Do i=1,NE

			N1=NodeNum(i,1)
			N2=NodeNum(i,2)
			N3=NodeNum(i,3)			

c	if (((Bound(N1)==0).or.(Bound(N2))==0).and.(Neib(i,1)==0)) then
		
			Do J=1,NE

			IF (I.NE.J) THEN

			if ((N1==NodeNum(j,1)).and.(N2==NodeNum(j,3)))then
			Neib(i,1)=j
			Neib(j,3)=i
			else if ((N1==NodeNum(j,2)).and.(N2==NodeNum(j,1)))then
			Neib(i,1)=j
			Neib(j,1)=i
			else if ((N1==NodeNum(j,3)).and.(N2==NodeNum(j,2)))then
			Neib(i,1)=j
			Neib(j,2)=i		
			end if
			end if
			end do
C			end if

C	if (((Bound(N2)==0).or.(Bound(N3))==0).and.(Neib(i,2)==0)) then

			Do J=1,NE

			IF (I.NE.J) THEN

			if ((N2==NodeNum(j,1)).and.(N3==NodeNum(j,3)))then
			Neib(i,2)=j
			Neib(j,3)=i
			else if ((N2==NodeNum(j,2)).and.(N3==NodeNum(j,1)))then
			Neib(i,2)=j
			Neib(j,1)=i
			else if ((N2==NodeNum(j,3)).and.(N3==NodeNum(j,2)))then
			Neib(i,2)=j
			Neib(j,2)=i

			end if
			end if
			end do
C			end if

C	if (((Bound(N3)==0).or.(Bound(N4))==0).and.(Neib(i,3)==0)) then
		   
		   Do J=1,NE

		   IF (I.NE.J) THEN

		   if ((N3==NodeNum(j,1)).and.(N1==NodeNum(j,3)))then
		   Neib(i,3)=j
		   Neib(j,3)=i
		   else if ((N3==NodeNum(j,2)).and.(N1==NodeNum(j,1)))then
		   Neib(i,3)=j
		   Neib(j,1)=i
		   else if ((N3==NodeNum(j,3)).and.(N1==NodeNum(j,2)))then
		   Neib(i,3)=j
		   Neib(j,2)=i
					
		   end if
		   end if
		   end do
C		   end if

	       END DO		

	end subroutine SetNeighbours
!****************************************************************************
      Subroutine FLUX       

	EADD=EPS
	
	DO I=1,NE

	FUO(I)=RHOO(I)*UO(I)

	FVO(I)=RHOO(I)*VO(I)

	HO(I)=EO(I)+PO(I)/RHOO(I)

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
		
C	DO J=1,MM

	DO I=1,NE

	NE1=NEIB(I,1)
	NE2=NEIB(I,2)
	NE3=NEIB(I,3)

	IF (NE1.EQ.ZERO)       THEN 
	
	NE1=NE1+1000
      
	RHOO(NE1)=RHOO(I) 
	
	HO(NE1)=HO(I) 

      QLO(NE1,3)=QLO(I,1) ! PLUS FOR INFINITY

      RLO(NE1,3)=RLO(I,1)
	
	PO(NE1)=PO(I)

	TOLD(NE1)=TOLD(I)

      DX(NE1,3)=DX(I,1)

	DY(NE1,3)=DY(I,1)

	DS(NE1,3)=DS(I,1)
	
	END IF  

	IF (NE2.EQ.ZERO)      THEN
	
	NE2=NE2+1000
      
	RHOO(NE2)=RHOO(I)

	HO(NE2)=HO(I) 

      QLO(NE2,1)=QLO(I,2) 

      RLO(NE2,1)=RLO(I,2)

	PO(NE2)=PO(I)

	TOLD(NE2)=TOLD(I)

	DX(NE2,1)=DX(I,2)

	DY(NE2,1)=DY(I,2)

	DS(NE2,1)=DS(I,2)

	END IF

	IF (NE3.EQ.ZERO)     THEN
	
	NE3=NE3+1000	

	RHOO(NE3)=RHOO(I)
	
	HO(NE3)=HO(I) 	

	QLO(NE3,2)=QLO(I,3)	

	RLO(NE3,2)=RLO(I,3)    

      PO(NE3)=PO(I)

	TOLD(NE3)=TOLD(I)

	DX(NE3,2)=DX(I,3)

	DY(NE3,2)=DY(I,3)

	DS(NE3,2)=DS(I,3)

	END IF

	! ROE'S AVERAGING

      RHOFH1=SQRT(RHOO(I)*RHOO(NE1))

	RHOFH2=SQRT(RHOO(I)*RHOO(NE2))

	RHOFH3=SQRT(RHOO(I)*RHOO(NE3))


	RHOS1=(SQRT(RHOO(I))+SQRT(RHOO(NE1)))

	RHOS2=(SQRT(RHOO(I))+SQRT(RHOO(NE2)))

	RHOS3=(SQRT(RHOO(I))+SQRT(RHOO(NE3)))


	PFHAT1=(PO(I)*SQRT(RHOO(I))+PO(NE1)*SQRT(RHOO(NE1)))/RHOS1

	PFHAT2=(PO(I)*SQRT(RHOO(I))+PO(NE2)*SQRT(RHOO(NE2)))/RHOS2

	PFHAT3=(PO(I)*SQRT(RHOO(I))+PO(NE3)*SQRT(RHOO(NE3)))/RHOS3


	CHAT1=(PFHAT1*GAMA/RHOFH1)**0.5D0

	CHAT2=(PFHAT2*GAMA/RHOFH2)**0.5D0

      CHAT3=(PFHAT3*GAMA/RHOFH3)**0.5D0

	
      UFHAT1=(UO(I)*SQRT(RHOO(I))+UO(NE1)*SQRT(RHOO(NE1)))/RHOS1

	UFHAT2=(UO(I)*SQRT(RHOO(I))+UO(NE2)*SQRT(RHOO(NE2)))/RHOS2

      UFHAT3=(UO(I)*SQRT(RHOO(I))+UO(NE3)*SQRT(RHOO(NE3)))/RHOS3

	
	VFHAT1=(VO(I)*SQRT(RHOO(I))+VO(NE1)*SQRT(RHOO(NE1)))/RHOS1

	VFHAT2=(VO(I)*SQRT(RHOO(I))+VO(NE2)*SQRT(RHOO(NE2)))/RHOS2

	VFHAT3=(VO(I)*SQRT(RHOO(I))+VO(NE3)*SQRT(RHOO(NE3)))/RHOS3

	
	QFHAT1=(QLO(I,1)*SQRT(RHOO(I))+QLO(NE1,1)*SQRT(RHOO(NE1)))/RHOS1
	
	!(UFHAT1*DY(I,1)-VFHAT1*DX(I,1))/DS(I,1)

	QFHAT2=(QLO(I,2)*SQRT(RHOO(I))+QLO(NE2,2)*SQRT(RHOO(NE2)))/RHOS2
	
	!(UFHAT2*DY(I,2)-VFHAT2*DX(I,2))/DS(I,2)

      QFHAT3=(QLO(I,3)*SQRT(RHOO(I))+QLO(NE3,3)*SQRT(RHOO(NE3)))/RHOS3
	
	!(UFHAT3*DY(I,3)-VFHAT3*DX(I,3))/DS(I,3)

	
	RFHAT1=(RLO(I,1)*SQRT(RHOO(I))+RLO(NE1,1)*SQRT(RHOO(NE1)))/RHOS1

	!(UFHAT1*DX(I,1)+VFHAT1*DY(I,1))/DS(I,1)

	RFHAT2=(RLO(I,2)*SQRT(RHOO(I))+RLO(NE2,2)*SQRT(RHOO(NE2)))/RHOS2

	!(UFHAT2*DX(I,2)+VFHAT2*DY(I,2))/DS(I,2)

      RFHAT3=(RLO(I,3)*SQRT(RHOO(I))+RLO(NE3,3)*SQRT(RHOO(NE3)))/RHOS3

	!(UFHAT3*DX(I,3)+VFHAT3*DY(I,3))/DS(I,3)


      HFHAT1=(HO(I)*SQRT(RHOO(I))+HO(NE1)*SQRT(RHOO(NE1)))/RHOS1

      HFHAT2=(HO(I)*SQRT(RHOO(I))+HO(NE2)*SQRT(RHOO(NE2)))/RHOS2

	HFHAT3=(HO(I)*SQRT(RHOO(I))+HO(NE3)*SQRT(RHOO(NE3)))/RHOS3


	TFHAT1=(TOLD(I)*SQRT(RHOO(I))+TOLD(NE1)*SQRT(RHOO(NE1)))/RHOS1

	TFHAT2=(TOLD(I)*SQRT(RHOO(I))+TOLD(NE2)*SQRT(RHOO(NE2)))/RHOS2

	TFHAT3=(TOLD(I)*SQRT(RHOO(I))+TOLD(NE3)*SQRT(RHOO(NE3)))/RHOS3


      ! 1st SIDE

	DP1=PO(NE1)-PO(I)

	DQ1=QLO(NE1,3)-QLO(I,1)

	DRHO1=RHOO(NE1)-RHOO(I)

	DR1=RLO(NE1,3)-RLO(I,1)

	
	V1=(DP1-RHOFH1*CHAT1*DQ1)/(2*CHAT1**2)

	COF1=ABS(QFHAT1-CHAT1)*V1


	V2=RHOFH1*DR1/CHAT1

	COF2=ABS(QFHAT1)*V2


	V3=DRHO1-DP1/CHAT1**2

	COF3=ABS(QFHAT1)*V3


	V4=(DP1+RHOFH1*CHAT1*DQ1)/(2*CHAT1**2)

	COF4=ABS(QFHAT1+CHAT1)*V4

      
	! 2nd SIDE

	DP2=PO(NE2)-PO(I)
	
	DQ2=QLO(NE2,1)-QLO(I,2)

	DRHO2=RHOO(NE2)-RHOO(I)

	DR2=RLO(NE2,1)-RLO(I,2)



	V5=(DP2-RHOFH2*CHAT2*DQ2)/(2*CHAT2**2)

	COF5=ABS(QFHAT2-CHAT2)*V5

      V6=RHOFH2*DR2/CHAT2  

	COF6=ABS(QFHAT2)*V6

	V7=DRHO2-DP2/CHAT2**2
	
	COF7=ABS(QFHAT2)*V7

	V8=(DP2+RHOFH2*CHAT2*DQ2)/(2*CHAT2**2)

	COF8=ABS(QFHAT2+CHAT2)*V8

 
      !3rd SIDE

	DP3=PO(NE3)-PO(I)
	
	DQ3=QLO(NE3,2)-QLO(I,3)

	DRHO3=RHOO(NE3)-RHOO(I)

	DR3=RLO(NE3,2)-RLO(I,3)



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

	PHI1=1./2.*(RHOO(I)*QLO(I,1)+RHOO(NE1)*QLO(NE1,1))

     &	-1./2.*(F11+F12+F13+F14)

	PHI2=1./2.*(RHOO(I)*QLO(I,2)+RHOO(NE2)*QLO(NE2,2))

     &	-1./2.*(F21+F22+F23+F24)

	PHI3=1./2.*(RHOO(I)*QLO(I,3)+RHOO(NE3)*QLO(NE3,3))

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
	
     &  +RHOO(NE1)*QLO(NE1,1)*UO(NE1)+PO(NE1)*DY(NE1,1)/DS(NE1,1))

     &  -1./2.*(F41+F42+F43+F44)

	SX2=1./2.*(RHOO(I)*QLO(I,2)*UO(I)+PO(I)*DY(I,2)/DS(I,2)
	
     &   +RHOO(NE2)*QLO(NE2,2)*UO(NE2)+PO(NE2)*DY(NE2,2)/DS(NE2,2))

     &   -1./2.*(F51+F52+F53+F54) 

	SX3=1./2.*(RHOO(I)*QLO(I,3)*UO(I)+PO(I)*DY(I,3)/DS(I,3)
	
     &   +RHOO(NE3)*QLO(NE3,3)*UO(NE3)+PO(NE3)*DY(NE3,3)/DS(NE3,3))

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
	
     &   +RHOO(NE1)*QLO(NE1,1)*VO(NE1)-PO(NE1)*DX(NE1,1)/DS(NE1,1))
	
     &	-1./2.*(F71+F72+F73+F74)

	DR2=1./2.*(RHOO(I)*QLO(I,2)*VO(I)-PO(I)*DX(I,2)/DS(I,2)
	
     &   +RHOO(NE2)*QLO(NE2,2)*VO(NE2)-PO(NE2)*DX(NE2,2)/DS(NE2,2))
     
     &    -1./2.*(F81+F82+F83+F84)

	DR3=1./2.*(RHOO(I)*QLO(I,3)*VO(I)-PO(I)*DX(I,3)/DS(I,3)
	
     &   +RHOO(NE3)*QLO(NE3,3)*VO(NE3)-PO(NE3)*DX(NE3,3)/DS(NE3,3))
     
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

	F113=COF7*1./2.*(VFHAT2**2+UFHAT2**2)

	F114=COF8*(HFHAT2+UFHAT2*CHAT2)


	F121=COF9*(HFHAT3-UFHAT3*CHAT3)

	F122=COF10*RFHAT3*CHAT3

	F123=COF11*1./2.*(VFHAT3**2+UFHAT3**2)

	F124=COF12*(HFHAT3+QFHAT3*CHAT3)	
		

	DEP1=1./2.*(RHOO(I)*QLO(I,1)*HO(I)+RHOO(NE1)*QLO(NE1,1)*HO(NE1))
	
     &	-1./2.*(F101+F102+F103+F104)

	DEP2=1./2.*(RHOO(I)*QLO(I,2)*HO(I)+RHOO(NE2)*QLO(NE2,2)*HO(NE2))
	
     &	-1./2.*(F111+F112+F113+F114)

	DEP3=1./2.*(RHOO(I)*QLO(I,3)*HO(I)+RHOO(NE3)*QLO(NE3,3)*HO(NE3))

     &	-1./2.*(F121+F122+F123+F124)
     		
	FINVE(I)=DEP1*DS(I,1)+DEP2*DS(I,2)+DEP3*DS(I,3)
      
      CONTINUE

	!! UPDATE ENERGY 

	RESE1(I)=-FINVE(I)/AREA(I)     

	! UPDATE ENERGY

      FE(I)=FEO(I)+DT*RESE1(I)

	E(I)=FE(I)/RHO(I)  

	END DO

C	END DO
      	
	END SUBROUTINE
!****************************************************************************
!     CALCULATE PARAMETERS
!
!****************************************************************************      
	subroutine UPDATE	
      
	! CALCULATE THE SUM OF THE RESIDUALS

	SUM1=0
	SUM2=0
	SUM3=0
	SUM4=0
		
	DO I=1,NE
	      
	SUM1(I)=1. !RHOO(I)
	SUM2(I)=1. !FTO(I)
	SUM3(I)=1. !FRO(I)
	SUM4(I)=1. !FEO(I)
	  
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
	
      P(I)=(E(I)-(U(I)**2+V(I)**2)/2.)*((GAMA-1)*RHO(I))

	PO(I)=P(I)

	T(I)=(E(I)-(U(I)**2+V(I)**2)/2.)/CV*(UI**2/TI)	 	

	TOLD(I)=T(I)

	C(I)=SQRT(GAMA*R*T(I)*TI)/UI 

	CO(I)=C(I)

	XM(I)=U(I)/C(I)	
	
	H(I)=E(I)+P(I)/RHO(I) !P(I)/(R*RHO(I))*UI**2/TI 

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
	
	RESCMAX=RESC(3)
	
	RESEMAX=RESET(3)

	RESMRMAX=RESMRT(3)

	RESMTMAX=RESMTT(3)
	
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

      E11=SQRT(E11)/(NE-2)
	E2=SQRT(E2)/(NE-2)
	E3=SQRT(E3)/(NE-2)
	E4=SQRT(E4)/(NE-2)

C	RESIDUAL=MAX(RESCMAX,RESEMAX,RESMRMAX,RESMTMAX)	

      RESIDUAL=MAX(E11,E2,E3,E4)	
	
	WRITE (*,*)  RESIDUAL 

	WRITE (2,*)  MC,RESIDUAL  
	
	MC=MC+1 
	
	TIME=TIME+DT

	END SUBROUTINE
!****************************************************************************
	Subroutine Output
		integer :: i,j
		
		open(5,file='output.plt')

		write(5,*),'VARIABLES="X" "Y" "T" "U" "V" "RHO" "P"'
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
      
	RHOO(NE1)=RHOO(I) 
	
	HO(NE1)=HO(I) 

      QLO(NE1,1)=QLO(I,1) ! PLUS FOR INFINITY

      RLO(NE1,1)=RLO(I,1)
	
	PO(NE1)=PO(I)

	TOLD(NE1)=TOLD(I)

      DX(NE1,1)=DX(I,1)

	DY(NE1,1)=DY(I,1)

	DS(NE1,1)=DS(I,1)
	
	END IF  

	IF (NE2.EQ.ZERO)      THEN
	
	NE2=NE2+1000
      
	RHOO(NE2)=RHOO(I)

	HO(NE2)=HO(I) 

      QLO(NE2,2)=QLO(I,2) 

      RLO(NE2,2)=RLO(I,2)

	PO(NE2)=PO(I)

	TOLD(NE2)=TOLD(I)

	DX(NE2,2)=DX(I,2)

	DY(NE2,2)=DY(I,2)

	DS(NE2,2)=DS(I,2)

	END IF

	IF (NE3.EQ.ZERO)     THEN
	
	NE3=NE3+1000	

	RHOO(NE3)=RHOO(I)
	
	HO(NE3)=HO(I) 	

	QLO(NE3,3)=QLO(I,3)	

	RLO(NE3,3)=RLO(I,3)    

      PO(NE3)=PO(I)

	TOLD(NE3)=TOLD(I)

	DX(NE3,3)=DX(I,3)

	DY(NE3,3)=DY(I,3)

	DS(NE3,3)=DS(I,3)

	END IF

	  TT(K1)=((T(I)+T(NE1))/2.)-((T(I)+T(NE2))/2.)+((T(I)+T(NE3))/2.)
	   
	  UU(K1)=((U(I)+U(NE1))/2.)-((U(I)+U(NE2))/2.)+((U(I)+U(NE3))/2.)

	  VV(K1)=((V(I)+V(NE1))/2.)-((V(I)+V(NE2))/2.)+((V(I)+V(NE3))/2.)

	  RR(K1)=((RHO(I)+RHO(NE1))/2.)-((RHO(I)+RHO(NE2))/2.)+
     &		((RHO(I)+RHO(NE3))/2.)

	  PP(K1)=((P(I)+P(NE1))/2.)-((P(I)+P(NE2))/2.)+((P(I)+P(NE3))/2.)


		
        TT(K2)=((T(I)+T(NE1))/2.)+((T(I)+T(NE2))/2.)-((T(I)+T(NE3))/2.)
	    
	  UU(K2)=((U(I)+U(NE1))/2.)+((U(I)+U(NE2))/2.)-((U(I)+U(NE3))/2.)

	  VV(K2)=((V(I)+V(NE1))/2.)+((V(I)+V(NE2))/2.)-((V(I)+V(NE3))/2.)

	  RR(K2)=((RHO(I)+RHO(NE1))/2.)+((RHO(I)+RHO(NE2))/2.)-
     &		((RHO(I)+RHO(NE3))/2.)

        PP(K2)=((P(I)+P(NE1))/2.)+((P(I)+P(NE2))/2.)-((P(I)+P(NE3))/2.)


		
	 TT(K3)=-((T(I)+T(NE1))/2.)+((T(I)+T(NE2))/2.)+((T(I)+T(NE3))/2.)
	    
	 UU(K3)=-((U(I)+U(NE1))/2.)+((U(I)+U(NE2))/2.)+((U(I)+U(NE3))/2.)

	 VV(K3)=-((V(I)+V(NE1))/2.)+((V(I)+V(NE2))/2.)+((V(I)+V(NE3))/2.)

	 RR(K3)=-((RHO(I)+RHO(NE1))/2.)+((RHO(I)+RHO(NE2))/2.)+
     &		((RHO(I)+RHO(NE3))/2.)

	  P(K3)=-((P(I)+P(NE1))/2.)+((P(I)+P(NE2))/2.)+((P(I)+P(NE3))/2.)
		  

		END DO

		
		do I=1,NN
		write(5,*) Xnode(I),Ynode(I),TT(I),UU(I),VV(I),RR(I),PP(I)
		end do

		do I=1,NE
		write(5,*) NodeNum(I,1),NodeNum(I,2),NodeNum(I,3)
		end do

	end Subroutine Output
!****************************************************************************
	End program 

