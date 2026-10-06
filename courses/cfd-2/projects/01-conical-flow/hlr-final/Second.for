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

      integer,parameter :: N=240	
  
	DoublePrecision :: TETA(N),TETAF(N),RHO(N),RHOO(N),RHOF(N),P(N)	
     &	,PO(N),PF(N),U(N) 

      ! TETA: ANGLE
	! TETAF: FACE ANGLE
	! RHO: DENSITY
	! RHOF: DENSITY AT FACE
	! RHOO: OLD DENSITY
	! P: PRESSURE
	! POLD: OLD PRESSURE
	! PF: FACE PRESSURE

	DoublePrecision :: UR(N),URF(N),UT(N),UTF(N),URO(N),UTO(N),FT(N)    
     &                   ,FTF(N),FR(N),FRF(N),FTO(N),FRO(N) 

	! UR:  RADIAL VELOCITY
	! UT:  TENGENTIAL VELOCITY
	! URO: OLD RADIAL VELOCITY
	! UTO: OLD TANGANTIAL VELOCITY
	! FT:  TANGANTIAL FLUX (RHO * U_ TETA)
	! FR:  RADIAL FLUX: (RHO * U_RADIAL)
	! FTO: OLD TANGANTIAL FLUX  
	! FRO: OLD RADIAL FLUX

	DoublePrecision :: SRT(N),SRR(N),SPP(N),STT(N)

	! SRT: TAU_ (R,TETA)
	! SRR: TAU_ (R,R)
	! SPP: TAU_(PHI,PHI)
	! STT: TAU_(TETA,TETA)
     
	DoublePrecision :: E(N),EF(N),EO(N),H(N),HF(N),HO(N),FE(N),FEO(N)	
     &              	,T(N),TF(N),TOLD(N)
	     
	! E: INTERNAL ENERGY
	! EO: OLD INTERNAL ENERGY
	! FE: FLUX OF INTERNAL ENERGY (RHO * E)
	! FEO: OLD FLUX OF INTERNAL ENERGY
	! H: ENTHALPY
	! HO: OLD ENTHALPY
	! T: TEMPERATURE
	! TOLD: OLD TEMP

	DoublePrecision ::RESC1(N),RESE1(N),RESMR1(N),RESMT1(N),XMU(N)
	DoublePrecision ::RESC2(N),RESE2(N),RESMR2(N),RESMT2(N)
      DoublePrecision ::RESC(N),RESE(N),RESMR(N),RESMT(N)	
	DoublePrecision :: DUTDT1(N),DUTDT2(N),DURDT1(N),DURDT2(N),
     &UR1(N),UR2(N),UT1(N),UT2(N),STT1(N),STT2(N),SRR1(N),SRR2(N),     
     &SPP1(N),SPP2(N),SRT1(N),SRT2(N),DUTDTETA(N),DURDTETA(N),DF(N)
     &,FVISCT(N),FVISCR(N),FVISCE(N),FMTV(N),FMRV(N),FEV(N),TIME    
	 
	

      DoublePrecision :: FINVT(N),FINVR(N),FINVE(N),SCI(N),SMTI(N),
     &SCV(N),SMIT(N),SMTV(N),SMRI(N),SMRV(N),SEI(N),SEV(N),TR(N),TL(N),     
     &XMUR(N),XMUL(N),XK(N),DFP(N),DFM(N),DRP(N),DRM(N),DEP(N),DEM(N)
      DoublePrecision :: Q(N),XM(N),URHAT(N),C(N),CO(N),DT(N)
      DoublePrecision :: DTETA(N),TETAC,UI,TI,URHAT1,UTHAT1,URHAT2,EADD	
     &,UTHAT2,GAMA,UHAT1,UHAT2,RE,R,XMUZ,CP,PR,RHOI,PIN,EI,CI,PRA,CV
	Integer :: JJ,JJJ,I,J,N1,N2,N3,K,N11,N22,N33,METHOD,BCTYPE,MM,ALFA 
     &, HANEL,IPLOT,SECOND
	DoublePrecision :: TFHAT1,TFHAT2,URFHAT1,URFHAT2,UTFHAT1,UTFHAT2,	 
     &RESET(N),RESMRT(N),RESMTT(N),SUM1(N),SUM2(N),SUM3(N),SUM4(N),QF(N)
     &,DTDT1(N),DTDT2(N),XMIN

      CALL SYSTEM ("del *.obj ")
      CALL SYSTEM ("del *.plg ")
	CALL SYSTEM ("del *.opt ")
	CALL SYSTEM ("del *.dsp ")
	CALL SYSTEM ("del *.dsw ")
	CALL SYSTEM ("del *.txt ")

	! Body of the Code
	
		CALL INPUT				
		CALL GEOM
		CALL INITIALIZE		
10		CALL BC
	    CALL FLUX
          CALL VISCUS
		CALL BC
          CALL UPDATE		
		CALL CALC
	    CALL BC
	    CALL CONVERGENCE 	
		If (RESIDUAL.GT.EPS) GOTO 10							   
		CALL OUTPUT 	

	CONTAINS
!****************************************************************************
!      INPUT DATA
!
!****************************************************************************
	Subroutine Input

	    Open(5,FILE='Output.plt')          		
	    Open(6,FILE='RHO.DAT')
	    Open(1,FILE='U.DAT')	
		Open(2,FILE='RES.DAT')		   
   
	    METHOD=3   ! 1:ROE, 2: VAN LEER 3: H-L-R

	    BCTYPE=3   ! 1:UNIFORM, 2:INVISCID  3:VISCOUS

	    HANEL=1

	    SECOND=0
		
		ALFA=0     ! SWITCH FOR VISCOUSITY

		IF (BCTYPE.EQ.3) ALFA=1  !SWITCH ON

      	GAMA=14.D-1

	    R=287.D0

	    CP=1004.5D0

	    CV=717.5D0

	    PRA=72.D-2      ! PRANDTL

	    MM=1

	    TIME=0

		RE=42.D4

	    TZERO=775.56D0

		XMIN=7.95D0 !4 !11

	    TI=TZERO/(1+0.2D0*XMIN**2)      !56.86D0

	    TCTI=(1+0.2D0*SQRT(PRA)*XMIN**2) !WALL TEMP

	    CI=SQRT(GAMA*R*TI)
	   
	    UI=XMIN*CI !1202.			       	      

	    TETAC=10.*3.14159265359d0/180.

	    DT=7.D-5

C	    DT=4.D-4 ! ROE 

	    EPS=1.D-4	    

	    XM(N)=-206.14/CI      !MACHT	    

	    RHOI=1.D-1

		PIN=RHOI*R*TI

	    XMUZ=286.D-6


	End subroutine 
!****************************************************************************
!    GEOMETRY
!
!****************************************************************************
	Subroutine GEOM		
		   
	    DO I=1,N
		
		DTETA(I)=((2.5*4.8/N))*(3.14159265359D0)/180. 	!2.5 
		TETA(I)=TETAC+(I-1.5)*DTETA(I)
		TETAF(I)=TETAC+(I-2)*DTETA(I)
	    END DO

	End subroutine 
!****************************************************************************
!    INITIALIZE
!
!****************************************************************************
	Subroutine Initialize
	
	    Open(10,FILE='UT.DAT')	
	    Open(20,FILE='UR.DAT')
		OPEN (30,FILE='RHOO.DAT') 
	    OPEN (40,FILE='P.DAT') 
	    OPEN (42,FILE='T.DAT') 

		DO I=1,N

	   GOTO 30
	   READ (10,*) UT(I)
	   READ (20,*) UR(I)
	   READ (30,*) RHO(I)
	   READ (40,*) P(I)
	   READ (42,*) T(I)

C	    GOTO 40

30		UR(I)=COS(TETA(I))*.10  
	    UT(I)=-SIN(TETA(I))*.10 
	    RHO(I)=1.  
	    P(I)=PIN/(RHOI*UI**2)
	    T(I)=1.   
40	    E(I)=CV*TI/UI**2+(UR(I)**2+UT(I)**2)/2.
	    H(I)=E(I)+P(I)/RHO(I)
	    C(I)=SQRT(GAMA*R*T(I)*TI)/UI
	    XM(I)=UT(I)/C(I)

	    END DO
		
		URO=UR	    
	    UTO=UT	   
	    PO=P	    
	    HO=H	    
	    RHOO=RHO  
	    CO=C	    
	    EO=E
		TOLD=T	    
				
	end subroutine 
!****************************************************************************
!    BOUNDARY CONDITIONS
!
!****************************************************************************
	Subroutine BC
	                 
	   SELECT CASE (BCTYPE) ! SELECT BC AT THE FIRST CELL

	   CASE (1)  ! UNIFORM FLOW  

         UR(1)=COS(TETA(1))
	   UT(1)=-SIN(TETA(1))
	   RHO(1)=1.
	   T(1)=1.
         E(1)=CV*TI/UI**2+(UR(1)**2+UT(1)**2)/2.
	   P(1)=PIN/(RHOI*UI**2) 	   
	   H(1)=E(1)+P(1)/RHO(1)
         C(1)=SQRT(GAMA*R*T(1)*TI)/UI	   

         CASE (2)  ! INVISCID FLOW, WALL 
	   UR(1)=UR(2)  
	   UT(1)=-UT(2)   
	   RHO(1)=RHO(2)  
	   E(1)=E(2)
	   H(1)=H(2)
	   C(1)=C(2)
	   P(1)=P(2)
         
	   CASE (3) ! VISCOUS FLOW, WALL

	   UR(1)=-UR(2)  
	   UT(1)=-UT(2) 
	   RHO(1)=RHO(2)  
	   E(1)=E(2)
	   H(1)=H(2)
	   C(1)=C(2)
	   P(1)=P(2)
	   T(1)=T(2)

	   END SELECT

         URO(1:2)=UR(1:2)	    
	   UTO(1:2)=UT(1:2)	  
	   RHOO(1:2)=RHO(1:2)
	   PO(1:2)=P(1:2)  
	   EO(1:2)=E(1:2)
	   HO(1:2)=H(1:2)
	   CO(1:2)=C(1:2) 
	   TOLD(1:2)=T(1:2)
	 
         ! FREE STREAM

         UR(N)=COS(TETA(N))
	   UT(N)=-SIN(TETA(N))
	   RHO(N)=1.
	   T(N)=1.
	   E(N)=CV*TI/UI**2+(UR(N)**2+UT(N)**2)/2. 
	   P(N)=PIN/(RHOI*UI**2)
	   H(N)=E(N)+P(N)/RHO(N)
         C(N)=SQRT(GAMA*R*TI)/UI	   

	   URO(N)=UR(N)	    
	   UTO(N)=UT(N)	  
	   RHOO(N)=RHO(N)
	   PO(N)=P(N)    
	   EO(N)=E(N)
	   EO(N)=E(N)
	   HO(N)=H(N)
	   CO(N)=C(N) 	     
	
	   DO I=1,N	

	   XNU=(ABS(UT(I))+C(I))/DTETA(I)
	   SIGMA=4./3.*XMUZ/(RHO(I)*DTETA(I))**2
C	   DT(I)=1D-3*MIN(1./XNU,1./SIGMA,(1./XNU+1./SIGMA))
C	   DT(I)=1.D-1/(DTETA(I)/(ABS(UTO(I)+CO(I)))
C     &        +4./3.*XMUZ/(RHOI*DTETA(I)**2)) 	
	   END DO  
   
	End subroutine 
!****************************************************************************
!     CALCULATE FLUX
!
!****************************************************************************
	Subroutine FLUX
	
	EADD=5.D-2*(1.D-1*DTETA(1))**2	

      FT(1)=RHO(1)*UT(1)

	FR(1)=RHO(1)*UR(1)

	H(1)=E(1)+P(1)/RHO(1)

	FE(1)=RHO(1)*E(1)

	FT(N)=RHO(N)*UT(N)

	FR(N)=RHO(N)*UR(N)

	H(N)=E(N)+P(N)/RHO(N)

	FE(N)=RHO(N)*E(N)
		
	DO I=1,N

	FTO(I)=RHOO(I)*UTO(I)

	FRO(I)=RHOO(I)*URO(I)

	HO(I)=EO(I)+PO(I)/RHOO(I)

	FEO(I)=RHOO(I)*EO(I)

	END DO

	!!! CONTINUITY

	DO I=2,N-1
	
	! ROE'S AVERAGING

      RHOFH1=SQRT(RHOO(I)*RHOO(I+1))

	RHOFH2=SQRT(RHOO(I)*RHOO(I-1))

	RHOS1=(SQRT(RHOO(I))+SQRT(RHOO(I+1)))

	RHOS2=(SQRT(RHOO(I))+SQRT(RHOO(I-1)))

	PFHAT1=(PO(I)*SQRT(RHOO(I))+PO(I+1)*SQRT(RHOO(I+1)))/RHOS1

	PFHAT2=(PO(I)*SQRT(RHOO(I))+PO(I-1)*SQRT(RHOO(I-1)))/RHOS2

	CHAT1=(CO(I)*SQRT(RHOO(I))+CO(I+1)*SQRT(RHOO(I+1)))/RHOS1  

	CHAT2=(CO(I)*SQRT(RHOO(I))+CO(I-1)*SQRT(RHOO(I-1)))/RHOS2

	URFHAT1=(URO(I)*SQRT(RHOO(I))+URO(I+1)*SQRT(RHOO(I+1)))/RHOS1

	URFHAT2=(URO(I)*SQRT(RHOO(I))+URO(I-1)*SQRT(RHOO(I-1)))/RHOS2

      UTFHAT1=(UTO(I)*SQRT(RHOO(I))+UTO(I+1)*SQRT(RHOO(I+1)))/RHOS1

	UTFHAT2=(UTO(I)*SQRT(RHOO(I))+UTO(I-1)*SQRT(RHOO(I-1)))/RHOS2

      HFHAT1=(HO(I)*SQRT(RHOO(I))+HO(I+1)*SQRT(RHOO(I+1)))/RHOS1

      HFHAT2=(HO(I)*SQRT(RHOO(I))+HO(I-1)*SQRT(RHOO(I-1)))/RHOS2

	TFHAT1=(TOLD(I)*SQRT(RHOO(I))+TOLD(I+1)*SQRT(RHOO(I+1)))/RHOS1

	TFHAT2=(TOLD(I)*SQRT(RHOO(I))+TOLD(I-1)*SQRT(RHOO(I-1)))/RHOS2

      ! I+0.5 

	DP1=(PO(I+1)-PO(I))

	DU1=(UTO(I+1)-UTO(I))

	DRHO1=(RHOO(I+1)-RHOO(I))

	DV1=(URO(I+1)-URO(I))

	! SECOND ORDER FLUXES

	IF (I.GT.3.AND.I.LT.N-2) THEN

	IF (SECOND.EQ.1) THEN

	RHOR1=(3./2.*RHOO(I+1)-1./2.*RHOO(I+2))

	RHOL1=(3./2.*RHOO(I)-1./2.*RHOO(I-1))

	PR1=(3./2.*PO(I+1)-1./2.*PO(I+2))

	PL1=(3./2.*PO(I)-1./2.*PO(I-1))

	CR1=(3./2.*CO(I+1)-1./2.*CO(I+2))

	CL1=(3./2.*CO(I)-1./2.*CO(I-1))

	UTR1=(3./2.*UTO(I+1)-1./2.*UTO(I+2))

	UTL1=(3./2.*UTO(I)-1./2.*UTO(I-1))

	VR1=(3./2.*URO(I+1)-1./2.*URO(I+2))

	VL1=(3./2.*URO(I)-1./2.*URO(I-1))

	HR1=(3./2.*HO(I+1)-1./2.*HO(I+2))

	HL1=(3./2.*HO(I)-1./2.*HO(I-1))

	TR1=(3./2.*TOLD(I+1)-1./2.*TOLD(I+2))

	TL1=(3./2.*TOLD(I)-1./2.*TOLD(I-1))

	RHOR2=(3./2.*RHOO(I)-1./2.*RHOO(I+1))

      RHOL2=(3./2.*RHOO(I-1)-1./2.*RHOO(I-2))

	PR2=(3./2.*PO(I)-1./2.*PO(I+1))
	
	PL2=(3./2.*PO(I-1)-1./2.*PO(I-2))

	CR2=(3./2.*CO(I)-1./2.*CO(I+1))

	CL2=(3./2.*CO(I-1)-1./2.*CO(I-2))

	UTR2=(3./2.*UTO(I)-1./2.*UTO(I+1))

	UTL2=(3./2.*UTO(I-1)-1./2.*UTO(I-2))
	
	VR2=(3./2.*URO(I)-1./2.*URO(I+1))

	VL2=(3./2.*URO(I-1)-1./2.*URO(I-2))	
	
	HR2=(3./2.*HO(I)-1./2.*HO(I+1))

	HL2=(3./2.*HO(I-1)-1./2.*HO(I-2))

	TR2=(3./2.*TOLD(I)-1./2.*TOLD(I+1))

	TL2=(3./2.*TOLD(I-1)-1./2.*TOLD(I-2))

	RHOFH1=SQRT(RHOR1*RHOL1)

	RHOFH2=SQRT(RHOR2*RHOL2)

	RHOS1=(SQRT(RHOR1)+SQRT(RHOL1))

	RHOS2=(SQRT(RHOR2)+SQRT(RHOL2))

	PFHAT1=(PL1*SQRT(RHOL1)+PR1*SQRT(RHOR1))/RHOS1

	PFHAT2=(PR2*SQRT(RHOR2)+PL2*SQRT(RHOL2))/RHOS2

	CHAT1=(CL1*SQRT(RHOL1)+CR1*SQRT(RHOR1))/RHOS1  

	CHAT2=(CR2*SQRT(RHOR2)+CL2*SQRT(RHOL2))/RHOS2

	URFHAT1=(VL1*SQRT(RHOL1)+VR1*SQRT(RHOR1))/RHOS1

	URFHAT2=(VR2*SQRT(RHOR2)+VL2*SQRT(RHOL2))/RHOS2

      UTFHAT1=(UTL1*SQRT(RHOL1)+UTR1*SQRT(RHOR1))/RHOS1

	UTFHAT2=(UTR2*SQRT(RHOR2)+UTL2*SQRT(RHOL2))/RHOS2

      HFHAT1=(HL1*SQRT(RHOL1)+HR1*SQRT(RHOR1))/RHOS1

      HFHAT2=(HR2*SQRT(RHOR2)+HL2*SQRT(RHOL2))/RHOS2

	TFHAT1=(TL1*SQRT(RHOL1)+TR1*SQRT(RHOR1))/RHOS1

	TFHAT2=(TR2*SQRT(RHOR2)+TL2*SQRT(RHOL2))/RHOS2

      DP1=PR1-PL1

	DU1=UTR1-UTL1 

	DRHO1=RHOR1-RHOL1     

	DV1=VR1-VL1

	END IF
      END IF 
	
	V1=(DP1-RHOFH1*CHAT1*DU1)/(2*CHAT1**2)

	COF1=ABS(UTFHAT1-CHAT1)*V1

	V2=-(DP1-(CHAT1**2)*DRHO1)/(CHAT1**2)

	COF2=ABS(UTFHAT1)*V2

	V3=RHOFH1*(DV1)/(URFHAT1+EADD)

	COF3=ABS(UTFHAT1)*V3

	V4=(DP1+RHOFH1*CHAT1*DU1)/(2*CHAT1**2)

	COF4=ABS(UTFHAT1+CHAT1)*V4

	! FIRST ORDER

	DP2=PO(I)-PO(I-1)
	
	DU2=(UTO(I)-UTO(I-1))

	DV2=(URO(I)-URO(I-1))

	DRHO2=(RHOO(I)-RHOO(I-1))

	IF (I.GT.3.AND.I.LT.N-2) THEN
	IF (SECOND.EQ.1) THEN

	DP2=PR2-PL2 

	DU2=UTR2-UTL2

	DRHO2=RHOR2-RHOL2 

	DV2=VR2-VL2 

	END IF
	END IF

	V5=(DP2-RHOFH2*CHAT2*DU2)/(2*CHAT2**2)

	COF5=ABS(UTFHAT2-CHAT2)*V5

      V6=-(DP2-(CHAT2**2)*DRHO2)/(CHAT2**2)  

	COF6=ABS(UTFHAT2)*V6

	V7=RHOFH2*DV2/(URFHAT2+EADD)
	
	COF7=ABS(UTFHAT2)*V7

	V8=(DP2+RHOFH2*CHAT2*DU2)/(2*CHAT2**2)

	COF8=ABS(UTFHAT2+CHAT2)*V8

      F11=COF1

	F12=COF2
     
	F13=0.

	F14=COF4

	! I-0.5

      F21=COF5

	F22=COF6

	F23=0.

	F24=COF8

	! SECOND ORDER FLUXES

	IF (I.GT.3.AND.I.LT.N-2) THEN

	WP1=(3./2.*FTO(I+1)-1./2.*FTO(I+2))  
      WP2=(3./2.*FTO(I)-1./2.*FTO(I-1))

	WM1=(3./2.*FTO(I)-1./2.*FTO(I+1))
	WM2=(3./2.*FTO(I-1)-1./2.*FTO(I-2))

	END IF

	! VAN LEER FLUX (SUBSONIC)

	FMASS1=RHOO(I)*CO(I)*(XM(I)+1)**2/4.        ! F PLUS J

	FMASS2=-RHOO(I+1)*CO(I+1)*(XM(I+1)-1)**2/4. ! F MINUS J+1

	FMASS3=RHOO(I-1)*CO(I-1)*(XM(I-1)+1)**2/4.  ! F PLUS J-1

	FMASS4=-RHOO(I)*CO(I)*(XM(I)-1)**2/4.       ! F MINUS J

       ! SECOND ORDER VAN LEER FLUX (SUBSONIC)

	IF (I.GT.3.AND.I.LT.N-2) THEN
	IF (SECOND.EQ.1) THEN

      FMASS1=RHOL1*CL1*((UTL1/CL1)+1)**2/4.        ! F PLUS J

	FMASS2=-RHOR1*CR1*((UTR1/CR1)-1)**2/4. ! F MINUS J+1

	FMASS3=RHOL2*CL2*((UTL2/CL2)+1)**2/4.  ! F PLUS J-1

	FMASS4=-RHOR2*CR2*((UTR2/CR2)-1)**2/4.       ! F MINUS J
	
	END IF
	END IF
	

	SELECT CASE (METHOD)

	CASE (1)  ! ROE

	DFP(I)=1./2.*(FTO(I)+FTO(I+1))-1./2.*(F11+F12+F13+F14)

	DFM(I)=1./2.*(FTO(I)+FTO(I-1))-1./2.*(F21+F22+F23+F24) 
!------------------------------------------------------------------------    
	IF (I.GT.3.AND.I.LT.N-2) THEN
	IF (SECOND.EQ.1) THEN	

      DFP(I)=1./2.*(WP1+WP2)-1./2.*(F11+F12+F13+F14)

	DFM(I)=1./2.*(WM1+WM2)-1./2.*(F21+F22+F23+F24) 

	END IF
	END IF
!------------------------------------------------------------------------    
	! CALCULATE FLUX RESIDUAL

	DF(I)=DFP(I)-DFM(I)
      
	CONTINUE
!------------------------------------------------------------------------    
	CASE(3) !H-L-R
      
	XLS1=V1**2+V2**2+V3**2+V4**2+EADD
	
	XLS2=V5**2+V6**2+V7**2+V8**2+EADD

	XUS1=V1**2*ABS(UTFHAT1-CHAT1)+V2**2*ABS(UTFHAT1)
     &	+V3**2*ABS(UTFHAT1)+V4**2*ABS(UTFHAT1+CHAT1)

	XUS2=V5**2*ABS(UTFHAT2-CHAT2)+V6**2*ABS(UTFHAT2)
     &	+V7**2*ABS(UTFHAT2)+V8**2*ABS(UTFHAT2+CHAT2)

      XHLR1=ABS(XUS1/XLS1)

	XHLR2=ABS(XUS2/XLS2)

      DFP(I)=1./2.*(FTO(I)+FTO(I+1))-1./2.*XHLR1*(RHOO(I+1)-RHOO(I))

	DFM(I)=1./2.*(FTO(I)+FTO(I-1))-1./2.*XHLR2*(RHOO(I)-RHOO(I-1))

	DF(I)=DFP(I)-DFM(I)
!------------------------------------------------------------------------    
	CASE(2)  ! VAN LEER	

	IF (XM(I).GE.1)	      DF(I)=FTO(I)-FTO(I-1)

	IF (XM(I).LE.-1)      DF(I)=FTO(I+1)-FTO(I)
	
	IF (ABS(XM(I)).LT.1) DF(I)=(FMASS1+FMASS2)-(FMASS3+FMASS4)

	IF (SECOND.EQ.1) THEN

      IF (XM(I).GE.1)	      DF(I)=WP2-WM2

	IF (XM(I).LE.-1) 	  DF(I)=WP1-WM1

	IF (ABS(XM(I)).LT.1)  DF(I)=(FMASS1+FMASS2)-(FMASS3+FMASS4)

	END IF		

	END SELECT
!-----------------------------------------------------------------------
	! CALCULATE FLUX RESIDUAL
	
	RESC1(I)=-DF(I)/DTETA(I)
!-----------------------------------------------------------------------
	!INVISCID CONTINUITY SOURCE
	
	X1=(2.*FRO(I-1)+FTO(I-1)/TAN(TETA(I-1)))

      X2=2.*(2*FRO(I)+FTO(I)/TAN(TETA(I))) 
     
      X21=(2.*FRO(I+1)+FTO(I+1)/TAN(TETA(I+1)))
	
C	SCI(I)=-(X1+X2+X21)/4.  	 
	
	! ROE AVERAGING FOR SOURCE

	SRC1=2*RHOFH1*URFHAT1+RHOFH1*UTFHAT1/TAN(TETAF(I+1))

	SRC2=2*RHOFH2*URFHAT2+RHOFH2*UTFHAT2/TAN(TETAF(I))

	SCI(I)=-(SRC1+SRC2)/2.

	SCV(I)=0.

      RESC2(I)=(SCI(I)+SCV(I))

	RESC(I)=RESC1(I)+RESC2(I)

	!! UPDATE CONTINUITY

	RHO(I)=RHOO(I)+DT(I)*RESC(I)
!------------------------------------------------------------------------
	!!!!!!!! TETA MOMENTUN

	FRT1=((GAMA-1)*UTO(I)+2.*CO(I))/GAMA
	
	FRT2=((GAMA-1)*UTO(I+1)-2.*CO(I+1))/GAMA

	FRT3=((GAMA-1)*UTO(I-1)+2.*CO(I-1))/GAMA

	FRT4=((GAMA-1)*UTO(I)-2.*CO(I))/GAMA

	IF (I.GT.3.AND.I.LT.N-2) THEN
	IF (SECOND.EQ.1) THEN

      FRT1=((GAMA-1)*UTL1+2.*CL1)/GAMA
	
	FRT2=((GAMA-1)*UTR1-2.*CR1)/GAMA

	FRT3=((GAMA-1)*UTL2+2.*CL2)/GAMA

	FRT4=((GAMA-1)*UTR2-2.*CR2)/GAMA
	
	END IF
	END IF
		
	F31=COF1*(UTFHAT1-CHAT1)
	F32=COF2*UTFHAT1
	F33=0.
	F34=COF4*(UTFHAT1+CHAT1)

	! I-0.5

      F41=COF5*(UTFHAT2-CHAT2)
	F42=COF6*UTFHAT2
	F43=0.
	F44=COF8*(UTFHAT2+CHAT2)

	! SECOND ORDER FLUXES
	IF (I.GT.3.AND.I.LT.N-2) THEN

	TP1=3./2.*(FTO(I+1)*UTO(I+1)+PO(I+1))
     &   -1./2.*(FTO(I+2)*UTO(I+2)+PO(I+2))
	
	TP2=3./2.*(FTO(I)*UTO(I)+PO(I))
     &   -1./2.*(FTO(I-1)*UTO(I-1)+PO(I-1))
      
	TM1=3./2.*(FTO(I)*UTO(I)+PO(I))
     &   -1./2.*(FTO(I+1)*UTO(I+1)+PO(I+1))
	
	TM2=3./2.*(FTO(I-1)*UTO(I-1)+PO(I-1))
     &   -1./2.*(FTO(I-2)*UTO(I-2)+PO(I-2))
      
	END IF  

	SELECT CASE (METHOD)

	CASE (1)  ! ROE'S AVERAGING

	! FIRST ORDER
    	
	SX1=1./2.*(FTO(I)*UTO(I)+PO(I)+FTO(I+1)*UTO(I+1)+PO(I+1))
     &-1./2.*(F31+F32+F33+F34)

	SX2=1./2.*(FTO(I)*UTO(I)+PO(I)+FTO(I-1)*UTO(I-1)+PO(I-1))
     &-1./2.*(F41+F42+F43+F44) 
	
	FINVT(I)=SX1-SX2
!------------------------------------------------------------------------    
	IF (I.GT.3.AND.I.LT.N-2) THEN
	IF (SECOND.EQ.1) THEN

	SY1=1./2.*(TP1+TP2)-1./2.*(F31+F32+F33+F34)

	SY2=1./2.*(TM1+TM2)-1./2.*(F31+F32+F33+F34)
	
	FINVT(I)=SY1-SY2 

	END IF
      END IF
	CONTINUE
!------------------------------------------------------------------------    
	CASE(3)

	SXHLR1=1./2.*(FTO(I)*UTO(I)+PO(I)+FTO(I+1)*UTO(I+1)+PO(I+1))
     &-1./2.*XHLR1*(FTO(I+1)-FTO(I))

	SXHLR2=1./2.*(FTO(I)*UTO(I)+PO(I)+FTO(I-1)*UTO(I-1)+PO(I-1))
     &-1./2.*XHLR2*(FTO(I)-FTO(I-1)) 
	
	FINVT(I)=SXHLR1-SXHLR2		 
!------------------------------------------------------------------------
	CASE(2)  ! VAN LEER	
	
	IF (XM(I).GE.1) THEN
	
	FINVT(I)=(FTO(I)*UTO(I)+PO(I)-FTO(I-1)*UTO(I-1)-PO(I-1))

	ELSE IF (XM(I).LE.-1) THEN

	FINVT(I)=(FTO(I+1)*UTO(I+1)+PO(I+1)-FTO(I)*UTO(I)-PO(I))

	ELSE IF (ABS(XM(I)).LT.1) THEN
	
	FINVT(I)=(FMASS1*FRT1+FMASS2*FRT2)-(FMASS3*FRT3+FMASS4*FRT4)

      END IF

      IF (SECOND.EQ.1) THEN
      
	IF (XM(I).GE.1)         FINVT(I)=TP2-TM2

	IF (XM(I).LE.-1) 	    FINVT(I)=TP1-TM1
  
	END IF
      
      END SELECT
!------------------------------------------------------------------------
	! CALCULATE FLUX RESIDUAL

      RESMT1(I)=-FINVT(I)/DTETA(I)
!------------------------------------------------------------------------
      !INVISCID T-MOM SOURCE
	
	X3=(3*FRO(I)*UTO(I)+FTO(I)*UTO(I)/TAN(TETA(I)))+
	
     &(3*FRO(I+1)*UTO(I+1)+FTO(I+1)*UTO(I+1)/TAN(TETA(I+1)))	 	 

	X4=(3*FRO(I)*UTO(I)+FTO(I)*UTO(I)/TAN(TETA(I)))+
	
     &(3*FRO(I-1)*UTO(I-1)+FTO(I-1)*UTO(I-1)/TAN(TETA(I-1)))
     
C    	IF (I.EQ.2) X4=0.
	
C	SMTI(I)=-(X3+X4)/4.	 

      ! ROE AVERAGING FOR SOURCE TERM 
	
	SRMT1=3*RHOFH1*URFHAT1*UTFHAT1+RHOFH1*UTFHAT1**2/TAN(TETAF(I+1))
	
	SRMT2=3*RHOFH2*URFHAT2*UTFHAT2+RHOFH2*UTFHAT2**2/TAN(TETAF(I))	 

	SMTI(I)=-(SRMT1+SRMT2)/2.	
	
	RESMT2(I)=SMTI(I)

	! SOURCE RESIDUAL

	RESMT(I)=RESMT1(I)+RESMT2(I)

	!! UPDATE T-MOM

C	FT(I)=FTO(I)+DT(I)*RESMT(I)
C	UT(I)=FT(I)/RHO(I)
!------------------------------------------------------------------------
	!!!!!!!! R MOMENTUN
	FRR1=URO(I)
	
	FRR2=URO(I+1)

	FRR3=URO(I-1)

	FRR4=URO(I)	

	IF (I.GT.3.AND.I.LT.N-2) THEN
	IF (SECOND.EQ.1) THEN

	FRR1=VL1
	
	FRR2=VR1

	FRR3=VL2

	FRR4=VR2
      	
	END IF
	END IF

	F51=COF1*URFHAT1

	F52=COF2*URFHAT1

	F53=COF3*URFHAT1

	F54=COF4*URFHAT1

	! I-0.5

      F61=COF5*URFHAT2

	F62=COF6*URFHAT2

	F63=COF7*URFHAT2

	F64=COF8*URFHAT2

	! SECOND ORDER FLUXES
	
	IF (I.GT.3.AND.I.LT.N-2) THEN

	RP1=3./2.*(FRO(I+1)*UTO(I+1))-1./2.*(FRO(I+2)*UTO(I+2))
	
	RP2=3./2.*(FRO(I)*UTO(I))-1./2.*(FRO(I-1)*UTO(I-1))
      
	RM1=3./2.*(FRO(I)*UTO(I))-1./2.*(FRO(I+1)*UTO(I+1))
	
	RM2=3./2.*(FRO(I-1)*UTO(I-1))-1./2.*(FRO(I-2)*UTO(I-2))

	END IF

	SELECT CASE (METHOD)

	CASE (1)  ! ROE'S AVERAGING

	DRP(I)=1./2.*(FRO(I)*UTO(I)+FRO(I+1)*UTO(I+1))
     &	  -1./2.*(F51+F52+F53+F54)

	DRM(I)=1./2.*(FRO(I)*UTO(I)+FRO(I-1)*UTO(I-1))
     &	  -1./2.*(F61+F62+F63+F64)
     		
	FINVR(I)=DRP(I)-DRM(I)  
!------------------------------------------------------------------------    
	IF (I.GT.3.AND.I.LT.N-2) THEN
	IF (SECOND.EQ.1) THEN

	SR1=1./2.*(RP1+RP2)-1./2.*(F51+F52+F53+F54)

	SR2=1./2.*(RM1+RM2)-1./2.*(F61+F62+F63+F64)
	
	FINVT(I)=SR1-SR2 

	END IF
	END IF
    	CONTINUE
!------------------------------------------------------------------------    
	CASE(3)

	DRP(I)=1./2.*(FRO(I)*UTO(I)+FRO(I+1)*UTO(I+1))
     &      -1./2.*XHLR1*(FRO(I+1)-FRO(I))

	DRM(I)=1./2.*(FRO(I)*UTO(I)+FRO(I-1)*UTO(I-1))
     &      -1./2.*XHLR2*(FRO(I)-FRO(I-1))
     		
	FINVR(I)=DRP(I)-DRM(I)      
!------------------------------------------------------------------------
	CASE(2)  ! VAN LEER	

	IF (XM(I).GE.1) THEN
	
	FINVR(I)=(FRO(I)*UTO(I)-FRO(I-1)*UTO(I-1))

	ELSE IF (XM(I).LE.-1) THEN
	
	FINVR(I)=(FRO(I+1)*UTO(I+1)-FRO(I)*UTO(I))
    
      ELSE IF (ABS(XM(I)).LT.1) THEN

	FINVR(I)=(FMASS1*FRR1+FMASS2*FRR2)-(FMASS3*FRR3+FMASS4*FRR4)

	END IF

	IF (SECOND.EQ.1) THEN
      
	IF (XM(I).GE.1)         FINVR(I)=RP2-RM2

	IF (XM(I).LE.-1) 	    FINVR(I)=RP1-RM1
  
	END IF

	END SELECT
!------------------------------------------------------------------------
	! CALCULATE FLUX RESIDUAL

	RESMR1(I)=-FINVR(I)/DTETA(I)
!------------------------------------------------------------------------
      !INVISCID R-MOM SOURCE
   
	X5=( 2*FRO(I)*URO(I)-FTO(I)*UTO(I)+FRO(I)*UTO(I)/TAN(TETA(I)) )
	
     &+  ( 2*FRO(I+1)*URO(I+1)-FTO(I+1)*UTO(I+1)+FRO(I+1)*UTO(I+1)
     
     &/TAN(TETA(I+1)))	 	 

	X6=(2*FRO(I)*URO(I)-FTO(I)*UTO(I)+FRO(I)*UTO(I)/TAN(TETA(I)))
	
     &  +( 2*FRO(I-1)*URO(I-1)-FTO(I-1)*UTO(I-1)+FRO(I-1)*UTO(I-1)
     
     &/TAN(TETA(I-1)))
     
C      SMRI(I)=-(X5+X6)/4.	 
	
	 ! ROE AVERAGING FOR SOURCE TERM 

      SRMR1=2*RHOFH1*URFHAT1**2-RHOFH1*UTFHAT1**2
	
     &	+RHOFH1*URFHAT1*UTFHAT1/TAN(TETAF(I+1))

	SRMR2=2*RHOFH2*URFHAT2**2-RHOFH2*UTFHAT2**2
	
     &	+RHOFH2*URFHAT2*UTFHAT2/TAN(TETAF(I))

	SMRI(I)=-(SRMR1+SRMR2)/2.

      RESMR2(I)=SMRI(I)

	RESMR(I)=RESMR1(I)+RESMR2(I)
	
	!! UPDATE R-MOM

C	FR(I)=FRO(I)+DT(I)*RESMR(I)
C	UR(I)=FR(I)/RHO(I)
!------------------------------------------------------------------------
	!!!!!!!! ENERGY EQUATION

	DOM=2.*(GAMA**2-1.)

	FRE1=(FRT1)**2/DOM+(URO(I)**2)/2.
	
	FRE2=(FRT2)**2/DOM+(URO(I+1)**2)/2.

	FRE3=(FRT3)**2/DOM+(URO(I-1)**2)/2.

	FRE4=(FRT4)**2/DOM+(URO(I)**2)/2.
	
	IF (I.GT.3.AND.I.LT.N-2) THEN
	IF (SECOND.EQ.1) THEN

	FRE1=(FRT1)**2/DOM+(VL1**2)/2.
	
	FRE2=(FRT2)**2/DOM+(VR1**2)/2.

	FRE3=(FRT3)**2/DOM+(VL2**2)/2.

	FRE4=(FRT4)**2/DOM+(VR2**2)/2.
		   	
	END IF
	END IF			

	! I+0.5 

      F71=COF1*(HFHAT1-UTFHAT1*CHAT1)

	F72=COF2*1./2.*(UTFHAT1**2+URFHAT1**2)

	F73=COF3*URFHAT1**2

	F74=COF4*(HFHAT1+UTFHAT1*CHAT1)

	! I-0.5

      F81=COF5*(HFHAT2-UTFHAT2*CHAT2)

	F82=COF6*1./2.*(URFHAT2**2+UTFHAT2**2)

	F83=COF7*URFHAT2**2

	F84=COF8*(HFHAT2+UTFHAT2*CHAT2)

	! SECOND ORDER FLUXES
	
	IF (I.GT.3.AND.I.LT.N-2) THEN

	EP1=3./2.*(FTO(I+1)*HO(I+1))-1./2.*(FTO(I+2)*HO(I+2))
	
	EP2=3./2.*(FTO(I)*HO(I))-1./2.*(FTO(I-1)*HO(I-1))
      
	EM1=3./2.*(FTO(I)*HO(I))-1./2.*(FTO(I+1)*HO(I+1))
	
	EM2=3./2.*(FTO(I-1)*HO(I-1))-1./2.*(FTO(I-2)*HO(I-2))
	
	END IF
	
	SELECT CASE (METHOD)	
	
	CASE (1)  ! ROE'S AVERAGING

	DEP(I)=1./2.*(FTO(I)*HO(I)+FTO(I+1)*HO(I+1))-1./2.*(F71+F72+F73+F74)

	DEM(I)=1./2.*(FTO(I)*HO(I)+FTO(I-1)*HO(I-1))-1./2.*(F81+F82+F83+F84)
     		
	FINVE(I)=DEP(I)-DEM(I)
!------------------------------------------------------------------------    
	IF (I.GT.3.AND.I.LT.N-2) THEN
	IF (SECOND.EQ.1) THEN

	SE1=1./2.*(EP1+EP2)-1./2.*(F51+F52+F53+F54)

	SE2=1./2.*(EM1+EM2)-1./2.*(F61+F62+F63+F64)
	
	FINVT(I)=SE1-SE2 

	END IF
      END IF

      CONTINUE
!------------------------------------------------------------------------    
	CASE(3)
	
	DEP(I)=1./2.*(FTO(I)*HO(I)+FTO(I+1)*HO(I+1))
     &	  -1./2.*XHLR1*(FEO(I+1)-FEO(I))

	DEM(I)=1./2.*(FTO(I)*HO(I)+FTO(I-1)*HO(I-1))
     &      -1./2.*XHLR2*(FEO(I)-FEO(I-1))
     		
	FINVE(I)=DEP(I)-DEM(I)   
!------------------------------------------------------------------------	
	CASE(2)  ! VAN LEER	
	
	IF (XM(I).GE.1) THEN
	
	FINVE(I)=(FTO(I)*HO(I)-FTO(I-1)*HO(I-1))

	ELSE IF (XM(I).LE.-1) THEN

      FINVE(I)=(FTO(I+1)*HO(I+1)-FTO(I)*HO(I))
	
	ELSE IF (ABS(XM(I)).LT.1) THEN
	
	FINVE(I)=(FMASS1*FRE1+FMASS2*FRE2)-(FMASS3*FRE3+FMASS4*FRE4)

	IF (HANEL.EQ.1) FINVE(I)=(FMASS1*HO(I)+FMASS2*HO(I+1))
     &	                    -(FMASS3*HO(I-1)+FMASS4*HO(I))

	END IF

	IF (SECOND.EQ.1) THEN
      
	IF (XM(I).GE.1)         FINVE(I)=EP2-EM2

	IF (XM(I).LE.-1) 	    FINVE(I)=EP1-EM1
  
	END IF

      END SELECT

	!! UPDATE ENERGY 

	RESE1(I)=-FINVE(I)/DTETA(I)      !!DVEL-FINVE(I)/DTETA(I)
!------------------------------------------------------------------------- 
	!! INVISCID ENERGY SOURCE

	X7=(2.*FRO(I)*HO(I)+FTO(I)*HO(I)/TAN(TETA(I)))
	
     &+(2.*FRO(I+1)*HO(I+1)+FTO(I+1)*HO(I+1)/TAN(TETA(I+1)))	 	 

	X8=(2.*FRO(I)*HO(I)+FTO(I)*HO(I)/TAN(TETA(I)))+
	
     &(2.*FRO(I-1)*HO(I-1)+FTO(I-1)*HO(I-1)/TAN(TETA(I-1))) 
     
C      SEI(I)=-(X7+X8)/4.

	! ROE AVERAGING FOR SOURCE TERM 

	SRE1=2*RHOFH1*URFHAT1*HFHAT1+RHOFH1*UTFHAT1*HFHAT1/TAN(TETAF(I+1))

	SRE2=2*RHOFH2*URFHAT2*HFHAT2+RHOFH2*UTFHAT2*HFHAT2/TAN(TETAF(I))

	SEI(I)=-(SRE1+SRE2)/2.

      RESE2(I)=SEI(I)

	RESE(I)=RESE1(I)+RESE2(I)

	! UPDATE ENERGY

C      FE(I)=FEO(I)+DT(I)*RESE(I)

C	E(I)=FE(I)/RHO(I)  

	END DO

	End subroutine
!****************************************************************************
!     CALCULATE VICSOUS TERMS
!
!****************************************************************************
      Subroutine VISCUS

	DO I=2,N-1

	RHOFH1=SQRT(RHOO(I)*RHOO(I+1))

	RHOFH2=SQRT(RHOO(I)*RHOO(I-1))
	
	RHOS1=(SQRT(RHOO(I))+SQRT(RHOO(I+1)))

	RHOS2=(SQRT(RHOO(I))+SQRT(RHOO(I-1)))

      PFHAT1=(PO(I)*SQRT(RHOO(I))+PO(I+1)*SQRT(RHOO(I+1)))/RHOS1

	PFHAT2=(PO(I)*SQRT(RHOO(I))+PO(I-1)*SQRT(RHOO(I-1)))/RHOS2	

	CHAT1=(CO(I)*SQRT(RHOO(I))+CO(I+1)*SQRT(RHOO(I+1)))/RHOS1  

	CHAT2=(CO(I)*SQRT(RHOO(I))+CO(I-1)*SQRT(RHOO(I-1)))/RHOS2

      URFHAT1=(URO(I)*SQRT(RHOO(I))+URO(I+1)*SQRT(RHOO(I+1)))/RHOS1

	URFHAT2=(URO(I)*SQRT(RHOO(I))+URO(I-1)*SQRT(RHOO(I-1)))/RHOS2

      UTFHAT1=(UTO(I)*SQRT(RHOO(I))+UTO(I+1)*SQRT(RHOO(I+1)))/RHOS1

	UTFHAT2=(UTO(I)*SQRT(RHOO(I))+UTO(I-1)*SQRT(RHOO(I-1)))/RHOS2

      HFHAT1=(HO(I)*SQRT(RHOO(I))+HO(I+1)*SQRT(RHOO(I+1)))/RHOS1

      HFHAT2=(HO(I)*SQRT(RHOO(I))+HO(I-1)*SQRT(RHOO(I-1)))/RHOS2

	TFHAT1=(TOLD(I)*SQRT(RHOO(I))+TOLD(I+1)*SQRT(RHOO(I+1)))/RHOS1

	TFHAT2=(TOLD(I)*SQRT(RHOO(I))+TOLD(I-1)*SQRT(RHOO(I-1)))/RHOS2

	IF (I.GT.3.AND.I.LT.N-2) THEN

	IF (SECOND.EQ.1) THEN

	RHOR1=(3./2.*RHOO(I+1)-1./2.*RHOO(I+2))

	RHOL1=(3./2.*RHOO(I)-1./2.*RHOO(I-1))

	PR1=(3./2.*PO(I+1)-1./2.*PO(I+2))

	PL1=(3./2.*PO(I)-1./2.*PO(I-1))

	CR1=(3./2.*CO(I+1)-1./2.*CO(I+2))

	CL1=(3./2.*CO(I)-1./2.*CO(I-1))

	UTR1=(3./2.*UTO(I+1)-1./2.*UTO(I+2))

	UTL1=(3./2.*UTO(I)-1./2.*UTO(I-1))

	VR1=(3./2.*URO(I+1)-1./2.*URO(I+2))

	VL1=(3./2.*URO(I)-1./2.*URO(I-1))

	HR1=(3./2.*HO(I+1)-1./2.*HO(I+2))

	HL1=(3./2.*HO(I)-1./2.*HO(I-1))

	TR1=(3./2.*TOLD(I+1)-1./2.*TOLD(I+2))

	TL1=(3./2.*TOLD(I)-1./2.*TOLD(I-1))

	RHOR2=(3./2.*RHOO(I)-1./2.*RHOO(I+1))

      RHOL2=(3./2.*RHOO(I-1)-1./2.*RHOO(I-2))

	PR2=(3./2.*PO(I)-1./2.*PO(I+1))
	
	PL2=(3./2.*PO(I-1)-1./2.*PO(I-2))

	CR2=(3./2.*CO(I)-1./2.*CO(I+1))

	CL2=(3./2.*CO(I-1)-1./2.*CO(I-2))

	UTR2=(3./2.*UTO(I)-1./2.*UTO(I+1))

	UTL2=(3./2.*UTO(I-1)-1./2.*UTO(I-2))
	
	VR2=(3./2.*URO(I)-1./2.*URO(I+1))

	VL2=(3./2.*URO(I-1)-1./2.*URO(I-2))	
	
	HR2=(3./2.*HO(I)-1./2.*HO(I+1))

	HL2=(3./2.*HO(I-1)-1./2.*HO(I-2))

	TR2=(3./2.*TOLD(I)-1./2.*TOLD(I+1))

	TL2=(3./2.*TOLD(I-1)-1./2.*TOLD(I-2))

	RHOFH1=SQRT(RHOR1*RHOL1)

	RHOFH2=SQRT(RHOR2*RHOL2)

	RHOS1=(SQRT(RHOR1)+SQRT(RHOL1))

	RHOS2=(SQRT(RHOR2)+SQRT(RHOL2))

	PFHAT1=(PL1*SQRT(RHOL1)+PR1*SQRT(RHOR1))/RHOS1

	PFHAT2=(PR2*SQRT(RHOR2)+PL2*SQRT(RHOL2))/RHOS2

	CHAT1=(CL1*SQRT(RHOL1)+CR1*SQRT(RHOR1))/RHOS1  

	CHAT2=(CR2*SQRT(RHOR2)+CL2*SQRT(RHOL2))/RHOS2

	URFHAT1=(VL1*SQRT(RHOL1)+VR1*SQRT(RHOR1))/RHOS1

	URFHAT2=(VR2*SQRT(RHOR2)+VL2*SQRT(RHOL2))/RHOS2

      UTFHAT1=(UTL1*SQRT(RHOL1)+UTR1*SQRT(RHOR1))/RHOS1

	UTFHAT2=(UTR2*SQRT(RHOR2)+UTL2*SQRT(RHOL2))/RHOS2

      HFHAT1=(HL1*SQRT(RHOL1)+HR1*SQRT(RHOR1))/RHOS1

      HFHAT2=(HR2*SQRT(RHOR2)+HL2*SQRT(RHOL2))/RHOS2

	TFHAT1=(TL1*SQRT(RHOL1)+TR1*SQRT(RHOR1))/RHOS1

	TFHAT2=(TR2*SQRT(RHOR2)+TL2*SQRT(RHOL2))/RHOS2

     	END IF
      END IF 

	TR(I)=TFHAT1  !(TOLD(I)+TOLD(I+1))/2.

	TL(I)=TFHAT2  !(TOLD(I)+TOLD(I-1))/2.

	XMUR(I)=ALFA*1.*(TI+1104.D-1)/(TR(I)*TI+1104.D-1)*(TR(I))**(3./2.)

	XMUL(I)=ALFA*1.*(TI+1104.D-1)/(TL(I)*TI+1104.D-1)*(TL(I))**(3./2.)
	
	DUTDT1(I)=(UTO(I+1)-UTO(I))/DTETA(I)

	DUTDT2(I)=(UTO(I)-UTO(I-1))/DTETA(I)

	DURDT1(I)=(URO(I+1)-URO(I))/DTETA(I)

	DURDT2(I)=(URO(I)-URO(I-1))/DTETA(I)

	UR1(I)=URFHAT1 !(URO(I)+URO(I+1))/2.

	UR2(I)=URFHAT2 !(URO(I-1)+URO(I))/2.

	UT1(I)=UTFHAT1 !(UTO(I)+UTO(I+1))/2.

	UT2(I)=UTFHAT2 !(UTO(I-1)+UTO(I))/2.

	TEF1=TAN(TETAF(I+1))

	TEF2=TAN(TETAF(I))
	
	STT1(I)=2./3.*XMUR(I)/RE*(2.*DUTDT1(I)+UR1(I)-UT1(I)/TEF1)

	STT2(I)=2./3.*XMUL(I)/RE*(2.*DUTDT2(I)+UR2(I)-UT2(I)/TEF2)

	SRR1(I)=-2./3.*XMUR(I)/RE*(DUTDT1(I)+2.*UR1(I)+UT1(I)/TEF1)

	SRR2(I)=-2./3.*XMUR(I)/RE*(DUTDT2(I)+2.*UR2(I)+UT2(I)/TEF2)

	SPP1(I)=2./3.*XMUR(I)/RE*(-DUTDT1(I)+UR1(I)+2.*UT1(I)/TEF1)

	SPP2(I)=2./3.*XMUR(I)/RE*(-DUTDT2(1)+UR2(I)+2.*UT2(I)/TEF2)

	SRT1(I)=XMUR(I)/RE*(DURDT1(I)-UT1(I))

	SRT2(I)=XMUL(I)/RE*(DURDT2(I)-UT2(I))

	! AVERAGING FOR FACES

	SRTT=(SRT1(I)+SRT2(I))/2.

	STTT=(STT1(I)+STT2(I))/2.

	SPPT=(SPP1(I)+SPP2(I))/2.

	SRRT=(SRR1(I)+SRR2(I))/2.
!------------------------------------------------------------------------- 
	!! UPDATE T-MOM

	!! ADD VISCOUS TERMS
	
	FMTV(I)=(STT1(I)-STT2(I))/DTETA(I)

	! VISCOUS SOURCE	

	SMTV(I)=2*SRTT+(STTT-SPPT)/TAN(TETA(I))

      RESMTT(I)=FMTV(I)+SMTV(I)+RESMT(I)

	FT(I)=FTO(I)+DT(I)*RESMTT(I)

	UT(I)=FT(I)/RHO(I)
!------------------------------------------------------------------------- 
	!! UPDATE R-MOM

	FMRV(I)=(SRT1(I)-SRT2(I))/DTETA(I)

	! VISCOUS R-MOM SOURCE

      SMRV(I)=SRRT-STTT-SPPT+SRTT/TAN(TETA(I))

      RESMRT(I)=FMRV(I)+SMRV(I)+RESMR(I)

      FR(I)=FRO(I)+DT(I)*RESMRT(I)

	UR(I)=FR(I)/RHO(I)
!------------------------------------------------------------------------- 
	!! UPDATE ENERGY

	DTDT1(I)=(TOLD(I+1)-TOLD(I))/DTETA(I)

	DTDT2(I)=(TOLD(I)-TOLD(I-1))/DTETA(I)

	Q1=-XMUR(I)/((GAMA-1.)*XMIN**2*RE*PRA)*DTDT1(I)

	Q2=-XMUL(I)/((GAMA-1.)*XMIN**2*RE*PRA)*DTDT2(I)

	Q(I)=Q1-Q2

	QF(I)=(Q1+Q2)/2.

	FEV1=STT1(I)*UT1(I)+SRT1(I)*UR1(I)

	FEV2=STT2(I)*UT2(I)+SRT2(I)*UR2(I)
	
	FEV(I)=(FEV1-FEV2-Q(I))/DTETA(I)

	! VISCOS ENERGY SOURCE

	UTS=(UTFHAT1+UTFHAT2)/2. !(2.*UTO(I)+UTO(I+1)+UTO(I-1))/4.

	URS=(URFHAT1+URFHAT2)/2. !(2.*URO(I)+URO(I+1)+URO(I-1))/4.

	TEX=TAN(TETA(I))
      
	SEV(I)=(SRTT+STTT/TEX)*UTS + (SRRT+SRTT/TEX)*URS-(QF(I)/TEX)
	
      RESET(I)=FEV(I)+SEV(I)+RESE(I)

      FE(I)=FEO(I)+DT(I)*RESET(I)
	
	E(I)=FE(I)/RHO(I)  

	END DO
	
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
		
	DO I=2,N-1
	      
	SUM1(I)=1. !RHOO(I)
	SUM2(I)=1. !FTO(I)
	SUM3(I)=1. !FRO(I)
	SUM4(I)=1. !FEO(I)
	  
	END DO

	! PUT THE VALUES IN THE LAST TIME STEP SOLUTION

	DO I=2,N-1
      
	URO(I)=UR(I) 	      

	UTO(I)=UT(I)       

      RHOO(I)=RHO(I)

	FTO(I)=FT(I)

	FRO(I)=FR(I)	

	FEO(I)=FE(I)
	
	EO(I)=E(I)
	
	END DO

	END SUBROUTINE	
!****************************************************************************
!     CALCULATE PARAMETERS
!
!****************************************************************************
      subroutine CALC
      	
	DO I=2,N-1
	
      P(I)=(E(I)-(UR(I)**2+UT(I)**2)/2.)*((GAMA-1)*RHO(I))

	PO(I)=P(I)

	T(I)=(E(I)-(UR(I)**2+UT(I)**2)/2.)/CV*(UI**2/TI)	 	

	TOLD(I)=T(I)

	C(I)=SQRT(GAMA*R*T(I)*TI)/UI 

	CO(I)=C(I)

	XM(I)=UT(I)/C(I)	
	
	H(I)=E(I)+P(I)/RHO(I) !P(I)/(R*RHO(I))*UI**2/TI 

	HO(I)=H(I)      
      
	END DO		         

	END SUBROUTINE     
!****************************************************************************
!     CHECK CONVERGENCE
!
!****************************************************************************
      subroutine CONVERGENCE 

	DO I=2,N-1
	
	RESC(I)=ABS(RESC(I)/SUM1(I))

	RESMTT(I)=ABS(RESMTT(I)/SUM2(I))

      RESMRT(I)=ABS(RESMRT(I)/SUM3(I))

      RESET(I)=ABS(RESET(I)/SUM4(I))

	END DO
	
	RESCMAX=RESC(3)
	
	RESEMAX=RESET(3)

	RESMRMAX=RESMRT(3)

	RESMTMAX=RESMTT(3)
	
	DO I=2,N-1
       
      IF (RESC(I).GT.RESCMAX)   RESCMAX=RESC(I)

	IF (RESET(I).GT.RESEMAX)   RESEMAX=RESET(I)

	IF (RESMRT(I).GT.RESMRMAX) RESMRMAX=RESMRT(I)

	IF (RESMTT(I).GT.RESMTMAX) RESMTMAX=RESMTT(I)

	END DO

	E1=0.
	E2=0.
	E3=0.
	E4=0.

	DO I=2,N-1
      
	E1=E1+ABS(RESC(I))**2
	E2=E2+ABS(RESMTT(I))**2
	E3=E3+ABS(RESMRT(I))**2
	E4=E4+ABS(RESET(I))**2

	END DO

      E1=SQRT(E1)/(N-2)
	E2=SQRT(E2)/(N-2)
	E3=SQRT(E3)/(N-2)
	E4=SQRT(E4)/(N-2)

C	RESIDUAL=MAX(RESCMAX,RESEMAX,RESMRMAX,RESMTMAX)	

      RESIDUAL=MAX(E1,E2,E3,E4)	
	
	WRITE (*,*)  RESIDUAL 

	WRITE (2,*)  MM,RESIDUAL !,E1,E2,E3,E4 !RESCMAX,RESEMAX !,RESMRMAX !,RESMTMAX 
	
	MM=MM+1 
	
	TIME=TIME+DT(1)
!---------------------------------------------------------------------------
	OPEN (80,FILE='TIM.DAT')
	
      IF ( (MM.EQ.1000).OR.(MM.EQ.1.D4).OR.(MM.EQ.2.D4).OR.(MM.EQ.3.D4)
     &.OR.(MM.EQ.4.D4).OR.(MM.EQ.5.D4).OR.(MM.EQ.6.D4).OR.(MM.EQ.5.D3)
     &.OR.(MM.EQ.15.D3).OR.(MM.EQ.25.D3).OR.(MM.EQ.35.D3)) THEN
C	DO J=1,50,10	
      
C	WRITE (80,*) ,'VARIABLES= "TETA" "T" "P" ' 

	DO I=2,N-1

	WRITE(80,*) TIME,(TETA(I)-TETAC)*180./(3.14159265359d0),T(I)
C     &,P(I)*(RHOI*UI**2)/PIN
      
	END DO
	END IF
C	END DO

	END SUBROUTINE
!****************************************************************************
!     CALCULATE OUTPUT
!
!****************************************************************************
 	Subroutine Output 	  

	write(5,*),'VARIABLES="TETA" "MACH" "UR" "UT" "MACHT" "RHO" "P"
     &	 "T" "ANGLE"'

      IF (METHOD==1) THEN 

	DO I=1,N
	WRITE(10,*) UT(I)     !Open(10,FILE='UT.DAT')	
	WRITE(20,*) UR(I)  !Open(20,FILE='UR.DAT')
	WRITE(30,*)	RHO(I) !OPEN (30,FILE='RHOO.DAT') 
	WRITE(40,*) P(I)  !OPEN (40,FILE='P.DAT') 
	WRITE(42,*) T(I)  !OPEN (42,FILE='T.DAT') 
	END DO
	
	END IF

	!(UR(I)+UR(I-1))/2.,(UT(I)+UT(I-1))/2.
	DO I=2,N-1
	U(I)=SQRT ( UR(I)**2+UT(I)**2)	
	WRITE(1,*)  UR(I),UT(I),XM(I) 
	WRITE(6,*)  RHO(I),P(I),T(I)
	WRITE(5,*)  (TETA(I)-TETAC)*180./(3.14159265359d0),U(I)*XMIN
     &,UR(I),UT(I),XM(I),RHO(I)
     &,P(I)*(RHOI*UI**2)/PIN,T(I),(ATAN(UT(I)/UR(I))+TETA(I))
     &*180./(3.14159265359d0)
	END DO          

	End Subroutine Output
!*****************************************************************************
	End program 

