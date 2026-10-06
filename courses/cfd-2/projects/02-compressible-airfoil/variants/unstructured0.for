!****************************************************************************
!
!  PROGRAM: unstructured
!
!
!****************************************************************************

	program unstructured

	implicit none

	! Variables
	Integer :: NNode,NElement,NBound,Nt
	DoublePrecision :: T(4000)
	DoublePrecision :: Xnode(5000),Ynode(5000),Xc(4000),Yc(4000),Area(3912)
	DoublePrecision :: Tn(4000),dTdx(4000),dTdy(4000)
	Integer  :: NodeNum(4000,4),Bound(5000),Neib(4000,4),ElementRank(4000)
	DoublePrecision :: Kx,Ky,dt,Err,gen,ro,cp
	DoublePrecision :: Tbc1,Tbc2,qbc1,qbc2,h1,Tinf1


	! Body of unstructured
	
		Call Input

		Call Set_Patitions

		Call Initialize

		Call SetNeighbours

		Nt=0
		Do while (Nt<4000) !(Err>1d-5)

			Call CalcGradient

			Call Solve

			Err=maxval(abs(Tn-T))
			T=Tn
			Nt=Nt+1
			write(*,*) Nt,Nt*dt,Err
		end do

		Call Output

	contains

!****************************************************************************
	subroutine Input
	
		Integer :: i,Trash
		character*20 :: string

		Open(1,file='Mesh.node')

		i=1
		DO while(i<5000)
  			READ(1,*,END=101)  XNode(i), YNode(i)
			i=i+1
		END DO

101		NNode=i-1

		close(1)

		Open(2,file='Mesh.element')

		i=1
		DO  while(i<4000)
			READ(2,*,END=102) NodeNum(i,1), NodeNum(i,2)
     &			              ,NodeNum(i,3), NodeNum(i,4)
			i=i+1
		END DO

102		Nelement=i-1

		close(2)

	end subroutine Input
!****************************************************************************
	Subroutine Set_Patitions
		Integer :: i

		Open(3,file='Mesh.epart')

		DO i=1,Nelement
			READ(3,*) elementRank(i)
		END DO	
		close(3)

		!*********************************** Boundary Conditions
		
		Bound=0
		DO  i=1, NNode
  			if (XNode(i)==-0.035d0) then
				Bound(i)=1
			else if (XNode(i)==+0.035d0) then
				Bound(i)=2
			else if (YNode(i)==-0.025d0) then
				Bound(i)=11
			else if (YNode(i)==+0.025d0) then
				Bound(i)=21
  			else if ((XNode(i)**2+YNode(i)**2)<(0.01001)**2) then
				Bound(i)=12
			end if
			
		END DO
				

	end subroutine Set_Patitions
!****************************************************************************

	subroutine Initialize
		integer :: i,N1,N2,N3,N4
		double precision ::s
		
		Kx=237d0		!w/mK
		Ky=2d0*Kx
		ro=903d0			!kg/m3
		cp=2700d0			!J/kgK
		
		gen=5d0			!w/m3

		!*********************************** Boundary Conditions
		Tbc1=20d0
		Tbc2=50d0

		qbc1=0d0	!w/m2
		qbc2=5d0

		h1=50d0		!w/m2K
		Tinf1=200d0
		
		T=20d0
		Tn=0d0

		Err=1

		Xc=0d0
		Yc=0d0

		s=0
		Do i=1,NElement
			N1=NodeNum(i,1)
			N2=NodeNum(i,2)
			N3=NodeNum(i,3)
			N4=NodeNum(i,4)

			Area(i)=0.5d0*abs(
     &		(XNode(N3)-XNode(N1))*(YNode(N2)-YNode(N4))
     &		-(YNode(N3)-YNode(N1))*(XNode(N2)-XNode(N4)))

			Xc(i)=(XNode(N1)+XNode(N2)+XNode(N3)+XNode(N4))/4d0
			Yc(i)=(YNode(N1)+YNode(N2)+YNode(N3)+YNode(N4))/4d0
			s=s+Area(i)
		end do
		write(*,*) s
		
		dt=1d0*minval(Area)*ro*cp/max(kx,ky)

	end subroutine Initialize
!****************************************************************************
	subroutine SetNeighbours
		integer :: i,j,N1,N2,N3,N4

!		allocate(Neib(NElement,4))

		Neib=0

		Do i=1,NElement
			N1=NodeNum(i,1)
			N2=NodeNum(i,2)
			N3=NodeNum(i,3)
			N4=NodeNum(i,4)

	if (((Bound(N1)==0).or.(Bound(N2))==0).and.(Neib(i,1)==0)) then
			do j=1,NElement
				if (.not.(i==j)) then
					if ((N1==NodeNum(j,1)).and.(N2==NodeNum(j,4)))then
						Neib(i,1)=j
						Neib(j,4)=i
			else if ((N1==NodeNum(j,2)).and.(N2==NodeNum(j,1)))then
						Neib(i,1)=j
						Neib(j,1)=i
			else if ((N1==NodeNum(j,3)).and.(N2==NodeNum(j,2)))then
						Neib(i,1)=j
						Neib(j,2)=i
			else if ((N1==NodeNum(j,4)).and.(N2==NodeNum(j,3)))then
						Neib(i,1)=j
						Neib(j,3)=i
					end if
				end if
			end do
			end if

			if (((Bound(N2)==0).or.(Bound(N3))==0).and.(Neib(i,2)==0)) then
			do j=1,NElement
				if (.not.(i==j)) then
					if ((N2==NodeNum(j,1)).and.(N3==NodeNum(j,4)))then
						Neib(i,2)=j
						Neib(j,4)=i
					else if ((N2==NodeNum(j,2)).and.(N3==NodeNum(j,1)))then
						Neib(i,2)=j
						Neib(j,1)=i
					else if ((N2==NodeNum(j,3)).and.(N3==NodeNum(j,2)))then
						Neib(i,2)=j
						Neib(j,2)=i
					else if ((N2==NodeNum(j,4)).and.(N3==NodeNum(j,3)))then
						Neib(i,2)=j
						Neib(j,3)=i
					end if
				end if
			end do
			end if

			if (((Bound(N3)==0).or.(Bound(N4))==0).and.(Neib(i,3)==0)) then
			do j=1,NElement
				if (.not.(i==j)) then
					if ((N3==NodeNum(j,1)).and.(N4==NodeNum(j,4)))then
						Neib(i,3)=j
						Neib(j,4)=i
					else if ((N3==NodeNum(j,2)).and.(N4==NodeNum(j,1)))then
						Neib(i,3)=j
						Neib(j,1)=i
					else if ((N3==NodeNum(j,3)).and.(N4==NodeNum(j,2)))then
						Neib(i,3)=j
						Neib(j,2)=i
					else if ((N3==NodeNum(j,4)).and.(N4==NodeNum(j,3)))then
						Neib(i,3)=j
						Neib(j,3)=i
					end if
				end if
			end do
			end if

			if (((Bound(N4)==0).or.(Bound(N1))==0).and.(Neib(i,4)==0)) then
			do j=1,NElement
				if (.not.(i==j)) then
					if ((N4==NodeNum(j,1)).and.(N1==NodeNum(j,4)))then
						Neib(i,4)=j
						Neib(j,4)=i
					else if ((N4==NodeNum(j,2)).and.(N1==NodeNum(j,1)))then
						Neib(i,4)=j
						Neib(j,1)=i
					else if ((N4==NodeNum(j,3)).and.(N1==NodeNum(j,2)))then
						Neib(i,4)=j
						Neib(j,2)=i
					else if ((N4==NodeNum(j,4)).and.(N1==NodeNum(j,3)))then
						Neib(i,4)=j
						Neib(j,3)=i
					end if
				end if
			end do
			end if

		end do
		

	end subroutine SetNeighbours
!****************************************************************************
	subroutine CalcGradient

		integer :: i,j,n
		

		do i=1,NElement


			dTdx(i)=0d0
			dTdy(i)=0d0

			do j=1,4
				if (Neib(i,j)>0) then
					dTdx(i)=dTdx(i)-(T(i)+T(Neib(i,j)))/2d0		
     &		*(YNode(NodeNum(i,mod(j,4)+1))-YNode(NodeNum(i,j)))
					dTdy(i)=dTdy(i)+(T(i)+T(Neib(i,j)))/2d0		
     &		*(XNode(NodeNum(i,mod(j,4)+1))-XNode(NodeNum(i,j)))
				else

			!*********************************** Boundary Conditions
			if (Bound(NodeNum(i,j))==1) then

				dTdx(i)=dTdx(i)-Tbc1		
     &				*(YNode(NodeNum(i,mod(j,4)+1))-YNode(NodeNum(i,j)))
				dTdy(i)=dTdy(i)+Tbc1
     &				*(XNode(NodeNum(i,mod(j,4)+1))-XNode(NodeNum(i,j)))

			else if (Bound(NodeNum(i,j))==2) then

				dTdx(i)=dTdx(i)-Tbc2
     &				*(YNode(NodeNum(i,mod(j,4)+1))-YNode(NodeNum(i,j)))
				dTdy(i)=dTdy(i)+Tbc2
     &				*(XNode(NodeNum(i,mod(j,4)+1))-XNode(NodeNum(i,j)))

			else if ((Bound(NodeNum(i,j))==11)
     &					.or.(Bound(NodeNum(i,j))==12)) then

				dTdx(i)=dTdx(i)-T(i)
     &				*(YNode(NodeNum(i,mod(j,4)+1))-YNode(NodeNum(i,j)))
				dTdy(i)=dTdy(i)+T(i)
     &				*(XNode(NodeNum(i,mod(j,4)+1))-XNode(NodeNum(i,j)))

			else if (Bound(NodeNum(i,j))==21) then

				dTdx(i)=dTdx(i)-T(i)
     &				*(YNode(NodeNum(i,mod(j,4)+1))-YNode(NodeNum(i,j)))
				dTdy(i)=dTdy(i)+T(i)
     &				*(XNode(NodeNum(i,mod(j,4)+1))-XNode(NodeNum(i,j)))

			end if
		end if
				
			end do
			dTdx(i)=dTdx(i)/Area(i)
			dTdy(i)=dTdy(i)/Area(i)


		end do


	end Subroutine CalcGradient

!****************************************************************************
	subroutine 	Solve

		integer :: i,j
		double precision :: Fd

		do i=1,NElement	
	
			Fd=0
	
		do j=1,4
			if (Neib(i,j)>0) then
			Fd=Fd-Kx*(dTdx(i)+dTdx(Neib(i,j)))/2d0
     &				*(YNode(NodeNum(i,mod(j,4)+1))-YNode(NodeNum(i,j)))
     &			 +Ky*(dTdy(i)+dTdy(Neib(i,j)))/2d0
     &				*(XNode(NodeNum(i,mod(j,4)+1))-XNode(NodeNum(i,j)))
			else
				
			!*********************************** Boundary Conditions
			if (Bound(NodeNum(i,j))==1) then

			Tn(i)=Tbc1+dTdx(i)*((XNode(NodeNum(i,mod(j,4)+1))
     &								+XNode(NodeNum(i,j)))/2.0-Xc(i))
     &				 +dTdy(i)*((YNode(NodeNum(i,mod(j,4)+1))
     &								+YNode(NodeNum(i,j)))/2.0-Yc(i))

			goto 1001
			else if (Bound(NodeNum(i,j))==2) then

			Tn(i)=Tbc2+dTdx(i)*((XNode(NodeNum(i,mod(j,4)+1))
     &								+XNode(NodeNum(i,j)))/2.0-Xc(i))
     &				+dTdy(i)*((YNode(NodeNum(i,mod(j,4)+1))
     &								+YNode(NodeNum(i,j)))/2.0-Yc(i))

			goto 1001
			else if (Bound(NodeNum(i,j))==11) then

			Fd=Fd+qbc1*sqrt(
     &				(YNode(NodeNum(i,mod(j,4)+1))-YNode(NodeNum(i,j)))**2
     &				+(XNode(NodeNum(i,mod(j,4)+1))-XNode(NodeNum(i,j)))**2)

			else if (Bound(NodeNum(i,j))==12) then

			Fd=Fd+qbc2*sqrt(
     &				(YNode(NodeNum(i,mod(j,4)+1))-YNode(NodeNum(i,j)))**2
     &				+(XNode(NodeNum(i,mod(j,4)+1))-XNode(NodeNum(i,j)))**2)

			else if (Bound(NodeNum(i,j))==21) then

				Fd=Fd+h1*(Tinf1-T(i))*sqrt(
     &				(YNode(NodeNum(i,mod(j,4)+1))-YNode(NodeNum(i,j)))**2
     &				+(XNode(NodeNum(i,mod(j,4)+1))-XNode(NodeNum(i,j)))**2)

			end if

			end if

		end do

		Tn(i)=T(i)+dt/(ro*cp)*(Fd/Area(i)+gen)



1001		end do



	end subroutine Solve
!****************************************************************************

	Subroutine Output
		integer :: i,j,nn(NNode),rank(NNode)
		Double precision :: Tnode(NNode)

		
		nn=0
		Tnode=0d0
		do i=1,Nelement
			do j=1,4
				Tnode(NodeNum(i,j))=Tnode(NodeNum(i,j))+T(i)
				nn(NodeNum(i,j))=nn(NodeNum(i,j))+1
				rank(NodeNum(i,j))=ElementRank(i)
			end do
		end do
		Tnode=Tnode/nn

		do i=1,NNode
			if (Bound(i)==1) then
				Tnode(i)=Tbc1
			else if (Bound(i)==2) then
				Tnode(i)=Tbc2
			end if
		end do

		open(2,file='output.plt')
		write(2,*),'VARIABLES="X" "Y" "T" "part"'
		write(2,*),'ZONE F=FEPOINT,ET=quadrilateral,N='
     &							, NNode, ',E=', NElement


		do i=1,NNode
			write(2,*) Xnode(i),Ynode(i),Tnode(i),rank(i)
		end do

		do i=1,Nelement
			write(2,*) NodeNum(i,1),NodeNum(i,2),NodeNum(i,3),NodeNum(i,4)
		end do

		close(2)


	end Subroutine Output
!****************************************************************************
	
	


	end program unstructured

