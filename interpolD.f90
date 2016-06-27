module interpol
  implicit none

  integer, private, parameter :: N_input=212! number of data points for total energies for the coupling elements
  integer, private, parameter :: N_states=2 ! number of states included in calculation
  integer, private, parameter :: N_coupl=1  ! number of coupling elements
  double precision :: R_grid(N_input)
  double precision :: Vdia(N_input,N_states),Vdia_NEW(N_input,N_states) ! total energy curves in the diabatic basis
  double precision :: coupl(N_input,N_coupl),coupl_NEW(N_input,N_coupl)

contains

!---------- reads input files for total energies
  subroutine read_input
    implicit none

    double precision :: dydR1,dydR2
    integer :: i,j

    !read in total diabatic energie surfaces
    open(98,file='Vdiabatic.inp',status='old')
    open(99,file='Coupldiabatic.inp',status='old')

    do i=1,N_input
       read(98,*) R_grid(i),Vdia(i,:)
       read(99,*) R_grid(i),coupl(i,:)
    end do

    close(98)
    close(99)

    !preprocess the total energy curves in the diabatic basis
    do i=1,N_states
       dydR1=(Vdia(2,i)-Vdia(1,i))/(R_grid(2)-R_grid(1))
       dydR2=(Vdia(N_input,i)-Vdia(N_input-1,i))/(R_grid(N_input)-R_grid(N_input-1))
       call spline(R_grid,Vdia(:,i),N_input,dydR1,dydR2,Vdia_NEW(:,i))
    end do

    !preprocess the diabatic coupling terms
    do i=1,N_coupl
       dydR1=(coupl(2,i)-coupl(1,i))/(R_grid(2)-R_grid(1))
       dydR2=(coupl(N_input,i)-coupl(N_input-1,i))/(R_grid(N_input)-R_grid(N_input-1))
       call spline(R_grid,coupl(:,i),N_input,dydR1,dydR2,coupl_NEW(:,i))
    end do

  end subroutine read_input


!---------- interpolates the diabatic coupling
  function dia_coupling(r,state)
    implicit none
    
    integer :: state
    double precision :: r,dia_coupling,help

    if (r.lt.R_grid(1)) then
!       print*, "WARNING: r is lower than smallest point in R_grid"
       dia_coupling=coupl(1,state)
       return
!       stop
    end if

    if (r.gt.R_grid(N_input)) then
       dia_coupling=coupl(N_input,state)
       return
    end if

    call splint(R_grid,coupl(:,state),coupl_NEW(:,state),N_input,r,help)
    dia_coupling=help

  end function dia_coupling


!---------- interpolates the total energy surfaces for the state "state" in the DIABATIC basis
  function energy_D(r,state)
    implicit none
    
    integer :: state
    double precision :: r,energy_D,help

    if (r.lt.R_grid(1)) then
!       print*, "WARNING: r is lower than smallest point in R_grid"
       energy_D=Vdia(1,state)
       return
!       stop
    end if

    if (r.gt.R_grid(N_input)) then
       energy_D=Vdia(N_input,state)
       return
    end if

    call splint(R_grid,Vdia(:,state),Vdia_NEW(:,state),N_input,r,help)
    energy_D=help

  end function energy_D


!---------- preprocesses the input data for interpolation with subroutine splint
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      DOUBLE PRECISION p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      end do
      return
      end SUBROUTINE spline


!---------- interpolation routine
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      DOUBLE PRECISION a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.d0) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
      return
      end SUBROUTINE splint


end module interpol
