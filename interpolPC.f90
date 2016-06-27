module interpol
  implicit none

  integer, private, parameter :: N_input=60 ! number of data points for total energies and the nonadiabatic coupling vector
  integer, private, parameter :: N_states=2 ! number of states included in calculation
  double precision :: R_grid(N_input)
  double precision :: Vad(N_input,N_states),Vad_NEW(N_input,N_states) ! total energy curves in the adiabatic basis
  double precision :: Vdia(N_input,N_states),Vdia_NEW(N_input,N_states) ! total energy curves in the diabatic basis

contains

!---------- reads input files for total energies
  subroutine read_input
    implicit none

    double precision :: dydR1,dydR2
    integer :: i,j

    !read in total energie curves and nonadiabatic coupling
    open(99,file='Vadiabatic.inp',status='old')
    open(98,file='Vdiabatic.inp',status='old')

    do i=1,N_input
       read(99,*) R_grid(i),Vad(i,:)
       read(98,*) R_grid(i),Vdia(i,:)
    end do

    close(99)

    !preprocess the total energy curves in the adiabatic basis
    do i=1,N_states
       dydR1=(Vad(2,i)-Vad(1,i))/(R_grid(2)-R_grid(1))
       dydR2=(Vad(N_input,i)-Vad(N_input-1,i))/(R_grid(N_input)-R_grid(N_input-1))
       call spline(R_grid,Vad(:,i),N_input,dydR1,dydR2,Vad_NEW(:,i))
    end do

    !preprocess the total energy curves in the diabatic basis
    do i=1,N_states
       dydR1=(Vdia(2,i)-Vdia(1,i))/(R_grid(2)-R_grid(1))
       dydR2=(Vdia(N_input,i)-Vdia(N_input-1,i))/(R_grid(N_input)-R_grid(N_input-1))
       call spline(R_grid,Vdia(:,i),N_input,dydR1,dydR2,Vdia_NEW(:,i))
    end do

  end subroutine read_input


!---------- interpolates the total energy curve for the state "state" in the ADIABATIC basis
  function energy_A(r,state)
    implicit none
    
    integer :: state
    double precision :: r,energy_A,help

    if (r.lt.R_grid(1)) then
!       print*, "WARNING: r is lower than smallest point in R_grid"
       energy_A=Vad(1,state)
       return
!       stop
    end if

    if (r.gt.R_grid(N_input)) then
       energy_A=Vad(N_input,state)
       return
    end if

    call splint(R_grid,Vad(:,state),Vad_NEW(:,state),N_input,r,help)
    energy_A=help

  end function energy_A


!---------- interpolates the total energy curve for the state "state" in the DIABATIC basis
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


!---------- calculates the off diagonal matrix elements of the hamiltonian in the diabatic basis
  function off_diag(r)
    implicit none

    double precision :: r,off_diag,help

    help = ( (energy_A(r,2)-(energy_D(r,1)+energy_D(r,2))/2.d0)**2 - (energy_D(r,1)-energy_D(r,2))**2/4.d0 )

    if (abs(help).lt.5.d-6) help = 0.d0
    if (help.lt.0.d0) help = 0.d0 ! using only help = abs(help) leads to a 3% change of the charge exchange probability

    off_diag = dsqrt(help)

  end function off_diag


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
