module data_grid
  integer:: Nt          			! Zeitschritte
  integer:: Nstates                             ! Anzahl der elektron. Zustände
  integer:: NRX,NRY                             ! Anzahl der Gridpunkte im Ortsraum
  double precision:: RX_init,RY_init		! Zentrum des Wellenpaketes im Ortsraum
  double precision:: PX_init,PY_init            ! Zentrum des Wellenpaketes im Impulsraum
  double precision:: kappaX,kappaY		! FWHM des Wellenpaketes im Impulsraum
  double precision:: dt				! Zeitschritt
  double precision:: dpr, dr			! Auflösung des Impuls- und Ortsraumgrids
  double precision:: R0                         ! starting grid parameter
  double precision:: Rend                       ! ending grid parameter
  double precision, allocatable, dimension(:):: prx,pry 	! Impulsgrid
  double precision, allocatable, dimension(:,:,:):: Pot         ! Potentialgrid
end module data_grid

module data_au
  double precision,parameter:: pi=3.141592653589793d0        ! that's just pi
  double precision,parameter:: m1=7296.3299d0,m2=51195.823d0 ! He and Si mass in a.u.
  double precision,parameter:: Mass=m1*m2/(m1+m2)
  complex*16,parameter:: im=(0.d0,1.d0)   
end module data_au


program paul_dynamik

use data_grid
use data_au

 implicit none
 integer:: I
 
 double precision R
 double precision, allocatable, dimension(:,:):: mu	! Kopplungsmatrixelemente
 
  call input				! Einlesen von Parametern
  allocate(prx(NRX), pry(NRY), Pot(NRX,NRY,3))
  call p_grid				! Initialisiert das Impulsraumgrid

  write(*,*)
  write(*,*) 'INITIAL KINETIC ENERGY [a.u.] : ', (PX_init**2+PY_init**2)/Mass*0.5d0
  write(*,*) 'REDUCED MASS [a.u.] : ', Mass
  write(*,*)

  
  allocate(mu(nrx,nry))

  call potential(mu)			! Einlesen / Berechnen der Potentiale

  call propagation_2D(mu)	        ! Propagation

  deallocate (mu)
  deallocate (prx, pry, Pot)
        
  stop	  
end program paul_dynamik


! _______________ Subroutines __________________________________________________


subroutine input

use data_grid
use data_au

implicit none

 open(10,file='input',status='old')
 
  read(10,*) dt			  ! dt = time step (a.u.).
  read(10,*) Nt			  ! Nt = number of time steps.	
  read(10,*) RX_init,RY_init      ! RI = center of initial Gaussian (a.u.).
  read(10,*) kappaX,kappaY        ! kappa = FWHM of initial Gaussian in momentum space (a.u.).
  read(10,*) PX_init,PY_init      ! PI = center of inial Gaussian in momentum space (a.u.)
  read(10,*) NRX,NRY              ! NRX,NRY = number of real space grid points
  read(10,*) R0                   ! R0 = starting point of grid (a.u.)
  read(10,*) Rend                 ! Rend = end point of grid (a.u.)

  dR = (Rend - R0) / (NRX - 1)
 
  dpr = (2.d0 * pi) / (dR * NRX)  

  print*,'_________________________'
  print*
  print*,'Parameters'
  print*,'_________________________'
  print*
  print*,'dt = ', SNGL(dt), 'a.u.'
  print*,'dpR = ', SNGL(dpR), 'a.u.'
  print*,'dR = ', SNGL(dR), 'a.u.'
  print*
  print*,'__________________________'
  print*		  
                 
 end subroutine
 
!...................... Impulsgrid......................
 
subroutine p_grid

use data_grid
use data_au
 implicit none 
 integer:: I 
      
  dpr = (2.d0 * pi) / (dR * NRX)  

  do I = 1, NRX
    if (I.le.(NRX / 2)) then    
    PRX(I) = (I - 1) * dpR    
    else    
    PRX(I) = - (NRX + 1 - I) * dpR
    end if    
  end do
       
  do I = 1, NRY
    if (I.le.(NRY / 2)) then    
    PRY(I) = (I - 1) * dpR
    else    
    PRY(I) = - (NRY + 1 - I) * dpR    
    end if    
  end do
    
  return
end subroutine  
