  subroutine potential(mu)
         
  use data_grid
  use data_au
  use interpol
 
  implicit none 
  
  integer:: I, J
  double precision:: RX, RY, GM, dist
  double precision,intent(out):: mu(nrx,nry)  ! Kopplungselemente
 
  call read_input

  do i = 1, NRX
     do j = 1, NRY
        rx = r0 + (i-1) *dr
        ry = r0 + (j-1) *dr
        dist=sqrt(rx**2+ry**2)
        pot(i,j,1) = energy_D(dist,1)
        pot(i,j,2) = energy_D(dist,2)
     end do
  end do

 do i = 1, nrx
    do j = 1, nry
       rx = r0 + (i-1) *dr
       ry = r0 + (j-1) *dr
       mu(i,j) = dia_coupling(dsqrt(rx**2+ry**2),1)  ! Kopplung 1-2
    end do
  end do

!  GM = energy_D(24.d0,1)


  GM = 1.d6 
  
  do I = 1, NRX			! Wird auf 0 geshiftet (nicht notwendig)
     do J = 1, NRY
        if(pot(i,j,1).le.GM) then     
           GM = pot(I,J,1)	
        end if
     end do
   end do  	   
  
!   print*,'global minimum:', sngl(GM *au2eV)
   
  pot = pot - gm
 


!  open(10,file='potentiale.out',status='unknown')
!  open(23,file='kopplung.out',status='unknown')
!
!  do I = 1, NRX
!     do J = 1, NRY
!        RX = R0 + (I-1) * dR     			  
!        RY = R0 + (J-1) * dR
!        write(10,*) sngl(RX), sngl(RY), sngl(pot(i,j,:))
!        write(23,*) sngl(RX), sngl(RY), sngl(mu(i,j))  
!     end do
!     write(10,*)
!     write(23,*)
!  end do 
!  
!  close(10,status='keep')   
!  close(23,status='keep') 
  
return
end subroutine
        
            
