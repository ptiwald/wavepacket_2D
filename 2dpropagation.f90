subroutine propagation_2D(mu)
	
use data_grid
use data_au

 implicit none   
 include "/usr/include/fftw3.f"	 
 
 integer I, J, K, TIMEK
 integer L, M, N, FNAME
 integer*8 planF, planB

 double precision:: time		! Die Zeit
 double precision:: RX,RY		! Kernabstand
 double precision:: cpm			! Abschneideparameter d. cut-off Funktion
 double precision:: verw(2), terw(2), perwx(2), perwy(2) ! Erwartungswerte
 double precision:: norm(2), evRX(2), evRY(2)	         ! Population, Ortserwartungswert
 double precision:: in_norm(2),out_norm(2),outnorm_sum   ! Hilfsgrößen um die auslaufende Norm zu berechnen
 double precision,intent(in):: mu(nrx,nry)             ! Kopplungen, Schwingungswellen.

 double precision, allocatable, dimension(:,:):: cof	   ! Abschneidefunktion
 double precision, allocatable, dimension(:,:,:):: pdens   ! Impulsraumdichten

 complex*16:: tout(2,2)					! Diagonalisierte Kopplungsmatrix
 complex*16, allocatable, dimension(:,:,:):: psi_ges	! Gesamtwellenfunktion in el. Zustaenden
 complex*16, allocatable, dimension(:,:):: psi, kprop	! Hilfsfunktion, kinetischer Propagator
 
 open(98,file='cut_off.out',status='unknown')
 open(100,file='psi0.out',status='unknown')
 open(101,file='psi0_imp.out',status='unknown')	 
! open(200,file='dens_1.out',status='unknown')
! open(201,file='dens_2.out',status='unknown') 
! open(202,file='dens_3.out',status='unknown') 
! open(203,file='pdens_1.out',status='unknown')
! open(204,file='pdens_2.out',status='unknown') 	
! open(205,file='pdens_3.out',status='unknown') 	
 open(800,file='R_erw.out',status='unknown')
 open(801,file='pot_erw.out',status='unknown')
 open(802,file='kin_erw.out',status='unknown')
 open(804,file='imp_erw.out',status='unknown')
 open(908,file='norm.out',status='unknown')
 open(909,file='outnorm.out',status='unknown')


 allocate(psi(NRX,NRY),kprop(NRX,NRY),psi_ges(NRX,NRY,2),cof(NRX,NRY))      
 allocate(pdens(NRX,NRY,2))

 print*
 print*,'Tuning FFTW...'

call dfftw_plan_dft_2d(planF, NRX, NRY, psi, psi, FFTW_FORWARD,FFTW_MEASURE)
call dfftw_plan_dft_2d(planB, NRX, NRY, psi, psi, FFTW_BACKWARD,FFTW_MEASURE)
  
             
 print*,'Done.'
 print*	  

 FNAME = 200
 psi_ges = (0.d0,0.d0)  
 out_norm = 0.d0

 do I = 1, NRX
    do J = 1, NRY
       RX = R0 + (I-1)* dR
       RY = R0 + (J-1)* dR
       kprop(I,J) = exp(-im *dt * (PRX(I)**2+PRY(J)**2)/(4.d0*Mass))  ! Propagator, kin. 
       psi_ges(I,J,2) = exp(- (1.3863d0/kappaX**2) * (RX - RX_init)**2)&
                      * exp(- (1.3863d0/kappaY**2) * (RY - RY_init)**2)&
                      * exp(im * PX_init * (RX-RX_init)) * exp(im * PY_init * (RY-RY_init))
    end do
 end do


! cpm = 3.d0  !................ Definition einer cut_off-Funktion (kreisförmig um den ursprung!!!)
! 
! do I = 1, NRX
!    do J = 1, NRY
!       RX = R0 + (I - 1) * dR
!       RY = R0 + (J - 1) * dR
!       if ((sqrt(RX**2).lt.(Rend - cpm)).AND.(sqrt(RY**2).lt.(Rend - cpm))) then
!          cof(I,J) = 1.d0
!       else
!          if ((sqrt(RX**2).gt.(Rend - cpm)).AND.(sqrt(RY**2).lt.(Rend - cpm))) &
!               cof(I,J) = cos(((sqrt(RX**2) - Rend + cpm) / -cpm) * (0.5d0 * pi))
!          if ((sqrt(RY**2).gt.(Rend - cpm)).AND.(sqrt(RX**2).lt.(Rend - cpm))) &
!               cof(I,J) = cos(((sqrt(RY**2) - Rend + cpm) / -cpm) * (0.5d0 * pi))
!          if ((sqrt(RX**2).gt.(Rend - cpm)).AND.(sqrt(RY**2).gt.(Rend - cpm))) &
!               cof(I,J) = cos(((sqrt(RX**2) - Rend + cpm) / -cpm) * (0.5d0 * pi))&
!                         *cos(((sqrt(RY**2) - Rend + cpm) / -cpm) * (0.5d0 * pi))
!          cof(I,J) = cof(I,J)**2
!       end if
!    end do
! end do


!!! SMOOTHER CUT-OFF FUNCTION FOR LOW ENERGY PARTICLES
 cpm = 5.d0  !................ Definition einer cut_off-Funktion
 
 do I = 1, NRX
    do J = 1, NRY
       RX = R0 + (I - 1) * dR
       RY = R0 + (J - 1) * dR
       if ((sqrt(RX**2).lt.(Rend - cpm)).AND.(sqrt(RY**2).lt.(Rend - cpm))) then
          cof(I,J) = 1.d0
       else
          if ((sqrt(RX**2).gt.(Rend - cpm)).AND.(sqrt(RY**2).lt.(Rend - cpm))) then
             cof(I,J) = (sin(((sqrt(RX**2) - Rend + cpm) / -cpm) * (0.5d0 * pi)))**4
             cof(I,J) = 1.d0 - cof(I,J)
          endif
          if ((sqrt(RY**2).gt.(Rend - cpm)).AND.(sqrt(RX**2).lt.(Rend - cpm))) then
             cof(I,J) = (sin(((sqrt(RY**2) - Rend + cpm) / -cpm) * (0.5d0 * pi)))**4
             cof(I,J) = 1.d0 - cof(I,J)
          end if
          if ((sqrt(RX**2).gt.(Rend - cpm)).AND.(sqrt(RY**2).gt.(Rend - cpm))) then
             cof(I,J) = (1.d0 - sin(((sqrt(RX**2) - Rend + cpm) / -cpm) * (0.5d0 * pi))**4)&
                       *(1.d0 - sin(((sqrt(RY**2) - Rend + cpm) / -cpm) * (0.5d0 * pi))**4)
          end if
       end if
    end do
 end do
  

 call integ(psi_ges, norm)	 !... Normierung
 print*,'norm =', sngl(norm) 
 	 
 psi_ges(:,:,2) = psi_ges(:,:,2) / sqrt(norm(2))
 
 call integ(psi_ges, norm)
 print*,'norm =', sngl(norm) 


 do I = 1, NRX
    do J = 1, NRY
       RX = R0 + (I-1) * dR
       RY = R0 + (J-1) * dR
       write(100,*) sngl(RX),sngl(RY),sngl(abs(psi_ges(I,J,2))**2) ! Anfangszustand  
       write(98,*)  sngl(RX),sngl(RY),sngl(cof(i,j))		   ! Abschneidefunktion
    end do
    write(100,*)
    write(98,*)
 end do

    do I = 1,NRX
       do J = 1,NRY
          psi(I,J) = psi_ges(I,J,2)   ! Hilfsgroesse, Kopie in 2D Groesse
       end do
    end do
    call dfftw_execute(planF,psi,psi)	 ! psi(r) --> psi(p)
    psi = psi / sqrt(dble(NRX*NRY))      ! Normierung nach FT

    do I = 1, NRX
       do J = 1, NRY
          write(101,*) sngl(PRX(I)),sngl(PRY(J)),sngl(abs(psi(I,J))**2) ! Anfangszustand  
       end do
       write(101,*)
    end do


!______________________________________________________________________
!                                                                     
!                   Propagation Loop                         
!______________________________________________________________________
                             

 print*
 print*,'2D propagation...'
 print*
	

timeloop: do TIMEK = 1, Nt

    time = TIMEK * dt

    if(mod(TIMEK,200).eq.0) then
      print*,'time:', sngl(time)	     
    end if


      terw = 0.d0
      perwx = 0.d0 
      perwy = 0.d0
      evRX = 0.d0
      evRY = 0.d0
      verw = 0.d0
 
  
    
!----------------- Propagation via Split Operator 3rd Ordnung ---------
!
!   
! ............... kinetische Propagation, 1. Teil ...................

    do K = 1,2	
      do I = 1,NRX
         do J = 1,NRY
            psi(I,J) = psi_ges(I,J,K)   ! Hilfsgroesse, Kopie in 1D Groesse
         end do
      end do 	

      call dfftw_execute(planF,psi,psi)	 ! psi(r) --> psi(p)

      psi = psi * kprop            ! kin. Propagation

      call dfftw_execute(planB,psi,psi)	 ! psi(p) --> psi(r)
      psi = psi / dble(NRX*NRY)	         ! Normieren, wg. numerischer FT
	 
      do I = 1,NRX
         do J = 1,NRY
            psi_ges(I,J,K) = psi(I,J)   ! Hilfsgroesse, zurueck kopieren
         end do
      end do      
    end do

! ............... potentielle Propagation ..........................
  
    
   do i = 1, NRX				
      do j = 1, NRY                             ! Kopplung zwischen den el. Zustaenden,
         call pulse2(tout, mu(i,j)) 		! geht ein in die off-Diagonalen -- hier
         psi_ges(i,j,1:2) = matmul(tout(1:2,1:2),psi_ges(i,j,1:2))  ! wird diagonalisiert	
      end do
   end do


   do k = 1, 2   
      do i = 1, NRX 
         do j = 1, NRY
            psi_ges(i,j,k) = psi_ges(i,j,k) * exp(-im * dt * pot(i,j,k)) ! potentielle Propagation    
         end do
      end do
   end do
  
! ............... kinetische Propagation, 2. Teil ...................  


    do K = 1,2	
      do I = 1, NRX
         do J = 1, NRY
            psi(I,J) = psi_ges(I,J,K)  	! Hilfsgroesse
         end do
      end do 
      call dfftw_execute(planF,psi,psi) ! psi(r) --> psi(p)
      psi = psi * kprop		        ! kin. Propagation
      psi = psi / sqrt(dble(NRX*NRY))	! Normieren, wg. numerischer FT

       do I = 1, NRX
          do J = 1, NRY
             terw(k) = terw(k) + abs(psi(i,j))**2 * ((PRX(I)**2+PRY(J)**2)/(2.d0*Mass))! Erwartungswerte, kin. Energie
             perwx(k) = perwx(k) + abs(psi(i,j))**2 * PRX(I)	! Impulserwartungswerte X
             perwy(k) = perwy(k) + abs(psi(i,j))**2 * PRY(J)	! Impulserwartungswerte Y
             pdens(i,j,k) = abs(psi(i,j))**2			! Impulsraumdichte
          end do
       end do

      call dfftw_execute(planB,psi,psi) ! psi(p) --> psi(r)                  
      psi = psi / sqrt(dble(NRX*NRY))	! Normieren, wg. numerischer FT

      do I = 1, NRX
         do J = 1, NRY
            psi_ges(I,J,K) = psi(I,J)		! Hilfsgroesse, zurueck kopieren
         end do
      end do      
    end do     
 
!-------------------- Erwartungswerte und Output -------------------------------

 
    do k = 1, 2
     do i = 1, NRX
        do j = 1, NRY
           RX = R0 + (i-1) *dR
           RY = R0 + (j-1) *dR
           evRX(k) = evRX(k) + abs(psi_ges(i,j,k))**2 * RX  ! Erwartungswerte, Kernabstand X
           evRY(k) = evRY(k) + abs(psi_ges(i,j,k))**2 * RY  ! Erwartungswerte, Kernabstand X
           verw(k) = verw(k) + abs(psi_ges(i,j,k))**2 * pot(i,j,k) ! Erwartungswerte, pot. Energie
        end do
     end do   
    end do
    evRX = evRX * dR**2
    evRY = evRY * dR**2
    verw = verw *dR**2
    ! * dR**2 wegen fftw !!!
    terw = terw *dR**2
    perwx = perwx *dR**2
    perwy = perwy *dR**2
  
 
    call integ(psi_ges, norm)			! Norm in den elektr. Zustaenden


    do j = 1, 2
	if (norm(j).ge.1.d-10) then
	evRX(j) = evRX(j) / norm(j)
	evRY(j) = evRY(j) / norm(j)
        verw(j) = verw(j) / norm(j)
        terw(j) = terw(j) / norm(j)
        perwx(j) = perwx(j) / norm(j)
        perwy(j) = perwy(j) / norm(j)
	end if
    end do
             
     
    write(800,*) sngl(time), sngl(evRX ), sngl(evRY )
    write(908,*) sngl(time), sngl(norm)	  
    write(801,*) sngl(time), sngl(verw)
    write(802,*) sngl(time), sngl(terw)
    write(804,*) sngl(time), sngl(perwx), sngl(perwy)

      	
	
 if(mod(TIMEK,20).eq.0) then
! if ( time.eq.14100.d0 ) then

   if (fname.gt.300) exit 

  open(fname, status='unknown')
  do I = 1, NRX, 20  		! Dichten in den elektronischen Zustaenden
     do J = 1, NRY, 20
        RX = R0 + (I-1) *dR 	           
        RY = R0 + (J-1) *dR
        write(FNAME,*) sngl(RX), sngl(RY), sngl(abs(psi_ges(I,J,1)**2)),sngl(abs(psi_ges(I,J,2)**2))
     end do
     write(FNAME,*)
  end do
  FNAME = FNAME + 1

!  open(fname, status='unknown')
!   do I = NRX/2+1, NRX, 2       ! Impulsraumdichten in den elektr. Zustaenden
!      do J = NRY/2+1, NRY, 2
!         write(FNAME,*) sngl(prx(i)), sngl(pry(j)), sngl(pdens(i,j,1)), sngl(pdens(i,j,2))
!      end do
!      write(FNAME,*)
!   end do
!  close(fname, status='keep')
!  FNAME = FNAME + 1  
!
!  open(fname, status='unknown')
!   do I = NRX/2+1, NRX, 2 !I = 1, NRX/2, 1
!      do J = 1, NRY/2, 2        
!         write(FNAME,*) sngl(prx(i)), sngl(pry(j)), sngl(pdens(i,j,1)), sngl(pdens(i,j,2))  
!      end do
!      write(FNAME,*)
!   end do
!  close(fname, status='keep')
!  FNAME = FNAME + 1  
!
!  open(fname, status='unknown')
!   do I = 1, NRX/2, 2 !I = 1, NRX/2, 1
!      do J = 1, NRY/2, 2        
!         write(FNAME,*) sngl(prx(i)), sngl(pry(j)), sngl(pdens(i,j,1)), sngl(pdens(i,j,2))  
!      end do
!      write(FNAME,*)
!   end do
!  close(fname, status='keep')
!  FNAME = FNAME + 1  

!   do I = nr/2+1, Nr, 4 	! Impulsraumdichten in den elektr. Zustaenden
!    write(203,*) sngl(time), sngl(pr(i)), sngl(pdens(i,1))
!    write(204,*) sngl(time), sngl(pr(i)), sngl(pdens(i,2))
!    write(205,*) sngl(time), sngl(pr(i)), sngl(pdens(i,3))
!   end do
!   do I = 1, Nr/2, 4
!    write(203,*) sngl(time), sngl(pr(i)), sngl(pdens(i,1))  
!    write(204,*) sngl(time), sngl(pr(i)), sngl(pdens(i,2))
!    write(205,*) sngl(time), sngl(pr(i)), sngl(pdens(i,3))
!   end do
!   write(203,*)
!   write(204,*)
!   write(205,*)
!
 end if


 call integ(psi_ges,norm)

 do i = 1, NRX
    do j = 1, NRY
       psi_ges(i,j,1) = psi_ges(i,j,1) * cof(i,j)	! Abschneidefunktion; die 
       psi_ges(i,j,2) = psi_ges(i,j,2) * cof(i,j)	! rauslaufende Funktion wird abgeschnitten
    end do
 end do

 call integ(psi_ges,in_norm)

 out_norm = out_norm + (norm - in_norm)

 write(909,*) sngl(time), sngl(out_norm)
 
 !------------------------------------------------------
! check sum of out_norm and end program if norm_sum ~ 1
!------------------------------------------------------

 outnorm_sum = 0.d0
 do i = 1, 2
    outnorm_sum = outnorm_sum + out_norm(i)
 end do

 if (outnorm_sum.gt.0.99d0) then
    print*, 'out_norm reached 1'
    exit
 end if

end do timeloop

 close(98, status='keep')	    
 close(100, status='keep')	    
 close(101, status='keep')	    
 close(200, status='keep')
 close(201, status='keep')
 close(202, status='keep')
 close(203, status='keep')
 close(800, status='keep')
 close(801, status='keep')
 close(802, status='keep')
 close(804, status='keep')
 close(805, status='keep')
 close(806, status='keep')
 close(908, status='keep')
 close(909, status='keep')
 
	 
                                   
 call dfftw_destroy_plan(planF)
 call dfftw_destroy_plan(planB)
 	    
 deallocate(psi, kprop, psi_ges, cof, pdens)
 
return
end subroutine


!_________________________________________________________


subroutine integ(psi, norm)
      
use data_grid
 implicit none
 integer I, J, K
      
 double precision,intent(out):: norm(2)
 complex*16,intent(in):: psi(NRX,NRY,2)
      
 norm = 0.d0
 
 do K = 1, 2 ! Norm in den 2 Zustaenden
  do I = 1, NRX
     do J = 1, NRY
        norm(K)= norm(K) + abs(psi(I,J,K))**2     				  
     end do
  end do
 end do
 
   norm = norm * dR * dR
     
     
return 
end subroutine integ
!______________________________________________________________
subroutine integ_onepsi(psi, norm)
      
use data_grid
 implicit none
 integer I, J
      
 double precision,intent(out):: norm
 complex*16,intent(in):: psi(NRX,NRY)
      
 norm = 0.d0
 
 do I = 1, NRX
    do J = 1, NRY
       norm= norm + abs(psi(I,J))**2     				  
    end do
 end do
 
 norm = norm * dR * dR
     
     
return 
end subroutine integ_onepsi
!______________________________________________________________
subroutine pulse2(tout, coupl)

use data_grid, only:dt
use data_au, only:im

implicit none

 integer:: J,info
 double precision:: w(2,2), u(2,2), d(2)
 double precision, intent(in)::coupl
 
 complex*16:: b(2,2), z(2,2)
 complex*16, intent(out):: tout(2,2)

!  ........... Definition der Kopplungsmatrix
   
u(1,1) = 0.d0
u(1,2) = coupl

u(2,1) = coupl
u(2,2) = 0.d0

! call jacobi(u,2,d) 	! Diagonalisierung, numerisch
 call jacobi2(u,2,2,d,w,info)

b= (0.d0,0.d0)
 
do J = 1,2
  b(J,J) = cmplx(dcos(dt*d(J)),-dsin(dt*d(J)))
! b(J,J) = exp (-im * dt * d(J))
end do

!jacobi
!z = matmul(u,b)
!tout = matmul(z,transpose(u))

! jacobi2
z = matmul(w,b)
tout = matmul(z,transpose(w))



return
end subroutine


!!-------------------------


      subroutine jacobi (mat,dim,ewerte)

         implicit none

         double precision genau
         parameter (genau=1.d-15)

         integer Jmax,mmax
         parameter (Jmax=15,mmax=18)
         integer matdim
         parameter (matdim=3)

         double precision mat(matdim,matdim)
         integer dim
         double precision ewerte(matdim)

         double precision s(matdim,matdim)
         integer ca,cb,p,q
         double precision c1,c2,t1,t2,t3,v1,v2,v3
         double precision tmp,l,n,t,m1,w,m
         logical flag

	s= 0.d0
         !!!!call fillmt(s,dim,dim,0.d0,matdim,matdim)

         do 1 ca=1,dim,1
            s(ca,ca)=1.d0
1           continue

         l=0.d0
         do 2 ca=2,dim,1
            do 12 cb=1,dim,1
               tmp=mat(ca,cb)
               l=l+2.d0*tmp*tmp
12          continue
2        continue

         n=dsqrt(l)
         m=genau*n/dim
         t=n

3        t=t/dim
4           do 6 q=2,dim,1
               do 16 p=1,q-1,1
                  flag=.false.
                  if (dabs(mat(p,q)).gt.t) then
                     flag=.true.
                     v1=mat(p,p)
                     v2=mat(p,q)
                     v3=mat(q,q)
                     m1=(v1-v3)/2.d0
                     if (m1.eq.0.d0) then
                           w=-1.d0
                        else
                           if (m1.gt.0.d0) then
                                 w=-v2/(dsqrt(v2*v2+m1*m1))
                              else
                                 w=v2/(dsqrt(v2*v2+m1*m1))
                              endif
                        endif

                     t1=w/dsqrt(2.d0*(1+dsqrt(1.d0-w/2.d0)))
                     t2=t1*t1
                     c1=dsqrt(1.d0-t2)
                     c2=c1*c1
                     t3=t1*c1

                     do 7 ca=1,dim,1
                        l=mat(ca,p)*c1-mat(ca,q)*t1
                        mat(ca,q)=mat(ca,p)*t1+mat(ca,q)*c1
                        mat(ca,p)=l
                        l=s(ca,p)*c1-s(ca,q)*t1
                        s(ca,q)=s(ca,p)*t1+s(ca,q)*c1
                        s(ca,p)=l
7                       continue
                     do 8 ca=1,dim,1
                        mat(p,ca)=mat(ca,p)
                        mat(q,ca)=mat(ca,q)
8                       continue
                     mat(p,p)=v1*c2+v3*t2-2*v2*t3
                     mat(q,q)=v1*t2+v3*c2+2*v2*t3
                     tmp=(v1-v3)*t3+v2*(c2-t2)
                     mat(p,q)=tmp
                     mat(q,p)=tmp
                     end if
16                continue
6              continue
               if (flag) go to 4
            if (m.lt.t) go to 3
            ewerte=0.d0
         !!!call fillvc(ewerte,dim,0.d0)
         do 9 ca=1,dim,1
            ewerte(ca)=mat(ca,ca)
9           continue
         do 10 ca=1,dim,1
            do 11 cb=1,dim,1
               mat(ca,cb)=s(ca,cb)
11          continue
10       continue

         return
         end


!______________________________________________________

      SUBROUTINE jacobi2(a,n,np,d,v,nrot)
      INTEGER n,np,nrot,NMAX
      DOUBLE PRECISION a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=1000)
      INTEGER i,ip,iq,j
      DOUBLE PRECISION c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.d0
11      continue
        v(ip,ip)=1.d0
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.d0
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.d0
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))
14        continue
15      continue
        if(sm.eq.0.d0)return
        if(i.lt.4)then
          tresh=0.2d0*sm/n**2
        else
          tresh=0.d0
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.d0*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.d0
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5d0*h/a(ip,iq)
                t=1.d0/(abs(theta)+sqrt(1.d0+theta**2))
                if(theta.lt.0.d0)t=-t
              endif
              c=1.d0/sqrt(1+t**2)
              s=t*c
              tau=s/(1.d0+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.d0
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.d0
23      continue
24    continue
      pause 'too many iterations in jacobi'
      return
      END

