program mainprogram

  use ws_param
  use ws_auxiliary
  use densita
  use fields
  use rspace


  implicit none

  integer :: lmax,it,ix,iter,itz,i,j
  integer :: Ncount(0:1),n,Proton,Neutron
  !
  integer, parameter :: itermax=1
  real (pr), parameter :: Emax=100.0_pr
  integer, parameter :: ialeroutine=0
  !
  real (pr) :: Rbox,diff, conv, cstso,cstcoul
  real (pr) :: Vdepth,Adiff,Radius,Naux, Vdn, Vdp
  real (pr) :: total_energy,xmu,total_energy_pre
  character*4 :: forceread
  character*2 :: pot
  logical (pr) :: w_so, w_coul, w_twodepth
  double precision, allocatable ::  PSaux(:,:)
  namelist /input/ Neutron, Proton, Rbox, forceread, lmax, Vdepth,Adiff,Radius, w_so,w_coul,pot,cstso,cstcoul, &
                    & w_twodepth, Vdp,Vdn


  ! ====== Initial parameters 

  diff = 0.0_pr
  Rbox=20.0_pr

  

  Vdepth=-51.d0 !-44.0192307692_pr   !51 - 33 (N-Z)/(N+Z) Pb208
  Adiff=0.67_pr
  Radius=3.200199d0  ! 7.52474001375_pr    !r0= 1.27; R = r0 A^1/3
  xmu=0.2_pr
  lmax= 6
  Proton=8
  Neutron=8
  forceread="t0t3"
  total_energy_pre=-9999.0_pr
  conv = 1.e-8
  w_so = .true.
  w_coul = .true.
  w_twodepth = .false.
  Vdn = Vdepth
  Vdp = Vdepth
  cstcoul = 1.0
  cstso = 1.0
!  pot = "WS"
  pot = "IW"

  
  open( unit=10, file ="input",status = "old" )
  read(10,nml=input)
  close(10)
  
   force=forceread
  Npoints=Nint(Rbox/mesh)


  ! zeroing vectors and allocation
  call zero
  call alloc

  call allocatedensity

  ! setting the  Skyrme interaction
  call set_force

  write(*,*)'Fixing the WS potential'
  call setWSpotential(Vdepth,Adiff,Radius,Proton,w_so,w_coul,pot,cstcoul,cstso,&
                  &   w_twodepth,Vdn,Vdp)

  write(*,*)'Solving Schroedinger equation'

  !write(*,*)
  !write(*,*)'-----------------------------------------------------'
  !write(*,*)'         IT       Etot                 Difference'
  !write(*,*)
  do iter=1,itermax
     ! I solve Schoreding equation to get w.functions

     if(ialeroutine.eq.0)then
        ! here I use the original Numerov form Milano
        allocate(PSaux(Nmaxstate,Npoints))
        itz=0
        call  boundary(mesh,Lmax,Emax, &
             Npoints,PSaux,Ncount(itz),Jused(:,itz),Lused(:,itz),EE(:,itz),itz)

        do n=1,Ncount(itz)
           do ix=1,Npoints
              PSi(ix,n,itz)=PSaux(n,ix)
           enddo
        enddo

        itz=1
        call  boundary(mesh,Lmax,Emax, &
             Npoints,PSaux,Ncount(itz),Jused(:,itz),Lused(:,itz),EE(:,itz),itz)
        do n=1,Ncount(itz)
           do ix=1,Npoints
              PSi(ix,n,itz)=PSaux(n,ix)
           enddo
        enddo

        deallocate(PSaux)
     else
        ! my routine
        call HFbasis(lmax,Ncount,Emax)
     endif

     ! derivative of the basis
     call derivative(Ncount,mesh)

     ! .. I sort them...
     call ordina(Ncount)

     ! .. I fill the degenracy

     call occupation(Proton,Neutron,Ncount,pot)

     ! -- i build the density 

     call builddensity(ncount,mesh)

     ! and I plot it
     call  plotdensity(mesh)
     do it=0,1
     Naux=0.0_pr
     do n=1,ncount(it)
        Naux=Naux+(jused(n,it)+1)*deg(n,it)
	!Rough estimate of the total energy by considering the sum or single-particles multiplied by degeneracy
	total_energy = total_energy+ ee(n,it)*(jused(n,it)+1)*deg(n,it)
     end do
     end do
!    if (pot == 'WS') print *, 'Estimate of the nucleus total energy:', total_energy
     if (pot == 'IW') then
     print *, 'Checking analytical formula for the eigenvalues'
     print *, '          n    ', 'Analitycal E_n' , '            Obtained  E_n'
     do i=1,10
     	print *,i, 20.7355300000d0*pi**2*i**2/Rbox**2, ee(i,0)
     end do 
     end if
     if (pot == 'FW') then
     print *, 'The occupied neutron state energies are:'
     print *,'Energy of single-particle (MeV)','     L ', '          2*J   ', 'Occupation'
     do it=0,0
     Naux=0.0_pr
     do n=1,ncount(0)
        Naux=Naux+(jused(n,it)+1)*deg(n,it)
        !if the state is occupied, print the energy
        if(deg(n,it).ne.0) then
        print *, ee(n,it),Lused(n,it),jused(n,it),Naux
        end if
!        do ix=1,Npoints
!           if(deg(n,it).ne.0)  write(999,*)ix,PSi(ix,n,it)
!        enddo
     enddo
     enddo
     end if
     if (pot == 'WS') then
     do it=0,1
     if (it==0) print *, 'The occupied neutron state energies are:'
     if (it==1) print *, 'The occupied proton state energies are:'
     print *,'Energy of single-particle (MeV)','     L ', '          2*J   ', 'Occupation'
     Naux=0.0_pr
     do n=1,ncount(it)
        Naux=Naux+(jused(n,it)+1)*deg(n,it)
        !if the state is occupied, print the energy
        if(deg(n,it).ne.0) then
        print *, ee(n,it),Lused(n,it),jused(n,it),Naux
        end if
!        do ix=1,Npoints
!           if(deg(n,it).ne.0)  write(999,*)ix,PSi(ix,n,it)
!        enddo
     enddo
     enddo
     end if
     open(unit=2000,file='wfs_test.dat')
     do j=1,6
     do i=1,Npoints
     	write(2000,*) i*mesh, sqrt(2.d0/Rbox)*sin(j*pi/Rbox*i*mesh)
     end do
        write(2000,*)
        write(2000,*)
      end do
     close(2000)

     ! create fields
     !call HFBpotentials(Npoints,mesh,xmu,Proton,Neutron,total_energy)
     !diff=abs(total_energy-total_energy_pre)

     !total_energy_pre=total_energy
     !write(*,*)iter,total_energy,diff
     !if(diff.lt.conv)exit
  enddo


  !=== writing on file

  open(unit=1000,file='neutron_singleparticles.dat')
  open(unit=1001,file='proton_singleparticles.dat')
  open(unit=2000,file='neutron_wfs.dat')
  open(unit=2001,file='protons_wfs.dat')
  write(1000+it,*) '#','Energy of single-particle (MeV)','Angular momentum L', 'J=L+S', 'Occupation'
  write(2000+it,*) '#','r (fm)','Wave function corresponding to energies in neutron_singleparticles.dat'
  do it=0,1
     Naux=0.0_pr
     do n=1,ncount(it)
        Naux=Naux+(jused(n,it)+1)*deg(n,it)
        !if the state is occupied, print the energy
        if(deg(n,it).ne.0) then
        write(1000+it,*)ee(n,it),Lused(n,it),jused(n,it),Naux
        !if the state is occupied, print the wave-funtion
        do ix=1,Npoints
        write(2000+it,*) ix*mesh,PSI(ix,n,it)
        end do
        write(2000+it,*)
        write(2000+it,*)
        end if
!        do ix=1,Npoints
!           if(deg(n,it).ne.0)  write(999,*)ix,PSi(ix,n,it)
!        enddo
     enddo
  enddo
  close(1000)
  close(1001)
  close(2000)
  close(2001)



  call dealloc
end program mainprogram


!
!!
!
