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
  !
  real (pr) :: Rbox,diff, conv
  real (pr) :: Vdepth,Adiff,Radius,Naux
  real (pr) :: total_energy,xmu,total_energy_pre
  character*4 :: forceread
  character*2 :: pot
  logical (pr) :: w_so, w_coul
  double precision, allocatable ::  PSaux(:,:)
  namelist /input/ Neutron, Proton, Rbox, forceread, lmax, Vdepth,Adiff,Radius, w_so,pot


  ! ====== Initial parameters 
  ! Do not touch here, only in input!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  diff = 0.0_pr
  Rbox=20.0_pr

  Npoints=Nint(Rbox/mesh)

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
!  pot = "WS"
  pot = "IW"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  open( unit=10, file ="input",status = "old" )
  read(10,nml=input)
  close(10)

   force=forceread



  ! zeroing vectors and allocation
call zero
call alloc

call allocatedensity

  ! setting the  Skyrme interaction
  ! not used here
!call set_force

write(*,*)'Setting the potential'
call setWSpotential(Vdepth,Adiff,Radius,Proton,w_so,w_coul,pot)

write(*,*)'Solving Schrodinger equation'
allocate(PSaux(Nmaxstate,Npoints)) ! This will be the wave-functions
	do itz = 0,1               ! This loop will be useful later to do protons and neutrons separately. For now, we do twice the same.
		!We use the Numerov algorithm
		call  boundary(mesh,Lmax,Emax, &
			Npoints,PSaux,Ncount(itz),Jused(:,itz), &
			Lused(:,itz),EE(:,itz),itz)
		
        	do n=1,Ncount(itz) !Ncount is the total number of states we keep after Numerov
			!We store the wave-functions for neutrons and protons           		
			do ix=1,Npoints 
				PSi(ix,n,itz)=PSaux(n,ix)
           		enddo
        	enddo
        enddo

        deallocate(PSaux)

! We do the derivative of the basis
call derivative(Ncount,mesh)

!We sort the energies from the smaller to the larger 
call ordina(Ncount)

!We fill the shells, for this exercise, we don't use it
call occupation(Proton,Neutron,Ncount,pot)

! We build the density of particles by considering the square of the wave-functions
call builddensity(ncount,mesh)

! This plots the density in the files DensityN.dat / DensityP.dat
call  plotdensity(mesh)
if (pot == 'IW') then
	print *, 'Checking analytical formula for the eigenvalues'
	print *, '          n    ', 'Analitycal E_n' , '            Obtained  E_n'
	!Here, you can program the analytical energies for an infinite well to check you have the same results.
	!Do it
!	do i=1,10
		print *, 'Compairison between analytical and obtained eigenvalues to fill'
!	end do 
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


  !=== writing on file
  !Here, we plot all the interesting quantities
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
