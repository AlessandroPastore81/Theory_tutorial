module fields

  use ws_param
  use densita


  real (pr) :: t0,x0,t1,t2,x1,x2,t3,x3
  real (pr) :: W,gamma,dmshb0
  integer :: ifecm,j2terms,ixtls

  character*4 force
  !
  !!
  !
  real (pr) :: a0,b0,a1,b1
  real (pr) :: a2,b2,a3,b3

  !
  real (pr) :: af1,af2,af3,af4,af5,af6
  real (pr) :: af7,af8,af9,af10,af11 
  real (pr) :: bf1,bf2,bf3,bf4,bf5,bf6,bf7,bf8
  real (pr) :: bf9,bf10,bf11,bf14

contains

  subroutine buildfield
    implicit none
  end subroutine buildfield


  !------------------------------------------------------------------
  subroutine set_force

    !this routine fix the parameters of the Skyrme force
    ! input and output are passed with the module skyrmepar

    implicit none

    gamma = 1.0d0
    dmshb0 = 1.0d0 / 20.7355300000d0
    w = 0.0d0
    !
    !! These three parameters are for testing purposes only...
    !

    !NOTA: if ifecm=0 it means that we make the one body c.o.m correction
    ! dmshb = dmshb0 / ( 1- 1 /A )
    !

    select case (trim(force))
       !
       
    case ("t0t3")
    
      j2terms = 0
       t0 = -1132.400
       x0 =     0.d0
       t1 =     0.d0
       x1 =     0.000000d0
       t2 =     0.d0
       x2 =     0.000000d0
       t3 = 23610.40d0
       x3 =     0.000000d0
       w  =    0.000000d0
       gamma =   1.000000d0
       !
       ixtls = 0
       ifecm = 0
      
      
    case ("SIII")
       !
       !! SIII SET OF PARAMETERS
       !
       j2terms = 0
       !
       t0 = -1128.750000d0
       x0 =     0.450000d0
       t1 =   395.000000d0
       x1 =     0.000000d0
       t2 =   -95.000000d0
       x2 =     0.000000d0
       t3 = 14000.000000d0
       x3 =     1.000000d0
       w  =   120.000000d0
       gamma =   1.000000d0
       !
       ixtls = 0
       ifecm = 0

    case ("SLY4")
       !
       !! SLY4 SET OF PARAMETERS
       !
       j2terms = 0
       !
       t0 = -2488.913000d0
       x0 =     0.834000d0
       t1 =   486.818000d0
       x1 =    -0.344000d0
       t2 =  -546.395000d0
       x2 =    -1.000000d0
       t3 = 13777.000000d0
       x3 =     1.354000d0
       w  =   123.000000d0
       gamma =   1.0d0 / 6.0d0
       !
       dmshb0 = 1.0d0 / 20.7355300000d0
       !
       ixtls = 0
       ifecm = 0

    case ("SLY5")
       !
       !! SLY5 SET OF PARAMETERS
       !
       j2terms = 1
       !
       t0 = -2483.450000d0
       x0 =     0.776000d0
       t1 =   484.230000d0
       x1 =    -0.317000d0
       t2 =  -556.690000d0
       x2 =    -1.000000d0
       t3 = 13757.000000d0
       x3 =     1.263000d0
       w  =   125.000000d0
       !
       gamma =   1.0d0 / 6.0d0
       !
       dmshb0 = 1.0d0 / 20.7355300000d0
       !
       ixtls = 0
       ifecm = 0
       !

    case ("SKM*")
       !
       j2terms = 0
       !
       t0 = -2645.00d0
       x0 =     0.09d0
       t1 =   410.00d0
       x1 =     0.00d0
       t2 =  -135.00d0
       x2 =     0.00d0
       t3 = 15595.00d0
       x3 =     0.00d0
       w  =   130.00d0
       gamma = 1.0d0 / 6.0d0
       !
       ixtls = 0
       ifecm = 0
       !
       case ("WSso")
       !
       j2terms = 0
       t0 = -1132.400
       x0 =     0.d0
       t1 =     0.d0
       x1 =     0.000000d0
       t2 =     0.d0
       x2 =     0.000000d0
       t3 = 23610.40d0
       x3 =     0.000000d0
       w  =    0.000000d0
       gamma =   1.000000d0
       !
    case default
       !
       print *, "Unknown force..."
       stop
       !
    end select
    !
    call paramaux

    return
  end subroutine set_force

  subroutine paramaux
    ! this routine is used to initialize some parameters used with Skyrme
    ! no input required, the parameters are initialized in the module 
    ! skyrmepar

    implicit none



    a0 = t0 * ( 2.d0 + x0 ) / 4.d0
    b0 = - t0 * ( 2.d0 * x0 + 1.d0 ) / 4.d0
    a1 = ( t1 * ( 2.d0 + x1 ) + t2 * ( 2.d0 + x2 ) ) / 8.d0
    b1 = - ( t1 * ( 2.d0 * x1 + 1.d0 ) - t2 * ( 2.d0 * x2 + 1.d0 ) ) / 8.d0
    a2 = ( 3.d0 * t1 * ( 2.d0 + x1 ) - t2 * ( 2.d0 + x2 ) ) / 32.d0
    b2 = - ( 3.d0 * t1 * ( 2.d0 * x1 + 1.d0 ) + t2 * ( 2.d0 * x2 + 1.d0 ) ) / 32.d0
    a3 = t3 * ( 2.d0 + x3 ) / 24.d0
    b3 = - t3 * ( 2.d0 * x3 + 1.d0 ) / 24.d0

    !
    af1 = t0 * ( 2.d0 + x0 ) / 4.d0
    af2 = - t0 * ( 2.d0 * x0 + 1.d0 ) / 4.d0
    af3 = t3 * ( 2.d0 + x3 ) / 24.d0
    af4 = - t3 * ( 2.d0 * x3 + 1.d0 ) / 24.d0
    af5 = ( t1 * ( 2.d0 + x1 ) + t2 * ( 2.d0 + x2 ) ) / 8.d0
    af6 = - ( 3.d0 * t1 * ( 2.d0 + x1 ) - t2 * ( 2.d0 + x2 ) ) / 32.d0
    af7 = - ( t1 * ( 2.d0 * x1 + 1.d0 ) - t2 * ( 2.d0 * x2 + 1.d0 ) ) / 8.d0
    af8 = ( 3.d0 * t1 * ( 2.d0 * x1 + 1.d0 ) + t2 * ( 2.d0 * x2 + 1.d0 ) ) / 32.d0
    !
    if ( j2terms == 1 ) then
       af9 = - ( t1 * x1 + t2 * x2 ) / 16.d0
       af10 = ( t1 - t2 ) / 16.d0
    else
       af9 = 0.d0
       af10 = 0.d0
    end if
    !
    !!
    !
    af11 = - w / 2.d0
    !
    !!
    !
    bf1 = 2.d0 * af1
    bf2 = 2.d0 * af2
    bf3 = ( 2.d0 + gamma ) * af3
    bf4 = 2.d0 * af4
    bf5 = af5
    bf6 = 2.d0 * af6
    bf7 = af7
    bf8 = 2.d0 * af8
    bf9 = 2.d0 * af9
    bf10 = 2.d0 * af10
    bf11 = - af11
    bf14 = gamma * af4

    return
  end subroutine paramaux
  !
  !!
  !
  subroutine HFBpotentials(Ngrid1,del1,xmu,Nprot,Neutr,total_energy)
    !=========================================================
    !abstract: routine that calculates the HF potentials
    ! based on a skyrme interaction, it calculates the total
    ! energy using the functional
    !
    !input
    !Ngrid1,del1: box dimension and mesh
    !xmu: mixing among the old and new fields
    !Nprot,Neutr: N,Z
    !output
    !total_energy: total energy using the functional
    !========================================================

    implicit none

    integer :: isospin,ngrid1,i,Amass,Nprot,Neutr
    real (pr) :: del1,xmu,eps,r,coef,sum1,sum2,t13,t43
    real (pr) :: hb,ymu

    real (pr), allocatable :: rho_tot(:),rho_gamma(:)
    real (pr), allocatable :: tau_tot(:),cur_tot(:)
    real (pr), allocatable :: drho(:,:),dtau(:,:),dcur(:,:)
    real (pr), allocatable :: d2rho(:,:),d2tau(:,:),d2cur(:,:)
    real (pr), allocatable :: tmp(:),drho_tot(:),d2rho_tot(:)
    real (pr), allocatable :: dcur_tot(:)
    real (pr), allocatable :: VCx(:),HMENx(:,:),VNx(:,:),VSONx(:,:)
    real (pr), allocatable :: dhmenx(:,:),d2hmenx(:,:)

    real (pr) :: ex_coulomb_energy,coulomb_energy ! for the energy
    real (pr) :: kinetic_energy(0:1),kinetic_energy_tot
    real (pr) :: spin_orbit_energy,r2,dmshb
    real (pr) :: h_rho,h_rhotau,h_drho,h_gamma,h_so
    real (pr) :: field_energy, total_energy

    !Nota: i prefer to use new vectors
    eps = spacing( 1.0d0 )

    t13 = 1.d0/3.d0
    t43 = 4.d0/3.d0
    coef = - 0.75d0 * ( 3.d0 / pi )**t13 * echarg ! Slater coefficient Coulomb potential
    total_energy=0.d0
    ymu= 1.d0-xmu
    amass=Nprot+Neutr
    if ( ifecm == 0 ) then ! center of mass correction
       dmshb = dmshb0 / ( 1.0d0 - 1.0d0 / amass )
    else
       dmshb = dmshb0
    end if
    hb = 1.0d0 / dmshb

    allocate(rho_tot(ngrid1),rho_gamma(ngrid1), &
         tau_tot(ngrid1),cur_tot(ngrid1),                          &
         drho(ngrid1,0:1),dtau(ngrid1,0:1),dcur(ngrid1,0:1),       &
         d2rho(ngrid1,0:1),d2tau(ngrid1,0:1),d2cur(ngrid1,0:1),    &
         tmp(-1:ngrid1),drho_tot(ngrid1),d2rho_tot(ngrid1)         &
         ,dcur_tot(ngrid1))

    allocate(VCx(ngrid1),HMENx(ngrid1,0:1),VNx(ngrid1,0:1),VSONx(ngrid1,0:1), &
         dHMENx(ngrid1,0:1),d2HMENx(ngrid1,0:1))


    do i=1,ngrid1
       rho_tot(i)=rho(i,0)+rho(i,1)
       rho_gamma(i)=rho_tot(i)**gamma
       tau_tot(i)=tau(i,0)+tau(i,1)
       cur_tot(i)=cur(i,0)+cur(i,1)
    enddo

    !---- calculation of the derivatives ----------------------------
    !
    !! The densities rho and \tilde\rho are extrapolated to 0 using
    !! a forth order polynomial P(x) with the assumptions:
    !!   P'(0) = 0
    !!   P(-h) = P(h)
    !! Outside the box, they are set to 0.
    !
    !! For J and \tilde J, the extrapolation is:
    !!   P(0) = 0
    !!   P(-h) = -P(h)
    !
    !! The derivatives of the densities nabla(rho), nabla(\tilde\rho),
    !! Delta(rho) and Delta(\tilde\rho)
    !! are computed using a 5 points formula.
    !
    do isospin=0,1
       do i=1,ngrid1
          tmp(i) = rho(i,isospin)
       enddo
       tmp(-1) = tmp(1)
       tmp(0) = ( 15.d0 * tmp(1) - 6.d0 * tmp(2) + tmp(3) ) / 10.d0
       call d1_and_d2(ngrid1,del1, tmp, drho(:,isospin), d2rho(:,isospin))
       !
       do i=1,ngrid1
          tmp(i) = tau(i,isospin)
       enddo
       tmp(-1) = tmp(1)
       tmp(0) = ( 15.d0 * tmp(1) - 6.d0 * tmp(2) + tmp(3) ) / 10.d0
       call d1_and_d2(ngrid1,del1, tmp, dtau(:,isospin), d2tau(:,isospin))
       !
       do i=1,ngrid1
          tmp(i) = cur(i,isospin)
       enddo
       tmp(-1) = - tmp(1)
       tmp(0) = 0.0d0
       call d1_and_d2(ngrid1,del1, tmp, dcur(:,isospin),d2cur(:,isospin))
    enddo


    do i=1,ngrid1
       drho_tot(i)=drho(i,0)+drho(i,1)
       d2rho_tot(i)=d2rho(i,0)+d2rho(i,1)
       dcur_tot(i)=dcur(i,0)+dcur(i,1)
    enddo

    !--------------------------------------------------------------------
    ! Coulomb (direct potential)

    vcx=0.d0
    sum1 = 0.0d0  ! the funny integral is taken in Skalski paper PRC63 2001
    sum2 = 0.0d0
    do i = 1, ngrid1
       r=i*del1
       sum1 = sum1 + fourpi * r*r * rho(i,1)
       sum2 = sum2 + fourpi * r * rho(i,1)
       vcx(i) = sum1 / r - sum2
    end do

    do i=1,ngrid1
       vcx(i) = echarg * del1 * ( vcx(i) + sum2 )
    enddo

    !----- here i calculate the energy using the energy functional  

    !
    !!............................................. Energy
    !
    if ( minval(rho_tot) < -1.E-10 ) then
       print '(" Negative density in subroutine ENERGY in points:")'
       do i = 1, ngrid1
          if ( rho_tot(i) < -1.E-10  ) print *, i
       end do
    end if
    !
    !!............................Exchange and  Direct Coulomb energy
    !
    coulomb_energy = 0.d0
    ex_coulomb_energy = 0.d0

    do i=1,ngrid1
       r=i*del1
       r2=r*r
       coulomb_energy =coulomb_energy +vcx(i) * rho(i,1) * r2  / 2.0d0
       ex_coulomb_energy = ex_coulomb_energy+ rho(i,1)**t43 * r2  * coef
    enddo

    !
    !!............................................. Kinetic energy
    kinetic_energy(0)=0.d0
    kinetic_energy(1)=0.d0
    do i=1,ngrid1
       r=i*del1
       r2=r*r
       kinetic_energy(0)=kinetic_energy(0)+tau(i,0)*hb*r2! neutrons
       kinetic_energy(1)=kinetic_energy(1)+tau(i,1)*hb*r2! protons
    enddo
    kinetic_energy_tot=kinetic_energy(0) + kinetic_energy(1)

    !
    !!............................................. Spin orbit energy
    !
    if ( ixtls == 0 ) then
       spin_orbit_energy =0.d0
       do i=1,ngrid1
          r=del1*i
          r2=r*r
          spin_orbit_energy = spin_orbit_energy + r2 * af11 * ( cur_tot(i) * drho_tot(i) &
               + cur(i,0) * drho(i,0) + cur(i,1) * drho(i,1) )
       enddo
    else
       stop 'errore'
    end if

    !
    !!............................................. Field energy
    !
    field_energy = 0.d0
    do i=1,ngrid1
       r=del1*i
       r2=r*r
       h_rho = a0 * rho_tot(i)**2 + b0 * ( rho(i,0)**2 + rho(i,1)**2 )
       !
       h_rhotau = a1 * rho_tot(i) * tau_tot(i)  &
            + b1 * ( rho(i,0) * tau(i,0) + rho(i,1) * tau(i,1) )
       !
       h_drho = a2 * drho_tot(i)**2 + b2 * ( drho(i,0)**2 + drho(i,1)**2 )
       !
       h_gamma = a3 * rho_gamma(i) * rho_tot(i)**2  &
            + b3 * rho_gamma(i) * ( rho(i,0)**2 + rho(i,1)**2 )
       !
       h_so = af9 * cur_tot(i)**2 + af10 * ( cur(i,0)**2 + cur(i,1)**2 )
       field_energy=field_energy+( r2 * ( h_rho + h_rhotau + h_drho   &
            + h_gamma + h_so ) )
    enddo

    !
    !!............................................. Rearrangement energy
    !

    !
    !!
    !
    kinetic_energy_tot     = kinetic_energy_tot     * del1 * fourpi ! 
    field_energy           = field_energy           * del1 * fourpi ! 
    spin_orbit_energy      = spin_orbit_energy      * del1 * fourpi ! 
    coulomb_energy         = coulomb_energy         * del1 * fourpi ! 
    ex_coulomb_energy      = ex_coulomb_energy      * del1 * fourpi ! 
    !

    total_energy = kinetic_energy_tot + field_energy   &
         + spin_orbit_energy + coulomb_energy + ex_coulomb_energy 
    !

    !
    !-- calculation of the fields ----------------------------------------
    !
    do isospin = 0,1
       do i=1,ngrid1 
          r=i*del1

          !
          !!.................................................... Central field
          !
          vnx(i,isospin) = rho_tot(i)* &  ! central field
               ( bf1 + bf3 * rho_gamma(i))      &
               + rho(i,isospin)                                          &
               * ( bf2 + bf4 * rho_gamma(i))      &
               + rho_gamma(i) * ( bf14 * ( rho(i,0)**2 + rho(i,1)**2 )   &
               ) / ( rho_tot(i) + eps )                                  &
               + bf5 * tau_tot(i)                                        &
               + bf6 * ( d2rho_tot(i) + 2.d0*drho_tot(i) / r )           &
               + bf7 * tau(i,isospin)                                    &
               + bf8 * ( d2rho(i,isospin) + 2.d0 * drho(i,isospin) / r ) & 
               + isospin * (vcx(i) + t43 * coef * rho(i,1)**t13 )

          !
          !!................................................. Effective mass
          !
          Hmenx(i,isospin) = bf5 * rho_tot(i) + bf7 * rho(i,isospin) + hb
          !
          !!................................................. Spin orbit
          !
          
          if ( ixtls == 0 ) then
             vnx(i,isospin) = vnx(i,isospin)               &
                  + bf11 * ( dcur_tot(i) + dcur(i,isospin) &
                  + 2.d0 * ( cur_tot(i) + cur(i,isospin) ) / r )
             vsonx(i,isospin) = bf9 * cur_tot(i) + bf10 * cur(i,isospin) &
                  - bf11 * ( drho_tot(i) + drho(i,isospin) )
          else
             stop 'error'
          end if
          ! to make compatible the spin orbit potential with the one used
          ! by the code in the routine boundary i have to multiply for the factor
          !                 -2/r
          vsonx(i,isospin)=-2.d0*vsonx(i,isospin)/r

          dhmenx(i,isospin)=bf5 * drho_tot(i)  + bf7 * drho(i,isospin)
          d2hmenx(i,isospin)=bf5 * d2rho_tot(i)  + bf7 * d2rho(i,isospin)

       enddo
       do i=1,ngrid1
          r=i*del1
          vnx(i,isospin)=vnx(i,isospin)+ dhmenx(i,isospin)/r
       enddo
    end do

    !putting together Direct Coulomb and exchange
    do i=1,ngrid1
       vcx(i)=vcx(i) + t43 * coef * rho(i,1)**t13 
    enddo
    !--- mixing the old filed with the new one

!    rewind(1000)
!    rewind(1001)
    do isospin=0,1
       do i=1,ngrid1
          r=i*del1
          !if(ifixws.eq.1)then
             Vpot(i,isospin) = xmu*vnx(i,isospin)+ymu*vpot(i,isospin)
             vso(i,isospin) = xmu*vsonx(i,isospin)+ymu*vso(i,isospin)
             hb2m(i,isospin) = xmu*hmenx(i,isospin)+ymu*hb2m(i,isospin)
             dhmen(i,isospin) = xmu*dhmenx(i,isospin)+ymu*dhmen(i,isospin)
             d2hmen(i,isospin) = xmu*d2hmenx(i,isospin)+ymu*d2hmen(i,isospin)
!             write(1000+isospin,*)r,Vnx(i,isospin),vsonx(i,isospin),hmenx(i,isospin),dhmenx(i,isospin),d2hmenx(i,isospin)
           
             if(isospin.eq.1)then
                vc(i)=xmu*vcx(i)+ymu*vc(i)
             endif
          !endif
       enddo
    enddo
    !11  format(1x,I6,3(1x,F20.15))

    !---- deallocating the memory ------------------------------
    deallocate(VCx,HMENx,VNx,VSONx,dHmenx,d2Hmenx)

    deallocate(rho_tot,rho_gamma,tau_tot,cur_tot,   &
         drho,dtau,dcur,d2rho,d2tau,d2cur,tmp      &
         ,drho_tot,d2rho_tot,dcur_tot)
    return
  end subroutine HFBpotentials

!
!!
!
end module fields
