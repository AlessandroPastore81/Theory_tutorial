module ws_param


  integer, parameter :: pr = selected_real_kind( p = 12 )
  real (pr), parameter :: hbar=197.326968_pr ! [MeVfm]
  real (pr), parameter :: massnucleon =(938.2720_pr+939.5653_pr)/2.0_pr
  real (pr), parameter :: mesh=0.1_pr
  integer, parameter :: Nmaxstate=10000

  real (pr), allocatable :: Vpot(:,:),Psi(:,:,:),Psi1(:,:,:)
  real (pr), allocatable :: hb2m(:,:),Vso(:,:),dhmen(:,:)
  real (pr), allocatable :: VC(:),d2hmen(:,:)
  integer :: Npoints

  real (pr), parameter :: pi=2.0_pr*asin(1.0_pr)
  real (pr), parameter :: fourpi =4.0_pr*pi

  real (pr), parameter :: echarg = 1.4399784_pr 

  real (pr) :: EE(Nmaxstate,0:1),deg(Nmaxstate,0:1)
  integer :: Lused(Nmaxstate,0:1),Jused(Nmaxstate,0:1)


  !


contains


  subroutine zero

    implicit none

    EE=0.0_pr
    Lused=0; Jused=0

    return

  end subroutine zero

  !
  !!
  !

  subroutine alloc

    implicit none

    ! potential
    allocate(Vpot(Npoints,0:1),hb2m(Npoints,0:1),Vso(Npoints,0:1))
    allocate(dhmen(Npoints,0:1), VC(Npoints),d2hmen(Npoints,0:1))
    ! w.f.
    allocate(Psi(Npoints,Nmaxstate,0:1))
    allocate(Psi1(Npoints,Nmaxstate,0:1))


    return
  end subroutine alloc

  !
  !!
  !
  subroutine dealloc
    implicit none
  
    deallocate(dhmen,VC,d2hmen)
    deallocate(Vpot,Psi,Psi1,hb2m,Vso)
    !
   
    return
  end subroutine dealloc

  !
  !!
  !
  subroutine ordina(Ncount)
    implicit none

    integer :: Ncount(0:1),it,i,j
    integer :: laux,jaux,ix
    real (pr) :: t,px

    do it=0,1

       do  i=1,Ncount(it)-1         
          do j=i,1,-1
             if (EE(j+1,it).lt.EE(j,it)) then

                t=EE(j,it)
                EE(j,it)=EE(j+1,it)
                EE(j+1,it)=t

                laux=Lused(j,it)
                Lused(j,it)=Lused(j+1,it)
                Lused(j+1,it)=laux

                jaux=Jused(j,it)
                Jused(j,it)=Jused(j+1,it)
                Jused(j+1,it)=jaux

                do ix=1,Npoints
                   px=psi(ix,j,it)
                   psi(ix,j,it)=psi(ix,j+1,it)
                   psi(ix,j+1,it)=px
                enddo

                do ix=1,Npoints
                   px=psi1(ix,j,it)
                   psi1(ix,j,it)=psi1(ix,j+1,it)
                   psi1(ix,j+1,it)=px
                enddo
             end if
          enddo
       enddo

    enddo

    return
  end subroutine ordina

!
!!
!

  subroutine derivative(ncount,h)
    implicit none

    integer :: it,ncount(0:1),lla
    integer :: ircm1,n
    real (pr) :: h,u(npoints),up(npoints)

    PSi1=0.0_pr

    do it=0,1
       do n=1,ncount(it)

          lla=Lused(n,it)
          u=0.d0
          up=0.d0
          do ircm1=1,Npoints
             u(ircm1)=Psi(ircm1,n,it)
          enddo

          !call derivative_sca( u, up, Npoints ,h,lla)!( f, df,n,h, l )
          call derivative_sca2( u, up, Npoints ,h)!( f, df,n,h, l )

          do ircm1=1,Npoints
             PSi1(ircm1,n,it)=up(ircm1)
            !write(888+it,*)ircm1,up(ircm1),Psi(ircm1,n,it)
          enddo
       enddo
    enddo
    return
  end subroutine derivative




subroutine derivative_sca( f, df,n,h, l )
  ! subroutine to calculate the first derivative of the single particle
  ! wave functions
  !
  !input
  ! f : vector to be derived
  ! n: : number of points of the grid
  ! h: step of integration/derivation
  ! l: angular momentum of the single particle wavefunction
  !output
  !df : first derivative of f
  implicit none
  integer:: n,l
  double precision :: f(n), df(n), sig
  double precision :: h, h_12,h_60
  integer :: i
  !
  df=0.d0
  h_12 = 12.d0 * h 
  h_60 = 60.d0 * h
  sig = ( modulo(l,2) - 0.5d0 ) * 2
  df(1) = ( 8.0d0 * f(2) - f(3) + sig * f(1) ) / h_12
  df(2) = ( 45.d0 * ( f(3) - f(1) ) - 9.d0 * f(4)  &
       + f(5) - sig * f(1) ) / h_60
  df(3) = ( 45.d0 * ( f(4) - f(2) ) - 9.d0 * ( f(5) - f(1) )  &
       + f(6) ) / h_60


  df(n) = ( -8.d0 * f(n-1) + f(n) + f(n-2) ) / h_12
  df(n-1) = ( 45.d0 * ( f(n) - f(n-2) ) + 9.d0 * f(n-3) &
       - f(n) - f(n-4) ) / h_60
  df(n-2) = ( 45.d0 * ( f(n-1) - f(n-3) ) &
       - 9.d0 * ( f(n) - f(n-4) ) - f(n-5) ) / h_60

  do i = 4, n - 3
     df(i) = ( 45.d0 * ( f(i+1) - f(i-1) )  &
          - 9.d0 * ( f(i+2) - f(i-2) ) &
          + f(i+3) - f(i-3) ) / h_60
  end do
  !
end subroutine derivative_sca


subroutine derivative_sca2( f, df,n,h)
  ! subroutine to calculate the first derivative of the single particle
  ! wave functions
  !
  !input
  ! f : vector to be derived
  ! n: : number of points of the grid
  ! h: step of integration/derivation
  !output
  !df : first derivative of f
  implicit none
  integer:: n,j,k,jj,i
  double precision :: f(n), df(n)
  double precision :: h, h_12,h_60
  double precision:: A(5,5),EMFACT,sum
!!$
  DATA A(1,1),A(1,2),A(1,3),A(1,4),A(1,5) /-50.d0,96.d0,-72.d0,32.d0,-6.d0/
  DATA A(2,1),A(2,2),A(2,3),A(2,4),A(2,5) /-6.d0,-20.d0,36.d0,-12.d0,2.d0/
  DATA A(3,1),A(3,2),A(3,3),A(3,4),A(3,5) /2.d0,-16.d0,0.d0,16.d0,-2.d0/
  DATA A(4,1),A(4,2),A(4,3),A(4,4),A(4,5) /-2.d0,12.d0,-36.d0,20.d0,6.d0/
  DATA A(5,1),A(5,2),A(5,3),A(5,4),A(5,5) /6.d0,-32.d0,72.d0,-96.d0,50.d0/
  DATA EMFACT/24.d0/
!!$  !
  !
  df=0.d0
  h_12 = 12.d0 * h
  h_60 = 60.d0 * h

  DO  J=1,N  ! i use this routine taken from French H.F. Code only for second
     K=3     ! derivative, i tested that Bennaceur's routine and this one give
     IF(J.LT.3)K=J   ! the same results on first derivatives.
     IF(J.GT.N-2)K=J-N+5
     SUM=0.d0
     DO I=1,5
        JJ=J+I-K
        SUM=SUM+A(K,I)*F(JJ)
     ENDDO
     DF(J)=SUM/(H*EMFACT)
  ENDDO
  ! NOTA: notice that comparing to HFBRAD the first derivative at the
  ! edge is not done making assumption on the w.f. outsie the box, but
  ! interpolating it with a 4 point formula, this is maybe a little less
  ! precise but it makes a better choice when we change boundary conditions
  ! (u/r)'=0 it is more complicated
  !
  return
end subroutine derivative_sca2


!
!!
!
subroutine d1_and_d2( n,h,f, df, d2f )
  ! subrotuine to do the first and (second) derivatve
  ! of the densities 

  implicit none
  !
  !! This subroutine computes the first and second derivative of
  !! of function evaluated on the meshpoints 1,...,npt.
  !! The input is the function f with extrapolated values in -1, 0.
  !
  integer, intent(in) :: n
  double precision, intent(in) :: f(-1:n),h
  double precision, intent(inout) :: df(1:n)
  double precision, intent(inout) :: d2f(1:n)
  !    double precision, intent(inout), optional :: d2f(1:n)

  integer :: i
  double precision h_12,hh_12
  !
  h_12 = 12.d0 * h         ! \
  hh_12 = 12.d0 * h * h
  do i = 1, n - 2
     df(i) = ( 8.d0 * ( f(i+1) - f(i-1) ) - f(i+2) + f(i-2) ) / h_12
     !       if ( present(d2f) ) &
     d2f(i) = ( - 30.d0 * f(i) + 16.d0 * ( f(i+1) + f(i-1) ) &
          - f(i+2) - f(i-2) ) / hh_12
  end do
  !
  df(n-1) = ( - f(n-4) + 6.d0 * f(n-3) - 18.d0 * f(n-2) &
       + 10.d0 * f(n-1) + 3.d0 * f(n) ) / h_12
  df(n) = ( 3.d0 * f(n-4) - 16.d0 * f(n-3) + 36.d0 * f(n-2) &
       - 48.d0 * f(n-1) + 25.d0 * f(n) ) / h_12
  !    if ( present(d2f) ) then
  d2f(n-1) = ( - f(n-4) + 4.d0 * f(n-3) + 6.d0 * f(n-2) &
       - 20.d0 * f(n-1) + 11.d0 * f(n) ) / hh_12
  d2f(n) = ( 11.d0 * f(n-4) - 56.d0 * f(n-3) + 114.d0 * f(n-2) &
       - 104.d0 * f(n-1) + 35.d0 * f(n) ) / hh_12
  !    end if
  !
end subroutine d1_and_d2
!
!!
!
end module ws_param
