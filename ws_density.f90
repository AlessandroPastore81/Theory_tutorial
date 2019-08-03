module densita

  use ws_param

  implicit none

  real (pr), allocatable :: rho(:,:),tau(:,:),cur(:,:) 

contains


  subroutine allocatedensity

    implicit none

    allocate(rho(Npoints,0:1),tau(Npoints,0:1),cur(Npoints,0:1))

    return
  end subroutine allocatedensity

!
!!
!

  subroutine deallocatedensity

    implicit none

    deallocate(rho,tau,cur)

    return

  end subroutine deallocatedensity

!
!!
!
  subroutine occupation(Z,N,Ncount,pot)

    implicit none

    integer :: it,N,Z
    integer :: istate,istate2
    integer :: Ncount(0:1)
    real (pr) :: Naux,Ncheck(0:1),tol
    real (pr) :: frac,diff
    character*2 :: pot

    tol=0.000001_pr
    deg=0.0_pr
    Ncheck(0)=N
    Ncheck(1)=Z

    if (pot == 'WS') then  
    do it=0,1
       do istate=1,ncount(it)

          deg(istate,it)=1.0_pr

          Naux=0.0_pr
          do istate2=1,istate
             Naux=Naux+deg(istate2,it)*(jused(istate2,it)+1.0_pr)
 !               write(*,*)Naux,deg(istate2,it)*(jused(istate2,it)+1.0_pr),it,istate
          enddo
          ! write(*,*)'ciao'

          if(abs(Naux-Ncheck(it)).lt.tol) go to 100

          if(Naux.gt.Ncheck(it))then
             diff=Naux-Ncheck(it)
             frac=(jused(istate,it)+1.0_pr-diff)/(jused(istate,it)+1.0_pr)

             deg(istate,it)=frac
             !write(*,*) frac,Naux,Ncheck(it)
             if(frac.lt.0) stop 'Error in occupation'
             go to 100
          endif

       enddo

100    continue
    enddo
    else
    deg(:,:)=1.0_pr
    end if





    return
  end subroutine occupation

  !
  !!
  !
  subroutine builddensity(ncount,mesh)
    implicit none

    integer :: it,istate,ix
    integer :: Ncount(0:1)
    real (pr) :: r,mesh,ccc,lhat,eso,j12

    rho=0.0_pr
    tau=0.0_pr
    cur=0.0_pr

    do it=0,1
       do istate=1,ncount(it)
          lhat=lused(istate,it)*(lused(istate,it)+1.0_pr)
          j12=jused(istate,it)/2.0_pr
          eso=-(j12*(j12+1.0_pr)-lhat-3.0_pr/4.0_pr)
          ccc=deg(istate,it)*(jused(istate,it)+1.0_pr)
          do ix=1,Npoints
             r=ix*mesh
             rho(ix,it)=rho(ix,it)+ccc*psi(ix,istate,it)**2/(r**2*fourpi)
             !
             tau(ix,it)=tau(ix,it)+ccc*((psi1(ix,istate,it)-psi(ix,istate,it)/r)**2 &
                  +lhat*(psi(ix,istate,it)/r)**2)/(r**2*fourpi)
             !
            ! if(ccc.ne.0)then
            !  write(888+it,*)ix,((psi1(ix,istate,it)-psi(ix,istate,it)/r)**2 &
            !      +lhat*(psi(ix,istate,it)/r)**2)/(r**2*fourpi)
            !  endif 
             cur(ix,it)=cur(ix,it)+ccc*eso*psi(ix,istate,it)**2/(r**3*fourpi)
          enddo

       enddo
    enddo

    return
  end subroutine builddensity

!
!!
!
  subroutine plotdensity(mesh)

    implicit none

    integer :: ix,it
    real (pr) :: mesh,r,som(0:1)


    open(unit=10,file='DensityN.dat')
    open(unit=11,file='DensityP.dat')

    som=0.0_pr
    
    do it=0,1
       write(10+it,*)'#','r (fm)', 'Rho particle density (fm^-3)','Tau kinetic density (fm^-5)','J spin-orbit density (fm^-4)'
       do ix=1,Npoints
          r=ix*mesh
          write(10+it,*)r,rho(ix,it),tau(ix,it),cur(ix,it)
    !Check to see if the correct number of particle is obtained by integrating the particle density
          som(it)=som(it)+rho(ix,it)*fourpi*r**2*mesh
       enddo
    enddo
    !write(*,*)'Check =',som(0),som(1)

    close(10)
    close(11)
    return
  end subroutine plotdensity

end module densita
