module ws_auxiliary
  use ws_param

contains


  subroutine setWSpotential(Vdepth,Adiff,Radius,Nprot,w_so,w_coul,pot,cstcoul,cstso,w_twodepth,Vdn,Vdp)

    implicit none

    integer :: i,itz,Nprot
    real (pr) :: Vdepth,Adiff,Radius, Vdp, Vdn
    real (pr) :: ff,r,cstcoul,cstso
    logical (pr) :: w_so, w_coul, w_twodepth
    character*2 :: pot



    dhmen=0.0_pr
    d2hmen=0.0_pr
    VC=0.0_pr
    Vso=0.0_pr 
    
    open(unit=1000,file='neutron_potential.dat')
    open(unit=1001,file='proton_potential.dat')
    
    do itz=0,1
       if (itz==0) write(1000,*)'#','r (fm)','Potential','Spin-orbit potential'
       if (itz==1) write(1001,*)'#','r (fm)','Potential','Spin-orbit potential', 'Coulomb'
       do i=1,npoints
          select case(pot)
          case("WS")
          	r=i*mesh
          	ff=exp(-(r-Radius)/Adiff)
          	Vpot(i,itz)=Vdepth*ff/(1.d0+ff)
                !Having a proton and a neutron depth
		if (w_twodepth) then
                if (itz==0) Vpot(i,itz)=Vdn*ff/(1.d0+ff)
                if (itz==1) Vpot(i,itz)=Vdp*ff/(1.d0+ff) 
                end if
          	hb2m(i,itz)= 20.7355300000d0!hbar**2/(2*massnucleon):
		if (w_so) then
          		Vso(i,itz)=cstso*Vpot(i,itz)/(r*Adiff*(1.d0+ff))
		end if

          	if(itz.eq.1)then ! defining the Coulomb potential
	     		if(w_coul) then
             			if(r.le.Radius) then
                			VC(i)=Nprot*echarg/Radius*(3.d0-(r/Radius)**2)*0.5d0
             			else
                			VC(i)=Nprot*echarg/r
            			endif
            		else
                		VC(i)=0.0d0
             		end if
             	Vpot(i,itz)=Vpot(i,itz)+cstcoul*VC(i)
          	endif
	  case("IW")
	     r=i*mesh		
             Vpot(i,itz)= 0.0d0
             hb2m(i,itz)= 20.7355300000d0
	  case("FW") 
	     r=i*mesh
	     if (r.lt.Radius) then		
             Vpot(i,itz)= Vdepth
             hb2m(i,itz)= 20.7355300000d0
             else
             Vpot(i,itz)= 0.0d0
             hb2m(i,itz)= 20.7355300000d0
             end if
	  !case("OH") 
          !   r=i*mesh
          !   Vpot(i,itz)= -r**2*m*w**2 / 2.d0
          case default
             print *, 'No potential available'
             stop
          end select
	
          if (itz==0) write(1000,*)r,Vpot(i,itz),Vso(i,itz)
          if (itz==1) write(1001,*)r,Vpot(i,itz),Vso(i,itz),VC(i)
       enddo

    enddo
    close(1000)
    close(1001)
  end subroutine setWSpotential

  !
  !!
  !

  subroutine HFbasis(lmax,Ncount,Emax)
    implicit none

    integer :: l,lmax,jj,spin,ix
    integer :: Ncount(0:1),iok,N,it

    real (pr) :: u(Npoints),Eout,Emax


    Ncount=0
    do it=0,1
       do l=0,lmax
          do spin=-1,1,2
             jj=2*l+spin
             if(jj.lt.0)go to 333

             do N=1,Nmaxstate


                call numerov(it,N,l,jj,Eout,Emax,U,iok)
                if(iok.eq.1.and.Eout.lt.Emax)then
                   Ncount(it)=Ncount(it)+1
                   EE(Ncount(it),it)=Eout
                   Lused(Ncount(it),it)=l
                   Jused(Ncount(it),it)=jj
                   do ix=1,Npoints
                      Psi(ix,Ncount(it),it)=U(ix)

                     ! write(999+it,*)ix,            psi(ix,Ncount(it),it)
                   enddo

                else
                   go to 333
                endif

             enddo
333          continue
          enddo

       enddo
    enddo

  end subroutine HFbasis

  !
  !!
  !

  subroutine numerov(it,N,l,jj,Eout,Ecut,Uout,iok)
    implicit none

    integer :: N,l,jj,iok,ie,Node
    integer :: ix,Nodecount,it
    integer  :: kstar,flagstar,ii,KL,k4,k3,k
    !
    real (pr) :: uout(Npoints),Eout,Ecut
    real (pr) :: utrial(Npoints),ff(Npoints)
    real (pr) :: vx,vm,vp,a1,a2,a3
    real (pr) :: gg(Npoints),rr
    real (pr) :: ESO,GI(Npoints),GIprimo(Npoints)
    real (pr) :: Eup,Edown,Etrial,CC
    real (pr) :: Fstar, Fstarprec,Wtemp,FAC
    !
    integer, parameter :: imax=2000
    real (pr), parameter :: dtol=10.0_pr**(-10.0_pr)



    Node=N-1
    utrial=0.0_pr
    iok=1 

    ESO=(JJ*(JJ+2.0_pr)-4.0_pr*L*(L+1)-3.0_pr)/8.0_pr
    if(L.eq.0)ESO=0.0_pr
    do ix=1,Npoints
       rr=mesh*ix
       ff(ix)=Vpot(ix,it) & !+dhmen(ix,it)/rr &
            +hb2m(ix,it)*L*(L+1)/rr**2+ESO*VSO(ix,it)!/rr
       gg(ix)=hb2m(ix,it) 
       GI(ix)=DHMEN(ix,it)/Hb2m(ix,it)
       GIprimo(ix)=hb2m(ix,it)*D2HMEN(ix,it)-(DHMEN(ix,it))**2
       GIprimo(ix)= GIprimo(ix)/(Hb2m(ix,it))**2      
    enddo
    Edown=Minval(ff)
    Eup=Ecut+10.0_pr


    !=== prepare spin-orbit 

    do ie=1,imax
       Etrial=0.5_pr*(Eup+Edown)    
       !write(*,*)Etrial

       flagstar=0
       do ix=1,npoints
          Fstar=ff(ix)-Etrial
          !c           if(L.gt.10)then
          if(flagstar.eq.0.and.Fstar.eq.0.d0)then
             kstar=ix
             flagstar=1             
          elseif(ix.ge.2)then
             if(flagstar.eq.0.and.Fstarprec*Fstar.lt.0.d0)then
                kstar=ix
                flagstar=1
             end if
          end if
          !c           end if
          Fstarprec=Fstar
       end do

       if(L.le.1)kstar=1

       utrial(1)=mesh**l

       do ix=1,Npoints-1
          if(ix.eq.1)then
             vp=-(ff(ix+1)-Etrial)/gg(ix+1)-0.5_pr*GI(ix+1)-0.25*GIprimo(ix+1)
             vx=-(ff(ix)-Etrial)/gg(ix)-0.5_pr*GI(ix)-0.25*GIprimo(ix)
             a1=2*(1.0_pr-5.0_pr*mesh**2/12.0_pr*vx)
             a3=(1.0_pr+mesh**2/12.0_pr*vp)
             utrial(ix+1)= (a1*utrial(ix))/a3
          else
             vm=-(ff(ix-1)-Etrial)/gg(ix-1)-0.5_pr*GI(ix-1)-0.25*GIprimo(ix-1)
             vp=-(ff(ix+1)-Etrial)/gg(ix+1)-0.5_pr*GI(ix+1)-0.25*GIprimo(ix+1)
             vx=-(ff(ix)-Etrial)/gg(ix)-0.5_pr*GI(ix)-0.25*GIprimo(ix)
             a1=2*(1.0_pr-5.0_pr*mesh**2/12.0_pr*vx)
             a2=(1.0_pr+mesh**2/12.0_pr*vm)
             a3=(1.0_pr+mesh**2/12.0_pr*vp)
             !write(*,*)a1,a2,a3,vm,vp,vx
             utrial(ix+1)= (a1*utrial(ix)-a2*utrial(ix-1))/a3
          endif
       enddo

       Nodecount=0

       do  ix=kstar,Npoints-1
          if(utrial(ix-1)*utrial(ix).lt.0.0_pr)Nodecount=Nodecount+1        
          !write(90,*)ix,utrial(ix)
       enddo


       if(Nodecount.gt.Node)then
          Eup=Etrial
       elseif(Nodecount.le.Node)then
          Edown=Etrial
       endif

       if(abs(Eup-Edown).lt.dtol)exit

    enddo

    if(Etrial.gt.Ecut)iok=-1

    !======== cleaning the tail

    do ix=npoints,kstar,-1
       if(Utrial(ix-1)*Utrial(ix).lt.0.d0)then
          ii=ix-1
          GO TO 70
       elseif(Utrial(ix-1)*Utrial(ix).eq.0.d0)then
          ii=ix-2
          GO TO 70
       end if
    end do
70  CONTINUE


!!$    if(Etrial.lt.0.d0)then
!!$
!!$       !searching for the position (=imax) of the first maximum/minimum 
!!$       !and computing the average point (=KL) between this max/min and
!!$       !the last node (=ii). This point (=KL) is found in the tail of
!!$       !last oscillation of the wave function, before of the (expected)
!!$       !exponentially decaying behaviour.
!!$       do ix=ii,kstar,-1   
!!$          if(dabs(Utrial(ix-1))-dabs(Utrial(ix)).le.0.d0)then
!!$             KL=(ix+ii)/2
!!$             GO TO 71
!!$          end if
!!$       end do
!!$71     CONTINUE
!!$
!!$       Wtemp=Utrial(KL)           !memorizes the actual value at KL
!!$       Utrial(npoints)=0.d0
!!$       Utrial(npoints-1)=1.E-10
!!$       k4=npoints-KL           !this is the number of points between Kl and the edge of the box
!!$       do k3=npoints-2,KL
!!$          k=k3        !k runs from nmaxt-1 to nmaxt-(nmaxt-KL)+1=KL+1
!!$          !---- propagation of the solution from the right to k-1=KL
!!$             vm=-(ff(k-1)-Etrial)/gg(k-1)-0.5_pr*GI(k-1)-0.25*GIprimo(k-1)
!!$             vp=-(ff(k+1)-Etrial)/gg(k+1)-0.5_pr*GI(k+1)-0.25*GIprimo(k+1)
!!$             vx=-(ff(k)-Etrial)/gg(k)-0.5_pr*GI(k)-0.25*GIprimo(k)
!!$             a1=2*(1.0_pr-5.0_pr*mesh**2/12.0_pr*vx)
!!$             a2=(1.0_pr+mesh**2/12.0_pr*vm)
!!$             a3=(1.0_pr+mesh**2/12.0_pr*vp)
!!$             !write(*,*)a1,a2,a3,vm,vp,vx
!!$             utrial(ix-1)= (a1*utrial(ix)-a3*utrial(ix+1))/a2
!!$       end do
!!$       !---- rescaling the right part in order to match the left part at KL
!!$       FAC=Wtemp/Utrial(KL)
!!$       do k=KL,npoints
!!$          Utrial(k)=Utrial(k)*FAC
!!$       end do
!!$    elseif(Etrial.ge.0.d0)then
!!$       do k=ii,npoints
!!$          Utrial(k)=Utrial(k)-dfloat(k-ii)*Utrial(npoints)/dfloat(Npoints-ii)
!!$       end do
!!$    end if


    !---- computing U(i)=W(i)/sqrt(HME(i))
    do ix=1,npoints
       Utrial(ix)=Utrial(ix)/(dsqrt(Hb2m(ix,it)))
    end do


    CC=0.0_pr!tiny(1.0_pr)

    do ix=1,Npoints
       CC=CC+utrial(ix)**2*mesh
    enddo

    utrial(:)=utrial(:)/sqrt(CC)
    do ix=1,Npoints
       Uout(ix)=utrial(ix)
    enddo
    Eout=Etrial


    return


  end subroutine numerov



end module ws_auxiliary
