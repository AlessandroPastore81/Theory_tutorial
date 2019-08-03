module rspace

contains

  subroutine boundary(h,Lmax,Ecut, &
       Nmaxt,PS,NNST,JJP,LP,EP,isospin)

    !========================================================================
    ! subroutine that calculates the single particle levels
    ! +++++++    INPUT   +++++++++++++++++++++++++++++++++++++++++++++
    ! the fields are passed through the module sp_potentials
    !
    ! Lmax: maximum value of angular momentum
    ! Ecut: cut off energy in the single particle spectrum
    ! nmaxt: number of points of the box
    ! nnp: maximum dimension of vectors
    ! isospin: isospin 0 neutrons,1 protons
    ! iwrite_sp: flag to write (1) or not (0) the single particle levels 
    !            on a file
    !
    ! +++++++    OUPUT   ++++++++++++++++++++++++++++++++++++++++++++
    ! PS: vector containing single particle wavefunctions
    ! nnst: number of states found within the cutoff
    ! JJP,LP,EP,LPS : vectors contiaing l,j,energy and a counter 
    !                 for every single particle state
    !
    !========================================================================

    use ws_param ,only : dhmen,d2hmen,hb2m,vpot,vso,Nmaxstate,PI

    implicit none 

    integer S,Nmax,j
    integer Lmin,Lmax,Nmaxt,nnp,i,ii,i1,i2
    Integer LP(Nmaxstate),JJP(Nmaxstate) 
    integer nnst,nh,ij,n,kx,l,no,iwrite_sp
    integer Amass,Nprot,Neutr,isospin,Igram
    integer nh_pre,nh_post,ncount,jj
    PARAMETER (Nmax=50000)

    double precision U(Nmaxt),GI(Nmaxt),GIprimo(Nmaxt)
    double precision ecut,emax0,h,r,e
    double precision PS(Nmaxstate,Nmaxt),EP(Nmaxstate)
    double precision EE(Nmaxstate),temp,ex


    Lmin=0
    igram=1
    emax0=2.d0*ecut+5.d0
    NNP=Nmaxt

    !---- azzero i vettori -----------------------
    U=0.d0
    PS=0.d0
    EP=0.d0     
    EE=0.d0
    !
    LP=0
    JJP=0

    nh=1

    do L=Lmin,Lmax

       do S=1,0,-1
          if(s.eq.1)ij=1
          if(s.eq.0)ij=2

          J=2*L+2*S-1
          if(L.eq.0.and.S.eq.0)GO TO 300
          nh_pre=nh

          do N=1,Nmax,1        

             do kx=1,nmaxt
                GI(kx)=DHMEN(kx,isospin)/Hb2m(kx,isospin)
                GIprimo(kx)=Hb2m(kx,isospin)*D2HMEN(kx,isospin)-(DHMEN(kx,isospin))**2
                GIprimo(kx)= GIprimo(kx)/(Hb2m(kx,isospin))**2                
             end do

             CALL NUMEROV  (h,nmaxt,Emax0,N,L,J,E,U,Vpot(:,isospin),  &
                  VSO(:,isospin),Hb2m(:,isospin),GI,GIprimo,no)

             if(no.eq.-1.and.e.lt.ecut) then
                write(*,*)'Probabile errore per single particles!!'
                write(*,*)e,n,l,j
                !stop
             endif
             if(no.ne.-1.and.e.lt.ecut)then
                EE(n)=e
                if(n.gt.1)then
                   if(EE(n).eq.EE(n-1))then
                      write(*,*)'Errore grossolano'
                   endif
                endif

               ! write(*,*)E,N,L,J


                if(nh.gt.Nmaxstate) then
                   write(*,*)'increase imaxlev!!'
                   stop
                endif
                LP(NH)=L
                EP(Nh)=E
                if(nh.gt.Nmax)stop 'Aumentare Nmax in boundary'
                JJP(Nh)=j !2*l-2*s+1

                do ii=1,nmaxt 
                   PS(NH,ii)=U(ii)
                enddo
                nh=nh+1
             else
                go to 300
             endif

          end do !N           

300       CONTINUE
          Nh_post=nh-1

       end do

    end do

    NNST=NH-1


    return
  END SUBROUTINE boundary




!0000000000000000000000000000000000000000000000000000000000000000000000

SUBROUTINE NUMEROV(h,Nrmax,Emax0,N,L,J,EFINAL,U,V,VSO,&
     HME,GI,GIprimo,no)

  ! use parameters, only:  Nrmax

  implicit none

  integer kt,nd, kstar,flagstar,i,ND1,ND2,ND3
  integer imax,ii,NDfinal,k3,k4,kl,k,L,J,no,niter
  integer Nmaxt,iabswf,iprogress,iout
  integer Node,N,Nrmax

  Double precision MA,POT 
  Double precision U(Nrmax),V(Nrmax),HME(Nrmax),VSO(Nrmax)
  Double precision GI(Nrmax),GIprimo(Nrmax),EFFE,KAPPAQ(Nrmax),&
       W(Nrmax),WAUX(Nrmax)
  Double precision DE,r,E1,WE1,WE3,EFINAL,Wtemp,sum
  Double precision Fstar,Fstarprec,Etrial,Emin,FAC
  Double precision Sig,POTmin,H,ESO,ABSWF,E2,E3
  Double precision EMAX,WE2,EMAX0

  no=0
  Nmaxt=Nrmax
  !-------- Description of this subroutine --------

  !1 - setting Emin equal to the minimum of the potential
  !2 - first search for the eigenfunction. It stops when the number of
  !    node is equal to N-1. (iprogress=1 for this step)
  !3 - second search for the eigenfunction. It stops when the number of
  !    node is equal to N. (iprogress=2 for this step). ETRIAL starts from
  !    the last value in the "iprogress=1 phase".
  !4 - checking for errors on node numbers (ND2-ND1 should be equal to 1)
  !5 - finer search for the eigenstate. The shooting method from here after is governed
  !    looking at the sign of the wave function at the edge of the box.
  !    The search stops when the difference between EMIN and EMAX is less than
  !    a fixed parameter DE, assigned at the beginning of the subroutine.
  !6 - cleaning the tail of the solutions, both for bound and unbound states.
  !    For bound states the exponentially decaying behaviour for r-->infinity
  !    is built by propagation of the solution from the right, and matching with the
  !    left part of the wave function. For unbound states, an addition of a linear function
  !    is performed in order to shift the last node to the edge of the box.
  !7 - going back to the physical wave function, by dividing by a factor
  !    sqrt((hbar**2/2m)*mbare/mstar) and normalizing this wave function.
  !    The outcoming wave function is called U(i). Moreover, U(i) is multiplied
  !    by the sign of U(kstar+1)

  !---- defining the accuracy for the energy (in MeV) ----
  kstar=0
  DE=0.5E-12 ! precision on the energy
  niter=5000
  NODE=N-1  !this is the search number of nodes         
  EMAX=EMAX0!ecut+emax0


  iabswf=-1
  if(nmaxt.gt.Nrmax)then
     write(*,*)'problem in Numerov'
     stop
  endif

  iprogress=1  !this index indicates the progress of the calculation 

  ESO=(J*(J+2.d0)-4.d0*L*(L+1)-3.d0)/8.d0
  if(L.eq.0)ESO=0.d0

56 CONTINUE


  !---- setting Emin equal to the minimum of the potential ----
  POTmin=0.d0
  do i=1,nmaxt
     R=H*I
     POT=(V(i)+ESO*VSO(i))+L*(L+1)*HME(i)/R**2.d0
     if(POT.lt.POTmin)POTmin=POT
  end do
  EMIN=POTmin

59 CONTINUE

  DO KT=1,NITER

     ETRIAL = (EMIN+EMAX)/2.d0
     W(1)=1.E-10
     ND=0
     IOUT=0

     !---- determining the point (kstar) before which the wave function
     !     cannot have an oscillating behaviour.
     flagstar=0
     do i=1,nmaxt
        R=H*I
        POT=(V(i)+ESO*VSO(i))+L*(L+1)*HME(i)/R**2.d0
        Fstar=POT-Etrial
        !c           if(L.gt.10)then
        if(flagstar.eq.0.and.Fstar.eq.0.d0)then
           kstar=i
           flagstar=1             
        elseif(i.ge.2)then
           if(flagstar.eq.0.and.Fstarprec*Fstar.lt.0.d0)then
              kstar=i
              flagstar=1
           end if
        end if
        !c           end if
        Fstarprec=Fstar

        EFFE=(ETRIAL-POT)/HME(i)
        KAPPAQ(i)=-0.5d0*GIprimo(i)-0.25d0*(GI(i))**2.d0+EFFE
     end do

     if(L.le.1)kstar=1

     !---- propagation of the solution till U(nmaxt)       
     DO I = 1,NMAXT-1!NMXT-2
        R=H*I

        !---- to compute W(2) it uses W(1) and W(0)=0.0
        if(i.eq.1)then
           W(i+1)=W(i)*2.d0*(1.d0-(5.d0*H**2.d0/12.d0)*KAPPAQ(i))
           W(i+1)=W(i+1)/(1.d0+(H**2.d0/12.d0)*KAPPAQ(i+1)) 
        else
           W(i+1)=W(i)*2.d0*(1.d0-(5.d0*H**2.d0/12.d0)*KAPPAQ(i)) &
                -W(i-1)*(1.d0+(H**2.d0/12.d0)*KAPPAQ(i-1))
           W(i+1)=W(i+1)/(1.d0+(H**2.d0/12.d0)*KAPPAQ(i+1))
        end if

        !---- counting the nodes beyond kstar, to avoid counting
        !     non-physical nodes
        if(i.ge.kstar)then
           IF(W(I+1)*W(I).LT.0)then
              ND=ND+2
           elseif(W(I+1)*W(I).eq.0)then
              ND=ND+1
           end if
        end if  !it needs to divide ND by 2 (made later)


     END DO

46   CONTINUE
     ND=ND/2


     !---- shooting method (the node at the edge of the box does not count)
     !     It varies the energy (using bisection) until the number of nodes
     !     is equal to:
     !                  N-1 for iprogress=1
     !                  N   for iprogress=2
     if(iprogress.eq.1)then
        if(ND.eq.NODE)then
           GO TO 61    
        elseif(ND.lt.NODE)then
           EMIN=ETRIAL
        elseif(ND.gt.NODE)then
           EMAX=ETRIAL
        end if
     elseif(iprogress.eq.2)then
        if(ND.eq.NODE+1)then
           GO TO 61    
        elseif(ND.lt.NODE+1)then
           EMIN=ETRIAL
        elseif(ND.gt.NODE+1)then
           EMAX=ETRIAL
        end if
     elseif(iprogress.eq.3)then
        E3=ETRIAL
        WE3=W(nmaxt)*W(kstar)
        if(w(kstar).eq.0.d0)then
           write(*,*)'ERROR for N,L,J=',N,L,J
           write(*,*)'w(kstar)=0.0'
           NO=-1 
           GO TO 100
        end if
        ND3=ND
        if(WE1*WE3.gt.0.d0)then     !solution 3 has W(nmaxt) of the same sign as for
           !solution 1. That is, their energies are both HIGHER
           !than the eigenstate.
           ND1=ND3
           WE1=WE3
           E1=E3
           GO TO 36  !when iprogress=3 the shooting method decisions are taken
           !outside the cycle on KT. 

        elseif(WE1*WE3.lt.0.d0)then !solution 3 has W(nmaxt) of the same sign as for
           !solution 2. That is, their energies are both LOWER
           !than the eigenstate.

           ND2=ND3
           WE2=WE3
           E2=E3
           do i=1,nmaxt
              Waux(i)=W(i)  !memorizing solution 2 (ND=N+1)
           end do
           GO TO 36  !when iprogress=3 the shooting method decisions are taken
           !outside the cycle on KT.
        elseif(WE1*WE3.eq.0.d0)then !since WE1 is NOT equal to zero, this means that
           !WE3=0.d0. In other words, solution 3 is the
           !eigenstate we are searching for.
           EFINAL=E3
           NDfinal=ND3
           do i=1,nmaxt
              Waux(i)=W(i)  !memorizing solution 2 (ND=N+1)
           end do
           GO TO 101
        end if
     end if


60   CONTINUE          
  end do  !KT

61 CONTINUE

  if(KT.eq.NITER+1)then  !KT ran 1 through NITER and a solution
     !has not been found. The energy of the
     !searched eigenstate is outside the energy range
     ![EMIN,EMAX]
     !write(*,*)'eccoci',kt,niter+1
     NO=-1 
     GO TO 100
  end if

  if(iprogress.eq.1)then
     ND1=ND            !memorizing number of nodes ND,wave function at the edge, energy
     WE1=W(nmaxt)*W(kstar)
     if(w(kstar).eq.0.d0)then
        write(*,*)'ERROR for N,L,J=',N,L,J
        write(*,*)'w(kstar)=0.0'
        NO=-1 
        GO TO 100
     end if
     E1=ETRIAL

     iprogress=iprogress+1

     GO TO 59
  elseif(iprogress.eq.2)then
     ND2=ND            !memorizing number of nodes ND,wave function at the edge, energy
     WE2=W(nmaxt)*W(kstar)
     if(w(kstar).eq.0.d0)then
        write(*,*)'ERROR for N,L,J=',N,L,J
        write(*,*)'w(kstar)=0.0'
        NO=-1 
        GO TO 100
     end if
     E2=ETRIAL

     do i=1,nmaxt
        Waux(i)=W(i)  !memorizing solution 2
     end do
  end if

  !---- checking for errors
  !     ND2 should be ND1+1
36 if(ND2-ND1.gt.1.or.ND2-ND1.lt.1)then
     NO=-1 
     GO TO 100
  end if

  !---- stepping to iprogress=3
  !     The shooting method from here after is governed
  !     looking at the sign of the wave function at the edge of the box.
  iprogress=3

  if(WE1*WE2.gt.0.d0)then
     NO=-1 
     GO TO 100         
  elseif(WE1*WE2.eq.0.d0.and.WE1.eq.0.d0)then
     NO=-1 
     GO TO 100
  elseif(WE1*WE2.eq.0.d0.and.WE1.ne.0.d0)then  !solution found in iprogress=2
     !is the right solution
     EFINAL=E2 
     NDfinal=ND2
     GO TO 101
  elseif(WE1*WE2.lt.0.d0)then  !we have to keep on playing the shooting method game
     !until the requested accuracy on the energy is required.         
     if(EMAX-EMIN.le.DE)then
        EFINAL=E2                                 
        NDfinal=ND2
        GO TO 101
     else 
        EMIN=E1                  !the solution is in between solution 1 (with ND=N) 
        EMAX=E2                  !and solution 2 (with ND=N+1)
        GO TO 59  !let's propagate the solution in the same way used for first and second phases.
        !A difference occurs: the shooting method is driven by the sign of the
        !wave function at the edge of the box.

     end if
  end if

  !---- last part of the game: the cleaning of the tail of the wave function
101 CONTINUE

  do i=1,nmaxt
     W(i)=Waux(i)  !solution 2 (=the one with ND=ND+1), which is memorized
     !in Waux, is now recorded in W(i)
  end do

  !---- cleaning the tail of wave functions
  ! A   finding the last node of the wave function (coordinate=ii+1 or ii+2)
  do i=nmaxt,kstar,-1
     if(W(i-1)*W(i).lt.0.d0)then
        ii=i-1
        GO TO 70
     elseif(W(i-1)*W(i).eq.0.d0)then
        ii=i-2
        GO TO 70
     end if
  end do
70 CONTINUE

  ! B   for bound states: propagation of the solution from the edge of the box for i>=ii
  !                       and subsequent matching to the left part (i<ii)
  !     for unbound states: addition of a linear function for i>=ii, in order
  !                         to shift the last node to the edge of the box.
  if(EFINAL.lt.0.d0)then

     !searching for the position (=imax) of the first maximum/minimum 
     !and computing the average point (=KL) between this max/min and
     !the last node (=ii). This point (=KL) is found in the tail of
     !last oscillation of the wave function, before of the (expected)
     !exponentially decaying behaviour.
     do imax=ii,kstar,-1   
        if(dabs(W(imax-1))-dabs(W(imax)).le.0.d0)then
           KL=(imax+ii)/2
           GO TO 71
        end if
     end do
71   CONTINUE

     Wtemp=W(KL)           !memorizes the actual value at KL
     W(nmaxt)=0.d0
     W(nmaxt-1)=1.E-10
     k4=nmaxt-KL           !this is the number of points between Kl and the edge of the box
     do k3=2,k4
        k=nmaxt-k3+1        !k runs from nmaxt-1 to nmaxt-(nmaxt-KL)+1=KL+1
        R=H*I
        !---- propagation of the solution from the right to k-1=KL
        W(k-1)=W(k)*2.d0*(1.d0-(5.d0*H**2.d0/12.d0)*KAPPAQ(k)) &
             -W(k+1)*(1.d0+(H**2.d0/12.d0)*KAPPAQ(k+1))
        W(k-1)=W(k-1)/(1.d0+(H**2.d0/12.d0)*KAPPAQ(k-1))
     end do
     !---- rescaling the right part in order to match the left part at KL
     FAC=Wtemp/W(KL)
     do k=KL,nmaxt
        W(k)=W(k)*FAC
     end do
  elseif(EFINAL.ge.0.d0)then
     do k=ii,nmaxt
        W(k)=W(k)-dfloat(k-ii)*W(nmaxt)/dfloat(nmaxt-ii)
     end do
  end if

  !---- computing U(i)=W(i)/sqrt(HME(i))
  do i=1,nmaxt
     U(i)=W(i)/(dsqrt(HME(i)))
  end do

  !---- normalization of the radial wave function
  !     (employing the extended Simpson's rule (eq.4.1.13 Fortran Numerical Recipes))
!!$  SUM=0.0d0
!!$  do i=1,nmaxt-1,2 !odd
!!$     SUM=SUM+(4.d0/3.d0)*U(i)**2.d0
!!$  end do
!!$  do i=2,nmaxt-1,2 !even
!!$     SUM=SUM+(2.d0/3.d0)*U(i)**2.d0
!!$  end do
!!$
!!$  SUM=SUM+(1.d0/3.d0)*U(nmaxt)**2.d0 !final point

 SUM=0.0d0
!!$  do i=1,nmaxt
!!$     SUM=SUM+U(i)**2.d0
!!$  end do
 do i=1,nmaxt-1
     SUM=SUM+U(i)**2.d0
  end do
  SUM=SUM+U(nmaxt)**2.d0/2.d0

  SUM=DSQRT(H*SUM)

  !---- giving a sign to the wave function, multiplying
  !     by the sign of U(kstar+1)
  SIG=U(kstar+1)/dabs(U(kstar+1))
  do i=1,nmaxt
     U(i)=U(i)*SIG/SUM
  end do

  !---- checking for divergent behaviour of the wave function
  !     at the edge of the box
  ABSWF=0.d0
  do i = 1,nmaxt
     if(EFINAL.lt.0.d0.and.dabs(U(i)).gt.ABSWF)then
        IABSWF=i
        ABSWF=dabs(U(i))
     end if
  end do

  IF (IABSWF.eq.nmaxt.and.EFINAL.lt.0.d0)then
     GO TO 100
  end if

  RETURN  !exiting without errors

  !---- exiting with error: use NO outside the subroutine
  !     to check exiting status
100 NO=-1
  RETURN

END SUBROUTINE NUMEROV


end module rspace
