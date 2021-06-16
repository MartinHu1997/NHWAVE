!-----------------------------------------------------------------------------------------------------------
!    This subroutine is used to evaluate Passive Tracers
!    Called by
!       main
!    Last update: May 2021, Lihan HU
!
!    Note that 
!
!    Change log:
!
!----------------------------------------------------------------------------------------------------------

!-----------------------------------------
!---------------- added by HU, 20210519
# if defined (TRACERS)
    subroutine eval_tracer(ISTEP)
!---------------------------------------------------------------------
!   Update Passive Tracer(s) concentration(s)
!   Last update: Lihan HU, 19/05/2021
!---------------------------------------------------------------------
    use global
    implicit none
    integer, intent(in) :: ISTEP
    real(SP), dimension(:,:,:), allocatable :: R5,DUfs,DVfs,Wfs
    real(SP), dimension(:), allocatable :: Acoef,Bcoef,Ccoef,Xsol,Rhs0
    real(SP) :: SchC
    integer :: i,j,k,Nlen,n

    allocate(R5(Mloc,Nloc,Kloc))
    allocate(DUfs(Mloc1,Nloc,Kloc))
    allocate(DVfs(Mloc,Nloc1,Kloc))
    allocate(Wfs(Mloc,Nloc,Kloc1))

    Nlen = Kend-Kbeg+1
    allocate(Acoef(Nlen))
    allocate(Bcoef(Nlen))
    allocate(Ccoef(Nlen))
    allocate(Xsol(Nlen))
    allocate(Rhs0(Nlen))

    DUfs = Ex
    DVfs = Ey
    Wfs = Omega

    do n = 1,N_TRACERS
    ! advection and diffusion 
    call adv_tracer_hlpa(DUfs,DVfs,Wfs,Trac,R5,n)

    if(VISCOUS_FLOW) then
      SchC = 1.0
      do i = Ibeg,Iend
      do j = Jbeg,Jend
        if(Mask(i,j)==0) cycle

        Nlen = 0
        do k = Kbeg,Kend
          Nlen = Nlen+1
          if(k==Kbeg) then
            Acoef(Nlen) = 0.0
          else
            Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
                 0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/SchC)/(0.5*dsig(k)*(dsig(k)+dsig(k-1)))
          endif

          if(k==Kend) then
            Ccoef(Nlen) = 0.0
          else
            Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1))+  &
                 0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/SchC)/(0.5*dsig(k)*(dsig(k)+dsig(k+1)))
          endif
        
          Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)

          Rhs0(Nlen) = DTrac(n,i,j,k)+dt*R5(i,j,k)
        enddo
      
        call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

        Nlen = 0
        do k = Kbeg,Kend
          Nlen = Nlen+1
          DTrac(n,i,j,k) = Xsol(Nlen)
        enddo
      enddo
      enddo

      ! update Tracers
      do k = Kbeg,Kend
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        DTrac(n,i,j,k) = ALPHA(ISTEP)*DTrac0(n,i,j,k)+BETA(ISTEP)*DTrac(n,i,j,k)
        if(Mask(i,j)==0) DTrac(n,i,j,k) = Zero
      enddo
      enddo
      enddo
    else
      ! update Tracers
      do k = Kbeg,Kend
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        DTrac(n,i,j,k) = ALPHA(ISTEP)*DTrac0(n,i,j,k)+BETA(ISTEP)*(DTrac(n,i,j,k)+dt*R5(i,j,k))
        if(Mask(i,j)==0) DTrac(n,i,j,k) = Zero
      enddo
      enddo
      enddo
    endif

    ! update Tracers
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      Trac(n,i,j,k) = DTrac(n,i,j,k)/D(i,j)
!      Trac_t(i,j,k) = Trac(n,i,j,k)
    enddo
    enddo
    enddo

    ! boundary condition and ghost cells
    call tracer_bc(n)

# if defined (PARALLEL)
    call phi_3D_exch(Trac(n,:,:,:))
!    call phi_3D_exch(Trac_t)
!    do k = Kbeg,Kend
!    do j = Jbeg,Jend
!    do i = Ibeg,Iend
!      Trac(n,i,j,k) = Trac_t(i,j,k)
!    enddo
!    enddo
!    enddo
# endif
    enddo

    deallocate(R5)
    deallocate(DUfs)
    deallocate(DVfs)
    deallocate(Wfs)
    deallocate(Acoef)
    deallocate(Bcoef)
    deallocate(Ccoef)
    deallocate(Xsol)
    deallocate(Rhs0)

    return
    end subroutine eval_tracer


    subroutine tracer_bc(N_C)
!-------------------------------------------------------------------
!
!   Boundary conditions for Tracers
!
!------------------------------------------------------------------
    use global
    implicit none
    integer, intent(in) :: N_C
    integer :: i,j,k

# if defined (PARALLEL)
    if(n_west.eq.MPI_PROC_NULL) then
# endif
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      do i = 1,Nghost
        if(Bc_X0==1.or.Bc_X0==2.or.Bc_X0==5) then ! added by Cheng for wall friction
          Trac(N_C,Ibeg-i,j,k) = Trac(N_C,Ibeg+i-1,j,k)
        elseif(Bc_X0==3) then
          Trac(N_C,Ibeg-i,j,k) = 2.0*Trin_X0(N_C,j,k)-Trac(N_C,Ibeg+i-1,j,k)
        elseif(Bc_X0==4) then
          Trac(N_C,Ibeg-i,j,k) = Trac(N_C,Ibeg+i-1,j,k)
        elseif(Bc_X0==8) then
          Trac(N_C,Ibeg-i,j,k) = Trin_X0(N_C,j,k)
        endif
      enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if(n_east.eq.MPI_PROC_NULL) then
# endif
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      do i = 1,Nghost
        if(Bc_Xn==1.or.Bc_Xn==2.or.Bc_Xn==5) then ! added by Cheng for wall friction
          Trac(N_C,Iend+i,j,k) = Trac(N_C,Iend-i+1,j,k)
        elseif(Bc_Xn==3) then
          Trac(N_C,Iend+i,j,k) = 2.0*Trin_Xn(N_C,j,k)-Trac(N_C,Iend-i+1,j,k)
        elseif(Bc_Xn==4) then
          Trac(N_C,Iend+i,j,k) = Trac(N_C,Iend-i+1,j,k)
        elseif(Bc_Xn==8) then
          Trac(N_C,Iend+i,j,k) = Trin_Xn(N_C,j,k)
        endif
      enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if(n_suth.eq.MPI_PROC_NULL) then
# endif
    do i = Ibeg,Iend
    do k = Kbeg,Kend
      do j = 1,Nghost
        if(Bc_Y0==1.or.Bc_Y0==2.or.Bc_Y0==5) then ! added by Cheng for wall friction
          Trac(N_C,i,Jbeg-j,k) = Trac(N_C,i,Jbeg+j-1,k)
        elseif(Bc_Y0==4) then
          Trac(N_C,i,Jbeg-j,k) = Trac(N_C,i,Jbeg+j-1,k)
        endif
      enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif 

# if defined (PARALLEL)
    if(n_nrth.eq.MPI_PROC_NULL) then
# endif
    do i = Ibeg,Iend
    do k = Kbeg,Kend
      do j = 1,Nghost
        if(Bc_Yn==1.or.Bc_Yn==2.or.Bc_Yn==5) then ! added by Cheng for wall friction
          Trac(N_C,i,Jend+j,k) = Trac(N_C,i,Jend-j+1,k)
        elseif(Bc_Yn==4) then
          Trac(N_C,i,Jend+j,k) = Trac(N_C,i,Jend-j+1,k)
        endif
      enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif

    do i = Ibeg,Iend
    do j = Jbeg,Jend
      do k = 1,Nghost
        Trac(N_C,i,j,Kbeg-k) = Trac(N_C,i,j,Kbeg+k-1)
      enddo
      do k = 1,Nghost
        Trac(N_C,i,j,Kend+k) = Trac(N_C,i,j,Kend-k+1)
      enddo
    enddo
    enddo

    end subroutine tracer_bc

    subroutine adv_tracer_hlpa(Flx,Fly,Flz,Phi,R5,N_C)
!--------------------------------------------------------
!   Subroutine for scalar convection and horizontal diffusion  
!   Last update: Gangfeng Ma, 04/04/2012
!-------------------------------------------------------  
    use global
    implicit none
    integer, intent(in) :: N_C
    real(SP), dimension(N_TRACERS,Mloc,Nloc,Kloc),  intent(in) :: Phi
    real(SP), dimension(Mloc1,Nloc,Kloc),  intent(in) :: Flx
    real(SP), dimension(Mloc,Nloc1,Kloc),  intent(in) :: Fly
    real(SP), dimension(Mloc,Nloc,Kloc1),  intent(in) :: Flz
    real(SP), dimension(Mloc,Nloc,Kloc), intent(inout) :: R5
    real(SP), dimension(:,:,:), allocatable :: Scalx,Scaly,Scalz,Sdiffx,Sdiffy
    real(SP) :: DUfs,DVfs,Wfs,Fww,Fw,Fp,Fe,hlpa,SchtH
    real(SP) :: L1top,L2top,L1bot,L2bot,alpha_c,beta_c,gamma_c, &
                dsigck,dsigck1,nuH_top,nuH_bot
    integer :: i,j,k

    allocate(Scalx(Mloc1,Nloc,Kloc))
    allocate(Scaly(Mloc,Nloc1,Kloc))
    allocate(Scalz(Mloc,Nloc,Kloc1))
    allocate(Sdiffx(MLoc,Nloc,Kloc))
    allocate(Sdiffy(Mloc,Nloc,Kloc))

    ! advection in x direction
    Scalx = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend+1
      DUfs = Flx(i,j,k)
      Fww = Phi(N_C,i-2,j,k)
      Fw  = Phi(N_C,i-1,j,k)
      Fp  = Phi(N_C,i,j,k)
      Fe  = Phi(N_C,i+1,j,k)
      Scalx(i,j,k) = DUfs*hlpa(DUfs,Fww,Fw,Fp,Fe)
    enddo
    enddo
    enddo

    ! advection in y direction
    Scaly = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend+1
    do i = Ibeg,Iend      
      DVfs = Fly(i,j,k)
      Fww = Phi(N_C,i,j-2,k)
      Fw  = Phi(N_C,i,j-1,k)
      Fp  = Phi(N_C,i,j,k)
      Fe  = Phi(N_C,i,j+1,k)
      Scaly(i,j,k) = DVfs*hlpa(DVfs,Fww,Fw,Fp,Fe)
    enddo
    enddo
    enddo

    ! advection in z direction
    Scalz = Zero
    do k = Kbeg+1,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      Wfs = Flz(i,j,k)
      Fww = Phi(N_C,i,j,k-2)
      Fw  = Phi(N_C,i,j,k-1)
      Fp  = Phi(N_C,i,j,k)
      Fe  = Phi(N_C,i,j,k+1)
      Scalz(i,j,k) = Wfs*hlpa(Wfs,Fww,Fw,Fp,Fe)
    enddo
    enddo
    enddo

    ! at boundaries
    call flux_tracer_bc(N_C,Scalx,Scaly,Scalz)

    ! Schmidt number
    SchtH = 1.0


    ! estimate horizontal diffusion
    Sdiffx = Zero; Sdiffy = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==0) cycle
      Sdiffx(i,j,k) = 0.5*((0.5*(Cmu(i+1,j,k)+Cmu(i,j,k))+  &
           0.5*(CmuHt(i+1,j,k)+CmuHt(i,j,k))/SchtH)*  &
           (D(i+1,j)+D(i,j))*(Phi(N_C,i+1,j,k)-Phi(N_C,i,j,k))-  &
           (0.5*(Cmu(i,j,k)+Cmu(i-1,j,k))+  &
           0.5*(CmuHt(i,j,k)+CmuHt(i-1,j,k))/SchtH)*  &
           (D(i,j)+D(i-1,j))*(Phi(N_C,i,j,k)-Phi(N_C,i-1,j,k)))/dx**2
      Sdiffy(i,j,k) = 0.5*((0.5*(Cmu(i,j+1,k)+Cmu(i,j,k))+  &
           0.5*(CmuHt(i,j+1,k)+CmuHt(i,j,k))/SchtH)*  &
           (D(i,j+1)+D(i,j))*(Phi(N_C,i,j+1,k)-Phi(N_C,i,j,k))-  &
           (0.5*(Cmu(i,j,k)+Cmu(i,j-1,k))+  &
           0.5*(CmuHt(i,j,k)+CmuHt(i,j-1,k))/SchtH)*  &
           (D(i,j)+D(i,j-1))*(Phi(N_C,i,j,k)-Phi(N_C,i,j-1,k)))/dy**2
    enddo
    enddo
    enddo

!    Sdiffx = Zero
!    Sdiffy = Zero
!    do k = Kbeg,Kend
!       !these are used for vertical gradient in non-uniform grid
!       !see Derakhti et.al 2015a Appendix B
!       dsigck  = (dsig(k)+dsig(k+1))/2.0 
!       dsigck1 = (dsig(k-1)+dsig(k))/2.0 
!       alpha_c = -dsigck/(dsigck+dsigck1)/dsigck1
!       beta_c = (dsigck-dsigck1)/(dsigck*dsigck1)
!       gamma_c = dsigck1/(dsigck+dsigck1)/dsigck
!    do j = Jbeg,Jend
!    do i = Ibeg,Iend
!      if(Mask(i,j)==0) cycle
!       Sdiffx(i,j,k) = (  (0.5*(Cmu(i+1,j,k)+Cmu(i,j,k))+0.5*(CmuHt(i+1,j,k)+CmuHt(i,j,k))/SchtH)! &
!                         *(Phi(N_C,i+1,j,k)*D(i+1,j)-Phi(N_C,i,j,k)*D(i,j))/dx &
!                        -(0.5*(Cmu(i,j,k)+Cmu(i-1,j,k))+0.5*(CmuHt(i,j,k)+CmuHt(i-1,j,k))/SchtH) !&
!                         *(Phi(N_C,i,j,k)*D(i,j)-Phi(N_C,i-1,j,k)*D(i-1,j))/dx &
!                       ) /dx & 
!                     + ( (Cmu(i+1,j,k)+CmuHt(i+1,j,k)/SchtH) &
!                         *(alpha_c*D(i+1,j)*Phi(N_C,i+1,j,k-1)*DelxSc(i+1,j,k-1)&
!                          +beta_c *D(i+1,j)*Phi(N_C,i+1,j,k  )*DelxSc(i+1,j,k  )&
!                          +gamma_c*D(i+1,j)*Phi(N_C,i+1,j,k+1)*DelxSc(i+1,j,k+1)) & 
!                        -(Cmu(i-1,j,k)+CmuHt(i-1,j,k)/SchtH) &
!                         *(alpha_c*D(i-1,j)*Phi(N_C,i-1,j,k-1)*DelxSc(i-1,j,k-1)&
!                          +beta_c *D(i-1,j)*Phi(N_C,i-1,j,k  )*DelxSc(i-1,j,k  )&
!                          +gamma_c*D(i-1,j)*Phi(N_C,i-1,j,k+1)*DelxSc(i-1,j,k+1)) & 
!                       ) /2.0/dx & 
!                     + (alpha_c*(Cmu(i,j,k-1)+CmuHt(i,j,k-1)/SchtH)*DelxSc(i,j,k-1)&
!                         *(D(i+1,j)*Phi(N_C,i+1,j,k-1)-D(i-1,j)*Phi(N_C,i-1,j,k-1))/2.0/dx &
!                        +beta_c*(Cmu(i,j,k)+CmuHt(i,j,k)/SchtH)*DelxSc(i,j,k)&
!                         *(D(i+1,j)*Phi(N_C,i+1,j,k  )-D(i-1,j)*Phi(N_C,i-1,j,k  ))/2.0/dx &
!                        +gamma_c*(Cmu(i,j,k+1)+CmuHt(i,j,k+1)/SchtH)*DelxSc(i,j,k+1)&
!                         *(D(i+1,j)*Phi(N_C,i+1,j,k+1)-D(i-1,j)*Phi(N_C,i-1,j,k+1))/2.0/dx &
!                        ) 
!
!      Sdiffy(i,j,k) = (  (0.5*(Cmu(i,j+1,k)+Cmu(i,j,k))+0.5*(CmuHt(i,j+1,k)+CmuHt(i,j,k))/SchtH) !&
!                         *(Phi(N_C,i,j+1,k)*D(i,j+1)-Phi(N_C,i,j,k)*D(i,j))/dy &
!                        -(0.5*(Cmu(i,j,k)+Cmu(i,j-1,k))+0.5*(CmuHt(i,j,k)+CmuHt(i,j-1,k))/SchtH) !&
!                         *(Phi(N_C,i,j,k)*D(i,j)-Phi(N_C,i,j-1,k)*D(i,j-1))/dy &
!                       ) /dy & 
!                     + ( (Cmu(i,j+1,k)+CmuHt(i,j+1,k)/SchtH) &
!                         *(alpha_c*D(i,j+1)*Phi(N_C,i,j+1,k-1)*DelySc(i,j+1,k-1)&
!                          +beta_c *D(i,j+1)*Phi(N_C,i,j+1,k  )*DelySc(i,j+1,k  )&
!                          +gamma_c*D(i,j+1)*Phi(N_C,i,j+1,k+1)*DelySc(i,j+1,k+1)) & 
!                        -(Cmu(i,j-1,k)+CmuHt(i,j-1,k)/SchtH) &
!                         *(alpha_c*D(i,j-1)*Phi(N_C,i,j-1,k-1)*DelySc(i,j-1,k-1)&
!                          +beta_c *D(i,j-1)*Phi(N_C,i,j-1,k  )*DelySc(i,j-1,k  )&
!                          +gamma_c*D(i,j-1)*Phi(N_C,i,j-1,k+1)*DelySc(i,j-1,k+1)) & 
!                       ) /2.0/dy & 
!                     + (alpha_c*(Cmu(i,j,k-1)+CmuHt(i,j,k-1)/SchtH)*DelySc(i,j,k-1)&
!                         *(D(i,j+1)*Phi(N_C,i,j+1,k-1)-D(i,j-1)*Phi(N_C,i,j-1,k-1))/2.0/dy &
!                        +beta_c*(Cmu(i,j,k)+CmuHt(i,j,k)/SchtH)*DelySc(i,j,k)&
!                         *(D(i,j+1)*Phi(N_C,i,j+1,k  )-D(i,j-1)*Phi(N_C,i,j-1,k  ))/2.0/dy &
!                        +gamma_c*(Cmu(i,j,k+1)+CmuHt(i,j,k+1)/SchtH)*DelySc(i,j,k+1)&
!                         *(D(i,j+1)*Phi(N_C,i,j+1,k+1)-D(i,j-1)*Phi(N_C,i,j-1,k+1))/2.0/dy &
!                        )  
!    enddo
!    enddo
!    enddo

    R5 = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==0) cycle
      R5(i,j,k) = -1.0/dx*(Scalx(i+1,j,k)-Scalx(i,j,k))  &
                  -1.0/dy*(Scaly(i,j+1,k)-Scaly(i,j,k))  &
                  -1.0/dsig(k)*(Scalz(i,j,k+1)-Scalz(i,j,k)) &
                  +Sdiffx(i,j,k)+Sdiffy(i,j,k)
    enddo
    enddo
    enddo

    deallocate(Scalx)
    deallocate(Scaly)
    deallocate(Scalz)
    deallocate(Sdiffx)
    deallocate(Sdiffy)

    end subroutine adv_tracer_hlpa


    subroutine flux_tracer_bc(N_C,Scalx,Scaly,Scalz)
!--------------------------------------------------------
!   Specify boundary conditions for scalar convection
!   Last update: Gangfeng Ma, 09/02/2011
!-------------------------------------------------------
    use global
    implicit none
    integer, intent(in) :: N_C
    real(SP), dimension(Mloc1,Nloc,Kloc), intent(inout) :: Scalx
    real(SP), dimension(Mloc,Nloc1,Kloc), intent(inout) :: Scaly
    real(SP), dimension(Mloc,Nloc,Kloc1), intent(inout) :: Scalz
    real(SP), dimension(Nloc,Kloc) :: Scal_X0,Scal_Xn
    integer :: i,j,k

    ! temporarily set it here
	! added by Cheng for initialization
	Scal_X0 = Zero
	Scal_Xn = Zero
      Scal_X0 = Trin_X0(N_C,:,:)
      Scal_Xn = Trin_Xn(N_C,:,:)
!    do k = Kbeg,Kend
!    do j = Jbeg,Jend
!      Scal_X0(j,k) = Trin_X0(N_C,j,k)
!      Scal_Xn(j,k) = Trin_Xn(N_C,j,k)
!    enddo
!    enddo

    ! left and right side
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
! added by Cheng for nesting. Please search for others with (COUPLING) in this subroutine
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(Bc_X0==1.or.Bc_X0==2.or.Bc_X0==5) then ! added by Cheng for wall friction
         Scalx(Ibeg,j,k) = Zero
       elseif(Bc_X0==3) then
         Scalx(Ibeg,j,k) = Ex(Ibeg,j,k)*Scal_X0(j,k)
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
!       if(Bc_Xn==1.or.Bc_Xn==2.or.Bc_X0==5) then ! added by Cheng for wall friction
         Scalx(Iend1,j,k) = Zero
!       elseif(Bc_Xn==3) then
!         Scalx(Iend1,j,k) = Din_Xn(j)*Uin_Xn(j,k)*Scal_Xn(j,k)
!       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif      

     ! front and back side
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       if(Bc_Y0==1.or.Bc_Y0==2.or.Bc_Y0==5) then ! added by Cheng for wall friction
         Scaly(i,Jbeg,k) = Zero
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif


# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       if(Bc_Yn==1.or.Bc_Yn==2.or.Bc_Yn==5) then ! added by Cheng for wall friction
         Scaly(i,Jend1,k) = Zero
       endif
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==0) then
        Scalx(i,j,k) = Zero
        Scalx(i+1,j,k) = Zero
        Scaly(i,j,k) = Zero
        Scaly(i,j+1,k) = Zero
      endif
    enddo
    enddo
    enddo

    do j = Jbeg,Jend
    do i = Ibeg,Iend
      Scalz(i,j,Kbeg) = Zero
      Scalz(i,j,Kend1) = Zero
    enddo
    enddo

    return
    end subroutine flux_tracer_bc
# endif
!-----------------------------------------
!-----------------------------------------




!-------------------------------------------
!---------------- added by HU, 20210517
!# if defined (TRACERS)
!    subroutine eval_tracer_old(ISTEP)
!!---------------------------------------------------------------------
!!   Update Passive Tracer(s) concentration(s)
!!   Last update: Lihan HU, 16/05/2021
!!---------------------------------------------------------------------
!    use global
!    implicit none
!    integer, intent(in) :: ISTEP
!    real(SP), dimension(:,:,:), allocatable :: R5,DUfs,DVfs,Wfs
!    real(SP), dimension(:), allocatable :: Acoef,Bcoef,Ccoef,Xsol,Rhs0
!    real(SP) :: SchC
!    integer :: i,j,k,IVAR,Nlen
!
!   allocate(R5(Mloc,Nloc,Kloc))
!    allocate(DUfs(Mloc1,Nloc,Kloc))
!    allocate(DVfs(Mloc,Nloc1,Kloc))
!    allocate(Wfs(Mloc,Nloc,Kloc1))
!
!    Nlen = Kend-Kbeg+1
!    allocate(Acoef(Nlen))
!    allocate(Bcoef(Nlen))
!    allocate(Ccoef(Nlen))
!    allocate(Xsol(Nlen))
!    allocate(Rhs0(Nlen))
!
!    DUfs = Ex
!    DVfs = Ey
!    Wfs = Omega
!
!    ! advection and diffusion 
!    IVAR = 7 
!    call adv_scalar_hlpa(DUfs,DVfs,Wfs,Trac,R5,IVAR)
!
!    if(VISCOUS_FLOW) then
!      SchC = 1.0
!      do i = Ibeg,Iend
!      do j = Jbeg,Jend
!        if(Mask(i,j)==0) cycle
!
!        Nlen = 0
!        do k = Kbeg,Kend
!          Nlen = Nlen+1
!          if(k==Kbeg) then
!            Acoef(Nlen) = 0.0
!          else
!            Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
!                 0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/SchC)/(0.5*dsig(k)*(dsig(k)+dsig(k-1)))
!          endif
!
!          if(k==Kend) then
!            Ccoef(Nlen) = 0.0
!          else
!            Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1))+  &
!                 0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/SchC)/(0.5*dsig(k)*(dsig(k)+dsig(k+1)))
!          endif
!        
!          Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)
!
!          Rhs0(Nlen) = DTrac(i,j,k)+dt*R5(i,j,k)
!        enddo
!      
!        call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)
!
!        Nlen = 0
!        do k = Kbeg,Kend
!          Nlen = Nlen+1
!          DTrac(i,j,k) = Xsol(Nlen)
!        enddo
!      enddo
!      enddo
!
!      ! update Tracers
!      do k = Kbeg,Kend
!      do j = Jbeg,Jend
!      do i = Ibeg,Iend
!        DTrac(i,j,k) = ALPHA(ISTEP)*DTrac0(i,j,k)+BETA(ISTEP)*DTrac(i,j,k)
!        if(Mask(i,j)==0) DTrac(i,j,k) = Zero
!      enddo
!      enddo
!      enddo
!    else
!      ! update Tracers
!      do k = Kbeg,Kend
!      do j = Jbeg,Jend
!      do i = Ibeg,Iend
!        DTrac(i,j,k) = ALPHA(ISTEP)*DTrac0(i,j,k)+BETA(ISTEP)*(DTrac(i,j,k)+dt*R5(i,j,k))
!        if(Mask(i,j)==0) DTrac(i,j,k) = Zero
!      enddo
!      enddo
!      enddo
!    endif
!
!    ! update Tracers
!    do k = Kbeg,Kend
!    do j = Jbeg,Jend
!    do i = Ibeg,Iend
!      Trac(i,j,k) = DTrac(i,j,k)/D(i,j)
!    enddo
!    enddo
!    enddo
!
!    ! boundary condition and ghost cells
!    call tracer_bc
!# if defined (PARALLEL)
!    call phi_3D_exch(Trac)
!# endif          
!
!    deallocate(R5)
!    deallocate(DUfs)
!    deallocate(DVfs)
!    deallocate(Wfs)
!    deallocate(Acoef)
!    deallocate(Bcoef)
!    deallocate(Ccoef)
!    deallocate(Xsol)
!    deallocate(Rhs0)
!
!    return
!    end subroutine eval_tracer_old
!
!
!    subroutine tracer_bc_old
!!-------------------------------------------------------------------
!!
!!   Boundary conditions for Tracers
!!
!!------------------------------------------------------------------
!    use global
!    implicit none
!    integer :: i,j,k
!
!# if defined (PARALLEL)
!    if(n_west.eq.MPI_PROC_NULL) then
!# endif
!    do j = Jbeg,Jend
!    do k = Kbeg,Kend
!      do i = 1,Nghost
!        if(Bc_X0==1.or.Bc_X0==2.or.Bc_X0==5) then ! added by Cheng for wall friction
!          Trac(Ibeg-i,j,k) = Trac(Ibeg+i-1,j,k)
!        elseif(Bc_X0==3) then
!          Trac(Ibeg-i,j,k) = 2.0*Trin_X0(j,k)-Trac(Ibeg+i-1,j,k)
!        elseif(Bc_X0==4) then
!          Trac(Ibeg-i,j,k) = Trac(Ibeg+i-1,j,k)
!        elseif(Bc_X0==8) then
!          Trac(Ibeg-i,j,k) = Trin_X0(j,k)
!        endif
!      enddo
!    enddo
!    enddo
!# if defined (PARALLEL)
!    endif
!# endif
!
!# if defined (PARALLEL)
!    if(n_east.eq.MPI_PROC_NULL) then
!# endif
!    do j = Jbeg,Jend
!    do k = Kbeg,Kend
!      do i = 1,Nghost
!        if(Bc_Xn==1.or.Bc_Xn==2.or.Bc_Xn==5) then ! added by Cheng for wall friction
!          Trac(Iend+i,j,k) = Trac(Iend-i+1,j,k)
!        elseif(Bc_Xn==3) then
!          Trac(Iend+i,j,k) = 2.0*Trin_Xn(j,k)-Trac(Iend-i+1,j,k)
!        elseif(Bc_Xn==4) then
!          Trac(Iend+i,j,k) = Trac(Iend-i+1,j,k)
!        elseif(Bc_Xn==8) then
!          Trac(Iend+i,j,k) = Trin_Xn(j,k)
!        endif
!      enddo
!    enddo
!    enddo
!# if defined (PARALLEL)
!    endif
!# endif
!
!# if defined (PARALLEL)
!    if(n_suth.eq.MPI_PROC_NULL) then
!# endif
!    do i = Ibeg,Iend
!    do k = Kbeg,Kend
!      do j = 1,Nghost
!        if(Bc_Y0==1.or.Bc_Y0==2.or.Bc_Y0==5) then ! added by Cheng for wall friction
!          Trac(i,Jbeg-j,k) = Trac(i,Jbeg+j-1,k)
!        elseif(Bc_Y0==4) then
!          Trac(i,Jbeg-j,k) = Trac(i,Jbeg+j-1,k)
!        endif
!      enddo
!    enddo
!    enddo
!# if defined (PARALLEL)
!    endif
!# endif 
!
!# if defined (PARALLEL)
!    if(n_nrth.eq.MPI_PROC_NULL) then
!# endif
!    do i = Ibeg,Iend
!    do k = Kbeg,Kend
!      do j = 1,Nghost
!        if(Bc_Yn==1.or.Bc_Yn==2.or.Bc_Yn==5) then ! added by Cheng for wall friction
!          Trac(i,Jend+j,k) = Trac(i,Jend-j+1,k)
!        elseif(Bc_Yn==4) then
!          Trac(i,Jend+j,k) = Trac(i,Jend-j+1,k)
!        endif
!      enddo
!    enddo
!    enddo
!# if defined (PARALLEL)
!    endif
!# endif
!
!    do i = Ibeg,Iend
!    do j = Jbeg,Jend
!      do k = 1,Nghost
!        Trac(i,j,Kbeg-k) = Trac(i,j,Kbeg+k-1)
!      enddo
!      do k = 1,Nghost
!        Trac(i,j,Kend+k) = Trac(i,j,Kend-k+1)
!      enddo
!    enddo
!    enddo
!
!    end subroutine tracer_bc_old
!# endif
!-----------------------------------------
!-----------------------------------------