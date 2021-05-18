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

# if defined (TRACERS)
    subroutine eval_tracer(ISTEP)
!---------------------------------------------------------------------
!   Update Passive Tracer(s) concentration(s)
!   Last update: Lihan HU, 16/05/2021
!---------------------------------------------------------------------
    use global
    implicit none
    integer, intent(in) :: ISTEP
    real(SP), dimension(:,:,:), allocatable :: R5,DUfs,DVfs,Wfs
    real(SP), dimension(:), allocatable :: Acoef,Bcoef,Ccoef,Xsol,Rhs0
    real(SP) :: SchC
    integer :: i,j,k,IVAR,Nlen

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

    ! advection and diffusion 
    IVAR = 7 
    call adv_scalar_hlpa(DUfs,DVfs,Wfs,Trac,R5,IVAR)

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

          Rhs0(Nlen) = DTrac(i,j,k)+dt*R5(i,j,k)
        enddo
      
        call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

        Nlen = 0
        do k = Kbeg,Kend
          Nlen = Nlen+1
          DTrac(i,j,k) = Xsol(Nlen)
        enddo
      enddo
      enddo

      ! update Tracers
      do k = Kbeg,Kend
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        DTrac(i,j,k) = ALPHA(ISTEP)*DTrac0(i,j,k)+BETA(ISTEP)*DTrac(i,j,k)
        if(Mask(i,j)==0) DTrac(i,j,k) = Zero
      enddo
      enddo
      enddo
    else
      ! update Tracers
      do k = Kbeg,Kend
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        DTrac(i,j,k) = ALPHA(ISTEP)*DTrac0(i,j,k)+BETA(ISTEP)*(DTrac(i,j,k)+dt*R5(i,j,k))
        if(Mask(i,j)==0) DTrac(i,j,k) = Zero
      enddo
      enddo
      enddo
    endif

    ! update Tracers
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      Trac(i,j,k) = DTrac(i,j,k)/D(i,j)
    enddo
    enddo
    enddo

    ! boundary condition and ghost cells
    call tracer_bc
# if defined (PARALLEL)
    call phi_3D_exch(Trac)
# endif          

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


    subroutine tracer_bc
!-------------------------------------------------------------------
!
!   Boundary conditions for Tracers
!
!------------------------------------------------------------------
    use global
    implicit none
    integer :: i,j,k

# if defined (PARALLEL)
    if(n_west.eq.MPI_PROC_NULL) then
# endif
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      do i = 1,Nghost
        if(Bc_X0==1.or.Bc_X0==2.or.Bc_X0==5) then ! added by Cheng for wall friction
          Trac(Ibeg-i,j,k) = Trac(Ibeg+i-1,j,k)
        elseif(Bc_X0==3) then
          Trac(Ibeg-i,j,k) = 2.0*Trin_X0(j,k)-Trac(Ibeg+i-1,j,k)
        elseif(Bc_X0==4) then
          Trac(Ibeg-i,j,k) = Trac(Ibeg+i-1,j,k)
        elseif(Bc_X0==8) then
          Trac(Ibeg-i,j,k) = Trin_X0(j,k)
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
          Trac(Iend+i,j,k) = Trac(Iend-i+1,j,k)
        elseif(Bc_Xn==3) then
          Trac(Iend+i,j,k) = 2.0*Trin_Xn(j,k)-Trac(Iend-i+1,j,k)
        elseif(Bc_Xn==4) then
          Trac(Iend+i,j,k) = Trac(Iend-i+1,j,k)
        elseif(Bc_Xn==8) then
          Trac(Iend+i,j,k) = Trin_Xn(j,k)
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
          Trac(i,Jbeg-j,k) = Trac(i,Jbeg+j-1,k)
        elseif(Bc_Y0==4) then
          Trac(i,Jbeg-j,k) = Trac(i,Jbeg+j-1,k)
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
          Trac(i,Jend+j,k) = Trac(i,Jend-j+1,k)
        elseif(Bc_Yn==4) then
          Trac(i,Jend+j,k) = Trac(i,Jend-j+1,k)
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
        Trac(i,j,Kbeg-k) = Trac(i,j,Kbeg+k-1)
      enddo
      do k = 1,Nghost
        Trac(i,j,Kend+k) = Trac(i,j,Kend-k+1)
      enddo
    enddo
    enddo

    end subroutine tracer_bc
# endif