!------------------------------------------------------------------------------------------------
!
!   pressure_PDI.F
!
!   Subroutines related to Poisson Solver with PDI Method
!
!    Lihan Hu, 08/08/2022
!
!----------------------------------------------------------------------------------------------
!
!    Subroutines included in this file:
!-----------------------------------------------------------------------------------------
!
!    (1) poisson_solver_pdi
!    (2) hypre_pres_solver_pdi
!    (3) generate_coef_rhs_pdi
!    (4) interpolate_pressure_pdi
!
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!
!    (1) Subroutine poisson_solver_pdi
!
!    Solve poisson equation for dynamic pressure (PDI Method)
!
!    Called by:  eval_duvw
!
!    Last update: 08/08/2022, Lihan Hu
!
!----------------------------------------------------------------------------------------
!
# if defined (PRESSUREACC)
     subroutine poisson_solver_pdi
!
     use global
     implicit none
     integer :: i,j,k,imask,kk
# if !defined (PARALLEL)
     ! variables for serial computation
     real(SP), dimension(:), allocatable :: Wksp
     integer,  dimension(:), allocatable :: IWksp
     real(SP), dimension(:), allocatable :: Phi
     real(SP) :: RPARM(30),Pbar(1)
     integer :: IPARM(30),S(1),IS(1),nwksp,inwksp,Ndim,Mdim,N,Maxnz,ierr,neq,Nlen
     external :: MIC3,IC3,SOR3,GMRES,CG,BCGS
# endif

     ! generate coefficient matrix and rhs
     call generate_coef_rhs_pdi
# if defined (PARALLEL)
     ! use HYPRE package for parallel computation
     call hypre_pres_solver_pdi
# else
     ! use NSPCG package for serial computation
     call dfault(IPARM,RPARM)

     ! reset default values
     IPARM(2) = itmax
     IPARM(3) = 3
     IPARM(4) = 33
     RPARM(1) = tol

      ! PDI_INT means the numerbers of interval subgrids of pressure grids
!      neqns2 = neqns / PDI_INT
      
     Ndim = 5*neqns2
     Mdim = 5*15
     N = neqns2
     Maxnz = 15
     nwksp = 30*neqns2
     inwksp = 10*neqns2

     allocate(Wksp(nwksp))
     allocate(Iwksp(inwksp))
      allocate(Phi(neqns2))
   
      !initial guess
      neq = 0
      do k = Kbeg,Kend,PDI_INT
      do j = Jbeg,Jend
      do i = Ibeg,Iend
          neq = neq+1
          Phi(neq) = P(i,j,k)
      enddo
      enddo
      enddo
!
!    solve Poisson equation
!
!        isolver = 1:
!        isolver = 2:
!        isolver = 3:
!  using preconditioned CG or GMRES
!
     if(isolver==1) then
       call nspcg(MIC3,CG,Ndim,Mdim,N,Maxnz,Coef2,JCoef2,S,IS,  &
           Phi,Pbar,Rhs2,Wksp,IWksp,nwksp,inwksp,IPARM,RPARM,ierr)
     elseif(isolver==2) then
       call nspcg(IC3,GMRES,Ndim,Mdim,N,Maxnz,Coef2,JCoef2,S,IS,  &
           Phi,Pbar,Rhs2,Wksp,IWksp,nwksp,inwksp,IPARM,RPARM,ierr)  
     elseif(isolver==3) then
       call nspcg(SOR3,GMRES,Ndim,Mdim,N,Maxnz,Coef2,JCoef2,S,IS,  &
           Phi,Pbar,Rhs2,Wksp,IWksp,nwksp,inwksp,IPARM,RPARM,ierr)
     endif   

!
!    Interpolation for dynamic pressure in original (velocity) grid
!
!    update the values in pressure grid
     neq = 0
     do k = Kbeg,Kend,PDI_INT
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       neq = neq+1
       P(i,j,k) = Phi(neq) 
     enddo
     enddo
     enddo
!    ==== interpolation ====
      if(PDI_INT.ne.1) call interpolate_pressure_pdi
      
     deallocate(Wksp)
     deallocate(Iwksp)
      deallocate(Phi)

# endif
!
!   fyshi gave boundary condition for dry cells
!   set zero for dry set is inaccurate
!   dry cells  (so what was actually done??)
!
!     do k = Kbeg,Kend
!     do j = Jbeg,Jend
!     do i = Ibeg,Iend
!       if(Mask(i,j)==0) then
!         P(i,j,k) = Zero
!         
!         ! south boundary 
!         if(Mask(i,j+1)==1)then
!           do imask=1,Nghost
!             P(i,j-imask+1,k)=P(i,j+imask,k)
!           enddo
!         ! north boundary
!         elseif(Mask(i,j-1)==1)then
!           do imask=1,Nghost
!             P(i,j+imask-1,k)=P(i,j-imask,k)
!           enddo
!         ! west boundary
!         elseif(Mask(i+1,j)==1)then
!           do imask=1,Nghost
!             P(i-imask+1,j,k)=P(i+imask,j,k)
!           enddo
!         ! east boundary
!         elseif(Mask(i-1,j)==1)then
!           do imask=1,Nghost
!             P(i+imask-1,j,k)=P(i-imask,j,k)
!           enddo
!         endif
!       endif 
!     enddo
!     enddo
!     enddo

!   collect into ghost cells
!
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       do i = 1,Nghost
         P(Ibeg-i,j,k) = P(Ibeg+i-1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do k = Kbeg,Kend
     do j = Jbeg,Jend 
       do i = 1,Nghost     
         P(Iend+i,j,k) = P(Iend-i+1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       do j = 1,Nghost
         P(i,Jbeg-j,k) = P(i,Jbeg+j-1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       do j = 1,Nghost
         P(i,Jend+j,k) = P(i,Jend-j+1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     call phi_3D_exch(P)
# endif

     end subroutine poisson_solver_pdi
# endif

!
!-----------------------------------------------------------------------------------------------
!
!    (2) Subroutine hypre_pres_solver_pdi
!
!    solve for dynamic pressure using HYPRE package
!
!    called by:  poisson_solver
!
!    Last update: 08/08/2022, Lihan Hu
!
!----------------------------------------------------------------------------------------------
!
# if defined (PRESSUREACC)
     subroutine hypre_pres_solver_pdi
!
# if defined (PARALLEL)
     use global
     implicit none
     integer, parameter :: ndim=3
     integer, parameter :: nentries=15
     integer :: i,j,k,n,ivalues,nvalues,neq,ientry,num_iterations,  &
                precond_id,n_pre,n_post,ierr
     integer*8 :: grid,stencil,matrix,vec_b,vec_x,solver,precond
     integer :: i_glob(Mloc),j_glob(Nloc),k_glob(Kloc)
!     integer, dimension(:) :: k_glob
     integer :: ilower(ndim),iupper(ndim),offsets(nentries,ndim),stencil_indices(nentries), &
                periodic_shift(ndim)
     real(SP) :: final_res_norm
     real(SP), dimension(:), allocatable :: values,Phi
     integer, dimension(:,:,:), allocatable :: indx 
     data ((offsets(i,j),j=1,ndim),i=1,nentries)/0,0,0,1,0,0,0,1,0,0,-1,1,-1,0,1,  &
             0,0,1,1,0,1,0,1,1,-1,0,0,0,-1,0,  &
             0,1,-1,1,0,-1,0,0,-1,-1,0,-1,0,-1,-1/
!
!    set up a three dimensional grid
!
     call HYPRE_StructGridCreate(MPI_COMM_WORLD,ndim,grid,ierr)
!
!    global indices
!
!      allocate(k_glob(Kglob/PDI_INT+2*Nghost))
      
!     do k = Kbeg,Kend
!     do j = Jbeg,Jend
!     do i = Ibeg,Iend
!       i_glob(i) = npx*(Iend-Ibeg+1)+i-Nghost
!       j_glob(j) = npy*(Jend-Jbeg+1)+j-Nghost
!       k_glob(k) = k-Nghost
!     enddo
!     enddo
!     enddo
      
      do k = Kbeg,Kend
          do j = Jbeg,Jend
          do i = Ibeg,Iend
              i_glob(i) = npx*(Iend-Ibeg+1)+i-Nghost
              j_glob(j) = npy*(Jend-Jbeg+1)+j-Nghost
              k_glob(k) = ((k-Kbeg)/PDI_INT) + 1
          enddo
          enddo
      enddo
      
     ilower(1) = i_glob(Ibeg)
     ilower(2) = j_glob(Jbeg)
     ilower(3) = k_glob(Kbeg)
!      ilower(3) = 1
     iupper(1) = i_glob(Iend)
     iupper(2) = j_glob(Jend)
     iupper(3) = k_glob(Kend)
!      iupper(3) = Kglob/PDI_INT

     call HYPRE_StructGridSetExtents(grid,ilower,iupper,ierr)
     if(PERIODIC_X.or.PERIODIC_Y) then
       if(PERIODIC_X) then
         periodic_shift(1) = Mglob
       else
         periodic_shift(1) = 0
       endif
       if(PERIODIC_Y) then
         periodic_shift(2) = Nglob
       else
         periodic_shift(2) = 0
       endif
       periodic_shift(3) = 0
       call HYPRE_StructGridSetPeriodic(grid,periodic_shift,ierr)
     endif
     call HYPRE_StructGridAssemble(grid,ierr)
!
!    define the discretization stencil
!
     call HYPRE_StructStencilCreate(ndim,nentries,stencil,ierr)
     do ientry = 1,nentries
       call HYPRE_StructStencilSetElement(stencil,(ientry-1),offsets(ientry,:),ierr)
     enddo

! create matrix object
     call HYPRE_StructMatrixCreate(MPI_COMM_WORLD,grid,stencil,matrix,ierr)
     call HYPRE_StructMatrixInitialize(matrix,ierr)

!    set the matrix coefficient
!
     do i = 1,nentries
       stencil_indices(i) = i-1
     enddo

     allocate(indx(Mloc,Nloc,Kloc))
      indx = Zero
     neq = 0
     do k = Kbeg,Kend,PDI_INT
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       neq = neq+1
       indx(i,j,k) = neq
     enddo
     enddo
     enddo
    
     nvalues = (Iend-Ibeg+1)*(Jend-Jbeg+1)*(Kend-Kbeg+1)*nentries
      nvalues = nvalues / PDI_INT
     allocate(values(nvalues))

     ivalues = 0
     do k = Kbeg,Kend,PDI_INT
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       do n = 1,nentries
         ivalues = ivalues+1
         values(ivalues) = Coef2(indx(i,j,k),n)
       enddo
     enddo
     enddo
     enddo

     call HYPRE_StructMatrixSetBoxValues(matrix,ilower,iupper,nentries,  &
                                  stencil_indices,values,ierr) 
     call HYPRE_StructMatrixAssemble(matrix,ierr)
     !call HYPRE_StructMatrixPrint(matrix,zero,ierr)
!
!    set up struct vectors for b and x
!
     call HYPRE_StructVectorCreate(MPI_COMM_WORLD,grid,vec_b,ierr)
     call HYPRE_StructVectorCreate(MPI_COMM_WORLD,grid,vec_x,ierr)

     call HYPRE_StructVectorInitialize(vec_b,ierr)
     call HYPRE_StructVectorInitialize(vec_x,ierr)
!
! set the vector coefficients
     call HYPRE_StructVectorSetBoxValues(vec_b,ilower,iupper,Rhs2,ierr)   
     call HYPRE_StructVectorAssemble(vec_b,ierr)     
     !call HYPRE_StructVectorPrint(vec_b,zero,ierr)

! initial guess
     allocate(Phi(neqns2))
     do k = Kbeg,Kend,PDI_INT
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Phi(indx(i,j,k)) = P(i,j,k)
     enddo
     enddo
     enddo

     call HYPRE_StructVectorSetBoxValues(vec_x,ilower,iupper,Phi,ierr)
     call HYPRE_StructVectorAssemble(vec_x,ierr)
      !call HYPRE_StructVectorPrint(vec_x,zero,ierr)

! set up and use a solver
     call HYPRE_StructGMRESCreate(MPI_COMM_WORLD,solver,ierr)
     call HYPRE_StructGMRESSetMaxIter(solver,itmax,ierr)
     call HYPRE_StructGMRESSetTol(solver,tol,ierr)
     call HYPRE_StructGMRESSetPrintLevel(solver,0,ierr)
     call HYPRE_StructGMRESSetLogging(solver,0,ierr)

! use symmetric SMG as preconditioner
     n_pre = 1; n_post = 1
     call HYPRE_StructSMGCreate(MPI_COMM_WORLD,precond,ierr)
     call HYPRE_StructSMGSetMemoryUse(precond,0,ierr)
     call HYPRE_StructSMGSetMaxIter(precond,1,ierr)
     call HYPRE_StructSMGSetTol(precond,0.0,ierr)
     call HYPRE_StructSMGSetNumPreRelax(precond,n_pre,ierr)
     call HYPRE_StructSMGSetNumPostRelax(precond,n_post,ierr)
     call HYPRE_StructSMGSetLogging(precond,0,ierr)

! set up preconditioner
     precond_id = 0
     call HYPRE_StructGMRESSetPrecond(solver,precond_id,precond,ierr)
     
! do the setup
     call HYPRE_StructGMRESSetup(solver,matrix,vec_b,vec_x,ierr)
 
! do the solve
     call HYPRE_StructGMRESSolve(solver,matrix,vec_b,vec_x,ierr)

! get results
     call HYPRE_StructVectorGetBoxValues(vec_x,ilower,iupper,Phi,ierr)

!     do k = Kbeg,Kend
!     do j = Jbeg,Jend
!     do i = Ibeg,Iend
!       P(i,j,k) = Phi(indx(i,j,k))
!     enddo
!     enddo
!     enddo
   
      ! update the values in pressure grid
     do k = Kbeg,Kend,PDI_INT
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       P(i,j,k) = Phi(indx(i,j,k)) 
     enddo
     enddo
     enddo

  !    ==== interpolation ====
      if(PDI_INT.ne.1) call interpolate_pressure_pdi

     ! get some info
     !call HYPRE_StructGMRESGetFinalRelati(solver,final_res_norm,ierr)
     !call HYPRE_StructGMRESGetNumIteratio(solver,num_iterations,ierr);
     !
     !if(myid.eq.0) then
     !  write(*,*)'Iterations = ',num_iterations
     !  write(*,*)'Final Relative Residual Norm = ',final_res_norm
     !endif

     ! free memory
     call HYPRE_StructGridDestroy(grid,ierr)
     call HYPRE_StructStencilDestroy(stencil,ierr)
     call HYPRE_StructMatrixDestroy(matrix,ierr)
     call HYPRE_StructVectorDestroy(vec_b,ierr)
     call HYPRE_StructVectorDestroy(vec_x,ierr)
     call HYPRE_StructGMRESDestroy(solver,ierr)
     call HYPRE_StructSMGDestroy(precond,ierr)

     deallocate(indx)
     deallocate(values)
     deallocate(Phi)

# endif
     return
     end subroutine hypre_pres_solver_pdi
# endif
!
!-----------------------------------------------------------------------------------------------
!
!    (3) Subroutine generate_coef_rhs_pdi
!
!    Generate coefficient matrix and rhs for pressure grid (PDI method)
!
!    called by:  poisson_solver_pdi
!
!    Last update: 08/08/2022, Lihan Hu
!
!----------------------------------------------------------------------------------------------
!
# if defined (PRESSUREACC)
     subroutine generate_coef_rhs_pdi
!
     use global
     implicit none
     integer :: i,j,k,neq,n,ic,kk
     real(SP) :: dsigku, dsigkd
     real(SP), dimension(:,:,:), allocatable :: DelxS,DelyS,DelzS,A1
      real(SP), dimension(:,:,:), allocatable :: U_p,V_p,W_p
      real(SP), dimension(:), allocatable :: dsig_p
     integer,  dimension(:,:,:), allocatable :: indx

     allocate(DelxS(Mloc,Nloc,Kloc1))
     allocate(DelyS(Mloc,Nloc,Kloc1))
     allocate(DelzS(Mloc,Nloc,Kloc1))
     allocate(A1(Mloc,Nloc,Kloc1))
     allocate(indx(Mloc,Nloc,Kloc))
      allocate(U_p(Mloc,Nloc,0:(Kglob/PDI_INT)))
      allocate(V_p(Mloc,Nloc,0:(Kglob/PDI_INT)))
      allocate(W_p(Mloc,Nloc,0:(Kglob/PDI_INT)))
      allocate(dsig_p(0:(Kglob/PDI_INT)))

     DelxS = Zero
     DelyS = Zero
     DelzS = Zero
     A1 = Zero

     do k = Kbeg,Kend1
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       DelxS(i,j,k) = (1.-sig(k))/D(i,j)*DelxH(i,j)-sig(k)/D(i,j)*DelxEta(i,j)
       DelyS(i,j,k) = (1.-sig(k))/D(i,j)*DelyH(i,j)-sig(k)/D(i,j)*DelyEta(i,j) 
       DelzS(i,j,k) = 1./D(i,j)

       A1(i,j,k) = DelxS(i,j,k)*DelxS(i,j,k)+DelyS(i,j,k)*DelyS(i,j,k)+  &
            DelzS(i,j,k)*DelzS(i,j,k)
     enddo
     enddo
     enddo
            
      ! generate the average U,V,W for pressure grid    
      U_p = Zero
      V_p = Zero
      W_p = Zero
      dsig_p = Zero
      
      if(PDI_INT==1) then
          dsig_p(0) = dsig(Kbeg-1)
          do j = Jbeg,Jend
          do i = Ibeg,Iend
              U_p(i,j,0) = U(i,j,Kbeg-1)
              V_p(i,j,0) = V(i,j,Kbeg-1)
              W_p(i,j,0) = W(i,j,Kbeg-1)
          enddo
          enddo         
      else
          dsig_p(0) = sig(Kbeg) - sig(Kbeg-Nghost)
          do j = Jbeg,Jend
          do i = Ibeg,Iend
          do kk = 1,Nghost
              U_p(i,j,0) = U_p(i,j,0) + U(i,j,(Kbeg-kk))*dsig(k-kk) / dsig_p(0)
              V_p(i,j,0) = V_p(i,j,0) + V(i,j,(Kbeg-kk))*dsig(k-kk) / dsig_p(0)
              W_p(i,j,0) = W_p(i,j,0) + W(i,j,(Kbeg-kk))*dsig(k-kk) / dsig_p(0)
          enddo
          enddo
          enddo
      endif
            
      neq = 0
      do k = Kbeg,Kend,PDI_INT
          neq = neq + 1
          dsig_p(neq) = sig(k+PDI_INT) - sig(k)
          do j = Jbeg,Jend
          do i = Ibeg,Iend
              do kk = 0,(PDI_INT-1)
                  U_p(i,j,neq) = U_p(i,j,neq) + U(i,j,(k+kk))*dsig(k+kk) / dsig_p(neq)
                  V_p(i,j,neq) = V_p(i,j,neq) + V(i,j,(k+kk))*dsig(k+kk) / dsig_p(neq)
                  W_p(i,j,neq) = W_p(i,j,neq) + W(i,j,(k+kk))*dsig(k+kk) / dsig_p(neq)
              enddo
          enddo
          enddo
      enddo
      
      !generate index for Possion vector (for solving Possion Equation)
      indx = Zero
      neq = 0
      do k = Kbeg,Kend,PDI_INT
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        neq = neq+1
        indx(i,j,k) = neq
      enddo
      enddo 
      enddo
       
     ! generate coefficient matrix and source terms
     Coef2 = Zero
     Rhs2 = Zero
      neq = 0
     do k = Kbeg,Kend,PDI_INT
      neq = neq + 1
     do j = Jbeg,Jend
     do i = Ibeg,Iend
!       if(indx(i,j,k)==0) cycle
       dsigku = dsig_p(neq)
       dsigkd = dsig_p(neq-1)

# if defined (POROUSMEDIA)
          Rhs2(indx(i,j,k)) = -((Uf(i+1,j,k)-Uf(i-1,j,k))/(2.0*dx)+(U_p(i,j,neq)-U_p(i,j,neq-1))/(0.5*(dsigku+dsigkd))*  &
              DelxS(i,j,k)+(Vf(i,j+1,k)-Vf(i,j-1,k))/(2.0*dy)+(V_p(i,j,neq)-V_p(i,j,neq-1))/(0.5*(dsigku+dsigkd))*  &
              DelyS(i,j,k)+(W_p(i,j,neq)-W_p(i,j,neq-1))/(0.5*(dsigku+dsigkd))*DelzS(i,j,k)-SourceC(i,j))*Rho0/dt* &
              (1+Cp_Por(i,j,k))
# else
          Rhs2(indx(i,j,k)) = -((Uf(i+1,j,k)-Uf(i-1,j,k))/(2.0*dx)+(U_p(i,j,neq)-U_p(i,j,neq-1))/(0.5*(dsigku+dsigkd))*  &
              DelxS(i,j,k)+(Vf(i,j+1,k)-Vf(i,j-1,k))/(2.0*dy)+(V_p(i,j,neq)-V_p(i,j,neq-1))/(0.5*(dsigku+dsigkd))*  &
              DelyS(i,j,k)+(W_p(i,j,neq)-W_p(i,j,neq-1))/(0.5*(dsigku+dsigkd))*DelzS(i,j,k)-SourceC(i,j))*Rho0/dt
# endif
     
!# if defined (POROUSMEDIA)
!       Rhs2(indx(i,j,k)) = -((Uf(i+1,j,k)-Uf(i-1,j,k))/(2.0*dx)+(U(i,j,k)-U(i,j,k-1))/(0.5*(dsig(k)+dsig(k-1)))*  &
!              DelxS(i,j,k)+(Vf(i,j+1,k)-Vf(i,j-1,k))/(2.0*dy)+(V(i,j,k)-V(i,j,k-1))/(0.5*(dsig(k)+dsig(k-1)))*  &
!              DelyS(i,j,k)+(W(i,j,k)-W(i,j,k-1))/(0.5*(dsig(k)+dsig(k-1)))*DelzS(i,j,k)-SourceC(i,j))*Rho0/dt* &
!              (1+Cp_Por(i,j,k))
!# else
!       Rhs2(indx(i,j,k)) = -((Uf(i+1,j,k)-Uf(i-1,j,k))/(2.0*dx)+(U(i,j,k)-U(i,j,k-1))/(0.5*(dsig(k)+dsig(k-1)))*  &
!              DelxS(i,j,k)+(Vf(i,j+1,k)-Vf(i,j-1,k))/(2.0*dy)+(V(i,j,k)-V(i,j,k-1))/(0.5*(dsig(k)+dsig(k-1)))*  &
!              DelyS(i,j,k)+(W(i,j,k)-W(i,j,k-1))/(0.5*(dsig(k)+dsig(k-1)))*DelzS(i,j,k)-SourceC(i,j))*Rho0/dt
!# endif

       ! Coef2[1~15] = [a8     / a9       a10      / a11         a12       / a13     / a14        a15        
       !               / a7        a6      / a5          a4        / a3      / a2          a1        ]
       ! the points =  [(i,j,k)/(i+1,j,k) (i,j+1,k)/(i,j-1,k+1) (i-1,j,k+1)/(i,j,k+1)/(i+1,j,k+1) (i,j+1,k+1)
       !                (i-1,j,k) (i,j-1,k)/(i,j+1,k-1) (i+1,j,k-1)/(i,j,k-1)/(i-1,j,k-1) (i,j-1,k-1)]

       Coef2(indx(i,j,k),1) = (2./(dx*dx)+2./(dy*dy)+A1(i,j,k)/(0.5*(dsigku+dsigkd)*dsigku)+  &
                               A1(i,j,k)/(0.5*(dsigku+dsigkd)*dsigkd))
       Coef2(indx(i,j,k),2) = -1./(dx*dx)
       Coef2(indx(i,j,k),3) = -1./(dy*dy)
       Coef2(indx(i,j,k),4) = (DelyS(i,j-1,k)/(2.*dy*(dsigku+dsigkd))+DelyS(i,j,k)/(2.*dy*(dsigku+dsigkd)))   
       Coef2(indx(i,j,k),5) = (DelxS(i-1,j,k)/(2.*dx*(dsigku+dsigkd))+DelxS(i,j,k)/(2.*dx*(dsigku+dsigkd)))
       Coef2(indx(i,j,k),6) = -A1(i,j,k)/(0.5*(dsigku+dsigkd)*dsigku)
       Coef2(indx(i,j,k),7) = -(DelxS(i+1,j,k)/(2.*dx*(dsigku+dsigkd))+DelxS(i,j,k)/(2.*dx*(dsigku+dsigkd)))
       Coef2(indx(i,j,k),8) = -(DelyS(i,j+1,k)/(2.*dy*(dsigku+dsigkd))+DelyS(i,j,k)/(2.*dy*(dsigku+dsigkd)))
       Coef2(indx(i,j,k),9) = -1./(dx*dx)
       Coef2(indx(i,j,k),10) = -1./(dy*dy)
       Coef2(indx(i,j,k),11) = (DelyS(i,j+1,k)/(2.*dy*(dsigku+dsigkd))+DelyS(i,j,k)/(2.*dy*(dsigku+dsigkd)))
       Coef2(indx(i,j,k),12) = (DelxS(i+1,j,k)/(2.*dx*(dsigku+dsigkd))+DelxS(i,j,k)/(2.*dx*(dsigku+dsigkd)))
       Coef2(indx(i,j,k),13) = -A1(i,j,k)/(0.5*(dsigku+dsigkd)*dsigkd)
       Coef2(indx(i,j,k),14) = -(DelxS(i-1,j,k)/(2.*dx*(dsigku+dsigkd))+DelxS(i,j,k)/(2.*dx*(dsigku+dsigkd)))
       Coef2(indx(i,j,k),15) = -(DelyS(i,j-1,k)/(2.*dy*(dsigku+dsigkd))+DelyS(i,j,k)/(2.*dy*(dsigku+dsigkd)))

!       Coef2(indx(i,j,k),1) = (2./(dx*dx)+2./(dy*dy)+A1(i,j,k)/(0.5*(dsig(k)+dsig(k-1))*dsig(k))+  &
!                              A1(i,j,k)/(0.5*(dsig(k)+dsig(k-1))*dsig(k-1)))
!       Coef2(indx(i,j,k),2) = -1./(dx*dx)
!       Coef2(indx(i,j,k),3) = -1./(dy*dy)
!       Coef2(indx(i,j,k),4) = (DelyS(i,j-1,k)/(2.*dy*(dsig(k)+dsig(k-1)))+DelyS(i,j,k)/(2.*dy*(dsig(k)+dsig(k-1))))   
!       Coef2(indx(i,j,k),5) = (DelxS(i-1,j,k)/(2.*dx*(dsig(k)+dsig(k-1)))+DelxS(i,j,k)/(2.*dx*(dsig(k)+dsig(k-1))))
!       Coef2(indx(i,j,k),6) = -A1(i,j,k)/(0.5*(dsig(k)+dsig(k-1))*dsig(k))
!       Coef2(indx(i,j,k),7) = -(DelxS(i+1,j,k)/(2.*dx*(dsig(k)+dsig(k-1)))+DelxS(i,j,k)/(2.*dx*(dsig(k)+dsig(k-1))))
!       Coef2(indx(i,j,k),8) = -(DelyS(i,j+1,k)/(2.*dy*(dsig(k)+dsig(k-1)))+DelyS(i,j,k)/(2.*dy*(dsig(k)+dsig(k-1))))
!       Coef2(indx(i,j,k),9) = -1./(dx*dx)
!       Coef2(indx(i,j,k),10) = -1./(dy*dy)
!       Coef2(indx(i,j,k),11) = (DelyS(i,j+1,k)/(2.*dy*(dsig(k)+dsig(k-1)))+DelyS(i,j,k)/(2.*dy*(dsig(k)+dsig(k-1))))
!       Coef2(indx(i,j,k),12) = (DelxS(i+1,j,k)/(2.*dx*(dsig(k)+dsig(k-1)))+DelxS(i,j,k)/(2.*dx*(dsig(k)+dsig(k-1))))
!       Coef2(indx(i,j,k),13) = -A1(i,j,k)/(0.5*(dsig(k)+dsig(k-1))*dsig(k-1))
!       Coef2(indx(i,j,k),14) = -(DelxS(i-1,j,k)/(2.*dx*(dsig(k)+dsig(k-1)))+DelxS(i,j,k)/(2.*dx*(dsig(k)+dsig(k-1))))
!       Coef2(indx(i,j,k),15) = -(DelyS(i,j-1,k)/(2.*dy*(dsig(k)+dsig(k-1)))+DelyS(i,j,k)/(2.*dy*(dsig(k)+dsig(k-1))))

     enddo
     enddo
     enddo

     ! fyshi added boundary conditions at masks face 02/15/2013
!     do i = Ibeg+1,Iend-1
!     do j = Jbeg+1,Jend-1
!     do k = Kbeg,Kend
!       if(mask(i,j)==0) then
!         ! left 
!         if(mask(i+1,j)==1) then
!           ic = indx(I+1,j,k)
!           Coef2(ic,1) = Coef2(ic,1)+Coef2(ic,9)
!           Coef2(ic,6) = Coef2(ic,6)+Coef2(ic,5)
!           Coef2(ic,13) = Coef2(ic,13)+Coef2(ic,14)
!           Coef2(ic,9) = Zero
!           Coef2(ic,5) = Zero
!           Coef2(ic,14) = Zero
!         ! right 
!         elseif(mask(i-1,j)==1) then
!           ic = indx(I-1,j,k)
!           Coef2(ic,1) = Coef2(ic,1)+Coef2(ic,2)
!           Coef2(ic,6) = Coef2(ic,6)+Coef2(ic,7)
!           Coef2(ic,13) = Coef2(ic,13)+Coef2(ic,12)
!           Coef2(ic,2) = Zero
!           Coef2(ic,7) = Zero
!           Coef2(ic,12) = Zero
!         ! south
!         elseif(mask(i,j+1)==1) then
!           ic = indx(i,J+1,k)
!           Coef2(ic,1) = Coef2(ic,1)+Coef2(ic,10)
!           Coef2(ic,6) = Coef2(ic,6)+Coef2(ic,4)
!           Coef2(ic,13) = Coef2(ic,13)+Coef2(ic,15)
!           Coef2(ic,10) = Zero
!           Coef2(ic,4) = Zero
!           Coef2(ic,15) = Zero
!         ! north
!         elseif(mask(i,j-1)==1) then
!           ic = indx(i,J-1,k)
!           Coef2(ic,1) = Coef2(ic,1)+Coef2(ic,3)
!           Coef2(ic,6) = Coef2(ic,6)+Coef2(ic,8)
!           Coef2(ic,13) = Coef2(ic,13)+Coef2(ic,11)
!           Coef2(ic,3) = Zero
!           Coef2(ic,8) = Zero
!           Coef2(ic,11) = Zero
!         endif ! end mask+1=1 
!       endif ! end mask=0
!     enddo
!     enddo
!     enddo

# if defined (OBSTACLE)
     do i = Ibeg+1,Iend-1
     do j = Jbeg+1,Jend-1
     do k = Kbeg+1,Kend-1
       if(set_flag(i,j,k)==1) then
        ! left 
         if(set_flag(i+1,j,k)==0) then
           ic = indx(I+1,j,k)
           if(ic==0) cycle
           Coef2(ic,1) = Coef2(ic,1)+Coef2(ic,9)
           Coef2(ic,9) = Zero
         ! right 
         elseif(set_flag(i-1,j,k)==0) then
           ic = indx(I-1,j,k)
           if(ic==0) cycle
           Coef2(ic,1) = Coef2(ic,1)+Coef2(ic,2)
           Coef2(ic,2) = Zero
         ! south
         elseif(set_flag(i,j+1,k)==0) then
           ic = indx(i,J+1,k)
           if(ic==0) cycle
           Coef2(ic,1) = Coef2(ic,1)+Coef2(ic,10)
           Coef2(ic,10) = Zero
         ! north
         elseif(set_flag(i,j-1,k)==0) then
           ic = indx(i,J-1,k)
           if(ic==0) cycle
           Coef2(ic,1) = Coef2(ic,1)+Coef2(ic,3)
           Coef2(ic,3) = Zero
         ! bottom
         elseif(set_flag(i,j,k+1)==0) then
           ic = indx(i,j,k+1)
           if(ic==0) cycle
           Coef2(ic,1) = Coef2(ic,1)+Coef2(ic,13)
           Coef2(ic,13) = Zero
         ! upper
         elseif(set_flag(i,j,k-1)==0) then
           ic = indx(i,j,k-1)
           if(ic==0) cycle
           Coef2(ic,1) = Coef2(ic,1)+Coef2(ic,6)
           Coef2(ic,6) = Zero
         endif  
       endif
     enddo
     enddo
     enddo
# endif
     !! boundary conditions
     ! left side (Neumann boundary condition: normal pressure gradient is zero)
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
! added by cheng for nesting, search (COUPLING) to find rest in this subroutine
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif
     i = Ibeg
     do k = Kbeg,Kend,PDI_INT
     do j = Jbeg,Jend
       ic = indx(i,j,k)
       Coef2(ic,1) = Coef2(ic,1) + Coef2(ic,9)
       Coef2(ic,6) = Coef2(ic,6) + Coef2(ic,5)
       Coef2(ic,13) = Coef2(ic,13) + Coef2(ic,14)
       Coef2(ic,9) = Zero
       Coef2(ic,5) = Zero
       Coef2(ic,14) = Zero
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

     ! right side (Neumann boundary condition: normal pressure gradient is zero)
# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
     i = Iend
     do k = Kbeg,Kend,PDI_INT
     do j = Jbeg,Jend
       ic = indx(i,j,k)
       Coef2(ic,1) = Coef2(ic,1)+Coef2(ic,2)
       Coef2(ic,6) = Coef2(ic,6)+Coef2(ic,7)
       Coef2(ic,13) = Coef2(ic,13)+Coef2(ic,12)
       Coef2(ic,2) = Zero
       Coef2(ic,7) = Zero
       Coef2(ic,12) = Zero
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

     ! front side (Neumann boundary condition: normal pressure gradient is zero)
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
     j = Jbeg
     do k = Kbeg,Kend,PDI_INT
     do i = Ibeg,Iend
       ic = indx(i,j,k)
       Coef2(ic,1) = Coef2(ic,1)+Coef2(ic,10)
       Coef2(ic,6) = Coef2(ic,6)+Coef2(ic,4)
       Coef2(ic,13) = Coef2(ic,13)+Coef2(ic,15)
       Coef2(ic,10) = Zero
       Coef2(ic,4) = Zero
       Coef2(ic,15) = Zero
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif

     ! back side (Neumann boundary condition: normal pressure gradient is zero)
# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
     j = Jend
     do k = Kbeg,Kend,PDI_INT
     do i = Ibeg,Iend
       ic = indx(i,j,k)
       Coef2(ic,1) = Coef2(ic,1)+Coef2(ic,3)
       Coef2(ic,6) = Coef2(ic,6)+Coef2(ic,8)
       Coef2(ic,13) = Coef2(ic,13)+Coef2(ic,11)
       Coef2(ic,3) = Zero
       Coef2(ic,8) = Zero
       Coef2(ic,11) = Zero
     enddo
     enddo
# if defined (COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
     endif
# endif
     ! bottom side (Neumann boundary condition: normal pressure gradient is zero 
      !             ... Shi Jian's Doctoral Thesis: the derivative of dynamic 
      !             ... pressure value along the ¦Ò direction of the bottom is zero 
      !             ... without taking the topographical change (over time) into account)
     k = Kbeg
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       ic = indx(i,j,k)
       dsigku = dsig_p(1)
       dsigkd = dsig_p(0)
# if defined (LANDSLIDE)
       if(SlideType(1:5)=='RIGID') then ! modified by Cheng to identify 2d and 3d landslide
         Rhs2(ic) = Rhs2(ic)+Rho0*(dsigku+dsigkd)*(Coef2(ic,13)*D(i,j)*Delt2H(i,j)+ &
            Coef2(ic,12)*D(i+1,j)*Delt2H(i+1,j)+Coef2(ic,11)*D(i,j+1)*Delt2H(i,j+1)+ &
            Coef2(ic,14)*D(i-1,j)*Delt2H(i-1,j)+Coef2(ic,15)*D(i,j-1)*Delt2H(i,j-1))
       endif
# endif

!added by Cheng for fluid slide and deformable slide
# if defined (FLUIDSLIDE) || defined (DEFORMABLESLIDE)
       Rhs2(ic) = Rhs2(ic)+Rho0*(dsigku+dsigkd)*(Coef2(ic,13)*D(i,j)*Delt2H(i,j)+ &
            Coef2(ic,12)*D(i+1,j)*Delt2H(i+1,j)+Coef2(ic,11)*D(i,j+1)*Delt2H(i,j+1)+ &
            Coef2(ic,14)*D(i-1,j)*Delt2H(i-1,j)+Coef2(ic,15)*D(i,j-1)*Delt2H(i,j-1))
# endif

# if defined (LANDSLIDE_COMPREHENSIVE)
       Rhs2(ic) = Rhs2(ic)+Rho0*(dsigku+dsigkd)*(Coef2(ic,13)*D(i,j)*Delt2H(i,j)+ &
            Coef2(ic,12)*D(i+1,j)*Delt2H(i+1,j)+Coef2(ic,11)*D(i,j+1)*Delt2H(i,j+1)+ &
            Coef2(ic,14)*D(i-1,j)*Delt2H(i-1,j)+Coef2(ic,15)*D(i,j-1)*Delt2H(i,j-1))
# endif
! end landslide comprehensive

# if defined (TWOLAYERSLIDE)
       if(D(i,j)>0.05) then
         Rhs2(ic) = Rhs2(ic)+Rho0*(dsigku+dsigkd)*(Coef2(ic,13)*D(i,j)*Delt2H(i,j)+ &             
            Coef2(ic,12)*D(i+1,j)*Delt2H(i+1,j)+Coef2(ic,11)*D(i,j+1)*Delt2H(i,j+1)+ &
            Coef2(ic,14)*D(i-1,j)*Delt2H(i-1,j)+Coef2(ic,15)*D(i,j-1)*Delt2H(i,j-1))
       endif
# endif

       Coef2(ic,6) = Coef2(ic,6)+Coef2(ic,13)
       Coef2(ic,7) = Coef2(ic,7)+Coef2(ic,12)
       Coef2(ic,8) = Coef2(ic,8)+Coef2(ic,11)
       Coef2(ic,5) = Coef2(ic,5)+Coef2(ic,14)
       Coef2(ic,4) = Coef2(ic,4)+Coef2(ic,15)
       Coef2(ic,13) = Zero
       Coef2(ic,12) = Zero
       Coef2(ic,11) = Zero
       Coef2(ic,14) = Zero
       Coef2(ic,15) = Zero
     enddo
     enddo

     ! top side (Dirichlet boundary)
     k = Kend1 - PDI_INT
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       ic = indx(i,j,k)
       Coef2(ic,4) = Zero
       Coef2(ic,5) = Zero
       Coef2(ic,6) = Zero
       Coef2(ic,7) = Zero
       Coef2(ic,8) = Zero
     enddo
     enddo
       
     ! take (i=2,j=2,k=2) to obtain the diagonal information (what is the JCoef??)
     JCoef2(1) = indx(Ibeg+1,Jbeg+1,Kbeg+PDI_INT)-indx(Ibeg+1,Jbeg+1,Kbeg+PDI_INT)  ! (i,j,k)
     JCoef2(2) = indx(Ibeg+2,Jbeg+1,Kbeg+PDI_INT)-indx(Ibeg+1,Jbeg+1,Kbeg+PDI_INT)  ! (i+1,j,k) 
     JCoef2(3) = indx(Ibeg+1,Jbeg+2,Kbeg+PDI_INT)-indx(Ibeg+1,Jbeg+1,Kbeg+PDI_INT)  ! (i,j+1,k)
     JCoef2(4) = indx(Ibeg+1,Jbeg,Kbeg+2*PDI_INT)-indx(Ibeg+1,Jbeg+1,Kbeg+PDI_INT)    ! (i,j-1,k+1)
     JCoef2(5) = indx(Ibeg,Jbeg+1,Kbeg+2*PDI_INT)-indx(Ibeg+1,Jbeg+1,Kbeg+PDI_INT)    ! (i-1,j,k+1)
     JCoef2(6) = indx(Ibeg+1,Jbeg+1,Kbeg+2*PDI_INT)-indx(Ibeg+1,Jbeg+1,Kbeg+PDI_INT)  ! (i,j,k+1)
     JCoef2(7) = indx(Ibeg+2,Jbeg+1,Kbeg+2*PDI_INT)-indx(Ibeg+1,Jbeg+1,Kbeg+PDI_INT)  ! (i+1,j,k+1)
     JCoef2(8) = indx(Ibeg+1,Jbeg+2,Kbeg+2*PDI_INT)-indx(Ibeg+1,Jbeg+1,Kbeg+PDI_INT)  ! (i,j+1,k+1)
     JCoef2(9) = indx(Ibeg,Jbeg+1,Kbeg+PDI_INT)-indx(Ibeg+1,Jbeg+1,Kbeg+PDI_INT)    ! (i-1,j,k)
     JCoef2(10) = indx(Ibeg+1,Jbeg,Kbeg+PDI_INT)-indx(Ibeg+1,Jbeg+1,Kbeg+PDI_INT)   ! (i,j-1,k)
     JCoef2(11) = indx(Ibeg+1,Jbeg+2,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+PDI_INT)   ! (i,j+1,k-1)
     JCoef2(12) = indx(Ibeg+2,Jbeg+1,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+PDI_INT)   ! (i+1,j,k-1)
     JCoef2(13) = indx(Ibeg+1,Jbeg+1,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+PDI_INT)   ! (i,j,k-1)
     JCoef2(14) = indx(Ibeg,Jbeg+1,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+PDI_INT)     ! (i-1,j,k-1)
     JCoef2(15) = indx(Ibeg+1,Jbeg,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+PDI_INT)     ! (i,j-1,k-1)

     deallocate(DelxS)
     deallocate(DelyS)
     deallocate(DelzS)
     deallocate(A1) 
     deallocate(indx)
      deallocate(U_p)
      deallocate(V_p)
      deallocate(W_p)
      deallocate(dsig_p)

     return
     end subroutine generate_coef_rhs_pdi
# endif
!
!-----------------------------------------------------------------------------------------------
!
!    (4) Subroutine interpolate_pressure_pdi
!
!    Interpolate the Pressure at PDI-grid to that in fined grids (PDI method)
!
!    called by:  poisson_solver_pdi
!
!    Last update: 10/08/2022, Lihan Hu
!
!----------------------------------------------------------------------------------------------
!
# if defined (PRESSUREACC)
subroutine interpolate_pressure_pdi
      use global
      implicit none
      integer :: i,j,k,imask,kk,neq,neq2,Nlen
      real(SP), dimension(:), allocatable :: Acoef,Bcoef,Ccoef,Xsol,Rhs0,hsig
      real(SP) :: akk,bkk,sigkk1,sigkk2,Pkk1,Pkk2
      logical  :: Linear_Interp=.False.

      Nlen = Kglob/PDI_INT
      allocate(Acoef(Nlen+1))
      allocate(Bcoef(Nlen+1))
      allocate(Ccoef(Nlen+1))
      allocate(Xsol(Nlen+1))
      allocate(Rhs0(Nlen+1))
      allocate(hsig(Nlen))
      
      neq = 0
      do k = Kbeg,Kend,PDI_INT
          neq = neq + 1
          hsig(neq) = sig(k+PDI_INT) - sig(k)
      enddo
      
!!    ==== interpolation ====
     if(Linear_Interp)then
      ! ---- linear interpolation ----
      !       p_kk(sig) = a_kk + b_kk*(sig-sigkk)
      do j = Jbeg,Jend
      do i = Ibeg,Iend
          P(i,j,Kend1) = 0   ! boundary condition at the free surface (P = 0)
          do k = Kbeg,Kend,PDI_INT
              akk = P(i,j,k)
              bkk = (P(i,j,k+PDI_INT)-P(i,j,k))/(sig(k+PDI_INT)-sig(k))
              do kk = 1,(PDI_INT-1)
                  P(i,j,k+kk) = akk + bkk*(sig(k+kk) - sig(k))
              enddo
          enddo
      enddo
      enddo
     else
      !! ---- cubic spline / to construct the third-order polynomials ----
      do j = Jbeg,Jend
      do i = Ibeg,Iend
       P(i,j,Kend1) = 0   ! boundary condition at the free surface (P = 0)       
       Acoef = Zero
       Bcoef = Zero
       Ccoef = Zero
       Rhs0  = Zero
       
       neq = 0
       do k = Kbeg,Kend1,PDI_INT
           neq = neq + 1
           if(k == Kbeg)then
              Acoef(neq) = 0.0
              Bcoef(neq) = 2.0
              Ccoef(neq) = 1.0
              Rhs0(neq)  = 6.0/(hsig(neq))*( (P(i,j,k+PDI_INT)-P(i,j,k))/hsig(neq) )
           elseif(k == Kend1)then
              Acoef(neq) = Small
              Bcoef(neq) = 2.0
              Ccoef(neq) = 0.0
              Rhs0(neq)  = 0.0
           else
              Acoef(neq) = (hsig(neq-1))/(hsig(neq-1)+hsig(neq))
              Bcoef(neq) = 2.0
              Ccoef(neq) = (hsig(neq))/(hsig(neq-1)+hsig(neq))
              Rhs0(neq)  = 6.0/(hsig(neq-1)+hsig(neq)) * &
                           ((P(i,j,k+PDI_INT)-P(i,j,k))/hsig(neq)-(P(i,j,k)-P(i,j,k-PDI_INT))/hsig(neq-1))
           endif      
       enddo
       
!       neq = 0
!       neq2 = Nlen + 1
!       do k = Kend1,Kbeg,(-PDI_INT)
!           neq = neq + 1
!           neq2 = neq2 - 1
!           if(k == Kbeg)then  ! here, neq = Nlen+1, neq2 = 0
!              Acoef(neq) = 1.0
!              Bcoef(neq) = 2.0
!              Ccoef(neq) = 0.0
!              Rhs0(neq)  = 6.0/(hsig(1))*( 0.0 - (P(i,j,k)-P(i,j,k+PDI_INT))/hsig(1) )
!           elseif(k == Kend1)then  ! here, neq = 1
!              Acoef(neq) = 0.0
!              Bcoef(neq) = 2.0
!              Ccoef(neq) = Small
!              Rhs0(neq)  = 0.0
!           else
!              Acoef(neq) = (hsig(neq2+1))/(hsig(neq2+1)+hsig(neq2))
!              Bcoef(neq) = 2.0
!              Ccoef(neq) = (hsig(neq2))/(hsig(neq2+1)+hsig(neq2))
!              Rhs0(neq)  = 6.0/(hsig(neq2+1)+hsig(neq2)) * &
!                           ((P(i,j,k-PDI_INT)-P(i,j,k))/hsig(neq2)-(P(i,j,k)-P(i,j,k+PDI_INT))/hsig(neq2+1))
!           endif      
!       enddo
       
      Xsol = Zero
      call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)
      
      neq = 0
      do k = Kbeg,Kend,PDI_INT
          neq = neq + 1
          sigkk1 = sig(k)
          sigkk2 = sig(k+PDI_INT)
          Pkk1 = P(i,j,k)
          Pkk2 = P(i,j,k+PDI_INT)
          do kk = 1,(PDI_INT-1)
              P(i,j,k+kk) = Xsol(neq)/(6.0*hsig(neq))* (sigkk2-sig(k+kk))**(3.0) &
                          + Xsol(neq+1)/(6.0*hsig(neq))* (sig(k+kk)-sigkk1)**(3.0) &
                          + (Pkk1 - Xsol(neq)/6.0*hsig(neq)*hsig(neq))*(sigkk2-sig(k+kk))/hsig(neq) &
                          + (Pkk2 - Xsol(neq+1)/6.0*hsig(neq)*hsig(neq))*(sig(k+kk)-sigkk1)/hsig(neq)
          enddo         
      enddo     
      
!      neq = 0
!      neq2 = Nlen + 1
!      do k = Kend1,(Kbeg+1),(-PDI_INT)
!          neq = neq + 1
!          neq2 = neq2 - 1
!          sigkk1 = sig(k)
!          sigkk2 = sig(k-PDI_INT)
!          Pkk1 = P(i,j,k)
!          Pkk2 = P(i,j,k-PDI_INT)
!          do kk = 1,(PDI_INT-1)
!              P(i,j,k-kk) = Xsol(neq)/(6.0*hsig(neq2))* (sigkk2-sig(k-kk))**(3.0) &
!                          + Xsol(neq+1)/(6.0*hsig(neq2))* (sig(k-kk)-sigkk1)**(3.0) &
!                          + (Pkk1 - Xsol(neq)/6.0*hsig(neq2)*hsig(neq2))*(sigkk2-sig(k-kk))/hsig(neq2) &
!                          + (Pkk2 - Xsol(neq+1)/6.0*hsig(neq2)*hsig(neq2))*(sig(k-kk)-sigkk1)/hsig(neq2)
!          enddo         
!      enddo   
      
      enddo
      enddo
      
     endif
      
      deallocate(Acoef)
      deallocate(Bcoef)
      deallocate(Ccoef)
      deallocate(Xsol)
      deallocate(Rhs0)

end subroutine interpolate_pressure_pdi
# endif
!
!
!    End of pressure_pdi.F
!
!--------------------------------------------------------------------------------------------------------
