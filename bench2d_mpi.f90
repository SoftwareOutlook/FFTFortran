PROGRAM commandline
  use omp_lib
  use, intrinsic :: iso_c_binding
  use mkl_cdft 
  use mpi
  !  include '/usr/include/fftw3.f03'                                            
  include '/opt/cray/fftw/default/ivybridge/include/fftw3-mpi.f03'

  Integer, Parameter :: wp = Selected_Real_Kind(15,307)  ! double real

  Integer, Parameter :: ip4 = Selected_Int_Kind(4)  ! integer*4

  integer   mpi_kind
  parameter (mpi_kind = wp)



  INTEGER :: nargs, n1, n2, n3, nq, fftlib ! Input arguments
  integer :: stat ! Allocat/deallocate stat 
  integer :: flag ! Error flag
  integer :: i, j, k, qq, m ! indices
  integer :: comm, ierr, my_id, no_proc ! mpi variables
  real(kind=wp) :: xo, yo, zo, a1, b1, c1, r  ! Used in definition of ellipsoid
  CHARACTER(LEN=100) :: option1, option2, option3, option4, option5 
  ! For reading inputs
  real (kind=wp), allocatable :: A(:,:,:) ! Array A
  real (kind=wp), allocatable :: B(:,:,:,:) ! B(:,:,:,i) is cuboid B_i
  real (kind=wp), allocatable :: C(:,:,:) ! C(:,:,:) is cuboid C_i
  real (kind=wp) :: s1,s2 ! used to define B
  real (kind=wp) :: tm1, tm2, tm_fft_init, tm_fft, tm_ifft_init, tm_ifft

  real (kind=wp) :: tm_fft_init_tot, tm_fft_tot, tm_ifft_init_tot, &
       tm_ifft_tot,  tm_ifft_init_max, tm_fft_init_max
  logical :: check

  call mpi_init(ierr)

  comm = mpi_comm_world
  ! find out process id and total no. processes
  call mpi_comm_rank(comm,my_id,ierr)
  call mpi_comm_size(comm,no_proc,ierr)
  write(*,*) 'my rank is ', my_id, ' out of ', no_proc, ' processes'

  !Read input from the command line
  !nargs = IARGC()
  nargs = COMMAND_ARGUMENT_COUNT() 
  IF (nargs .ne. 5) THEN
     goto 10
  ELSE
     CALL GET_COMMAND_ARGUMENT(1,option1) !Grab the first command line argument
     ! and store it in temp variable 'option1'
     CALL GET_COMMAND_ARGUMENT(2,option2) !Grab the 2nd command line argument
     ! and store it in temp variable 'option2'
     CALL GET_COMMAND_ARGUMENT(3,option3) !Grab the 3rd command line argument
     ! and store it in temp variable 'option3'
     CALL GET_COMMAND_ARGUMENT(4,option4) !Grab the 4th command line argument
     ! and store it in temp variable 'option4'
     !  1: FFTE  NOT included due to bugs in real-complex fft
     !  2: FFTW
     !  3: MKL
     !  4: P3DFFT NOT included due to fussy installation scripts
     !  5: P3DFFT++ NOT included due to poor installation instructions
     CALL GET_COMMAND_ARGUMENT(5,option5) !Grab the 5th command line argument
     ! and store it in temp variable 'option5'


     read(option1,*) n1 !Now convert string to integer
     read(option2,*) n2 !Now convert string to integer
     read(option3,*) n3 !Now convert string to integer
     read(option4,*) nq !Now convert string to integer
     read(option5,*) fftlib !Now convert string to integer
     write(*,'(a,i8)') "Variable n1 = ", n1
     write(*,'(a,i8)') "Variable n2 = ", n2
     write(*,'(a,i8)') "Variable n3 = ", n3
     write(*,'(a,i8)') "Variable nq = ", nq
     write(*,'(a,i8)') "Variable fftlib = ", fftlib

     if (n1.lt.1 .or. n2.lt.1 .or. n3.lt.1 .or. nq.lt.1 .or. fftlib .lt. 2 &
          .or. fftlib .gt. 5) then
        goto 10
     endif

  ENDIF

  !  nthreads = omp_get_num_threads()
  !  write(*,*) 'nthreads',nthreads

  ! Allocate arrays
  allocate(A(n1,n2,n3),stat=stat)
  if (stat .ne. 0) then
     write(*, '(a)') "Error allocating A"
     goto 100
  end if

  allocate(B(n1,n2,n3,nq),stat=stat)
  if (stat .ne. 0) then
     write(*, '(a)') "Error allocating B"
     deallocate(A)
     goto 100
  end if
  allocate(C(n1,n2,n3),stat=stat)
  if (stat .ne. 0) then   
     write(*, '(a)') "Error allocating C"
     deallocate(A,B)
     goto 100
  end if

  ! Initialise total times
  tm_fft_init_tot=0.0_wp
  tm_fft_tot=0.0_wp 
  tm_ifft_init_tot=0.0_wp 
  tm_ifft_tot=0.0_wp
  tm_ifft_init_max=0.0_wp
  tm_fft_init_max=0.0_wp


  ! Set A
  if (n1.eq.1) then
     xo = real(n1,wp)
  else
     xo = 0.6*real(n1,wp)
  end if
  if (n2.eq.1) then
     yo = real(n2,kind=wp)
  else 
     yo = 0.4*real(n2,wp)
  end if
  if (n3.eq.1) then
     zo = real(n3,wp)
  else
     zo = real(n3,wp)/3.0
  end if
  a1 = 0.3*real(n1,wp)
  b1 = 0.35*real(n2,wp)
  c1 = real(n3,wp)/3.0_wp
  !$ tm1=omp_get_wtime()
  do i=1,n1
     do j=1,n2
        do k=1,n3
           r = ((real(i,wp)-xo)/a1)**2 + ((real(j,wp)-yo)/b1)**2 + &
                ((real(k,wp)-zo)/c1)**2
           if (r .le. 1.0_wp) then
              A(i,j,k) = r + 0.5_wp
           else
              A(i,j,k) = 0.5_wp
           end if
        end do
     end do
  end do
  !$ tm2=omp_get_wtime()
  write(*,*) 'Set-up A time=', tm2-tm1

  ! Set B  
  m = 0
  !$ tm1=omp_get_wtime()

  !$OMP PARALLEL PRIVATE (i,j,k,qq) SHARED (n1,n2,n3,B,nq)
  !$OMP DO COLLAPSE(4)
  do i=1,n1
     do j=1,n2
        do k=1,n3
           do qq=1,nq
              B(i,j,k,qq) = (real(i*qq,kind=wp)/real(j*k,kind=wp))
           end do
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL    

  !$ tm2=omp_get_wtime()
  write(*,*) 'Set-up B time (no norm)=', tm2-tm1


  ! Set-up each 2D slice and perform FFT
  ! Each slice formed in C(:,:,:) by performing element-wise multiplaction of 
  ! A with B(:,:,:,qq)
  do qq=1,nq
     !$  tm1=omp_get_wtime()
     !$OMP PARALLEL DO PRIVATE (j,k,s1,s2)
     do i=1,n1
        do j=1,n2
           do k=1,n3
              s1 = A(i,j,k)
              s2 = B(i,j,k,qq)
              C(i,j,k) = s1*s2

           end do
        end do
     end do
     !$OMP END PARALLEL DO


     !$ tm2=omp_get_wtime()
     write(*,*) 'Set-up Cq time=', tm2-tm1

     !  write(*,*) 'C(n1/2,n2/2,n3/2)', C


     ! Perform FFT on each 2D slice
     check=.true.
     call fft_bench(comm,n1,n2,n3,C,fftlib,check,flag,tm_fft_init,tm_fft,&
          tm_ifft_init,tm_ifft)
     write(*,'(a10,i8,4e10.3e2)') "Matrix",qq,tm_fft_init,tm_fft,tm_ifft_init,tm_ifft
     tm_fft_init_tot = tm_fft_init_tot + tm_fft_init
     tm_fft_tot = tm_fft_tot +tm_fft
     tm_ifft_init_tot = tm_ifft_init_tot +tm_ifft_init
     tm_ifft_tot = tm_ifft_tot +tm_ifft
     tm_fft_init_max = max(tm_fft_init_max,tm_fft_init)
     tm_ifft_init_max = max(tm_ifft_init_max,tm_ifft_init)


  end do

  i = 1
  !$  i = omp_get_max_threads() 
  tm1 = real(nq*n3,kind=wp)
  tm2 = real(nq,kind=wp)

  if (my_id .eq. 0) then
     write(*,'(a8,7i8,8e10.3e2)') "Average",fftlib,no_proc,i,n1,n2,n3,nq,&
          tm_fft_init_tot/tm2,&
          tm_fft_init_max,tm_ifft_init_tot/tm2,tm_ifft_init_max,&
          tm_fft_tot,tm_fft_tot/tm1,tm_ifft_tot,tm_ifft_tot/tm1
  end if




  ! Deallocate arrays
  deallocate(A,B,C, stat=stat)
  if (stat .ne. 0) then
     write(*, '(a)') "Error deallocating arrays"
     goto 100
  end if

  goto 100

10 write(*,'(a)') "usage: ./bench2d.exe n1 n2 n3 nq fftlib"
  write(*,'(a)') " n1=positive integer : first dimension of cuboid"
  write(*,'(a)') " n2=positive integer : second dimenstion of cuboid"
  write(*,'(a)') " n3=positive integer : third dimenstion of cuboid"
  write(*,'(a)') " nq=positive integer : number of multiplying cuboids"
  write(*,'(a)') " fftlib=positive integer less than 6: FFT library to use"
  write(*,'(a)') "   fftlib=1: FFTE"
  write(*,'(a)') "   fftlib=2: FFTW"
  write(*,'(a)') "   fftlib=3: MKL"

100 continue

  call mpi_finalize(ierr)

contains

  subroutine fft_bench(comm,n1,n2,n3,C,fftlib,check,flag,tm_fft_init,tm_fft,&
       tm_ifft_init,tm_ifft)
    integer, intent(inout) :: comm ! mpi communicator
    integer, intent(in) :: n1,n2,n3 ! Array dimensions
    real (kind=wp), intent(in) :: C(n1,n2,n3) ! Input array 
    integer, intent(in) :: fftlib ! fft library to use
    !  1: FFTE
    !  2: FFTW
    !  3: MKL
    !  4: P3DFFT
    !  5: P3DFFT++
    logical, intent(in) :: check ! Additionally, perform inverse, element-wise
    ! division by Bi and compare with A
    integer, intent(out) :: flag ! 0 : all fine
    ! -1: error: check is true but A or Bi missing
    real(kind=wp), intent(out) :: tm_fft_init ! total initialisation time fft
    real(kind=wp), intent(out) :: tm_fft ! total time fft
    real(kind=wp), intent(out) :: tm_ifft_init ! total initialisation time ifft
    real(kind=wp), intent(out) :: tm_ifft ! total time ifft

    ! Local variables and arrays
    complex(kind=wp), allocatable :: Dk(:,:), work(:,:), workmkl(:),local(:)
    complex(kind=wp), allocatable :: X(:,:)
      REAL(kind=wp), POINTER :: X_IN(:,:)
      COMPLEX(kind=wp), POINTER :: X_OUT(:,:)
 

   real(kind=wp) :: nrm,tm1,tm2, t, s1
    integer :: stat, k, i, j, ntemp, nthreads,m_padded

    type(DFTI_DESCRIPTOR_DM), POINTER :: My_Desc_Handle!, My_Desc_Handle_Inv


    integer :: Status, L(2)
    integer :: strides_in(3)
    integer :: strides_out(3)

    integer ::  nx,nx_out,start_x,start_x_out,size

    integer   elementsize, rootrank
    parameter (elementsize = 16)


    flag = 0
    tm_fft_init = 0.0_wp
    tm_fft = 0.0_wp
    tm_ifft_init = 0.0_wp
    tm_ifft = 0.0_wp

    !     write(*,*) 'a',A
    !     write(*,*) 'B',B
    !     write(*,*) 'C',C

    select case (fftlib)

    case (1)
       ! FFTE
       ! Check that n1 and n2 factorise into powers of 2, 3 and 5
       do i=1,2
          if (i.eq.1) then
             ntemp = n1
          else
             ntemp = n2
          end if
          do while (mod(ntemp,2).eq.0)
             ntemp = ntemp/2
          end do
          do while (mod(ntemp,3).eq.0)
             ntemp = ntemp/3
          end do
          do while (mod(ntemp,5).eq.0)
             ntemp = ntemp/5
          end do
          if (ntemp .ne. 1) then
             flag = -4
             goto 20
          end if
       end do


       allocate(Dk(n1,n2),stat=stat)
       if (stat .ne. 0) then
          flag = -2
          goto 20
       end if
       allocate(work(n1/2+1,n2),stat=stat)
       if (stat .ne. 0) then
          flag = -2
          goto 20
       end if

       do k=1,n3
          ! Copy each slice into Dk
          do i=1,n1
             do j=1,n2
                !     write(*,*) 'c',i,j,k,C(i,j,k)
                Dk(i,j) = cmplx(C(i,j,k),kind=wp)
                !     write(*,*) i,j,k,Dk(i,j)
             end do

          end do
          if (k.eq.1) then
             !$          tm1 = omp_get_wtime()
             call DZFFT2D(Dk,n1,n2,0,work)
             !$          tm2 = omp_get_wtime()
             tm_fft_init = tm_fft_init + tm2 - tm1
          end if
          !$      tm1 = omp_get_wtime()
          call DZFFT2D(Dk,n1,n2,-1,work)
          !$      tm2 = omp_get_wtime()
          tm_fft = tm_fft + tm2 - tm1
          !        write(*,*) 'fft time=', tm2-tm1


          !    do i=1,n1/2+1
          !      do j=1,n2
          !        write(*,*) i,j,k,Dk(i,j)
          !      end do
          !    end do
          if (check) then
             if (k.eq.1) then
                !$          tm1 = omp_get_wtime()
                call ZDFFT2D(Dk,n1,n2,0,work)
                !$          tm2 = omp_get_wtime()
                tm_ifft_init = tm_ifft_init + tm2 - tm1
             end if
             !$        tm1 = omp_get_wtime()
             call ZDFFT2D(Dk,n1,n2,1,work)
             !$        tm2 = omp_get_wtime()
             tm_ifft = tm_ifft + tm2 - tm1

             if (k.eq.1) then
                nrm = 0.0_wp
             end if
             call check_error(n1,n2,C(:,:,k),Dk,nrm)

             if (k.eq.n3) then
                !nrm = sqrt(nrm)
                write(*,*) 'k, nrm^2:',k,nrm
             end if
          end if
       end do


       deallocate(Dk, stat=stat)
       if (stat .ne. 0) then
          flag = -3
          goto 20
       end if
       deallocate(work, stat=stat)
       if (stat .ne. 0) then
          flag = -3
          goto 20
       end if

    case (2)
       ! FFTW

       call fft_bench_fftw(comm,n1,n2,n3,C,check,flag,tm_fft_init,tm_fft,&
            tm_ifft_init,tm_ifft)




    case (3) ! MKL

       nthreads = 1
       !$    nthreads=omp_get_max_threads()
       call mkl_domain_set_num_threads(nthreads, MKL_DOMAIN_FFT)
       write(*,'(a14,i5)') "MKL threads=",nthreads

       allocate(X(n1,n2),stat=stat)
       if (stat .ne. 0) then
          flag = -2
          goto 20
       end if

       X(:,:) = 0.0_wp

       L(1) = n1
       L(2) = n2
       m_padded = n1/2+1


       call mpi_barrier(comm,ierr)
       !$          tm1 = omp_get_wtime()
       Status = DftiCreateDescriptorDM( comm,My_Desc_Handle, DFTI_DOUBLE,&
            DFTI_REAL, 2, L )
       !  write(*,*) 'Status1', Status
       if (status .ne. 0) then
          if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
             write(*,*) 'Error: ', DftiErrorMessage(status)
          endif
       endif


       !     3. obtain some values of configuration parameters by calls to
       !        dftigetvaluedm
       !
       status = dftigetvaluedm(my_desc_handle,cdft_local_size,size)
       if (status .ne. 0) then
          if (.not. dftierrorclass(status,dfti_no_error)) then
             write(*,*) 'error: ', dftierrormessage(status)
          endif
       endif


       status = dftigetvaluedm(my_desc_handle,cdft_local_nx,nx)
       if (status .ne. 0) then
          if (.not. dftierrorclass(status,dfti_no_error)) then
             write(*,*) 'error: ', dftierrormessage(status)
          endif
       endif


       status = dftigetvaluedm(my_desc_handle,cdft_local_x_start,start_x)
       if (status .ne. 0) then
          if (.not. dftierrorclass(status,dfti_no_error)) then
             write(*,*) 'error: ', dftierrormessage(status)
          endif
       endif


       status = dftigetvaluedm(my_desc_handle,cdft_local_out_nx,nx_out)
       if (status .ne. 0) then
          if (.not. dftierrorclass(status,dfti_no_error)) then
             write(*,*) 'error: ', dftierrormessage(status)
          endif
       endif


       status = dftigetvaluedm(my_desc_handle,cdft_local_out_x_start,start_x_out)
       if (status .ne. 0) then
          if (.not. dftierrorclass(status,dfti_no_error)) then
             write(*,*) 'error: ', dftierrormessage(status)
          endif
       endif



       allocate(local(size), workmkl(size), stat=status)
       if (stat .ne. 0) then
          flag = -2
          goto 20
       end if



      CALL C_F_POINTER (C_LOC(LOCAL), X_IN, [2*M_PADDED,NX])
      CALL C_F_POINTER (C_LOC(LOCAL), X_OUT, [M_PADDED,NX_OUT])




       !
       !     4. specify a value(s) of configuration parameters by a call(s) to
       !        dftisetvaluedm
       !
!       status = dftisetvaluedm(my_desc_handle,cdft_workspace,workmkl)
!       if (status .ne. 0) then
!          if (.not. dftierrorclass(status,dfti_no_error)) then
!             write(*,*) 'error: ', dftierrormessage(status)
!          endif
!       endif



       Status = DftiCommitDescriptorDM( My_Desc_Handle)
       !  write(*,*) 'Status5', Status

       if (status .ne. 0) then
          if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
             write(*,*) 'Error: ', DftiErrorMessage(status)
          endif

       endif
       call mpi_barrier(comm,ierr)
       !$          tm2 = omp_get_wtime()

       tm_fft_init = tm_fft_init + tm2 - tm1

       if (check) then
          call mpi_barrier(comm,ierr)
          !$      tm1 = omp_get_wtime()
          call mpi_barrier(comm,ierr)
          !$      tm2 = omp_get_wtime()
          tm_ifft_init = tm_ifft_init + tm2 - tm1
          allocate(Dk(n1,n2),stat=stat)
          if (stat .ne. 0) then
             flag = -2
             goto 20
          end if

          Dk(:,:) = 0.0_wp
          nrm = 0.0_wp


       end if


       do k=1,n3     

          ! Copy slice from C(:,:,k)
          do i=1,n1
             do j=1,nx
                t = 0.1**12
                if (C(i,j+start_x-1,k) .lt. t) then
                   s1=0.0_wp
                else
                   s1=C(i,j+start_x-1,k)
                end if
                X_in(i,j) = cmplx(s1,kind=wp)
             end do
          end do

          

          !         write(*,*) X
          call mpi_barrier(comm,ierr)
          !$      tm1 = omp_get_wtime()


          call mpi_barrier(comm,ierr)


          Status = DftiComputeForwardDM( My_Desc_Handle,local,workmkl)

          if (status .ne. 0) then
             if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
                write(*,*) 'Error: ', DftiErrorMessage(status)
             endif
          endif

          !      write(*,*) 'Status4', Status

          call mpi_barrier(comm,ierr)



          !$      tm2 = omp_get_wtime()
          !     write(*,*) 'fft time=', tm2-tm1
          tm_fft = tm_fft + tm2 - tm1

          !         write(*,*) X

          if (check) then




          call mpi_barrier(comm,ierr)


             Status = DftiComputeBackwardDM( My_Desc_Handle, local, workmkl)

             if (status .ne. 0) then
                if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
                   write(*,*) 'Error: ', DftiErrorMessage(status)
                endif
             endif

          call mpi_barrier(comm,ierr)

             !$          tm2 = omp_get_wtime()

             tm_ifft = tm_ifft + tm2 - tm1
             !        write(*,*) 'ifft time=', tm2-tm1
             !           write(*,*) X
             ! Copy slice from X to Dk
             do i=1,n1
                do j=1,nx

                   Dk(i,j+start_x-1) = cmplx(X_in(i,j)/real(n1*n2,kind=wp),kind=wp)
                end do
             end do
 

   call check_error(n1,nx,C(1:n1,start_x:start_x+nx-1,k),Dk(1:n1,start_x:start_x+nx-1),nrm)
             !      write(*,*) 'nrm',nrm,k

             if (k.eq.n3) then
                write(*,*) 'k,nrm,my_id',k,nrm,my_id

                deallocate(Dk,stat=stat)
                if (stat .ne. 0) then
                   flag = -3
                   goto 20
                end if


             end if

          end if


          if (k.eq.n3) then
             Status = DftiFreeDescriptorDM(My_Desc_Handle)

             if (status .ne. 0) then
                if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
                   write(*,*) 'Error: ', DftiErrorMessage(status)
                endif
             endif


             !         if (k.eq.n3) then
             deallocate(X,local,workmkl,stat=stat)
             if (stat .ne. 0) then
                flag = -3
                goto 20
             end if

             !         deallocate(X_2D,stat=stat)
             !         if (stat .ne. 0) then
             !           flag = -3
             !           goto 20
          end if


       end do



    end select

    call mpi_barrier(comm,ierr)
    return

20  select case (flag)
    case (-1)
       write(*,'(a)') "Error check requested  but either A or Bi missing"
       ! should never be possible to reach this error
    case (-2)
       write(*,'(a)') "Allocation error"
    case (-3)
       write(*,'(a)') "Deallocation error"
    case (-4)
       write(*,'(a)') "n1 and n2 must be factorisable into powers of 2, 3 and 5"


    end select


  end subroutine fft_bench


  subroutine fft_bench_fftw(comm,n1,n2,n3,C,check,flag,tm_fft_init,tm_fft,&
       tm_ifft_init,tm_ifft)
    integer, intent(inout) :: comm ! mpi communicator
    integer, intent(in) :: n1,n2,n3 ! Array dimensions
    real (kind=wp), intent(in) :: C(n1,n2,n3) ! Input array 
    logical, intent(in) :: check ! Additionally, perform inverse, element-wise
    ! division by Bi and compare with A
    integer, intent(out) :: flag ! 0 : all fine
    ! -1: error: check is true but A or Bi missing
    real(kind=wp), intent(out) :: tm_fft_init ! total initialisation time fft
    real(kind=wp), intent(out) :: tm_fft ! total time fft
    real(kind=wp), intent(out) :: tm_ifft_init ! total initialisation time ifft
    real(kind=wp), intent(out) :: tm_ifft ! total time ifft

    ! Local variables and arrays
    !complex(kind=wp), allocatable :: Dk(:,:)
    real(kind=wp) :: nrm,tm1,tm2, n1n2
    integer :: stat, k, i, j, nthreads
    integer :: my_id,ierr,local_n1
    integer(c_int) :: flags
    integer(c_intptr_t) :: n1_4,n2_4
    integer(kind=4) :: nthreads_4

    integer(c_intptr_t)   :: alloc_local, local_n2, local_j_offset


    !    type(C_PTR) :: plan, iplan

    !    real(C_DOUBLE), dimension(n1,n2) :: in, iout
    !    real(C_DOUBLE), dimension(n1,n2) :: out, iin

    flag = 0
    tm_fft_init = 0.0_wp
    tm_fft = 0.0_wp
    tm_ifft_init = 0.0_wp
    tm_ifft = 0.0_wp

    n1_4 = int(n1,kind=c_intptr_t)
    n2_4 = int(n2,kind=c_intptr_t)
    flags = int(0,kind=c_int) 

    nthreads = 1
    !$  nthreads=omp_get_max_threads()                                                                                
    nthreads_4 = int(nthreads,kind=ip4)

    call mpi_barrier(comm,ierr)
    !$          tm1 = omp_get_wtime() 

    stat = fftw_init_threads()
    call fftw_mpi_init()
    call fftw_plan_with_nthreads(nthreads_4)
    alloc_local = fftw_mpi_local_size_2d(n2_4,n1_4, comm, &
         local_n2, local_j_offset)

    !   plan = fftw_plan_r2r_2d(n2_4,n1_4, in,out,FFTW_R2HC,FFTW_R2HC,flags)

    call mpi_barrier(comm,ierr)
    !$          tm2 = omp_get_wtime()
    tm_fft_init = tm_fft_init + tm2 - tm1

    if (check) then
       call mpi_barrier(comm,ierr)
       !$   tm1 = omp_get_wtime()
       ! iplan = fftw_plan_r2r_2d(n2_4,n1_4, iin,iout,FFTW_HC2R,&
       !      FFTW_HC2R,flags)

       call mpi_barrier(comm,ierr)
       !$   tm2 = omp_get_wtime()
       tm_ifft_init = tm_ifft_init + tm2 - tm1

    end if

    call fft_bench_fftw_mpi(comm,n1,n2,n3,C,check,flag,tm_fft_init,tm_fft,&
         tm_ifft_init,tm_ifft,local_n2,local_j_offset)

  end subroutine fft_bench_fftw

  subroutine fft_bench_fftw_mpi(comm,n1,n2,n3,C,check,flag,tm_fft_init,tm_fft,&
       tm_ifft_init,tm_ifft,local_n2,local_j_offset)
    integer, intent(inout) :: comm ! mpi communicator
    integer, intent(in) :: n1,n2,n3 ! Array dimensions
    real (kind=wp), intent(in) :: C(n1,n2,n3) ! Input array 
    logical, intent(in) :: check ! Additionally, perform inverse, element-wise
    ! division by Bi and compare with A
    integer, intent(out) :: flag ! 0 : all fine
    ! -1: error: check is true but A or Bi missing
    real(kind=wp), intent(inout) :: tm_fft_init ! total initialisation time fft
    real(kind=wp), intent(inout) :: tm_fft ! total time fft
    real(kind=wp), intent(inout) :: tm_ifft_init ! total initialisation time ifft
    real(kind=wp), intent(inout) :: tm_ifft ! total time ifft

    integer(c_intptr_t), intent(in)   :: local_n2, local_j_offset

    ! Local variables and arrays
    complex(kind=wp), allocatable :: Dk(:,:)
    real(kind=wp) :: nrm,tm1,tm2, n1n2
    integer :: stat, k, i, j, nthreads
    integer :: my_id,ierr,local_n1
    integer(c_int) :: flags
    integer(c_intptr_t) :: n1_4,n2_4
    integer(kind=4) :: nthreads_4
    type(C_PTR) :: plan, iplan

    real(C_DOUBLE), dimension(n1,local_n2) :: in, iout
    real(C_DOUBLE), dimension(n1,local_n2) :: out, iin

    flag = 0
    !    tm_fft_init = 0.0_wp
    !    tm_fft = 0.0_wp
    !    tm_ifft_init = 0.0_wp
    !    tm_ifft = 0.0_wp

    n1_4 = int(n1,kind=c_intptr_t)
    n2_4 = int(n2,kind=c_intptr_t)
    flags = int(0,kind=c_int)
    !$          tm1 = omp_get_wtime()
    nthreads = 1
    !$  nthreads=omp_get_max_threads()                                                                                
    nthreads_4 = int(nthreads,kind=ip4)
    ! stat = fftw_init_threads()
    ! call fftw_plan_with_nthreads(nthreads_4)

    call mpi_barrier(comm,ierr)
    !$          tm1 = omp_get_wtime() 
    plan = fftw_mpi_plan_r2r_2d(n2_4,n1_4, in,out,comm,FFTW_R2HC,FFTW_R2HC,flags)
    call mpi_barrier(comm,ierr)
    !$          tm2 = omp_get_wtime()
    tm_fft_init = tm_fft_init + tm2 - tm1

    if (check) then
       call mpi_barrier(comm,ierr)
       !$   tm1 = omp_get_wtime()
       iplan = fftw_mpi_plan_r2r_2d(n2_4,n1_4, iin,iout,comm,FFTW_HC2R,&
            FFTW_HC2R,flags)
       call mpi_barrier(comm,ierr)
       !$   tm2 = omp_get_wtime()
       tm_ifft_init = tm_ifft_init + tm2 - tm1

       allocate(Dk(n1,local_n2),stat=stat)
       if (stat .ne. 0) then
          flag = -2
          goto 20
       end if



    end if


    do k=1,n3


       ! Copy each slice into in
       do i=1,n1
          do j=1,local_n2
             !     write(*,*) 'c',i,j,k,C(i,j,k)
             in(i,j) = C(i,j+local_j_offset,k)
             !     write(*,*) i,j,k,Dk(i,j)
          end do
       end do
       call mpi_barrier(comm,ierr)
       !$   tm1 = omp_get_wtime()
       call fftw_mpi_execute_r2r(plan, in, out)
       call mpi_barrier(comm,ierr)
       !$   tm2 = omp_get_wtime()
       !       write(*,*) 'fft time=', tm2-tm1
       tm_fft = tm_fft + tm2 - tm1

       if (check) then

          ! Copy out into iin
          do i=1,n1
             do j=1,local_n2
                !     write(*,*) 'c',i,j,k,C(i,j,k)
                iin(i,j) = out(i,j)
                !     write(*,*) i,j,k,Dk(i,j)
             end do

          end do
          call mpi_barrier(comm,ierr)
          !$      tm1 = omp_get_wtime()
          call fftw_mpi_execute_r2r(iplan, iin, iout)
          call mpi_barrier(comm,ierr)
          !$      tm2 = omp_get_wtime()
          !      write(*,*) 'ifft time=', tm2-tm1
          tm_ifft = tm_ifft + tm2 - tm1

          if (k.eq.1) then
             nrm = 0.0_wp
          end if

          n1n2 = real(n1*n2,kind=wp)
          !$OMP PARALLEL DO PRIVATE(j)
          do i=1,n1
             do j=1,local_n2
                Dk(i,j) = real(iout(i,j),kind=wp)/n1n2
             end do
          end do
          !$OMP END PARALLEL DO

          !          write(*,*) iout(n1/2,n2/2), Dk(n1/2,n2/2), C(n1/2,n2/2,k)
          call check_error(n1,int(local_n2),C(:,1+local_j_offset:local_n2+local_j_offset,k),Dk,nrm)

          if (k.eq.n3) then

             write(*,*) 'k, nrm^2:',k,nrm

             call mpi_barrier(comm,ierr)
             call fftw_destroy_plan(iplan)
             call mpi_barrier(comm,ierr)
             deallocate(Dk, stat=stat)
             if (stat .ne. 0) then
                flag = -3
                goto 20
             end if
          end if
       end if

       if (k.eq.n3) then
          call mpi_barrier(comm,ierr)
          call fftw_destroy_plan(plan)
          call mpi_barrier(comm,ierr)
          call fftw_mpi_cleanup()

          call fftw_cleanup_threads()

       end if
    end do

    return

20  select case (flag)
    case (-2)
       write(*,'(a)') "Allocation error"
    case (-3)
       write(*,'(a)') "Deallocation error"
    case (-4)
       write(*,'(a)') "n1 and n2 must be factorisable into powers of 2, 3 and 5"


    end select

  end subroutine fft_bench_fftw_mpi


  subroutine check_error(n1,n2,A,C,nrm)
    integer, intent(in) :: n1,n2 ! Array dimensions
    real(kind=wp), intent(in) :: A(n1,n2) ! Input array A
    complex(kind=wp), intent(in) :: C(n1,n2) ! Input array C
    real(kind=wp), intent(inout) :: nrm ! 2-norm of A-C

    ! local variables
    integer :: i,j,k
    complex(kind=wp) :: s, t

    write(*,*) 'nrm in',nrm

    !$OMP PARALLEL DO REDUCTION(+:nrm) PRIVATE(i,j,k,s,t) COLLAPSE(2)
    do i = 1,n1
       do j= 1,n2
          s= cmplx(A(i,j),kind=wp)-C(i,j)
          t = s*conjg(s)
          nrm = nrm + real(t,kind=wp)
!                    write(*,*) A(i,j), C(i,j),s,t,nrm
       end do
    end do
    !$OMP END PARALLEL DO

    write(*,*) 'nrm out',nrm



  end subroutine check_error

  integer function mkl_cdft_data_d(comm,rootrank,elementsize,dim,   &
       &                                 lengths,global,nx,start_x,local, &
       &                                 flag)

    ! include 'mpif.h'

    !     mpi related integer should have the kind mpi expect
    integer, allocatable :: counts(:),displs(:),buf(:)
    integer:: i,tmp(2),req,stat(mpi_status_size)
    integer:: comm,rootrank,nproc,nrank,mpi_err


    integer elementsize,dim,lengths(*),nx,start_x,flag,status
    complex(8) global(*),local(*)
    intent(in) comm,rootrank,elementsize,dim,lengths,nx,start_x,flag

    integer fd

    call mpi_comm_rank(comm,nrank,mpi_err)
    if (mpi_err/=mpi_success) goto 100

    if (nrank==rootrank) then
       call mpi_comm_size(comm,nproc,mpi_err)
       if (mpi_err/=mpi_success) goto 100
       allocate(counts(nproc),displs(nproc),buf(2*nproc), stat=status)
       if(status/=0) goto 100
    end if

    fd=1
    do i=1,dim-1
       fd=fd*lengths(i)
    end do

    tmp(1)=nx*fd*elementsize
    tmp(2)=(start_x-1)*fd

    call mpi_gather(tmp,2_mpi_kind,mpi_integer,buf,2_mpi_kind,        &
         &                mpi_integer,rootrank,comm,mpi_err)
    if (mpi_err/=mpi_success) goto 100

    if (nrank==rootrank) then
       counts=buf(1:2*nproc-1:2)
       displs=buf(2:2*nproc:2)
    end if

    if (flag==0) then

       call mpi_irecv(local,tmp(1),mpi_byte,rootrank,123_mpi_kind,    &
            &                  comm,req,mpi_err)
       if (mpi_err/=mpi_success) goto 100

       if (nrank==rootrank) then
          do i=0,nproc-1
             call mpi_send(global(displs(i+1)+1),counts(i+1),mpi_byte,&
                  &                       i,123_mpi_kind,comm,mpi_err)
             if (mpi_err/=mpi_success) goto 100
          end do
       end if

       call mpi_wait(req,stat,mpi_err)
       if (mpi_err/=mpi_success) goto 100
    endif

    if (flag==1) then

       call mpi_isend(local,tmp(1),mpi_byte,rootrank,222_mpi_kind,    &
            &                  comm,req,mpi_err)
       if (mpi_err/=mpi_success) goto 100

       if (nrank==rootrank) then
          do i=0,nproc-1
             call mpi_recv(global(displs(i+1)+1),counts(i+1),mpi_byte,&
                  &                       i,222_mpi_kind,comm,stat,mpi_err)
             if (mpi_err/=mpi_success) goto 100
          end do
       end if
       call mpi_wait(req,stat,mpi_err)
       if (mpi_err/=mpi_success) goto 100

    end if

    if (nrank==rootrank) deallocate(counts,displs,buf)

100 mkl_cdft_data_d=mpi_err

  end function mkl_cdft_data_d


  integer function mkl_cdft_scatterdata_d(comm,rootrank,elementsize,&
       &                                        dim,lengths,global_in,nx, &
       &                                        start_x,local_in)

    integer::  rootrank, comm
    integer :: elementsize,dim,lengths(*),nx,start_x
    complex(8) global_in(*),local_in(*)

    mkl_cdft_scatterdata_d=mkl_cdft_data_d(comm,rootrank,elementsize, &
         &                                       dim,lengths,global_in,nx,  &
         &                                       start_x,local_in,0)

  end function mkl_cdft_scatterdata_d


  integer function mkl_cdft_gatherdata_d(comm,rootrank,elementsize, &
       &                                       dim,lengths,global_out,nx, &
       &                                       start_x,local_out)

    integer:: rootrank, comm
    integer:: elementsize,dim,lengths(*),nx,start_x
    complex(8):: global_out(*),local_out(*)

    mkl_cdft_gatherdata_d=mkl_cdft_data_d(comm,rootrank,elementsize,  &
         &                                      dim,lengths,global_out,nx,  &
         &                                      start_x,local_out,1)

  end function mkl_cdft_gatherdata_d



END PROGRAM commandline
