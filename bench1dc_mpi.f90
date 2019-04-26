program commandline
  !$ use omp_lib
  use, intrinsic :: iso_c_binding
  use mkl_cdft 
  use mpi
  implicit none
  !  include '/usr/include/fftw3.f03'
  include '/opt/cray/fftw/default/ivybridge/include/fftw3-mpi.f03'

  integer, parameter :: wp = selected_real_kind(15,307)  ! double real

  integer, parameter :: ip4 = selected_int_kind(4)  ! integer*4
  integer   mpi_kind
  parameter (mpi_kind = wp)


  integer :: nargs, n1, n2, nq, fftlib ! input arguments
  integer :: stat ! allocat/deallocate stat 
  integer :: flag ! error flag
  integer :: i, j, k, qq, m ! indices
  integer :: comm, ierr, my_id, no_proc ! mpi variables
  real(kind=wp) :: xo, yo, a1, b1, r  ! used in definition of ellipsoid
  character(len=100) :: option1, option2, option3, option4 
  ! for reading inputs
  complex (kind=wp), allocatable :: a(:,:) ! array a
  complex (kind=wp), allocatable :: b(:,:,:) ! b(:,:,i) is matrix b_i
  complex (kind=wp), allocatable :: c(:,:) ! c(:,:) is matrix c_i
  real (kind=wp) :: s1,s2 ! used to define b
  complex(kind=wp) :: t1,t2 ! used to define c
  real (kind=wp) :: tm1, tm2, tm_fft_init, tm_fft, tm_ifft_init, tm_ifft

  real (kind=wp) :: tm_fft_init_tot, tm_fft_tot, tm_ifft_init_tot,&
       tm_ifft_tot, tm_fft_init_max, tm_ifft_init_max

  logical :: init, check

  call mpi_init(ierr)

  comm = mpi_comm_world

  ! find out process id and total no. processes
  call mpi_comm_rank(comm,my_id,ierr)
  call mpi_comm_size(comm,no_proc,ierr)
  write(*,*) 'my rank is ', my_id, ' out of ', no_proc, ' processes'

  !read input from the command line
  !nargs = iargc()
  nargs = command_argument_count() 
  if (nargs .ne. 4) then
     goto 10
  else
     call get_command_argument(1,option1) !grab the first command line argument
     ! and store it in temp variable 'option1'
     call get_command_argument(2,option2) !grab the 2nd command line argument
     ! and store it in temp variable 'option2'
     call get_command_argument(3,option3) !grab the 3rd command line argument
     ! and store it in temp variable 'option3'
     !  1: ffte
     !  2: fftw
     !  3: mkl
     !  4: p3dfft not included due to fussy installation scripts
     !  5: p3dfft++ not included due to poor installation instructions
     call get_command_argument(4,option4) !grab the 4th command line argument
     ! and store it in temp variable 'option4'


     read(option1,*) n1 !now convert string to integer
     read(option2,*) n2 !now convert string to integer
     read(option3,*) nq !now convert string to integer
     read(option4,*) fftlib !now convert string to integer
     if (my_id .eq. 0) then
        write(*,'(a,i8)') "variable n1 = ", n1
        write(*,'(a,i8)') "variable n2 = ", n2
        write(*,'(a,i8)') "variable nq = ", nq
        write(*,'(a,i8)') "variable fftlib = ", fftlib
     end if
     if (n1.lt.1 .or. n2.lt.1 .or. nq.lt.1 .or. fftlib .lt. 1 &
          .or. fftlib .gt. 3) then
        goto 10
     endif

  endif

  ! allocate arrays
  allocate(a(n1,n2),stat=stat)
  if (stat .ne. 0) then
     write(*, '(a)') "error allocating a"
     goto 100
  end if

  allocate(b(n1,n2,nq),stat=stat)
  if (stat .ne. 0) then
     write(*, '(a)') "error allocating b"
     deallocate(a)
     goto 100
  end if
  allocate(c(n1,n2),stat=stat)
  if (stat .ne. 0) then
     write(*, '(a)') "error allocating c"
     deallocate(a,b)
     goto 100
  end if

  ! initialise total times
  tm_fft_init_tot=0.0_wp
  tm_fft_tot=0.0_wp
  tm_ifft_init_tot=0.0_wp
  tm_ifft_tot=0.0_wp
  tm_fft_init_max=0.0_wp
  tm_ifft_init_max=0.0_wp



  ! set a

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

  a1 = 0.3*real(n1,wp)
  b1 = 0.35*real(n2,wp)
  !  c1 = real(n3,wp)/3.0_wp
  !$ tm1=omp_get_wtime()
  do i=1,n1
     do j=1,n2
        r = ((real(i,wp)-xo)/a1)**2 + ((real(j,wp)-yo)/b1)**2
        if (r .le. 1.0_wp) then
           a(i,j) = cmplx(r + 0.5_wp,-2.0_wp*r+0.5_wp,kind=wp)
        else
           a(i,j) = cmplx(0.5_wp,-1.5_wp,kind=wp)
        end if
     end do
  end do
  !$ tm2=omp_get_wtime()
  write(*,*) 'set-up a time=', tm2-tm1, "rank", my_id
  !  write(*,*) a
  ! set b


  m = 0
  !$ tm1=omp_get_wtime()

  !$omp parallel do private (j,k,qq) shared(n1,n2,nq,b) collapse(3)
  do i=1,n1
     do j=1,n2
        do qq=1,nq
           s1 = real(i*qq,kind=wp)/real(j,kind=wp)
           b(i,j,qq) = cmplx(s1,-0.5_wp*s1,kind=wp)          
        end do
     end do
  end do
  !$omp end parallel do


  !$ tm2=omp_get_wtime()
  write(*,*) 'set-up b time (no norm)=', tm2-tm1, "rank", my_id



  !  write(*,*) 'b(n1/2,n2/2,n3/2,q)', b(n1/2,n2/2,n3/3,:) 

  ! normalise norm(b(i,j,k,:),2) to equal 1 


!!$ tm1=omp_get_wtime()

  ! !$omp parallel do private(j,s1)

  !  do i=1,n1
  !    do j=1,n2
  !        s1=0.0_wp
  !        do qq=1,nq
  !          s1 = s1 + (b(i,j,qq))**2
  !        end do 
  !        s1 = s1**0.5
  !        if (s1 .ge. 0.00000001_wp) then
  !          do qq=1,nq
  !            b(i,j,qq) = b(i,j,qq)/s1
  !          end do 
  !        else
  !            b(i,j,:) = 1.0_wp/(real(nq,kind=wp)**0.5)
  !        end if
  !    end do
  !  end do
  !  !$omp end parallel do



!!$ tm2=omp_get_wtime()
  !   write(*,*) 'set-up b time=', tm2-tm1

  !   write(*,*) b

  !   write(*,*) 'b(n1/2,n2/2,n3/2,q)', b(n1/2,n2/2,n3/3,:)


  ! set init to true so that fft initilisation performed first
  init = .true.

  ! set-up each 2d slice and perform fft
  ! each slice formed in c(:,:,:) by performing element-wise multiplaction of 
  ! a with b(:,:,:,qq)
  do qq=1,nq
     !$  tm1=omp_get_wtime()
     !$omp parallel do private (j,s1,s2) collapse(2)
     do i=1,n1
        do j=1,n2
           t1 = a(i,j)
           t2 = b(i,j,qq)
           c(i,j) = t1*t2
        end do
     end do
     !$omp end parallel do

     !  write(*,*) 'qq',qq
     !  write(*,*) c



     !$ tm2=omp_get_wtime()
     write(*,*) 'set-up cq time=', tm2-tm1, "rank", my_id

     ! write(*,*) 'c(n1/2,n2/2,n3/2)', c(n1/2,n2/2,n3/3)


     ! perform fft on each 2d slice
     check=.true.
     call fft_bench(comm,n1,n2,c,fftlib,check,flag,tm_fft_init,tm_fft,&
          tm_ifft_init,tm_ifft)

     write(*,'(a10,i8,4e10.3e2)') "matrix",qq,tm_fft_init,tm_fft,tm_ifft_init,tm_ifft
     tm_fft_init_tot = tm_fft_init_tot + tm_fft_init
     tm_fft_tot = tm_fft_tot +tm_fft
     tm_ifft_init_tot = tm_ifft_init_tot +tm_ifft_init
     tm_ifft_tot = tm_ifft_tot +tm_ifft
     tm_fft_init_max = max(tm_fft_init_max,tm_fft_init)
     tm_ifft_init_max = max(tm_ifft_init_max,tm_ifft_init)


  end do
  i = 1
  !$  i = omp_get_max_threads()                                                                                     
  tm1 = real(nq*n2,kind=wp)
  tm2 = real(nq,kind=wp)
  if (my_id .eq. 0) then
     write(*,'(a8,6i8,8e10.3e2)') "average",fftlib,no_proc,i,n1,n2,nq,&
          tm_fft_init_tot/tm2,&
          tm_fft_init_max,tm_ifft_init_tot/tm2,tm_ifft_init_max,&
          tm_fft_tot,tm_fft_tot/tm1,tm_ifft_tot,tm_ifft_tot/tm1
  end if




  ! deallocate arrays
  deallocate(a,b,c, stat=stat)
  if (stat .ne. 0) then
     write(*, '(a)') "error deallocating arrays"
     goto 100
  end if

  goto 100

10 write(*,'(a)') "usage: ./bench2d.exe n1 n2 n3 nq fftlib"
  write(*,'(a)') " n1=positive integer : first dimension of cuboid"
  write(*,'(a)') " n2=positive integer : second dimension of cuboid"
  write(*,'(a)') " nq=positive integer : number of multiplying cuboids"
  write(*,'(a)') " fftlib=positive integer less than 4: fft library to use"
  write(*,'(a)') "   fftlib=1: ffte"
  write(*,'(a)') "   fftlib=2: fftw"
  write(*,'(a)') "   fftlib=3: mkl"
  ! write(*,'(a)') "   fftlib=4: p3dfft"
  ! write(*,'(a)') "   fftlib=5: p3dfft++"


100 continue
  call mpi_finalize(ierr)

contains

  subroutine fft_bench(comm,n1,n2,c,fftlib,check,flag,tm_fft_init,tm_fft,&
       tm_ifft_init,tm_ifft)
    integer, intent(inout) :: comm ! mpi communicator
    integer, intent(in) :: n1,n2 ! array dimensions
    complex (kind=wp), intent(in) :: c(n1,n2) ! input array 
    integer, intent(in) :: fftlib ! fft library to use
    !  1: ffte
    !  2: fftw
    !  3: mkl
    !  4: p3dfft
    !  5: p3dfft++
    logical, intent(in) :: check ! additionally, perform inverse, element-wise
    ! division by bi and compare with a
    integer, intent(out) :: flag ! 0 : all fine
    ! -1: error: check is true but a or bi missing
    real(kind=wp), intent(out) :: tm_fft_init ! total initialisation time fft
    real(kind=wp), intent(out) :: tm_fft ! total time fft
    real(kind=wp), intent(out) :: tm_ifft_init ! total initialisation time ifft
    real(kind=wp), intent(out) :: tm_ifft ! total time ifft

    ! local variables and arrays
    complex(kind=wp), allocatable :: dk(:,:), dkffte(:,:,:), work(:,:,:), work1(:,:)
    complex(kind=wp), allocatable :: x(:), workmkl(:)

    complex(kind=wp), allocatable :: local(:)
    real(kind=wp) :: nrm,tm1,tm2, t
    complex(kind=wp) :: t1,s
    integer :: stat, k, i, ntemp, nthreads, nn, n1_local
    integer :: my_id, npu,ierr ! mpi variables
    integer   elementsize, rootrank
    parameter (elementsize = 16)

    integer   nx,nx_out,start_x,start_x_out,size

    type(dfti_descriptor_dm), pointer :: my_desc_handle!, my_desc_handle_inv
    integer :: status

    call mpi_comm_rank(comm, my_id, ierr)
    call mpi_comm_size(comm, npu, ierr)
    write(*,*) 'fft_bench, rank ', my_id

    flag = 0

    tm_fft_init = 0.0_wp
    tm_fft = 0.0_wp
    tm_ifft_init = 0.0_wp
    tm_ifft = 0.0_wp



    !   write(*,*) 'benchfft'
    !write(*,*) a

    !   write(*,*) b

    !   write(*,*) c



    select case (fftlib)

    case (1)
       ! ffte
       ! check that n1 factorises into powers of 2, 3 and 5


       ntemp = n1
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

       nn = ceiling(real(n1)/real(npu))
       n1_local = min(n1,(my_id+1)*nn) - my_id*nn
       call mpi_barrier(comm,ierr)
       write(*,*) 'n1_local,n1,npu,nn,my_id', n1_local, n1,npu, nn,my_id

       ntemp = n1_local
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



       allocate(dkffte(1,1,n1_local),stat=stat)
       if (stat .ne. 0) then
          flag = -2
          goto 20
       end if
       allocate(work(1,1,n1_local),stat=stat)
       if (stat .ne. 0) then
          flag = -2
          goto 20
       end if
       allocate(work1(n1_local,1),stat=stat)
       if (stat .ne. 0) then
          flag = -2
          goto 20
       end if


       do k=1,n2
          ! copy each slice into dk



          !$omp parallel do 
          do i=1,n1_local
             !     write(*,*) 'c',i,j,k,c(i,j,k)
             dkffte(1,1,i) = c(i+my_id*nn,k)
             write(*,*) my_id,i,dkffte(1,1,i)!, c(i+my_id*nn,k)


          end do
          !$omp end parallel do
          write(*,*) 'dk filled', my_id
          call mpi_barrier(comm,ierr)

          if (k.eq.1) then
             !$          tm1 = omp_get_wtime()
             call pzfft3d(dkffte,work,1,1,n1,comm,npu,0)
             !$          tm2 = omp_get_wtime()
             tm_fft_init = tm_fft_init + tm2 - tm1

          end if
          !$   tm1 = omp_get_wtime()
          call pzfft3d(dkffte,work,1,1,n1,comm,npu,-1)
          !$   tm2 = omp_get_wtime()
          tm_fft = tm_fft + tm2 - tm1
          !        write(*,*) 'fft time=', tm2-tm1


          do i=1,n1_local
             !      do j=1,n2
             write(*,*) 'out',my_id,i,work(1,1,i)
             !      end do
          end do
          if (check) then
             if (k.eq.1) then
                !$          tm1 = omp_get_wtime()
                !  call zfft1d(dk,work,work1,n1,comm,my_id,npu,0)
                !$          tm2 = omp_get_wtime()
                tm_ifft_init = tm_ifft_init + tm2 - tm1
             end if
             !$        tm1 = omp_get_wtime()
             call pzfft3d(work,dkffte,1,1,n1,comm,npu,1)
             !$        tm2 = omp_get_wtime()
             tm_ifft = tm_ifft + tm2 - tm1

             do i=1,n1_local
                !      do j=1,n2                                                                                                                                                           
                !write(*,*) 'iout',my_id,i,dk(1,i), c(i+my_id*nn,k)
                ! if (npu .eq. 1) then
                work1(i,1) = dkffte(1,1,i)
                ! else
                !   work1(i,1) = dk(1,i)*n1
                ! end if

                write(*,*) 'iout',my_id,i,dkffte(1,1,i), c(i+my_id*nn,k)


                !      end do                                                                                                                                                              
             end do


             if (k.eq.1) then
                nrm = 0.0_wp
             end if
             call check_error(n1_local,c(1+my_id*nn:my_id*nn+n1_local,k),work1,nrm)
             write(*,*) 'nrm', nrm, my_id
             if (k.eq.n2) then
                !nrm = sqrt(nrm)
                write(*,*) 'rank, k, nrm^2:',my_id,k,nrm
             end if
          end if

       end do

       deallocate(dkffte, stat=stat)
       if (stat .ne. 0) then
          flag = -3
          goto 20
       end if
       deallocate(work,work1, stat=stat)
       if (stat .ne. 0) then
          flag = -3
          goto 20
       end if

    case (2)
       ! fftw

       call fft_bench_fftw(comm,n1,n2,c,check,flag,tm_fft_init,tm_fft,&
            tm_ifft_init,tm_ifft)


    case (3) ! mkl
       nthreads = 1
       !$    nthreads=omp_get_max_threads()
       call mkl_set_num_threads(nthreads)
       write(*,'(a14,i5)') "mkl threads=",nthreads



       do k=1,n2     
          if (k.eq.1) then
             !allocate(x_2d(2*(n1/2+1),n2),stat=stat)
             !if (stat .ne. 0) then
             !  flag = -2
             !  goto 20
             !end if
             allocate(x(n1),stat=stat)
             if (stat .ne. 0) then
                flag = -2
                goto 20
             end if
             !             equivalence(x_2d,x)
             !          l(1) = n1
             !          l(2) = 1
             call mpi_barrier(comm,ierr)

             !$          tm1 = omp_get_wtime()
             status = dfticreatedescriptordm(comm, my_desc_handle, dfti_double,&
                  dfti_complex, 1, n1 )
             !  write(*,*) 'status1', status
             if (status .ne. 0) then
                if (.not. dftierrorclass(status,dfti_no_error)) then
                   write(*,*) 'error: ', dftierrormessage(status)
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




             !
             !     4. specify a value(s) of configuration parameters by a call(s) to
             !        dftisetvaluedm
             !
             status = dftisetvaluedm(my_desc_handle,cdft_workspace,workmkl)
             if (status .ne. 0) then
                if (.not. dftierrorclass(status,dfti_no_error)) then
                   write(*,*) 'error: ', dftierrormessage(status)
                endif
             endif



             status = dfticommitdescriptordm( my_desc_handle)
             !  write(*,*) 'status5', status

             if (status .ne. 0) then
                if (.not. dftierrorclass(status,dfti_no_error)) then
                   write(*,*) 'error: ', dftierrormessage(status)
                endif
             endif


             call mpi_barrier(comm,ierr)
             !$          tm2 = omp_get_wtime()
             tm_fft_init = tm_fft_init + tm2 - tm1


          end if

          ! copy slice from c(:,:,k)
          do i=1,n1
             t = 0.1**12
             t1 = c(i,k)
             if (real(t1*conjg(t1),kind=wp).lt. t) then
                s = 0.0_wp
             else
                s = t1
             end if
             x(i) = s
          end do

          call mpi_barrier(comm,ierr)
          !$      tm1 = omp_get_wtime()
          !  write(*,*) x

          rootrank=0
          status = mkl_cdft_scatterdata_d(comm,rootrank,elementsize,1,      &
               [n1],x,nx,start_x,local)
          if (status .ne. 0) then
             if (.not. dftierrorclass(status,dfti_no_error)) then
                write(*,*) 'error: ', dftierrormessage(status)
             endif
          endif

          status = dfticomputeforwarddm( my_desc_handle, local )

          ! write(*,*) x

          !  write(*,*) 'stat',status
          if (status .ne. 0) then
             if (.not. dftierrorclass(status,dfti_no_error)) then
                write(*,*) 'error: ', dftierrormessage(status)
             endif
          endif

          !
          !     gather data among processors
          !
          status = mkl_cdft_gatherdata_d(comm,rootrank,elementsize,         &
               1,[n1],x,nx,start_x,local)


          call mpi_barrier(comm,ierr)
          !$      tm2 = omp_get_wtime()
          tm_fft = tm_fft + tm2 - tm1
          !        write(*,*) 'fft time=', tm2-tm1

          if (check) then

             if (k.eq.1) then
                call mpi_barrier(comm,ierr)

                !$      tm1 = omp_get_wtime()
                !            status = dfticreatedescriptor( my_desc_handle_inv, dfti_double,&
                !              dfti_complex, 1, l(1) )
                !          if (status .ne. 0) then
                !            if (.not. dftierrorclass(status,dfti_no_error)) then
                !                write(*,*) 'error: ', dftierrormessage(status)
                !            endif
                !endif


                !            status = dfticommitdescriptor( my_desc_handle_inv)
                !          if (status .ne. 0) then
                !            if (.not. dftierrorclass(status,dfti_no_error)) then
                !                write(*,*) 'error: ', dftierrormessage(status)
                !            endif
                !endif
                call mpi_barrier(comm,ierr)

                !$      tm2 = omp_get_wtime()
                tm_ifft_init = tm_ifft_init + tm2 - tm1

                allocate(dk(n1,1),stat=stat)
                if (stat .ne. 0) then
                   flag = -2
                   goto 20
                end if
                dk(:,:) = 0.0_wp
                nrm = 0.0_wp

             end if

             call mpi_barrier(comm,ierr)
             !$          tm1 = omp_get_wtime()

             status = mkl_cdft_scatterdata_d(comm,rootrank,elementsize,1,      &
                  [n1],x,nx,start_x,local)
             if (status .ne. 0) then
                if (.not. dftierrorclass(status,dfti_no_error)) then
                   write(*,*) 'error: ', dftierrormessage(status)
                endif
             endif


             status = dfticomputebackwarddm( my_desc_handle, local )
             !                                                                                                                                                                                    
             !     gather data among processors                                                                                                                                                   
             !                                                                                                                                                                                    
             status = mkl_cdft_gatherdata_d(comm,rootrank,elementsize,         &
                  1,[n1],x,nx,start_x,local)


             call mpi_barrier(comm,ierr)

             !$          tm2 = omp_get_wtime()
             !            write(*,*) 'ifft time=', tm2-tm1

             tm_ifft = tm_ifft + tm2 - tm1

             !    write(*,*) x



             ! copy slice from x to dk
             s = real(n1,kind=wp)
             do i=1,n1

                dk(i,1) = x(i)/s
             end do
             ! write(*,*) dk
             ! write(*,*) c(:,k)
             if (my_id .eq. rootrank) then 
                call check_error(n1,c(:,k),dk,nrm)
             end if
             !       write(*,*) 'nrm',nrm,k

             if (k.eq.n2 .and. my_id .eq. rootrank) then

                write(*,*) 'k,nrm',k,nrm

                !  status = dftifreedescriptor(my_desc_handle_inv)

             end if

          end if

          if (k.eq.n2) then
             status = dftifreedescriptordm(my_desc_handle)

             deallocate(x,dk,local,workmkl,stat=stat)
             if (stat .ne. 0) then
                flag = -3
                goto 20
             end if

             !         deallocate(x_2d,stat=stat)
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
       write(*,'(a)') "error check requested  but either a or bi missing"
       ! should never be possible to reach this error
    case (-2)
       write(*,'(a)') "allocation error"
    case (-3)
       write(*,'(a)') "deallocation error"
    case (-4)
       write(*,'(a)') "n1 and n1_local must be factorisable into powers of 2, 3 and 5"
       !     call mpi_abort(comm,flag)

    end select
    call mpi_abort(comm,flag,ierr)


  end subroutine fft_bench

  subroutine fft_bench_fftw(comm,n1,n2,c,check,flag,tm_fft_init,tm_fft,&
       tm_ifft_init,tm_ifft)
    integer, intent(inout) :: comm ! mpi communicator  
    integer, intent(in) :: n1,n2 ! array dimensions
    complex (kind=wp), intent(in) :: c(n1,n2) ! input array 
    logical, intent(in) :: check ! additionally, perform inverse, element-wise
    ! division by bi and compare with a
    integer, intent(out) :: flag ! 0 : all fine
    ! -1: error: check is true but a or bi missing
    real(kind=wp), intent(out) :: tm_fft_init ! total initialisation time fft
    real(kind=wp), intent(out) :: tm_fft !   total time fft
    real(kind=wp), intent(out) :: tm_ifft_init ! total initialisation time ifft
    real(kind=wp), intent(out) :: tm_ifft ! total time ifft

    ! local variables and arrays
    complex(kind=wp), allocatable :: dk(:,:)
    real(kind=wp) :: nrm,tm1,tm2
    integer :: stat, k, i,nthreads
    integer(c_int) :: flags
    integer(c_intptr_t) :: n1_4
    integer(kind=4) :: nthreads_4
    integer :: my_id,ierr,local_n1

    type(c_ptr) :: plan, iplan, cin, ciout, cout, ciin

    integer(c_intptr_t)   :: alloc_local, local_ni, local_i_start, &
         local_no, local_o_start
    integer(c_intptr_t)   :: ialloc_local, ilocal_ni, ilocal_i_start, &
         ilocal_no, ilocal_o_start


    complex(c_double_complex),pointer :: in(:), iout(:)
    complex(c_double_complex),pointer :: out(:), iin(:)

    call mpi_comm_rank(comm, my_id, ierr );
    flag = 0

    tm_fft_init = 0.0_wp
    tm_fft = 0.0_wp
    tm_ifft_init = 0.0_wp
    tm_ifft = 0.0_wp
    n1_4 = int(n1,kind=c_intptr_t)
    flags = int(0,kind=c_int)
    nthreads=1
    !$  nthreads = omp_get_max_threads()
    if (my_id .eq.0) then
       write(*,*) 'nthreads',nthreads
    end if
    nthreads_4 = int(nthreads,kind=ip4)

    call mpi_barrier(comm,ierr)
    !$    tm1 = omp_get_wtime()
    !stat=fftw_init_threads()
    ! if (stat .eq. 0) then
    !    write(*,*) 'fftw_init_threads stat', stat        
    ! end if

    stat = fftw_init_threads()
    call fftw_mpi_init()
    call fftw_plan_with_nthreads(nthreads_4)

    !   get local data size and allocate
    alloc_local = fftw_mpi_local_size_1d(n1_4, comm, &
         fftw_forward,flags, local_ni, local_i_start, local_no, local_o_start               )
    cin = fftw_alloc_complex(alloc_local)
    cout = fftw_alloc_complex(alloc_local)

    call c_f_pointer(cin, in, [local_ni])
    call c_f_pointer(cout, out, [local_no])
    ! call fftw_plan_with_nthreads(nthreads_4)

    plan = fftw_mpi_plan_dft_1d(n1_4, in,out,comm,fftw_forward,flags)
    call mpi_barrier(comm,ierr)

    !$    tm2 = omp_get_wtime()
    tm_fft_init = tm_fft_init + tm2 - tm1

    !    write(*,'(a5,i5,a10,i8,a10,i8)') 'rank', my_id,'local_ni',local_ni,'local_i_start', local_i_start
    !    write(*,'(a5,i5,a10,i8,a10,i8)') 'rank', my_id,'local_no', local_no, 'local_o_start', local_o_start

    if (check) then
       call mpi_barrier(comm,ierr)

       !$   tm1 = omp_get_wtime()
       !   get local data size and allocate                                                                                                                                                 
       ialloc_local = fftw_mpi_local_size_1d(n1_4, comm, &
            fftw_backward,flags, ilocal_ni, ilocal_i_start, ilocal_no, ilocal_o_start               )
       ciin = fftw_alloc_complex(ialloc_local)
       ciout = fftw_alloc_complex(ialloc_local)

       call c_f_pointer(ciin, iin, [ilocal_ni])
       call c_f_pointer(ciout, iout, [ilocal_no])


       iplan = fftw_mpi_plan_dft_1d(n1_4, iin,iout,comm,fftw_backward,&
            flags)
       call mpi_barrier(comm,ierr)

       !$   tm2 = omp_get_wtime()
       tm_ifft_init = tm_ifft_init + tm2 - tm1

       !       write(*,'(a5,i5,a10,i8,a10,i8)') 'rank', my_id,'ilocal_ni',ilocal_ni,'ilocal_i_start', ilocal_i_start
       !       write(*,'(a5,i5,a10,i8,a10,i8)') 'rank', my_id,'ilocal_no', ilocal_no, 'ilocal_o_start', ilocal_o_start

       allocate(dk(n1,1),stat=stat)
       if (stat .ne. 0) then
          flag = -2
          goto 20
       end if
    end if

    do k=1,n2

       ! copy each slice into in
       do i=1,local_ni
          !     write(*,*) 'c',i,j,k,c(i,j,k)
          in(i) = c(i+local_i_start,k)
          !     write(*,*) i,j,k,dk(i,j)

       end do
       call mpi_barrier(comm,ierr)

       !$   tm1 = omp_get_wtime()
       call fftw_mpi_execute_dft(plan, in, out)
       call mpi_barrier(comm,ierr)

       !$   tm2 = omp_get_wtime()
       tm_fft = tm_fft + tm2 - tm1

       !        write(*,*) 'fft time=', tm2-tm1

       if (check) then

          ! copy out into iin
          do i=1,local_no
             !             write(*,*) 'c',i,in(i),out(i)
             iin(i) = out(i)
             !     write(*,*) i,j,k,dk(i,j)

          end do
          call mpi_barrier(comm,ierr)

          !$      tm1 = omp_get_wtime()
          call fftw_mpi_execute_dft(iplan, iin, iout)
          call mpi_barrier(comm,ierr)
          !$      tm2 = omp_get_wtime()




          !        write(*,*) 'ifft time=', tm2-tm1
          tm_ifft = tm_ifft + tm2 - tm1

          if (k.eq.1) then
             nrm = 0.0_wp
          end if

          !$omp parallel do private(i)
          do i=1,ilocal_no
             dk(i,1) = cmplx(iout(i),kind=wp)/real(n1,kind=wp)
          end do
          !$omp end parallel do

          !          write(*,*) iout(n1/2,n2/2), dk(n1/2,n2/2), c(n1/2,n2/2,k)
          local_n1 = int(ilocal_no)
          call check_error(local_n1,c(ilocal_o_start+1:ilocal_o_start+ilocal_no,k),&
               dk(1:local_n1,1),nrm)


          if (k.eq.n2) then

             write(*,*) 'k, nrm^2:',k,nrm, my_id

             call fftw_destroy_plan(iplan)
             deallocate(dk, stat=stat)
             if (stat .ne. 0) then
                flag = -3
                goto 20
             end if

          end if


       end if



       if (k.eq.n2) then
          call fftw_destroy_plan(plan)

          call mpi_barrier(comm,ierr)

          !         call fftw_free(cin)
          !         call fftw_free(cout)
          !call mpi_barrier(comm,ierr)
          call fftw_mpi_cleanup()

          call fftw_cleanup_threads()

       end if



    end do



    return

20  select case (flag)
    case (-1)
       write(*,'(a)') "error check requested  but either a or bi missing"
       ! should never be possible to reach this error
    case (-2)
       write(*,'(a)') "allocation error"
    case (-3)
       write(*,'(a)') "deallocation error"
    case (-4)
       write(*,'(a)') "n1 and n2 must be factorisable into powers of 2, 3 and 5"


    end select


  end subroutine fft_bench_fftw


  subroutine check_error(n1,a,c,nrm)
    integer, intent(in) :: n1 ! array dimensions
    complex(kind=wp), intent(in) :: a(n1,1) ! input array a
    complex(kind=wp), intent(in) :: c(n1,1) ! input array c
    real(kind=wp), intent(inout) :: nrm ! 2-norm of a-c

    ! local variables
    integer :: i
    complex(kind=wp) :: s, t
    !   write(*,*) 'nrm',nrm

    !$omp parallel do private(s,t) reduction(+:nrm)
    do i = 1,n1
       s= a(i,1)-c(i,1)
       t = s*conjg(s)
       !                 write(*,*) 's,t',a(i,1), c(i,1),s,t, real(t,kind=wp),my_id 
       nrm = nrm + real(t,kind=wp)
    end do
    !$omp end parallel do
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

end program commandline
