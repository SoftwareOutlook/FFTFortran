PROGRAM commandline
!$ use omp_lib
  use, intrinsic :: iso_c_binding
  use mkl_dfti 
  implicit none
  include '/opt/cray/fftw/default/ivybridge/include/fftw3-mpi.f03'

  Integer, Parameter :: wp = Selected_Real_Kind(15,307)  ! double real

  Integer, Parameter :: ip4 = Selected_Int_Kind(4)  ! integer*4



  INTEGER :: nargs, n1, n2, nq, fftlib ! Input arguments
  integer :: stat ! Allocat/deallocate stat 
  integer :: flag ! Error flag
  integer :: i, j, k, qq, m ! indices
  real(kind=wp) :: xo, yo, a1, b1, r  ! Used in definition of ellipsoid
  CHARACTER(LEN=100) :: option1, option2, option3, option4 
  ! For reading inputs
  real (kind=wp), allocatable :: A(:,:) ! Array A
  real (kind=wp), allocatable :: B(:,:,:) ! B(:,:,i) is cube B_i
  real (kind=wp), allocatable :: C(:,:) ! C(:,:,:) is cube C_i
  real (kind=wp) :: s1,s2 ! used to define B
  real (kind=wp) :: tm1, tm2, tm_fft_init, tm_fft, tm_ifft_init, tm_ifft

  real (kind=wp) :: tm_fft_init_tot, tm_fft_tot, tm_ifft_init_tot, &
                   tm_ifft_tot, tm_ifft_init_max, tm_fft_init_max

  logical :: init, check

  !Read input from the command line
  !nargs = IARGC()
  nargs = COMMAND_ARGUMENT_COUNT() 
  IF (nargs .ne. 4) THEN
    goto 10
  ELSE
    CALL GET_COMMAND_ARGUMENT(1,option1) !Grab the first command line argument
    ! and store it in temp variable 'option1'
    CALL GET_COMMAND_ARGUMENT(2,option2) !Grab the 2nd command line argument
    ! and store it in temp variable 'option2'
    CALL GET_COMMAND_ARGUMENT(3,option3) !Grab the 3rd command line argument
    ! and store it in temp variable 'option3'
         !  1: FFTE NOT included due to bugs in real-complex fft
         !  2: FFTW
         !  3: MKL
         !  4: P3DFFT
         !  5: P3DFFT++ NOT included due to poor installation instructions
    CALL GET_COMMAND_ARGUMENT(4,option4) !Grab the 4th command line argument
    ! and store it in temp variable 'option4'


    read(option1,*) n1 !Now convert string to integer
    read(option2,*) n2 !Now convert string to integer
    read(option3,*) nq !Now convert string to integer
    read(option4,*) fftlib !Now convert string to integer
    write(*,'(a,i8)') "Variable n1 = ", n1
    write(*,'(a,i8)') "Variable n2 = ", n2
    write(*,'(a,i8)') "Variable nq = ", nq
    write(*,'(a,i8)') "Variable fftlib = ", fftlib

    if (n1.lt.1 .or. n2.lt.1 .or. nq.lt.1 .or. fftlib .lt. 2 &
        .or. fftlib .gt. 4) then
      goto 10
    endif

  ENDIF

  ! Allocate arrays
  allocate(A(n1,n2),stat=stat)
  if (stat .ne. 0) then
      write(*, '(a)') "Error allocating A"
      goto 100
  end if

  allocate(B(n1,n2,nq),stat=stat)
  if (stat .ne. 0) then
      write(*, '(a)') "Error allocating B"
      deallocate(A)
      goto 100
  end if
  allocate(C(n1,n2),stat=stat)
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
  tm_fft_init_max = 0.0_wp
  tm_ifft_init_max = 0.0_wp


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


  a1 = 0.3*real(n1,wp)
  b1 = 0.35*real(n2,wp)

!$ tm1=omp_get_wtime()
  do i=1,n1
    do j=1,n2
         r = ((real(i,wp)-xo)/a1)**2 + ((real(j,wp)-yo)/b1)**2
         if (r .le. 1.0_wp) then
            A(i,j) = r + 0.5_wp
         else
            A(i,j) = 0.5_wp
         end if
       !  write(*,*) i,j,k, A(i,j,k)
    end do
  end do
!$ tm2=omp_get_wtime()
   write(*,*) 'Set-up A time=', tm2-tm1
 !  write(*,*) A
  ! Set B

  
  m = 0
!$ tm1=omp_get_wtime()

!$OMP PARALLEL DO PRIVATE (j,k,qq) SHARED(n1,n2,nq,B) COLLAPSE(3)
    do i=1,n1
      do j=1,n2
        do qq=1,nq
          B(i,j,qq) = real(i*qq,kind=wp)/real(j,kind=wp)
        end do
      end do
    end do
!$OMP END PARALLEL DO


!$ tm2=omp_get_wtime()
   write(*,*) 'Set-up B time (no norm)=', tm2-tm1

  ! Set init to true so that fft initilisation performed first
  init = .true.

  ! Set-up each 2D slice and perform FFT
  ! Each slice formed in C(:,:,:) by performing element-wise multiplaction of 
  ! A with B(:,:,:,qq)
  do qq=1,nq
!$  tm1=omp_get_wtime()
!$OMP PARALLEL DO PRIVATE (j,s1,s2) COLLAPSE(2)
    do i=1,n1
      do j=1,n2
          s1 = A(i,j)
          s2 = B(i,j,qq)
          C(i,j) = s1*s2
      end do
    end do
!$OMP END PARALLEL DO

!$ tm2=omp_get_wtime()
   write(*,*) 'Set-up Cq time=', tm2-tm1
    
  ! Perform FFT on each 2D slice
    check=.true.
    call fft_bench(n1,n2,C,fftlib,check,flag,tm_fft_init,tm_fft,&
       tm_ifft_init,tm_ifft)

    write(*,'(a10,i8,4e10.3e2)') "Matrix",qq,tm_fft_init,tm_fft,tm_ifft_init,tm_ifft
    tm_fft_init_tot = tm_fft_init_tot + tm_fft_init
    tm_fft_tot = tm_fft_tot +tm_fft
    tm_ifft_init_tot = tm_ifft_init_tot +tm_ifft_init
    tm_ifft_tot = tm_ifft_tot +tm_ifft
    tm_fft_init_max = max(tm_fft_init_max,tm_fft_init)
    tm_ifft_init_max = max(tm_ifft_init_max,tm_ifft_init)


  end do
    j = 1
    i = 1
!$  i = omp_get_max_threads()   
    write(*,*) 'i',i
    tm1 = real(nq*n2,kind=wp)
    tm2 = real(nq,kind=wp)
    write(*,'(a8,6i8,8e10.3e2)') "Average",fftlib,j,i,n1,n2,nq,&
       tm_fft_init_tot/tm2,&
       tm_fft_init_max,tm_ifft_init_tot/tm2,tm_ifft_init_max,&
       tm_fft_tot,tm_fft_tot/tm1,tm_ifft_tot,tm_ifft_tot/tm1


  ! Deallocate arrays
  deallocate(A,B,C, stat=stat)
  if (stat .ne. 0) then
      write(*, '(a)') "Error deallocating arrays"
      goto 100
  end if

  goto 100

10  write(*,'(a)') "usage: ./bench2d.exe n1 n2 n3 nq fftlib"
    write(*,'(a)') " n1=positive integer : first dimension of cuboid"
    write(*,'(a)') " n2=positive integer : second dimension of cuboid"
    write(*,'(a)') " nq=positive integer : number of multiplying cuboids"
    write(*,'(a)') " fftlib=positive integer less than 4: FFT library to use"
    write(*,'(a)') "   fftlib=1: FFTE"
    write(*,'(a)') "   fftlib=2: FFTW"
    write(*,'(a)') "   fftlib=3: MKL"
   ! write(*,'(a)') "   fftlib=4: P3DFFT"
   ! write(*,'(a)') "   fftlib=5: P3DFFT++"


100 continue

 contains
   
  subroutine fft_bench(n1,n2,C,fftlib,check,flag,tm_fft_init,tm_fft,&
       tm_ifft_init,tm_ifft)
    integer, intent(in) :: n1,n2 ! Array dimensions
    real (kind=wp), intent(in) :: C(n1,n2) ! Input array 
    integer, intent(in) :: fftlib ! fft library to use
         !  1: FFTE NOTE real->complex ffte has bugs and gives wrong results
         !  2: FFTW
         !  3: MKL
         !  4: P3DFFT NOT include
         !  5: P3DFFT++ NOT included
    logical, intent(in) :: check ! Additionally, perform inverse, element-wise
                                 ! division by Bi and compare with A
    integer, intent(out) :: flag ! 0 : all fine
                                 ! -1: error: check is true but A or Bi missing
    real(kind=wp), intent(out) :: tm_fft_init ! total initialisation time fft
    real(kind=wp), intent(out) :: tm_fft ! total time fft
    real(kind=wp), intent(out) :: tm_ifft_init ! total initialisation time ifft
    real(kind=wp), intent(out) :: tm_ifft ! total time ifft

    ! Local variables and arrays
    complex(kind=wp), allocatable :: Dk(:,:), work(:,:)
    real(kind=wp), allocatable :: X(:)
    real(kind=wp) :: nrm,tm1,tm2,s,t
    integer :: stat, k, i, j, ntemp, nthreads

    type(DFTI_DESCRIPTOR), POINTER :: My_Desc_Handle, My_Desc_Handle_Inv
    integer :: Status, L(2)
    integer :: strides_in(3)
    integer :: strides_out(3)


    flag = 0

    tm_fft_init = 0.0_wp
    tm_fft = 0.0_wp
    tm_ifft_init = 0.0_wp
    tm_ifft = 0.0_wp


    
!   write(*,*) 'benchfft'
   !write(*,*) A

!   write(*,*) B

!   write(*,*) C

 

    select case (fftlib)

    case (1)
      ! FFTE
      ! Check that n1 factorises into powers of 2, 3 and 5
    
    
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
    


      allocate(Dk(n1,1),stat=stat)
      if (stat .ne. 0) then
       flag = -2
       goto 20
      end if
      allocate(work(n1/2+1,1),stat=stat)
      if (stat .ne. 0) then
       flag = -2
       goto 20
      end if


      do k=1,n2
        ! Copy each slice into Dk
!$OMP PARALLEL DO 
        do i=1,n1
       !     write(*,*) 'c',i,j,k,C(i,j,k)
            Dk(i,1) = cmplx(C(i,k),kind=wp)
       !     write(*,*) i,j,k,Dk(i,j)
       

       end do
!$OMP END PARALLEL DO

        if (k.eq.1) then
!$          tm1 = omp_get_wtime()
          call DZFFT2D(Dk,n1,1,0,work)
!$          tm2 = omp_get_wtime()
            tm_fft_init = tm_fft_init + tm2 - tm1

        end if 
!$   tm1 = omp_get_wtime()
        call DZFFT2D(Dk,n1,1,-1,work)
!$   tm2 = omp_get_wtime()
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
            call ZDFFT2D(Dk,n1,1,0,work)
!$          tm2 = omp_get_wtime()
            tm_ifft_init = tm_ifft_init + tm2 - tm1
          end if
!$        tm1 = omp_get_wtime()
          call ZDFFT2D(Dk,n1,1,1,work)
!$        tm2 = omp_get_wtime()
          tm_ifft = tm_ifft + tm2 - tm1

          if (k.eq.1) then
            nrm = 0.0_wp
          end if
          call check_error(n1,C(:,k),Dk,nrm)

          if (k.eq.n2) then
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

      call fft_bench_fftw(n1,n2,C,check,flag,tm_fft_init,tm_fft,&
       tm_ifft_init,tm_ifft)


    case (3) ! MKL

    nthreads = 1
!$    nthreads=omp_get_max_threads()  
     call mkl_set_num_threads(nthreads)      
     write(*,'(a14,i5)') "MKL threads=",nthreads
     do k=1,n2     
        if (k.eq.1) then
          !allocate(X_2D(2*(n1/2+1),n2),stat=stat)
          !if (stat .ne. 0) then
          !  flag = -2
          !  goto 20
          !end if
          allocate(X(2*(n1/2+1)),stat=stat)
          if (stat .ne. 0) then
            flag = -2
            goto 20
          end if
!             equivalence(X_2D,X)
          L(1) = n1
          L(2) = 1

          strides_in(1) = 0
          strides_in(2) = 1
          strides_in(3) = 2*(n1/2+1)
          strides_out(1) = 0
          strides_out(2) = 1
          strides_out(3) = n1/2+1

!$          tm1 = omp_get_wtime()
          Status = DftiCreateDescriptor( My_Desc_Handle, DFTI_DOUBLE,&
            DFTI_REAL, 1, L(1) )
  !        write(*,*) 'Status1', Status
          if (status .ne. 0) then
            if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
                write(*,*) 'Error: ', DftiErrorMessage(status)
            endif
           endif 


          Status = DftiSetValue(My_Desc_Handle, DFTI_CONJUGATE_EVEN_STORAGE,&
            DFTI_COMPLEX_COMPLEX)
!          write(*,*) 'Status2', Status

          if (status .ne. 0) then
            if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
                write(*,*) 'Error: ', DftiErrorMessage(status)
            endif
endif

          !Status = DftiSetValue(My_Desc_Handle, DFTI_INPUT_STRIDES, strides_in)
          !write(*,*) 'Status3', Status

          !Status = DftiSetValue(My_Desc_Handle, DFTI_OUTPUT_STRIDES, &
          !  strides_out)
          !write(*,*) 'Status4', Status

          Status = DftiCommitDescriptor( My_Desc_Handle)
!          write(*,*) 'Status5', Status

          if (status .ne. 0) then
            if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
                write(*,*) 'Error: ', DftiErrorMessage(status)
            endif
endif

!$          tm2 = omp_get_wtime()
            tm_fft_init = tm_fft_init + tm2 - tm1


        end if

        ! Copy slice from C(:,:,k)
        do i=1,n1
          do j=1,1
            t = 0.1**12
            if (C(i,k).lt. t) then
               s = 0.0_wp
            else
               s = C(i,k)
            end if
             X(i+(j-1)*2*(n1/2+1)) = s
          end do
        end do


!$      tm1 = omp_get_wtime()


        Status = DftiComputeForward( My_Desc_Handle, X )
          if (status .ne. 0) then
            if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
                write(*,*) 'Error: ', DftiErrorMessage(status)
            endif
endif

!$      tm2 = omp_get_wtime()
            tm_fft = tm_fft + tm2 - tm1
!        write(*,*) 'fft time=', tm2-tm1

        if (check) then
     
          if (k.eq.1) then

!$      tm1 = omp_get_wtime()
            Status = DftiCreateDescriptor( My_Desc_Handle_Inv, DFTI_DOUBLE,&
              DFTI_REAL, 1, L(1) )
          if (status .ne. 0) then
            if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
                write(*,*) 'Error: ', DftiErrorMessage(status)
            endif
endif

            Status = DftiSetValue(My_Desc_Handle_Inv,&
              DFTI_CONJUGATE_EVEN_STORAGE,&
              DFTI_COMPLEX_COMPLEX)
          if (status .ne. 0) then
            if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
                write(*,*) 'Error: ', DftiErrorMessage(status)
            endif
endif

           ! Status = DftiSetValue(My_Desc_Handle_Inv, DFTI_INPUT_STRIDES,&
           !   strides_out)
           ! Status = DftiSetValue(My_Desc_Handle_Inv, DFTI_OUTPUT_STRIDES, &
           !   strides_in)
            Status = DftiCommitDescriptor( My_Desc_Handle_Inv)
          if (status .ne. 0) then
            if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
                write(*,*) 'Error: ', DftiErrorMessage(status)
            endif
endif

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


!$          tm1 = omp_get_wtime()

          Status = DftiComputeBackward( My_Desc_Handle_Inv, X )

!$          tm2 = omp_get_wtime()
!            write(*,*) 'ifft time=', tm2-tm1

            tm_ifft = tm_ifft + tm2 - tm1



        ! Copy slice from X to Dk
        do i=1,n1
          do j=1,1

Dk(i,j) = X(i+(j-1)*2*(n1/2+1))/real(n1,kind=wp)
          end do
        end do
call check_error(n1,C(:,k),Dk,nrm)
 !       write(*,*) 'nrm',nrm,k

          if (k.eq.n2) then

            write(*,*) 'k,nrm',k,nrm

            Status = DftiFreeDescriptor(My_Desc_Handle_Inv)

          end if

        end if

        if (k.eq.n2) then
          Status = DftiFreeDescriptor(My_Desc_Handle)

          deallocate(X,Dk,stat=stat)
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


  end subroutine

  subroutine fft_bench_fftw(n1,n2,C,check,flag,tm_fft_init,tm_fft,&
       tm_ifft_init,tm_ifft)
    integer, intent(in) :: n1,n2 ! Array dimensions
    real (kind=wp), intent(in) :: C(n1,n2) ! Input array 
    logical, intent(in) :: check ! Additionally, perform inverse, element-wise
                                 ! division by Bi and compare with A
    integer, intent(out) :: flag ! 0 : all fine
                                 ! -1: error: check is true but A or Bi missing
    real(kind=wp), intent(out) :: tm_fft_init ! total initialisation time fft
    real(kind=wp), intent(out) :: tm_fft ! total time fft
    real(kind=wp), intent(out) :: tm_ifft_init ! total initialisation time ifft
    real(kind=wp), intent(out) :: tm_ifft ! total time ifft

    ! Local variables and arrays
    complex(kind=wp), allocatable :: Dk(:,:)
    real(kind=wp) :: nrm,tm1,tm2
    integer :: stat, k, i, nthreads
    integer(kind=4) :: n1_4, flags,nthreads_4


    type(C_PTR) :: plan, iplan

    real(C_DOUBLE), dimension(n1) :: in, iout
    real(C_DOUBLE), dimension(n1) :: out, iin

    flag = 0

    tm_fft_init = 0.0_wp
    tm_fft = 0.0_wp
    tm_ifft_init = 0.0_wp
    tm_ifft = 0.0_wp
    nthreads = 1
!$  nthreads=omp_get_max_threads()  
    nthreads_4 = int(nthreads,kind=ip4)

           n1_4 = int(n1,kind=ip4)
           flags = int(0,kind=ip4)
!$          tm1 = omp_get_wtime()
           stat = fftw_init_threads()
           call fftw_plan_with_nthreads(nthreads_4)
           plan = fftw_plan_r2r_1d(n1_4, in,out,FFTW_R2HC,flags)
!$          tm2 = omp_get_wtime()
            tm_fft_init = tm_fft_init + tm2 - tm1

       if (check) then

!$   tm1 = omp_get_wtime()
           iplan = fftw_plan_r2r_1d(n1_4, iin,iout,FFTW_HC2R,&
                   flags)
!$   tm2 = omp_get_wtime()
    
           tm_ifft_init = tm_ifft_init + tm2 - tm1

           allocate(Dk(n1,1),stat=stat)
           if (stat .ne. 0) then
             flag = -2
             goto 20
           end if


       end if

       do k=1,n2


        ! Copy each slice into in
        do i=1,n1
       !     write(*,*) 'c',i,j,k,C(i,j,k)
            in(i) = C(i,k)
       !     write(*,*) i,j,k,Dk(i,j)

        end do

!$   tm1 = omp_get_wtime()
        call fftw_execute_r2r(plan, in, out)
!$   tm2 = omp_get_wtime()
        tm_fft = tm_fft + tm2 - tm1

!        write(*,*) 'fft time=', tm2-tm1

        if (check) then

         ! Copy out into iin
         do i=1,n1
        !     write(*,*) 'c',i,j,k,C(i,j,k)
             iin(i) = out(i)
        !     write(*,*) i,j,k,Dk(i,j)

        end do

!$      tm1 = omp_get_wtime()
        call fftw_execute_r2r(iplan, iin, iout)
!$      tm2 = omp_get_wtime()
!        write(*,*) 'ifft time=', tm2-tm1
        tm_ifft = tm_ifft + tm2 - tm1

          if (k.eq.1) then
            nrm = 0.0_wp
          end if

!$OMP PARALLEL DO PRIVATE(i)
          do i=1,n1
             Dk(i,1) = real(iout(i),kind=wp)/real(n1,kind=wp)
          end do
!$OMP END PARALLEL DO

!          write(*,*) iout(n1/2,n2/2), Dk(n1/2,n2/2), C(n1/2,n2/2,k)

          call check_error(n1,C(:,k),Dk,nrm)


         if (k.eq.n2) then

            write(*,*) 'k, nrm^2:',k,nrm

           call fftw_destroy_plan(iplan)
           deallocate(Dk, stat=stat)
           if (stat .ne. 0) then
             flag = -3
             goto 20
           end if
         
         end if


        end if


   
        if (k.eq.n2) then
           call fftw_destroy_plan(plan)
           call fftw_cleanup_threads()
        end if

      

         end do



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


  end subroutine


  subroutine check_error(n1,A,C,nrm)
    integer, intent(in) :: n1 ! Array dimensions
    real(kind=wp), intent(in) :: A(n1,1) ! Input array A
    complex(kind=wp), intent(in) :: C(n1,1) ! Input array C
    real(kind=wp), intent(inout) :: nrm ! 2-norm of A-C

    ! local variables
    integer :: i
    complex(kind=wp) :: s, t
 !   write(*,*) 'nrm',nrm

!$OMP PARALLEL DO PRIVATE(s,t) REDUCTION(+:nrm)
    do i = 1,n1

          s= cmplx(A(i,1),kind=wp)-C(i,1)
          t = s*conjg(s)
!          write(*,*) 's,t',s,t, real(t,kind=wp)
          nrm = nrm + real(t,kind=wp)
    end do
!$OMP END PARALLEL DO
  end subroutine

END PROGRAM commandline
