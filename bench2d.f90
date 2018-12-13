PROGRAM commandline
  IMPLICIT NONE

  Integer, Parameter :: wp = Selected_Real_Kind(15,307)  ! double real

  INTEGER :: nargs, n1, n2, n3, nq, fftlib ! Input arguments
  integer :: stat ! Allocat/deallocate stat 
  integer :: flag ! Error flag
  integer :: i, j, k, l, qq, p, m ! indices
  real(kind=wp) :: xo, yo, zo, a1, b1, c1, r  ! Used in definition of ellipsoid
  CHARACTER(LEN=100) :: option1, option2, option3, option4, option5 
  ! For reading inputs
  real (kind=wp), allocatable :: A(:,:,:) ! Array A
  real (kind=wp), allocatable :: B(:,:,:,:) ! B(:,:,:,i) is cuboid B_i
  real (kind=wp), allocatable :: C(:,:,:) ! C(:,:,:) is cuboid C_i
  real (kind=wp) :: s1,s2,s3,t1,t2 ! used to define B
  logical :: init, check

  !Read input from the command line
  nargs = IARGC()
  IF (nargs .ne. 5) THEN
    goto 10
  ELSE
    CALL GETARG(1,option1) !Grab the first command line argument
    ! and store it in temp variable 'option1'
    CALL GETARG(2,option2) !Grab the 2nd command line argument
    ! and store it in temp variable 'option2'
    CALL GETARG(3,option3) !Grab the 3rd command line argument
    ! and store it in temp variable 'option3'
    CALL GETARG(4,option4) !Grab the 4th command line argument
    ! and store it in temp variable 'option4'
         !  1: FFTE
         !  2: FFTW
         !  3: MKL
         !  4: P3DFFT
         !  5: P3DFFT++
    CALL GETARG(5,option5) !Grab the 5th command line argument
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

    if (n1.lt.1 .or. n2.lt.1 .or. n3.lt.1 .or. nq.lt.1 .or. fftlib .lt. 1 &
        .or. fftlib .gt. 5) then
      goto 10
    endif

  ENDIF

  ! Allocate arrays
  allocate(A(n1,n2,n3),stat=stat)
  if (stat .ne. 0) then
      write(*, '(a)') "Error allocating A"
      return
  end if

  allocate(B(n1,n2,n3,nq),stat=stat)
  if (stat .ne. 0) then
      write(*, '(a)') "Error allocating B"
      deallocate(A)
      return
  end if
  allocate(C(n1,n2,n3),stat=stat)
  if (stat .ne. 0) then
      write(*, '(a)') "Error allocating C"
      deallocate(A,B)
      return
  end if

  ! Set A
  xo = 0.6*real(n1,wp)
  yo = 0.4*real(n2,wp)
  zo = real(n3,wp)/3.0
  a1 = 0.3*real(n1,wp)
  b1 = 0.35*real(n2,wp)
  c1 = real(n3,wp)/3.0_wp
  do i=1,n1
    do j=1,n2
      do k=1,n3
         r = ((real(i,wp)-xo)/a1)**2 + ((real(j,wp)-yo)/b1)**2 + &
             ((real(k,wp)-zo)/c1)**2
         if (r .le. 1) then
            A(i,j,k) = r
         else
            A(i,j,k) = 0.0_wp
         end if
       !  write(*,*) i,j,k, A(i,j,k)
      end do
    end do
  end do

  ! Set B

  
  m = 0
  do while (3*m+1.le.nq)
    do i=1,n1
      do j=1,n2
        do k=1,n3
         s1=1.0_wp
         s2=1.0_wp
         s3=1.0_wp
         t1=1.0_wp
         t2=1.0_wp
         do p=1,m+1
           s1 = s1*(real(i*(m+2),wp)/real(p*n1,wp) - 1.0_wp )
           s2 = s2*(real(j*(m+2),wp)/real(p*n2,wp) - 1.0_wp )
           s3 = s3*(real(k*(m+2),wp)/real(p*n3,wp) - 1.0_wp )
         end do
         do qq=1,m
           t1 = t1*(real(j*(m+1),wp)/real(qq*n2,wp) - 1.0_wp )
           t2 = t2*(real(k*(m+1),wp)/real(qq*n3,wp) - 1.0_wp )
         end do
         B(i,j,k,3*m+1) = s1*t1*t2
         if (3*m+2 .le. nq) then
           B(i,j,k,3*m+2) = s1*s2*t2
         end if
         if (3*m+3 .le. nq) then
           B(i,j,k,3*m+3) = s1*s2*s3
         end if
        end do
      end do
    end do
    m=m+1
  end do

  ! Normalise norm(B(i,j,k,:),2) to equal 1 
  do i=1,n1
    do j=1,n2
      do k=1,n3
        s1=0.0_wp
        do qq=1,nq
          s1 = s1 + (B(i,j,k,qq))**2
        end do 
        s1 = s1**0.5
        do qq=1,nq
          B(i,j,k,qq) = B(i,j,k,qq)/s1
        end do 
      end do
    end do
  end do

  ! Set init to true so that fft initilisation performed first
  init = .true.

  ! Set-up each 2D slice and perform FFT
  ! Each slice formed in C(:,:,:) by performing element-wise multiplaction of 
  ! A with B(:,:,:,qq)
  do qq=1,nq
    do i=1,n1
      do j=1,n2
        do k=1,n3
          C(i,j,k) = A(i,j,k)*B(i,j,k,qq)
        end do
      end do
    end do
    
  ! Perform FFT on each 2D slice
    check=.true.
    call fft_bench(n1,n2,n3,C,fftlib,init,check,flag,A=A,Bi=B(:,:,:,qq))

  end do
  

  


  ! Deallocate arrays
  deallocate(A,B,C, stat=stat)
  if (stat .ne. 0) then
      write(*, '(a)') "Error deallocating arrays"
      return
  end if

  return

10  write(*,'(a)') "usage: ./bench2d.exe n1 n2 n3 nq fftlib"
    write(*,'(a)') " n1=positive integer : first dimension of cuboid"
    write(*,'(a)') " n2=positive integer : second dimenstion of cuboid"
    write(*,'(a)') " n3=positive integer : third dimenstion of cuboid"
    write(*,'(a)') " nq=positive integer : number of multiplying cuboids"
    write(*,'(a)') " fftlib=positive integer less than 6: FFT library to use"
    write(*,'(a)') "   fftlib=1: FFTE"
    write(*,'(a)') "   fftlib=2: FFTW"
    write(*,'(a)') "   fftlib=3: MKL"
    write(*,'(a)') "   fftlib=4: P3DFFT"
    write(*,'(a)') "   fftlib=5: P3DFFT++"
    return


 contains
   
  subroutine fft_bench(n1,n2,n3,C,fftlib,init,check,flag,A,Bi)
    integer, intent(in) :: n1,n2,n3 ! Array dimensions
    real (kind=wp), intent(in) :: C(n1,n2,n3) ! Input array 
    integer, intent(in) :: fftlib ! fft library to use
         !  1: FFTE
         !  2: FFTW
         !  3: MKL
         !  4: P3DFFT
         !  5: P3DFFT++
    logical, intent(inout) :: init  ! Has fft routine been initialise?
    logical, intent(in) :: check ! Additionally, perform inverse, element-wise
                                 ! division by Bi and compare with A
    integer, intent(out) :: flag ! 0 : all fine
                                 ! -1: error: check is true but A or Bi missing
    real(kind=wp), intent(in), optional :: A(n1,n2,n3) ! Input array A
    real(kind=wp), intent(in), optional :: Bi(n1,n2,n3) ! Input array Bi

    ! Local variables and arrays
    complex(kind=wp), allocatable :: Dk(:,:), work(:,:)
    real(kind=wp) :: nrm
    integer :: stat, k, i, j, iopt, ntemp

    flag = 0
    if (check .and. ((.not. present(A)) .or. (.not. present(Bi)))) then
      flag = -1 ! Should not be possible to reach this error flag
      goto 20
    end if
 

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
        if (init) then
          call DZFFT2D(Dk,n1,n2,0,work)
        end if 
        call DZFFT2D(Dk,n1,n2,-1,work)
    !    do i=1,n1/2+1
    !      do j=1,n2
    !        write(*,*) i,j,k,Dk(i,j)
    !      end do
    !    end do
        if (check) then
          if (init) then
            call ZDFFT2D(Dk,n1,n2,0,work)
          end if
          call ZDFFT2D(Dk,n1,n2,1,work)

          if (k.eq.1) then
            nrm = 0.0_wp
          end if
          call check_error(n1,n2,C(:,:,k),Dk,nrm)

          if (k.eq.n3) then
            !nrm = sqrt(nrm)
            write(*,*) 'k, nrm^2:',k,nrm
          end if
        end if
        
        if (init) then
          init = .false.
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

    ! case (2)



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


  subroutine check_error(n1,n2,A,C,nrm)
    integer, intent(in) :: n1,n2 ! Array dimensions
    real(kind=wp), intent(in) :: A(n1,n2) ! Input array A
    complex(kind=wp), intent(in) :: C(n1,n2) ! Input array C
    real(kind=wp), intent(inout) :: nrm ! 2-norm of A-C

    ! local variables
    integer :: i,j,k
    complex(kind=wp) :: s, t

    do i = 1,n1
      do j= 1,n2
        do k = 1,n3
          s= cmplx(A(i,j))-C(i,j)
          t = s*conjg(s)
          nrm = nrm + real(t,kind=wp)
        end do
      end do
    end do
    


  end subroutine

END PROGRAM commandline
