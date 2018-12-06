PROGRAM commandline
  IMPLICIT NONE
  INTEGER :: nargs, n1, n2, n3, nq, fftlib ! Input arguments
  integer :: stat ! Allocat/deallocate stat 
  CHARACTER(LEN=100) :: option1, option2, option3, option4, option5
  double precision, allocatable :: A(:,:,:) ! Array A
  double precision, allocatable :: Q(:,:,:,:) ! Q(:,:,:,i) is cuboid Q_i

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
    CALL GETARG(4,option5) !Grab the 5th command line argument
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
  allocate(A(n1,n2,n3,stat=stat)
  if (stat .ne. 0) then
      write(*, '(a)') "Error allocating A"
      return
  end if

  allocate(Q(n1,n2,n3,nq,stat=stat)
  if (stat .ne. 0) then
      write(*, '(a)') "Error allocating Q"
      deallocate(A)
      return
  end if

  ! Set A




  ! Deallocate arrays
  deallocate(A,Q, stat=stat)
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
    return

END PROGRAM commandline
