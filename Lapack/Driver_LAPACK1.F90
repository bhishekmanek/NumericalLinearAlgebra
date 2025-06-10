! Project3/Question2/Lapack-part1
! Main program for calling LAPACK routines to calculate e-values and evectors
! for any input matrix

   program Driver_lapack

   use linal

   implicit none

   character*100 filename1,JOBZ,UPLO,SYM_CHECK,JOBVR,JOBVL
   INTEGER N,LWORK,LDA,msize,nsize,INFO,i,j,k,SDIM,LDVL,LDVR
   real, dimension(:,:), allocatable, save :: Asave,mat,A,VL,VR,rem,a_v
   real, dimension(:,:), allocatable, save :: e_value,a_v1,rem1,q_mat,matV
   real, dimension(:), allocatable :: W,WR,WI,WORK
   logical, dimension(:), allocatable :: BWORK
   real, dimension(:),allocatable, save :: evalue
   real :: sum_eigen, trace_mat

   if (iargc() .ne. 1) then
       write(*,*)"Problem with your filename ! Enter correctly !"
       stop
   endif

   call getarg(1,filename1)

   call readandallocatemat(mat,msize,nsize,filename1)

! Checking whether matrix is symmetric or not

   if (frob_norm(mat - transpose(mat),msize,msize) .le. epsilon(0.)) then
       SYM_CHECK = 'YES'
   else
       SYM_CHECK = 'NO'
   endif

   allocate (W(msize),rem(msize,msize),a_v(msize,msize),e_value(msize,msize))
   allocate (Asave(msize,msize),rem1(msize,msize),a_v1(msize,msize))

! For Symmetric Matrices - ssyev
   JOBZ = 'V'
   UPLO = 'U'
   N = msize
   A = mat
   LDA = msize
   LWORK = 50
   Asave = mat

   allocate (WORK(LWORK))

! For Non-symmetrix Matrices - sgeev

   JOBVL = 'V'
   JOBVR = 'V'
   LDVL = msize
   LDVR = msize
   INFO = 0

   allocate (WR(N),WI(N),VL(LDVL,msize),VR(LDVR,msize))

  if (SYM_CHECK == 'YES') then

     trace_mat = trace(Asave,msize)
     call ssyev (JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO)
          a_v = matmul(mat,A)

          do i = 1,msize
             rem(:,i) = a_v(:,i) - W(i)*A(:,i)
          enddo 

  elseif (SYM_CHECK == 'NO') then

     trace_mat = trace(Asave,msize)
     call sgeev (JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
          a_v = matmul(mat,VR)

     do i = 1,msize
        rem(:,i) = a_v(:,i) - WR(i)*VR(:,i)
     enddo

  endif


! Printing everything in terminal

   write(*,*)"The input matrix is:"
   call matrix_print(mat,msize,msize)
   write(*,*)"Eigenvalues of given matrix:"

   if (SYM_CHECK == 'YES') then
      do i=1,N
         write(*,*)W(i)
      enddo

      sum_eigen = sum(W)
      write(*,*)"  "
      write(*,*)"The Eigenvectors are:"
      call matrix_print(A,msize,msize)
      write(*,*)" "
      write(*,*)"Av - lambda*v is"
      do i = 1,msize
         print*,rem(:,i)
      enddo
      write(*,*)" "

      print*,"The trace of input matrix is:",trace_mat
      write(*,*)" "
      print*,"The sum of Eigenvalues is:",sum_eigen

   elseif (SYM_CHECK == 'NO') then

      write(*,*)" "

      do i=1,N
         write(*,*)WR(i)
      enddo

      sum_eigen = sum(WR)

      write(*,*)" "
      print*,"The trace of input matrix is:",trace_mat
      write(*,*)" "
      print*,"The sum of Eigenvalues is:",sum_eigen

   endif

   deallocate (A,mat,W,WORK,WR,WI,Asave,VL,VR,rem,a_v)
   deallocate (e_value,a_v1,rem1)

   end program Driver_lapack
