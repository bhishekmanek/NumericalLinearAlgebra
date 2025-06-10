! Project3/Question3
! Main Program for calculating e-values fo a Lehmer matrix
! Takes input from user of dimension of Lehmer matrix 

   program Driver_Lehmer

   use linal

   implicit none

   character*100 JOBZ,UPLO,choice,filename
   INTEGER N,LWORK,LDA,msize,nsize,INFO,i,j
   real, dimension(:,:), allocatable, save :: mat,A,rem,a_v,e_value
   real, dimension(:), allocatable :: W,WORK
   real, dimension(:,:), allocatable, save :: q_mat,matV
   real, dimension(:),allocatable, save :: evalue
   real :: i1,j1,temp
   real :: start,finish

   call cpu_time(start)

   write(*,*)"Enter the size of square Matrix: "
   read(*,*)msize
   nsize = msize

   allocate (rem(msize,msize),a_v(msize,msize),e_value(msize,msize))

! Creates Lehmer matrix of input size
   call create_lehmer(mat,msize)

   write(*,*)"What method you want to use for Eigenvalue calculation?"
   write(*,*)"Lapack routine (L) or my routine(M)?"
   read(*,*)choice

! For Symmetric Matrices - ssyev
   JOBZ = 'N'
   UPLO = 'U'
   N = msize
   A = mat
   LDA = max(1,msize)
   LWORK = 1000

   allocate (WORK(LWORK),W(N))

   if (choice == 'L') then

   open(10,file='EV_L.dat',status='unknown')

      call ssyev (JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO)

      write(*,*)"The input matrix is:"
      call matrix_print(mat,msize,msize)

      write(*,*)"Eigenvalues of given matrix:"
      do i=1,N
         write(*,*)W(i)
         write(10,*)i,W(i)
      enddo

   elseif (choice == 'M') then

   open(11,file='EV_M.dat',status='unknown')

      write(*,*)"The input matrix is:"
      call matrix_print(mat,msize,msize)
   
      call qr_ev(mat,msize,nsize,q_mat,matV,evalue)

      e_value = matmul(matmul(transpose(q_mat),mat),q_mat)
      write(*,*)"Finally converged Eigenvalues:"
      do i=1,msize
         evalue(i) = e_value(i,i)
      enddo

      do i=1,msize-1
         do j=i+1,msize
         if (evalue(i) .gt. evalue(j)) then
            temp = evalue(i)
            evalue(i) = evalue(j)
            evalue(j) = temp
         else
           continue
         endif
      enddo
     enddo  

     do i=1,msize
        print*,evalue(i)
        write(11,*)i,evalue(i)
     enddo

   else 
 
      write(*,*)"Wrong choice! Enter L or M"
   
   endif

   deallocate (A,mat,W,WORK)

   call cpu_time(finish)
   print*,'Time taken for the run is',finish-start, 'seconds'

  end program Driver_Lehmer
