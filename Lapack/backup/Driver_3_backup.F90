   program Driver_lapack

   use linal

   implicit none

   character*100 filename1,JOBZ,UPLO,SYM_CHECK,JOBVR,JOBVL
   INTEGER N,LWORK,LDA,msize,nsize,INFO,i,j,k,SDIM,LDVL,LDVR
   real, dimension(:,:), allocatable, save :: Asave,mat,A,VL,VR,rem,a_v,e_value,a_v1,rem1
   real, dimension(:), allocatable :: W,WR,WI
   real, dimension(:), allocatable :: WORK
   logical, dimension(:), allocatable :: BWORK
   real, dimension(:,:), allocatable, save :: q_mat,matV
   real, dimension(:),allocatable, save :: evalue

   if (iargc() .ne. 1) then
       write(*,*)"Problem with your filename ! Enter correctly !"
       stop
   endif

   call getarg(1,filename1)

   call readandallocatemat(mat,msize,nsize,filename1)

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

   call ssyev (JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO)

   write(*,*)"The input matrix is:"
   call matrix_print(mat,msize,msize)

   write(*,*)"Eigenvalues of given matrix:"
   do i=1,N
   write(*,*)W(i)
   enddo

   write(*,*)"  "
   write(*,*)"The Eigenvectors are:"
   call matrix_print(A,msize,msize)

   a_v = matmul(mat,A)

   write(*,*)"Av - lambda*v is"
   do i = 1,msize
   rem(:,i) = a_v(:,i) - W(i)*A(:,i)
   print*,rem(:,i)
   enddo


  elseif (SYM_CHECK == 'NO') then

   call sgeev (JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)

   write(*,*)"The input matrix is:"
   call matrix_print(mat,msize,msize)
   write(*,*)"Eigenvalues of given matrix:"
   do i=1,N
   write(*,*)WR(i)
   enddo
   write(*,*)" " 
   write(*,*)"The Eigenvectos are:"
   do i=1,msize
      write(*,*)(VR(i,j),j=1,msize)
   enddo

   a_v = matmul(mat,VR)

   write(*,*)"  "
   write(*,*)"Av - lambda*v is"
   do i = 1,msize
   rem(:,i) = a_v(:,i) - WR(i)*VR(:,i)
   print*,rem(:,i)
   enddo

  endif
   
   write(*,*)" " 
   call qr_ev(mat,msize,nsize,q_mat,matV,evalue)
   write(*,*)"Finally converged eigenvectors"
   do i=1,msize
      print*,(matV(i,j),j=1,msize)
   enddo

   e_value = matmul(matmul(transpose(q_mat),mat),q_mat)
   write(*,*)"Finally converged Eigenvalues:"
   do i=1,msize
      evalue(i)=e_value(i,i)
      print*,e_value(i,i)
   enddo

   write(*,*)" "
   call matrix_print(Asave,msize,msize)
   write(*,*)" " 
   call matrix_print(matV,msize,msize)


   a_v1 = matmul(Asave,matV)

   call matrix_print(a_v1,msize,msize)


   write(*,*)"  dang"
   write(*,*)"Av - lambda*v is"
   do i = 1,msize
   rem1(:,i) =  a_v1(:,i)- evalue(i)*matV(:,i)
   print*,rem1(:,i)
   enddo



   deallocate (A,mat,W,WORK)

  end program Driver_lapack
