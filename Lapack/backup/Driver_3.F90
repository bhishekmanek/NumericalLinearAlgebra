   program Driver_lapack

   use linal

   implicit none

   character*100 filename1,SYM_CHECK
   INTEGER msize,nsize,i,j,k
   real, dimension(:,:), allocatable, save :: Asave,mat,A,rem,a_v,e_value,a_v1,rem1
   real, dimension(:,:), allocatable, save :: q_mat,matV
   real, dimension(:),allocatable, save :: evalue
   real, save :: trace_mat,sum_eigen

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

   allocate (matV(msize,msize),rem(msize,msize),a_v(msize,msize),e_value(msize,msize))
   allocate (Asave(msize,msize),rem1(msize,msize),a_v1(msize,msize),evalue(msize))

   Asave = mat

  if (SYM_CHECK == 'YES') then

   trace_mat = trace(Asave,msize)
   call qr_ev(mat,msize,nsize,q_mat,matV,evalue)
   e_value = matmul(matmul(transpose(q_mat),mat),q_mat)
   a_v1 = matmul(Asave,matV)
   do i=1,msize
      evalue(i)=e_value(i,i)
   enddo
   do i = 1,msize
   rem1(:,i) = a_v1(:,i) - evalue(i)*matV(:,i)
   enddo

  elseif (SYM_CHECK == 'NO') then

   write(*,*)"NOTE - The input is real but non-symmetric."
   write(*,*)"Continuing with the QR algorithm ..."
   write(*,*)" "
   trace_mat = trace(Asave,msize)
   call qr_ev(mat,msize,nsize,q_mat,matV,evalue)
   e_value = matmul(matmul(transpose(q_mat),mat),q_mat)
   a_v1 = matmul(Asave,matV)
   do i=1,msize
      evalue(i)=e_value(i,i)
   enddo
   do i = 1,msize
   rem1(:,i) = a_v1(:,i) - evalue(i)*matV(:,i)
   enddo

  endif
   

   write(*,*)"The input matrix is:"
   call matrix_print(Asave,msize,msize)
!   write(*,*)" " 
   write(*,*)"Finally converged eigenvectors"
   do i=1,msize
      print*,(matV(i,j),j=1,msize)
   enddo
   write(*,*)" " 
  write(*,*)"Corresponding Eigenvalues:"
   do i=1,msize
      print*,e_value(i,i)
   enddo
   write(*,*)" " 
   write(*,*)'Trace of input Matrix is:',trace_mat
   write(*,*)"  "
   sum_eigen = sum(evalue)
   write(*,*)'Sum of Eigenvalues is : ',sum_eigen
   write(*,*)"  "
   write(*,*)"Av - lambda*v is"

   do i = 1,msize
   print*,rem1(:,i)
   enddo









   deallocate (mat)

  end program Driver_lapack
