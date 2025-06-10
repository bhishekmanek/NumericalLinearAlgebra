program Driver_iterative

use linal

implicit none

character*100 SYM_CHECK
real, dimension(:,:), save, allocatable :: matD
real, dimension(:), save, allocatable :: b,x
real :: diag,accuracy
integer :: msize,i,choice

! Input from the user
write(*,*)'We will be using either Gauss-Jacobi, Gauss-Seidel method or Conjugate gradient (CG)/Smart CG method to solve Ax=b'
write(*,*)" "
write(*,*)'Enter the diagonal element D for the matrix: '
read(*,*)diag
write(*,*)'Enter the size of Matrix'
read(*,*)msize
write(*,*)'Enter the choice of method you want to use:'
write(*,*)'1. Gauss-Jacobi or 2. Gauss-Seidel or 3. Conjugate gradient or 4. Smart Conjugate gradient'
read(*,*)choice

allocate (matD(msize,msize),b(msize),x(msize))

accuracy = 1E-6
matD = 1.
write(*,*)'Accuracy demanded is',accuracy
do i=1,msize
   matD(i,i) = msize + 1. - real(i)
!   matD(i,i) = diag
   b(i) = real(i)
enddo

write(*,*)'The input matrix A is:'
call matrix_print(matD,msize,msize)
write(*,*)" "
write(*,*)'The RHS vector b is:'
do i=1,msize
print*,b(i)
enddo
write(*,*)" "
x = 0.
if (choice == 1) then
call gauss_jacobi(matD,b,msize,accuracy,x)
elseif (choice == 2) then
call gauss_seidel(matD,b,msize,accuracy,x)
elseif (choice == 3) then

   if (frob_norm(matD - transpose(matD),msize,msize) .le. epsilon(0.)) then
       SYM_CHECK = 'YES'
       write(*,*)'The input matrix is symmetric'
   else
       SYM_CHECK = 'NO'
       write(*,*)'The input matrix is not symmetric, stopping ...'
       stop
   endif


call conj_grad(matD,b,msize,accuracy,x)
elseif (choice == 4) then
call smart_cg(matD,b,msize,accuracy,x)
else 
write(*,*)'Err ! Wrong choice. Enter 1, 2, 3 or 4 for your choices'
endif

write(*,*)" "
write(*,*)'The solution vector x for Ax=b is:'
do i=1,msize
write(*,*)x(i)
enddo

end program Driver_iterative
