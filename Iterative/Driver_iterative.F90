! /ams213a-bhishek/Project4/Question1
! Dirver routine for Iterative methods to solve for Ax=b
! Following driver can call the following methods- 
! 1. Gauss-Jacobi
! 2. Gauss-Seidel 
! 3. Conjugate gradient 
! 4. Smart conjugate gradient
! 5. Diagonal conditioning applied to conjugate gradient 

! Input - matD,initial guess x,accuracy,diagonal element of matD, choice of method
! Ouput - x, Error in calculation at different iteration step

program Driver_iterative

use linal

implicit none

character*100 SYM_CHECK
real, dimension(:,:), save, allocatable :: matD
real, dimension(:), save, allocatable :: b,x
real :: diag,accuracy
integer :: msize,i,choice

! Input from the user
write(*,*)'We will be using either Gauss-Jacobi, Gauss-Seidel method or'    
write(*,*)'Conjugate gradient (CG)/Smart CG or Diagonally preconditioned conjugate gradient method to solve Ax=b'
write(*,*)" "
write(*,*)'Enter the diagonal element D for the matrix: '
read(*,*)diag
write(*,*)'Enter the size of Matrix'
read(*,*)msize
write(*,*)'Enter the choice of method you want to use:'
write(*,*)'1. Gauss-Jacobi'
write(*,*)'2. Gauss-Seidel'
write(*,*)'3. Conjugate gradient'
write(*,*)'4. Smart Conjugate gradient'
write(*,*)'5. Diagonally pre-conditioned CG method'
read(*,*)choice

allocate (matD(msize,msize),b(msize),x(msize))

accuracy = 1E-6
matD = 1.

write(*,*)'Accuracy demanded is',accuracy

! Defining the diagonal elements of the matrix. Comment out according to your requirement
do i=1,msize
!   matD(i,i) = real(i)
   matD(i,i) = diag
   b(i) = real(i)
enddo

! Inital guess solution x
x = 0.

if (choice == 1) then
   call gauss_jacobi(matD,b,msize,accuracy,x)
elseif (choice == 2) then
   call gauss_seidel(matD,b,msize,accuracy,x)
elseif (choice == 3) then

! Check whether the matrix is symmetric or not.
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
elseif (choice == 5) then
   call diag_pre_cond(matD,b,msize,accuracy,x)
else 
   write(*,*)'Err ! Wrong choice. Enter 1, 2, 3 or 4 for your choices'
endif

write(*,*)" "
write(*,*)'The solution vector x for Ax=b is:'

! Uncomment if you want solutions to be printed
!do i=1,msize
!write(*,*)x(i)
!enddo

deallocate (matD,b,x)

end program Driver_iterative
