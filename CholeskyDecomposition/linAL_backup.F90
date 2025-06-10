! /ams213a-bhishek/Project1/LU-Decomposition

! This module contains Linear Algebra routines as well as basic matrix properties calculation routines

! List of Subroutines/Functions:
! 1. readandallocatemat - Reads matrix from an input file and allocates
! 2. trace - This function calculates trace of a given input square matrix
! 3. euclid_norm - Calculates the euclidean norm of a given input vector
! 4. matrix_print - Prints any given input matrix in readable form
! 5. eliminate_gauss - Performs Gaussian elimination with partial pivoting
! 6. gauss_backsub - Backsubstitution to calculate X for the problem of type UX = Y
! 7. lu_decomp - Performs LU decomposition with partial pivoting
! 8. lu_backsub - Performs backsubtitution to solve LUX = PY
! 9. error_mat - Calculates error matrix E = AX - B 

module LinAl

implicit none

contains

!********************************************************

subroutine readandallocatemat(mat,msize,nsize,filename)

character*100 filename
real, dimension(:,:), allocatable, intent(in out) :: mat
integer, intent(OUT) :: msize, nsize

integer :: i,j

! Reads a file containing the matrix A 
! Sample file:
! 3  4
! 1.2     1.3     -1.4    3.31
! 31.1    0.1     5.411   -1.23
! -5.4    7.42    10      -17.4
! Note that the first 2 lines are the matrix dimensions, 
! then the next msize lines are the matrix entries
! Note that entries must be separated by a tab.
! Then allocates an array of size (msize,nsize), populates the matrix,
! and returns the array. 

! This routine takes as INPUTS:
! filename = a character string with the name of file to read
! This routine returns as OUTPUTS:
! msize = first dimension of read matrix
! nsize = second dimension of read matrix
! mat = read matrix. Note that mat type is real.

open(10,file=filename)

! Read the matrix dimensions
read(10,*) msize,nsize

! Allocate matrix
allocate(mat(msize,nsize))

! Read matrix
do i=1,msize
   read(10,*) ( mat(i,j), j=1,nsize )
enddo

close(10)

end subroutine readandallocatemat

!***********************************************************************

! This function calculates trace of input matrix of dimension dim_mat

real function trace(matrix,dim_mat)

implicit none

! I/O Variables
real, dimension(dim_mat,dim_mat),intent(IN) :: matrix
integer, intent (IN) :: dim_mat 

! Local variables
integer :: i

trace = 0.

do i=1,dim_mat
trace = trace + matrix(i,i)
enddo

end function trace

!*********************************************************************

! This function calculates the Euclidean norm for a given input vector
! of dimension dim_vec


real function euclid_norm(vector,dim_vec)

implicit none

! I/O variables
real, dimension(dim_vec), intent(IN) :: vector
integer, intent(IN) :: dim_vec

! Local variables
integer i

euclid_norm = 0.

do i=1,dim_vec
   euclid_norm = euclid_norm + (vector(i)**2)
enddo

euclid_norm = sqrt(euclid_norm)

end function euclid_norm

!*******************************************************************

! This routines prints an input matrix of dimension m by n in 
! readable form


subroutine matrix_print(matrix,n,m)

implicit none

! I/O variables
real, dimension(m,n), intent(IN), allocatable :: matrix(:,:)
integer, intent(IN) :: n,m

! Local variables
integer i,j

do i=1,n,1
   write(*,*)(matrix(i,j),j=1,m)
enddo

   write(*,*) " "

end subroutine matrix_print

!*********************************************************
! Gaussian elimination routine - 
! Input - Matrix A and B, and their dimensions
! Output - Upper triangular matrix A 
!        - Matrix B modified similar to operations applied on A
!        - Logical ising which indicated if the matrix A 
!          is singular or not, TRUE if det(A)=0
!          FALSE is det(A) /= 0

subroutine eliminate_gauss(matA,matB,nsize,msize,nsize1,ising)

implicit none

! I/O variables
real, intent(INOUT), allocatable :: matB(:,:)
real, intent(INOUT), allocatable :: matA(:,:)
logical, intent(OUT) :: ising

integer, intent(IN) :: nsize,msize,nsize1

! Local variables
real :: pivot
integer, dimension(1) :: pivot_loc
real :: mult,det_A
real, dimension(:,:),allocatable :: tempa,tempb
integer :: i,j,k

allocate (tempa(msize,msize),tempb(msize,msize))

do k=1,msize-1

pivot = maxval(abs(matA(k:msize,k))) ! Finds Pivot
pivot_loc = maxloc(abs(matA(k:msize,k))) + k-1 ! Finds Pivot Location

if (pivot_loc(1) .ne. k) then
! Swapping of rows in matA
    tempa(k,1:msize) = matA(k,1:msize)
    matA(k,1:msize) = matA(pivot_loc(1),1:msize)
    matA(pivot_loc(1),1:msize) = tempa(k,1:msize)

    tempb(k,1:nsize1) = matB(k,1:nsize1)
    matB(k,1:nsize1) = matB(pivot_loc(1),1:nsize1)
    matB(pivot_loc(1),1:nsize1) = tempb(k,1:nsize1)
endif

     if (matA(k,k) == 0.) then
        write(*,*)'Diagonal Element 0 !!!'
        stop
     endif


   do i=k+1,msize
! Calculating multiplying factors
     mult = matA(i,k)/matA(k,k)
      matA(i,k) = 0.
!     write(*,*)matA(k,k)
!       do j=k+1,msize
          matA(i,k+1:msize) = matA(i,k+1:msize) - mult*matA(k,k+1:msize)
          matB(i,1:nsize1) = matB(i,1:nsize1) - mult*matB(k,1:nsize1)
!        enddo
   enddo     
enddo

det_A = 1.0

do i=1,nsize
   do j=1,msize
      det_A = det_A * matA(i,i)
   enddo
enddo

if (det_A == 0.) then
   ising = .TRUE.
else
   ising = .FALSE.
endif

end subroutine eliminate_gauss

!***********************************************************

subroutine gauss_backsub(matA,matB,msize,nsize1,x,ising)

implicit none

real, intent(IN), allocatable :: matA(:,:),matB(:,:)
real, intent(OUT), allocatable :: x(:,:)
integer, intent(IN) :: msize,nsize1
logical, intent(IN) :: ising
integer :: i,j,k
real :: summ

allocate (x(msize,nsize1))
if (ising .eqv. .TRUE.) then
      write(*,*)"Matrix is singular! Doing backsubstitution is crazy so quitting!"
   stop
else

do i=1,msize
!write(*,*)(matA(i,j),j=1,msize)
enddo

do i=msize,1,-1
   if (matA(i,i) == 0.) then
     write(*,*)'Dividing by zero is never a good thing! Try some other matrix!'
      stop
   endif

!summ=0.0d0
do k=1,nsize1
summ=0.0
do j=i+1,msize
summ = summ + matA(i,j)*x(j,k)
enddo
x(i,k) = (matB(i,k)-summ)/matA(i,i)
enddo
enddo
endif
end subroutine gauss_backsub

!***********************************************************

! LU decomposition matrix
! Inputs are matrix A, its dimension, and ising (logical)
! Outputs are - 
! A = LU (A written in form such that the lower diagonal elements are from L)
!         and the Upper diagonal elements

subroutine lu_decomp(matA,msize,ising,s)

implicit none

! I/O variables
real, dimension(:,:), allocatable, intent(INOUT) :: matA
integer, intent(IN) :: msize
logical, intent(IN) :: ising
integer, dimension(:), allocatable, intent(OUT) :: s

! Local Variables
real, dimension(:,:), allocatable :: tempa
real, dimension(:), allocatable :: temps
integer n,i,j,k
real :: pivot
integer, dimension(1) :: pivot_loc
real :: mult

! Allocating the temporary local matrices used for swapping and also permutation vector
allocate (s(msize), tempa(msize,msize), temps(msize))

! Initializing Permutation Vector
do j=1,msize
   s(j) = j
enddo

if (ising .eqv. .FALSE.) then

! Loop over column
do j=1,msize
pivot = maxval(abs(matA(j:msize,j))) ! Finds Pivot
pivot_loc = maxloc(abs(matA(j:msize,j))) + j-1 ! Finds Pivot Location
if (pivot_loc(1) .ne. j) then
! Swapping of rows in matA
    tempa(j,1:msize) = matA(j,1:msize)
    matA(j,1:msize) = matA(pivot_loc(1),1:msize)
    matA(pivot_loc(1),1:msize) = tempa(j,1:msize)

! Swapping of Permutation vector
    temps(j) = s(j)
    s(j) = s(pivot_loc(1))
    s(pivot_loc(1)) = temps(j)

endif

if (matA(j,j) == 0) then
write(*,*)'Diagonal Element 0 !!!'
stop
endif

! Creating L and storing in matA
do i=j+1,msize
matA(i,j) = matA(i,j)/matA(j,j) ! Multiplication factor
   do k=j+1,msize
! Creates Lower matrix and stores in matA instead of creating L and U separately
      matA(i,k) = matA(i,k) - matA(i,j)*matA(j,k)
   enddo
enddo

enddo

else
write(*,*)'Matrix is Singluar ! Cannot proceed with LU decomposition'
endif
end subroutine lu_decomp

!*********************************************************************************

! LU Backsubstitution 
! Input : matrix A and B, permutation vector s, dimensions of A, B, s
! Input & Output : matrix X after doing backsubstitution

subroutine lu_backsub(matA,msize,nsize1,matB,s,matX)

implicit none

! I/O variables
real, dimension(:,:), allocatable, intent(IN) :: matA, matB
real, dimension(:,:), allocatable, intent(INOUT) :: matX
integer, dimension(:), allocatable, intent(IN) :: s
integer, intent(IN) :: msize,nsize1

! Local variables
real, dimension(:,:),allocatable :: temp
real :: summ
integer :: i,j,k
real, dimension(msize) :: Bc,Xc
real, dimension(msize,nsize1) :: y

allocate (matX(msize,nsize1), temp(msize,nsize1))

temp = matB

! Calculating Pb and storing it in y
do j=1,msize
   y(j,:) = temp(s(j),:)
enddo

! Forward substitution to calculate y = L_inv*P*b
do j=1,msize-1
   do i=j+1,msize
      y(i,1:nsize1) = y(i,1:nsize1) - y(j,1:nsize1)*matA(i,j)
   enddo
enddo

! Backward substitution to get X, UX = Y
do i=msize,1,-1
   if (matA(i,i) == 0.) then
      write(*,*)'Diagonal element is zero! No Backward substitution'
      stop
    endif

do j=1,nsize1
summ = 0.0
   do k=i+1,msize
      summ = summ + matA(i,k)*matX(k,j)
   enddo
   matX(i,j) = (y(i,j)-summ)/matA(i,i)
enddo
enddo

end subroutine lu_backsub

!************************************************************

subroutine cholesky_fact(matA,msize,ising)

implicit none

real, intent(INOUT), allocatable :: matA(:,:)
integer, intent(IN) :: msize
logical, intent(IN) :: ising

integer :: i,j,k

if (ising .eqv. .FALSE.) then
do j=1,msize
   do k=1,j-1
     matA(j,j) = matA(j,j) - matA(j,k)*matA(j,k)
   enddo
   matA(j,j) = sqrt(matA(j,j))


do i=j+1,msize
   do k=1,j-1
      matA(i,j) = matA(i,j) - matA(i,k)*matA(j,k)
   enddo
   matA(i,j) = matA(i,j)/matA(j,j)
enddo
enddo

else 
write(*,*)'Matrix is singular or not positive definite!'
stop
endif

end subroutine cholesky_fact

!***********************************************************

subroutine cholesky_back(matA,b,msize,x,y)

implicit none

real, intent(IN), allocatable :: matA(:,:),b(:)
integer, intent(IN) :: msize

real, intent(OUT), allocatable :: x(:),y(:)

integer :: i,j,k
real :: summ

allocate (x(msize),y(msize))


do i=1,msize
   summ = b(i)
   do j=1,i-1
      summ = summ - y(j)*matA(i,j)
   enddo
   y(i) = summ/matA(i,i)
enddo

do i=msize,1,-1
   if (matA(i,i) == 0.) then
      stop
   endif
   do k=i+1,msize
      y(i) = y(i)- matA(k,i)*x(k)
   enddo
   x(i) = y(i)/matA(i,i)
enddo


end subroutine cholesky_back

!**********************************************************

subroutine rms_error(x,msize,error)

real, dimension(msize), intent(IN) :: x
integer, intent(IN) :: msize

real, intent(OUT) :: error

real :: summ,fit_num,xc,xval
!real, dimension(:), allocatable :: xc,xval
integer :: i,j,num_points,nsize

open(19,file='atkinson.dat',status='unknown')
read(19,*)num_points,nsize
!allocate (xc(num_points),xval(num_points))

summ = 0.0
do i=1,num_points
read(19,*)xc,xval
fit_num = 0.0
do j=msize,1,-1
fit_num = fit_num + (x(j)*(xc**(j-1)))
enddo
!write(*,*)fit_num,xval
summ = summ + (fit_num-xval)**2
enddo
close(19)
error = sqrt(summ/num_points)

end subroutine rms_error


!************************************************************
! Calculates error matrix
! Input : matrix A, B and X
! Output : matrix E
subroutine error_mat(matA,matB,matX,matE,msize,nsize1)

implicit none

real, intent(IN), allocatable :: matA(:,:), matB(:,:), matX(:,:)
integer, intent(IN) :: msize,nsize1

real, intent(OUT), allocatable :: matE(:,:)

allocate(matE(msize,nsize1))

matE = matmul(matA,matX) - matB

end subroutine error_mat


end module LinAl
