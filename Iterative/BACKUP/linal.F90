! /ams213a-bhishek/Project1/LU-Decomposition

! This module contains Linear Algebra routines as well as basic matrix properties calculation routines

! List of Subroutines/Functions:
! 1. readandallocatemat - Reads matrix from an input file and allocates
! 2. trace - This function calculates trace of a given input square matrix
! 3. euclid_norm - Calculates the euclidean norm of a given input vector
! 4. frob_norm - This function calculates the Frobenius norm of an input matrix
! 5. matrix_print - Prints any given input matrix in readable form
! 6. eliminate_gauss - Performs Gaussian elimination with partial pivoting
! 7. gauss_backsub - Backsubstitution to calculate X for the problem of type UX = Y
! 8. lu_decomp - Performs LU decomposition with partial pivoting
! 9. lu_backsub - Performs backsubtitution to solve LUX = PY
! 10. cholesky_fact - Performs Cholesky decomposition such that A=LL*
! 11. cholesky_back - Performs forward and back susbtitution on LL*x = b
! 12. qr_decomp - Performs QR decomposition on A to return orthogonal Q and upper triangular R
! 13. qr_ev - Calculates e-values and e-vectors of a symmetric matrix A using QR decomposition
! 14. rms_error - Calculates rms error for the fitted curve and input data 
! 15. create_lehmer - Creates Lehmer matrix of input dimension
! 16. plot_dat - Stores the fitted polynomial in a file "plot.dat" 
! 17. error_mat - Calculates error matrix E = AX - B 

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

!************************************************************************

! This function calculates the Forbenius norm of a given input matrix
! of dimension msize by nsize

real function frob_norm(matrix,msize,nsize)

implicit none

real, dimension(msize,nsize), intent(IN) :: matrix
integer, intent(IN) :: msize,nsize

integer i,j

frob_norm = 0.

do i=1,msize
   do j=1,nsize
      frob_norm = frob_norm + (abs(matrix(i,j)))**2
   enddo
enddo

frob_norm = sqrt(frob_norm)

end function frob_norm

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

   if (abs(matA(k,k)) <= epsilon(0.)) then
      write(*,*)'Diagonal Element 0 !!!'
      stop
   endif


   do i=k+1,msize
! Calculating multiplying factors
      mult = matA(i,k)/matA(k,k)
      matA(i,k) = 0.
      matA(i,k+1:msize) = matA(i,k+1:msize) - mult*matA(k,k+1:msize)
      matB(i,1:nsize1) = matB(i,1:nsize1) - mult*matB(k,1:nsize1)
   enddo     
enddo

det_A = 1.0

do i=1,nsize
   do j=1,msize
      det_A = det_A * matA(i,i)
   enddo
enddo

if (abs(det_A) <= epsilon(0.)) then
   ising = .TRUE.
else
   ising = .FALSE.
endif

end subroutine eliminate_gauss

!******************************************************************
! This subroutine performs backsubsitution to solve for x in AX=B
! Input - A, B, msize, nsize1 (number of RHS vectors), ising
!         ising indicates if the matri A is singular or not
! Output - X (solution matrix or vector depending on RHS B)


subroutine gauss_backsub(matA,matB,msize,nsize1,x,ising)

implicit none

real, intent(IN), dimension(msize,msize) :: matA
real, intent(IN), dimension(msize,nsize1) :: matB
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

   do i=msize,1,-1
      if (abs(matA(i,i)) <= epsilon(1.)) then
         write(*,*)'Dividing by zero is never a good thing! Try some other matrix!'
         stop
      endif

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
real, dimension(msize,msize), intent(INOUT) :: matA
integer, intent(IN) :: msize
logical, intent(IN) :: ising
integer, dimension(:), allocatable, intent(OUT) :: s

! Local Variables
real, dimension(:,:), allocatable :: tempa
integer, dimension(:), allocatable :: temps
integer i,j,k
real :: pivot
integer, dimension(1) :: pivot_loc

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


   if (abs(matA(j,j)) <= epsilon(0.)) then
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
   if (abs(matA(i,i)) <= epsilon(0.)) then
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

! Cholesky Decomposition
! Input & output : matrix A given as input and returns A again with LL* stored
! L* is stored as upper diagonal elements of A
! Input : Dimesnion of Matrix A, Ising: Logical tells whether matrix singular or not
! ising = FALSE means Matrix non-singular, ising = TRUE means Matrix singular

subroutine cholesky_fact(matA,msize,ising)

implicit none

! I/O Variables
real, intent(INOUT), allocatable :: matA(:,:)
integer, intent(IN) :: msize
logical, intent(IN) :: ising

! Local variables
integer :: i,j,k

if (ising .eqv. .FALSE.) then
   do j=1,msize
! New diagonal elements calculation
      do k=1,j-1
         matA(j,j) = matA(j,j) - matA(j,k)*matA(j,k)
      enddo

   if (matA(j,j) < epsilon(0.)) then
      write(*,*)'Matrix is not positive definite! Diagonal element negative'
      stop
   else
      matA(j,j) = sqrt(matA(j,j))
   endif

! Below diagonal elements calcuations
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

! Cholesky Backsubstitution
! Input : Matrix A, Vector b, dimension of matrix A and length of vector b 
! Output : Vector x and y

subroutine cholesky_back(matA,b,msize,x,y)

implicit none

! I/O variables
real, intent(IN), allocatable :: matA(:,:),b(:)
integer, intent(IN) :: msize
real, intent(OUT), allocatable :: x(:),y(:)

! Local variables
integer :: i,j,k
real :: summ

allocate (x(msize),y(msize))

! Fowrawd substitution solving Ly = b
do i=1,msize
   summ = b(i)
   do j=1,i-1
      summ = summ - y(j)*matA(i,j)
   enddo
   y(i) = summ/matA(i,i)
enddo

! Backward substitution solving L*x = y
do i=msize,1,-1
   if (abs(matA(i,i)) <= epsilon(0.)) then
       write(*,*)"Errrr! Diagonal element zero!"
       stop
   endif
   do k=i+1,msize
      y(i) = y(i)- matA(k,i)*x(k)
   enddo
   x(i) = y(i)/matA(i,i)
enddo

end subroutine cholesky_back

!**********************************************************
! This subroutine performs QR factorization on A
! Q is orthogonal matrix and R is upper triangular

! Input/Output - matrix A as input and returns R stored in A
! Input - msize, nsize dimensions of A
! Output - Q stored in q_mat

subroutine qr_decomp(matA,q_mat,msize,nsize)

! I/O Variables
real, dimension(msize,nsize), intent(INOUT) :: matA
real, intent(OUT), allocatable :: q_mat(:,:)
integer, intent(IN) :: msize, nsize

! Local variables
integer :: i,j
real :: summ,vec_norm
real, allocatable :: mult(:,:),s(:),v(:,:)
real, allocatable :: iden(:,:),house_holde(:,:),temp_h(:,:)

allocate (s(nsize),v(msize,nsize),mult(msize,msize))
allocate (temp_h(msize,msize))
allocate (iden(msize,msize),house_holde(msize,msize),q_mat(msize,msize))

! Initializing
   v = 0.
   s = 0.
   iden = 0.
   house_holde = 0.
   q_mat = 0.
   temp_h = 0.

do i=1,msize
   iden(i,i) = 1.
   temp_h(i,i) = 1.
enddo

do j=1,nsize
   mult = 0.

   summ = euclid_norm(matA(j:msize,j),msize-j+1)
   s(j) = sign(1.,matA(j,j))*summ
   v(j,j) = matA(j,j) + s(j)
   v(j+1:msize,j) = matA(j+1:msize,j)
  
! Vector v
   vec_norm =  euclid_norm(v(j:msize,j),msize-j+1)
   v(j:msize,j) = v(j:msize,j)/vec_norm

! Calculating V*Vt
   mult(j:msize,j:msize) = matmul(v(j:msize,j:nsize),transpose(v(j:msize,j:nsize)))

! Calculation of Q matrix
   if (j <= nsize) then
 
       house_holde = iden - 2.*mult
       q_mat = matmul(temp_h,house_holde)
       temp_h = q_mat
              
   endif

! Calculation of R stored in A
   mult(j:msize,j:nsize) = 2.*matmul(mult(j:msize,j:msize),matA(j:msize,j:nsize))
   matA(j:msize,j:nsize) = matA(j:msize,j:nsize) - mult(j:msize,j:nsize)

enddo

  deallocate (s,v,mult,temp_h,iden)

end subroutine qr_decomp

!*******************************************************************
! This matrix calcualted the eigenvalues and eigenvectors of an input symmetric matrix matA
! Input - matA, msize, nsize
! Output - matA, matV, q_mat, evalue

subroutine qr_ev(matA,msize,nsize,q_mat,matV,evalue)

! I/O Variables
real, dimension(msize,nsize), intent(INOUT) :: matA
real, allocatable, dimension(:,:), intent(OUT) :: matV,q_mat
real, dimension(:), allocatable, intent(OUT) :: evalue
integer, intent(IN) :: msize,nsize

! Local Variables
real :: tolerance,error
real, dimension(:), allocatable :: temp_e
integer :: i,j,iter_count
real, dimension(:,:), allocatable :: e_value1,Asave

allocate (e_value1(msize,nsize),evalue(msize),matV(msize,nsize),q_mat(msize,nsize))
allocate (Asave(msize,msize))

tolerance = 1e-5

do i=1,msize
   do j=1,nsize
      if (i == j) then
         matV(i,j) = 1.
      else
         matV(i,j) = 0.
      endif
   enddo
enddo

do i=1,msize
   evalue(i) = 0.
enddo

Asave = matA
error = 10.0
iter_count = 0

do while (error > tolerance)

   iter_count = iter_count + 1
   temp_e = evalue
   call qr_decomp(matA,q_mat,msize,nsize)
   matA = matmul(matA,q_mat)
   matV = matmul(matV,q_mat)
      e_value1 = matmul(matmul(transpose(q_mat),Asave),q_mat)

   do i=1,msize,1
      evalue(i) = e_value1(i,i)
   enddo 

   error = euclid_norm(evalue-temp_e,msize)

enddo

end subroutine qr_ev

!**********************************************************

! RMS error - (sqrt(summ(x-xi)**2)/num_points)
! Input : Vector a and its length, Vector a is just the coefficients of fitting function, filename for input data
! Output : rms error
! NOTE : a are the input coefficients of fitting function and xc read locally
!        are the points at which function is to be evaluated and further rms error calculated

subroutine rms_error(a,msize,error,filename1)

! I/O Variables
character*100, intent(IN) :: filename1
real, dimension(msize), intent(IN) :: a
integer, intent(IN) :: msize
real, intent(OUT) :: error

! Local Variables
real :: summ,fit_num,xc,xval
integer :: i,j,num_points,nsize

! Opening file which contains the points where fitting function to be evaluated
open(19,file=filename1,status='unknown')
read(19,*)num_points,nsize

! Reading points
summ = 0.0
do i=1,num_points
   read(19,*)xc,xval
   fit_num = 0.0

! Calculating value of fitting function at read points xc
   do j=msize,1,-1
      fit_num = fit_num + (a(j)*(xc**(j-1)))
   enddo
   summ = summ + (fit_num-xval)**2
enddo

close(19)

error = sqrt(summ/num_points)

end subroutine rms_error

!**********************************************************
! Creates Lehmer matrix
! Input - Dimension of matrix
! Output - Lehmer matrix

subroutine create_lehmer(mat,msize)

real, intent(OUT), allocatable :: mat(:,:)
integer, intent(IN) :: msize

integer i,j
real i1,j1

allocate (mat(msize,msize))

   do i=1,msize,1
      do j=1,msize,1
         i1 = real(i)
         j1 = real(j)
         mat(i,j) = min(i1,j1)/max(i1,j1)
      enddo
   enddo

end subroutine create_lehmer

!***********************************************************
! Print data for fitted curve in file
! Input - Coefficients of fitted polynomial and degree of polynomial
! Output - Created a file plot.dat containing the fitted curve data in [0,1)
subroutine plot_data(a,msize)

! I/O variable
real, dimension(msize), intent(IN) :: a
integer, intent(IN) :: msize

! Local variables
integer :: i,j
real :: x,y

open(11,file='plot.dat',status='unknown')

do i=1,201,1
   y = 0.
   x = (i-1)*0.005
   do j=msize,1,-1
      y = y + (a(j)*(x**(j-1)))
   enddo
   write(11,*)x,y
enddo

close(11)
end subroutine plot_data

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

!*************************************************************

subroutine gauss_jacobi(A,b,msize,accuracy,x)

implicit none

real, dimension(msize,msize), intent(IN) :: A
real, dimension(msize), intent(IN) :: b
real, intent(IN) :: accuracy
integer, intent(IN) :: msize

real, dimension(msize), intent(INOUT) :: x(:)

real, allocatable :: r(:,:)
real, allocatable :: err_vec(:),y(:),w1(:)

real :: error
integer :: i,j,k

open(1,file='Error_Iteration_J.dat',status='unknown')

allocate (r(msize,msize),w1(msize),err_vec(msize),y(msize))

k=0

r = A

do i=1,msize
          r(i,i) = 0.
enddo

do i=1,msize
   err_vec(i) =  b(i) - dot_product(A(i,:),x(:))
enddo

error = euclid_norm(err_vec,msize)

do while (error > accuracy)
k = k + 1
   do i=1,msize
      y(i) = (b(i) - dot_product(r(i,:),x(:)))/A(i,i)
   enddo

   err_vec = y - x
   x = y
   error = euclid_norm(err_vec,msize)
   write(*,*)'Error at Iteration step : ', k, 'is', error

   write(1,*)k,error

enddo

   close(1)

end subroutine gauss_jacobi

!***************************************************************

subroutine gauss_seidel(A,b,msize,accuracy,x)

implicit none

real, dimension(msize,msize), intent(IN) :: A
real, dimension(msize), intent(IN) :: b
real, intent(IN) :: accuracy
integer, intent(IN) :: msize

real, dimension(msize), intent(INOUT) :: x(:)

real, allocatable :: r(:,:)
real, allocatable :: err_vec(:),y(:),w1(:)

real :: error
integer :: i,j,k

open(3,file='Error_Iteration_S.dat',status='unknown')

allocate (r(msize,msize),w1(msize),err_vec(msize),y(msize))

k=0

r = A

do i=1,msize
          r(i,i) = 0.
enddo

do i=1,msize
   err_vec(i) =  b(i) - dot_product(A(i,:),x(:))
enddo

error = euclid_norm(err_vec,msize)

do while (error > accuracy)
 
  y = x

k = k + 1
   do i=1,msize
      x(i) = (b(i) - dot_product(r(i,:),x(:)))/A(i,i)
   enddo

   err_vec = (x - y)

   error = euclid_norm(err_vec,msize)
   write(*,*)'Error at Iteration step : ', k, 'is', error
   write(3,*)k,error

enddo

   close(3)

end subroutine gauss_seidel

!************************************************************

subroutine conj_grad(A,b,msize,accuracy,x)

real, dimension(msize,msize), intent(IN) :: A
real, dimension(msize), intent(IN) :: b
real, dimension(msize), intent(INOUT) :: x
real, intent(IN) :: accuracy

integer, intent(IN) :: msize

real, dimension(:), allocatable :: y,p,err_vec
real :: error,alpha,beta
integer :: i,j,k

k = 0

allocate (p(msize),err_vec(msize),y(msize))

open(5,file='Error_Iteration_CJ.dat',status='unknown')

p = b

do i=1,msize
   err_vec(i) =  b(i) - dot_product(A(i,:),x(:))
enddo

error = euclid_norm(err_vec,msize)

do while (error > accuracy)
k = k + 1
do i=1,msize
y(i) = dot_product(A(i,:),p(:))
enddo
alpha = (dot_product(p,err_vec))/(dot_product(p,y))
x = x + alpha*p
err_vec = err_vec - alpha*y
error = euclid_norm(err_vec,msize)
beta = -(dot_product(err_vec,y))/(dot_product(p,y))
p = err_vec + beta*p

write(*,*)'Error at Iteration step : ', k, 'is', error
write(5,*)k,error

enddo

end subroutine conj_grad

!**************************************************************

subroutine smart_cg(A,b,msize,accuracy,x)

real, dimension(msize,msize), intent(IN) :: A
real, dimension(msize), intent(IN) :: b
real, dimension(msize), intent(INOUT) :: x
real, intent(IN) :: accuracy

integer, intent(IN) :: msize

real, dimension(:), allocatable :: y,p,err_vec
real :: error,alpha,beta,tempe
integer :: i,j,k

k = 0

allocate (p(msize),err_vec(msize),y(msize))

open(5,file='Error_Iteration_CJ.dat',status='unknown')

p = b

do i=1,msize
   err_vec(i) =  b(i) - dot_product(A(i,:),x(:))
enddo

error = euclid_norm(err_vec,msize)
tempe = error
do while (error > accuracy)
k = k + 1
do i=1,msize
y(i) = dot_product(A(i,:),p(:))
enddo
alpha = (error**2)/(dot_product(p,y))
x = x + alpha*p
err_vec = err_vec - alpha*y
error = euclid_norm(err_vec,msize)
beta = (error**2)/(tempe**2)
p = err_vec + beta*p
tempe = error

write(*,*)'Error at Iteration step : ', k, 'is', error
write(5,*)k,error

enddo

end subroutine smart_cg

end module LinAl
