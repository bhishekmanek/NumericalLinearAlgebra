! Driver program for Cholesky decomposition

! Makes matrix A and vector b using given data from any input file
! Calls Linear Algebra routine to perform decomposition 
! Returns LL* stored in A and does backsubstituion to get solution vector
! Calculates error for solution and also rms error between input and fitted curve

program Driver_Cholesky

use linal

implicit none

character*100 filename1
real, dimension(:,:), save, allocatable :: a, matA, Asave, matB, Bsave, y
real, dimension(:), save, allocatable :: x,b, error,z
real, save :: error1,error_norm
integer, save :: msize,nsize
logical :: ising
integer i,j,num_points
real :: val

ising = .FALSE.

if (iargc() .ne. 1) then
   write(*,*)"Problem with your filename ! Enter correctly !"
   stop
endif

call getarg(1,filename1)

write(*,*)'Order of polynomial fitting to be done:'
read(*,*)msize
msize = msize + 1

! Reading data from input file

open(10,file=filename1,status='unknown')
read(10,*)num_points,nsize

allocate (x(msize),y(num_points,1),a(num_points,msize),matA(msize,msize))
allocate (matB(msize,1),b(msize),error(msize))

do i=1,num_points
   read(10,*)val,y(i,1)
   do j=1,msize
      a(i,j) = val**(j-1)
   enddo
enddo

close(10)

! Calculating normal equations
matA = matmul(transpose(a),a)
matB = matmul(transpose(a),y)

Asave = matA
Bsave = matB

b = matB(:,1)

call cholesky_fact(matA,msize,ising)
call cholesky_back(matA,b,msize,x,z)
call rms_error(x,msize,error1,filename1)
error = b - matmul(Asave,x)
error_norm = euclid_norm(error,msize)
call plot_data(x,msize)

! Printitng everything in terminal screen

write(*,*)'The normal matrices are'
do i=1,msize,1
   write(*,*)(Asave(i,j),j=1,msize),"|",b(i)
enddo
write(*,*)" "

write(*,*)'After Cholesky decomposition the Nomral matrix written in the form of a lower triangular matrix.'
write(*,*)'The upper diagonal elements are kept same as the elements of Normal matrix'
call matrix_print(matA,msize,msize)
write(*,*)" "

write(*,*)'Solving for x using Cholesky decomposition and backsubstitution...'
do i=1,msize
   write(*,*)x(i)
enddo
write(*,*)" "

write(*,*)'The error b-Ax is:'
do i=1,msize
   write(*,*)error(i)
enddo
write(*,*)" "
write(*,*)"||b-Ax|| - Euclid norm of error is ",error_norm
write(*,*)" "
write(*,*)'The rms error between fitted curve and data is',error1
write(*,*)" "
write(*,*)'The machine precision is ', epsilon(0.)

deallocate(x,y,a,matA,matB,b,error)

end program Driver_Cholesky
