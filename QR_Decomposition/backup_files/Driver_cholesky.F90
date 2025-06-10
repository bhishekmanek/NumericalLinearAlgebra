program Driver

use linal

implicit none

character*100 filename1,filename2
real, dimension(:,:), save, allocatable :: a, matA, Asave, matB, Bsave, matE, y
real, dimension(:), save, allocatable :: vec,norm,x,b, error,z
real, save :: trace_mat,error1
integer, save :: msize,nsize,nsize1
logical :: ising
integer i,j,k,num_points
real :: start, finish

!allocate (error(msize))

ising = .FALSE.

write(*,*)'Order of polynomial fitting to be done:'
read(*,*)msize
msize = msize + 1
open(10,file='atkinson.dat',status='unknown')
read(10,*)num_points,nsize

allocate (x(num_points),y(num_points,1),a(num_points,msize),matA(msize,msize))
allocate (matB(msize,1),b(msize),error(msize))

do i=1,num_points
read(10,*)x(i),y(i,1)
do j=1,msize
a(i,j) = x(i)**(j-1)
enddo
enddo
close(10)
do i=1,num_points
!write(*,*)(a(i,j),j=1,msize)
enddo

! Calculating normal equations
matA = matmul(transpose(a),a)
matB = matmul(transpose(a),y)

open(11,file='matA_atk.dat',status='unknown')
open(12,file='matB_atk.dat',status='unknown')

write(11,*)msize,msize
write(12,*)msize,1
do i=1,msize,1
write(11,*)(matA(i,j),j=1,msize)
write(12,*)matB(i,1)
enddo

Asave = matA
Bsave = matB

filename1 = 'matA_atk.dat'
filename2 = 'matB_atk.dat'

!call readandallocatemat(matA,msize,nsize,filename1)
!Asave = matA
!call readandallocatemat(matB,msize,nsize1,filename2)
!Bsave = matB

b = matB(:,1)

call cholesky_fact(matA,msize,ising)
call cholesky_back(matA,b,msize,x,z)
call rms_error(x,msize,error1)
!error = matmul(Asave,x) - b

do i=1,msize
!write(*,*)error(i)
enddo

write(*,*)'The normal matrices are'
do i=1,msize,1
write(*,*)(Asave(i,j),j=1,msize),"|",b(i)
enddo
write(*,*)" "

write(*,*)'After Cholesky decomposition the Nomral matrix written in the form of a lower triangular matrix.'
write(*,*)'The upper diagonal elements are kept same as the elements of Normal matrix'
do i=1,msize
write(*,*)(matA(i,j),j=1,msize)
enddo
write(*,*)" "

write(*,*)'Solving for x using Cholesky decomposition and backsubstitution...'
do i=1,msize
write(*,*)x(i)
enddo
write(*,*)" "

error = matmul(Asave,x) - b
write(*,*)'The error matrix is:'
do i=1,msize
write(*,*)error(i)
enddo
write(*,*)" "
write(*,*)'The rms error is',error1
write(*,*)" "
write(*,*)'The machine precision is ', epsilon(0.)


end program Driver
