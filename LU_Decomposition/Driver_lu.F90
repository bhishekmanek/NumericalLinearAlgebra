! Main program to perfom LU decompositiong

program Driver_lu

use LinAl

implicit none

character*100 filename1,filename2
real, dimension(:,:), allocatable, save :: matA, Asave, matB1
real, dimension(:), save, allocatable :: vec,norm
integer, dimension(:), allocatable, save :: s 
real, dimension(:,:), allocatable, save :: matB, Bsave, matX, matE
integer, save :: msize, nsize, nsize1
real :: start,finish
logical :: ising

integer j,i

allocate(norm(nsize1))

call cpu_time(start)

filename1 = 'Amat.dat'
filename2 = 'Bmat.dat'
ising = .FALSE.

call readandallocatemat(matA,msize,nsize,filename1)
Asave = matA
call readandallocatemat(matB,msize,nsize1,filename2)
Bsave = matB
call lu_decomp(matA,msize,ising,s)
call lu_backsub(matA,msize,nsize1,matB,s,matX)
call error_mat(Asave,matB,matX,matE,msize,nsize1)
do i=1,nsize1

	norm(i) = euclid_norm(matE(:,i),msize)

enddo

!********************************************************************
! Priniting everything out
write(*,*)'The input matrix A is'
call matrix_print(Asave,msize,msize)
write(*,*)'The input matrix B is'
call matrix_print(Bsave,msize,nsize1)
write(*,*)'After performing LU decomposition, the LU matrix is:'
call matrix_print(matA,msize,msize)
write(*,*)'The permutation vector is:'
do i=1,msize

	write(*,*)s(i)

enddo
write(*,*)'Solving for X ...'
write(*,*)'X = '
call matrix_print(matX,msize,nsize1)
write(*,*)'The error matrix is:'
call matrix_print(matE,msize,nsize1)
write(*,*)'Norms of column vectors of the error matrix '
write(*,*)(norm(i),i=1,nsize1)
write(*,*)'Machine precision is',epsilon(0.)

call cpu_time(finish)
print*,'Time taken for the run is',finish-start, 'seconds'

stop
end
