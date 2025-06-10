program Driver_matrix

use linal

implicit none

character*100 filename1,filename2
real, dimension(:,:), save, allocatable :: x, matA, Asave, matB, Bsave, matE
real, dimension(:), save, allocatable :: vec,norm
real, save :: trace_mat
integer, save :: msize,nsize,nsize1
logical :: ising
integer i,j,k
real :: start, finish

call cpu_time(start)

allocate(norm(nsize1))

filename1 = 'Amat.dat'
filename2 = 'bmat.dat'

call readandallocatemat(matA,msize,nsize,filename1)
Asave = matA
call readandallocatemat(matB,msize,nsize1,filename2)
Bsave = matB
!*************************************************
trace_mat = trace(matA,msize)
write(*,*)'Input matrix is'
call matrix_print(Asave,msize,msize)
write(*,*)" "
write(*,*)'Trace of input matrix is ',trace_mat
!*************************************************

!************************************************
!vec = matA(:,1)
norm(1) = euclid_norm(matA(:,1),msize)
write(*,*)'Input vector is'
do i=1,msize
write(*,*)matA(i,1)
enddo
write(*,*)" "
write(*,*)"Euclidean norm of vector is ",norm(1)
!**********************************************

!call eliminate_gauss(matA,matB,nsize,msize,nsize1,ising)
!call gauss_backsub(matA,matB,msize,nsize1,x,ising)
!call error_mat(matA,matB,x,matE,msize,nsize1)
!do i=1,nsize1
!norm(i) = euclid_norm(matE(:,i),msize)
!enddo
!write(*,*)'The input matrix A is'
!call matrix_print(Asave,msize,msize)
!write(*,*)'The input matrix B is'
!call matrix_print(Bsave,msize,nsize1)
!write(*,*)'After pivoting, the Upper traingular matrix (U) is'
!call matrix_print(matA,msize,msize)
!write(*,*)'Solving for X ...'
!write(*,*)'X = '
!call matrix_print(x,msize,nsize1)
!write(*,*)'Matrix Y is'
!call matrix_print(matB,msize,nsize1)
!write(*,*)'Error matrix E = AX - B is'
!call matrix_print(matE,msize,nsize1)
!write(*,*)'Norms of column vectors of the error matrix '
!write(*,*)(norm(i),i=1,nsize1)
!write(*,*)'Machine precision is',epsilon(0.)

!call cpu_time(finish)
!print*,'Time taken for the run is',finish-start, 'seconds'

end program Driver_matrix
