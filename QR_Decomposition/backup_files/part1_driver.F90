program Driver_QR

use linal

implicit none

character*100 filename1,filename2
real, dimension(:,:), save, allocatable :: matA,Asave,v,q_mat
real, dimension(:), save, allocatable :: s
integer, save :: msize, nsize

integer :: i

fileName1 = 'matA.dat'
fileName2 = 'matB.dat'

call readandallocatemat(matA,msize,nsize,filename1)
Asave = matA
write(*,*)"Input matrix A:"
call matrix_print(Asave,msize,nsize)
!print*,msize,nsize
allocate (q_mat(msize,msize))
!call readandallocatemat(matB,msize,nsize1,filename2)
!Bsave = matB

call qr_decomp(matA,q_mat,msize,nsize)
do i=1,nsize
!print*,s(i)
enddo

end program Driver_QR
