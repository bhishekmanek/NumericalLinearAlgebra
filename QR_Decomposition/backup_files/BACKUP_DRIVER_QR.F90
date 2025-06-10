program Driver_QR

use linal

implicit none

logical :: ising
character*100 filename1,filename2
real, dimension(:,:), save, allocatable :: tempi,Bsave, matA,Asave,v,q_mat,acc_mat,ortho_mat,matB,a,y
real, dimension(:,:), allocatable :: q_matt,val1,x_val,matA_new,matB_new
real, dimension(:), save, allocatable :: s,error,b,x,b_new
real, save :: fr_norm_acc,fr_norm_ortho,val2
integer, save :: msize, nsize,num_points

integer :: i,j,k,nsize1

!******************
write(*,*)'Order of polynomial fitting to be done:'
read(*,*)nsize
nsize = nsize + 1
open(10,file='atkinson.dat',status='unknown')
read(10,*)msize,nsize1
allocate (x(msize),y(msize,1),a(msize,nsize),matA(msize,nsize))
allocate (matB(msize,1),b(msize),error(msize),x_val(nsize,1),matA_new(nsize,nsize),matB_new(nsize,1))
allocate (acc_mat(msize,nsize),ortho_mat(msize,msize),tempi(msize,msize),Bsave(msize,1))
allocate (q_matt(msize,msize),val1(msize,msize),b_new(msize),Asave(msize,nsize))

do i=1,msize
read(10,*)x(i),y(i,1)
do j=1,nsize
a(i,j) = x(i)**(j-1)
enddo
enddo
close(10)

matA = a
matB = y

open(11,file='matA_atk.dat',status='unknown')
open(12,file='matB_atk.dat',status='unknown')

write(11,*)msize,msize
write(12,*)msize,1
do i=1,msize,1
write(11,*)(matA(i,j),j=1,nsize)
write(12,*)matB(i,1)
enddo
!*****************

Asave = matA
Bsave = matB

filename1 = 'matA_atk.dat'
filename2 = 'matB_atk.dat'

b = matB(:,1)

!call readandallocatemat(matA,msize,nsize,filename1)
!Asave = matA
!call readandallocatemat(matB,msize,nsize1,filename2)
!Bsave = matB

call qr_decomp(matA,q_mat,msize,nsize)

!write(*,*)"The input matrix A is "
!call matrix_print(Asave,msize,nsize)
write(*,*)"The matrix B is "
!call matrix_print(matB,msize,1)
write(*,*)"Working on QR decomposition ..." 
write(*,*)" "
!write(*,*)"After decomposition Q matrix is "
!call matrix_print(q_mat,msize,msize)
!write(*,*)"So Q transpos is"
q_matt = transpose(q_mat)
!call matrix_print(q_matt,msize,msize)
!write(*,*)"R matrix is "
!call matrix_print(matA,msize,nsize)
write(*,*)"A - QR matrix"
acc_mat = Asave - matmul(q_mat,matA)
call matrix_print(acc_mat,msize,nsize)
fr_norm_acc = frob_norm(acc_mat,msize,nsize)
write(*,*)"Frobenius norm of A-QR is ",fr_norm_acc

do i=1,msize
   do j=1,msize
      if (i .eq. j) then
         tempi(i,j) = 1.
      else
         tempi(i,j) = 0.
      endif
   enddo
enddo

ortho_mat = matmul(q_matt,q_mat) - tempi(:,:)
write(*,*)"Qt*Q - I matrix is "
call matrix_print(ortho_mat,msize,msize)
fr_norm_ortho = frob_norm(ortho_mat,msize,msize)
write(*,*)"Frobenius notm of Qt*Q - I is",fr_norm_ortho

matB(:,1) = matmul(q_matt,b)
!write(*,*)"Printing B:"
!call matrix_print(matB,msize,1)
ising = .FALSE.
matA_new = matA(1:nsize,1:nsize)
matB_new(1:nsize,1) = matB(1:nsize,1)
call gauss_backsub(matA_new,matB_new,nsize,1,x_val,ising)
write(*,*)" "
write(*,*)" Solution Vector x is: "
do i=1,nsize
print*,x_val(i,1)
enddo

end program Driver_QR
