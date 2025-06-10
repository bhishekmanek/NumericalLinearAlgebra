! Driver routine for QR decomposition

! This program creates the matrix A from given data points and calls
! qr_decomp routine in module linal to perform QR decomposition

! Input - A, b
! Output - Q, R, x

   program Driver_QR

   use linal

   implicit none

! Variables to be saved and used locally

   logical :: ising
   character*100 filename1
   real, dimension(:,:), save, allocatable :: Bsave,matA,Asave,q_mat,acc_mat,ortho_mat,matB
   real, dimension(:,:), allocatable :: x_val, tempi
   real, dimension(:), allocatable :: error,b
   real, save :: fr_norm_acc,fr_norm_ortho,error1,error_norm
   integer, save :: msize, nsize
   real :: x
   integer :: i,j,nsize1

     write(*,*)'Order of polynomial fitting to be done:'
     read(*,*)nsize
     nsize = nsize + 1

     if (iargc() .ne. 1) then
        write(*,*)"Problem with your filename ! Enter correctly !"
        stop
     endif

     call getarg(1,filename1)

     open(10,file=filename1,status='unknown')
     read(10,*)msize,nsize1

     allocate (matA(msize,nsize),b(nsize))
     allocate (matB(msize,1),x_val(nsize,1))
     allocate (acc_mat(msize,nsize),ortho_mat(msize,msize),tempi(msize,msize),Bsave(msize,1))
     allocate (Asave(msize,nsize),error(nsize))

     do i=1,msize
        read(10,*)x,matB(i,1)
        do j=1,nsize
           matA(i,j) = x**(j-1)
        enddo
     enddo

     close(10)

     Asave = matA
     Bsave = matB

! Defining identity matrix to be used for QtQ - I
     do i=1,msize
        do j=1,msize
           if (i .eq. j) then
              tempi(i,j) = 1.
           else
              tempi(i,j) = 0.
           endif
        enddo
     enddo

     call qr_decomp(matA,q_mat,msize,nsize)
     acc_mat = Asave - matmul(q_mat,matA)
     fr_norm_acc = frob_norm(acc_mat,msize,nsize)
     ortho_mat = matmul(transpose(q_mat),q_mat) - tempi(:,:)
     fr_norm_ortho = frob_norm(ortho_mat,msize,msize)
     matB(:,1) = matmul(transpose(q_mat),matB(:,1))
     ising = .FALSE.
     call gauss_backsub(matA(1:nsize,1:nsize),matB(1:nsize,1),nsize,1,x_val,ising)
     b = matB(1:nsize,1)
     error = b - matmul(matA(1:nsize,1:nsize),x_val(1:nsize,1))
     error_norm = euclid_norm(error,nsize)
     call rms_error(x_val,nsize,error1,filename1)
     call plot_data(x_val,nsize)

! Printing to terminal screen starts from here
     write(*,*)"Working on QR decomposition ..." 
     write(*,*)" "
     write(*,*)"A - QR matrix"
     call matrix_print(acc_mat,msize,nsize)
     write(*,*)"Frobenius norm of A-QR is ",fr_norm_acc
     write(*,*)"Qt*Q - I matrix is "
     call matrix_print(ortho_mat,msize,msize)
     write(*,*)"Frobenius norm of Qt*Q - I is",fr_norm_ortho
     write(*,*)" "
     write(*,*)" Solution Vector x is: "
     call matrix_print(x_val,nsize,1)   
     write(*,*)'The error b-Ax is:'
     do i=1,nsize
        write(*,*)error(i)
     enddo 
     write(*,*)" "
     write(*,*)"||b-Ax|| - Euclid norm of error is ",error_norm
     write(*,*)" "
     write(*,*)'The rms error between fitted curve and data is',error1
     write(*,*)" "
     write(*,*)'The machine precision is ', epsilon(0.)

     deallocate (matA,matB,Asave,Bsave,ortho_mat,acc_mat,q_mat,tempi,error)

   end program Driver_QR
