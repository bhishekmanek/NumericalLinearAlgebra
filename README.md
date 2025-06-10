! This folder contains Linear Algebra routines as well as basic matrix properties calculation routines
! The main file is linal.F90 and all individual routines have been tested in separate folders with dedicated makefiles.

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
