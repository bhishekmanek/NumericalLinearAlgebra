program rms_error

real :: summ,x(21),rval(21),fit_value(21),fit_value5(21),error
integer :: i,msize,nsize

open(10,file='atkinson.dat',status='unknown')
read(10,*)msize,nsize
summ = 0.0
do i=1,21
read(10,*)x(i),rval(i)
fit_value(i) = 7.66867352*(x(i)**3)-11.1281977*(x(i)**2)+4.72584343*(x(i))+0.574662268
fit_value5(i)= 17.7211895*(x(i)**5)-46.5010490*(x(i)**4)+51.0233231*(x(i)**3)-28.0456734*(x(i)**2)+7.14999056*x(i)+0.511144936
summ = summ + (rval(i)-fit_value(i))**2
enddo

error=sqrt(summ/21)
print*,error
end program rms_error
