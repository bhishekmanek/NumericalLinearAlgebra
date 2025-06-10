program error

real valuem(101),valuel(101),summ
integer m,i

open(1,file='EV_M100.dat',status='unknown')
open(2,file='EV_L100.dat',status='unknown')

do i=1,100
read(1,*)m,valuem(i)
read(2,*)m,valuel(i)
enddo

summ=0.

do i=1,100
summ = summ + (abs(valuem(i)-valuel(i)))**2
enddo

write(*,*)'Error is', sqrt(summ/100.)

end program error
