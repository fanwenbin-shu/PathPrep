program main
implicit none
integer Natoms, i, j, k, window
real(8), allocatable :: coord(:,:), coordNew(:,:), eleMass(:)
character(4), allocatable :: eleName(:)
logical alive

! read coordinates from TS.xyz ! angstrom
open(999, file='TS.xyz', status='old', action='read') ! xyz -> 999
read(999, *) Natoms
allocate(coord(3, Natoms), eleName(Natoms), eleMass(Natoms))
read(999, *)
do i = 1, Natoms
read(999, *) eleName(i), coord(:, i)
!# write(*,'(3F11.6)') coord(:,i)
!# write(*,*) eleName(i)
end do
close(999)
! Get mass
call eleName2Mass(Natoms, eleName, eleMass)

! PES initialization
call init

inquire(file='spline.xyz', exist=alive)
if(.not. alive) then
write(*,*) 'Perform restricted optimization...'
call opt(Natoms, coord, eleName, eleMass)
end if

inquire(file='umbrella.xyz', exist=alive)
if(.not. alive) then
write(*,*) 'Perform interpolation...'
call intp(Natoms)
end if

!call Gibbs(Natoms, eleMass)

end program

subroutine Gibbs(Natoms, eleMass)
implicit none
integer Natoms, i, j, k, l, fileStatus
real(8) coord(3, Natoms), hessian(3*Natoms,3*Natoms), freq(3*Natoms)
real(8) xi, eleMass(Natoms), ZPE, S(15) ! 100K - 1500K
real(8) energy
character(4) :: eleName(Natoms)

open(427, file='gibbs.txt', status='unknown', action='write') ! Gibbs -> gbs -> 427

open(862, file='umbrella.xyz', status='unknown', action='read') !spline -> spl -> 775
do while (.true.)
read(862, *, iostat=fileStatus) Natoms
if (fileStatus /= 0) exit
read(862, *) xi ! title line
do i = 1, Natoms
read(862, *) eleName(i), coord(:, i)
end do

! ZPE
call freqCalc(Natoms, coord, eleMass, freq)
ZPE = sum(freq) * 0.00285914353799542d0 * 0.5d0 ! cm-1 to kcal/mol
do k = 1, 15
call SCalc(Natoms, freq*100d0*2.99792458d8, k*100d0, S(k))
end do
call potReal(Natoms, coord, energy)
write(427, '(20F11.6)') xi, energy, energy + S + ZPE

end do
close(775)

end subroutine

subroutine SCalc(Natoms, freq, T, S)
implicit none
integer Natoms, i
real(8) freq(Natoms * 3), T, S
real(8), parameter :: h = 6.62607015d-34, kB = 1.380649d-23

S = 0d0

do i = 1, 3*Natoms
if (freq(i) .gt. 1) then
S = S + &
& (h * freq(i)) / (kB * T * (dexp(h * freq(i) / (kB * T)) - 1)) - &
& dlog(1 - dexp(-h * freq(i) / (kB * T)))
!#write(*,*) dlog(1 - dexp(-h * freq(i) / (kB * T)))
end if
end do
S = S * 8.314d0 / 1000d0 * 4.184d0
! * R, J to kJ, kJ to kcal

end subroutine

! Interpolation for xi
subroutine intp(Natoms)
implicit none
integer Natoms, Nframe, i, j, k, l, fileStatus
real(8), allocatable :: coord(:,:,:), xiGiven(:)
real(8), allocatable :: coeff(:,:,:), xiList(:), coordOut(:,:,:)
real(8) xi, energy
integer, parameter :: Nwindow = 221
character(4) eleName(Natoms)
real(8) coord1D(Natoms*3), rr, mc(3), R, x1, x2, y1, y2, gam, test(3)

! read all geoms for interpolation
write(*,*) 'Reading geometries from `spline.xyz`...'
open(775, file='spline.xyz', status='old', action='read')
Nframe = 0
do while (.true.)
read(775, *, iostat=fileStatus) Natoms
if (fileStatus /= 0) exit
read(775, *)
do i = 1, Natoms
read(775, *)
end do
Nframe = Nframe + 1
end do
close(775)

allocate(coord(3, Natoms, Nframe), xiGiven(Nframe))
open(775, file='spline.xyz', status='old', action='read')
do i = 1, Nframe
read(775, *) Natoms
read(775, *) xiGiven(i)
do j = 1, Natoms
read(775, *) eleName(j), coord(:,j,i)
end do
end do
close(775)
write(*,'(A,I6,A)') ' Read ', Nframe, ' frames! '

! spline
! fit first
allocate(coeff(3, Natoms, Nframe))
do i = 1, Natoms
do j = 1, 3
call splineFit(xiGiven, coord(j, i, :), Nframe, coeff(j, i, :))
end do
end do
write(*,*) 'Spline fitting completed! '

! generate xi list
allocate(xiList(Nwindow)) ! -0.05 ~ 1.05, 111 points in sum
xiList(1) = -0.05d0
do i = 2, Nwindow
xiList(i) = xiList(i-1) + 0.005d0
end do
!#write(*,*) xiList

! interpolate new points
allocate(coordOut(3, Natoms, Nwindow))
do k = 1, Nwindow
do i = 1, Natoms
do j = 1, 3
call splineIntp(xiGiven, coord(j, i, :), coeff(j, i, :), Nframe, xiList(k), coordOut(j, i, k))
end do
end do
end do
write(*,*) 'Interpolation completed! '

! Output interpolated coordinates
open(688, file='umberlla.xyz', status='unknown', action='write') ! out->688
open(736, file='out_ene-xi.txt', status='unknown', action='write') ! spline energy -> sen -> 736
open(263, file='umbrella_configurations', status='unknown', action='write') !configuration -> cof -> 263
open(348, file='jacobi.txt', status='unknown', action='write') ! distance -> dis -> 348
do k = 1, Nwindow
! write .xyz
write(688, *) Natoms
call getXi(Natoms, coordOut(:,:,k), xi)
call potReal(Natoms, coordOut(:,:,k), energy)
write(688, '(F10.5)') xi
do j = 1, Natoms
write(688,'(A,3F15.8)') eleName(j), coordOut(:,j,k)
end do
! write xi and energy
write(736, '(3E18.8)') xi, energy * 627.509474063056d0, xi - xiList(k) ! to kcal/mol, difference
! write umbrella configuration in RPMDrate format
write(263, *) 
write(263, '(A5,F7.4)') 'xi = ', xi
do j = 1, Natoms
write(263,'(A5,3F12.6)') eleName(j), coordOut(:,j,k) / 0.529177210903d0
end do

end do
close(688)
close(736)
close(263)

end subroutine

subroutine opt(Natoms, coord, eleName, eleMass)
implicit none
integer Natoms, i, j, k, l
real(8) coord(3, Natoms), coordBkup(3, Natoms), grad(3, Natoms)
real(8) xi, energy, anchor, eleMass(Natoms)
character(4) eleName(Natoms)
real(8), parameter :: space = 0.005d0, maxXi = 1.1d0, minXi = -0.1d0, alpha = 0.05d0
integer, parameter :: step = 1000
integer window, windowR, windowL, c ! count
real(8), allocatable :: xiHis(:,:), energyHis(:,:), diffHis(:,:)
real(8), allocatable :: optHis(:,:,:,:), coordHis(:,:,:), xiList(:)

window = ceiling((maxXi - minXi) / space)
windowR = ceiling((maxXi - 1d0) / space)
!#write(*,*) 'windowR', window
allocate(xiHis(step, window), energyHis(step, window), diffHis(step, window))
allocate(optHis(3, Natoms, step, window), coordHis(3, Natoms, window), xiList(window))

coordBkup = coord

! forward or right
c = window - windowR + 1 ! record energy etc. from here
do anchor = 1d0, maxXi, space
!#write(*,*) 'c', c
    do i = 1, step
    call gradCalc(Natoms, coord, grad, anchor)
    coord = coord - alpha * grad
    !call writeXYZ(Natoms, eleName, coord, 688, xi, energy)
    call getXi(Natoms, coord, xi)
    call potReal(Natoms, coord, energy)
    optHis(:,:,i,c) = coord
    xiHis(i,c) = xi
    energyHis(i,c) = energy
    diffHis(i,c) = xi - anchor
    xiList(c) = anchor
    end do
!call writeXYZ(Natoms, eleName, coord, 775, xi, energy)
coordHis(:,:,c) = coord
write(*,'(3(A,E18.8))') 'Anchor xi: ', anchor, &
                       &',  xi:', xi,&
                       &',  diff:', (xi - anchor) !/ 0.005d0
c = c + 1 ! count + 1

end do

! backword or left
coord = coordBkup ! restore the initial geom
c = window - windowR
do anchor = 1d0-space, minXi-0.001d0, -space ! 0.001 means the lower bound small than minXi
!#write(*,*) 'c', c
    do i = 1, step
    call gradCalc(Natoms, coord, grad, anchor)
    coord = coord - alpha * grad
    !call writeXYZ(Natoms, eleName, coord, 688, xi, energy)
    call getXi(Natoms, coord, xi)
    call potReal(Natoms, coord, energy)
    optHis(:,:,i,c) = coord
    xiHis(i,c) = xi
    energyHis(i,c) = energy
    diffHis(i,c) = xi - anchor
    xiList(c) = anchor
    end do
!call writeXYZ(Natoms, eleName, coord, 775, xi, energy)
coordHis(:,:,c) = coord
write(*,'(3(A,E18.8))') 'Anchor xi: ', anchor, &
                       &',  xi:', xi,&
                       &',  diff:', (xi - anchor) !/ 0.005d0
c = c - 1
end do

!#write(*,'(6F10.4)') xiList

! write all intermediate information
! all geometries in optimization
write(*,*) 'Writing all geometries in optimization...'
open(688, file='opt.xyz', status='unknown', action='write') !opt -> 688
do i = 1, window
do j = 1, step
!call writeXYZ(Natoms, eleName, optHis(:,:,j,i), 688)
end do
end do
close(688)

! last geom
write(*,*) 'Writing last geometries in each window...'
open(775, file='spline.xyz', status='unknown', action='write') !spl -> 775
open(636, file='opt_ene-xi.txt', status='unknown', action='write') ! opt energy -> oen -> 636
do k = 1, window
write(775, *) Natoms
call getXi(Natoms, coordHis(:,:,k), xi)
call potReal(Natoms, coordHis(:,:,k), energy)
write(775, '(F10.5)') xi
write(636, '(3E18.8)') xi, energy * 627.509474063056d0 ! to kcal/mol, difference
do j = 1, Natoms
write(775,'(A,3F15.8)') eleName(j), coordHis(:,j,k)
end do
end do
close(775)

! xi, energy, diff in opt
write(*,*) 'Writing evolution of xi, energy in optimization...'
open(94,  file='opt_xi.txt', status='unknown', action='write') ! xi -> 94
open(363, file='opt_energy.txt', status='unknown', action='write') ! ene -> 363
open(343, file='opt_difference.txt', status='unknown', action='write') ! dif -> 343
do i = 1, step
write(94, '(I5, 10000E18.8)') i, xiHis(i, :)
write(363,'(I5, 10000E18.8)') i, energyHis(i, :)
write(343,'(I5, 10000E18.8)') i, diffHis(i, :)
end do
close(94)
close(363)
close(343)
open(94,  file='opt_conv_xi.txt', status='unknown', action='write') ! xi -> 94
open(363, file='opt_conv_energy.txt', status='unknown', action='write') ! ene -> 363
open(343, file='opt_conv_difference.txt', status='unknown', action='write') ! dif -> 343
do i = 2, step
write(94, '(I5, 10000E18.8)') i, abs(xiHis(i, :) - xiHis(i-1, :))
write(363,'(I5, 10000E18.8)') i, abs(energyHis(i, :) - energyHis(i-1, :))
write(343,'(I5, 10000E18.8)') i, abs(diffHis(i, :) - diffHis(i-1, :))
end do
close(94)
close(363)
close(343)

end subroutine

! Export coordinates to .xyz file
subroutine writeXYZ(Natoms, eleName, coord, fileID)!, xi, energy)
implicit none
integer Natoms, fileID, i
character(4) eleName(Natoms)
character(64) info
real(8) coord(3, Natoms), xi, energy

call getXi(3, coord, xi)
!call potReal(3, coord, energy)

write(fileID,'(I4)') Natoms
write(fileID,'(2F18.8)') xi!, energy
do i = 1, Natoms
write(fileID,'(A,3F15.8)') eleName(i), coord(:,i)
end do
end subroutine

! Calculate the gradient
subroutine gradCalc(Natoms, coord, grad, xit)
implicit none
integer Natoms
real(8) coord(3, Natoms), xit
real(8) grad(3, Natoms)
real(8) eneFor, eneBak, coordBackup ! forward, backward, backup
real(8) :: delta = 1d-5 ! displacement
integer i, j

do j = 1, Natoms
do i = 1, 3 ! x,y,z
! finite displacement
coordBackup = coord(i,j)
coord(i,j) = coord(i,j) + delta
call pot(Natoms, coord, eneFor, xit)
coord(i,j) = coord(i,j) - 2d0 * delta
call pot(Natoms, coord, eneBak, xit)
coord(i,j) = coordBackup
! calculate gradient
grad(i,j) = (eneFor - eneBak) / (2d0 * delta)
end do
!#write(*,'(3F11.6)') grad(:,j)
end do
end subroutine

! Get each mass from element name
subroutine eleName2Mass(Natoms, eleName, eleMass)
implicit none
integer Natoms
character(4) :: allName(20) = &
(/'H', 'He',   'Li',   'Be',   'B', &
 'C',  'N',    'O',    'F',    'Ne', &
 'Na', 'Mg',   'Al',   'Si',   'P', &
 'S',  'Cl',   'Ar',   'K',    'Ca'/)
real(8) :: allMass(20) = &
(/1.008d0,      4.002602d0, 6.94d0,       9.0121831d0,    10.81d0, &
 12.011d0,      14.007d0,   15.999d0,     18.998403163d0, 20.1797d0, &
 22.98976928d0, 24.305d0,   26.9815385d0, 28.085d0,       30.973761998d0, &
 32.06d0,       35.45d0,    39.948d0,     39.0983d0,      40.078d0/)
character(4) eleName(Natoms)
real(8), intent(out) :: eleMass(Natoms)
integer i, j

do i = 1, Natoms
do j = 1, 20
!# write(*,*) i, j, eleName(i), allName(j)
if (trim(eleName(i)) == trim(allName(j))) then
eleMass(i) = allMass(j)
exit
!else
!write(*, *) '[ERROR] Element not supported: ', eleName(i)
end if
end do
end do
!# write(*,*) 'Element mass: '
!# write(*,*) eleMass
end subroutine

! Shift centroid to origin
subroutine cleanGeom(Natoms, eleMass, coord)
integer Natoms
real(8) coord(3, Natoms), coord3(3, Natoms), cm(3) ! center of mass
real(8) eleMass(Natoms), sumMass(3), totalMass

!coord3 = reshape(coord, (/3, Natoms/))
!# coord3(1,:) = coord3(1,:) + 1d0
!# write(*,*) 'Coordintaes'
!# write(*,'(3F11.6)') coord3

! Find center of mass
sumMass = 0d0
do j = 1, Natoms
do i = 1, 3
sumMass(i) = sumMass(i) + coord3(i, j) * eleMass(j)
!# write(*,'(3F11.6)') coord3(i, j) * eleMass(j)
end do
end do
!# write(*,*) 'Sum of mass * r'
!# write(*,'(3F11.6)') sumMass


! total mass
totalMass = sum(eleMass)
!# write(*,*) 'Total mass: '
!# write(*,'(3F11.6)') totalMass

! center of mass
cm = sumMass / totalMass
!# write(*,*) 'Center of mass: '
!# write(*,'(3F11.6)') cm

! shift CM to origin
do j = 1, Natoms
coord3(:, j) = coord3(:, j) - cm
end do
!# write(*,*) 'Coordintaes shifted '
!# write(*,'(3F11.6)') coord3

!coord = reshape(coord3, (/Natoms/))

end subroutine

! Below two subroutines copied from `Fortran Receipts` Chapter 3. ! excerpted on 2019-Feb-03 from book
! Usage: 
!   call splineFit(x, y, <amount of points>, <returned coefficient list>) only once
!   call splineIntp(x, y, <coefficient list>, <amount of points>, <new x>, <returned y>)
SUBROUTINE splineFit(x,y,n,y2)
INTEGER n, NMAX
real(8) x(n), y(n), y2(n)
PARAMETER (NMAX=1000) 
! 给定矩阵 x(1:n) 和 y(1:n) 包括函数列表, 如 y_i = f(x_i), x1 < x2 < ... < xN, 
! 并且给定 1 点和 n 点的一阶导数值 yp1 和 ypn. 这个程序返回长度为 n 的矩阵 y2(1:n), 
! 包括插值函数在列表中 x_i 点的二阶导数. 
! 参数: NMAX 是参与的最大 n 值.
INTEGER i,k
REAL(8) p,qn,sig,un,u(NMAX)
!write(*,*) '下边界条件为自然' !2019-01-30 10:35:49 Fan, WENBIN 
y2(1)=0.
u(1)=0.
do i=2,n-1 ! 这是三对角算法的分解循环. y2 和 u 用来临时储存分解因子
    !write(*,*) x(i-1), x(i), x(i+1)
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    !write(*,*) 'sig', sig !2019-01-30 10:40:10 Fan, WENBIN
    p=sig*y2(i-1)+2.
    y2(i)=(sig-1.)/p
    u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    !write(*,*) 'u(i)', u(i) !2019-01-30 10:38:00 Fan, WENBIN
enddo

!write(*,*) '上边界条件为自然' !2019-01-30 10:36:06 Fan, WENBIN
qn=0.
un=0.
y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
do k=n-1,1,-1 ! 这是三对角算法的向后替换循环
    y2(k)=y2(k)*y2(k+1)+u(k)
enddo
return
END

SUBROUTINE splineIntp(xa,ya,y2a,n,x,y)
INTEGER n
real(8) x,y,xa(n),y2a(n),ya(n)
! 给定长度为 n 的矩阵 xa(1:n) 和 ya(1:n), 表示的是顺序的函数列表. 给定矩阵 y2a(1:n), 它是上面 spline 的输出. 给定 x, 这个程序返回三次样条插值 y.
INTEGER k,khi,klo
real(8) a,b,h
klo=1 ! 我们将使用二分法找到正确的位置
khi=n ! 如果随机的连续调用此程序, 那这样是最好的. 如果是有序的连续调用, 并且间隔很小, 那可以先存储之前 klo 和 khi 的值, 并测试它们是否适合在下一次调用时.
1	if (khi-klo.gt.1) then
    k=(khi+klo)/2
    if(xa(k).gt.x)then
        khi=k
    else
        klo=k
    endif
goto 1
endif ! 现在 klo 和 khi 包括输入值 x
h=xa(khi)-xa(klo)
if (h.eq.0.) pause 'bad xa input in splint' ! xa 必须不同
a=(xa(khi)-x)/h ! 现在计算三次样条多项式
b=(x-xa(klo))/h
y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
return
END

! Below subroutine copied from my code `ExactTS` ! written on 2020-Mar-10

! Vibrational analysis
subroutine vibAna(Natoms, eleMass, hessian, freq)
use lapack95, only: syev
implicit none
integer Natoms
real(8) hessian(3 * Natoms, 3 * Natoms), hessianMW(3 * Natoms, 3 * Natoms) ! mass-weighted
real(8) eleMass(Natoms), eigenvalue(3 * Natoms)
real(8) freq(3 * Natoms), freqReduced(3*Natoms - 6), wavenumber(3 * Natoms)
real(8) tmp ! temp value for sorting freq
integer i, j, m, n
real(8) amu2kg, a2m, au2J
real(8), parameter :: pi = 4.d0*atan(1.d0)
!integer fileID

! Convert to mass-weighted Hessian
do j = 0, Natoms - 1
do i = 0, Natoms - 1
! 9 elements for each atom ! 3 * 3 matrix, ie. xx, xy, xz, yx, yy, yz, zx, zy, zz
do n = 1, 3
do m = 1, 3
hessianMW(i*3+m, j*3+n) = hessian(i*3+m, j*3+n) / dsqrt(eleMass(i+1) * eleMass(j+1))
end do
end do
end do
end do

! Eigenvalue
call syev(hessianMW, eigenvalue)
!# write(*,'(6F11.6)') eigenvalue

! To frequency ! Below 13 lines were referred Sobereva's `Hess2freq`. 
amu2kg=1.66053878D-27
a2m=1D-10
au2J=4.35974434D-18
eigenvalue=eigenvalue*au2J/a2m**2/amu2kg !convert force constant from a.u. to SI
do i = 1, Natoms * 3
if (eigenvalue(i)<0) then
freq(i) = -dsqrt(abs(eigenvalue(i))) / (2 * pi)
else
freq(i) = dsqrt(eigenvalue(i)) / (2 * pi)
end if
end do

freq = freq/2.99792458D10
!write(fileID, '(300F11.2)') freq

! sort modes by absolute value
do i = 1, Natoms * 3
do j = i + 1, Natoms * 3
if (abs(freq(i))>abs(freq(j))) then
tmp=freq(i)
freq(i)=freq(j)
freq(j)=tmp
end if
end do
end do

freqReduced = freq(7:)

! ! sort freq by value
! do i = 1, Natoms * 3 - 6
! do j = i + 1, Natoms * 3 - 6
! if (freqReduced(i) .gt. freqReduced(j)) then
! tmp = freqReduced(i)
! freqReduced(i) = freqReduced(j)
! freqReduced(j) = tmp
! end if
! end do
! end do

! ! output freq
! write(*,'(A)') 'Freq(cm-1): '
! write(*,'(6F12.2)') freq(7:)

end subroutine

! Calculate the Hessian matrix
subroutine hessianCalc(Natoms, coord1D, hessian)
implicit none
integer Natoms
real(8) coord1D(3 * Natoms)
real(8) hessian(Natoms * 3, Natoms * 3)
real(8) energy2nd

integer i, j

do j = 1, Natoms * 3
do i = 1, Natoms * 3
call der2nd(Natoms, coord1D, i, j, energy2nd)
hessian(i, j) = energy2nd
!# write(*,*) hessian(i,j)
end do
end do

!# ! check Hessian: $H_{i,j} = H_{j,i}$ ! check passed! 
!# do j = 1, Natoms * 3
!# do i = 1, Natoms * 3
!# write(*,'(3F11.6)') hessian(i, j), hessian(j, i), hessian(i, j) - hessian(j, i)
!# end do
!# end do

end subroutine

! Calculate the 2nd derivation
subroutine der2nd(Natoms, coord1D, i, j, energy)
implicit none
integer Natoms, i, j
real(8) coord1D(3 * Natoms), energy
real(8) :: delta = 1d-5, delta4squared = 4d-10
real(8) e1, e2, e3, e4 ! See ref: William, Fortran Reciept, Page 182
real(8) bkup1, bkup2 ! backup energy

! e1: x + delta, y + delta
bkup1 = coord1D(i)
bkup2 = coord1D(j)
coord1D(i) = coord1D(i) + delta
coord1D(j) = coord1D(j) + delta
call TSenergy(Natoms, reshape(coord1D, (/3, Natoms/)), e1)
coord1D(i) = bkup1
coord1D(j) = bkup2

! e2: x + delta, y - delta
bkup1 = coord1D(i)
bkup2 = coord1D(j)
coord1D(i) = coord1D(i) + delta
coord1D(j) = coord1D(j) - delta
call TSenergy(Natoms, reshape(coord1D, (/3, Natoms/)), e2)
coord1D(i) = bkup1
coord1D(j) = bkup2

! e3: x - delta, y + delta
bkup1 = coord1D(i)
bkup2 = coord1D(j)
coord1D(i) = coord1D(i) - delta
coord1D(j) = coord1D(j) + delta
call TSenergy(Natoms, reshape(coord1D, (/3, Natoms/)), e3)
coord1D(i) = bkup1
coord1D(j) = bkup2

! e4: x - delta, y - delta
bkup1 = coord1D(i)
bkup2 = coord1D(j)
coord1D(i) = coord1D(i) - delta
coord1D(j) = coord1D(j) - delta
call TSenergy(Natoms, reshape(coord1D, (/3, Natoms/)), e4)
coord1D(i) = bkup1
coord1D(j) = bkup2

energy = (e1 - e2 - e3 + e4) / delta4squared

end subroutine

! inverse of any square matrix
subroutine pinv(N, A, pinvA, info)
use lapack95, only: gesvd
implicit none
integer N
real(8), intent(in) :: A(N,N)
real(8), intent(out) :: pinvA(N,N)
real(8) U(N,N), S(N), VT(N,N)!, WW(3)
real(8) S1(N,N), V(N,N), VS1(N,N)
real(8) work(N*N), lwork
integer, intent(out) :: info
integer i

call dgesvd('A', 'A', N, N, A, N, S, U, N, VT, N, work, N*N, info)

!#write(*,*) 'A'
!#write(*,'(6F11.4)') A

! pseudo inversion for S
S1 = 0d0
do i = 1, N
if (S(i) .gt. 1d-10) S1(i,i) = 1/S(i)
end do

! A-1 = V * S-1 * U^T
VS1 = matmul(transpose(VT), S1)
pinvA = matmul(VS1, transpose(U))

end subroutine

subroutine freqCalc(Natoms, coord, eleMass, freq)
implicit none
integer Natoms
real(8) coord(3, Natoms), coord1D(3*Natoms), freq(3*Natoms)
real(8) hessian(3*Natoms,3*Natoms), eleMass(Natoms)

call hessianCalc(Natoms, reshape(coord, (/3*Natoms/)), hessian)
call vibAna(Natoms, eleMass, hessian, freq)

end subroutine
