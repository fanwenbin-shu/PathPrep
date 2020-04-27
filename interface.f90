! 20200322 13:04:33 Wenbin, FAN @ SHU
! 
! An interface connecting your PES
! 
! Initializing your potential energy surface
subroutine init
call prepot
end subroutine

! Hock your PES here
! Unit: Hartree and Angstrom, please! 
subroutine potReal(Natoms, coord, energy)
implicit none
integer, intent(in) :: Natoms
real(8), intent(in) :: coord(3, Natoms)
real(8), intent(out) :: energy
real(8) CART, PENGYGS
INTEGER NATOM
PARAMETER (NATOM=25)
COMMON/USROCM/ PENGYGS
COMMON/USRICM/ CART(NATOM,3)

cart = transpose(coord) / 0.5291772d0
call FH2POT
energy = PENGYGS

end subroutine

! add potential energy and quadra potential
! `xit`: expected reaction coordinate
subroutine pot(Natoms, coord, energy, xit)
implicit none
integer, intent(in) :: Natoms
real(8), intent(in) :: coord(3, Natoms), xit ! expected xi
real(8), intent(out) :: energy
real(8) xi, quadra

call potReal(Natoms, coord, energy)
call getXi(Natoms, coord, xi)

quadra = 10d0 * (xi - xit) ** 2d0 ! force constant here
energy = energy + quadra
!write(*,*) energy, quadra

end subroutine

! define your reaction coordinate `xi` here. 
! here is the example for `FH2` system. 
subroutine getXi(Natoms, coord, xi)
implicit none
integer Natoms
real(8), intent(in) :: coord(3, Natoms)
real(8) dist
real(8), parameter :: rinf = 8d0
real(8) :: R, Rtmp, s0
real(8) r32, r31, r21, s(2), s1
real(8), intent(out) :: xi
integer i

do i = 1, 3
Rtmp = (coord(i, 3) * 0.5d0 + coord(i, 2) * 0.5d0 - coord(i, 1))
R = R + Rtmp * Rtmp
end do
R = dsqrt(R)
s0 = rinf - R

r32 = dist(Natoms, coord, 3, 2)
r31 = dist(Natoms, coord, 3, 1)
r21 = dist(Natoms, coord, 2, 1)

s(1) = r32 - 0.777931d0 - (r31 - 1.472932d0)
s(2) = r32 - 0.777931d0 - (r21 - 1.472932d0)

s1 = maxval(s)

xi = s0 / (s0 - s1)

end subroutine

function dist(Natoms, coord, No1, No2)
implicit none
integer Natoms, No1, No2, i
real(8) coord(3, Natoms), tmp
real(8) dist

do i = 1, 3
tmp = coord(i, No1) - coord(i, No2)
dist = dist + tmp * tmp
end do

dist = dsqrt(dist)

return
end function

! energy interface for gradients and second derivates
subroutine TSenergy(Natoms, coord, energy)
implicit none
integer, intent(in) :: Natoms
real(8), intent(in) :: coord(3, Natoms)
real(8), intent(out) :: energy

call potReal(Natoms, coord, energy)

end subroutine
