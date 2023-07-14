module kindset

! Sets integer and real kind
implicit none

! Everything is really public, but the private scope is 
! imposed anyway.

private

public rk, ik

integer, parameter:: rk = selected_real_kind(15,307), &
					 ik = selected_int_kind(9)

end module kindset