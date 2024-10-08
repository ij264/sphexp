	SUBROUTINE gasdev(harvest)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: harvest
	REAL(SP) :: rsq,v1,v2
	REAL(SP), SAVE :: g
	LOGICAL, SAVE :: gaus_stored=.false.
	if (gaus_stored) then
		harvest=g
		gaus_stored=.false.
	else
		do
			call ran1(v1)
			call ran1(v2)
			v1=2.0_sp*v1-1.0_sp
			v2=2.0_sp*v2-1.0_sp
			rsq=v1**2+v2**2
			if (rsq > 0.0 .and. rsq < 1.0) exit
		end do
		rsq=sqrt(-2.0_sp*log(rsq)/rsq)
		harvest=v1*rsq
		g=v2*rsq
		gaus_stored=.true.
	end if
        END SUBROUTINE gasdev

