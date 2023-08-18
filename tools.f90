module tools

use kindset
use mod_splines
implicit none

private
public:: decisions, tauchen, stationary, aggregation, evolution, &
         linspace, logspace, gridlookup1, gridlookup, &
         profitfromk, ergodicdist

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                 !
!     decisions                   !
!                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine decisions(neps, nk, precision, wage, qrf, delta, alpha, nu, &
                     kgrid, evec, pie, kmax, gloc)

integer(ik):: neps, nk, ie, iik, iterationv, gloc(neps, nk)

real(rk):: distance, precision, &
           wage, qrf, delta, alpha, nu, &
           kgrid(nk), &
           evec(neps), pie(neps, neps), pievec(neps), epsilon, &
           kval, xval, &
           v(neps, nk), tv(neps, nk), &            
           ev(nk), vvec(neps), evmat(neps, nk), &
           kmax(neps, nk)
                          
intent(in):: neps, nk, precision, wage, qrf, delta, alpha, nu, kgrid, evec, pie
intent(out):: kmax, gloc


distance = 2.0_rk*precision
iterationv = 0_ik
v = 0_rk

do
   if (distance.lt.precision/(1.0_rk**(1.0_rk)).or.iterationv.ge.300_ik) then
       exit
   end if
        
   iterationv = iterationv + 1_ik
   
    do ie = 1, neps
    pievec(1:neps) = pie(ie, 1:neps)
        do iik = 1, nk
            vvec(1:neps) = v(1:neps, iik)
            evmat(ie, iik) = dot_product(pievec(1:neps), vvec(1:neps))
        end do
    end do
      
   do ie = 1, neps
    epsilon = evec(ie)
    ev(1:nk) = evmat(ie, 1:nk)

        do iik = 1, nk
            kval = kgrid(iik)
            xval = profitfromk(epsilon, alpha, nu, wage, kval) + (1-delta)*kval

            call maxgss(ev(1:nk), xval, nk, kgrid(1:nk), qrf, &
                        tv(ie, iik), kmax(ie, iik), gloc(ie, iik))
            
        end do
        
   end do

   distance = maxval(dabs(tv - v))
   if (iterationv.le.5_ik.or.modulo(iterationv,10).eq.0_ik) then
       
   write(*,'(1x,a,i3,a,f9.5 )') ' v iteration ', iterationv, '  normv = ', distance  
   
   end if
   
   v = tv
 
end do

end subroutine decisions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                       !
!   Golden section                                      !
!                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine maxgss(ev, xval, nk, kgrid, qrf, &
                  vcmax, kmax, gloc)
                      
type(spline) :: sp
integer(ik):: nk, ikv, ikpvec(1), gloc, i, n

real(rk):: ev(nk), qrf, &
           kgrid(nk), xval, &
           vcmax, kmax, &
           dval(nk), vc(nk), kfval, h, x, y

intent(in):: ev, xval, nk, kgrid, qrf
intent(out):: vcmax, kmax, gloc

kmax = 0.0_rk
vcmax = 0.0_rk
gloc = 0_ik

do ikv = 1, nk
    kfval = kgrid(ikv)
    dval(ikv) = xval - kfval        
    vc(ikv) = dval(ikv) + qrf*ev(ikv)
end do

!Cubic spline version
if (nk.gt.100) then
    !n = 50
    !h = (kgrid(size(kgrid))-kgrid(1))/(n-1)
    !print *, 'Cubic Spline Interpolation Demo'
    sp = spline(kgrid, vc)
    !print *, ""
    !print '(1x,a6,1x,a18,1x,a18)', "Index", "x", "y"
    !do i=1,n
    !    x = kgrid(1) + (i-1)*h
    !    y = sp%value(x)
    !    print '(1x,i6,1x,g18.11,1x,g18.6)', i, x, y
    !end do
    kmax = sp%extrema()
    gloc = sp%indexof(kmax)
    vcmax = sp%value(kmax)

    !print *, "Local Extrema"
    !print '(1x,a6,1x,a18,1x,a18)', "Index", "x", "y"
    !Print '(1x,i6,1x,g18.11,1x,g18.6)', gloc, kmax, vcmax

else
    ikpvec = maxloc(vc)
    gloc = ikpvec(1)

    vcmax = vc(gloc)
    kmax = kgrid(gloc)
    !Print '(1x,i6,1x,g18.11,1x,g18.6)', gloc, kmax, vcmax

endif

end subroutine maxgss


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                       !
!       Stationary Distribution         !
!                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine stationary(neps, evec, pie, pie0, pid, &
                      kin, mu0in, g, nkg, kagrid, precision, gloc, mu, &
                      tmu, producernum, kagg, kfagg)

integer(ik):: neps, nkg, kin, gloc(neps, nkg)
real(rk):: evec(neps), pie(neps, neps), pie0(neps), pid, &
           qval, zval, wage, delta, alpha, nu, &
           mu0in, g(neps, nkg), kagrid(nkg), precision, mu(neps, nkg), &
           tmu(neps, nkg), producernum

integer(ik):: iterationm
real(rk):: distancem, kagg, kfagg

intent(in)::neps, evec, pie, pie0, pid, &
            kin, mu0in, &
            g, nkg, kagrid, precision, gloc
intent(out):: tmu, producernum, kagg, kfagg
intent(inout):: mu

distancem = 2.0_rk*precision; iterationm = 0_ik

do
    if (distancem.lt.precision/(10.0_rk**(4.0_rk)).or.iterationm.gt.350_ik) then
        exit
    end if

    tmu = 0.0_rk

    call evolution(neps, nkg, kin, gloc, evec, pie, pie0, qval, zval, wage, delta, alpha, nu, kagrid, pid, mu0in, mu, g, &
                   tmu, producernum, kagg, kfagg)

    distancem = maxval(dabs(tmu - mu)); iterationm = iterationm + 1_ik
    
    mu = tmu
    
    write(*,'(1x,a, i4,4(a,f10.6))') '     stationary ', iterationm, '  entrants = ', mu0in, '    norm = ', distancem, '    kagg = ', kagg, &
                                         ' measure ', sum(tmu)
    
end do

end subroutine stationary


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                           !
!   Compute the distribution next period    !
!                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine evolution(neps, nkg, kin, gloc, evec, pie, pie0, qval, zval, wage, delta, alpha, nu, kagrid, pid, mu0in, mu0, ga, &
                     tmu, producernum, kagg, kfagg)
                 
integer(ik):: neps, nkg, kin, gloc(neps, nkg), &
              ie, ix, ikloc, ief
real(rk):: evec(neps), pie(neps, neps), pie0(neps), epsilon, &
           qval, zval, wage, &
           delta, alpha, nu, term, term1, term2, term3, &
           kstar(neps), &
           kagrid(nkg), &
           pid, mu0in, &
           mu0(neps, nkg), mu(neps, nkg), muval0, tmu(neps, nkg), &
           ga(neps, nkg), &
           kval, kvalf, &
           producernum, &
           kagg, kfagg

intent(in):: neps, nkg, kin, gloc, evec, pie, pie0, qval, zval, wage, delta, alpha, nu, kagrid, pid, mu0in, mu0, ga
intent(out):: tmu, producernum, kagg, kfagg

mu = mu0
producernum = 0.0_rk      
kagg = 0.0_rk
kfagg = 0.0_rk
term = 1 - alpha - nu

do ie = 1, neps
    epsilon = evec(ie)
    term1 = (((1/qval)-1+delta)/alpha)**((nu-1)/term)
    term2 = (zval*epsilon)**(1/term)
    term3 = (nu/wage)**(nu/term)
    kstar(ie) = term1*term2*term3
end do
do ie = 1, neps    
    epsilon = evec(ie)
    do ix = 1, nkg    

        muval0 = mu(ie, ix)
        if (muval0.gt.0.0_rk) then
        kval = kagrid(ix)
        kagg = kagg + kval*muval0
                        
        kvalf = ga(ie, ix)
        ikloc = gloc(ie, ix)

        producernum = producernum + muval0
            
        do ief = 1, neps         
            tmu(ief, ikloc) = tmu(ief,ikloc) + pie(ie, ief)*muval0
        end do
        
        kfagg = kfagg + kvalf*muval0
        end if
    end do          ! end of ix loop
end do              ! end of ie loop                   
                
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                               !
!       Exit and entrants       !
!                               !    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

tmu = (1.0_rk - pid)*tmu

kfagg = (1.0_rk - pid)*kfagg

ix = kin

do ie = 1, neps

    kvalf = ga(ie, ix)
    ikloc = gloc(ie, ix)
    
    do ief = 1, neps                
        tmu(ief, ikloc) = tmu(ief,ikloc) + pie(ie, ief)*pie0(ie)*mu0in
    end do
                            
    kfagg = kfagg + pie0(ie)*mu0in*kvalf
    
end do

end subroutine evolution


!!!!!!!!!!!!!!!!!!!!!!!!!
!                       !
!       profitfromk     !
!                       !
!!!!!!!!!!!!!!!!!!!!!!!!!

function profitfromk(epsilon, alpha, nu, wage, kp)

real(rk):: epsilon, alpha, nu, wage, kp
real(rk):: profitfromk, nuterm, wageterm, eval, kterm, output, profit

intent(in):: epsilon, alpha, nu, wage, kp

nuterm = 1.0_rk/(1.0_rk - nu)
wageterm = (nu/wage)**(nu*nuterm)
eval = epsilon**nuterm
kterm = kp**(alpha*nuterm)
output = eval*wageterm*kterm
profit = (1.0 - nu)*output

profitfromk = profit

end function profitfromk

!!!!!!!!!!!!!!!!!!!!!
!                   !
!   aggregation     !
!                   !
!!!!!!!!!!!!!!!!!!!!!

subroutine aggregation(neps, evec, alpha, nu, nkg, kagrid, mu, wage, &
                       output, kagg, nagg)

integer(ik):: neps, nkg, ie, ix

real(rk):: evec(neps), alpha, nu, kagrid(nkg), mu(neps, nkg), &
           wage, output, kagg, nagg
           
real(rk):: nuterm, wageterm, epsilon, eval, kval, mu0val, nval, yval

intent(in):: neps, evec, alpha, nu, nkg, kagrid, mu, wage
intent(out):: output, kagg, nagg
           
nuterm = 1.0_rk/(1.0_rk - nu)
wageterm = (nu/wage)**nuterm

nval = 0.0_rk
kval = 0.0_rk
yval = 0.0_rk
nagg = 0.0_rk
kagg = 0.0_rk
output = 0.0_rk
do ie = 1, neps
    epsilon = evec(ie)
    eval = epsilon**nuterm
        
    do ix = 1, nkg
        kval = kagrid(ix)
        
        nval = eval*wageterm*(kval**(alpha*nuterm))
        yval = epsilon*(kval**alpha)*(nval**nu)        
        mu0val = mu(ie, ix)
        nagg = nagg + nval*mu0val
        kagg = kagg + kval*mu0val
        output = output + yval*mu0val
        
    end do
end do
       
    
end subroutine aggregation



! ==========================================================
! ==========================================================
!                     FUNCTIONS
! ==========================================================
! ==========================================================


! ====================================================================
!         FUNCTION: gridlookup: find nearest (at or below) gridloc
! ====================================================================
subroutine gridlookup1(gridnum, xdist, xval, ixlow, iweight)

integer(ik):: gridnum, ixhigh, ixlow, ixplace 
real(rk):: xdist(gridnum), xval, iweight
intent(in):: gridnum, xdist, xval
intent(out):: ixlow, iweight

ixhigh = gridnum; ixlow = 1_ik

do

    if ((ixhigh	- ixlow).le.1_ik) then
	    exit
    end	if

    ixplace	= (ixhigh +	ixlow)/2_ik

    if (xdist(ixplace).ge.xval)	then
		    ixhigh = ixplace
    else
		    ixlow =	ixplace
    end	if

end	do

iweight = (xdist(ixlow + 1_ik) - xval)/(xdist(ixlow + 1_ik) - xdist(ixlow))

end subroutine gridlookup1

! ====================================================================
!         FUNCTION: gridlookup: find nearest (at or below) gridloc
! ====================================================================
pure function gridlookup(gridnum, xdist, xval)

integer(ik):: gridnum, ixhigh, ixlow, ixplace, gridlookup 
real(rk):: xdist(gridnum), xval
intent(in):: gridnum, xdist, xval

ixhigh = gridnum; ixlow = 1_ik

do

    if ((ixhigh	- ixlow).le.1_ik) then
	    exit
    end	if

    ixplace	= (ixhigh +	ixlow)/2_ik

    if (xdist(ixplace).ge.xval)	then
		    ixhigh = ixplace
    else
		    ixlow =	ixplace
    end	if

end	do

gridlookup = ixlow

end function gridlookup


! ==========================================================
! ==========================================================
!        FUNCTIONS THROWN IN FOR LATER USE - IN CASE...
! ==========================================================
! ==========================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!										!
!	TAUCHEN	and associated programs		!
!										!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! G. Tauchen (1986) 'Finite State Markov-Chain Approximations to 
! Univariate and Vector Autoregressions' Economics Leters 20: 177-181 

! Contains numerical integration of normal density to yield distribution.
subroutine tauchen(meaninnov, stdinnov, persistence, multiple, znum, zvec, pi)

integer(ik):: znum, j, gridsize, k
real(rk)::meaninnov, stdinnov, persistence, multiple, zvec(znum), pi(znum, znum), &
		  f1, f0, stdz, zlow, zhigh, meanz, w, z, lowerbound

intent(in):: meaninnov, stdinnov, persistence, multiple, znum
intent(out):: zvec, pi

! Set endpoints of stochastic variable grid at multiple m of the 
! standard deviation.  

stdz = stdinnov**2.0_rk
stdz = stdz/(1.0_rk - persistence**2.0_rk)
stdz = dsqrt(stdz)
meanz = meaninnov/(1.0_rk - persistence)
zlow = meanz - stdz*multiple
zhigh = meanz + stdz*multiple

lowerbound = meaninnov - stdinnov*dmax1(10.0_rk, 2.0_rk*multiple)
gridsize = 10000

call linspace(zlow, zhigh, znum, zvec)

pi = 0.0_rk

w = (zhigh - zlow)/dble(znum-1_ik)

! this is equations (3a) and (3b) from tauchen (1986)

do j = 1, znum
	
	z = zvec(1) - persistence*zvec(j)
	f1 = normal(z + w/2.0_rk, meaninnov, stdinnov, gridsize, lowerbound)
	pi(j,1) = f1
	
	do k = 2, znum - 1
		z = zvec(k) - persistence*zvec(j)
		f1 = normal(z + w/2.0_rk, meaninnov, stdinnov, gridsize, lowerbound)
		f0 = normal(z - w/2.0_rk, meaninnov, stdinnov, gridsize, lowerbound)
		pi(j,k) = f1 - f0
	end do

	z = zvec(znum) - persistence*zvec(j)
	f0 = normal(z - w/2.0_rk, meaninnov, stdinnov, gridsize, lowerbound)
	pi(j,znum) = 1.0_rk - f0

end do

end subroutine tauchen   

function normal(upperbound, mean, sd, gridsize, lowerbound)

! Program uses numerical integration with gridsize grid points to 
! compute the pr(lowerbound < x < upperbound) for a normal distri-
! bution with mean given and standard deviation of sd.  Numerical
! integration is based on a simple upper and lower darboux sum 
! method where the integrand is the normal density function.    

! Aubhik, bastion of stupidity, 18.04.2002

optional:: mean, sd, gridsize, lowerbound
integer(ik):: gridsize
real(rk):: lowerbound, upperbound, mean, sd, increment, normal
real(rk), allocatable:: x0dl(:), x1du(:), f0dl(:), f1du(:), f(:)

! Arguments for mean, sd, grid size for numerical integration and the 
! lower bound of the interval are optional.
if (present(mean)) then
	continue
else							! set mean to 0 if not specified	
	mean = 0.0_rk
end if

if(present(sd)) then
	continue
else							! set standard deviation to 1 is not specified
	sd = 1.0_rk
end if

if(present(lowerbound)) then 
	continue
else							! give probability from roughly - infinity	
	lowerbound = -10.0*sd + mean
end if


if(present(gridsize)) then
	continue					! use a 1000 points for numerical integration if not speficied
else
	gridsize = 100000_ik			
end if

! Grids for lower and upper darboux sums, remember Megan and math analysis at DRL.
increment = (upperbound - lowerbound)/dble(gridsize)

allocate(x0dl(gridsize), x1du(gridsize), f0dl(gridsize), f1du(gridsize), f(gridsize))

call linspace(lowerbound, upperbound - increment, gridsize, x0dl)
call linspace(lowerbound + increment, upperbound, gridsize, x1du)

! The normaldensity subroutine gives the value of the normal density at each grid point.
call normaldensity(gridsize, x0dl, mean, sd, f0dl)
call normaldensity(gridsize, x1du, mean, sd, f1du)

! Average the function values.
f = (f0dl + f1du)/2.0_rk

! Weight each one by the width of the rectangle of integration.
normal = sum(f*increment)

contains

subroutine normaldensity(numberofpoints, vectorofpoints, mean, sd, densities)

integer(ik):: numberofpoints
real(rk):: vectorofpoints(numberofpoints), densities(numberofpoints), mean, sd, &
		   variance, pi, coefficient, transcend(numberofpoints)

intent(in):: numberofpoints, vectorofpoints, mean, sd
intent(out):: densities

! From the matlab version of this program, the reference is apparently 
! page 244, chapter 4 of bopwerman and o'connell (1997) applied statistics

variance = sd**2.0_rk
pi = 2.0_rk*dasin(1.0_rk)
coefficient = variance*2.0_rk*pi
coefficient = 1.0_rk/dsqrt(coefficient)

! Note the vector-valued operations as a scalar, mean, is subtracted from each 
! element of a vector, then this resultant vector is divided by two scalars and
! finally used as a vector of exponents for the exponential function.

transcend = (vectorofpoints - mean)**2.0_rk
transcend = transcend/variance
transcend = transcend/2.0_rk
transcend = -1.0_rk*transcend
densities = coefficient*dexp(transcend)

end subroutine normaldensity

end function normal


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	   Ergodic distribution disteps using Markov Chain pieps        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ergodicdist(neps, pieps, precision, disteps)

! variables sent it
integer(ik)::neps
real(rk):: pieps(neps, neps), precision

! variables sent out
real(rk):: disteps(neps)

! variables used neither sent in nor out
real(rk):: diste0mat(neps,1), distemat(neps,1), piepst(neps,neps), distance

intent(in):: neps, pieps, precision
intent(out):: disteps

! Iterate to determine the ergodic distribution of shocks, diste
diste0mat = 0.0_rk; disteps = 0.0_rk; diste0mat = 1.0_rk/dble(neps); distemat = 0.0_rk
piepst = transpose(pieps)

distance = 2.0_rk*precision

do 
	distemat = matmul(piepst, diste0mat)
	distance = maxval(dabs(distemat - diste0mat))

	if (distance.lt.precision) then
		exit
	end if
	diste0mat = distemat
end do

disteps(1:neps) = distemat(1:neps,1)

end subroutine ergodicdist



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Sample random draws from a discrete distribution using E0Draw	!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine e0draw(ne, diste, plantsim, ei0)

integer(ik):: ne, plantsim, col, ei0(plantsim), plant
integer(ik), allocatable:: minplace(:)
real(rk):: diste(ne), pisum(ne), rv(plantsim), pid(ne), pisumtotal 

intent(in):: ne, diste, plantsim
intent(out):: ei0

pisum = 0.0_rk; pisum(1) = diste(1)

do col = 2, ne

	pisum(col) = pisum(col - 1) + diste(col)

end do

pisumtotal = pisum(ne) ! corrected from pisumtotal = sum(pisum)

call random_number(rv)
allocate(minplace(1))

do plant = 1, plantsim

	! subtract scalar from vector implies subtracting it from each element
	pid = pisum - rv(plant)*pisumtotal
	! Find least upper bound in pisum for rv, remember that minplace is an array!
	minplace = minloc(pid, mask = pid.gt.0.0_rk)
	ei0(plant) = max(minplace(1),1)
end do

deallocate(minplace)

end subroutine e0draw


! ====================================================================
!      FUNCTION:  stddev [compute sample standard deviation]
! ====================================================================
pure function stddev(vector)

! compute standard deviation of a vector, 
implicit none
integer(ik):: n
real(rk):: vector(:), mean, stddev, num
real(rk), allocatable:: vectoruse(:)

intent(in):: vector

n = size(vector); allocate(vectoruse(n))
num = dble(n)

mean = sum(vector)/num
vectoruse = (vector - mean)**2.0_rk

! num = dble(n - 1_ik) ! hogg and craig (1978), 4th edition, page 124.
stddev = sum(vectoruse)/num
stddev = sqrt(stddev)

end function stddev


!!!!!!!!!!!!!!!!!
!				!
!	COVARIANCE	!
!				!
!!!!!!!!!!!!!!!!!

function correlation(vec1, vec2, switch)

! Covariance or correlation between two vectors.  Note this program must recalculate the 
! standard deviations of vec1 and vec2 without using STDDEV.  Internal subroutines 
! cannot access other internal subroutines.  

optional:: switch
integer(ik):: n1, n2, switch
real(rk):: vec1(:), vec2(:), mean1, mean2, std1, std2, num, correlation !, sttdev
real(rk), allocatable:: vec(:)

intent(in):: vec1, vec2, switch

n1 = size(vec1)
n2 = size(vec2)

if (n1.ne.n2) then
	write(*,*) ' You bungler, I cannot calculate a correlation between different length series '	
	read(*,*)
	return
else 
	allocate(vec(n1))
end if

num = dble(n1)


mean1 = sum(vec1)/num
mean2 = sum(vec2)/num

vec = (vec1 - mean1)*(vec2 - mean2)

correlation = sum(vec)/num		! this is the covariance

if (present(switch)) then
	if (switch.eq.0_ik) then	! return covariance not correlation coefficient
		return
	end if	
end if

! Compute standard deviations
vec = (vec1 - mean1)**2.0_rk
std1 = sum(vec)/num
std1 = dsqrt(std1)

vec = (vec2 - mean2)**2.0_rk
std2 = sum(vec)/num
std2 = dsqrt(std2)

correlation = correlation/(std1*std2)	! correlation coefficient is the default

end function correlation


!!!!!!!!!!!!!!!!!
!	LinSpace	!
!!!!!!!!!!!!!!!!!

subroutine linspace(lb, ub, gridnum, x)

integer(ik):: gridnum, j1
real(rk):: lb, ub, x(gridnum), y(gridnum-2_ik)

intent(in):: lb, ub, gridnum
intent(out):: x

! This subroutine is written to reproduce the MATLAB 4.2 func-
! tion linspace.  It generates a vector which contains a grid 
! of n equally spaced elements between loweround and upper-
! bound.  Note this version produces a row vector.

do j1 = 1_ik, gridnum - 2_ik
	y(j1) = dble(j1)
end do

! Note that I am multiplying an n-2x1 vector, Y, by a scalar 
! (which implies multiplying each element of the vector by 
! this scalar) and then adding another scalar (which implies 
! adding this scalar to each element of the vector).

y = ((ub - lb)/(gridnum-1_ik))*y + lb

x(1_ik) = lb
x(2_ik:gridnum-1_ik) = dble(y)
x(gridnum) = ub

end subroutine linspace


!!!!!!!!!!!!!!!!!
!	LogSpace	!
!!!!!!!!!!!!!!!!!

subroutine logspace(lowerbound, upperbound, n, x)

integer:: n, j1
real(rk):: lowerbound, upperbound, term0, lb, ub
real(rk):: x(n), y(n-1)
intent(in):: lowerbound, upperbound, n
intent(out):: x

! This function departs from its MATLAB 5.3 counterpart in that 
! lower and upper bounds are the actual bounds on the logarithmically
! space vector.

term0 = dlog(10.0_rk)
lb = dlog(lowerbound)/term0
ub = dlog(upperbound)/term0

! this is lifted directly from the linspace subroutine. 
do j1 = 1, n - 2
	y(j1) = dble(j1)
end do

y = ((ub - lb)/(n-1))*y + lb

x(1) = lb
x(2:n-1) = y
x(n) = ub

! This is extremely cool, raising a scalar to a vector yields a vector.
x = 10.0_rk**x

end subroutine logspace


end module tools 
