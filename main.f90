program main

use kindset
use tools

implicit none

integer(kind=ik), parameter :: neps = 15
integer(kind=ik), parameter :: nk   = 100
integer(kind=ik), parameter :: nkg  = 100
integer(kind=ik), parameter :: kin = 5

integer(kind=ik) :: olgo, ie, iie, ix, iteration

integer(ik), allocatable:: gloc(:,:)

real(kind=rk), parameter    :: alpha = 0.265
real(kind=rk), parameter    :: nu = 0.6
real(kind=rk), parameter    :: delta = 0.067
real(kind=rk), parameter    :: betas = 0.96
real(kind=rk), parameter    :: qrf = betas
real(kind=rk), parameter    :: eta = 2.15
real(kind=rk), parameter    :: pid = 0.1
real(kind=rk), parameter    :: theta = 1.00
real(kind=rk), parameter    :: mu0in = pid
real(kind=rk), parameter    :: rhoeps = 0.653
real(kind=rk), parameter    :: stdinnoveps = 0.135
real(kind=rk), parameter    :: tauchene= 2.0
real(kind=rk), parameter    :: klow = 0.0
real(kind=rk), parameter    :: klow0 = 0.0
real(kind=rk), parameter    :: khigh = 6.0
real(kind=rk), parameter    :: khigh0 = 6.0
real(kind=rk), parameter    :: precision = 10.0**(-3.0)
real(kind=rk), parameter    :: p0 = 2.00

real(kind=rk):: wage, fval, flow, fhigh, kagg, kfagg, kagg2, output, iagg, cagg, nagg, tfp, producernum, &
                distance, starttime, endtime, ps, p, plow, phigh
           
real(kind=rk), allocatable  :: pie0(:), pie(:,:), evec(:), pievec(:), kgrid(:), kagrid(:), g(:,:), mu(:,:), tmu(:,:)
                               
character(30):: datafile, datafile1, datafile2, datestring, timestring

call cpu_time ( starttime )

datafile = 'workhorse'
call date_and_time(datestring, timestring)
olgo = len_trim(datafile)

datafile(olgo + 1:olgo+4) = timestring(1:4)
datafile1 = datafile(1:olgo+4)
datafile1(olgo+5:olgo+8) = '.txt'
olgo = len_trim(datafile)
datafile(olgo+1:olgo+3) = '.h5'
datafile2 = datafile(1:olgo)


! Idiosyncatic shock
allocate (evec(neps), pie(neps, neps), pievec(neps), pie0(neps))

call tauchen(0.0_rk, stdinnoveps, rhoeps, tauchene, neps, evec, pie)
evec = dexp(evec)    

call ergodicdist(neps, pie, 10.0_rk**(-6.0_rk), pie0)

! Grids on capital
allocate(kgrid(nk), kagrid(nkg))
call linspace(klow, khigh, nk, kgrid)
call linspace(klow0, khigh0, nkg, kagrid)

! Allocations must be before iterative steps
allocate(g(neps, nk), gloc(neps, nk), mu(neps, nk), tmu(neps, nk))

write(*,*) '+..........................................+'
write(*,*) '+   Starting Steady State Solution         +'
write(*,*) '+..........................................+'
write(*,*)

plow = p0*0.97_rk
phigh = p0*1.03_rk 

write(*,*) ' The plow evaluation '
write(*,*)

distance = 2.0*precision
iteration = 1
wage = eta/plow

call decisions(neps, nk, precision, wage, qrf, delta, alpha, nu, kgrid, evec, pie, g, gloc)

mu = 0.0
mu(1:neps, nkg/2) = pie0

call stationary(neps, evec, pie, pie0, pid, kin, mu0in, g, nkg, kagrid, precision, gloc, mu, &
                tmu, producernum, kagg, kfagg) 

call aggregation(neps, evec, alpha, nu, nkg, kagrid, tmu, wage, output, kagg2, nagg)

tfp = output/((kagg**alpha)*(nagg**nu))
iagg = kfagg  - (1.0 - delta)*kagg2
cagg = output - iagg 
ps = 1.0/cagg

fval = plow - ps
flow = fval

write(*,*)
write(*,*) '====================    results from plow evaluation    ==================== '
write(*,*)
write(*,'(1x, a, i4, a, 2f8.5, a, 4f8.4)') ' iteration ', iteration, ' (p,w) ', plow, wage, ' k, k2, n and output = ', &
    kagg, kagg2, nagg, output
write(*,'(1x,a,2f6.3,a,f8.5)') ' tfp, cagg = ', tfp, cagg, ' Tp = ', ps
write(*,'(1x,3(a,f9.4))') '    measure ', sum(tmu),  '   producers ', producernum
write(*,*) '====================                                       ==================== '
write(*,*); write(*,*)




write(*,*) ' The phigh evaluation '
write(*,*)

iteration = iteration + 1
wage = eta/phigh

call decisions(neps, nk, precision, wage, qrf, delta, alpha, nu, kgrid, evec, pie, g, gloc)

mu = 0.0
mu(1:neps, nkg/2) = pie0

call stationary(neps, evec, pie, pie0, pid, kin, mu0in, g, nkg, kagrid, precision, gloc, mu, &
                tmu, producernum, kagg, kfagg) 

call aggregation(neps, evec, alpha, nu, nkg, kagrid, tmu, wage, output, kagg2, nagg)

tfp = output/((kagg**alpha)*(nagg**nu))
iagg = kfagg  - (1.0 - delta)*kagg2
cagg = output - iagg 
ps = 1.0/cagg

fval = phigh - ps
fhigh = fval

write(*,*)
write(*,*) '====================    results from phigh evaluation    ==================== '
write(*,'(1x, a, i4, a, 2f8.4, a, 4f8.4)') ' iteration ', iteration, ' (p,w) ', phigh, wage, ' k, k2, n and output = ', &
    kagg, kagg2, nagg, output
write(*,'(1x,a,2f6.3,a,f8.5)') ' tfp, cagg = ', tfp, cagg, ' Tp = ', ps
write(*,'(1x,3(a,f9.4))') '    measure ', sum(tmu),  '   producers ', producernum
write(*,*) '====================                                       ==================== '
write(*,*); write(*,*)



write(*,*) 'Preparing for the p bisection'
write(*,*)

distance = 2.0*precision
do
    if (flow*fhigh.gt.0.0) then
        write(*,'(1x,a)') ' cannot bisect for steady state equilibrium p'
        write(*,'(1x,a,f8.4, f12.4)') ' plow  and flow = ', plow, flow
        write(*,'(1x,a,f8.4, f12.4)') ' phigh and fhigh = ', phigh, fhigh
        read(*,*)
        exit
    end if
    
    if (distance.lt.precision) then
        exit
    end if

iteration = iteration + 1    
p = (plow + phigh)/2.0
wage = eta/p

write(*,*); write(*,*)
write(*,'(1x, a, i4, a)') ' beginning p-iteration ', iteration, '.'
write(*,*)

g = 0.0

call decisions(neps, nk, precision, wage, qrf, delta, alpha, nu, kgrid, evec, pie, g, gloc)

mu = 0.0
mu(1:neps, nkg/2) = pie0

call stationary(neps, evec, pie, pie0, pid, kin, mu0in, g, nkg, kagrid, precision, gloc, mu, &
                tmu, producernum, kagg, kfagg)

call aggregation(neps, evec, alpha, nu, nkg, kagrid, tmu, wage, output, kagg2, nagg)

tfp = output/((kagg**alpha)*(nagg**nu))
iagg = kfagg - (1.0 - delta)*kagg2 
cagg = output - iagg
ps = 1.0/cagg

fval = p - ps
if (fval*fhigh.gt.0.0) then
    phigh = p
    fhigh = fval
else
    plow = p
    flow = fval
end if
distance = dabs(fval)

write(*,*)
write(*,*) '====================    results of current p evaluation    ==================== '
write(*,'(1x, a, i4, a, 2f8.5, a, 4f8.5)') ' iteration ', iteration, ' (p,wage) ', p, wage, ' k, k2, n and output = ', &
    kagg, kagg2, nagg, output
write(*,'(1x, a, f6.3, a, f6.3, a, f9.6)')  '    capital/output = ', kagg2/output, '    investment/capital = ', iagg/kagg2, '   excess demand: ', distance  
write(*,'(1x,a,2f6.3,a,f8.5)') ' tfp, cagg = ', tfp, cagg, ' Tp = ', ps
write(*,'(1x,3(a,f9.4))') '    measure ', sum(tmu),  '   producers ', producernum

write(*,*) '====================                                       ==================== '
write(*,*); write(*,*)


end do

call cpu_time ( endtime )
write(*,'(1x,a,f8.2,a)') '  elapsed time: ', endtime - starttime, ' seconds '

!!!!!!!!!!!!!!!!!!!!
! Store results    !
!!!!!!!!!!!!!!!!!!!!
!open(unit = 50, file = datafile2, status = 'unknown', action = 'write')
!
!do ix = 1, nkg
!    do ie = 1, neps
!        write(50,*) mu(ie, ix)
!    end do
!end do
!
!write(50,*) p, alpha, nu, delta, betas, theta, mu0in, kin, eta, pid, klow, khigh
!
!close (50)


write(*,*)
write(*,*)
write(*,*) '  Epsilon support, transition matrix, and ergodic dist. '
write(*,*)
do ie = 1, neps
    write(*,'(1x,f8.4)', advance = 'no') evec(ie)
end do
write(*,*)
write(*,*)

do ie = 1, neps
    do iie = 1, neps
        write(*,'(1x,f8.4)', advance = 'no') pie(ie, iie)
    end do
    write(*,'(1x,f8.4)') sum(pie(ie,1:neps))
end do
write(*,*)

do ie = 1, neps
    write(*,'(1x,f8.4)', advance = 'no') pie0(ie)
end do
write(*,'(1x,f8.4)') sum(pie0(1:neps))

write(*,*)
write(*,*)
write(*,*)
write(*,*) '  Capital choice. '
write(*,*)
do ie = 1, neps
    write(*,'(1x,f8.4)', advance = 'no') g(ie, nk)
end do
write(*,*); write(*,*)

write(*,*) '  Done! '
write(*,*); write(*,*)


end program main
