program main_sphere_alfven
use iso_fortran_env, only: real64
use mhd_variational_integrator
implicit none
integer, parameter :: dp=real64
type(sphere_mesh) :: mesh
type(dec_complex) :: dec
real(dp), allocatable :: u(:),A(:),m(:)
real(dp) :: dt,cs,mu0,time,ekin,emag,eint,mass,flux,exacterr,residual
integer :: nsteps,nout,k,it
call build_icosphere(mesh)
call build_dec(mesh,dec)
allocate(u(mesh%ne),A(mesh%ne),m(mesh%nf))
call initialize_state(mesh,dec,u,A,m)
dt=1.0e-3_dp
cs=0.2_dp
mu0=1.0_dp
nsteps=100
nout=10
time=0.0_dp
write(*,'(a)') 'step time residual iterations exact_B_error total_flux total_mass kinetic magnetic internal'
do k=0,nsteps
    if (mod(k,nout)==0) then
        call diagnostics(mesh,dec,u,A,m,mu0,ekin,emag,eint,mass,flux,exacterr)
        write(*,'(i6,1x,es13.5,1x,es13.5,1x,i5,1x,es13.5,1x,es13.5,1x,es13.5,1x,es13.5,1x,es13.5,1x,es13.5)') &
            k,time,residual,it,exacterr,flux,mass,ekin,emag,eint
    end if
    if (k==nsteps) exit
    call advance_step(mesh,dec,u,A,m,dt,cs,mu0,residual,it)
    time=time+dt
end do
end program
