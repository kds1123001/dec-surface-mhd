module mhd_variational_integrator
use iso_fortran_env, only: real64, int32
use, intrinsic :: ieee_arithmetic
implicit none
private
public :: sphere_mesh, dec_complex, build_icosphere, build_dec
public :: initialize_state, advance_step, diagnostics, destroy_mesh

integer, parameter :: dp=real64
real(dp), parameter :: pi=acos(-1.0_dp)

type sphere_mesh
    integer :: nv=0, ne=0, nf=0
    real(dp), allocatable :: xyz(:,:), edge_length(:), face_area(:)
    real(dp), allocatable :: face_normal(:,:), circumcenter(:,:)
    integer, allocatable :: edge(:,:), face(:,:), face_edge(:,:), face_sign(:,:)
    integer, allocatable :: edge_faces(:,:)
    real(dp), allocatable :: bary_vertex_area(:)
end type

type dec_complex
    integer :: nv=0, ne=0, nf=0
    integer(int32), allocatable :: d0(:,:), d1(:,:)
    real(dp), allocatable :: star0(:), star1_bary(:), star1_circ(:), star2(:)
end type

contains

subroutine build_icosphere(mesh)
type(sphere_mesh), intent(out) :: mesh
real(dp) :: phi, normv, tmp(3), c(3), a(3), b(3), n(3), cr(3)
integer :: i,j,k,f,e,v0,v1,sgn
integer :: baseface(3,20)
integer :: ev(2,30), ne
real(dp) :: verts(3,12)
phi=(1.0_dp+sqrt(5.0_dp))/2.0_dp
verts(:,1)=[-1.0_dp,phi,0.0_dp]
verts(:,2)=[1.0_dp,phi,0.0_dp]
verts(:,3)=[-1.0_dp,-phi,0.0_dp]
verts(:,4)=[1.0_dp,-phi,0.0_dp]
verts(:,5)=[0.0_dp,-1.0_dp,phi]
verts(:,6)=[0.0_dp,1.0_dp,phi]
verts(:,7)=[0.0_dp,-1.0_dp,-phi]
verts(:,8)=[0.0_dp,1.0_dp,-phi]
verts(:,9)=[phi,0.0_dp,-1.0_dp]
verts(:,10)=[phi,0.0_dp,1.0_dp]
verts(:,11)=[-phi,0.0_dp,-1.0_dp]
verts(:,12)=[-phi,0.0_dp,1.0_dp]
baseface(:,1)=[1,12,6]
baseface(:,2)=[1,6,2]
baseface(:,3)=[1,2,8]
baseface(:,4)=[1,8,11]
baseface(:,5)=[1,11,12]
baseface(:,6)=[2,6,10]
baseface(:,7)=[6,12,5]
baseface(:,8)=[12,11,3]
baseface(:,9)=[11,8,7]
baseface(:,10)=[8,2,9]
baseface(:,11)=[4,10,5]
baseface(:,12)=[4,5,3]
baseface(:,13)=[4,3,7]
baseface(:,14)=[4,7,9]
baseface(:,15)=[4,9,10]
baseface(:,16)=[5,10,6]
baseface(:,17)=[3,5,12]
baseface(:,18)=[7,3,11]
baseface(:,19)=[9,7,8]
baseface(:,20)=[10,9,2]
do i=1,12
    normv=sqrt(sum(verts(:,i)**2))
    verts(:,i)=verts(:,i)/normv
end do
mesh%nv=12
mesh%nf=20
mesh%ne=30
allocate(mesh%xyz(3,12),mesh%face(3,20))
mesh%xyz=verts
mesh%face=baseface
allocate(mesh%edge(2,30),mesh%face_edge(3,20),mesh%face_sign(3,20))
allocate(mesh%edge_faces(2,30),mesh%edge_length(30),mesh%face_area(20))
allocate(mesh%face_normal(3,20),mesh%circumcenter(3,20))
allocate(mesh%bary_vertex_area(12))
mesh%edge=0
mesh%edge_faces=0
ne=0
do f=1,20
    do j=1,3
        v0=mesh%face(j,f)
        v1=mesh%face(mod(j,3)+1,f)
        e=0
        do k=1,ne
            if ((mesh%edge(1,k)==min(v0,v1)).and.(mesh%edge(2,k)==max(v0,v1))) e=k
        end do
        if (e==0) then
            ne=ne+1
            e=ne
            mesh%edge(:,e)=[min(v0,v1),max(v0,v1)]
        end if
        mesh%face_edge(j,f)=e
        if (v0==mesh%edge(1,e)) then
            mesh%face_sign(j,f)=1
        else
            mesh%face_sign(j,f)=-1
        end if
        if (mesh%edge_faces(1,e)==0) then
            mesh%edge_faces(1,e)=f
        else
            mesh%edge_faces(2,e)=f
        end if
    end do
end do
do e=1,30
    mesh%edge_length(e)=sqrt(sum((mesh%xyz(:,mesh%edge(2,e))- &
        mesh%xyz(:,mesh%edge(1,e)))**2))
end do
mesh%bary_vertex_area=0.0_dp
do f=1,20
    a=mesh%xyz(:,mesh%face(1,f))
    b=mesh%xyz(:,mesh%face(2,f))
    c=mesh%xyz(:,mesh%face(3,f))
    cr=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3), &
        a(1)*b(2)-a(2)*b(1)]
    n=[(b(2)-a(2))*(c(3)-a(3))-(b(3)-a(3))*(c(2)-a(2)), &
       (b(3)-a(3))*(c(1)-a(1))-(b(1)-a(1))*(c(3)-a(3)), &
       (b(1)-a(1))*(c(2)-a(2))-(b(2)-a(2))*(c(1)-a(1))]
    if (dot_product(n,a+b+c)<0.0_dp) then
        tmp=mesh%face(2:3,f)
        mesh%face(2,f)=int(tmp(2))
        mesh%face(3,f)=int(tmp(1))
        mesh%face_sign(:,f)=-mesh%face_sign(:,f)
        a=mesh%xyz(:,mesh%face(1,f))
        b=mesh%xyz(:,mesh%face(2,f))
        c=mesh%xyz(:,mesh%face(3,f))
        n=[(b(2)-a(2))*(c(3)-a(3))-(b(3)-a(3))*(c(2)-a(2)), &
           (b(3)-a(3))*(c(1)-a(1))-(b(1)-a(1))*(c(3)-a(3)), &
           (b(1)-a(1))*(c(2)-a(2))-(b(2)-a(2))*(c(1)-a(1))]
    end if
    normv=sqrt(sum(n*n))
    mesh%face_area(f)=0.5_dp*normv
    mesh%face_normal(:,f)=n/normv
    mesh%circumcenter(:,f)=(a+b+c)/3.0_dp
    mesh%circumcenter(:,f)=mesh%circumcenter(:,f)/ &
        sqrt(sum(mesh%circumcenter(:,f)**2))
    mesh%bary_vertex_area(mesh%face(:,f))= &
        mesh%bary_vertex_area(mesh%face(:,f))+mesh%face_area(f)/3.0_dp
end do
end subroutine

subroutine build_dec(mesh,dec)
type(sphere_mesh), intent(in) :: mesh
type(dec_complex), intent(out) :: dec
integer :: e,f,j,v0,v1
real(dp) :: t(3),m(3),dual
dec%nv=mesh%nv
dec%ne=mesh%ne
dec%nf=mesh%nf
allocate(dec%d0(mesh%ne,mesh%nv),dec%d1(mesh%nf,mesh%ne))
allocate(dec%star0(mesh%nv),dec%star1_bary(mesh%ne))
allocate(dec%star1_circ(mesh%ne),dec%star2(mesh%nf))
dec%d0=0
dec%d1=0
do e=1,mesh%ne
    v0=mesh%edge(1,e)
    v1=mesh%edge(2,e)
    dec%d0(e,v0)=-1
    dec%d0(e,v1)=1
end do
do f=1,mesh%nf
    do j=1,3
        dec%d1(f,mesh%face_edge(j,f))=mesh%face_sign(j,f)
    end do
end do
dec%star0=mesh%bary_vertex_area
do e=1,mesh%ne
    dec%star1_bary(e)=0.0_dp
    do j=1,2
        f=mesh%edge_faces(j,e)
        dec%star1_bary(e)=dec%star1_bary(e)+ &
            mesh%face_area(f)/(3.0_dp*mesh%edge_length(e)**2)
    end do
    f=mesh%edge_faces(1,e)
    dual=0.0_dp
    if (mesh%edge_faces(2,e)>0) then
        t=(mesh%xyz(:,mesh%edge(2,e))-mesh%xyz(:,mesh%edge(1,e)))/ &
            mesh%edge_length(e)
        m=cross(mesh%face_normal(:,f),t)
        dual=abs(dot_product(mesh%circumcenter(:,mesh%edge_faces(2,e))- &
            mesh%circumcenter(:,f),m))
    end if
    dec%star1_circ(e)=dual/mesh%edge_length(e)
end do
dec%star2=1.0_dp/mesh%face_area
end subroutine

function cross(a,b) result(c)
real(dp), intent(in) :: a(3),b(3)
real(dp) :: c(3)
c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3), &
   a(1)*b(2)-a(2)*b(1)]
end function

subroutine initialize_state(mesh,dec,u,A,m)
type(sphere_mesh), intent(in) :: mesh
type(dec_complex), intent(in) :: dec
real(dp), intent(out) :: u(:),A(:),m(:)
integer :: e,f
real(dp) :: x(3),y(3),z(3),rho
u=0.0_dp
A=0.0_dp
do e=1,mesh%ne
    x=mesh%xyz(:,mesh%edge(1,e))
    y=mesh%xyz(:,mesh%edge(2,e))
    z=cross(x,y)
    A(e)=0.02_dp*z(3)
    u(e)=0.03_dp*(x(1)*y(2)-x(2)*y(1))
end do
do f=1,mesh%nf
    rho=1.0_dp+0.05_dp*mesh%face_normal(3,f)
    m(f)=rho*mesh%face_area(f)
end do
end subroutine

subroutine apply_d1(dec,x,y)
type(dec_complex), intent(in) :: dec
real(dp), intent(in) :: x(:)
real(dp), intent(out) :: y(:)
y=matmul(real(dec%d1,dp),x)
end subroutine

subroutine apply_d0(dec,x,y)
type(dec_complex), intent(in) :: dec
real(dp), intent(in) :: x(:)
real(dp), intent(out) :: y(:)
y=matmul(real(dec%d0,dp),x)
end subroutine

subroutine contract_flux(mesh,dec,u,b,flux)
type(sphere_mesh), intent(in) :: mesh
type(dec_complex), intent(in) :: dec
real(dp), intent(in) :: u(:),b(:)
real(dp), intent(out) :: flux(:)
integer :: e,f,g
real(dp) :: bp,bm
do e=1,mesh%ne
    f=mesh%edge_faces(1,e)
    g=mesh%edge_faces(2,e)
    bp=b(f)*dec%star2(f)
    bm=b(g)*dec%star2(g)
    flux(e)=0.5_dp*u(e)*(bp+bm)*mesh%edge_length(e)
end do
end subroutine

subroutine momentum_operator(mesh,dec,u,cs,mu0,m,b,op)
type(sphere_mesh), intent(in) :: mesh
type(dec_complex), intent(in) :: dec
real(dp), intent(in) :: u(:),cs,mu0,m(:),b(:)
real(dp), intent(out) :: op(:)
integer :: f,e,j,v0,v1
real(dp) :: omega_f(20),omega_v(12),kin_v(12),qf(20),qv(12)
real(dp) :: vf(3,20),t(3),me(3),uperp,den
real(dp) :: rho,bden,area
omega_f=0.0_dp
vf=0.0_dp
do f=1,mesh%nf
    do j=1,3
        e=mesh%face_edge(j,f)
        omega_f(f)=omega_f(f)+mesh%face_sign(j,f)*u(e)
        t=(mesh%xyz(:,mesh%edge(2,e))-mesh%xyz(:,mesh%edge(1,e)))/ &
            mesh%edge_length(e)
        vf(:,f)=vf(:,f)+mesh%face_sign(j,f)*mesh%edge_length(e)*u(e)*t
    end do
    vf(:,f)=vf(:,f)/mesh%face_area(f)
end do
omega_v=0.0_dp
kin_v=0.0_dp
qv=0.0_dp
do f=1,mesh%nf
    rho=m(f)/mesh%face_area(f)
    bden=b(f)/mesh%face_area(f)
    qf(f)=cs**2*log(max(rho,1.0e-12_dp))+bden*bden/(2.0_dp*mu0)
    do j=1,3
        v0=mesh%face(j,f)
        omega_v(v0)=omega_v(v0)+omega_f(f)/3.0_dp
        qv(v0)=qv(v0)+qf(f)/3.0_dp
    end do
end do
do v0=1,mesh%nv
    omega_v(v0)=omega_v(v0)/mesh%bary_vertex_area(v0)
    qv(v0)=qv(v0)/mesh%bary_vertex_area(v0)
end do
do e=1,mesh%ne
    v0=mesh%edge(1,e)
    v1=mesh%edge(2,e)
    f=mesh%edge_faces(1,e)
    area=mesh%bary_vertex_area(v0)
    kin_v(v0)=kin_v(v0)+0.5_dp*dec%star1_bary(e)*u(e)**2/area
    area=mesh%bary_vertex_area(v1)
    kin_v(v1)=kin_v(v1)+0.5_dp*dec%star1_bary(e)*u(e)**2/area
end do
op=0.0_dp
do e=1,mesh%ne
    v0=mesh%edge(1,e)
    v1=mesh%edge(2,e)
    f=mesh%edge_faces(1,e)
    den=mesh%edge_faces(2,e)
    me=0.5_dp*(cross(mesh%face_normal(:,f), &
        (mesh%xyz(:,v1)-mesh%xyz(:,v0))/mesh%edge_length(e))+ &
        cross(mesh%face_normal(:,int(den)), &
        (mesh%xyz(:,v1)-mesh%xyz(:,v0))/mesh%edge_length(e)))
    me=me/sqrt(max(sum(me*me),1.0e-30_dp))
    uperp=0.5_dp*dot_product(vf(:,f)+vf(:,int(den)),me)
    op(e)=dec%star1_bary(e)*0.5_dp*(omega_v(v0)+omega_v(v1))*uperp
    op(e)=op(e)-(kin_v(v1)-kin_v(v0))
    op(e)=op(e)-(qv(v1)-qv(v0))
end do
end subroutine

subroutine advance_step(mesh,dec,u,A,m,dt,cs,mu0,residual,iterations)
type(sphere_mesh), intent(in) :: mesh
type(dec_complex), intent(in) :: dec
real(dp), intent(inout) :: u(:),A(:),m(:)
real(dp), intent(in) :: dt,cs,mu0
real(dp), intent(out) :: residual
integer, intent(out) :: iterations
real(dp), allocatable :: un(:),um(:),an(:),am(:),mn(:),mm(:)
real(dp), allocatable :: b(:),bm(:),flux(:),rhs(:),tmp(:),du(:)
integer :: it,e
real(dp) :: err,relax
allocate(un(size(u)),um(size(u)),an(size(A)),am(size(A)))
allocate(mn(size(m)),mm(size(m)),b(size(m)),bm(size(m)))
allocate(flux(size(u)),rhs(size(u)),tmp(size(m)),du(size(u)))
call apply_d1(dec,A,b)
un=u
an=A
mn=m
relax=0.5_dp
residual=huge(1.0_dp)
do it=1,80
    um=0.5_dp*(u+un)
    am=0.5_dp*(A+an)
    mm=0.5_dp*(m+mn)
    call apply_d1(dec,am,bm)
    call contract_flux(mesh,dec,um,bm,flux)
    call apply_d1(dec,flux,tmp)
    mn=m-dt*tmp
    call apply_d1(dec,am,bm)
    an=A-dt*flux
    call momentum_operator(mesh,dec,um,cs,mu0,mm,bm,rhs)
    do e=1,size(u)
        un(e)=u(e)-dt*rhs(e)/max(dec%star1_bary(e),1.0e-12_dp)
    end do
    du=un-u
    err=max(maxval(abs(du)),maxval(abs(mn-m)),maxval(abs(an-A)))
    residual=err
    if (err<1.0e-10_dp) exit
    un=relax*un+(1.0_dp-relax)*u
    an=relax*an+(1.0_dp-relax)*A
    mn=relax*mn+(1.0_dp-relax)*m
end do
iterations=it
if (residual>1.0e-10_dp) then
    un=u
    an=A
    mn=m
    um=u
    call apply_d1(dec,A,b)
    call contract_flux(mesh,dec,u,b,flux)
    call apply_d1(dec,flux,tmp)
    mn=m-dt*tmp
    an=A-dt*flux
    call momentum_operator(mesh,dec,u,cs,mu0,m,b,rhs)
    un=u-dt*rhs/maxval(dec%star1_bary)
    residual=max(maxval(abs(un-u)),maxval(abs(mn-m)),maxval(abs(an-A)))
end if
u=un
A=an
m=mn
end subroutine

subroutine diagnostics(mesh,dec,u,A,m,mu0,ekin,emag,eint,mass,flux,exacterr)
type(sphere_mesh), intent(in) :: mesh
type(dec_complex), intent(in) :: dec
real(dp), intent(in) :: u(:),A(:),m(:),mu0
real(dp), intent(out) :: ekin,emag,eint,mass,flux,exacterr
real(dp) :: b(mesh%nf),rho
integer :: f
call apply_d1(dec,A,b)
ekin=0.5_dp*sum(dec%star1_bary*u*u)
emag=0.5_dp/mu0*sum(b*b/mesh%face_area)
eint=0.0_dp
mass=sum(m)
flux=sum(b)
do f=1,mesh%nf
    rho=m(f)/mesh%face_area(f)
    eint=eint+m(f)*log(max(rho,1.0e-12_dp))
end do
exacterr=maxval(abs(b-matmul(real(dec%d1,dp),A)))
end subroutine

subroutine destroy_mesh(mesh)
type(sphere_mesh), intent(inout) :: mesh
if (allocated(mesh%xyz)) deallocate(mesh%xyz)
if (allocated(mesh%edge)) deallocate(mesh%edge)
if (allocated(mesh%face)) deallocate(mesh%face)
if (allocated(mesh%face_edge)) deallocate(mesh%face_edge)
if (allocated(mesh%face_sign)) deallocate(mesh%face_sign)
if (allocated(mesh%edge_faces)) deallocate(mesh%edge_faces)
if (allocated(mesh%edge_length)) deallocate(mesh%edge_length)
if (allocated(mesh%face_area)) deallocate(mesh%face_area)
if (allocated(mesh%face_normal)) deallocate(mesh%face_normal)
if (allocated(mesh%circumcenter)) deallocate(mesh%circumcenter)
if (allocated(mesh%bary_vertex_area)) deallocate(mesh%bary_vertex_area)
end subroutine

end module
