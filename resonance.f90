program gridbasis
implicit none
integer, parameter :: nb=3000,nx=nb,ne=nb+nb**2,lwork=3*nb-1
real *8 pi,xo,x(nx),xst,xend,dx,rwork(2*nb),alpha
complex*16 ho(nb,nb),t(nb,nb),vv(nx),zi,eigval(nb),work(lwork),vl(nb,nb),to(nb,nb),&
vr(nb,nb),vi(nx),vo(nb,nb),h(nb,nb),vii(nb,nb),norm(nb,nb),y(ne),emat(nb),psi(nb,nb) ,cmat(nb,nb),psi1(nb,nb)
real*8 e(nb),v(nb,nb),c(nb,nb),lambda,hmat(nb,nb),dl,lambda1,theta
complex*16 vl1(nb,nb),vr1(nb,nb),eigval1(nb),vmat(nb,nb),v2i(nb,nb),xx(nx),v5(nb,nb)
integer m,n,k,i,j,info,ii,itheta,path,iii,ilamda
real *8:: lamda,xl,xr,mheu
complex*16 :: cap(nb,nb)
pi=dacos(-1.d0)
xst = 0.d0
xend = +30.d0
dx=(xend-xst)/nx
zi =(0.d0,1.d0)
alpha =0.05d0
lamda=1.5d0
path=8
xr=+22.d0
xl=-1000.d0
mheu=1.d0


do ilamda =1,20
   lamda=0.1*ilamda
   ! xr=ilamda
   ! xl=-1000.d0
   do m=1,nb
      do n=1,nb
         if ( m .eq. n ) then
            t(m,n)=(pi**2/(3.d0*2*dx**2))
         end if
         if ( m .ne. n) then
            t(m,n)=((-1)**(m-n)/(2.d0*dx**2))*2.d0/(m-n)**2
         end if
         ! write(*,*) m,n,v(m,n)
      end do
   end do

   to=t
   do i=1,nb
      do j=1,nb
         write(88,*) t(i,j),i,j
      enddo
   enddo
   vv=0.0
   do m=1,nx
      x(m)= xst+(m-1)*dx
!!! potential in use - piecewise-pot
!! Potential barrier on a half plane
      if( x(m) .ge. 3.d0 .and. x(m) .le. 4.d0) then
         vv(m)=3.d0
      endif
!! Parabolic Well Potential
!! if (abs(x(m)) .le. 4.5d0) vv(m)=x(m)**2/2.d0 !! Guassian potential

!! Double barrier potential
!! case I
! if( x(m) .ge. -4.5d0 .and. x(m) .le. -3.0d0) vv(m)=8.d0
! if( x(m) .ge. 3.d0 .and. x(m) .le. 4.5d0) vv(m)=8.d0
!! Case II
! if( x(m) .ge. -3.5d0 .and. x(m) .le. -3.0d0) vv(m)=8.d0
! if( x(m) .ge. 3.d0 .and. x(m) .le. 3.5d0) vv(m)=8.d0

!!!! addition of kinetic and potential part
      to(m,m)=to(m,m)+vv(m)
      write(10,*)x(m),real(vv(m)),aimag(vv(m))
   end do
   
   do itheta=0,100!200
      theta=0.002d0*itheta
      ! write(6,*)xl,xr,x,lamda,theta,mheu,nx,dx,path; stop
      call cal_SES_CAP_for_split_operator(cap,xl,xr,x,lamda,theta,mheu,nx,dx,path)
      !cap=0.d0
      ho=to+cap
      call ZGEEV( 'N', 'V', nb,ho, nb, eigval1, VL1, nb, VR1, nb, WORK, LWORK, RWORK, INFO)

      call sorting2(eigval1,vr1,nb)

      do i = 1,nb
         write(100+ilamda,*) real(eigval1(i)),aimag(eigval1(i));call flush(100+ilamda)
      end do

      norm=matmul(transpose(vr1),vr1)

      do i=1,nb
         vr1(:,i)=vr1(:,i)/sqrt(norm(i,i))  
      enddo

      if (itheta .eq. 200 ) then
         do ii=1,nx
            write(nb*100,*) x(ii),(real(vr1(ii,83)))+abs(eigval1(83)),&
                 (-real(vr1(ii,90)))+abs(eigval1(90)),(real(vr1(ii,98)))+abs(eigval1(98))!,(real(vr1(ii,32)))!,abs(real(vr1(ii,13)))
         enddo
         close(nb*100)
      end if
   enddo
enddo
stop
end program gridbasis
 

subroutine sorting2(a,b,n)
!integer, parameter :: dp = selected_real_kind( p=15, r=307 )
  integer n,i,j
  complex*16 a(n),temp,temp2(n),b(n,n)
  do i=1,n
   do  j=1,n-i
     if (real(a(j)) .gt. real(a(j+1)))then
        temp=a(j)   ;  temp2(:)=b(:,j)
        a(j)=a(j+1) ;  b(:,j)=b(:,j+1)
        a(j+1)=temp ;  b(:,j+1)=temp2(:)
     end if
  end do
  end do
! write(6,*)'values', a
end subroutine sorting2






subroutine cal_SES_CAP_for_split_operator(cap,xl,xr,x,lamda,theta,mheu,nx,dx,path)
  implicit none
  integer nx,path,l,j,ii,i
  complex*16 v3mat(nx),Kar(nx),cap(nx,nx)
  integer,save:: ic=1
  real*8 ke(nx,nx)
  real*8 dx,mheu,theta,x(nx),xl,xr,lamda,t1,t2
  
  call ke_dvr_ses(ke,nx,dx,mheu)

  kar=0.d0;v3mat=0.d0

  call SES_grid_basis_for_split_operator(kar,v3mat,nx,x,xl,xr,lamda,theta,mheu,path)
!  do i=10,nx
!     if (sum(abs(v3mat(i-10:i))) .lt. 1.d-10)then
!        ii=i; goto 10
!     endif
!  enddo





!10 continue
!  do i=1,nx
!     write(15,*)x(i),real(kar(i)),aimag(kar(i)),real(v3mat(i)),aimag(v3mat(i))
!  enddo
  
!  !ii=nx
!  write(6,*)ii,nx
  !call cpu_time(t1)
  ! write(6,*) 'time taken by SES grid',t2-t1
  cap=0.d0
  do l=1,nx
     do j=1,nx
        cap(l,j)=kar(l)*ke(l,j)*kar(j)-ke(l,j)
     enddo
     cap(l,l)=cap(l,l)+v3mat(l)
  enddo
  !cap=cap-ke
  !  cap=matmul(Kar,matmul(ke, Kar))    +v3mat-ke
  !call cpu_time(t2)
  !write(6,*) 'time taken by cap',t2-t1
  !stop
end subroutine cal_SES_CAP_for_split_operator



subroutine SES_grid_basis_for_split_operator(kar,v3mat,nx,x,xl,xr,lamda,theta,mheu,path)
  implicit none
  
  integer nx,path
  real*8 x(nx),xst,xl,xr,lamda,theta,mheu,l,pi,dx
  complex*16 g(nx),gp(nx),g2p(nx),v3mat(nx),Kar(nx)
  complex*16 zi,f(nx),fp(nx),f2p(nx),v3(nx),v0(nx),v1(nx),v2(nx)
  
  integer i,m,n
  zi=(0.d0,1.d0)
  pi=4.d0*atan(1.d0)
  dx=x(2)-x(1)
  
  
  call gdiff(g,gp,g2p,x,xl,xr,nx,lamda,theta,path)
  
  
  do i=1,nx
     write(12,*)x(i),real(g(i)),aimag(g(i)),real(gp(i)),real(g2p(i)),aimag(g2p(i))
  enddo
  
  close(12)
  f=1.d0+(exp(zi*theta)-1.d0)*g
  fp=(exp(zi*theta)-1.d0)*gp
  f2p=(exp(zi*theta)-1.d0)*g2p


  do i=1,nx
     write(19,*) x(i),real(f(i)),aimag(f(i))
  enddo
  
  v3=-1.d0/8.d0/mheu*(2.d0*f2p*f -3.d0*fp**2)/f**4
    
  kar=0 ;v3mat=0.d0 
  do m=1,nx
     Kar(m)  = 1.d0/f(m)
     v3mat(m)=v3(m) 
  enddo
  
end subroutine SES_GRID_BASIS_FOR_SPLIT_OPERATOR





!---------------------------------------------------------------
subroutine gdiff(g,gp,g2p,x,xl,xr,nx,lamda,theta,path)
  implicit none
  integer nx,i,path
  real*8 x(nx),lamda,theta,xl,xr
  complex*16  g(nx),gp(nx),g2p(nx),g1,g2,gfn
  
  real*8,parameter:: dx=0.001d0 !this dx is not same as main program dx
  
  do i=1,nx
     g(i)=gfn(x(i),xl,xr,lamda,theta,path)
     write(18,*)x(i),real(g(i)),aimag(g(i))
   enddo
  
  do i=1,nx
     g1=gfn(x(i)-dx,xl,xr,lamda,theta,path)
     g2=gfn(x(i)+dx,xl,xr,lamda,theta,path)
     gp(i)=(g2-g1)/2.d0/dx
     g2p(i)=(g2+g1-2.d0*g(i))/dx**2
     !write(18,*)x(i)-dx,g1
     !write(17,*)x(i)+dx,g2
  enddo
  
end subroutine gdiff



 
function gfn(x,xl,xr,lamda,theta,path)
  
  implicit none
  integer path
  real*8 x,lamda,theta,xl,xr
  complex*16 gfn,zi,ap,bp,dp,ep,thetax,aap
  real*8,parameter::lamdap=1.d0
  real*8,parameter:: x0p=0.d0
  
  zi=(0.d0,1.d0)
  
  if(path .eq. 0) then
     gfn=1.d0+0.5d0*(tanh(lamda*(x-xr))-tanh(lamda*(x-xl)))!moiseyev
      
  elseif(path .eq. 1) then
     gfn=(1.d0+exp((xr-x)*lamda))**(-1)+(1.d0+exp((-xl+x)*lamda))**(-1)!woodsaxon
     
  elseif(path .eq. 2) then
     
     gfn=(1.d0+exp((xr-x)*lamda))**(-1)+(1.d0+exp((-xl+x)*lamda))**(-1)+x*lamda*&
          (exp((xr-x)*lamda)*(1.d0+exp((xr-x)*lamda))**(-2)-exp((-xl+x)*lamda)*&
          (1.d0+exp((-xl+x)*lamda))**(-2))
     
  elseif(path .eq. 3) then
     if (x.gt. xl .and. x.le.xr) then
        gfn=0.d0
     else
        
        gfn=1.d0-exp(-lamda*(x-xr)**2)&
             +2.d0*lamda*(x-xr)**2*exp(-lamda*(x-xr)**2)&
             -exp(-lamda*(x-xl)**2)&
             +2.d0*lamda*(x-xl)**2*exp(-lamda*(x-xl)**2)  !elander
     endif
  elseif(path .eq. 4) then
     gfn=1.d0 !complex scalling
     
     !  elseif(path .eq. 5) then
     !     gfn=(1.d0+exp((xr-x)*lamda))**(-1.d0)+(1.d0+exp((-xl+x)*lamda))**(-1.d0) &
     !-p*zi*((1.d0+exp((x0p-x)*lamdap))**(-1.d0)+(1.d0+exp((x0p+x)*lamdap))**(-1.d0))
     !     gfn=(exp(zi*gfn*theta)-1.d0)/(exp(zi*theta)-1)&
     !          +zi*theta*x*exp(zi*gfn*theta)/(exp(zi*theta)-1)&
     !          *((exp((xr-x)*lamda)*(1.d0+exp((xr-x)*lamda))**(-2.d0)&
     !          -exp((-xl+x)*lamda)*(1.d0+exp((-xl+x)*lamda))**(-2.d0))*lamda &
     !-p*zi*((exp((x0p-x)*lamdap)*(1.d0+exp((x0p-x)*lamdap))**(-2.d0)&
     !          -exp((x0p+x)*lamdap)*(1.d0+exp((x0p+x)*lamdap))**(-2.d0))*lamdap))
     
     
     
  elseif(path .eq. 6) then
     ap=lamda*exp((xr-x)*lamda)*(1.d0+exp((xr-x)*lamda))**(-2)
     bp=-lamda*exp((-xl+x)*lamda)*(1.d0+exp((-xl+x)*lamda))**(-2)
     
     dp=lamda*exp((xr-x)*lamda)*(1.d0+exp((xr-x)*lamda))**(-2)
     
     dp=dp+2.d0*x*lamda**(2)*exp((xr-x)*lamda*2.d0)*(1.d0+exp((xr-x)*lamda))**(-3)
     
     dp=dp-x*lamda**(2)*exp((xr-x)*lamda)*(1.d0+exp((xr-x)*lamda))**(-2)
     
     ep=lamda*exp((-xl+x)*lamda)*(1.d0+exp((-xl+x)*lamda))**(-2)
     
     ep=ep-2.d0*x*lamda**(2)*exp((-xl+x)*lamda*2.d0)*(1.d0+exp((-xl+x)*lamda))**(-3)
     ep=ep+x*lamda**(2)*exp((-xl+x)*lamda)*(1.d0+exp((-xl+x)*lamda))**(-2)
     
     aap=ap+bp+dp-ep
     
     
     thetax=((1.d0+exp((xr-x)*lamda))**(-1.d0)+(1.d0+exp((-xl+x)*lamda))**(-1)+x*lamda*&
          (exp((xr-x)*lamda)*(1.d0+exp((xr-x)*lamda))**(-2)-exp((-xl+x)*lamda)*&
          (1.d0+exp((-xl+x)*lamda))**(-2)))*theta   !thetax is taylor form gfn multiply by theta
     
     gfn=(exp(zi*thetax)-1.d0)/(exp(zi*theta)-1.d0)&
          +(zi*theta*x*exp(zi*thetax)/(exp(zi*theta)-1.d0))*aap
     
     
     !write(73,*)x,real(thetax),real(bp),real(dp),real(ep),aimag(ap),aimag(bp),aimag(dp),aimag(ep)   
     
  elseif(path .eq. 7) then!elander space
     if (x.gt.xl.and.x.le.xr) then
        gfn=0.d0
     else
        ap=-(x-xr)*exp(-lamda*(x-xr)**2)
        bp=-2.d0*(x-xr)*exp(-lamda*(x-xr)**2)*(1.d0-lamda*(x-xr)**2)
        dp=-(x-xl)*exp(-lamda*(x-xl)**2)
        ep=-2.d0*(x-xl)*exp(-lamda*(x-xl)**2)*(1.d0-lamda*(x-xl)**2)
        aap=ap+bp+dp+ep
        thetax=(1.d0-exp(-lamda*(x-xr)**2)&
             +2.d0*lamda*(x-xr)**2*exp(-lamda*(x-xr)**2)&
             -exp(-lamda*(x-xl)**2)&
             +2.d0*lamda*(x-xl)**2*exp(-lamda*(x-xl)**2))*theta
        
        gfn=(exp(zi*thetax)-1.d0)/(exp(zi*theta)-1.d0)&
             -2.d0*lamda*(zi*theta*x*exp(zi*thetax)/(exp(zi*theta)-1.d0))*aap
        
     endif
     
  elseif(path .eq. 8) then
     if(theta.lt. 1e-7)then
        gfn=0
     else
        thetax=(1.d0+0.5d0*(tanh(lamda*(x-xr))-tanh(lamda*(x-xl))))*theta
        
     
        aap=0.5d0*lamda*(1.d0/(cosh(lamda*(x-xr)))**2-1.d0/(cosh(lamda*(x-xl)))**2)
      
     gfn=(exp(zi*thetax)-1.d0)/(exp(zi*theta)-1.d0)&
          +(zi*theta*x*exp(zi*thetax)/(exp(zi*theta)-1.d0))*aap
       write(17,*) x,real(gfn),aimag(gfn)
     endif
  elseif(path .eq. 9) then
     gfn=exp(-x**2)
  endif
  

end function gfn

subroutine ke_dvr_ses(ke,nx,dx,mheu)
  implicit none
  integer nx,i,j
  real*8 ke(nx,nx),mheu,dx,pi
  pi=datan(1.d0)*4.d0
  do i=1,nx
     do j=1,nx
        if(i.eq.j)then
           KE(i,j)=pi**2.d0/3.d0
        else
           KE(i,j)=2.d0/(i-j)**2*(-1.d0)**(i-j)
        endif
     enddo
  enddo
  ke=ke/2.d0/mheu/dx**2
end subroutine ke_dvr_ses


subroutine propagation_taylor_series(h_total,udt,nb,dt)!udt=short time evolution operator
  implicit none
  integer i,nb,j,ii,jj,n_taylor
  complex*16 h_total(nb,nb),udt(nb,nb),storage(nb,nb),zi
  real*8 aa,dt
  zi=(0.d0,1.d0)
  udt=0.d0
  n_taylor=5000

  storage=-zi*h_total*dt
  udt=storage
  do i=1,nb
   udt(i,i)=udt(i,i)+1.d0
enddo
aa=dt

do i=2,n_taylor
   storage=-zi*dt*matmul(h_total,storage)/i
   udt=udt+storage
   if(sum(abs(storage(:,:))) .lt. 1.0d-16) goto 10
   write(6,*) 'i-taylor,term',i,sum(abs(storage(:,:)));call flush(6)
enddo
10 write(6,*)'exp(-i*h*dt) done',i
end subroutine propagation_taylor_series



