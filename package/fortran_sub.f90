      subroutine DTE(nz,ny,nx,u0,v0,tk0,u1,v1,tk1,p,TE3d,TE2d)
      !f2py -c fortran_sub.f90 -m f90_sub
      !Input: u, v, tk are the difference of member and mean
      !       p is the pressure for certain member     
      !Example: te3d,te2d = dte(u,v,tk,p,[nz,ny,nx]) 
      !numpy.float64(variable)
      !Difference Total Energy  
      !Zhang et al. 2003, Palmer et al. 1998 eq.11, 
      !Nielsen and Schumacher 2016 eq(B3),(B4)
      !https://stackoverflow.com/questions/37592168/passing-arrays-from-python-to-fortran-with-f2py-results-error-2-shape-related
      implicit none
      integer i, j, k 
      integer, intent(in) :: nz, ny, nx
      real*4,dimension(nz,ny,nx), intent(in) :: u0, v0, tk0, &
                                             &  u1, v1, tk1, &
                                             &   p
      real*8 :: du, dv, dt, dpb
      real*8, dimension(nz,ny,nx), intent(out) :: TE3d
      real*8, dimension(ny,nx), intent(out) :: TE2d
      real*8, parameter :: kappa = 1004.9/270.0 !Cp/Tr
      
      TE3d(:,:,:) = 0.0; TE2d(:,:) = 0.0      
      do i = 1, nx  
         do j = 1, ny  
            dpb = dble(p(nz,j,i)-p(1,j,i))
            do k = 1, nz
               du = dble(u1(k,j,i)-u0(k,j,i))  
               dv = dble(v1(k,j,i)-v0(k,j,i))  
               dt = dble(tk1(k,j,i)-tk0(k,j,i))  
               TE3d(k,j,i) = 0.5*(du**2.0 + dv**2.0 + kappa*dt**2.0)
               if (k.lt.nz) then
                   TE2d(j,i) = TE2d(j,i) + &
                     & dble(p(k+1,j,i)-p(k,j,i))/dpb * TE3d(k,j,i) 
               end if
            end do 
         end do 
      end do
      end subroutine DTE
     
      subroutine TKE(nz,ny,nx,u,v,p,KE2d)
      !Input: u, v, p are the value of ensmeble mean  
      !Nielsen and Schumacher 2016 eq(B5)
      implicit none
      integer i, j, k 
      integer, intent(in) :: nz, ny, nx
      real*4, dimension(nz,ny,nx), intent(in) :: u, v, p
      real*8 :: dpb
      real*8, dimension(ny,nx), intent(out) :: KE2d
      
      KE2d(:,:) = 0.0      
      do i = 1, nx  
         do j = 1, ny  
            dpb = dble(p(nz,j,i)-p(1,j,i))
            do k = 1, nz-1
               KE2d(j,i) = KE2d(j,i) + dble(p(k+1,j,i)-p(k,j,i))/dpb &
                         & * 0.5*dble(u(k,j,i)**2.0 + v(k,j,i)**2.0)
            end do 
         end do 
      end do
      end subroutine TKE 

      subroutine DTE_wave(nz,nk,du,dv,dtk,p,TE1d)
      !input standard deviation u, v, t, avg(p)
      implicit none
      integer z, k
      integer, intent(in) :: nz,nk
      real*4,dimension(nz), intent(in) :: p
      real*4,dimension(nz,nk), intent(in) :: du, dv, dtk
      real*8 :: duu, dvv, dtt, dpb
      real*8, dimension(nk), intent(out) :: TE1d
      real*8, parameter :: kappa = 1004.9/270.0 !Cp/Tr

      TE1d(:) = 0.0
      dpb = dble(p(nz)-p(1))
      do k = 1,nk ! spectrum
         do z = 1, nz
            duu = dble(du(z,k))
            dvv = dble(dv(z,k))
            dtt = dble(dtk(z,k))
            if (z.lt.nz) then
                TE1d(k) = TE1d(k) + dble(p(z+1)-p(z))/dpb * &
                  & 0.5*(duu**2 + dvv**2 + kappa*dtt**2)
            end if
         end do
      end do
      end subroutine DTE_wave
       
      subroutine to1davg(nz,ny,nx,var3d,vary,varx,var2d)
      implicit none
      integer i, j, k
      integer, intent(in) :: nz, ny, nx
      real*8, dimension(nz,ny,nx), intent(in) :: var3d 
      real*8, intent(out) :: vary(ny), varx(nx), var2d(ny,nx)
      integer :: nj(ny), ni(nx), n2d(ny,nx)
    
      vary = 0.0; varx = 0.0; var2d = 0.0
      nj = 0; ni = 0; n2d = 0
      do i = 1, nz
         do j = 1, ny
            do k = 1, nx
               if (var3d(k,j,i) .gt. -1.d10) then
                   vary(j) = vary(j) + var3d(k,j,i)
                   nj(j) = nj(j) + 1
                   varx(i) = varx(i) + var3d(k,j,i)
                   ni(i) = ni(i) + 1
                   var2d(j,i) = var2d(j,i) + var3d(k,j,i)
                   n2d(j,i) = n2d(j,i) + 1
               end if
            end do
         end do
      end do
        
      do j = 1, ny
         !print*,"--",nj(j)
         vary(j) = vary(j)/dble(nj(j))
         do i = 1, nx
            var2d(j,i) = var2d(j,i)/dble(n2d(j,i))
         end do
      end do
      do i = 1, nx
         varx(i) = varx(i)/dble(ni(i))
      end do
   
      end subroutine to1davg    
   
      subroutine read_bin(ff,u,v,t)
      !double check the reading of binary file
      implicit none
      character(len=100), intent(in) :: ff
      integer, parameter :: nx = 150, ny = 150, nz = 150, nt = 61
      real*4, dimension(nx,ny,nz,nt), intent(out) :: u, v, t

      open(11, file=trim(ff),form="unformatted",access="direct", &
          & recl=4*nx*ny*nz*nt)
      
      read(11,rec=1) u
      read(11,rec=2) v
      read(11,rec=3) t
      close(11)
      end subroutine read_bin

      subroutine read_bin2(ff,nv,nz,ny,nx,t,var)
      !double check the reading of binary file
      implicit none
      character(len=100), intent(in) :: ff
      integer, intent(in) :: nv,nz,ny,nx,t
      real*4, dimension(nx,ny,nz,nv), intent(out) :: var

      open(11, file=trim(ff),form="unformatted",access="direct", &
          & recl=4*nx*ny*nz*nv)

      read(11,rec=t+1) var
      print*,'Fortran: (x,y,z,v)(4,3,2,1)',var(4,3,2,1)
      close(11)
      end subroutine read_bin2
 
