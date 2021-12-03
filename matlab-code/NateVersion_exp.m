% Vorticity advection equation with a semi-Lagrangian advection algorithm.
% This version includes background shear and dust as a second fluid.
% Dust can drift with respect to gas with terminal velocity.
% This version has full two-way coupling between gas and dust.
profile on 
Nx            = 256;
Ny            = 256;
Lx            = 4.0;
Ly            = 4.0;
dx            = Lx/Nx;
dy            = Ly/Ny
x             = (0:1:Nx-1)*dx-0.5*Lx;
y             = (0:1:Ny-1)*dy-0.5*Ly;
[x,y]         = ndgrid(x,y);
kx            = (0:1:Nx-1);
ky            = (0:1:Ny-1);
kx(Nx/2+2:Nx) = kx(Nx/2+2:Nx)-Nx;
ky(Ny/2+2:Ny) = ky(Ny/2+2:Ny)-Ny;
kx            = kx*(2*pi/Lx);
ky            = ky*(2*pi/Ly);
[kx,ky]       = ndgrid(kx,ky);
k2            = kx.^2+ky.^2;
k2(1,1)       = 1.0E64;

% Parameters
rho0          = 1.0;
omega         = 1.0;
shear         = -1.5*0;
dPdR          = -0.10;
tau           = 0.005;
dt            = 1.0/8.0;
nt            = 4096;
bufx          = Nx/16;
bufy          = Ny/16;
num_sl_disp_iter  = 4;
num_pressure_iter = 4;
hypviscpow    = 8;
hypvisc       = exp(-8.0*((k2/(pi/dx)^2).^(hypviscpow/2)));
hypvisc(1,1)  = 1.0;
interp_method1 = 'cubic';
interp_method2 = 'spline';

% Initialize arrays
wzq           = zeros(Nx,Ny,nt+1);
vx            = zeros(Nx,Ny);
vxb           = zeros(Nx,Ny);
vy            = zeros(Nx,Ny);
psi           = zeros(Nx,Ny);
delx          = zeros(Nx,Ny);
dely          = zeros(Nx,Ny);
xi            = zeros(Nx,Ny);
yi            = zeros(Nx,Ny);
vx_buf        = zeros(Nx+2*bufx,Ny+2*bufy);
vy_buf        = zeros(Nx+2*bufx,Ny+2*bufy);
wz_buf        = zeros(Nx+2*bufx,Ny+2*bufy);
x_buf         = zeros(Nx+2*bufx,Ny+2*bufy);
y_buf         = zeros(Nx+2*bufx,Ny+2*bufy);

rho           = zeros(Nx,Ny,nt+1);
vxw_x         = zeros(Nx,Ny);
vxw_y         = zeros(Nx,Ny);
nlxf          = zeros(Nx,Ny);
nlyf          = zeros(Nx,Ny);
hf            = zeros(Nx,Ny); 
dPdx          = zeros(Nx,Ny);
dPdy          = zeros(Nx,Ny);            
qx            = zeros(Nx,Ny);
qy            = zeros(Nx,Ny);
divq          = zeros(Nx,Ny);
crlq          = zeros(Nx,Ny);
rho_buf       = zeros(Nx+2*bufx,Ny+2*bufy);
divq_buf      = zeros(Nx+2*bufx,Ny+2*bufy);
crlq_buf      = zeros(Nx+2*bufx,Ny+2*bufy);

x_buf = add_buffer(x,bufx,bufy);
x_buf(1:bufx,:) = x_buf(1:bufx,:)-Lx;
x_buf(end-bufx+1:end,:) = x_buf(end-bufx+1:end,:)+Lx;
y_buf = add_buffer(y,bufx,bufy);
y_buf(:,1:bufy) = y_buf(:,1:bufy)-Ly;
y_buf(:,end-bufy+1:end) = y_buf(:,end-bufy+1:end)+Ly;

% initialize background shear (set to zero above via variable "shear")
vxb  = -shear*y.*(0.5*(tanh(8.0*(y+0.4*Ly))-tanh(8.0*(y-0.4*Ly))));

% initialize vorticity
wzq(:,:,1) = 0*x;

% initialize dust density
rho(:,:,1) = 0.1*(1.0+0.01*noise2d(sqrt(k2),pi/dx/64,pi/dx/2,1,1))

% find velocity from vorticity via streamfunction
psi = fft2d(wzq(:,:,1),-2)./k2;
vx  = fft2d(+1i*ky.*psi,+2);
vy  = fft2d(-1i*kx.*psi,+2);
vxw_x = -0.5*fft2d(vx.^2+vy.^2,-2);
vxw_y = 1i*ky.*vxw_x;
vxw_x = 1i*kx.*vxw_x;
vxw_x = vxw_x+fft2d(vy.*(wzq(:,:,1)+2.0*(omega+shear)),-2);
vxw_y = vxw_y-fft2d(vx.*(wzq(:,:,1)+2.0* omega       ),-2);

% compute gas pressure and dust drift velocity
fac1 = (rho(:,:,1)/rho0)./((1.0+rho(:,:,1)/rho0).^2+(2.0*omega*tau)^2);
for jj=1:num_pressure_iter
    nlxf = vxw_x+qx;
    nlyf = vxw_y+qy;
    hf   = -1i*(kx.*nlxf+ky.*nlyf)./k2;
    dPdx = fft2d(1i*kx.*hf,+2);
    dPdy = fft2d(1i*ky.*hf,+2)+dPdR;
    qx   = fft2d(fac1.*((1.0+rho(:,:,1)/rho0).*dPdx+2.0*omega*tau*dPdy),-2);
    qy   = fft2d(fac1.*((1.0+rho(:,:,1)/rho0).*dPdy-2.0*omega*tau*dPdx),-2);
end
divq = dt*fft2d(1i*(kx.*qx+ky.*qy),+2)*(rho0*tau)
crlq = dt*fft2d(1i*(kx.*qy-ky.*qx),+2);

% iterate to find displacements
vx_buf = add_buffer(vx+vxb,bufx,bufy);
vy_buf = add_buffer(vy    ,bufx,bufy);
delx          = dt*vx;
dely          = dt*vy;
xi            = x-delx;
yi            = y-dely;
xi(xi> Lx/2)  = xi(xi> Lx/2)-Lx;
xi(xi<-Lx/2)  = xi(xi<-Lx/2)+Lx;
yi(yi> Ly/2)  = yi(yi> Ly/2)-Ly;
yi(yi<-Ly/2)  = yi(yi<-Ly/2)+Ly;
for jj=1:num_sl_disp_iter-1
    delx          = dt*interpn(x_buf,y_buf,vx_buf,xi,yi,interp_method1);
    dely          = dt*interpn(x_buf,y_buf,vy_buf,xi,yi,interp_method1);
    xi            = x-delx;
    yi            = y-dely;
    xi(xi> Lx/2)  = xi(xi> Lx/2)-Lx;
    xi(xi<-Lx/2)  = xi(xi<-Lx/2)+Lx;
    yi(yi> Ly/2)  = yi(yi> Ly/2)-Ly;
    yi(yi<-Ly/2)  = yi(yi<-Ly/2)+Ly;
end

% advect forward one step
wz_buf = add_buffer(wzq(:,:,1),bufx,bufy);
rho_buf= add_buffer(rho(:,:,1),bufx,bufy);
wzq(:,:,2) = fft2d(hypvisc.*fft2d(interpn(x_buf,y_buf,wz_buf,xi,yi,interp_method2)  + crlq,-2),+2);
rho(:,:,2) = fft2d(hypvisc.*fft2d(interpn(x_buf,y_buf,rho_buf,xi,yi,interp_method2) - divq,-2),+2);

for tt=2:nt
tt

% find velocity from vorticity via streamfunction
psi = fft2d(wzq(:,:,tt),-2)./k2;
vx  = fft2d(+1i*ky.*psi,+2);
vy  = fft2d(-1i*kx.*psi,+2);
vxw_x = -0.5*fft2d(vx.^2+vy.^2,-2);
vxw_y = 1i*ky.*vxw_x;
vxw_x = 1i*kx.*vxw_x;
vxw_x = vxw_x+fft2d(vy.*(wzq(:,:,tt)+2.0*(omega+shear)),-2);
vxw_y = vxw_y-fft2d(vx.*(wzq(:,:,tt)+2.0* omega       ),-2);

% compute gas pressure and dust drift velocity
fac1 = (rho(:,:,tt)/rho0)./((1.0+rho(:,:,tt)/rho0).^2+(2.0*omega*tau)^2);
for jj=1:num_pressure_iter
    nlxf = vxw_x+qx;
    nlyf = vxw_y+qy;
    hf   = -1i*(kx.*nlxf+ky.*nlyf)./k2;
    dPdx = fft2d(1i*kx.*hf,+2);
    dPdy = fft2d(1i*ky.*hf,+2)+dPdR;
    qx   = fft2d(fac1.*((1.0+rho(:,:,tt)/rho0).*dPdx+2.0*omega*tau*dPdy),-2);
    qy   = fft2d(fac1.*((1.0+rho(:,:,tt)/rho0).*dPdy-2.0*omega*tau*dPdx),-2);
end
divq = dt*fft2d(1i*(kx.*qx+ky.*qy),+2)*(rho0*tau);
crlq = dt*fft2d(1i*(kx.*qy-ky.*qx),+2);

% iterate to find displacements
vx_buf = add_buffer(vx+vxb,bufx,bufy);
vy_buf = add_buffer(vy    ,bufx,bufy);
for jj=1:num_sl_disp_iter
    delx          = dt*interpn(x_buf,y_buf,vx_buf,xi,yi,interp_method1);
    dely          = dt*interpn(x_buf,y_buf,vy_buf,xi,yi,interp_method1);
    xi            = x-delx;
    yi            = y-dely;
    xi(xi> Lx/2)  = xi(xi> Lx/2)-Lx;
    xi(xi<-Lx/2)  = xi(xi<-Lx/2)+Lx;
    yi(yi> Ly/2)  = yi(yi> Ly/2)-Ly;
    yi(yi<-Ly/2)  = yi(yi<-Ly/2)+Ly;
end
xi2           = x-2.0*delx;
yi2           = y-2.0*dely;
xi2(xi2> Lx/2)= xi2(xi2> Lx/2)-Lx;
xi2(xi2<-Lx/2)= xi2(xi2<-Lx/2)+Lx;
yi2(yi2> Ly/2)= yi2(yi2> Ly/2)-Ly;
yi2(yi2<-Ly/2)= yi2(yi2<-Ly/2)+Ly;

% advect forward one more step
wz_buf = add_buffer(wzq(:,:,tt-1),bufx,bufy);
rho_buf= add_buffer(rho(:,:,tt-1),bufx,bufy);
divq_buf = add_buffer(divq(:,:),bufx,bufy);
crlq_buf = add_buffer(crlq(:,:),bufx,bufy);
wzq(:,:,tt+1) = fft2d(hypvisc.*fft2d(interpn(x_buf,y_buf,wz_buf,xi2,yi2,interp_method2)  + 2.0*interpn(x_buf,y_buf,crlq_buf,xi,yi,interp_method2),-2),+2);
rho(:,:,tt+1) = fft2d(hypvisc.*fft2d(interpn(x_buf,y_buf,rho_buf,xi2,yi2,interp_method2) - 2.0*interpn(x_buf,y_buf,divq_buf,xi,yi,interp_method2),-2),+2);
end

%save('rho_005.mat', 'rho', '-v7.3')
%save('wzq_005.mat', 'wzq', '-v7.3')

profile viewer

function [g]=fft2d(f,dir)
% If dir<0: real physical space     --> complex frequency space
% If dir>0: complex frequency space --> real physical space
    Nx = size(f,1);
    Ny = size(f,2); 
    if dir<0
        g = fft2(real(f));
        g(Nx/2+1,:) = 0.0;
        g(:,Ny/2+1) = 0.0;
    elseif dir>0
        f(Nx/2+1,:) = 0.0;
        f(:,Ny/2+1) = 0.0;
        g = real(ifft2(f));
    end
end

function [g]=add_buffer(f,bufx,bufy)
    Nx = size(f,1);
    Ny = size(f,2); 
    g = zeros(Nx+2*bufx,Ny+2*bufy);    
    g(1+bufx:Nx+bufx,1+bufy:Ny+bufy) = f;
    g(1:bufx,1+bufy:Ny+bufy) = f(Nx-bufx+1:Nx,:);
    g(1+bufx+Nx:Nx+2*bufx,1+bufy:Ny+bufy) = f(1:bufx,:);
    g(:,1:bufy) = g(:,Ny+1:Ny+bufy);
    g(:,Ny+bufy+1:Ny+2*bufy) = g(:,1+bufy:2*bufy);
end

%%%Initialize in Fourier space a grid of random noise 
function noise = noise2d(k,kmin,kmax,kpow,rms_noise) 
[imax,jmax] = size(k);
for jj=2:jmax
    for ii=2:imax
        if ((k(ii,jj)>=kmin)&&(k(ii,jj)<=kmax))
            noise(ii,jj) = exp(2i*pi*rand)/(k(ii,jj)^kpow);
        end
    end
end
noise = fft2d(noise,+2);
rms = sqrt(mean(mean(noise.^2)));
noise = noise/rms*rms_noise;
end