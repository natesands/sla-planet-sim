k = rand(4,4);
noise2d(k,.1,.6,1,1)



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