%k = rand(4,4)

%noise2d(k,.1,.6,1,1)
k =  [0.7749, 0.3998, 0.9106, 0.1361;
      0.8173, 0.2599, 0.1818, 0.8693;
      0.8687, 0.8001, 0.2638, 0.5797;
      0.0844, 0.4314, 0.1455, 0.5499]


kkft = fft2d(k,-2)
kc = [8.0732 + 0.0000i,  1.0436 + 0.2438i,   0.0208 + 0.0000i,   1.0436 - 0.2438i;
  -0.2909 - 0.9171i,  -0.2497 - 0.7399i,   1.3969 - 0.6213i,  -1.2315 - 0.6533i;
   1.3942 + 0.0000i,  -0.1052 - 1.2120i,   1.7838 + 0.0000i,  -0.1052 + 1.2120i;
  -0.2909 + 0.9171i,  -1.2315 + 0.6533i,   1.3969 + 0.6213i,  -0.2497 + 0.7399i]
kc
kk = fft2d(kc,2)
kk = fft2d(kkft,2)

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