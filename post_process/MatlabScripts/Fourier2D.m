%====================================================
%
% Compute 2D Fourier Transform of Function
%   z(1:ny,1:nx)
%
%   Input: F(1:ny,1:nx)
%	   y(1:ny)
%	   x(1:nx)
%
%   Output: TFF(1:ny,1:nx)
%	    ky(1:ny)
%	    kx(1:nx)
%
%====================================================

function [TFF,ky,kx] = Fourier2D(F0,y0,x0)

nx0 = length(x0);
nx  = 2*fix(nx0/2);
hnx = nx/2;
ny0 = length(y0);
ny  = 2*fix(ny0/2);
hny = ny/2;

x   = x0(1:nx);
y   = y0(1:ny);
F   = F0(1:ny,1:nx);

%--------------------------------------------
% Creation of the vector in the Fourier space
% ky(1:ny) kx(1:nx)

Lx           = x(nx) - x(1);
dx           = x(2) - x(1);
dkx          = 2*pi./(Lx+dx);
kx(1:hnx)    = -(hnx:-1:1) * dkx;
kx(hnx+1:nx) = (0:hnx-1) * dkx;

Ly           = y(ny) - y(1);
dy           = y(2) - y(1);
dky          = 2*pi./(Ly+dy);
ky(1:hny)    = -(hny:-1:1) * dky;
ky(hny+1:ny) = (0:hny-1) * dky;

%--------------------------------------------
TFF = zeros(ny,nx);

var = fft2(F)/nx/ny;
AA(:,1:hnx)     = var(:,hnx+1:nx);
AA(:,hnx+1:nx)  = var(:,1:hnx);
TFF(1:hny,:)    = AA(hny+1:ny,:);
TFF(hny+1:ny,:) = AA(1:hny,:);
