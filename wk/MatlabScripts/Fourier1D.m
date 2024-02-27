%====================================================
%
% Compute 1D Fourier Transform of Function F(1:nt,1:nx)
% along the coordinate X
%
%   Input: F(1:nt,1:nx),  !!!  "nx" has to be odd  !!!
%	   x(1:nx)
%
%   Output: TFF(1:nt,1:nx)
%	    kx(1:nx)
%
%====================================================

function [TFF,kx] = Fourier1D(F0,x0)

nx0 = length(x0);
nx  = 2*fix(nx0/2);
hnx = nx/2;
x   = x0(1:nx);
F   = F0(:,1:nx);

%--------------------------------------------
% Creation of the vector in the Fourier space
% kx(1:nx)

Lx           = x(nx) - x(1);
dx           = x(2) - x(1);
dkx          = 2*pi./(Lx+dx);
kx(1:hnx)    = -(hnx:-1:1) * dkx;
kx(hnx+1:nx) = (0:hnx-1) * dkx;

%--------------------------------------------
% Fourier transformation using Matlab fft
n1  = size(F,1);
TFF = zeros(n1,nx);

var = fft(F')'/nx;
TFF(:,1:hnx)    = var(:,hnx+1:nx);
TFF(:,hnx+1:nx) = var(:,1:hnx);
