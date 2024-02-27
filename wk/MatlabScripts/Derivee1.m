%==========================================================
%
%  Calcul de la derivee premiere d'une fonction
%  en differences finies a l'ordre o(dx^4)
%
%  Input:	F  = Fonction a deriver (vecteur 1D)
%			dx = pas de la variable selon laquelle deriver
%			Optionnel:	PERIODIQUE= 1 si F periodique
%						    0 sinon (par defaut)
%
%  Output:	dFdx = derivee premiere de F
%
%==========================================================

function [dFdx] = Derivee1(F,dx,varargin)

PERIODIQUE = 0;
sze  = size(varargin);
if (max(sze) == 1)
  PERIODIQUE = varargin{1};
end
nx = length(F);

c0 = 2/3;
dFdx(3:nx-2) = c0/dx * ( F(4:nx-1)-F(2:nx-3) - (F(5:nx)-F(1:nx-4))/8 );

c1 = 4/3;
c2 = 25/12;
c3 = 5/6;
if (PERIODIQUE==0)
   dFdx(1)    = (-F(5)/4 + c1*F(4) - 3*F(3) + 4*F(2) - c2*F(1))/dx;
   dFdx(nx)   = (F(nx-4)/4 - c1*F(nx-3) + 3*F(nx-2) - 4*F(nx-1) + c2*F(nx))/dx;
   dFdx(2)    = (F(5)/12 - F(4)/2 + F(3)/c0 - c3*F(2) - F(1)/4)/dx;
   dFdx(nx-1) = ( F(nx)/4 + c3*F(nx-1) - F(nx-2)/c0 + F(nx-3)/2 - F(nx-4)/12)/dx;
else
   dFdx(1)    = c0/dx * ( F(2)-F(nx) - (F(3)-F(nx-1))/8 );
   dFdx(2)    = c0/dx * ( F(3)-F(1) - (F(4)-F(nx))/8 );
   dFdx(nx)   = c0/dx * ( F(1)-F(nx-1) - (F(2)-F(nx-2))/8 );
   dFdx(nx-1) = c0/dx * ( F(nx)-F(nx-2) - (F(1)-F(nx-3))/8 );
end
