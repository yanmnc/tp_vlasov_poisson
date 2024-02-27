%==========================================================
%
%  Calcul de la derivee seconde d'une fonction
%  en differences finies a l'ordre o(dx^4)
%
%  Input:	F  = Fonction a deriver (vecteur 1D)
%			dx = pas de la variable selon laquelle deriver
%			Optionnel:	PERIODIQUE= 1 si F periodique
%				                    0 sinon (par defaut)
%
%  Output:	d2Fdx = derivee seconde de F
%
%==========================================================

function [d2Fdx] = Derivee2(F,dx,varargin)

PERIODIQUE = 0;
sze  = size(varargin);
if (max(sze) == 1)
  PERIODIQUE = varargin{1};
end
nx  = length(F);
dx2 = dx*dx;

d2Fdx(3:nx-2) = (-30*F(3:nx-2) + 16*(F(4:nx-1)+F(2:nx-3)) - ...
                (F(1:nx-4)+F(5:nx))) / (12*dx2);

c1 = 11/12;
c2 = 14/3;
c3 = 9.5;
c4 = 26/3;
c5 = 35/12;
c6 = 5/3;
c7 = 11/12;
if (PERIODIQUE==0)
   d2Fdx(1)  = (c1*F(5) - c2*F(4) + c3*F(3) - c4*F(2) + c5*F(1))/dx2;
   d2Fdx(nx) = (c1*F(nx-4) - c2*F(nx-3) + c3*F(nx-2) - c4*F(nx-1) + c5*F(nx))/dx2;
   d2Fdx(2)  = (-F(5)/12 + F(4)/3 + F(3)/2 - c6*F(2) + c7*F(1))/dx2;
   d2Fdx(nx-1) = (-F(nx-4)/12 + F(nx-3)/3 + F(nx-2)/2 - c6*F(nx-1) + c7*F(nx))/dx2;
else
   d2Fdx(1)    = (-30*F(1) + 16*(F(2)+F(nx)) - (F(nx-1)+F(3))) / (12*dx2);
   d2Fdx(2)    = (-30*F(2) + 16*(F(3)+F(1)) - (F(nx)+F(4))) / (12*dx2);
   d2Fdx(nx)   = (-30*F(nx) + 16*(F(1)+F(nx-1)) - (F(nx-2)+F(2))) / (12*dx2);
   d2Fdx(nx-1) = (-30*F(nx-1) + 16*(F(nx)+F(nx-2)) - (F(nx-3)+F(1))) / (12*dx2);
end
