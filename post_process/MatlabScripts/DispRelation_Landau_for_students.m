clear
No = 1000;
Ng = 801;

Disp_real = zeros(Ng,No);
Disp_imag = zeros(Ng,No);

gamma_max = 0.1; 
gamma_min = -3.5; 
omega_max = 3.5;
omega_min = 0.5;
pas_o = (omega_max-omega_min)/No;
pas_g = (gamma_max-gamma_min)/ (Ng-1);

omega = (omega_min+pas_o):pas_o:omega_max;
gamma = (gamma_min:pas_g:gamma_max);

%v   = -4:0.0001:4;
vmax = 7;
Nv   = 255;
dv   = 2*vmax/(Nv+2);
v    = -vmax:dv:vmax;
v2   = v.*v;
fv0  = exp(-0.5*v2)/sqrt(2*pi);

figure(1);
plot(v,fv0,'-r.',[0 0],[0 0.5],'k');axis([-vmax vmax 0 0.5]);grid
xlabel('v');ylabel('f(v)');
title(['distribution function'])

% Selection of the wave vector
kstring = input('Wave vector: [default: k=0.5] k = ','s');
if (~isempty(kstring))
  k = str2num(kstring);
else
  k=0.5;
end
k2 = k*k;

% Dispersion relation
z          = (ones(Ng,1)*omega + 1i.*(gamma')*ones(1,No))/(sqrt(2)*abs(k));
dispersion = k2 + 1 + z.*zetaf(z); 
D          = abs(dispersion);
Log_D      = -log(D);

Nmax = 40;
Nf1  = 15;
Nf2  = 5;
Nf3  = 10;
fscale = zeros(1,Nmax);
maxf = max(max(Log_D));
minf = 0.6;
pasf  = (maxf-minf)/(Nf1-1);
fscale((Nmax-Nf1+1):Nmax) = minf:pasf:maxf;
maxf = 0.5;
minf = -0.5;
pasf  = (maxf-minf)/(Nf2-1);
fscale((Nmax-Nf1-Nf2+1):(Nmax-Nf1)) = minf:pasf:maxf;
maxf =-2;
minf = min(min(Log_D));
pasf  = (maxf-minf)/(Nf3-1);
fscale(1:Nf3) = minf:pasf:maxf;

lwidth = 1;
figure(2);
set(newplot,'fontsize',16)
contour(omega,gamma,Log_D,fscale,'b','linewidth',lwidth);
xlabel('Real part of \omega');ylabel('Imaginary part of \omega');
title(['poles of the dispersion relation'])
hold on
  plot([omega_min omega_max],[0 0],'r')
hold off

[ig,io]=find(Log_D==max(max(Log_D)));
omega_sol = omega(io);
gamma_sol = gamma(ig);
display('=================================')
display('Solution of Dispersion Relation:')
display(['wave number k  = ',num2str(k)])
display(['real frequency = ',num2str(omega_sol)])
display(['growth rate    = ',num2str(gamma_sol)])
display(['phase velocity = ',num2str(omega_sol/k)])
display('=================================')
