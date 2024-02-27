clear
No = 1000;
Ng = 801;

Disp_real = zeros(Ng,No);
Disp_imag = zeros(Ng,No);

gamma_min = -0.8; 
gamma_max = 0.5; 
omega_min = 0.5;
omega_max = 1.5;
%gamma_min = -0.1; 
%gamma_max = 0.2; 
%omega_min = 0.9;
%omega_max = 1.3;
pas_o = (omega_max-omega_min)/No;
pas_g = (gamma_max-gamma_min)/ (Ng-1);

omega = (omega_min+pas_o):pas_o:omega_max;
gamma = (gamma_min:pas_g:gamma_max);

%================================================
% Characteristics of distribution functions
%================================================
disp('*****************')
disp('Characteristics of the equilibrium distribution functions')
disp('   F1 = (1-eps) * exp{-v^2/2} / sqrt(2*pi)')
disp('   F2 = eps * exp{-(v-v0)^2/(2*T0)} / sqrt(2*pi*T0)')
disp('*****************')
eps_string = input('Value of epsilon: [default: eps=0.1] eps = ','s');
if (~isempty(eps_string))
  eps = str2num(eps_string);
else
  eps = 0.1;
end
v0_string = input('Value of v0: [default: v0=3.8] v0 = ','s');
if (~isempty(v0_string))
  v0 = str2num(v0_string);
else
  v0 = 3.8;
end
T0_string = input('Value of T0: [default: T0=0.2] T0 = ','s');
if (~isempty(T0_string))
  T0 = str2num(T0_string);
else
  T0 = 0.2;
end

% Velocity
vmax  = 8;
v     = -vmax:1.e-4:vmax;
v2    = v.*v;
vmv02 = (v-v0).*(v-v0);

% Distribution functions
f1 = (1-eps) * exp(-0.5*v2) / sqrt(2*pi);
f2 = eps * exp(-vmv02/(2*T0)) / sqrt(2*pi*T0);
ftot = f1 + f2;

figure();
subplot(211)
plot(v,f1,'r.',v,f2,'b.',v,ftot,'k');grid
xlabel('v');ylabel('f(v)');
title(['distribution function'])
subplot(212)
semilogy(v,f1,'r.',v,f2,'b.',v,ftot,'k');grid
axis([-vmax vmax min(f1) 2*max(ftot)]);
xlabel('v');ylabel('f(v)');
title(['distribution function'])

%================================================
% Dispersion relation
%================================================
% Selection of the wave vector
kstring = input('Wave vector: [default: k=0.377] k = ','s');
if (~isempty(kstring))
  k = str2num(kstring);
else
  k=0.377;
end
k2 = k*k;

z  = ones(Ng,1)*omega + 1i.*(gamma')*ones(1,No);
z1 = z / (sqrt(2)*abs(k));
z2 = (z-k*v0) / (sqrt(2*T0)*abs(k));

dispersion1 = (1-eps) * (1 + z1.*zetaf(z1)); 
dispersion2 = eps/T0 * (1 + z2.*zetaf(z2)); 
dispersion  = k2 + dispersion1 + dispersion2;

D     = abs(dispersion);
Log_D = -log(D);

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
figure();
set(newplot,'fontsize',20)
contour(omega,gamma,Log_D,fscale,'b','linewidth',lwidth);
xlabel('Real part of \omega');ylabel('Imaginary part of \omega');
title(['poles of the dispersion function'])
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
display(['phase velocity = ',num2str(sqrt(omega_sol^2 + gamma_sol^2)/k)])
display('=================================')
