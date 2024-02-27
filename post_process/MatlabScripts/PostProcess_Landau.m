%-----------------------------------------------------
%  file : PsotProcess_Landau.m
%  date : 2019-03-01
%-----------------------------------------------------
fontsize = 18;

% Approximate solution for k<<1
omega_th = sqrt(1+sqrt(1+12*kx0^2))/sqrt(2);
gamma_th = -sqrt(pi/2)*omega_th/abs(kx0)*exp(-omega_th^2/(2*kx0^2))/...
  (2*kx0^2/omega_th^3*(1+6*kx0^2/omega_th^2));
vphi_th  = omega_th/kx0;


disp('*********************************')
disp(['Wave vector k=',num2str(kx0)])
disp('Using "DispRelation_Landau.m" to solve the dispersion relation, enter:');
omega_str = input('   Real frequency [default: limit k<<1]:      omega = ','s');
gamma_str = input('   Imaginary frequency [default: limit k<<1]: gamma = ','s');
disp('*********************************')
if (~isempty(omega_str)); omega_th = str2num(omega_str); end
if (~isempty(gamma_str)); gamma_th = str2num(gamma_str); end
vphi_th  = sqrt(omega_th^2+gamma_th^2)/kx0;

%-----------------------------------------------------
%*** plot the potential energy ***
Enpot_th = Enpot(1)*cos(omega_th*time).^2.*exp(2*gamma_th*time);

figure(1);
set(newplot,'fontsize',fontsize)
semilogy(time,Enpot,'r+-',time,Enpot_th,'b--'); grid;
hold on;
  plot(time,Enpot_th,'-k');
hold off
xlabel('time'); ylabel('Potential energy')

%return
%-----------------------------------------------------
%*** plot energy & entropy conservation ***
Entot      = Enkin+Enpot;
Entot_init = (Enkin(1)+Enpot(1))/2.;

figure(2);
set(newplot,'fontsize',fontsize)
subplot(211)
  plot(time,Enkin,'r+-',time,Enpot,'b-o');grid
  xlabel('time'); 
  hold on;
  plot(time,Entot-Entot_init,'g+-');
  hold off;
  legend('\delta E_{kin}','\delta E_{pot}','\delta E_{tot}');
subplot(212)
  semilogy(time,(entropy-entropy(1))./entropy,'r+-'); grid;
  xlabel('time'); ylabel('\delta Entropy / Entropy')

%-----------------------------------------------------
% real frequency - comparison to theoretical value
[FTPhi_t,omega] = Fourier1D(Phi1D_evol(1,:),time);

% index at which the spectrum is maximum
iomega_max = max(find(abs(FTPhi_t)==max(abs(FTPhi_t))));
omega_max  = omega(max(iomega_max));
disp('=========================')
disp(['    omega at max[FTPhi] = ',num2str(omega_max)])
disp('=========================')

figure(3)
set(newplot,'fontsize',fontsize)
subplot(211)
  plot(time,Phi1D_evol(1,:),'-r.');grid
  xlabel('time'); ylabel('phi(t,x=0)')
subplot(212)
  semilogy(omega,abs(FTPhi_t),'-r.');grid
  xlabel('frequency'); ylabel('|FT[phi(t,x=0)]|')
  hold on;
  plot(omega_max*[1 1],[min(abs(FTPhi_t)) max(abs(FTPhi_t))],'r--')
  plot(omega_th*[1 1],[min(abs(FTPhi_t)) max(abs(FTPhi_t))],'k--','linewidth',2)
  hold off
  legend('Spectrum','max of spectrum','theoretical frequency');


%-----------------------------------------------------
% Small scales in velocity space
dfv_evol = f1Dv_evol-f1Dv_evol(:,1)*ones(1,Ntime);

%figure(4)
for it=1:Ntime
  [FTdfv(it,:),kv] = Fourier1D(dfv_evol(:,it)',vg);
  %semilogy(kv,abs(FTdfv(it,:)),'-r.');grid;axis([kv(1) kv(end) 1e-10 1e-2])
  %xlabel('k_v');ylabel('|TF \deltaf|');
  %title(['it = ',num2str(it),' / ',num2str(Ntime)])
  %pause(0.05)
end
dkv = kv(2)-kv(1);

figure(4)
set(newplot,'fontsize',fontsize)
pcolor(kv-dkv/2,time,(abs(FTdfv)));shading('flat');colorbar
xlabel('k_v');ylabel('time')
hold on;
  plot(kx0*time+2*dkv,time,'k--')
  plot(-kx0*time-2*dkv,time,'k--')
hold off
title('Fourier Transform of \deltaf in velocity space')

%-----------------------------------------------------
if (f2D_evol_true==1) 
%-----------------------------------------------------
% f2D plot
fM   = exp(-vg.^2/2)/sqrt(2*pi);
fM2D = ones(Nx,1)*fM';

figure(5)
for itdiag=1:5:Ntime
  pcolor(xg,vg,(f2D_evol(:,:,itdiag)-f2D_evol(:,:,1))'); shading('flat'); colorbar();
  hold on;plot([xg(1) xg(end)],vphi_th*[1 1],'w');hold off
  xlabel('x coordinate'); ylabel('velocity')
  title(['\delta f(x,v) at time ',num2str(time(itdiag)),' / ',num2str(time(end))]) 
  pause(0.1)
end
%-----------------------------------------------------
end
%-----------------------------------------------------
