%-----------------------------------------------------
%  file : PostProcess_BumpOnTail.m
%  date : 2015-02-15
%-----------------------------------------------------

%-----------------------------------------------------
% Approximate solution for k<<1
omega_th = sqrt(1+sqrt(1+12*kx0^2))/sqrt(2);
gamma_th = 0;
vphi_th  = omega_th/kx0;


disp('*********************************')
disp(['Wave vector k=',num2str(kx0)])
disp('Using "DispRelation_BoT.m" to solve the dispersion relation, enter:');
omega_str = input('   Real frequency [default: limit k<<1]:      omega = ','s');
gamma_str = input('   Imaginary frequency [default: 0]: gamma = ','s');
disp('*********************************')
if (~isempty(omega_str)) 
  omega_th = str2num(omega_str); 
  vphi_th  = sqrt(omega_th^2+gamma_th^2)/kx0;
end
if (~isempty(gamma_str))
  gamma_th = str2num(gamma_str);
end

%-----------------------------------------------------
% Fourier transform in space => kx
for it=1:Ntime
  [FTPhi(it,:),kx] = Fourier1D(Phi1D_evol(:,it)',xg);
end

figure(1)
subplot(211)
for it=1:fix(Ntime/20):Ntime
  %[FTPhi_tmp,kx] = Fourier1D(Phi1D_evol(:,it)',xg);
  plot(kx,log10(abs(FTPhi(it,:))),'-r.');grid
  axis([kx(1) kx(end) -10 0])
  xlabel('kx'); ylabel('abs(FT[phi])')
  title(['it = ',num2str(it),' / ',num2str(Ntime)])
  pause(0.05)
end
hold on;plot(kx,log10(abs(FTPhi(1,:))),'-b.');hold off

subplot(212)
pcolor(kx,time,log10(abs(FTPhi)));shading('flat');colorbar
xlabel('wave vector kx'); ylabel('time');
title('Log10 |F[phi]|')
caxis([-10 max(max(log10(abs(FTPhi))))])


%-----------------------------------------------------
% Fourier transform in time for the "most unstable" mode

% wave number having the largest magnitude at t_end 
ikx_max = max(find(abs(FTPhi(end,:))==max(abs(FTPhi(end,:)))));
kx_max  = kx(ikx_max);
% Time evolution of this mode - comparison to theoretical growth rate
figure(2);
subplot(211)
semilogy(time,abs(FTPhi(:,ikx_max)),'-r.');grid
title(['Mode kx = ',num2str(kx_max)]);
xlabel('time');ylabel('abs(FT[phi])')
hold on
  plot(time,abs(FTPhi(1,ikx_max)*exp(gamma_th*time).*cos(omega_th*time)),'k')
  %plot(time,abs(FTPhi(1,ikx_max)*exp(gamma_th*time))/5,'k')
hold off
subplot(212)
[FTPhi_kmax,freq] = Fourier1D(real(FTPhi(1:min(100,Ntime),ikx_max)'),time(1:min(100,Ntime)));
semilogy(freq,abs(FTPhi_kmax),'-r.');grid
xlabel('freq');ylabel('abs(FT[phi_{kxmax}])')

% index at which the spectrum is maximum
iomega_max = max(find(abs(FTPhi_kmax)==max(abs(FTPhi_kmax))));
omega_max  = freq(max(iomega_max));
disp('=========================')
disp(['    k_max     = ',num2str(kx_max)])
disp(['    omega_max = ',num2str(omega_max)])
disp('=========================')

hold on
  plot([1 1]*omega_max,[min(abs(FTPhi_kmax)) max(abs(FTPhi_kmax))],'k')
hold off

%-----------------------------------------------------
%*** plot energy conservation ***
Entot      = Enkin+Enpot;
Entot_init = (Enkin(1)+Enpot(1))/2.;
figure(3);

subplot(211)
plot(time,Enkin,'r+-',time,Enpot,'b-o');grid
hold on;
plot(time,Entot-Entot_init,'g+-');
hold off;
legend('Enkin','Enpot','Entot');

%*** plot entropy conservation ***
subplot(212)
semilogy(time,abs(entropy-entropy(1))./entropy,'r+-'); grid;
xlabel('time'); ylabel('|delta Entropy| / Entropy')

%-----------------------------------------------------
% equilibrium distribution function

% phase velocity of kx_max mode
vphi_max = omega_max/kx_max;

figure(4)
for it=1:fix(Ntime/20):Ntime
  semilogy(vg,f1Dv_evol(:,it),'-r.');grid
  hold on;
  plot([1 1]*vphi_max,[1e-5 max(f1Dv_evol(:,1))],'k')
  hold off
  xlabel('velocity v'); ylabel('Distribution Function')
  title(['it = ',num2str(it),' / ',num2str(Ntime)])
  %axis([vg(1) vg(end) 0 max(f1Dv_evol(:,1))])
  axis([vg(1) vg(end) 1.e-6 2])
  pause(0.05)
end

%-----------------------------------------------------
% f2D plot
fM   = exp(-vg.^2/2)/sqrt(2*pi);
fM2D = ones(Nx,1)*fM';

figure(5)
for itdiag=1:fix(Ntime/10):Ntime
  pcolor(xg,vg,(f2D_evol(:,:,itdiag)-f2D_evol(:,:,1))'); shading('flat'); colorbar();
  hold on;plot([xg(1) xg(end)],vphi_th*[1 1],'w');hold off
  xlabel('x coordinate'); ylabel('velocity')
  title(['\delta f(x,v) at time ',num2str(time(itdiag)),' / ',num2str(time(end))]) 
  pause(0.05)
end

%-----------------------------------------------------
% Small scales in velocity space

dfv_evol = f1Dv_evol-f1Dv_evol(:,1)*ones(1,Ntime);
for it=1:Ntime
  [FTdfv(it,:),kv] = Fourier1D(dfv_evol(:,it)',vg);
end

figure(6)
subplot(211)
for it=1:fix(Ntime/10):Ntime
  %[FTdfv_tmp,kv] = Fourier1D(dfv_evol(:,it)',vg);
  semilogy(kv,abs(FTdfv(it,:)),'-r.');grid;axis([kv(1) kv(end) 1e-10 1e-1])
  xlabel('k_v');ylabel('|TF \deltaf|');
  title(['it = ',num2str(it),' / ',num2str(Ntime)])
  pause(0.05)
end
dkv = kv(2)-kv(1);

subplot(212)
pcolor(kv-dkv/2,time,log10(abs(FTdfv)));shading('flat');colorbar
xlabel('k_v');ylabel('time')
title('Fourier Transform of \deltaf in velocity space (log scale)')
hold on;
  plot(kx0*time+2*dkv,time,'k--')
  plot(-kx0*time-2*dkv,time,'k--')
hold off

return

%-----------------------------------------------------
% resonance broadening: width of island

% Width of the island (separatrix)
delta_v = 4*sqrt(abs(FTPhi(:,ikx_max)));

figure(7);
dfv_evol = f1Dv_evol-f1Dv_evol(:,1)*ones(1,Ntime);
pcolor(vg,time,dfv_evol');shading('flat');colorbar
hold on;
plot(time*0,time,'k')
plot(vphi_max*ones(1,Ntime),time,'w');
plot(vphi_max+delta_v/2,time,'w--',vphi_max-delta_v/2,time,'w--');
hold off
xlabel('velocity');ylabel('time');title('f(x_0,v,t)-f(x_0,v,0)')


