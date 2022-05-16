clc
clear
Fs = 1000;            % Sampling frequency                    
Ts = 1/Fs;             % Sampling period  
Tc=Ts*20000;
fc=1/Tc;
L = 2000000;             % Length of signal
t = (0:L-1)*Ts;        % Time vector


bits=5;
niveles=2^bits-1;   %-2^(bits-1)=-16 a  2^(bits-1)-1=15

uu = (2^bits-1)/2*sin(2*pi*fc*t); %señal  sinusoide de fc Hz 
u=ceil(uu);


 figure;plot(t(1:30000),uu(1:30000),'k.-')
 hold on
% for i=1:niveles
% nn(i,:)=i*ones(1,length(uu));
% plot(t,nn(i,:),'y--')
% end
 plot(t(1:30000),u(1:30000),'r.-')
% xlabel('t (milliseconds)')

%Analizo el espectro-------------------------
f = Fs*(0:(L/2))/L;
uu_fft=fft(uu);

P2 = abs(uu_fft/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);


u_fft=fft(u);
P2_d = abs(u_fft/L);
P1_d = P2_d(1:L/2+1);
P1_d(2:end-1) = 2*P1_d(2:end-1);
 figure;plot(f,20*log10(P1),'k.-')   %Grafica sin 20log10
 hold on;plot(f,20*log10(P1_d),'r.-') 
 xlabel('f [Hz]')
 ylabel('Potencia [dB]')


% ii=0;
% for ti=1:length(P1_d)
%     if P1_d(ti)~=0
%         ii=ii+1;
%        P1_d_sin_ceros(ii)=P1_d(ti);
%        ff_sin_ceros(ii)=f(ti);
%     end
% end


% figure;plot(ff_sin_ceros,20*log10(P1_d_sin_ceros),'r*-')   %señal digitalizada, le saco los ceros porq sino se me iba a INF cuando le hago log10
% hold on;plot(f,20*log10(P1),'ko-') 
% xlabel('f (Hz)')
 
M=3;  %cantidad de ADCs
AM=1;   %ADC extra para la aleatorizacion
MAM=M+AM;   %cantidad de ADCs totales
TTS=200;

Ag=randn(1,M+AM)*0.1;%error de amplitud de cada ADC
Ao=randn(1,M+AM)*0.1;%error de offset de cada ADC
At=floor(rand(1,M+AM)*TTS/8);%error de jitter de cada ADC, error en muestreo



sec7 = comm.PNSequence('Polynomial','x^3+x+1', ...
    'InitialConditions',[0 0 1],'SamplesPerFrame',1);
% sec15 = comm.PNSequence('Polynomial','x^4+x^3+1', ...
%     'InitialConditions',[0 0 0 1],'SamplesPerFrame',1);
% sec63 = comm.PNSequence('Polynomial','x^6+x^5+1', ...
%     'InitialConditions',[0 0 0 1 0 1],'SamplesPerFrame',1);
% sec63 = comm.PNSequence('Polynomial','x^7+x+1', ...
%     'InitialConditions',[0 1 0 0 1 0 1 ],'SamplesPerFrame',1);
% sec255 = comm.PNSequence('Polynomial','x^8+x^7+x^2+x+1', ...
%     'InitialConditions',[0 1 0 0 1 0 1 0],'SamplesPerFrame',1);
% sec511 = comm.PNSequence('Polynomial','x^9+x^5+1', ...
%     'InitialConditions', [0 0 0 1 0 0 0 1 0],'SamplesPerFrame',1);
% sec8191 = comm.PNSequence('Polynomial','x^13+x^12+x^10+x^9+1', ...
%     'InitialConditions', [0 0 1 0 0 0 0 1 0 0 0 1 0],'SamplesPerFrame',1);
secmucho = comm.PNSequence('Polynomial','x^40+x^5+x^4+x^3+1', ...
    'InitialConditions',[0 0 0 1 0 0 0 1 0 0 0 1 1 1 0 0 0 1 0 0 0 1 1 1 0 0 0 1 0 0 0 1 1 1 0 0 0 1 0 0],'SamplesPerFrame',1);


mem1_7=2;   %registro donde guarda el ADC que se eligio el ciclo anterior
mem2_7=3;    %registro donde guarda el ADC que se eligio 2 ciclos antes
elige0_7=4;  %registro donde guarda el ADC que se peude elegir
elige1_7=1;   %registro donde guarda el ADC que se peude elegir

% mem1_15=2;   %registro donde guarda el ADC que se eligio el ciclo anterior
% mem2_15=3;    %registro donde guarda el ADC que se eligio 2 ciclos antes
% elige0_15=4;  %registro donde guarda el ADC que se peude elegir
% elige1_15=1;   %registro donde guarda el ADC que se peude elegir
% 
% %repite todo para cada secuencia
% mem1_63=2;  %registro donde guarda el ADC que se eligio el ciclo anterior
% mem2_63=3;   %registro donde guarda el ADC que se eligio 2 ciclos antes
% elige0_63=4;  %registro donde guarda el ADC que se peude elegir
% elige1_63=1;   %registro donde guarda el ADC que se peude elegir
% 
% %repite todo para cada secuencia 511
% mem1_511=2;  %registro donde guarda el ADC que se eligio el ciclo anterior
% mem2_511=3;   %registro donde guarda el ADC que se eligio 2 ciclos antes
% elige0_511=4;  %registro donde guarda el ADC que se peude elegir
% elige1_511=1;   %registro donde guarda el ADC que se peude elegir
% 
% %repite todo para cada secuencia 8191
% mem1_8191=2;  %registro donde guarda el ADC que se eligio el ciclo anterior
% mem2_8191=3;   %registro donde guarda el ADC que se eligio 2 ciclos antes
% elige0_8191=4;  %registro donde guarda el ADC que se peude elegir
% elige1_8191=1;   %registro donde guarda el ADC que se peude elegir

%repite todo para cada secuencia 1.099511627775000e+12
mem1_m=2;  %registro donde guarda el ADC que se eligio el ciclo anterior
mem2_m=3;   %registro donde guarda el ADC que se eligio 2 ciclos antes
elige0_m=4;  %registro donde guarda el ADC que se peude elegir
elige1_m=1;   %registro donde guarda el ADC que se peude elegir

%repite todo para cada secuencia caotica
mem1_c=2;  %registro donde guarda el ADC que se eligio el ciclo anterior
mem2_c=3;   %registro donde guarda el ADC que se eligio 2 ciclos antes
elige0_c=4;  %registro donde guarda el ADC que se peude elegir
elige1_c=1;   %registro donde guarda el ADC que se peude elegir

ind=1;
indi=1;
indo=1;
indu=1;
inda=1;
inde=1;
seq_7 =sec7(); 
% seq_15 =sec15();    
% seq_63 =sec63();
% seq_511 =sec511();    
% seq_8191 =sec8191();
seq_m =secmucho();
ci=1;
seq_c =TWBM;
ss(ci)=seq_c;
ci=ci+1;
    if (seq_c<0.5)
        seq_c=0;
    else
        seq_c=1;
    end
%*****       %actualiza los registros usando sec7
        if (seq_7==0) 
                ind_7=elige0_7;
                elige0_7=mem2_7;
            else
                 ind_7=elige1_7;
              elige1_7=mem2_7;
            end
                mem2_7=mem1_7;
                mem1_7=ind_7;
                
%************   %actualiza los registros usando sec_m LSFR muiy largo
 
if (seq_m==0) 
    ind_m=elige0_m;
    elige0_m=mem2_m;
else
    ind_m=elige1_m;
    elige1_m=mem2_m;
end
mem2_m=mem1_m;
mem1_m=ind_m;

%************   %actualiza los registros usando sec_c caotico mapa
 
if (seq_c==0) 
    ind_c=elige0_c;
    elige0_c=mem2_c;
else
    ind_c=elige1_c;
    elige1_c=mem2_c;
end
mem2_c=mem1_c;
mem1_c=ind_c;

for k=1:1:L


 %muestreada con 1 solo ADC
if (k<TTS*indi)  %retengo el valor durante todo TTS
    yy(k)=u(indi*TTS);  %señal original con solo el ruido de cuantificacion
else
    yy(k)=u(indi*TTS); 
    indi=indi+1;
    seq_7 =sec7(); 
        %*****       %actualiza los registros usando sec7
 end    

 %muestreada con 4 ADCs secuenciales cada uno con su error de amplitud,
 %tiempo y fase que se repiten en la misma muestra
if (k<indo*TTS+At(mod(indo-1,MAM)+1))  %retengo el valor durante todo TTS
    y(k)=(1+Ag(mod(indo-1,MAM)+1))*u(indo*TTS+At(mod(indo-1,MAM)+1))+Ao(mod(indo-1,MAM)+1);  
else
    y(k)=(1+Ag(mod(indo-1,MAM)+1))*u(indo*TTS+At(mod(indo-1,MAM)+1))+Ao(mod(indo-1,MAM)+1); 
    indo=indo+1;
end  

 %muestreada con 4 ADCs seleccionados con secuencia LSR de periodo 7
if (k<indu*TTS+At(ind_7))  
    y_7(k)=(1+Ag(ind_7))*u(indu*TTS+At(ind_7))+Ao(ind_7);  
else
     y_7(k)=(1+Ag(ind_7))*u(indu*TTS+At(ind_7))+Ao(ind_7); 
    indu=indu+1;
    %*****       %actualiza los registros usando sec7
       seq_7 =sec7(); 
        if (seq_7==0) 
                ind_7=elige0_7;
                elige0_7=mem2_7;
            else
                 ind_7=elige1_7;
              elige1_7=mem2_7;
            end
         mem2_7=mem1_7;
         mem1_7=ind_7;
end  

 %muestreada con 4 ADCs seleccionados con secuencia LSR de periodo muy
 %largo
if (k<inda*TTS+At(ind_m))  
    y_m(k)=(1+Ag(ind_m))*u(inda*TTS+At(ind_m))+Ao(ind_m);  
else
     y_m(k)=(1+Ag(ind_m))*u(inda*TTS+At(ind_m))+Ao(ind_m); 
    inda=inda+1;

%************   %actualiza los registros usando secm
    seq_m =secmucho(); 
    if (seq_m==0) 
        ind_m=elige0_m;
         elige0_m=mem2_m;
    else
        ind_m=elige1_m;
        elige1_m=mem2_m;
    end
        mem2_m=mem1_m;
        mem1_m=ind_m;    
    
    
end  


 %muestreada con 4 ADCs seleccionados con secuencia caotica
 %largo
if (k<inde*TTS+At(ind_c))  
    y_c(k)=(1+Ag(ind_c))*u(inde*TTS+At(ind_c))+Ao(ind_c);  
else
     y_c(k)=(1+Ag(ind_c))*u(inde*TTS+At(ind_c))+Ao(ind_c); 
    inde=inde+1;

%************   %actualiza los registros usando secm
seq_c =rand();
ss(ci)=seq_c;
ci=ci+1;
    if (seq_c<0.5)
        seq_c=0;
    else
        seq_c=1;
    end
    
    
    if (seq_c==0) 
        ind_c=elige0_c;
         elige0_c=mem2_c;
    else
        ind_c=elige1_c;
        elige1_c=mem2_c;
    end
        mem2_c=mem1_c;
        mem1_c=ind_c;    
    
    
end  

% t_7(ind)=t(k+At(ind_7));
% y_15(ind)=(1+Ag(ind_15))*u(k+At(ind_15))+Ao(ind_15);   %señal muestreada eligiendo aleatoriamente con PRNG aleatorio,mas ruido de cuantificacion
% t_15(ind)=t(k+At(ind_15));
%  y_63(ind)=(1+Ag(ind_63))*u(k+At(ind_63))+Ao(ind_63);   %señal muestreada eligiendo aleatoriamente con PRNG caotico,mas ruido de cuantificacion
%  t_63(ind)=t(k+At(ind_63));
%  y_511(ind)=(1+Ag(ind_511))*u(k+At(ind_511))+Ao(ind_511);
%   t_63(ind)=t(k+At(ind_511));
%  y_8191(ind)=(1+Ag(ind_8191))*u(k+At(ind_8191))+Ao(ind_8191);
%   t_8191(ind)=t(k+At(ind_8191));
%   y_m(ind)=(1+Ag(ind_m))*u(k+At(ind_m))+Ao(ind_m);
%     t_m(ind)=t(k+At(ind_m));
%  
% ind=ind+TTS;
end

kk=30000;
figure;
plot(yy(1:kk),'.-');hold on   % 1 ADC
plot(y(1:kk),'r.-')          % 4 ADCS secuanciales
plot(y_7(1:kk),'g.-')        % 4 ADCs elige secuencia con LSR de periodo 7
plot(y_m(1:kk),'m.-')        % 4 ADCs elige secuencia con LSR de periodo muiy largo
plot(y_c(1:kk),'c.-')        % 4 ADCs elige secuencia con caotica
legend('1 ADC','4 ADC secuencial','4 ADC LSR 7','4 ADC LSR muy largo','4 ADC SEC CAOTICA')

figure;subplot(2,1,1);plot(yy(1:kk),'.-')
hold on;subplot(2,1,2);plot(y(1:kk),'r.-')
xlabel('tiempo [seg.]')
ylabel('tensión [V]')
subplot(2,1,1),ylabel('tensión [V]')

kk=60000;
figure;subplot(2,1,1);plot(y(1:kk),'r.-')
hold on;subplot(2,1,2);plot(y_m(1:kk),'m.-')
xlabel('tiempo [seg.]')
ylabel('tensión [V]')
subplot(2,1,1),ylabel('tensión [V]')


% figure; plot(yy,'.-')
% hold on; plot(y,'r.-')
% 
%  figure
%  plot(t_yy(1:100),yy(1:100),'k.-') %grafica en en negro señal cuantificada
%  hold on
%  plot(t_y(1:100),y(1:100),'o-') %grafica en azul la señal muestreada secuencialmente con cada ADC
%  plot(t_7(1:100),y_7(1:100),'r*-') %grafica en azul la señal muestreada secuencialmente con cada ADC
% % hold on
% plot(yy(1:1000),'r.-')  %grafica en rojo la señal original solo afectada con el ruido de cuantificacion
% plot(y_15(1:100),'k.-')   %grafica en negro la señal muestreada aleatoriamente con PRNG aleatorio
% plot(y_m(1:100),'m.-')   %grafica en celeste la señal muestreada aleatoriamente con PRNG caotico
% legend('señal muestreada secuencialmente con cada ADC','señal original solo afectada con el ruido de cuantificacion','señal muestreada con PRNG aleatorio','señal muestreada con PRNG caotico');
% 
figure
L=length(yy);
yy_fft=fft(yy);
P2_yy = abs(yy_fft/L);
P1_yy = P2_yy(1:L/2+1);
P1_yy(2:end-1) = 2*P1_yy(2:end-1);
hold on;plot(f(1:length(P1_yy)),20*log10(P1_yy),'k.-')  

y_fft=fft(y);  % 4 ADCS secuanciales
P2_y = abs(y_fft/L);
P1_y = P2_y(1:L/2+1);
P1_y(2:end-1) = 2*P1_y(2:end-1);
hold on;plot(f(1:length(P1_y)),20*log10(P1_y),'ro-')  


figure;subplot(2,1,1);plot(f(1:length(P1_yy)),20*log10(P1_yy),'k.-')  
subplot(2,1,1);axis([0.008 4.3 -42 25])
hold on;subplot(2,1,2);plot(f(1:length(P1_y)),20*log10(P1_y),'ro-')
subplot(2,1,2);axis([0.008 4.3 -42 25])
xlabel('frecuencia [Hz]')
ylabel('Potencia [dB]')
subplot(2,1,1),ylabel('Potencia [dB]')


y7_fft=fft(y_7);
P2_y7 = abs(y7_fft/L);
P1_y7 = P2_y7(1:L/2+1);
P1_y7(2:end-1) = 2*P1_y7(2:end-1);
hold on;plot(f(1:length(P1_y7)),20*log10(P1_y7),'g*-')  

ym_fft=fft(y_m);
P2_ym = abs(ym_fft/L);
P1_ym = P2_ym(1:L/2+1);
P1_ym(2:end-1) = 2*P1_ym(2:end-1);
hold on;plot(f(1:length(P1_ym)),20*log10(P1_ym),'m*-')  

yc_fft=fft(y_c);
P2_yc = abs(yc_fft/L);
P1_yc = P2_yc(1:L/2+1);
P1_yc(2:end-1) = 2*P1_yc(2:end-1);
hold on;plot(f(1:length(P1_yc)),20*log10(P1_yc),'c.-')  


figure;subplot(3,1,1);plot(f(1:length(P1_y)),20*log10(P1_y),'r.-')
subplot(3,1,1);axis([0.008 4.3 -42 25])
hold on;subplot(3,1,2);plot(f(1:length(P1_y7)),20*log10(P1_y7),'g.-')  
subplot(3,1,2);axis([0.008 4.3 -42 25])
hold on;subplot(3,1,3);plot(f(1:length(P1_ym)),20*log10(P1_ym),'m.-')  
subplot(3,1,3);axis([0.008 4.3 -42 25])
xlabel('frecuencia [Hz]')
ylabel('Potencia [dB]')
subplot(3,1,1),ylabel('Potencia [dB]')
subplot(3,1,2),ylabel('Potencia [dB]')

%---------------------------------------------------------------------------

Y_fft = abs(fft((y/NN)));  %FFT de la señal muestreada secuencialmente con cada ADC
P1=20*log10(Y_fft);

YY_fft = abs(fft((yy/NN)));
PP1=20*log10(YY_fft);

Y_7_fft = abs(fft((y_7/NN)));   %FFT de la señal muestreada aleatoriamente con PRNG aleatorio
P_7_1=20*log10(Y_7_fft);

Y_15_fft = abs(fft((y_15/NN)));   %FFT de la señal muestreada aleatoriamente con PRNG aleatorio
P_15_1=20*log10(Y_15_fft);

Y_63_fft = abs(fft((y_63/NN)));   %FFT de la señal muestreada aleatoriamente con PRNG aleatorio
P_63_1=20*log10(Y_63_fft);

Y_511_fft = abs(fft((y_511/NN)));   %FFT de la señal muestreada aleatoriamente con PRNG aleatorio
P_511_1=20*log10(Y_511_fft);

Y_8191_fft = abs(fft(double(y_8191/NN)));   %FFT de la señal muestreada aleatoriamente con PRNG aleatorio
P_8191_1=20*log10(Y_8191_fft);

Y_m_fft = abs(fft((y_m/NN)));   %FFT de la señal muestreada aleatoriamente con PRNG aleatorio
P_m_1=20*log10(Y_m_fft);

figure
plot(f(2:length(P1)/2),PP1(2:length(P1)/2),'ko-') %grafica FFT de la señal muestreada con ADC ideal, solo ruido cuantificacion
hold on
plot(f(2:length(P1)/2),P1(2:length(P1)/2),'*-') %grafica FFT de la señal muestreada secuencialmente con cada ADC
hold on
%plot(f,PP1,'r.-')
plot(f(2:length(P1)/2),P_7_1(2:length(P1)/2),'m>-') 
plot(f(2:length(P1)/2),P_15_1(2:length(P1)/2),'k*-')  %grafica FFT de la señal muestreada aleatoriamente con PRNG aleatorio
plot(f(2:length(P1)/2),P_63_1(2:length(P1)/2),'co-') %FFT de la señal muestreada aleatoriamente con PRNG caotico
plot(f(2:length(P1)/2),P_511_1(2:length(P1)/2),'g.-') %FFT de la señal muestreada aleatoriamente con PRNG caotico
plot(f(2:length(P1)/2),P_8191_1(2:length(P1)/2),'m.-') %FFT de la señal muestreada aleatoriamente con PRNG caotico
plot(f(2:length(P1)/2),P_m_1(2:length(P1)/2),'r.-') %FFT de la señal muestreada aleatoriamente con PRNG caotico
legend('ADC ideal','secuencial','LSR T=7','LSR T=15','LSR T=63','LSR T=511','LSR T=8191','LSR T=muy alto');
xlabel('f (Hz)')
ylabel('dB')

