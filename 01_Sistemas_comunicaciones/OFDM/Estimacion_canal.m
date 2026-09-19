clc;clear;close all;


%% Parametros
MQAM = 16;
N= 64;
Ncp= N/4;
Ni= 48; % portadoras datos
Np= 4; % Portadoras pilotos
Nu= 12; % Portadoras nulas
NSimb_ofdm= 100;
indices_nulos= [-32;-31;-30;-29;-28;-27;0; 27;28;29;30;31] + N*0.5 +1; % Bandas de guarda
indices_pilotos=[-21;-7;7;21] + N*0.5 +1;
hcanal = [0.2 1 -0.2 0.5 0.1];
SNR= 100;

[preT,Xpre_da,Xpre_pi] = ref_ofdm(N,2*Ncp,indices_nulos,indices_pilotos);

%---------------------Canal

EntradaCanal= preT;
senyalCanal= filter(hcanal,1,EntradaCanal);

senyalCanal = awgn(senyalCanal,SNR,"measured");

Salida_Canal= senyalCanal;

%------------------------Receptor

SenyalRx= Salida_Canal;

% primero localizo el preambulo
%Tiene tamanyo  2*Ncp + N +N en este caso 160

% busco los simbolos de entremaiento

Y_entrenamiento_uno= SenyalRx(2*Ncp+1:2*Ncp+N);
Y_entrenamiento_dos=SenyalRx(2*Ncp+1+N:2*Ncp+2*N);

Y_pro= (Y_entrenamiento_uno+Y_entrenamiento_dos)*0.5; % promediando nos quitamos el ruido



[Ypre_da,Ypre_pi] = ofdmdemod(Y_pro,N,0,0,indices_nulos,indices_pilotos); % Y_pro ya no tiene el prefijo ciclico  son los dos simbolos de entremaientos promediados

Hest_da = Ypre_da ./ Xpre_da; % simbolos de entremiento en frecuencia
Hest_pi = Ypre_pi ./ Xpre_pi; % simbolos de pilotos en frecuencia


Hcanal = fft(hcanal,N);

I_dat=setdiff(1:N,[indices_nulos; indices_pilotos]);

figure(10)
plot(1:64,fftshift(abs(Hcanal)),I_dat,abs(Hest_da),'*',indices_pilotos,abs(Hest_pi),'ok')






function [preT,Xpre_da,Xpre_pi] = ref_ofdm(Nfft,Ncp,I_nu,I_pi)
rng(35)
Npi = length(I_pi);
Nda = Nfft - Npi - length(I_nu);
Xpre_da = randi([0 1],Nda,1)*2 -1;
Xpre_pi = randi([0 1],Npi,1)*2 -1;
ref = ofdmmod(Xpre_da,Nfft,Ncp,I_nu,I_pi,Xpre_pi);
preT = [ref; ref(end-Nfft+1:end)];
end