% OFDM Basico
% Este archivo contiene:
%   1) Canal con ruido
%   2) Simulacion para varios S/N
%   3) Evalaucion BER /SNR

clc;clear;close all;

%% Parametros
MQAM = 16;
N= 64;
Ncp= N/4;
Ni= 48; % portadoras datos
Np= 4; % Portadoras pilotos
Nu= 12; % Portadoras nulas
NSimb_ofdm= 1000;
Nbits= log2(MQAM)*Ni*NSimb_ofdm;
SNR= 4:1:18; % evaluo de 1 hasta 30dB
BerSnr= zeros(length(SNR),1);%% Generacion de bits

Bits_T= randi([0 1], Nbits,1);
rng(35)

%% Mapeo

SimbQAM_T = qammod(Bits_T,MQAM,"bin","InputType","bit","UnitAveragePower",true);
SimbQAM_T1= reshape(SimbQAM_T,Ni,[]);% pongo asi porque esa sera generado con ofdmod

%% Generacion con ofdmmod

% Calculamos los indices

DC= N*0.5+1;% esto sale en 33
Nyquist= -32 + N*0.5+1; % esto sale a 1

indices_nulos= [-32;-31;-30;-29;-28;-27;0; 27;28;29;30;31] + N*0.5 +1; % Bandas de guarda
indices_pilotos=[-21;-7;7;21] + N*0.5 +1;
indice_datos= setdiff(1:N,[indices_pilotos;indices_nulos]);

% Pilotos
Pilotos = zeros(Np,NSimb_ofdm); % Aqui me guardo los pilotos

SimbOfdm_1 = ofdmmod(SimbQAM_T1,N,Ncp,indices_nulos,indices_pilotos,Pilotos);

%--------------------------Canal
for k= 1: length(SNR)
    EntradaCanal= SimbOfdm_1;
    senyalCanal = awgn(EntradaCanal,SNR(k),"measured");
    Salida_Canal= senyalCanal;

    %-------------------------- Receptor

    %% Recepcion

    senyalRx= Salida_Canal;

    [SimbOfdm_R,Pilotos_R]= ofdmdemod(senyalRx,N,Ncp,Ncp,indices_nulos,indices_pilotos);
    SimbOfdm_R=SimbOfdm_R(:);
    Bits_R =qamdemod(SimbOfdm_R,MQAM,"bin","OutputType","bit","UnitAveragePower",true);

    %--------- Evalaucion

    %% BER / SER
    [err, Ber]= biterr(Bits_T,Bits_R);
    BerSnr(k)=Ber;
end
%% Constelacion
scatterplot(SimbOfdm_R)

%% BER/ SNR
figure
semilogy(SNR,BerSnr,'o-')
grid on
xlabel('SNR [dB]')
ylabel('BER')
title('BER vs SNR - 16 QAM OFDM')