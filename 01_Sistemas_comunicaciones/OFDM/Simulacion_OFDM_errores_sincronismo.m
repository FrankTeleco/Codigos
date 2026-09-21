% OFDM Basico
% Este archivo contiene:
%   1) Canal con ruido
%   2) Simulacion para varios S/N
%   3) Evaluacion BER / SNR
%   4) Multicamino
%   5) Estimacion del canal y ecualizacion
%   6) Comparacion EQ ideal vs EQ con ruido
%   7) He añadido errores de sincronismo Temporal y offset en frecuencia
%      para ver sus efectos-
 

clc;
clear;
close all;

%% Parametros

MQAM = 16;
N = 64;
Ncp = N/4;

Ni = 48;          % Portadoras de datos
Np = 4;           % Portadoras piloto
Nu = 12;          % Portadoras nulas

NSimb_ofdm = 10;
Nbits = log2(MQAM)*Ni*NSimb_ofdm;

%SNR = 8:26;
SNR=200;
BerSnr_ruido = zeros(length(SNR),1);
BerSnr_ideal = zeros(length(SNR),1);

hcanal = [0.2 1 -0.2 0.5 0.1];

rng(10)

%% Generacion de bits

Bits_T = randi([0 1],Nbits,1);

%% Mapeo QAM

SimbQAM_T = qammod(Bits_T,MQAM,"bin","InputType","bit","UnitAveragePower",true);
SimbQAM_T1 = reshape(SimbQAM_T,Ni,[]);

%% Indices OFDM

DC = N*0.5 + 1;
Nyquist = -32 + N*0.5 + 1;

indices_nulos = [-32;-31;-30;-29;-28;-27;0;27;28;29;30;31] + N*0.5 + 1;
indices_pilotos = [-21;-7;7;21] + N*0.5 + 1;
indice_datos = setdiff(1:N,[indices_pilotos;indices_nulos]);

%% Pilotos

Pilotos = zeros(Np,NSimb_ofdm);

%% Generacion del preambulo

[preT,Xpre_da,Xpre_pi] = ref_ofdm(N,2*Ncp,indices_nulos,indices_pilotos);

%% Generacion de señal OFDM

SimbOfdm_1 = ofdmmod(SimbQAM_T1,N,Ncp,indices_nulos,indices_pilotos,Pilotos);

%% Estimacion ideal del canal

Entrenamiento_ideal = filter(hcanal,1,preT);

entrenamiento_ideal_1 = Entrenamiento_ideal(2*Ncp+1:2*Ncp+N);
entrenamiento_ideal_2 = Entrenamiento_ideal(2*Ncp+N+1:2*Ncp+2*N);

entrenamiento_ideal_pro = (entrenamiento_ideal_1 + entrenamiento_ideal_2)*0.5;

[Ypre_da_ideal,Ypre_pi_ideal] = ofdmdemod(entrenamiento_ideal_pro,N,0,0,indices_nulos,indices_pilotos);

Hest_da_ideal = Ypre_da_ideal ./ Xpre_da;
Hest_pi_ideal = Ypre_pi_ideal ./ Xpre_pi;

eq_ideal = 1 ./ Hest_da_ideal;
EQ_ideal = repmat(eq_ideal,1,NSimb_ofdm);

%% Bucle SNR


    %% Canal multicamino

    EntradaCanal = SimbOfdm_1;
    senyalCanal = filter(hcanal,1,EntradaCanal);
    senyalCanal = awgn(senyalCanal,SNR,"measured");
    senyalCanal= [senyalCanal(2:end);0];
    % Offset ebn freuencia
    CFO = 1/N*0.05;% 
    Ls = length(senyalCanal);

    senyalCanal = senyalCanal .* exp(1j*2*pi*CFO*(1:Ls)');
    Salida_Canal = senyalCanal;

    %% Preambulo con ruido

    Entrenamiento_canal = filter(hcanal,1,preT);
    Entrenamiento_canal = awgn(Entrenamiento_canal,SNR,"measured");
    %Entrenamiento_canal=[zeros(5,1);Entrenamiento_canal(1:end-5)];
    Entrenamiento_canal= [Entrenamiento_canal(2:end);0];

    %% Recepcion

    senyalRx = Salida_Canal;
    Senyal_entrenamiento_RX = Entrenamiento_canal;

    %% Estimacion del canal con ruido

    entrenamiento_RX_1 = Senyal_entrenamiento_RX(2*Ncp+1:2*Ncp+N);
    entrenamiento_RX_2 = Senyal_entrenamiento_RX(2*Ncp+N+1:2*Ncp+2*N);

    entrenamiento_pro = (entrenamiento_RX_1 + entrenamiento_RX_2)*0.5;

    [Ypre_da,Ypre_pi] = ofdmdemod(entrenamiento_pro,N,0,0,indices_nulos,indices_pilotos);

    Hest_da = Ypre_da ./ Xpre_da;
    Hest_pi = Ypre_pi ./ Xpre_pi;

    eq = 1 ./ Hest_da;
    EQ = repmat(eq,1,NSimb_ofdm);

    %% Demodulacion OFDM

    [SimbOfdm_R,Pilotos_R] = ofdmdemod(senyalRx,N,Ncp,Ncp,indices_nulos,indices_pilotos);

    %% Ecualizacion con estimacion ruidosa

    SimbOfdm_R_EQ_ruido = SimbOfdm_R .* EQ;

    Bits_R_ruido = qamdemod(SimbOfdm_R_EQ_ruido(:),MQAM,"bin","OutputType","bit","UnitAveragePower",true);

    [err_ruido,Ber_ruido] = biterr(Bits_T,Bits_R_ruido);









%% Constelacion

scatterplot(SimbOfdm_R_EQ_ruido(:))
title('Constelacion ecualizada')



%% Funcion preambulo

function [preT,Xpre_da,Xpre_pi] = ref_ofdm(Nfft,Ncp,I_nu,I_pi)

rng(35)

Npi = length(I_pi);
Nda = Nfft - Npi - length(I_nu);

Xpre_da = randi([0 1],Nda,1)*2 - 1;
Xpre_pi = randi([0 1],Npi,1)*2 - 1;

ref = ofdmmod(Xpre_da,Nfft,Ncp,I_nu,I_pi,Xpre_pi);

preT = [ref; ref(end-Nfft+1:end)];

end