% OFDM Basico
% Plantilla basica de Tx Rx con generacion manual y automatica
clc;clear;close all;

%% Parametros
MQAM = 16;
N= 64;
Ncp= N/4;
Ni= 48; % portadoras datos
Np= 4; % Portadoras pilotos
Nu= 12; % Portadoras nulas
NSimb_ofdm= 1;
Nbits= log2(MQAM)*Ni*NSimb_ofdm;

%% Generacion de bits

Bits_T= randi([0 1], Nbits,1);

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
Pilotos = ones(Np,NSimb_ofdm); % Aqui me guardo los pilotos

SimbOfdm_1 = ofdmmod(SimbQAM_T1,N,Ncp,indices_nulos,indices_pilotos,Pilotos);

%% Generacion Manual

SimbOfdm= zeros(N,NSimb_ofdm);

SimbOfdm(indice_datos,:)=SimbQAM_T1;
SimbOfdm(indices_pilotos,:)=Pilotos;
SimbOfdm(indices_nulos,:)=0;

SimbOfdm = ifft(ifftshift(SimbOfdm,1), N, 1);% ifftshift dezoplaza el vector N/2 posiciones para que quede DC positivas Negativas

SimbOfdm= [SimbOfdm(end-Ncp+1:end,:);SimbOfdm];

SimbOfdm= SimbOfdm(:);

%z=isequal(SimbOfdm, SimbOfdm_1); % comprobamos que ambas generaciones son iguales


%-------------------------- Receptor

%% Recepcion

senyalRx= SimbOfdm;

SimbOfdm_R= reshape(senyalRx,N+Ncp,[]);
SimbOfdm_R= SimbOfdm_R(Ncp+1:end,:);
SimbOfdm_R= fft(SimbOfdm_R,N);
SimbOfdm_R = fftshift(SimbOfdm_R,1); % vuelvo a la  misma posiciosn
SimbOfdm_R= SimbOfdm_R(indice_datos,:);

SimbOfdm_R=SimbOfdm_R(:);

Bits_R =qamdemod(SimbOfdm_R,MQAM,"bin","OutputType","bit","UnitAverage",true);

%% Evalaucion

[Ber, erro]= biterr(Bits_T,Bits_R)

