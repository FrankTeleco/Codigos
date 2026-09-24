% vamos a calibrar buscando el errore
clc;clear;close all;

fs= 200e3;
fc1= 20e3;
n= 0:fs-1;
Tono_1= exp(1j*2*pi*n*fc1/fs);

senyalTx= Tono_1(:);

senyalTx = 0.6*senyalTx/ max(abs(senyalTx)); % esto es para usar un 60% del rango del DAC

%% Setup pluto
fs_pluto= 200e3; % en Banda base
fc= 2.42e9; % frecuencia central de tranmision

tx= sdrtx("Pluto");

tx.CenterFrequency=fc;
tx.BasebandSampleRate=fs_pluto;
tx.Gain=0;
txRadioInfo= info(tx);


transmitRepeat(tx,senyalTx);

% bloque de recepcion

numero_muestras= 1*fs_pluto; % 5 segundos de señal debe ser impar y estar entre 2 y 16777216

rx= sdrrx("Pluto");
rx.CenterFrequency=fc;
rx.BasebandSampleRate= fs_pluto;
rx.SamplesPerFrame=numero_muestras;
rx.OutputDataType="double";

% recibir bloque

[senyalRx,~,overflow] = rx();

% comprobar overflow
if overflow
    warning("Se ha producido overflow en RX");
end


% procesamos lo que hemos recibido

figure(10)

plot(n(10:100)/fs,real(senyalRx(10:100)),'r')

hold on

plot(n(10:100)/fs,real(senyalTx(10:100)),'b')ue mas



[c,lags] = xcorr(senyalRx,senyalTx);

[~,idx] = max(abs(c));

retardo_muestras = lags(idx);

retardo_segundos = retardo_muestras/fs_pluto;

fprintf("Retardo = %d muestras\n",retardo_muestras);
fprintf("Retardo = %.3f us\n",retardo_segundos*1e6);