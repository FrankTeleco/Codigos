% vamos a calibrar buscando el errore
clc;clear;close all;

fs= 200e3;
fc1= 20e3;
fc2= 40e3;
fc3= 80e3;
n= 0:fs-1;
Tono_1= exp(1j*2*pi*n*fc1/fs);
Tono_2= exp(1j*2*pi*n*fc2/fs);
Tono_3= exp(1j*2*pi*n*fc3/fs);

senyal = Tono_1+Tono_2+Tono_3; senyal= senyal(:);

%pwelch(senyal,1024,[],[],fs,"centered"); % 3 tosnos

senyal = 0.6*senyal/ max(abs(senyal)); % esto es para usar un 60% del rango del DAC

%% Setup pluto
fs_pluto= 200e3; % en Banda base
fc= 2.42e9; % frecuencia central de tranmision

tx= sdrtx("Pluto");

tx.CenterFrequency=fc;
tx.BasebandSampleRate=fs_pluto;
tx.Gain=0;
txRadioInfo= info(tx);


transmitRepeat(tx,senyal);

% bloque de recepcion

numero_muestras= 5*fs_pluto; % 5 segundos de señal debe ser impar y estar entre 2 y 16777216

rx= sdrrx("Pluto");
rx.CenterFrequency=fc;
rx.BasebandSampleRate= fs_pluto;
rx.SamplesPerFrame=numero_muestras;
rx.ShowAdvancedProperties=true;
rx.FrequencyCorrection= -2.1018e-5;

% recibir bloque

[data,~,overflow] = rx();

% comprobar overflow
if overflow
    warning("Se ha producido overflow en RX");
end


% procesamos lo que hemos recibido

Nfft= 2^18;

% FFT
Xf = fftshift(fft(data,Nfft));
P = 20*log10(abs(Xf)/max(abs(Xf))); % potencia normalizada 
f= (-Nfft*0.5:Nfft*0.5-1)* fs_pluto/Nfft;

% visualizacion
figure;
plot(f/1e3,P);
grid on;
xlabel("Frecuencia en banda base [kHz]");
ylabel("Magnitud [dB]");
title("Espectro recibido con ADALM-PLUTO");
xlim([-fs_pluto/2 fs_pluto/2]/1e3);

% calculamos el error 

% buscamos los picvos

[picos,posicion]=findpeaks(abs(Xf),"SortStr","descend","NPeaks",3);

f_medidas = sort(f(posicion));

f_ideales = [20e3 40e3 80e3];

error_Hz = f_medidas - f_ideales;
error_Hz= mean(error_Hz);
ppm = error_Hz/fc * 1e6;