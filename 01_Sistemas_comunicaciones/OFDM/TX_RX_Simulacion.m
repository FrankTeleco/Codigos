clc;clear;close all;




%% Generador de bits


%% Parametros
MQAM = 16;
N= 64;
Ncp= N/4;
Ni= 48; % portadoras datos
Np= 4; % Portadoras pilotos
Nu= 12; % Portadoras nulas


%% Generacion de bits

archivo = 'LaBachata.mp3';
inicio= 30;
duracion= 10;
K= 3;
generadores=[7 5];

[bits_codificados, bits_originales, Fs, numMuestras, numCanales] = cancion_a_bits_conv(archivo,inicio,duracion,K,generadores);
Bits_T=bits_codificados(:);

NSimb_ofdm= ceil(length(Bits_T)/(Ni*log2(MQAM)));

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




% 
% audio_recuperado = bits_conv_a_cancion(bits_codificados,Fs,numMuestras,numCanales,K,generadores);
% 
% sound(audio_recuperado, Fs);




function [bits_codificados, bits_originales, Fs, numMuestras, numCanales] = ...
    cancion_a_bits_conv(archivo, inicio, duracion, K, generadores)

% Información del audio
info = audioinfo(archivo);
Fs = info.SampleRate;

% Seleccionar fragmento
muestraInicio = floor(inicio * Fs) + 1;
muestraFin = floor((inicio + duracion) * Fs);

muestraFin = min(muestraFin, info.TotalSamples);

[audio, Fs] = audioread(archivo, ...
    [muestraInicio muestraFin]);

% Guardar dimensiones
[numMuestras, numCanales] = size(audio);

% Limitar rango
audio = max(min(audio, 1), -1);

% Pasar a PCM de 16 bits
audio_int16 = int16(audio * 32767);

% Vectorizar
muestras = audio_int16(:);

% Reinterpretar como uint16
muestras_uint16 = typecast(muestras, 'uint16');

% Pasar cada muestra a 16 bits
bits_matriz = de2bi(muestras_uint16, ...
    16, ...
    'left-msb');

% Convertir a vector de bits
bits_originales = bits_matriz.';
bits_originales = bits_originales(:).';

% Crear código convolucional
trellis = poly2trellis(K, generadores);

% Añadir bits de terminación
memoria = K - 1;

bits_entrada = [bits_originales ...
    zeros(1, memoria)];

% Codificar
bits_codificados = convenc(bits_entrada, trellis);


end


function audio_recuperado = bits_conv_a_cancion( ...
    bits_codificados, Fs, numMuestras, numCanales, K, generadores)

%% 1. Crear el mismo código convolucional

trellis = poly2trellis(K, generadores);

%% 2. Decodificar con Viterbi

traceback = 5*(K-1);

bits_decodificados = vitdec( ...
    bits_codificados, ...
    trellis, ...
    traceback, ...
    'term', ...
    'hard');

%% 3. Quitar los bits de terminación

memoria = K - 1;

bits_recuperados = bits_decodificados(1:end-memoria);

%% 4. Convertir grupos de 16 bits en muestras

bits_matriz = reshape(bits_recuperados, 16, []).';

muestras_uint16 = bi2de(bits_matriz, 'left-msb');

%% 5. Reinterpretar uint16 como int16

muestras_int16 = typecast( ...
    uint16(muestras_uint16), ...
    'int16');

%% 6. Recuperar la forma original del audio

audio_int16 = reshape( ...
    muestras_int16, ...
    numMuestras, ...
    numCanales);

%% 7. Convertir a double

audio_recuperado = double(audio_int16) / 32767;

%% 8. Guardar canción recuperada

audiowrite('cancion_recuperada.wav', ...
    audio_recuperado, Fs);

fprintf('Cancion recuperada correctamente.\n');
fprintf('Guardada como cancion_recuperada.wav\n');

end