clc;clear;close all;
archivo = 'LaBachata.mp3';
inicio= 30;
duracion= 10;
K= 3;
generadores=[7 5];

[bits_codificados, bits_originales, Fs, numMuestras, numCanales] = cancion_a_bits_conv(archivo,inicio,duracion,K,generadores);

audio_recuperado = bits_conv_a_cancion(bits_codificados,Fs,numMuestras,numCanales,K,generadores);

sound(audio_recuperado, Fs);




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

fprintf('Bits originales: %d\n', length(bits_originales));
fprintf('Bits codificados: %d\n', length(bits_codificados));
fprintf('Muestras: %d\n', numMuestras);
fprintf('Canales: %d\n', numCanales);

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