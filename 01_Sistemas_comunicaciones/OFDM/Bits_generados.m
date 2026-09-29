clc; clear;close all;


archivo = 'LaBachata.mp3';

fid = fopen(archivo,'rb');
bytes = fread(fid,Inf,'*uint8');
fclose(fid);

Bits_T = int2bit(bytes,8,true);
Bits_T = Bits_T(:);
Bits_R= Bits_T;
% Bits_R = bits recibidos tras OFDM

bytes_R = bit2int(Bits_R,8,true);

fid = fopen('cancion_recuperada_2.mp3','wb');
fwrite(fid,bytes_R,'uint8');
fclose(fid);

[y,fs] = audioread('cancion_recuperada_2.mp3');
sound(y,fs);