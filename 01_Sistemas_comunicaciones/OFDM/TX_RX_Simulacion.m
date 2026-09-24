clc;clear;close all;


%% Parametros
MQAM = 16;
N= 64;
Ncp= N/4;
Ni= 48; % portadoras datos
Np= 4; % Portadoras pilotos
Nu= 12; % Portadoras nulas

%% Generacion de bits
%archivo = 'LaBachata.mp3';
archivo ='Karol G - Si antes te hubiera conocido (HQ).mp3';
inicio= 100;
duracion= 30;
Bits_T = cancion_a_bits_conv(archivo,inicio,duracion,Ni,MQAM);
NSimb_ofdm= round(length(Bits_T)/(Ni*log2(MQAM)));

%% Mapeo
SimbQAM_T = qammod(Bits_T,MQAM,"bin","InputType","bit","UnitAveragePower",true);
SimbQAM_T= reshape(SimbQAM_T,Ni,[]);% pongo asi porque esa sera generado con ofdmod

%% Generacion con ofdmmod
% Calculamos los indices
DC= N*0.5+1;% esto sale en 33
Nyquist= -32 + N*0.5+1; % esto sale a 1
indices_nulos= [-32;-31;-30;-29;-28;-27;0; 27;28;29;30;31] + N*0.5 +1; % Bandas de guarda
indices_pilotos=[-21;-7;7;21] + N*0.5 +1;
indice_datos= setdiff(1:N,[indices_pilotos;indices_nulos]);
% Pilotos
Pilotos = ones(Np,NSimb_ofdm); % Aqui me guardo los pilotos
SimbOfdm = ofdmmod(SimbQAM_T,N,Ncp,indices_nulos,indices_pilotos,Pilotos);
%% Genero los simbolos Pilotos y los añado
[preT,Xpre_da,Xpre_pi] = ref_ofdm(N,2*Ncp,indices_nulos,indices_pilotos);

SimbOfdm=[preT;SimbOfdm];

%% Canal

retardo = 100;
CFO = 0.15/N;
hcanal = [0.2 1 -0.2 0.5 0.1];
EntradaCanal= [zeros(retardo,1); SimbOfdm];
EntradaCanal= filter(hcanal,1,EntradaCanal);
n_canal= (0:numel(EntradaCanal)-1).';
EntradaCanal= EntradaCanal.*exp(1j*2*pi*CFO*n_canal);
%SenyalCanal= SimbOfdm;
SenyalCanal= awgn(EntradaCanal,30,"measured");

SalidaCanal= SenyalCanal;
%-------------------------- Receptor
%% Recepcion
senyalRx= SalidaCanal;

%% Sincronizacion Temporal

filtro_inversor= conj(preT(end:-1:1));

correlacion= filter(filtro_inversor,1,senyalRx);

[valor,posicion]= max(abs(correlacion));
% Adelantar la ventana dentro del CP para absorber el desplazamiento del pico por el canal
posicion= posicion-Ncp/2;

% Extraemos los simbolos de entrenamiento

Simbolo_Entrenamiento= senyalRx(posicion-2*N+1:posicion);


sal = sum(Simbolo_Entrenamiento(1:N).*conj(Simbolo_Entrenamiento( N+1:end)));
f_estimado= -angle(sal)/(2*pi*N);
n= 0: length(senyalRx)-1;n=n.';
% Aplicamos la correccion
senyalRx_corr= senyalRx.*exp(-1j*2*pi*n*f_estimado);

%% Estimacion del canal con el preambulo corregido
pre_rx= senyalRx_corr(posicion-2*N+1:posicion);
pre1= pre_rx(1:N);
pre2= pre_rx(N+1:2*N);
preT_pro= (pre1+pre2)/2;
[Ypre_da,Ypre_pi]= ofdmdemod(preT_pro,N,0,0,indices_nulos,indices_pilotos);
Hest_da= Ypre_da./Xpre_da;
Hest_pi= Ypre_pi./Xpre_pi;

% Extraemos los pilotos

% figure(10)
% plot(abs(senyalRx_corr))

senyalRx_corr_2= senyalRx_corr(posicion+(1:NSimb_ofdm*(N+Ncp)));

[SimbOfdm_R, pilotosR]= ofdmdemod(senyalRx_corr_2,N,Ncp,Ncp,indices_nulos,indices_pilotos);
% Vamos a comparar los pilotos

% Ecualizacion del canal estimado, incluida la fase por el ajuste temporal
EQ = repmat(1./Hest_da,1,NSimb_ofdm);
SimbOfdm_R = SimbOfdm_R .* EQ;

% Pilotos esperados teniendo en cuenta el canal
HkPk = Pilotos .* repmat(Hest_pi,1,NSimb_ofdm);

% Una fase por simbolo OFDM: tamaño 1 x NSimb_ofdm
CFO_seg = angle(sum(pilotosR .* conj(HkPk),1));

% Repetir la correccion para las 48 portadoras de datos
exp_CFO = repmat(exp(-1j*CFO_seg),Ni,1);

SimbOfdm_R = SimbOfdm_R .* exp_CFO;


SimbOfdm_R=SimbOfdm_R(:);

Bits_R =qamdemod(SimbOfdm_R,MQAM,"bin","OutputType","bit","UnitAveragePower",true);


%% Evalaucion
[erro, Ber]= biterr(Bits_T,Bits_R)

scatterplot(SimbOfdm_R)

[audio_recuperado, Fs] = bits_conv_a_cancion(Bits_R,Ni,MQAM);
sound(audio_recuperado, Fs);

function Bits_T =cancion_a_bits_conv(archivo, inicio, duracion, Ni, MQAM)
% Codigo convolucional fijo para transportar el audio
K=3;
generadores=[7 5];
% Información del audio
info = audioinfo(archivo);
Fs = info.SampleRate;
%assert(duracion==10,'El formato fijo requiere duracion=10 segundos.');
% Seleccionar fragmento
muestraInicio = floor(inicio * Fs) + 1;
muestraFin = floor((inicio + duracion) * Fs);
muestraFin = min(muestraFin, info.TotalSamples);
if muestraInicio<=muestraFin
    [audio, Fs] = audioread(archivo, [muestraInicio muestraFin]);
else
    audio=zeros(0,info.NumChannels);
end
% Ajustar siempre a 44,1 kHz, dos canales y 441000 muestras por canal
if Fs~=44100 && ~isempty(audio)
    audio=resample(audio,44100,Fs);
end
Fs=44100;
if size(audio,2)==1
    audio=[audio audio];
else
    audio=audio(:,1:2);
end
numMuestras=441000;
numCanales=2;
audio=audio(1:min(size(audio,1),numMuestras),:);
audio=[audio;zeros(numMuestras-size(audio,1),numCanales)];
% Limitar rango
audio = max(min(audio, 1), -1);
% Pasar a PCM de 16 bits
audio_int16 = int16(audio * 32767);
% Vectorizar
muestras = audio_int16(:);
% Reinterpretar como uint16
muestras_uint16 = typecast(muestras, 'uint16');
% Pasar cada muestra a 16 bits
bits_matriz = de2bi(muestras_uint16, 16, 'left-msb');
% Convertir a vector de bits
bits_originales = bits_matriz.';
bits_originales = bits_originales(:).';
% Crear código convolucional
trellis = poly2trellis(K, generadores);
% Añadir bits de terminación
memoria = K - 1;
bits_entrada = [bits_originales zeros(1, memoria)];
% Codificar
bits_codificados = convenc(bits_entrada, trellis);
% Devolver los bits listos para el mapeo QAM y OFDM
Bits_T=double(bits_codificados(:));
bitsPorOFDM=Ni*log2(MQAM);
resto=mod(-length(Bits_T),bitsPorOFDM);
Bits_T=[Bits_T;zeros(resto,1)];
end
function [audio_recuperado, Fs] = bits_conv_a_cancion(Bits_R, Ni, MQAM)
% Formato fijo conocido por el receptor
Fs=44100;
numMuestras=441000;
numCanales=2;
K=3;
generadores=[7 5];
% Retirar el relleno OFDM antes de decodificar
numBitsCodificados=(numMuestras*numCanales*16+K-1)*numel(generadores);
bitsPorOFDM=Ni*log2(MQAM);
numBitsTrama=ceil(numBitsCodificados/bitsPorOFDM)*bitsPorOFDM;
assert(numel(Bits_R)==numBitsTrama,'Bits_R debe contener una trama completa de audio.');
bits_codificados=double(Bits_R(1:numBitsCodificados));
bits_codificados=bits_codificados(:);
%% 1. Crear el mismo código convolucional
trellis = poly2trellis(K, generadores);
%% 2. Decodificar con Viterbi
traceback = 5*(K-1);
bits_decodificados = vitdec(bits_codificados,trellis,traceback,'term','hard');
%% 3. Quitar los bits de terminación
memoria = K - 1;
bits_recuperados = bits_decodificados(1:end-memoria);
%% 4. Convertir grupos de 16 bits en muestras
bits_matriz = reshape(bits_recuperados, 16, []).';
muestras_uint16 = bi2de(bits_matriz, 'left-msb');
%% 5. Reinterpretar uint16 como int16
muestras_int16 = typecast(uint16(muestras_uint16), 'int16');
%% 6. Recuperar la forma original del audio
audio_int16 = reshape(muestras_int16, numMuestras, numCanales);
%% 7. Convertir a double
audio_recuperado = double(audio_int16) / 32767;
%% 8. Guardar canción recuperada
audiowrite('cancion_recuperada.wav', audio_recuperado, Fs);
fprintf('Cancion recuperada correctamente.\n');
fprintf('Guardada como cancion_recuperada.wav\n');
end



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