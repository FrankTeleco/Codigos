clc;clear;close all;


%% Parametros
MQAM = 4;
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
SimbQAM_T= reshape(SimbQAM_T,Ni,[]);

%% Generacion con ofdmmod
DC= N*0.5+1;
Nyquist= -32 + N*0.5+1;

indices_nulos= [-32;-31;-30;-29;-28;-27;0;27;28;29;30;31] + N*0.5 +1;
indices_pilotos=[-21;-7;7;21] + N*0.5 +1;
indice_datos= setdiff(1:N,[indices_pilotos;indices_nulos]);

%% Pilotos
Pilotos = ones(Np,NSimb_ofdm);

SimbOfdm = ofdmmod( ...
    SimbQAM_T, ...
    N, ...
    Ncp, ...
    indices_nulos, ...
    indices_pilotos, ...
    Pilotos);

%% Preambulo
[preT,Xpre_da,Xpre_pi] = ref_ofdm( ...
    N, ...
    2*Ncp, ...
    indices_nulos, ...
    indices_pilotos);

SimbOfdm=[preT;SimbOfdm];

%% Canal
retardo = 100;
CFO = 0.15/N;
hcanal = [0.2 1 -0.2 0.5 0.1];

EntradaCanal = [zeros(retardo,1); ...
                SimbOfdm; ...
                zeros(length(hcanal)-1,1)];

EntradaCanal = filter(hcanal,1,EntradaCanal);

n_canal = (0:length(EntradaCanal)-1).';

EntradaCanal = EntradaCanal .* ...
    exp(1j*2*pi*CFO*n_canal);

SenyalCanal = awgn(EntradaCanal,30,"measured");

SalidaCanal = SenyalCanal;

%% Recepcion
senyalRx = SalidaCanal;

%% Sincronizacion Temporal
filtro_inversor = conj(preT(end:-1:1));

correlacion = filter(filtro_inversor,1,senyalRx);

[valor,posicion] = max(abs(correlacion));

% Primera muestra util del primer simbolo de entrenamiento
inicioEntrenamiento = posicion - 2*N + 1;

%% Extraemos los dos simbolos de entrenamiento
Simbolo_Entrenamiento = ...
    senyalRx(inicioEntrenamiento:inicioEntrenamiento+2*N-1);

%% Estimacion CFO
sal = sum( ...
    Simbolo_Entrenamiento(1:N) .* ...
    conj(Simbolo_Entrenamiento(N+1:end)) ...
    );

f_estimado = -angle(sal)/(2*pi*N);

n = (0:length(senyalRx)-1).';

%% Correccion CFO
senyalRx_corr = ...
    senyalRx .* exp(-1j*2*pi*n*f_estimado);

%% Estimacion del canal con el preambulo corregido
pre_rx = senyalRx_corr( ...
    inicioEntrenamiento:inicioEntrenamiento+2*N-1);

pre1 = pre_rx(1:N);
pre2 = pre_rx(N+1:2*N);

preT_pro = (pre1+pre2)/2;

[Ypre_da,Ypre_pi] = ofdmdemod( ...
    preT_pro, ...
    N, ...
    0, ...
    0, ...
    indices_nulos, ...
    indices_pilotos);

Hest_da = Ypre_da./Xpre_da;
Hest_pi = Ypre_pi./Xpre_pi;

%% Extraer exactamente los simbolos OFDM de datos
inicioDatos = inicioEntrenamiento + 2*N;

numMuestrasDatos = NSimb_ofdm*(N+Ncp);

finDatos = inicioDatos + numMuestrasDatos - 1;

senyalRx_corr_2 = ...
    senyalRx_corr(inicioDatos:finDatos);

%% Demodulacion OFDM
[SimbOfdm_R,pilotosR] = ofdmdemod( ...
    senyalRx_corr_2, ...
    N, ...
    Ncp, ...
    Ncp, ...
    indices_nulos, ...
    indices_pilotos);

%% Ecualizacion
EQ = repmat(1./Hest_da,1,NSimb_ofdm);

SimbOfdm_R = SimbOfdm_R .* EQ;

%% Pilotos esperados
HkPk = Pilotos .* repmat(Hest_pi,1,NSimb_ofdm);

%% Seguimiento CPE residual
CFO_seg = angle(sum( ...
    pilotosR .* conj(HkPk), ...
    1));

exp_CFO = repmat( ...
    exp(-1j*CFO_seg), ...
    Ni, ...
    1);

SimbOfdm_R = SimbOfdm_R .* exp_CFO;

SimbOfdm_R = SimbOfdm_R(:);

%% Demodulacion QAM
Bits_R = qamdemod( ...
    SimbOfdm_R, ...
    MQAM, ...
    "bin", ...
    "OutputType","bit", ...
    "UnitAveragePower",true);

%% Evaluacion
[erro,Ber] = biterr(Bits_T,Bits_R)

scatterplot(SimbOfdm_R)

[audio_recuperado,Fs] = ...
    bits_conv_a_cancion(Bits_R,Ni,MQAM);

sound(audio_recuperado,Fs);


%% Funcion audio -> bits
function Bits_T = cancion_a_bits_conv(archivo,inicio,duracion,Ni,MQAM)

K = 3;
generadores = [7 5];

info = audioinfo(archivo);
Fs = info.SampleRate;

muestraInicio = floor(inicio*Fs)+1;
muestraFin = floor((inicio+duracion)*Fs);
muestraFin = min(muestraFin,info.TotalSamples);

if muestraInicio <= muestraFin
    [audio,Fs] = audioread(archivo,[muestraInicio muestraFin]);
else
    audio = zeros(0,info.NumChannels);
end

if Fs ~= 44100 && ~isempty(audio)
    audio = resample(audio,44100,Fs);
end

Fs = 44100;

if size(audio,2) == 1
    audio = [audio audio];
else
    audio = audio(:,1:2);
end

numMuestras = 441000;
numCanales = 2;

audio = audio(1:min(size(audio,1),numMuestras),:);

audio = [audio; ...
    zeros(numMuestras-size(audio,1),numCanales)];

audio = max(min(audio,1),-1);

audio_int16 = int16(audio*32767);

muestras = audio_int16(:);

muestras_uint16 = typecast(muestras,'uint16');

bits_matriz = de2bi(muestras_uint16,16,'left-msb');

bits_originales = bits_matriz.';
bits_originales = bits_originales(:).';

trellis = poly2trellis(K,generadores);

memoria = K-1;

bits_entrada = [bits_originales zeros(1,memoria)];

bits_codificados = convenc(bits_entrada,trellis);

Bits_T = double(bits_codificados(:));

bitsPorOFDM = Ni*log2(MQAM);

resto = mod(-length(Bits_T),bitsPorOFDM);

Bits_T = [Bits_T;zeros(resto,1)];

end


%% Funcion bits -> audio
function [audio_recuperado,Fs] = ...
    bits_conv_a_cancion(Bits_R,Ni,MQAM)

Fs = 44100;
numMuestras = 441000;
numCanales = 2;

K = 3;
generadores = [7 5];

numBitsCodificados = ...
    (numMuestras*numCanales*16+K-1)*numel(generadores);

bitsPorOFDM = Ni*log2(MQAM);

numBitsTrama = ...
    ceil(numBitsCodificados/bitsPorOFDM)*bitsPorOFDM;

assert( ...
    numel(Bits_R)==numBitsTrama, ...
    'Bits_R debe contener una trama completa de audio.');

bits_codificados = double(Bits_R(1:numBitsCodificados));
bits_codificados = bits_codificados(:);

trellis = poly2trellis(K,generadores);

traceback = 5*(K-1);

bits_decodificados = vitdec( ...
    bits_codificados, ...
    trellis, ...
    traceback, ...
    'term', ...
    'hard');

memoria = K-1;

bits_recuperados = ...
    bits_decodificados(1:end-memoria);

bits_matriz = ...
    reshape(bits_recuperados,16,[]).';

muestras_uint16 = ...
    bi2de(bits_matriz,'left-msb');

muestras_int16 = ...
    typecast(uint16(muestras_uint16),'int16');

audio_int16 = ...
    reshape(muestras_int16,numMuestras,numCanales);

audio_recuperado = ...
    double(audio_int16)/32767;

audiowrite( ...
    'cancion_recuperada.wav', ...
    audio_recuperado, ...
    Fs);

fprintf('Cancion recuperada correctamente.\n');
fprintf('Guardada como cancion_recuperada.wav\n');

end


%% Funcion preambulo
function [preT,Xpre_da,Xpre_pi] = ...
    ref_ofdm(Nfft,Ncp,I_nu,I_pi)

rng(35)

Npi = length(I_pi);
Nda = Nfft-Npi-length(I_nu);

Xpre_da = ...
    randi([0 1],Nda,1)*2-1;

Xpre_pi = ...
    randi([0 1],Npi,1)*2-1;

ref = ofdmmod( ...
    Xpre_da, ...
    Nfft, ...
    Ncp, ...
    I_nu, ...
    I_pi, ...
    Xpre_pi);

preT = [ref;ref(end-Nfft+1:end)];

end