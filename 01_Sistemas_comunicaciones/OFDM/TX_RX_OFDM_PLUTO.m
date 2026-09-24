clc;clear;close all;


%% Parametros
MQAM = 16;
Nbits = 1e6; % Total de bits de informacion a transmitir
NbitsTrama = 100000; % Bits por trama para no exceder la memoria del Pluto
N= 64;
Ncp= N/4;
Ni= 48; % portadoras datos
Np= 4; % Portadoras pilotos
Nu= 12; % Portadoras nulas

%% setup pluto
fc = 1e9;
FsPluto = 1e6;
tx = sdrtx('Pluto','RadioID','usb:0','CenterFrequency',fc,'BasebandSampleRate',FsPluto,'Gain',-30);
rx = sdrrx('Pluto','RadioID','usb:0','CenterFrequency',fc,'BasebandSampleRate',FsPluto,'GainSource','Manual','Gain',10,'OutputDataType','double','EnableBurstMode',true,'NumFramesInBurst',1);

%% Indices y preambulo
% Calculamos los indices
DC= N*0.5+1;% esto sale en 33
Nyquist= -32 + N*0.5+1; % esto sale a 1
indices_nulos= [-32;-31;-30;-29;-28;-27;0; 27;28;29;30;31] + N*0.5 +1; % Bandas de guarda
indices_pilotos=[-21;-7;7;21] + N*0.5 +1;
indice_datos= setdiff(1:N,[indices_pilotos;indices_nulos]);
[preT,Xpre_da,Xpre_pi] = ref_ofdm(N,2*Ncp,indices_nulos,indices_pilotos);
filtro_inversor = conj(preT(end:-1:1));
bitsPorOFDM = Ni*log2(MQAM);
bitsEvaluados = 0;
erroTotal = 0;
numTramas = ceil(Nbits/NbitsTrama);

% Para interrumpir: Ctrl+C y despues release(tx); release(rx);
for iTrama = 1:numTramas
    %% Generacion de bits
    bitsEstaTrama = min(NbitsTrama,Nbits-bitsEvaluados);
    Bits_T = randi([0 1],bitsEstaTrama,1);
    NSimb_ofdm = ceil(bitsEstaTrama/bitsPorOFDM);
    numRelleno = NSimb_ofdm*bitsPorOFDM-bitsEstaTrama;
    Bits_trama = [Bits_T; zeros(numRelleno,1)];

    %% Mapeo y generacion OFDM
    SimbQAM_T = qammod(Bits_trama,MQAM,"bin","InputType","bit","UnitAveragePower",true);
    SimbQAM_T = reshape(SimbQAM_T,Ni,[]);
    Pilotos = ones(Np,NSimb_ofdm);
    SimbOfdm = ofdmmod(SimbQAM_T,N,Ncp,indices_nulos,indices_pilotos,Pilotos);
    SimbOfdm = [preT; SimbOfdm];

    %% Canal: Pluto real, sin canal simulado
    release(rx);
    release(tx);
    rx.SamplesPerFrame = 2*numel(SimbOfdm)+2*numel(preT);
    EntradaCanal = 0.8*SimbOfdm/max(abs(SimbOfdm));
    transmitRepeat(tx,EntradaCanal);

    %% Recepcion
    valido = false;
    while ~valido
        [senyalRx,valido,overflow] = rx();
    end
    % Una rafaga independiente: no concatenar capturas ni descartar el overflow inicial
    release(tx);
    release(rx);

    %% Sincronizacion temporal
    correlacion = filter(filtro_inversor,1,senyalRx);
    % Evitar seleccionar la ultima repeticion sin todos sus datos detras
    inicioBusqueda = numel(preT)+Ncp/2;
    finBusqueda = numel(senyalRx)-NSimb_ofdm*(N+Ncp)+Ncp/2;
    [valor,indice] = max(abs(correlacion(inicioBusqueda:finBusqueda)));
    posicion = inicioBusqueda+indice-1-Ncp/2;

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


    %% Evaluacion de los bits utiles de esta trama
    Bits_R = Bits_R(1:bitsEstaTrama);
    [erro,BerTrama] = biterr(Bits_T,Bits_R);
    erroTotal = erroTotal+erro;
    bitsEvaluados = bitsEvaluados+bitsEstaTrama;
    Ber = erroTotal/bitsEvaluados;
    fprintf('Trama %d/%d | Bits: %d | BER trama: %.6g | BER acumulada: %.6g\n',iTrama,numTramas,bitsEvaluados,BerTrama,Ber);
end

release(tx);
release(rx);
scatterplot(SimbOfdm_R);
fprintf('Bits evaluados: %d | Errores: %d | BER: %.6g\n',bitsEvaluados,erroTotal,Ber);

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