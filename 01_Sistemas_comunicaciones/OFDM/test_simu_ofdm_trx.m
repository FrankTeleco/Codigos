%% Transmisor OFDM
clear
close all

%% Condiciones transmision

SNRdB = 300;
hcanal = 1;     %[0.2  1  -0.2  0.5  0.1];[1];
CFO_sim = 0;    %0.2/64;

% Carga datos a transmitir
load("Bits_transmision_OFDM.mat")
Nbits_dat = length(Bits_trx);


%% Ejemplo WLAN

Mqam = 4;

% Parametros OFDM
Nfft = 64;
Ncp = 16;
Nda  = 48;
Npi  = 4;
PosDC = Nfft/2+1;
Nnu = Nfft - Nda - Npi;
Nba = (Nnu - 2)/2;
I_nu = [(1:Nba+1) PosDC (Nfft-Nba+1:Nfft)]';
I_pi = [-21 -7 7 21]' +Nfft/2+1;

% Numero simbolos para transmision 
NSofdm =  ceil(Nbits_dat / (Nda * log2(Mqam)) );


%% Preambulo sincronizacion y estima canal
preT= ref_ofdm(Nfft,Ncp*2,I_nu,I_pi);


%% Genera pilotos sincro CFO
rng(20)
BitsPi = randi([0 1],Npi*1,NSofdm);
Xpi = pskmod(BitsPi,2); %zeros(Npi,NSofdm);


%% Genera simbolos QAM y OFDM
Xda = qammod(Bits_trx,Mqam,InputType="bit",UnitAveragePower=true);
Xda = reshape(Xda,Nda,NSofdm);
sofdm = ofdmmod(Xda,Nfft,Ncp,I_nu,I_pi,Xpi);


%% Monta trama
trama_ofdm = [zeros(1e4,1); preT; sofdm; zeros(2e4,1)];


%% Canal ruido + multicamino + CFO
ysal = filter(hcanal,1,trama_ofdm);
ysal = ysal .*exp(1j*2*pi*CFO_sim*(1:length(ysal))');
ysal = awgn(ysal,SNRdB,"measured");

figure
plot(abs(ysal))



%% Funcion preambulo


function [preT,Xpre_da,Xpre_pi] = ref_ofdm(Nfft,Ncp,I_nu,I_pi)
rng(35)

Npi = length(I_pi);Nda = Nfft - Npi - length(I_nu);

Xpre_da = randi([0 1],Nda,1)*2 -1;
Xpre_pi = randi([0 1],Npi,1)*2 -1;
ref = ofdmmod(Xpre_da,Nfft,Ncp,I_nu,I_pi,Xpre_pi);
preT = [ref; ref(end-Nfft+1:end)];

end
