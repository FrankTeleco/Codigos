clc;clear;close all;

%% RECEPTOR

datos = load("Captura_RTL_ofdm.mat");
sr = datos.sig_rcx;
sr = sr(:);

load("Bits_transmision_OFDM.mat")
Nbits_dat = length(Bits_trx);

%% Parametros WLAN

Mqam = 4;

%% Parametros OFDM

Nfft = 64;
Ncp = 16;
Nda = 48;
Npi = 4;

PosDC = Nfft/2+1;
Nnu = Nfft-Nda-Npi;
Nba = (Nnu-2)/2;

I_nu = [(1:Nba+1) PosDC (Nfft-Nba+1:Nfft)]';
I_pi = [-21 -7 7 21]'+Nfft/2+1;

NSofdm = ceil(Nbits_dat/(Nda*log2(Mqam)));

%% Visualizacion señal capturada

figure
plot(abs(sr))
grid on
title("Señal capturada RTL-SDR")
xlabel("Muestra")
ylabel("|sr|")

%% Preambulo conocido

[preT,Xpre_da,Xpre_pi] = ref_ofdm(Nfft,Ncp*2,I_nu,I_pi);

%% Sincronizacion temporal

preCOR = conj(preT(end:-1:1));

salco = filter(preCOR,1,sr);

[ma,po] = max(abs(salco));

poIni = po-4;

figure
plot(abs(salco))
grid on
hold on
plot(po,ma,"o")
title("Correlacion con el preambulo")
xlabel("Muestra")
ylabel("|Correlacion|")

%% Estimacion CFO

poIni_pre = poIni-Ncp*2-Nfft*2;

pini1 = poIni_pre+Ncp*2;
pini2 = pini1+Nfft;

z = sum(sr(pini1+(1:Nfft)).*conj(sr(pini2+(1:Nfft))));

CFO_est = -angle(z)/(2*pi*Nfft);

disp("CFO estimado:")
disp(CFO_est)

%% Correccion CFO

sr_sin_cfo = sr.*exp(-1j*2*pi*CFO_est*(1:length(sr))');

%% Estimacion canal

preT_pro = 0.5*(sr_sin_cfo(pini1+(1:Nfft))+sr_sin_cfo(pini2+(1:Nfft)));

[Ypre_da,Ypre_pi] = ofdmdemod(preT_pro,Nfft,0,0,I_nu,I_pi);

Hest_da = Ypre_da./Xpre_da;
Hest_pi = Ypre_pi./Xpre_pi;

figure
Hpre = fftshift(fft(preT_pro,Nfft));
plot(-32:31,abs(Hpre))
grid on
title("Canal estimado")
xlabel("Subportadora")
ylabel("|H|")

%% Ecualizador

EQ_est = repmat(1./Hest_da,1,NSofdm);

%% Genera los mismos pilotos usados en TX

rng(20)

BitsPi = randi([0 1],Npi,NSofdm);

Xpi = pskmod(BitsPi,2);

HkPk = Xpi.*repmat(Hest_pi,1,NSofdm);

%% Extraccion simbolos OFDM

Lsofdm = NSofdm*(Ncp+Nfft);

sr_ofdm = sr_sin_cfo(poIni+(1:Lsofdm));

[Yda,Ypi] = ofdmdemod(sr_ofdm,Nfft,Ncp,Ncp,I_nu,I_pi);

%% Compensacion canal

Yda_eq = Yda.*EQ_est;

scatterplot(Yda_eq(:))
title("Despues de ecualizar canal")

%% CFO residual con pilotos

CFO_seg = angle(sum(Ypi.*conj(HkPk),1));

exp_CFO = repmat(exp(-1j*CFO_seg),Nda,1);

Yda_eq_cfo = Yda_eq.*exp_CFO;

figure
plot(CFO_seg)
grid on
title("Fase residual estimada")
xlabel("Simbolo OFDM")
ylabel("Fase [rad]")

scatterplot(Yda_eq_cfo(:))
title("Despues de corregir CFO residual")

%% Bits recibidos

Bits_rcx = qamdemod(Yda_eq_cfo,Mqam,OutputType="bit",UnitAveragePower=true);

Bits_rcx = Bits_rcx(:);

if length(Bits_rcx) > Nbits_dat
    Bits_rcx = Bits_rcx(1:Nbits_dat);
end

[ne,be] = biterr(Bits_trx(:),Bits_rcx);

disp("Numero de errores:")
disp(ne)

disp("BER:")
disp(be)

%% LLR

llr_rcx = qamdemod(Yda_eq_cfo,Mqam,OutputType="llr",UnitAveragePower=true);

llr_rcx = llr_rcx(:);

if length(llr_rcx) > Nbits_dat
    llr_rcx = llr_rcx(1:Nbits_dat);
end

%% Reconstruccion imagen

reconstruye_jpg(llr_rcx);

disp("Imagen reconstruida: imagen_rcx.jpg")

%% FUNCIONES

function [preT,Xpre_da,Xpre_pi] = ref_ofdm(Nfft,Ncp,I_nu,I_pi)

rng(35)

Npi = length(I_pi);
Nda = Nfft-Npi-length(I_nu);

Xpre_da = randi([0 1],Nda,1)*2-1;
Xpre_pi = randi([0 1],Npi,1)*2-1;

ref = ofdmmod(Xpre_da,Nfft,Ncp,I_nu,I_pi,Xpre_pi);

preT = [ref;ref(end-Nfft+1:end)];

end