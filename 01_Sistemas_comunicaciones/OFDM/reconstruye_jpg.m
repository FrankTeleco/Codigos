function  reconstruye_jpg(llr_rcx)
% 
% Reconstruye la imagen a partir de los bits recibidos
%
% Como entrada requiere los LLR a partir de los símbolos recibidos 
% (ecualizados y corregidos de CFO)
%

    Nbits_dat = 173664;
    
    % llr_rcx= qamdemod(Yr_eq_cfo,Mqam,OutputType="llr",UnitAveragePower=true);

    llr_rcx = llr_rcx(:);
    Lllr = length(llr_rcx);

    if Lllr ~= Nbits_dat
        error("El numero de bits recibidos es erroneo")
    end

    load("matriz_LDPC_648_R12.mat")
    [NumFil,NumCol] = size(Hs);
    NumBitsLDPC = NumCol-NumFil;
    Rate = NumBitsLDPC/NumCol;
    
    codcfg = ldpcEncoderConfig(Hs);
    deccfg = ldpcDecoderConfig(codcfg,"layered-bp");
    Niter = 30; %Numero de iteraciones el decodificador

    NumBlo = Nbits_dat/NumCol;

    %LLRrec = -(Bits_rcx*2-1);
    LLRrec = reshape(llr_rcx(:),NumCol,NumBlo);
    [BitsDec,Nitusadas] = ldpcDecode(LLRrec,deccfg,Niter);

    BitsDec = BitsDec(:);

    % Guarda jpg
    da_bi_da = bit2int(BitsDec(:),8);
    fe = fopen(['imagen_rcx.jpg'],'w');
    fwrite(fe,da_bi_da,"uint8");
    fclose(fe);

    
end


