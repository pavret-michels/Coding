clear all; %#ok<CLALL>
close all;
clc;
tic();%Tempo inicial

%Definindo a Pot�ncia de Transmiss�o
SNR_db = 0:5:20; %Pot�ncia relativa do sinal e do ru�do em Decib�is
SNR=10.^(SNR_db/10); %Pot�ncia relativa do sinal e do ru�do decimal
PS=1; %Pot�ncia do Sinal em Watts;

%Definindo a Transmiss�o
QAM = 16; %Modula��o de amplitude em quadratura
n= 10000;%Dimens�o do Vetor
Tx = 32; %N�mero de Antenas Transmissoras
Rx = Tx; %N�mero de Antenas Receptoras
PN = (Tx*PS)./SNR; %Pot�ncia do Ru�do
Pmed = (2/3)*(QAM-1); %Potencia M�dia para um Sinal QAM;
Vmed =  sqrt(Pmed); %Valor m�dio da Tens�o
Bits_Transmitido = Tx*n*log2(QAM);%Quantidade de Bits Transmitidos
k=sqrt(1/2);
fprintf('*** M�ltiplas Antenas com Multiplexa��o Espacial ***\n-> Antenas Transmissoras:%d\n-> Antenas Receptoras:%d\n',Tx,Rx);

%Definindo o Vetor de Erros (Acumulando os Dados de Erros)
Vetor_Erros_MRC=zeros(length(SNR_db),1);
Vetor_Erros_ZF=zeros(length(SNR_db),1);
Vetor_Erros_MMSE=zeros(length(SNR_db),1);

for i=1:length(SNR_db) %La�o de repeti��o FOR Referente a SNR_db
    %Medida de dispers�o em torno da m�dia de uma vari�vel aleat�ria.
    Desvio_Padrao = sqrt(PN(i)/2);
    %Iniciando as vari�veis de Erros dos Detectores
    Erro_MRC=0;
    Erro_ZF=0;
    Erro_MMSE=0;
    
    %Retorno o hor�rio que o Algoritmo executa o comando
    tempo = clock; 
    fprintf(['\n ',num2str(tempo(3)),'/',num2str(tempo(2)),' ',num2str(tempo(4)),'h',num2str(tempo(5)),'min',num2str(round(tempo(6))),'s'])
   
    for ii=1:n %La�o de repeti��o FOR referente aos "Slices" da Matriz
        %Gerando Pacote de Dados
        inf=randi([0 (QAM-1)],Tx,1);%Informa��o Inicial
        Pacote_Dados_Modulado = (1/Vmed)*qammod(inf,QAM);%Modula��o do Pacote de Dados
        Pacote_Dados_Binario = de2bi(inf,log2(QAM));%Convertendo para Bin�rio
    
        %Desvanecimento do Sinal Rayleigh 
        h = k*(1*randn(Rx,Tx) + 1i*randn(Rx,Tx));
        AWGN = Desvio_Padrao*(1*randn(Rx,1) + 1i*randn(Rx,1)); %Gerando o ru�do
        
        %Algoritmos de Detec��o
        MRC = ctranspose(h);%Detector Maximal Ration Combining
        ZF = pinv(h); %Detector Zero Forcing
        MMSE = pinv(ctranspose(h)*h + (1/SNR(i))*eye(Tx))*ctranspose(h); %Detector Mininum Mean Square Error
        %Fun��o eye cria a Matriz Identidade(I)
        
        %Transmiss�o do Sinal
        Sinal_Recebido = h*Pacote_Dados_Modulado + AWGN; %Transmitindo o Sinal com o ru�do
        
        %Sinal detectado 
        Sinal_Detectado_MRC = MRC*Sinal_Recebido;
        Sinal_Detectado_ZF = ZF*Sinal_Recebido;
        Sinal_Detectado_MMSE = MMSE*Sinal_Recebido;
                
        %Demodula��o do Sinal
        Pacote_Dados_Recebido_MRC = qamdemod(Vmed*Sinal_Detectado_MRC,QAM);
        Pacote_Dados_Recebido_ZF = qamdemod(Vmed*Sinal_Detectado_ZF,QAM);
        Pacote_Dados_Recebido_MMSE = qamdemod(Vmed*Sinal_Detectado_MMSE,QAM);
        
        %Convertendo o Sinal recebido para Bin�rio
        Pacote_Dados_Recebido_Binario_MRC = de2bi(Pacote_Dados_Recebido_MRC,log2(QAM));
        Pacote_Dados_Recebido_Binario_ZF = de2bi(Pacote_Dados_Recebido_ZF,log2(QAM));
        Pacote_Dados_Recebido_Binario_MMSE = de2bi(Pacote_Dados_Recebido_MMSE,log2(QAM));
     
        %Compara��o bit a bit dos Dados Enviado com os Dados Recebido
        Erro_MRC = Erro_MRC+symerr(Pacote_Dados_Binario,Pacote_Dados_Recebido_Binario_MRC);
        Erro_ZF = Erro_ZF+symerr(Pacote_Dados_Binario,Pacote_Dados_Recebido_Binario_ZF);
        Erro_MMSE = Erro_MMSE+symerr(Pacote_Dados_Binario,Pacote_Dados_Recebido_Binario_MMSE);
    end
    %Retorna o valor de SNR_db 
    fprintf('\nSimula��o com SNR_db: %ddb \n',SNR_db(i));
        
    %Acumulando o Total de Erros no Vetor
    Vetor_Erros_MRC(i) = Erro_MRC;
    Vetor_Erros_ZF(i) = Erro_ZF;
    Vetor_Erros_MMSE(i) = Erro_MMSE;
end
%Fun��o Bit Error Rate (Taxa de Erro Por Bit)
BERMRC = Vetor_Erros_MRC/(Bits_Transmitido);
BERZF = Vetor_Erros_ZF/(Bits_Transmitido);
BERMMSE = Vetor_Erros_MMSE/(Bits_Transmitido);

%Plota o Gr�fico de Fun��o em monolog sendo o Eixo y a Vari�vel SNR_db
semilogy(SNR_db,BERMRC,'r',SNR_db,BERZF,'b',SNR_db,BERMMSE,'k')
legend('MRC','ZF','MMSE');
title('Multiplexa��o Espacial'); %Define o T�tulo de Gr�fico
xlabel('SNR (db)')
ylabel('Taxa de Erro Por Bit')

%Retorno o tempo total gasto na execu��o do Algoritmo
tf=toc();
fprintf('\n Tempo Total de Simula��o: %.3f\n',tf)