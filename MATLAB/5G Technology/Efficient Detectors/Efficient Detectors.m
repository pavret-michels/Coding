clear all; %#ok<CLALL>
close all;
clc;
tic();%Tempo inicial

%Definindo a Potência de Transmissão
SNR_db = 0:5:20; %Potência relativa do sinal e do ruído em Decibéis
SNR=10.^(SNR_db/10); %Potência relativa do sinal e do ruído decimal
PS=1; %Potência do Sinal em Watts;

%Definindo a Transmissão
QAM = 16; %Modulação de amplitude em quadratura
n= 10000;%Dimensão do Vetor
Tx = 32; %Número de Antenas Transmissoras
Rx = Tx; %Número de Antenas Receptoras
PN = (Tx*PS)./SNR; %Potência do Ruído
Pmed = (2/3)*(QAM-1); %Potencia Média para um Sinal QAM;
Vmed =  sqrt(Pmed); %Valor médio da Tensão
Bits_Transmitido = Tx*n*log2(QAM);%Quantidade de Bits Transmitidos
k=sqrt(1/2);
fprintf('*** Múltiplas Antenas com Multiplexação Espacial ***\n-> Antenas Transmissoras:%d\n-> Antenas Receptoras:%d\n',Tx,Rx);

%Definindo o Vetor de Erros (Acumulando os Dados de Erros)
Vetor_Erros_MRC=zeros(length(SNR_db),1);
Vetor_Erros_ZF=zeros(length(SNR_db),1);
Vetor_Erros_MMSE=zeros(length(SNR_db),1);

for i=1:length(SNR_db) %Laço de repetição FOR Referente a SNR_db
    %Medida de dispersão em torno da média de uma variável aleatória.
    Desvio_Padrao = sqrt(PN(i)/2);
    %Iniciando as variáveis de Erros dos Detectores
    Erro_MRC=0;
    Erro_ZF=0;
    Erro_MMSE=0;
    
    %Retorno o horário que o Algoritmo executa o comando
    tempo = clock; 
    fprintf(['\n ',num2str(tempo(3)),'/',num2str(tempo(2)),' ',num2str(tempo(4)),'h',num2str(tempo(5)),'min',num2str(round(tempo(6))),'s'])
   
    for ii=1:n %Laço de repetição FOR referente aos "Slices" da Matriz
        %Gerando Pacote de Dados
        inf=randi([0 (QAM-1)],Tx,1);%Informação Inicial
        Pacote_Dados_Modulado = (1/Vmed)*qammod(inf,QAM);%Modulação do Pacote de Dados
        Pacote_Dados_Binario = de2bi(inf,log2(QAM));%Convertendo para Binário
    
        %Desvanecimento do Sinal Rayleigh 
        h = k*(1*randn(Rx,Tx) + 1i*randn(Rx,Tx));
        AWGN = Desvio_Padrao*(1*randn(Rx,1) + 1i*randn(Rx,1)); %Gerando o ruído
        
        %Algoritmos de Detecção
        MRC = ctranspose(h);%Detector Maximal Ration Combining
        ZF = pinv(h); %Detector Zero Forcing
        MMSE = pinv(ctranspose(h)*h + (1/SNR(i))*eye(Tx))*ctranspose(h); %Detector Mininum Mean Square Error
        %Função eye cria a Matriz Identidade(I)
        
        %Transmissão do Sinal
        Sinal_Recebido = h*Pacote_Dados_Modulado + AWGN; %Transmitindo o Sinal com o ruído
        
        %Sinal detectado 
        Sinal_Detectado_MRC = MRC*Sinal_Recebido;
        Sinal_Detectado_ZF = ZF*Sinal_Recebido;
        Sinal_Detectado_MMSE = MMSE*Sinal_Recebido;
                
        %Demodulação do Sinal
        Pacote_Dados_Recebido_MRC = qamdemod(Vmed*Sinal_Detectado_MRC,QAM);
        Pacote_Dados_Recebido_ZF = qamdemod(Vmed*Sinal_Detectado_ZF,QAM);
        Pacote_Dados_Recebido_MMSE = qamdemod(Vmed*Sinal_Detectado_MMSE,QAM);
        
        %Convertendo o Sinal recebido para Binário
        Pacote_Dados_Recebido_Binario_MRC = de2bi(Pacote_Dados_Recebido_MRC,log2(QAM));
        Pacote_Dados_Recebido_Binario_ZF = de2bi(Pacote_Dados_Recebido_ZF,log2(QAM));
        Pacote_Dados_Recebido_Binario_MMSE = de2bi(Pacote_Dados_Recebido_MMSE,log2(QAM));
     
        %Comparação bit a bit dos Dados Enviado com os Dados Recebido
        Erro_MRC = Erro_MRC+symerr(Pacote_Dados_Binario,Pacote_Dados_Recebido_Binario_MRC);
        Erro_ZF = Erro_ZF+symerr(Pacote_Dados_Binario,Pacote_Dados_Recebido_Binario_ZF);
        Erro_MMSE = Erro_MMSE+symerr(Pacote_Dados_Binario,Pacote_Dados_Recebido_Binario_MMSE);
    end
    %Retorna o valor de SNR_db 
    fprintf('\nSimulação com SNR_db: %ddb \n',SNR_db(i));
        
    %Acumulando o Total de Erros no Vetor
    Vetor_Erros_MRC(i) = Erro_MRC;
    Vetor_Erros_ZF(i) = Erro_ZF;
    Vetor_Erros_MMSE(i) = Erro_MMSE;
end
%Função Bit Error Rate (Taxa de Erro Por Bit)
BERMRC = Vetor_Erros_MRC/(Bits_Transmitido);
BERZF = Vetor_Erros_ZF/(Bits_Transmitido);
BERMMSE = Vetor_Erros_MMSE/(Bits_Transmitido);

%Plota o Gráfico de Função em monolog sendo o Eixo y a Variável SNR_db
semilogy(SNR_db,BERMRC,'r',SNR_db,BERZF,'b',SNR_db,BERMMSE,'k')
legend('MRC','ZF','MMSE');
title('Multiplexação Espacial'); %Define o Título de Gráfico
xlabel('SNR (db)')
ylabel('Taxa de Erro Por Bit')

%Retorno o tempo total gasto na execução do Algoritmo
tf=toc();
fprintf('\n Tempo Total de Simulação: %.3f\n',tf)