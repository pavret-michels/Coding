clear all; %#ok<CLALL>
close all; clc; tic();%Tempo inicial

%Definindo a Transmiss�o
fprintf('*** Efici�ncia Espectral do Sinal ***');
M = [50;100;150;200;250;300;350;400;450;500]; %Antena Transmissora
K = [5;10;15;20;25;30;35;40;45;50]; %Usu�rios
k=sqrt(1/2);
Monte_Carlo= 10;%N�mero de repeti��es

%Definindo o Vetor do C�lculo de Efici�ncia Energ�tica
%Perfect Channel State Information
Energy_Efficiency_MRC=zeros(length(K),1);
Energy_Efficiency_ZF=zeros(length(K),1);
Energy_Efficiency_MMSE=zeros(length(K),1);

%Definindo o Vetor do C�lculo de Efici�ncia Espectral
%Perfect Channel State Information
Spectral_Efficiency_MRC=zeros(length(K),1);
Spectral_Efficiency_ZF=zeros(length(K),1);
Spectral_Efficiency_MMSE=zeros(length(K),1);
PN = 1; %Pot�ncia do Ru�do

%Medida de dispers�o em torno da m�dia de uma vari�vel aleat�ria.    
Desvio_Padrao = sqrt(PN/2);

v_SNR_db=-20:1:50;
for I_SNR=1:length(v_SNR_db)
    %Definindo a Pot�ncia de Transmiss�o
    SNR_db = v_SNR_db(I_SNR); %Pot�ncia relativa do sinal e do ru�do em Decib�is
    SNR=10^(SNR_db/10); %Pot�ncia relativa do sinal e do ru�do em decimal
    %Pot�ncia M�dia do Usu�rio
    P_u = SNR./M; %DownLink
    PN = 1; %Pot�ncia do Ru�do

for A=1:length(K)%La�o referente ao n�mero de Antenas/Usu�rios 
        %Vari�veis auxiliares do La�o Monte-Carlo
        Spectral_Efficiency_MRC_Parcial=0;
        Spectral_Efficiency_ZF_Parcial=0;
        Spectral_Efficiency_MMSE_Parcial=0;
        Energy_Efficiency_MRC_Parcial=0;
        Energy_Efficiency_ZF_Parcial=0;
        Energy_Efficiency_MMSE_Parcial=0;
        
    for B=1:Monte_Carlo %La�o Monte-Carlo
        %Desvanecimento do Sinal Rayleigh 
        h = k*(1*randn(M(A),K(A)) + 1i*randn(M(A),K(A)));
        AWGN = Desvio_Padrao*(1*randn(M(A),1) + 1i*randn(M(A),1)); %Gerando o ru�do
        %Algoritmos de Detec��o 
        %Perfect Channel State Information
        MRC = h;%Detector Maximal Ration Combining
        ZF = h*(pinv(ctranspose(h)*h)); %Detector Zero Forcing
        MMSE = h*(pinv((ctranspose(h)*h+(1/SNR)*eye(K(A))))); %Detector Mininum Mean Square Error
                              
        for Kth = 1:K(A)
            %C�lculo do Numerador
            Numerator_MRC = P_u(A)*abs(ctranspose(MRC(:,Kth))*h(:,Kth))^2;
            Numerator_ZF = P_u(A)*abs(ctranspose(ZF(:,Kth))*h(:,Kth))^2;
            Numerator_MMSE = P_u(A)*abs(ctranspose(MMSE(:,Kth))*h(:,Kth))^2;
                        
            %M�todo por Itera��o
            Denominator_MRC = 0;
            Denominator_ZF = 0;
            Denominator_MMSE = 0;
                for MI = 1:K(A)
                      if MI ~= Kth
                       Denominator_MRC = Denominator_MRC + abs(ctranspose(MRC(:,Kth))*h(:,MI))^2;
                       Denominator_ZF = Denominator_ZF + abs(ctranspose(ZF(:,Kth))*h(:,MI))^2;
                       Denominator_MMSE = Denominator_MMSE + abs(ctranspose(MMSE(:,Kth))*h(:,MI))^2;
                      end
                end
                    %C�lculo da Efici�ncia Espectral Parcial
                    Spectral_Efficiency_MRC_Parcial = Spectral_Efficiency_MRC_Parcial + log2(1 + (Numerator_MRC/(P_u(A)*Denominator_MRC + norm(MRC(:,Kth))^2)));
                    Spectral_Efficiency_ZF_Parcial = Spectral_Efficiency_ZF_Parcial + log2(1 + (Numerator_ZF/(P_u(A)*Denominator_ZF + norm(ZF(:,Kth))^2)));
                    Spectral_Efficiency_MMSE_Parcial = Spectral_Efficiency_MMSE_Parcial + log2(1 + (Numerator_MMSE/(P_u(A)*Denominator_MMSE + norm(MMSE(:,Kth))^2)));
        end
    end
        Spectral_Efficiency_MRC(A) = Spectral_Efficiency_MRC_Parcial/Monte_Carlo;
        Spectral_Efficiency_ZF(A) = Spectral_Efficiency_ZF_Parcial/Monte_Carlo;
        Spectral_Efficiency_MMSE(A) = Spectral_Efficiency_MMSE_Parcial/Monte_Carlo;
        
        %Energy-Efficiency
        Energy_Efficiency_MRC(A) = (1/sum(P_u(A))).*Spectral_Efficiency_MRC(A);
        Energy_Efficiency_ZF(A) = (1/sum(P_u(A))).*Spectral_Efficiency_ZF(A);
        Energy_Efficiency_MMSE(A) = (1/sum(P_u(A))).*Spectral_Efficiency_MMSE(A);
    
end
    V_Spectral_Efficiency_MRC(I_SNR)=mean(Spectral_Efficiency_MRC);
    V_Energy_Efficiency_MRC(I_SNR)=mean(Energy_Efficiency_MRC);
    V_Spectral_Efficiency_MMSE(I_SNR)=mean(Spectral_Efficiency_MMSE);
    V_Energy_Efficiency_MMSE(I_SNR)=mean(Energy_Efficiency_MMSE);
    V_Spectral_Efficiency_ZF(I_SNR)=mean(Spectral_Efficiency_ZF);
    V_Energy_Efficiency_ZF(I_SNR)=mean(Energy_Efficiency_ZF);
end
 
%Gr�fico Spectral Efficiency 
%Perfect Channel State Information
figure(1)
handle=plot(M,Spectral_Efficiency_MRC,'r-x',M,Spectral_Efficiency_ZF,'b-x',M,Spectral_Efficiency_MMSE,'k--x');
set(handle,'LineWidth',1.5);
title(legend,'Algoritmo de Detec��o (A)');
legend('MRC','ZF','MMSE');
grid on
title('Spectral Efficiency - Perfect CSI'); %Define o T�tulo de Gr�fico
xlabel('Number of Base Station Antennas (M)')
%xlabel('Number of Users (K)')
ylabel('Spectral-Efficiency (bits/s/Hz)')

%Gr�fico Energy-Efficiency
figure(2)
handle=semilogy(Spectral_Efficiency_MRC,Energy_Efficiency_MRC,'r-x',Spectral_Efficiency_ZF,Energy_Efficiency_ZF,'b-x',Spectral_Efficiency_MMSE,Energy_Efficiency_MMSE,'k--x');
set(handle,'LineWidth',1.5);
title(legend,'Algoritmo de Detec��o (A)');
legend('MRC','ZF','MMSE');
grid on
title('Energy Efficiency - Perfect CSI'); %Define o T�tulo de Gr�fico
xlabel('Spectral-Efficiency (bits/s/Hz)')
ylabel('Efici�ncia Energ�tica (bits/s/J)')

figure(3)
handle=plot(v_SNR_db,V_Energy_Efficiency_MMSE,'k-',v_SNR_db,V_Spectral_Efficiency_MMSE,'k-x',v_SNR_db,V_Energy_Efficiency_MRC,'r-',v_SNR_db,V_Spectral_Efficiency_MRC,'r-x',v_SNR_db,V_Energy_Efficiency_ZF,'b--',v_SNR_db,V_Spectral_Efficiency_ZF,'b--x');
set(handle,'LineWidth',1.5);
title(legend,'Algoritmo de Detec��o (A)');
legend('MMSE_{Energy Efficiency}','MMSE_{Spectral Efficiency}','MRC_{Energy Efficiency}','MRC_{Spectral Efficiency}','ZF_{Energy Efficiency}','ZF_{Spectral Efficiency}');
xlabel('Rela��o Sinal-Ru�do (db)')
ylabel('Energy Efficiency x Spectral-Efficiency')

tf=toc();
fprintf('\n Tempo Total de Simula��o: %.2f(min)\n ',tf/60)