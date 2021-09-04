clear;
close all;
clc;
format long;
%% Dados (Artigo Sharma (2008))
lambda =  3.7e9;                           %Parâmetro de lamé para sólidos elásticos
mi     =  7.9e9;                             %Parâmetro de lamé para sólidos elásticos
beta   =  4224000;                            %Coeficiente de tensão térmica
k_c      =  170;                              %Condutividade térmica
T_0     =  300;                              %Temperatura de referência
t_rel   =  10^(-5);                              %Tempo de relaxação
rho    =  2216;                              %Densidade do agregado sólido
c_e     =  1040;                             %Calor específico relacionado à tensão constante


f_min = 1; 
f_max = 1.e8;
df_1  = 1;
df_2  = 1.e2;
num_1 = linspace(f_min,10000,df_1);
num_2 = linspace(10000,f_max,df_2);
f_1   = f_min+(0:num_1)*df_1;
f_2   = f_1(end) + (1:num_2)*df_2;
f     = [f_1,f_2];
omega = 2*pi*f;

% Alocacação de memória

D_2 = zeros(1,length(omega));
D_1 = zeros(1,length(omega));
k_1 = zeros(1,length(omega));
C_pf = zeros(1,length(omega));


for j=1:length(omega)
    
D_0 = k_c*(lambda+2*mi);
D_1(j) =  -(omega(j)^2)*rho*k_c + (-t_rel*(omega(j)^2) + 1i*omega(j))*((lambda+2*mi)*rho*c_e + (beta^2)*T_0);
D_2(j) =  -(rho^2)*(omega(j)^2)*c_e*(-t_rel*(omega(j)^2)+i*omega(j));
k_1(j) =   sqrt( (-D_1(j) -sqrt(D_1(j)^2 -4*D_0*D_2(j)))/(2*D_0));
C_pf(j) = omega(j)/real(k_1(j));
    
end

%% Gráficos
g1 = semilogx(f,C_pf,'b--','LineWidth',0.5);
set(g1,'LineWidth',2);
title('Velocidade da Onda P Rapida','Interpreter','latex','Fontsize',14);
xlabel('frequencia(Hz)','Interpreter','latex','Fontsize',14);
xlim([1e2 1e10]);
ylabel('Cpf(m/s)','Interpreter','latex','fontsize',14);
legend('Arenito - 300k');


