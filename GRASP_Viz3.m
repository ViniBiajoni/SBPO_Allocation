% Find critical ck-tuple through search algorithms
clear all
original_path=pwd;
format long
current_path=pwd; % get current path
%addpath(strcat(current_path,'\LinearSE'),path); % add Linear SE folder
addpath(strcat(current_path,'\medidas'),path); % add medidas folder
addpath(strcat(current_path,'\sistemas'),path); % add sistemas folder
cd(original_path);
clc

 global barras_atuais
 global UM
 global lote
 global H
 global nm
 global w_k
 global M
 global start_test
 global time_max
%##########################################################################
%% 03 bus TEST SYSTEM
%med_file='3Barras_tese.med';
%med_file='3Barras_tese(2).med';
%sis_file='3Barras_tese.cdf';
%##########################################################################
%% 06 bus TEST SYSTEM
%med_file='6BusCase01_TPS2015.med';
%med_file='6BusCase02_TPS2015.med';
%med_file='6BusCase05_TPS2015.med';
%med_file='6BusCase06_TPS2015.med';
%med_file='6BusCase07_TPS2015.med';
%med_file='6BusCase08_TPS2015.med';
%med_file='6Bus_v2Case03_TPS2015.med';
%med_file='6Bus_v2Case04_TPS2015.med';
%med_file='6BusCaseUM1.med';
%med_file='6BusCaseUM2.med';
%med_file='6BusCaseUM3.med';
%med_file='6BusCaseUM4.med';
%med_file='6BusCaseUM5.med';
%med_file='6BusCaseUM6.med';
%med_file='6Busv2CaseUM5_UM6.med';

%sis_file='6Bus.cdf';
%sis_file='6Busv2.cdf';
%##########################################################################
%% IEEE 06 bus TEST SYSTEM
%med_file='ieee6Case01.med';
%sis_file='ieee6.cdf';
%##########################################################################
%% IEEE 14 bus TEST SYSTEM
%med_file='ieee14_london.med';
%med_file='ieee14_london_teste.med';
%med_file='ieee14exampleTPS2017.med';

%med_file='ViniIEEE14BusFull_Meas.med';
%########14 red#######
% med_file='IEEE14Bus_ViniRedundancy0.4.med';
% med_file='IEEE14Bus_ViniRedundancy0.5.med';
% med_file='IEEE14Bus_ViniRedundancy0.6.med';
% med_file='IEEE14Bus_ViniRedundancy0.7.med';
% med_file='IEEE14Bus_ViniRedundancy0.8.med';
% med_file='IEEE14Bus_ViniRedundancy0.9.med';
%#####################
%sis_file='ieee14.cdf';
%##########################################################################
%% IEEE 24 bus TEST SYSTEM
%sis_file='ieee24.cdf';
%########24 red#######
%med_file='IEEE24Bus_Complete_ViniRedundancy1.5_1.med';
%med_file='IEEE24Bus_Complete_ViniRedundancy2_1.med';
%med_file='IEEE24Bus_Complete_ViniRedundancy2.5_1.med';
%#####################
%med_file='ViniIEEE24BusFull_Meas.med';

%med_file='ieee24Case01_SBSE2014.med';%Caso Base OK
%med_file='24BusUM10eUM15.med';%Caso base com UM10 e UM15 OK
%med_file='ieee24Case02(4)CBA2016.med';
%med_file='IEEE24Case02_TPS2015.med';
%med_file='ieee24Case03_TPS2015.med';
%med_file='ieee24Case04_TPS2015.med';
%med_file='ieee24Case05_TPS2015.med';
%med_file='ieee24Case06_TPS2015.med';
%med_file='ieee24Case07_TPS2015.med';
%med_file='ieee24Case07(1)_TPS2015.med';
%med_file='ieee24Case07(2)_TPS2015.med';
%med_file='ieee24Case07(3)_TPS2015.med';
%med_file='ieee24Case08_TPS2015.med';
%med_file='IEEE24Case08(2)_TPS2015.med';
%med_file='IEEE24Case08(3)_TPS2015.med';
%med_file='ieee24Case10_TPS2015.med';
%med_file='ieee24Case11_TPS2015.med';
%med_file='ieee24Case12_TPS2015.med';
%med_file='ieee24Case13_TPS2015.med';
%med_file='ieee24Case14_TPS2015.med';
%med_file='ieee24_gen_meet_2013_Obs_CASE01.med';
%med_file='ieee24_gen_meet_2013_Obs_CASE01.med';
%med_file='ieee24_gen_meet_2013_case02.med';
%med_file='ieee24_gen_meet_2013_case03.med';
%med_file='ieee24_gen_meet_2013_case04.med';
%med_file='ieee24_gen_meet_2013_case05.med';
%##########################################################################
%% IEEE 30 bus TEST SYSTEM
%med_file='ieee30_observability_PSCC_2014.med';%Caso Base 37 meas 20 MU
%med_file='S1_CS10_81M.med';
%med_file='S1_CS10_81M_62M_PMU26491012.med';
%med_file='ieee30_observability_full_2014.med';
%econk1= dlmread('IEEE30PSCCk1.txt');
%med_file='ieee30_observability_PSCC_2014B.med'; % Caso 42 meas
med_file='ieee30_observability_PSCC_2014B_43.med'; % Caso 43 meas

%med_file='ViniIEEE30BusFull_Meas.med';
%########30 red#######
%med_file='IEEE30Bus_Complete_ViniRedundancy1.5_1.med';
%med_file='IEEE30Bus_Complete_ViniRedundancy2_1.med';
%med_file='IEEE30Bus_Complete_ViniRedundancy2_1.5.med';
%#####################

sis_file='ieee30.cdf';
%##########################################################################
%% IEEE 118 bus TEST SYSTEM
%med_file='S1_CS5_352M.med'; %Plan1
%med_file='S1_CS5_352M_90M_PMU.med'; %Plan2
%econk1= dlmread('IEEE118k1.txt');
%med_file='118_teste.med';
%########118 red#######
% med_file='IEEE118Bus_Complete_ViniRedundancy1.5_1.med';
% med_file='IEEE118Bus_Complete_ViniRedundancy2_1.med';
% med_file='IEEE118Bus_Complete_ViniRedundancy2.5_!.med';
%########118 red#######
% med_file='IEEE118Bus_ViniRedundancy0.4.med';
% med_file='IEEE118Bus_ViniRedundancy0.5.med';
% med_file='IEEE118Bus_ViniRedundancy0.6.med';
% med_file='IEEE118Bus_ViniRedundancy0.7.med';
% med_file='IEEE118Bus_ViniRedundancy0.8.med';
% med_file='IEEE118Bus_ViniRedundancy0.9.med';
%#####################
%sis_file='ieee118.cdf';
%##########################################################################
% Convençoes:
% - Medidas de fluxo de potencia ativa:     tipo 1
% - Medidas de injecao de potencia ativa:   tipo 2
% - Medidas de angulo:                      tipo 3
% - Medidas de fluxo de potencia reativa:   tipo 4
% - Medidas de injecao de potencia reativa: tipo 5
% - Medidas de tensao:                      tipo 6
% - Medidas de corrente de ramo real        tipo 7
% - Medidas de corrente de ramo imag        tipo 8
% - Medidas de injecao de corrente real     tipo 9
% - Medidas de injecao de corrente imag     tipo 10
% - Medidas de fluxo de 1 => 2   ---> de=1 e para=2
% - Medidas de tensao na barra 3 ---> de=0 e para=3
%##########################################################################

[Caso, Barra, Ramo]=ler_sistema(sis_file); % read network data
[medidores]=ler_medidores(med_file); % read measurement data
%generate_full_linear_case(Caso,Ramo);
%######################Gen Diff Redundancies Sets##########################
% This part is now ina script named Cases_Generator

%gen_diff_Red(medidores,Caso.NB,Caso.NR,Barra,Ramo,Caso);
%generate_sets(Caso.NB,Caso.NR,redundancy,medidores,Barra,Ramo,Caso);
%generate_sets_UM_complete(Caso.NB,Caso.NR,nMedgen,medidores,Barra,Ramo,Caso,red);
fprintf('\nPower network with %d Buses, %d Branches and %d Measurements\n',Caso.NB, Caso.NR, length(medidores.num));
fprintf('\nMeasurement system design:\n\t1-From file:\tdefault\n\t2-Custom(complete measurement units):\tcustom\n\t3-Redund(Selected Redundancy):\tRed\n'); % selects the measurement design to be assessed
%Measurement_system_design=input('\nMeasurement_system_design:       ','s');
Measurement_system_design='default';
%Measurement_system_design='custom';

 
%% create branch to node incidence matrix
A=zeros(Caso.NR, Caso.NB);
MAdj= zeros(Caso.NB,Caso.NB);
for i=1:Caso.NR
    A(i,Ramo.de(i))=+1;
    A(i,Ramo.para(i))=-1;
    MAdj(Ramo.de(i),Ramo.para(i))=1;
    MAdj(Ramo.para(i),Ramo.de(i))=1; 
end

% create Jacobian matrix and measurement system with only real/active measurements
switch Measurement_system_design
    case 'custom'
        fprintf('\nNumber of system buses: %d\n',Caso.NB);
        %meas_location=zeros(1,nm); %reset meas_location
        mu_location=input('installed MUs (vector):   ');
        is_pmu=input('MU type of each MU: 0(conventional) or 1(phasor)');
        if(isempty(mu_location))
            mu_location=1:Caso.NB;
        end
        if(isempty(is_pmu))% all installed MUs are conventional
            is_pmu=zeros(1,Caso.NB);
        end
        % overwrite current measurement system
        [medidores]=generate_complete_munits(mu_location,is_pmu,A,Caso);
        [H , order] = LinJ(medidores,Ramo,Barra,Caso);
        medidores.num = medidores.num(order);
        medidores.de = medidores.de(order);
        medidores.para = medidores.para(order);
        medidores.tipo = medidores.tipo(order);
        medidores.PMU_num = medidores.PMU_num(order);
        medidores.leitura = medidores.leitura(order);
        [nm,ns]=size(H);
        
        %create Measurement unit and measurement location arrays
        meas_location=zeros(1,nm);
        for i=1:nm
            if(medidores.tipo(i)==1 || medidores.tipo(i)==7) % check if the ith measurement is a power flow/branch current
                pos=find(logical(mu_location==medidores.de(i)),1); %give the fort ocurrence of the bus in mu_location
                meas_location(i)=pos;
            else if (medidores.tipo(i)==2 || medidores.tipo(i)==3) % check if the ith measurement is a power injection or voltage angle
                    pos=find(logical(mu_location==medidores.para(i)),1); %give the fort ocurrence of the bus in mu_location
                    meas_location(i)=pos;
                end
            end
        end
        nmu=length(mu_location);
        
        
        %          nmu=length(mu_location);
        %          fprintf('\nNumber of measurements units: %d\n',length(mu_location));
        %          dim=length(mu_location);
        %          kmax=input('Maximum cardinality to be searched:   ');
        %          kmin=input('Minimum cardinality to be searched:   ');
        %          klim=kmax;
    case 'default'
        idz=logical( medidores.ok==0 & ( medidores.tipo==1 | medidores.tipo==2 | medidores.tipo==3 | medidores.tipo==7 | medidores.tipo==9) ); % medidores ativos de fluxo e injecao
        medidores.num = medidores.num(idz);
        medidores.de = medidores.de(idz);
        medidores.para = medidores.para(idz);
        medidores.tipo = medidores.tipo(idz);
        medidores.PMU_num = medidores.PMU_num(idz);
        medidores.leitura = medidores.leitura(idz);
        [H , order] = LinJ(medidores,Ramo,Barra,Caso);
        medidores.num = medidores.num(order);
        medidores.de = medidores.de(order);
        medidores.para = medidores.para(order);
        medidores.tipo = medidores.tipo(order);
        medidores.PMU_num = medidores.PMU_num(order);
        medidores.leitura = medidores.leitura(order);
        [nm,ns]=size(H);
        %create Measurement unit and measurement location arrays
        meas_location=zeros(1,nm);
        mu_location=zeros(1,Caso.NB);
        nmu=0;
        for i=1:nm
            if(medidores.tipo(i)==1 || medidores.tipo(i)==7) %check if it is a power flow/branch current measurement
                pos=find(logical(mu_location==medidores.de(i)),1); %get the first ocurrence of the bus in mu_location
                if(isempty(pos)) % check if the bus is present in mu_location
                    nmu=nmu+1; %increase the number of measurement units
                    mu_location(nmu)=medidores.de(i); %add pmu location
                    meas_location(i)=nmu;
                else
                    meas_location(i)=pos;
                end
            else if (medidores.tipo(i)==2 || medidores.tipo(i)==3 || medidores.tipo(i)==9) %check if it is a active power/ real current injection or a voltage phase measurement
                    pos=find(logical(mu_location==medidores.para(i)),1); %get the first ocurrence of the bus in mu_location
                    if(isempty(pos)) % check if the bus is present in mu_location
                        nmu=nmu+1;
                        mu_location(nmu)=medidores.para(i); %add pmu location
                        meas_location(i)=nmu;  %increase the number of measurement units
                    else
                        meas_location(i)=pos;
                    end
                end
            end
        end
        mu_location((nmu+1):end)=[];% shrink mu_location       
end

%% #####################################Pre-Processo###################################################################

[Barras_Livres,UMs_Livres]= UMs_complementares(mu_location,Caso.NB);
UM= repmat(struct('H_parc',[],'Barra',[],'Num_Medidas',[]),1,UMs_Livres);
UM= montaH(UM,MAdj,Barras_Livres,H);% Medidas Scada
disp(UM)
barras_atuais=Barras_Livres;
lote=3;%tamanho do lote
%w_k=[10 1];% pesos 118_base
%w_k=[1 5];% pesos 118_pesado
w_k=[10 5 2];%pesos 30
Caso_Otim=[num2str(Caso.NB) 'BUS_' 'lote_' num2str(lote) '_UMs_' 'GRASP_Viz3' ]; % Define Caso

G_inicial=H'*H;
E=eye(nm,nm)-((H)*(G_inicial^-1)*H');
dlmwrite('Covariancia.txt', E, 'delimiter', ' ', 'precision', '%.15f');
dlmwrite('nmed.txt',nm);
dos('Crit_find_CPU.exe');
num_crit=(dlmread('Criticalidades.txt'))';
best_solution_inicial=sum((w_k.*num_crit)); % valor

%% Avalia Solucao Pronta
% 
% solution=zeros(1,UMs_Livres);
% 
% solution=[1 0 0 1 0 0 0 0 0 0 1 0 ];
% %solution=[1 1 1 1 1 1 1 1 1 1 1 1 ];

% % Atualiza a best_solution
% if  best_fobj > fobj_atual
%     
%     best_fobj=fobj_atual;
%     best_crit= num_crit;
%     
% end

%% ######################################Otimizacao########################################################################### 
alpha_vec=[0.3 0.5 0.8];
index_alpha=1;
time_max=input('Digite o tempo max em min');
while index_alpha<=length(alpha_vec) % testa varios alfas
teste=1; % contador de testes
max_testes=3;
best_crit_testes=[];
best_sol_testes=[]; % tupla binaria
start_testes= tic; 
end_testes=[];
max_it=10;      
solutions_path=zeros(max_testes,max_it+1); % Caminho da otimizacao contando a condicao inicial do sistema
while teste<=max_testes
%% CONSTRUTIVO GRASP - baseado no numero de Medidas (Cria uma solucao avaliando pela quantidade de medidas disponivel pela barra)
    iter_max=1;
    %best_sol=1000000;
    best_solution= best_solution_inical;
    best_crit=[];
    best_sol=[];
    start_test=tic;     
    solutions_path(teste,iter_max)=best_solution; % insere a condicao inicial do sistema
    %% Inicia Tabela Hash/ Lista Tabu de Vizinhos - Memoria Global
    M = containers.Map('KeyType','char','ValueType','any'); 
    while iter_max<=max_it &&((toc(start_test)/60)<time_max)
        solution=zeros(1,UMs_Livres);
        sol_aux=ones(1,UMs_Livres);
        ordenacao1=(1:UMs_Livres);
        orden_aux=(1:UMs_Livres);%reduz de tamanho
        alpha=alpha_vec(index_alpha);
        LC=[];
        for i=1:length(orden_aux)
            qualidade(1,i)= UM(orden_aux(i)).Num_Medidas;
            LC=[LC;orden_aux(i)];
        end
        a=[qualidade;LC'];
        [idx,order]=sort(a(1,:),'descend');
        b= a(:,order);
        LC=b(2,:);
        orden_aux1=nonzeros(LC.*sol_aux);
        orden_aux2=orden_aux1;
        for i=1:length(orden_aux2)
            pos_orig(orden_aux2(i),1)= orden_aux2(i);
            pos_orig(orden_aux2(i),2)= i;
        end
        while sum(solution)<lote
            tam_LRC= ceil(max(1,length(orden_aux1)*alpha));
            LRC=orden_aux1(1:tam_LRC);
            random_pos=randi(length(LRC));
            solution(LRC(random_pos))= 1;
            sol_aux(pos_orig(LRC(random_pos),2))=0;
            orden_aux1=nonzeros(orden_aux2'.*sol_aux);
        end
        actual_fobj=0;
        actual_crit=[];
        [H_temp]= remonta_H(solution,barras_atuais,UM,H,lote);
        [fobj_atual,num_crit]= eval_solution(H_temp,H,w_k,nm);
        actual_fobj=fobj_atual; % fobj construtivo
        actual_crit=num_crit; 
        actual_solution=solution; % solution construtivo
        %tempo_total=0;
        %iter_max=0;
 
 %% processo iterativo - BUSCA LOCAL VIZ UNICA
    k=1;
%     %% Inicia Tabela Hash
%     M = containers.Map('KeyType','char','ValueType','any'); % Evita muitas repeticoes com calculo da fobj
    while (k==1)
        [sol_viz,fobj_viz,crit_viz] = est_viz3_altern(actual_solution);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Avalia Soluçao Vizinhaca %%%%%%%%%%%%%%%%%%%%%%%%%%%
       if (fobj_viz < actual_fobj)
            actual_fobj=fobj_viz; % valor da funcao objetivo
            actual_crit=crit_viz; % qtde de criticalidades
            actual_solution= sol_viz; % tupla binaria referencia da solucao 
       else
            k=k+1;
       end
       
    end
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Avalia Soluçao c/ Best %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (actual_fobj < best_solution)
        best_solution=actual_fobj; % valor da funcao objetivo
        best_crit=actual_crit; % qtde de criticalidades
        best_sol= actual_solution; % tupla binaria referencia da solucao 
    end
    
    solutions_path(teste,iter_max+1)=best_solution; % Caminho da otimizacao
    iter_max=iter_max+1;
% tic
% while (tempo_total < 5)&&(iter_max<150)
%    
%     [H_temp]= remonta_H(solution,barras_atuais,UM,H,lote);
%     [fobj_atual,num_crit]= eval_solution(H_temp,H,w_k,nm);
%     tempo_total=toc/60; 
%     iter_max=iter_max+1;
%     disp('ok')
%     
% end

 end
    end_GRASP= toc(start_test);  % count elapsed time for a complete GRASP
    best_crit_testes=[best_crit_testes;best_crit];
    best_sol_testes=[best_sol_testes;best_sol];
    end_testes=[end_testes;end_GRASP];
    teste=teste+1;
end

time_total_test=toc(start_testes);
%end_testes=toc(start_testes); % count elapsed time for all tests

%% PLOTAGENS E SALVAMENTOS
dlmwrite([Caso '_Criticalidades_' '_alpha_' num2str(alpha_vec(index_alpha)) '.txt'],best_crit_testes);
dlmwrite([Caso '_Tuplas_bin_' '_alpha_' num2str(alpha_vec(index_alpha)) '.txt'],best_sol_testes);
dlmwrite([Caso '_Fobj_values_' '_alpha_' num2str(alpha_vec(index_alpha)) '.txt'],solutions_path);
dlmwrite([Caso 'total_time_' '_alpha_' num2str(alpha_vec(index_alpha)) '.txt'],end_testes);

index_alpha=index_alpha+1; % testa outro alpha
end

disp('saiu')
disp('OK')