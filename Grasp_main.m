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
%% IEEE 57 bus TEST SYSTEM
med_file='Caso57b104m.txt'; %Plan1
%med_file='S1_CS5_352M_90M_PMU.med'; %Plan2

%#####################
sis_file='ieee57.cdf';
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

%Intilinprog2

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
UM= montaH(UM,MAdj,Barras_Livres,H,Caso.NB);% Medidas Scada
disp(UM)
barras_atuais=Barras_Livres;
lote=input('Digite o tamanho do lote');%tamanho do lote
% w_k=[20 1.5];% pesos 118_base(a)
w_k=[1 1];% pesos 118_(b)
%w_k=[10 3 1.5];%pesos 30
Caso_Otim=[num2str(Caso.NB) 'BUS_' 'lote_' num2str(lote) '_UMs_' 'GRASP_Completo' ]; % Define Caso

G_inicial=H'*H;
E=eye(nm,nm)-((H)*(G_inicial^-1)*H');
dlmwrite('Covariancia.txt', E, 'delimiter', ' ', 'precision', '%.15f');
dlmwrite('nmed.txt',nm);
tic
dos('Crit_find_CPU.exe');
toc
num_crit=(dlmread('Criticalidades.txt'))';
disp(num_crit)
best_solution_inical=sum((w_k.*num_crit)); % valor


%% Testa Guloso - Para Ver quao Perto a Metaheuristica Está Chegando
start_exaustivo= tic;
combos= combntns(1:UMs_Livres,lote);
best_solution=0;
best_tuple=[];
best_crit=[];
for i=1:length(combos(:,1))
    disp('Teste\n')
    disp(i)
    disp('\n----')
    sol=zeros(1,UMs_Livres);
   for j=1:length(combos(1,:))
      sol(combos(i,j))=1;
   end
    actual_fobj=0;
    actual_crit=[];
    [H_temp]= remonta_H(sol,barras_atuais,UM,H,lote);
    [fobj_atual,num_crit]= eval_solution(H_temp,H,w_k,nm);
    actual_fobj= fobj_atual;
    actual_crit= num_crit; 
    current_sol= sol;
    
    if (best_solution==0)||(best_solution>actual_fobj)
        best_solution= actual_fobj;
        best_tuple=current_sol;    
        best_crit=actual_crit;
    end
    toc(start_exaustivo)
    allocated_buses=nonzeros(best_tuple.*barras_atuais);
%     if (toc(start_exaustivo)/60)>=90
%         break;
%     end
end
toc(start_exaustivo)
allocated_buses=nonzeros(best_tuple.*barras_atuais);

%% ######################################Otimizacao########################################################################### 

%alpha_vec=[0.3 0.5 0.8];
%alpha_vec=[0.3 0.8];%118 novo
alpha_vec(1)=0.1;%118 novo teste 2
index_alpha=1;
time_max=input('Digite o tempo max em min');
max_testes=input('Digite o numero de testes');
max_it=input('Digite o numero de iterações');
while index_alpha<=length(alpha_vec) % testa varios alfas

teste=1; % contador de testes
best_crit_testes=[];
best_sol_testes=[]; % tupla binaria
start_testes= tic; 
end_testes=[];
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
    %% Inicia Tabela Hash/ Lista Tabu 
    M = containers.Map('KeyType','char','ValueType','any'); 
while (iter_max<= max_it)&&((toc(start_test)/60)<time_max)%&&(no_change<5)
    solution=zeros(1,UMs_Livres);
    sol_aux=ones(1,UMs_Livres);
    ordenacao1=(1:UMs_Livres);
    orden_aux=(1:UMs_Livres);%reduz de tamanho
    alpha=alpha_vec(index_alpha);
    %alpha=0.5;
    %alpha=0.8;
    %alpha=0.3;
    LC=[];
    for i=1:length(orden_aux)
        qualidade(1,i)= UM(orden_aux(i)).Num_Medidas;
        LC=[LC;orden_aux(i)];
    end
    a=[qualidade;LC'];
    [idx,order]=sort(a(1,:),'descend');%ordena decrescente
    %[idx,order]=sort(a(1,:));%ordena crescente
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
    
    %% Teste estruturas de vizinhanca
    %[s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz1_altern(actual_solution);
    %[s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz2_altern(actual_solution);
    %[s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz3_altern(actual_solution);
    %[s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz4_altern(actual_solution);
 
 %% processo iterativo - BUSCA LOCAL VND
    sol_VND=[];
    fobj_VND=0;
    crit_VND=[];
    [sol_VND,fobj_VND,crit_VND] =VND(actual_solution,actual_fobj,actual_crit,Caso.NB);
    %[sol_VND,fobj_VND,crit_VND]=RVND(actual_solution,actual_fobj,actual_crit,Caso.NB);%Teste
    %eventual com ordenacao aleatoria do VND
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Avalia Soluçao %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (fobj_VND < best_solution)
        best_solution=fobj_VND; % valor da funcao objetivo
        best_crit=crit_VND; % qtde de criticalidades
        best_sol= sol_VND; % tupla binaria referencia da solucao 
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
%disp(index_alpha);
%disp(iter_max);
end
    end_GRASP= toc(start_test);  % count elapsed time for a complete GRASP
    best_crit_testes=[best_crit_testes;best_crit];
    best_sol_testes=[best_sol_testes;best_sol];
    end_testes=[end_testes;end_GRASP];
    teste=teste+1;
end

time_total_test=toc(start_testes);
%disp('alpha:')
%disp(alpha_vec(index_alpha))
%disp('Time_total=')
%disp(time_total)
%end_testes=toc(start_testes); % count elapsed time for all tests

%% PLOTAGENS E SALVAMENTOS
dlmwrite([Caso_Otim '_Criticalidades_' '_alpha_' num2str(alpha_vec(index_alpha)) '.txt'],best_crit_testes);%tupla
dlmwrite([Caso_Otim '_Tuplas_bin_' '_alpha_' num2str(alpha_vec(index_alpha)) '.txt'],best_sol_testes);%valor fpbj
dlmwrite([Caso_Otim '_Fobj_values_' '_alpha_' num2str(alpha_vec(index_alpha)) '.txt'],solutions_path); %caminho da otimizacao 
dlmwrite([Caso_Otim 'total_time_' '_alpha_' num2str(alpha_vec(index_alpha)) '.txt'],end_testes);%tempo

index_alpha=index_alpha+1; % testa outro alpha
end

disp('saiu')
disp('ok')