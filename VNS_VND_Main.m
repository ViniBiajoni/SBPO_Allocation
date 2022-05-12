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

med_file='Caso57b150m.txt';

sis_file='ieee57.cdf';

w_k=[1 3.7 1 1.5];

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
UM= montaH(UM,MAdj,Barras_Livres,H,Caso.NB);% Medidas Scada
disp(UM)
barras_atuais=Barras_Livres;

lote=input('Digite o tamanho do lote');



%% Escolhe tipo de Sol Inicial
tipo_inicio=input('Digite 0 para sol inicial Aleatoria e 1 para Gulosa \n');

if tipo_inicio==0
    tipo_initial_sol='Aleatoria';
else
    tipo_initial_sol='Gulosa'; 
end

%% Inicializa Caso
Caso_Otim=[num2str(Caso.NB) 'BUS_' 'lote_' num2str(lote) '_UMs_' 'VNS_VND_' tipo_initial_sol ]; % Define Caso

G_inicial=H'*H;
E=eye(nm,nm)-((H)*(G_inicial^-1)*H');
dlmwrite('Covariancia.txt', E, 'delimiter', ' ', 'precision', '%.15f');
dlmwrite('nmed.txt',nm);
%tic
dos('CritFindGPU.exe');
%toc
num_crit=(dlmread('Criticalidades.txt'))';
best_solution_inicial=sum((w_k.*num_crit)); % valor


%% Testa Guloso - Para Ver quao Perto a Metaheuristica Está Chegando
% combos= combntns(1:UMs_Livres,lote);
% best_solution=0;
% best_tuple=[];
% best_crit=[];
% for i=1:length(combos(:,1))
%     sol=zeros(1,:);
%    for j=1:length(combos(1,:))
%       sol(combos(i,j))=1;
%    end
%     actual_fobj=0;
%     actual_crit=[];
%     [H_temp]= remonta_H(sol,barras_atuais,UM,H,lote);
%     [fobj_atual,num_crit]= eval_solution(H_temp,H,w_k,nm);
%     actual_fobj= fobj_atual;
%     actual_crit= num_crit; 
%     current_sol= solution;
%     
%     if (best_solution==0)||(best_solution>actual_fobj)
%         best_solution= actual_fobj;
%         best_tuple=current_sol;    
%         best_crit=actual_crit;
%     end
%     
% end

%% ######################################Otimizacao VNS-VND########################################################################### 
teste=1; % contador de testes
max_testes=input('Digite o maximo de testes');
best_crit_testes=[];
best_sol_testes=[]; % tupla binaria
start_testes= tic; 
end_testes=[];
max_it=input('Digite o maximo de iteracoes');    
time_max=input('Digite o tempo max em min');
solutions_path=zeros(max_testes,max_it+1); % Caminho da otimizacao contando a condicao inicial do sistema
while teste<=max_testes
    iter_max=1;
    %best_sol=1000000;
    best_crit=[];
    best_sol=[];
    start_test=tic;
    
%% solucao inicial
% Pegando as UMs com mais medidas disponíveis
if tipo_inicio==0
 
% Aleatoria
     complete =0;
     solution=zeros(1,UMs_Livres);
     while (complete ==0)
         sorteio=randi(UMs_Livres);
         if (solution(sorteio)==0)
             solution(sorteio)=1;
         end
 
         if (sum(solution)==lote)
             complete=1;
         end
     end
     
%Guloso
else
    solution=zeros(1,UMs_Livres);
    for i=1:UMs_Livres
        temp(i)= i; 
        barras(i)= UM(i).Barra;
        meds(i)= UM(i).Num_Medidas;
    end
%     prior=[30 33 40 54 70];
%     for j=1:length(prior)
%        solution(prior(j))=1;%coloca barras prioritarias
%        barras(j)=0;
%        meds(j)=0;
%        temp(j)=0;
%     end
%     barras=nonzeros(barras)';
%     meds=nonzeros(meds)';
%     temp=nonzeros(temp)';
    a=[temp;barras;meds];
    [idx,order]=sort(a(3,:),'descend'); % ordena decrescente
    %[idx,order]=sort(a(3,:)); % ordena crescente
    temp2=a(1,order);
    for i=1:lote%-length(prior))
   
        solution(temp2(1,i))=1;
    
    end

end

best_fobj=0;
best_crit=[];
[H_temp]= remonta_H(solution,barras_atuais,UM,H,lote);
[fobj_atual,num_crit]= eval_solution(H_temp,H,w_k,nm);
%% solucao inicial
best_fobj=fobj_atual;
best_crit=num_crit;
best_sol=solution; 

%% Inicializa Caminho das Solucoes
solutions_path(teste,iter_max)=best_fobj; % insere a condicao inicial do sistema
    
%% Inicia Tabela Hash/ Lista Tabu - Acelera as Iteraçoes Buscando Apenas Novos Vizinhos 
M = containers.Map('KeyType','char','ValueType','any'); 

while (iter_max<= max_it)&&((toc(start_test)/60)<time_max)
    
   k=1;
    %% Teste estruturas de vizinhanca
    %[s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz1_altern(actual_solution);
    %[s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz2_altern(actual_solution);
    %[s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz3_altern(actual_solution);
    %[s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz4_altern(actual_solution);
  while (k<=4)&&((toc(start_test)/60)<time_max)
 %% Gera Vizinho Qualquer na Vizinhanca k
    [actual_solution,actual_fobj,actual_crit]= gera_viz_VNS(best_sol,k); % gera uma vizinhanca aleatoria
 %% processo iterativo - BUSCA LOCAL VND
    sol_VND=[];
    fobj_VND=0;
    crit_VND=[];
    [sol_VND,fobj_VND,crit_VND] =VND(actual_solution,actual_fobj,actual_crit,Caso.NB);
    %[sol_VND,fobj_VND,crit_VND]=RVND(actual_solution,actual_fobj,actual_crit,Caso.NB);%Teste
    %eventual com ordenacao aleatoria do VND
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Avalia Soluçao %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (fobj_VND < best_fobj)
        best_fobj=fobj_VND; % valor da funcao objetivo
        best_crit=crit_VND; % qtde de criticalidades
        best_sol= sol_VND; % tupla binaria referencia da solucao 
        k=1;
    else
        k=k+1;
    
    end
    

  end

iter_max=iter_max+1;
solutions_path(teste,iter_max)=best_fobj; % Caminho da otimizacao
disp(iter_max)
end
    end_VNS= toc(start_test);  % count elapsed time for the VNS-VND iterations
    best_crit_testes=[best_crit_testes;best_crit];
    best_sol_testes=[best_sol_testes;best_sol];
    end_testes=[end_testes;end_VNS];
    teste=teste+1;
end
time_total_test=toc(start_testes);% count elapsed time for all tests

best_sol_testes = print_UM_Sol(max_testes,best_sol_testes,Barras_Livres,lote);

%end_testes=toc(start_testes); 
%% PLOTAGENS E SALVAMENTOS
dlmwrite([Caso_Otim '_Criticalidades_' '.txt'],best_crit_testes);%tupla
dlmwrite([Caso_Otim '_Tuplas_bin_' '.txt'],best_sol_testes);%valor fobj
dlmwrite([Caso_Otim '_Fobj_values_' '.txt'],solutions_path); %caminho da otimizacao 
dlmwrite([Caso_Otim 'total_time_' '.txt'],end_testes);%tempo por iteracao

disp('saiu')