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

%% IEEE 30 bus TEST SYSTEM
%med_file='ieee30_observability_PSCC_2014.med';%Caso Base 37 meas 20 MU
%med_file='S1_CS10_81M.med';
%med_file='S1_CS10_81M_62M_PMU26491012.med';
%med_file='ieee30_observability_full_2014.med';
%econk1= dlmread('IEEE30PSCCk1.txt');
%med_file='ieee30_observability_PSCC_2014B.med'; % Caso 42 meas
%med_file='ieee30_observability_PSCC_2014B_43.med'; % Caso 43 meas
%sis_file='ieee30.cdf';
%##########################################################################
%% IEEE 118 bus TEST SYSTEM
med_file='S1_CS5_352M.med'; %Plan1
%med_file='S1_CS5_352M_90M_PMU.med'; %Plan2
%econk1= dlmread('IEEE118k1.txt');
% med_file='118_teste.med';
%#####################
sis_file='ieee118.cdf';
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
w_k=[20 1.5];% pesos 118_base(a)
% w_k=[1 1];% pesos 118_(b)
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
% combos= (combntns(UMs_Livres,lote));
% num_waves = ;
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
    %% Print Parcial
    toc(start_exaustivo);
    allocated_buses=nonzeros(best_tuple.*barras_atuais);
%     fileID = fopen('Teste_busca_exaustiva_118a_lote6.txt','w');
    fileID = fopen('Teste_busca_exaustiva_118a_lote9.txt','w');
%     fileID = fopen('Teste_busca_exaustiva_118b_lote10.txt','w');
    time = toc(start_exaustivo);
    fprintf(fileID,'Time = %2.1f seconds\n',time);
    fprintf(fileID,'Fobj = %6.2f\n',best_solution);
    fclose(fileID);
    %Printa Alocacao
%     fileID_tuple = 'Alocacoes_busca_exaustiva_118a_lote6.txt';
      fileID_tuple = 'Alocacoes_busca_exaustiva_118a_lote9.txt';
%     fileID_tuple = 'Alocacoes_busca_exaustiva_118b_lote10.txt';
    dlmwrite(fileID_tuple,allocated_buses');
     if (time/60)>=600
          break;
     end
end
toc(start_exaustivo)
allocated_buses=nonzeros(best_tuple.*barras_atuais);