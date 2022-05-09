function [medidores]=generate_complete_munits(mu_location,is_pmu,A,Caso)
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

% Prof. Milton's convention
% Conventional meaurement unit : bus power injection (tipo 2) and power flows (tipo 1) at every incident branch
% Conventional meaurement unit : bus voltage phase (tipo 3) and real branch current (tipo 7) at every incident branch

medidores=struct('tipo',[],'de',[],'para',[],'num',[],'leitura',[], 'PMU_num',[]);
k=0;
for i=1:length(mu_location)
    k=k+1; % increase counter
    
    if (~is_pmu(i))% Measurement unit is not a phasor one
        medidores.tipo=[medidores.tipo 2]; %add active power injection
        medidores.PMU_num=[medidores.PMU_num,-i]; % PMU_num is negative for conventional MU
    else
        medidores.tipo=[medidores.tipo 3]; %add voltage phase      
        medidores.PMU_num=[medidores.PMU_num,+i]; % PMU_num is positive for Phasor MU
    end
    medidores.de=[medidores.de 0];
    medidores.para=[medidores.para mu_location(i)];
    medidores.leitura=[medidores.leitura k];
    medidores.num=[medidores.num k];
    
    for j=1:Caso.NR
        if(abs(A(j,mu_location(i)))==1) %add flow measurement
            k=k+1; % increase counter
            if (is_pmu(i)==0) % Measurement unit is not a phasor one
                medidores.tipo=[medidores.tipo 1]; %add a branch power flow
                medidores.PMU_num=[medidores.PMU_num,-i]; % PMU_num is negative for conventional MU
            else
                medidores.tipo=[medidores.tipo 7]; %add a real branch current
                medidores.PMU_num=[medidores.PMU_num,+i]; % PMU_num is positive for Phasor MU
            end
            medidores.de=[medidores.de mu_location(i)];
            if (A(j,mu_location(i))==+1) % locate bus to
                bus_to=find(logical((A(j,:)==-1)),1);
            else
                bus_to=find(logical((A(j,:)==+1)),1);
            end
            medidores.para=[medidores.para bus_to];
            medidores.leitura=[medidores.leitura k];
            medidores.num=[medidores.num k];
        end
    end
end
