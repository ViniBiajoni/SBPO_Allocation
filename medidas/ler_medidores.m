function [medidores, Nmed]=ler_medidores(med_file)
% Le arquivo de dados de medidores
%
% ler_medidores              -> Solicita o nome do arquivo de dados
% ler_medidores(med_file)    -> Utiliza o arquivo especificado
% 
% Estrutura de medidores
% num           nº do medidor
% de            nº da barra DE
% para          nº da barra PARA
% circ          nº do circuito
% tipo          tipo de medidor
% PMU_num       Medida associada a PMU n
% ok            Utiliza ou nao o medidor
% acc           exatidão do medidor
% fs            fundo de escala do medidor
% dp            desvio-padrão
% ref           valor de referencia
% leitura       valor medido
%
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
%
medidores.num = [];
medidores.de = [];
medidores.para = [];
medidores.circ = [];
medidores.tipo = [];
medidores.PMU_num = [];
medidores.ok = [];
medidores.acc = [];
medidores.fs = [];
medidores.dp = [];
medidores.ref = [];
medidores.leitura = [];
medidores.v_estimado = [];
medidores.residuo = [];
medidores.dpc = [];
medidores.rn =[];
if nargin < 1 %Avalia se a funcao recebeu argumentos ou nao
    disp('====================================================================')
    med_file = input('Ler medidores no arquivo:[medidores.dat] ', 's');
    if isempty(med_file); med_file = 'medidores.dat'; end
end
% Abre o arquivo para leitura
fid = fopen(med_file,'r');
% Ler o arquivo de medidores
[dados,count] = fscanf(fid, '%g', [12, inf]);%armazena as linhas em colunas
% Fecha o arquivo após a leitura
fclose(fid);
[linha, coluna] = size(dados);
% Insere os dados na estrutura de medidores
for i = 1 : coluna
    medidores.num(i)	= 		dados(1,i);
    medidores.de(i)		=		dados(2,i);
    medidores.para(i)	= 		dados(3,i);
    medidores.circ(i)	= 		dados(4,i);
    medidores.tipo(i)	= 		dados(5,i);
    medidores.PMU_num(i)=       dados(6,i);
    medidores.ok(i)		= 		dados(7,i);
    medidores.acc(i)	= 		dados(8,i);
    medidores.fs(i)		= 		dados(9,i);
    medidores.dp(i)		= 		dados(10,i);
    medidores.ref(i)	= 		dados(11,i);
    medidores.leitura(i)=       dados(12,i);
end
Nmed = length(medidores.num);
return

