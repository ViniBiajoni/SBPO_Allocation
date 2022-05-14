%% lendo o arquivo e montando matriz de conexão
clc
clear all
Conexao = dlmread('conex57.txt');
[numero_barras,c] = size(Conexao);
for i = 1:numero_barras
    for j = 1:numero_barras
        if i==j 
            Conexao(i,j)=1;
        end
    end
end

tamc = numero_barras;
b=ones(tamc,1);
f=ones(tamc,1);
b=b*-1;
Conexao=Conexao*-1;
intcon=(1:tamc);
lb=zeros(tamc,1);
ub=ones(tamc,1);
y=intlinprog(f,intcon,Conexao,b,[],[],lb,ub);
barras_alocadas = y.*(1:numero_barras)';
barras_alocadas = nonzeros(barras_alocadas);
dlmwrite('Barras_Alocadas.txt', barras_alocadas);