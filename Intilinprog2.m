%% lendo o arquivo e montando matriz de conexão
%     A= (dlmread('conex118.txt'));
%     A= dlmread('conex30.txt');
%     A= dlmread('conex24.txt');
    A = dlmread('conex57.txt');
    A=A(:,1:2);
    [m,n]=size(A);
       
    for i =1:m
        Conexao(A(i,1),A(i,2))=1;
        Conexao(A(i,2),A(i,1))=1;
    end
    
  
    for i = 1:tamc
        for j = 1:tamc
            if i==j 
                Conexao(i,j)=1;
            end
        end
    end


Conexao= MAdj;

b=ones(tamc,1);
f=ones(tamc,1);
b=b*-1;
Conexao=Conexao*-1;
intcon=(1:tamc);
lb=zeros(tamc,1);
ub=ones(tamc,1);
y=intlinprog(f,intcon,Conexao,b,[],[],lb,ub);

num_UMs=input('Digite o numero de UMs a serem alocadas');
[taml,tamc]=size(Conexao);
% caso_med=dlmread('ieee118CaseFULLPF.med');
%%Injecoes
% for i=1:Caso.NB
%    caso_med() 
% end
%% Fluxos

for cases=1:10
caso_med(:,7)=1;    
yaux=zeros(sum(y(:)),1);
barras_livres=zeros(118-length(yaux),1);
cont=1;
cont2=1;
%% Projeto

for i=1:tamc
    if (y(i)==1)
       yaux(cont,1)=i; 
       cont=cont+1;
    else
        barras_livres(cont2,1)=i;
        cont2=cont2+1;
    end
end

for sorteio=1:num_UMs
   
    value= randi(length(barras_livres));
    yaux=[yaux;barras_livres(value)];
    barras_livres(value)=0;
    barras_livres= nonzeros(barras_livres)';
end

for i=1:length(yaux)
    for j=1:length(caso_med(:,1))
        if (caso_med(j,2)==0 && caso_med(j,3)==yaux(i))||(caso_med(j,2)==yaux(i))
            caso_med(j,7)=0;
        end
    end
    
end
%% Name the file STAND ALONE
tam = num2str(tamc);
fileID = ['118_teste_' num2str(num_UMs) '_UMs' '_case_' num2str(cases) '.med'];
dlmwrite(fileID,caso_med,'delimiter',' ');
%name=[tam 'bus_begin.txt'];
%dlmwrite(name,yaux);
end

