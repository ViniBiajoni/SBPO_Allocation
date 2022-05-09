function [actual_solution,actual_fobj,actual_crit]= gera_viz_VNS(sol_atual,k)
global barras_atuais
global UM
global lote
global H
global nm
global w_k
s_tupla=sol_atual;
if k==1
    
    %% Troca uma Pos Genérica
    best_viz_fobj=0;
    best_viz_tupla=[];
    best_viz_crit=[];
    ativas= (1:length(s_tupla)).*(s_tupla);
    n_ativas= (1:length(s_tupla)).*(s_tupla(:)==0)';
    ativas=(nonzeros(ativas))';
    n_ativas=(nonzeros(n_ativas))';        
    vizinho=s_tupla; % vizinho a ser avaliado
    desliga=randi(length(ativas)); % posicao para desligar
    liga=randi(length(n_ativas)); % posicao para ligar
    vizinho(ativas(desliga))=0;
    vizinho(n_ativas(liga))=1;
         

end

if k==2
    
    %% Inicializa sol inicial vizinhanca
    best_viz_fobj=0;
    best_viz_tupla=[];
    best_viz_crit=[];
    ativas= (1:length(s_tupla)).*(s_tupla);
    n_ativas= (1:length(s_tupla)).*(s_tupla(:)==0)';
    ativas=(nonzeros(ativas))';
    n_ativas=(nonzeros(n_ativas))';
    vizinho= troca_2_posicoes(s_tupla,ativas,n_ativas);
               
end


if k==3
   
   %% Troca uma Pos Genérica e Rotaciona o Vetor
    best_viz_fobj=0;
    best_viz_tupla=[];
    best_viz_crit=[];
    ativas= (1:length(s_tupla)).*(s_tupla);
    n_ativas= (1:length(s_tupla)).*(s_tupla(:)==0)';
    ativas=(nonzeros(ativas))';
    n_ativas=(nonzeros(n_ativas))';        
    vizinho=s_tupla; % vizinho a ser avaliado
    desliga=randi(length(ativas)); % posicao para desligar
    liga=randi(length(n_ativas)); % posicao para ligar
    vizinho(ativas(desliga))=0;
    vizinho(n_ativas(liga))=1;
    vizinho=flip(vizinho);
     
end


if k==4
    
   %% Inicializa sol inicial vizinhanca-> Troca 2 e inverte
    best_viz_fobj=0;
    best_viz_tupla=[];
    best_viz_crit=[];
    ativas= (1:length(s_tupla)).*(s_tupla);
    n_ativas= (1:length(s_tupla)).*(s_tupla(:)==0)';
    ativas=(nonzeros(ativas))';
    n_ativas=(nonzeros(n_ativas))';
    vizinho= troca_2_posicoes(s_tupla,ativas,n_ativas);
    vizinho=flip(vizinho);
    
end


%% Avalia vizinho
actual_crit=[];
[H_temp]= remonta_H(vizinho,barras_atuais,UM,H,lote);
[fobj_atual,num_crit]= eval_solution(H_temp,H,w_k,nm);
actual_fobj=fobj_atual;
actual_crit=num_crit;
actual_solution=vizinho;

end

