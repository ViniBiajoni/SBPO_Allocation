function [s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz5(s_tupla)
    
% Vizinhanca 1 deterministica
% Troca uns por zeros

    %% Inicializa sol inicial vizinhanca
    best_viz_fobj=0;
    best_viz_tupla=[];
    best_viz_crit=[];
    
    %% Avalia a vizinhanca dentro de um numero maximo de iteracoes
    %time_max=tic;
    sol_base=s_tupla;
    ativas= (1:length(s_tupla)).*(s_tupla);
    n_ativas= (1:length(s_tupla)).*(s_tupla(:)==0)';
    ativas=(nonzeros(ativas))';
    n_ativas=(nonzeros(n_ativas))';
    
    for i=1:length(ativas)
       vizinho=sol_base; % vizinho a ser avaliado
       vizinho(ativas(i))=0; % desativo
       for j=1:length(n_ativas) 
          vizinho(n_ativas(i))=1; %ativo 
          %% Avalia vizinho
          actual_fobj=0;
          actual_crit=[];
          [H_temp]= remonta_H(vizinho,barras_atuais,UM,H,lote);
          [fobj_atual,num_crit]= eval_solution(H_temp,H,w_k,nm);
           viz_fobj=fobj_atual;
          if (best_viz_fobj == 0)||(best_viz_fobj>actual_fobj)  
                 best_viz_fobj=actual_fobj;
                 best_viz_tupla=vizinho;
                 best_viz_crit=actual_crit;
          end
       end         
    end
    
 %% Retorna Melhor vizinho
    s_viz_sol=best_viz_tupla;
    s_viz_fobj=best_viz_fobj;
    s_viz_crit=best_viz_crit;

end

