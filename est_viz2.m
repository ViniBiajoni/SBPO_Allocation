% Troca 2 posicoes
function [s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz2(s_tupla)
    
   %% Inicializa sol inicial vizinhanca
    best_viz_fobj=0;
    best_viz_tupla=[];
    best_viz_crit=[];
    
    %% Avalia a vizinhanca dentro de um numero maximo de iteracoes
    iter=0;
    iter_max=100;
    %time_max=tic;
    sol_base=s_tupla;
    ativas= (1:length(s_tupla)).*(s_tupla);
    n_ativas= (1:length(s_tupla)).*(s_tupla(:)==0)';
    ativas=(nonzeros(ativas))';
    n_ativas=(nonzeros(n_ativas))';
    %last_randi=0; % usar caso deseje criar um valor tabu
    
    while (iter<=iter_max) %&& ((toc(time_max)/60)<=1)
            
        vizinho=sol_base; % salva sol base vizinho a ser avaliado
        vizinho= troca_2_posicoes(vizinho,ativas,n_ativas);
                
        %% Avalia vizinho
        actual_crit=[];
        [H_temp]= remonta_H(vizinho,barras_atuais,UM,H,lote);
        [fobj_atual,num_crit]= eval_solution(H_temp,H,w_k,nm);
        viz_fobj=fobj_atual;
        actual_crit=num_crit;
        
        if (best_viz_fobj == 0)||(best_viz_fobj>viz_fobj)
            
                 best_viz_fobj=viz_fobj;
                 best_viz_tupla=vizinho;
                 best_viz_crit=actual_crit;
        end
        
        iter=iter+1;
    end
    
 %% Retorna Melhor vizinho
    s_viz_sol=best_viz_tupla;
    s_viz_fobj=best_viz_fobj;
    s_viz_crit=best_viz_crit;

end