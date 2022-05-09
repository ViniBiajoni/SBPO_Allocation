%Troca 1 pos e inverte vetor
function [s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz3_altern_hash_t(s_tupla)

 global barras_atuais
 global UM
 global lote
 global H
 global nm
 global w_k
 global M
 global start_test
 global time_max
 %% Inicializa sol inicial vizinhanca

 %% Avalia a solucao central
 fobj_central=0;
 crit_central=[];
%  [H_temp]= remonta_H(s_tupla,barras_atuais,UM,H,lote);
%  [fobj_central,crit_central]= eval_solution(H_temp,H,w_k,nm);

    %% Inicializa sol inicial vizinhanca
    best_viz_fobj=0;
    best_viz_tupla=[];
    best_viz_crit=[];
    %time_max=tic;
    sol_base=s_tupla;
    ativas= (1:length(s_tupla)).*(s_tupla);
    n_ativas= (1:length(s_tupla)).*(s_tupla(:)==0)';
    ativas=(nonzeros(ativas))';
    n_ativas=(nonzeros(n_ativas))';
    cont=0;
    %% Avalia o tamanho da instancia 
    comb_zeros=nchoosek(sum(s_tupla(:)==0),1);
    comb_ones=nchoosek(sum(s_tupla(:)==1),1);
    combinacoes=comb_zeros*comb_ones;
    abort=0;   
 if combinacoes<=50   
    %% Efetua Combinacoes
    for i=1:length(ativas)
       vizinho=sol_base; % vizinho a ser avaliado
       vizinho(ativas(i))=0; % desativo
       for j=1:length(n_ativas) 
          vizinho(n_ativas(j))=1; %ativo
          viz_temp=vizinho;
          vizinho=flip(vizinho); %gira o vizinho
          %% Avalia se ja foi analisado
          key=num2str(vizinho);
          if isKey(M,key)==0 % so testa se nao estiver na hash Table
            %% Avalia vizinho
            actual_fobj=0;
            actual_crit=[];
            [H_temp]= remonta_H(vizinho,barras_atuais,UM,H,lote);
            [actual_fobj,actual_crit]= eval_solution(H_temp,H,w_k,nm);
            M(key)=[actual_fobj actual_crit];%salva as criticalidades
                
            if (best_viz_fobj == 0)||(best_viz_fobj>actual_fobj)  
                best_viz_fobj=actual_fobj;
                best_viz_tupla=vizinho;
                best_viz_crit=actual_crit;
            end
            
          else
              vector_map= M(key);
              actual_fobj= vector_map(1);
              actual_crit= vector_map(2:end);
              if (best_viz_fobj == 0)||(best_viz_fobj>actual_fobj)  
                 best_viz_fobj=actual_fobj;
                 best_viz_tupla=vizinho;
                 best_viz_crit=actual_crit;
              end
          end
          
          viz_temp(n_ativas(j))=0; %desativa para andar para o proximo 
          vizinho=viz_temp;
          if (toc(start_test)/60)>=time_max
              abort=1;  
              break;    
          end
       end
       
       if abort==1
            break;
       end
    end
    
 else
        iter_limite=50;
        iter=0;
        done=0;
        vizinho=sol_base; % vizinho a ser avaliado
       while done==0
            desliga=randi(length(ativas)); % posicao para desligar
            pos_deslig=ativas(desliga);
            vizinho(ativas(desliga))=0; %desliga primeiro vizinho
            ativas(desliga)=0;
            ativas=(nonzeros(ativas))';%desabilita posicao para futuros sorteios
        for i=1:length(n_ativas)
          vizinho(n_ativas(i))=1;
          viz_temp=vizinho;
          vizinho=flip(vizinho); %gira o vizinho
          %% Avalia se ja foi analisado
          key=num2str(vizinho);
          if isKey(M,key)==0 % so testa se nao estiver na hash Table
            %% Avalia vizinho
            actual_fobj=0;
            actual_crit=[];
            [H_temp]= remonta_H(vizinho,barras_atuais,UM,H,lote);
            [actual_fobj,actual_crit]= eval_solution(H_temp,H,w_k,nm);
            M(key)=[actual_fobj actual_crit];%salva as criticalidades
           
            if (best_viz_fobj == 0)||(best_viz_fobj>actual_fobj)  
                 best_viz_fobj=actual_fobj;
                 best_viz_tupla=vizinho;
                 best_viz_crit=actual_crit;
            end
          else
              
              vector_map= M(key);
              actual_fobj= vector_map(1);
              actual_crit= vector_map(2:end);
              if (best_viz_fobj == 0)||(best_viz_fobj>actual_fobj)  
                 best_viz_fobj=actual_fobj;
                 best_viz_tupla=vizinho;
                 best_viz_crit=actual_crit;
              end
            
          end
          viz_temp(n_ativas(i))=0; %desativa para andar para o proximo 
          vizinho=viz_temp;
          iter=iter+1;
          if (iter>iter_limite)||(toc(start_test)/60)>=time_max
              done=1;
              break;
          end
          
        end
        vizinho(pos_deslig)=1; %religa primeiro vizinho    
       end
         
 end
    
       %% Retorna Melhor vizinho
        s_viz_sol=best_viz_tupla;
        s_viz_fobj=best_viz_fobj;
        s_viz_crit=best_viz_crit;  
        


end