% Troca duas posicoes 
% A determinacao se é deterministica ou nao depende do tamanho da instancia
% Se demorar muito terei que sair por tempo
function [s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz2_altern(s_tupla)

 global barras_atuais
 global UM
 global lote
 global H
 global nm
 global w_k
 global M
 global start_test
 global time_max
%% Avalia a solucao central
 fobj_central=0;
 crit_central=[];
 [H_temp]= remonta_H(s_tupla,barras_atuais,UM,H,lote);
 [fobj_central,crit_central]= eval_solution(H_temp,H,w_k,nm);
    
%% Inicializa sol inicial vizinhanca
    best_viz_fobj=0;
    best_viz_tupla=[];
    best_viz_crit=[];
    
    %% Avalia o tamanho da instancia
    
    comb_zeros=nchoosek(sum(s_tupla(:)==0),2);
    comb_ones=nchoosek(sum(s_tupla(:)==1),2);
    %comb_ones=500;
    
    
    %% Avalia a vizinhanca dentro de um numero maximo de iteracoes
    iter=0;
    %iter_max=100;
    %time_max=tic;
    sol_base=s_tupla;
    ativas= (1:length(s_tupla)).*(s_tupla);
    n_ativas= (1:length(s_tupla)).*(s_tupla(:)==0)';
    ativas=(nonzeros(ativas))';
    n_ativas=(nonzeros(n_ativas))';
    %last_randi=0; % usar caso deseje criar um valor tabu
    cont=0;
    if (comb_zeros*comb_ones <= 200)
        
        finish_now=0;% force finish 
        
        for i=1:length(ativas)
            vizinho=s_tupla;
            vizinho(ativas(i))=0;%desliga primeiro
            for j=i+1:length(ativas)
                vizinho(ativas(j))=0;%desliga segundo
                for n=1:length(n_ativas)
                    vizinho(n_ativas(n))=1;%liga primeiro
                    for m=n+1:length(n_ativas)
                        vizinho(n_ativas(m))=1;%liga segundo
                        key=num2str(vizinho);
                        if isKey(M,key)==0 % so testa se nao estiver na hash Table
                            cont=cont+1;
                            M(key)=vizinho;
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
                        end
                      vizinho(n_ativas(m))=0;% reinicia segundo ligado
                      if (toc(start_test)/60)>=time_max
                        finish_now=1;
                        break;
                      end
                    end
                      
                      if finish_now==1
                         break; 
                      end
                      vizinho(n_ativas(n))=0;% reinicia primeiro ligado
                end
                vizinho(ativas(j))=1;% reinicia segundo desligado
                if finish_now==1
                    break; 
                end
            end   
            
            if finish_now==1
               break; 
            end
            
        end
     
          
    else
      %%  Avalia parte dos vizinhos  de forma semi-aleatória
      
      complete=1;
      iter_total=0;
      total_time=tic;
      testa=0;
      while complete ~= 0 && (length(ativas)>1)
          pos=1;
          vizinho=s_tupla;
          desliga=randi(length(ativas)); % posicao para desligar
          vizinho(ativas(desliga))=0; %desliga primeiro vizinho
          ativas(desliga)=0;
          ativas=(nonzeros(ativas))';%desabilita posicao para futuros sorteios
          %vizinho_base=vizinho; 
        while pos<=length(ativas)
            
          vizinho(ativas(pos))=0; % desliga segundo vizinho
          pos_desliga_2=ativas(pos);
          n_ativ_aux = n_ativas;
          
          while length(n_ativ_aux)>1 && (complete~=0)
              
                liga=randi(length(n_ativ_aux)); % posicao para ligar primeira opcao
                vizinho(n_ativ_aux(liga))=1; %liga primeiro vizinho
                pos_liga_1=n_ativ_aux(liga);
                n_ativ_aux(liga)=0;
                n_ativ_aux=(nonzeros(n_ativ_aux))';%desabilita posicao para futuros sorteios
                %aux=n_ativ_aux;%guardar vetor auxiliar com as posicoes disponiveis ainda
                pos2=1;
                
                while (iter_total<=200) && (pos2<=length(n_ativ_aux))
                    
                    vizinho(n_ativ_aux(pos2))=1; %liga segundo vizinho
                    key=num2str(vizinho);
                    if isKey(M,key)==0 % so testa se nao estiver na hash Table
                        cont=cont+1;
                        M(key)=vizinho;
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
                    end
                    
                      iter_total=iter_total+1;%atualiza o total de iteracoes
                      if iter_total>200 || ((toc(start_test)/60)>time_max) % limita esse loop a 5 minutos        
                         complete=0;
                         break;
                      end
                      vizinho(n_ativ_aux(pos2))=0; % reinicia
                      pos2=pos2+1;
                 end
                 %n_ativ_aux=aux;
                 vizinho(pos_liga_1)=0;%reinica os desligados
               
           end
           vizinho(pos_desliga_2)=1; % reinicia segundo vizinho desligado
           pos=pos+1;
           if complete==0
              break; 
           end
        end    
        
     end
  end    
    
 if cont==0
        s_viz_sol=s_tupla;
        s_viz_fobj=fobj_central;
        s_viz_crit=crit_central;
    else
       %% Retorna Melhor vizinho
        s_viz_sol=best_viz_tupla;
        s_viz_fobj=best_viz_fobj;
        s_viz_crit=best_viz_crit;  
        
  end

end