% solucao inicial
 % Pegando as UMs com mais medidas disponíveis
%  solution=zeros(1,UMs_Livres);
%  for i=1:UMs_Livres
%     temp(i)= i; 
%     barras(i)= UM(i).Barra;
%     meds(i)= UM(i).Num_Medidas;
%  end
%  a=[temp;barras;meds];
%  [idx,order]=sort(a(3,:),'descend');
%  temp2=a(1,order);
% for i=1:lote
%    
%     solution(temp2(1,i))=1;
%     
% end
%  disp(solution);
%  best_fobj=0;
%  best_crit=[];
%  [H_temp]= remonta_H(solution,barras_atuais,UM,H,lote);
%  [fobj_atual,num_crit]= eval_solution(H_temp,H,w_k,nm);
%  best_fobj=fobj_atual;
%  best_crit=num_crit;
%  
 
%% Aleatoria
 %pop=[];
 %pop_max=0; 
 %while (pop_max<100)
 % Aleatoria
%     complete =0;
%     solution=zeros(1,UMs_Livres);
% 
%     while (complete ==0)
%         sorteio=randi(UMs_Livres);
%         if (solution(sorteio)==0)
%             solution(sorteio)=1;
%         end
% 
%         if (sum(solution)==lote)
%             complete=1;
%             %pop_max=pop_max+1;
%             %pop=[pop;solution];
%         end
%     end
%  %end
%  
% best_fobj=0;
% best_crit=[];
% [H_temp]= remonta_H(solution,barras_atuais,UM,H,lote);
% [fobj_atual,num_crit]= eval_solution(H_temp,H,w_k,nm);
% best_fobj=fobj_atual;
% best_crit=num_crit;
% best_solution=solution;

%% Testa fobj para Multiplas Solucoes
% fobj=[];
% crit=[];
% ordenacao=(1:pop_max);
% tic()
%  for i=1:pop_max
%      solution=pop(i,:);
%     [H_temp]= remonta_H(solution,barras_atuais,UM,H,lote);
%     [fobj_atual,num_crit]= eval_solution(H_temp,H,w_k,nm);
%     fobj=[fobj;fobj_atual];
%     crit=[crit;num_crit];
%     
%  end
% toc()
% a=[fobj';ordenacao];
% [idx,order]=sort(a(1,:));
% ordenacao= a(2,order);  % indexado das melhores para as piores solucoes