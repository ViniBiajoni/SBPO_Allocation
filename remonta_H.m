function [H]= remonta_H(solution,barras_atuais,UM,H,lote)
%teste= solution.*barras_atuais;
aux=(1:length(barras_atuais));
teste2=aux.*solution;
teste2=nonzeros(teste2);
% if length(teste2)==5
%     
%     disp('erro')
% end
%teste= nonzeros(teste)';
nm=0;

for i=1:lote
  
   H= [H;UM(teste2(i)).H_parc];
    
end
    
return
end

