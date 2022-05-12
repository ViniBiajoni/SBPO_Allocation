function barras_sol = print_UM_Sol(max_testes,best_sol_testes,Barras_Livres,lote)
barras_sol=zeros(max_testes,lote);
for k=1:max_testes
    j=1;
    for i=1:length(Barras_Livres)
        if(best_sol_testes(k,i)==1)
           barras_sol(k,j) = Barras_Livres(i);
           j=j+1;
        end
    end
end
end

