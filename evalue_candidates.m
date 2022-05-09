function  [melhor_atual,custo_atual,removed]= evalue_candidates1(D,C,list_compras,lojas_disp,origem)

c= zeros(1,length(lojas_disp)); %custo para ir ao candidato
lojas= zeros(1,length(lojas_disp)); % lojas associadas 
removed=[]; %lojas que nao temmais produtos validos
    for i=1:length(lojas_disp)
        
        lojas(i)= lojas_disp(i);
        
        if (min(C(lojas_disp(i)+1,:)) == Inf )
       
            c(i)= Inf;
        
        else
       
            for j=1:length(C(lojas_disp(i)+1,:))
                
                if ((list_compras(i)*C(lojas_disp(i)+1,j)) ~= Inf)
                    
                    c(i)= c(i) + list_compras(j)*C(lojas_disp(i)+1,j);
                
                end
                
            end
            
            if (c(i)==0)
                c(i)=Inf;
                removed=[removed i];
            else
                c(i)= c(i) + D(origem+1,lojas_disp(i)+1); % adiciona a distancia (custo da distancia)
            end
        end
    
    end
     
    [minimo,pos]= min(c);
    melhor_atual= lojas(pos); %melhor loja
    custo_atual= minimo;
    
    return

end
