function vizinho= troca_2_posicoes(vizinho,ativas,n_ativas)
        
        desliga=randi(length(ativas)); % posicao para desligar
        liga=randi(length(n_ativas)); % posicao para ligar
        vizinho(ativas(desliga))=0;
        vizinho(n_ativas(liga))=1;
        ativas(desliga)=0;
        n_ativas(liga)=0;
        
        %%  Atualiza as posicoes disponiveis
        ativas=(nonzeros(ativas))';
        n_ativas=(nonzeros(n_ativas))';
        
        %% Sorteia as segundas opções
        desliga=randi(length(ativas)); % posicao para desligar
        liga=randi(length(n_ativas)); % posicao para ligar
        vizinho(ativas(desliga))=0;
        vizinho(n_ativas(liga))=1;

end

