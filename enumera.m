

function combs = enumera(combsInWave,lote,UMs_livres,combs_geradas,Cn)
    
    nUMs = UMs_livres;
    card = lote;
    combs = zeros(combsInWave,card);
    for linha = 1:combsInWave %da pra mudar para um limite desejado 
        nZ = nUMs - card;
        nO = card; 
        n = linha + combs_geradas;
        for i=1:nUMs
            nZ = nZ - 1;
            
            if (nO <= 0)||(nUMs - i <= 0)
                zcomb = 0;
            else
                zcomb = Cn(nUMs - i,nO);
            end

            if zcomb < n && nO > 0
                combs(linha,(card - nO)+1) = i; %pego a linha do combs
%                 disp(combs(linha,(card - nO)+1))
                nO = nO - 1;
                nZ = nZ + 1;
                n = n - zcomb;
            end
        end
         dlmwrite('combs.txt',combs(linha,:),'-append');
%        disp('-------------------------------------')
    end
    combs = sortrows(combs);
    disp(combs);
end

