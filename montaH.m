function UM = montaH(UM,MAdj,free_bus,H,NB)

num_estados= length(H(1,:));
pmu=0;

    for j=1:length(free_bus)    
        inj= sum(MAdj(free_bus(j),:));
        num_estados2=0;
        num_flux= inj;
        if (num_estados==NB)
            H_parc=zeros(num_flux+1,num_estados);
            pmu=1;
        else
            H_parc=zeros(num_flux+1,num_estados+1);
            num_estados2=num_estados+1;
        end
        cont=2;
        H_parc(1,free_bus(j))=inj;
        for k=1:num_estados2
            if (MAdj(free_bus(j),k)==1)
                H_parc(1,k)= -1;
                H_parc(cont,free_bus(j))=1;
                H_parc(cont,k)=-1;
                cont=cont+1;
            end
        end
        if (pmu==0)
            
            H_parc=H_parc(:,2:num_estados2);

        end
        UM(j).H_parc= H_parc;
        UM(j).Barra= free_bus(j);
        UM(j).Num_Medidas= length(H_parc(:,1));
    end
        
     return  
end

