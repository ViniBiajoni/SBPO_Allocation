function [Barras_Livres,UMs_Livres]= UMs_complementares(mu_location,NB)

mu_location= nonzeros(mu_location);
Barras_Livres=[];
UMs_Livres=0;

for i=1:NB
   
    pos=find(mu_location==i);
    
    if (isempty(pos))
        Barras_Livres=[Barras_Livres i];
        UMs_Livres=UMs_Livres+1;
    end
    
end

return

end

