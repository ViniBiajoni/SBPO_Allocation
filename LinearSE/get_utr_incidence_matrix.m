function [AMU,mu,nmu]=get_utr_incidence_matrix(nm,medidores,Caso)
AMU=zeros(nm,Caso.NB);
mu=false(Caso.NB,1);
nmu=Caso.NB;
for imeas=1:nm
    if medidores.de(imeas)~=0 && medidores.para(imeas)~=0 %medida de fluxo
    AMU(imeas,medidores.de(imeas))=1;
    else
        AMU(imeas,medidores.para(imeas))=1;
    end   
end
for ibus=1:Caso.NB
    if nnz(AMU(:,ibus)) > 0
        mu(ibus)=true; % UTR at bus ibus has measurements
    end
return;
end