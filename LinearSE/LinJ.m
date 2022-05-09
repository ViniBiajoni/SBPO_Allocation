function [H,order]=LinJ(medidores,Ramo,Barra,Caso) % calcula o Jacobiano H

adjacente=[];
nb=Caso.NB;
nm=max(size(medidores.num));


phasemeas=find(medidores.tipo==3,1);
if(isempty(phasemeas))
    ns=nb-1;
    reference=Barra.num(1);
else
    ns=nb;
    reference=-1;
end
% encontrar ramos incidentes
for i=1:nb
    idx=logical( Barra.num(i)==Ramo.de | Barra.num(i)==Ramo.para ); %seleciona ramos adjacentes
    de=Ramo.de(idx);
    para=Ramo.para(idx);
    k=0;
    for j=1:size(de,1)
        k=k+1;
        if (de(j) == Barra.num(i)) % se barra de então para eh adjacente
            adjacente{i}(k)=para(j);
        else
            adjacente{i}(k)=de(j); % se barra para entao de eh adjacente
        end
    end
end



H=zeros(nm,ns);
Aux=1:nm;
%% medidas de angulo
idx=logical(medidores.tipo==3);
if ~isempty(find(idx,1))
    angulo.num=medidores.num(idx);
    angulo.leitura=medidores.leitura(idx);
    angulo.para=medidores.para(idx);
    nma=max(size(angulo.num));
    order(1:nma)=Aux(idx);
    for l=1:nma
        i=angulo.para(l);
        H(l,i)=1;
    end
else
    nma=0;
end

%%  medidas de injecao de potencia
idx=logical(medidores.tipo==2);
if ~isempty(find(idx,1))
    injecoes.num=medidores.num(idx);
    injecoes.leitura=medidores.leitura(idx);
    injecoes.para=medidores.para(idx);
    nmp=max(size(injecoes.num));
    order(nma+1:nma+nmp)=Aux(idx);
    for l=1:nmp
        i=injecoes.para(l);
        H(l+nma,adjacente{i}(:))=-1; % barras adjacentes recebem -1
        H(l+nma,i)=-sum(H(l+nma,adjacente{i}(:))); % barras com a injecção recebem o oposto da soma
    end
else
    nmp=0;
end
%% medidas de fluxo de potencia
idx=logical(medidores.tipo==1);
if ~isempty(find(idx,1))
    fluxos.num=medidores.num(idx);
    fluxos.leitura=medidores.leitura(idx);
    fluxos.para=medidores.para(idx);
    fluxos.de=medidores.de(idx);
    nmpij=max(size(fluxos.num));
    order(nma+nmp+1:nma+nmp+nmpij)=Aux(idx);

    for l=1:nmpij
        de=fluxos.de(l);
        para=fluxos.para(l);
        H(l+nma+nmp,de)=+1;
        H(l+nma+nmp,para)=-1;
    end
else
    nmpij=0;
end
%% medidas de Injeçao de corrente
idx=logical(medidores.tipo==9);
if ~isempty(find(idx,1))
    injecoes.num=medidores.num(idx);
    injecoes.leitura=medidores.leitura(idx);
    injecoes.para=medidores.para(idx);
    nmii=max(size(injecoes.num));
    order(nma+nmp+nmpij+1:nma+nmp+nmpij+nmii)=Aux(idx);
    for l=1:nmii
        i=injecoes.para(l);
        H(l+nma+nmp+nmpij,adjacente{i}(:))=-1; % barras adjacentes recebem -1
        H(l+nma+nmp+nmpij,i)=-sum(H(l+nma+nmp+nmpij,adjacente{i}(:))); % barras com a injecção recebem o oposto da soma
    end
else
    nmii=0;
end
%% medidas de corrente nos ramos
idx=logical(medidores.tipo==7);
if ~isempty(find(idx,1))
    fluxos.num=medidores.num(idx);
    fluxos.leitura=medidores.leitura(idx);
    fluxos.para=medidores.para(idx);
    fluxos.de=medidores.de(idx);
    nmiij=max(size(fluxos.num));
    order(nma+nmp+nmpij+nmii+1:nma+nmp+nmpij+nmii+nmiij)=Aux(idx);

    for l=1:nmiij
        de=fluxos.de(l);
        para=fluxos.para(l);
        H(l+nma+nmp+nmpij+nmii,de)=+1;
        H(l+nma+nmp+nmpij+nmii,para)=-1;
    end
%else
    %nmiij=0;
end
%%
ref_idx=logical(Barra.num~=reference);

H=H(:,ref_idx');



save Jacobiano H medidores

return