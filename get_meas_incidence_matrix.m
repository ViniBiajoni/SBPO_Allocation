function [M]=get_meas_incidence_matrix(medidores,nm,Ramo,Caso)
M=zeros(nm,Caso.NR);
for imeas=1:nm
    if medidores.de(imeas)~=0 && medidores.para(imeas)~=0 %medida de fluxo

        iramo=find(logical(medidores.de(imeas)==Ramo.de & medidores.para(imeas)==Ramo.para));
        if ~isempty(iramo)
            M(imeas,iramo)=1;
            continue;
        end
        iramo=find(logical(medidores.de(imeas)==Ramo.para & medidores.para(imeas)==Ramo.de));
        if ~isempty(iramo)
            M(imeas,iramo)=- 1;
            continue;
        end

    else % medida de injeção
        iramo=find(medidores.para(imeas)==Ramo.de);
        if ~isempty(iramo)
            M(imeas,iramo)=1;
        end
        iramo=find(medidores.para(imeas)==Ramo.para);
        if ~isempty(iramo)
            M(imeas,iramo)=- 1;
        end
    end
end
return;
end