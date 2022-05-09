function [A]=get_net_incidence_matrix(Ramo,Caso)
A=zeros(Caso.NR,Caso.NB);
for ilin=1:Caso.NR    
    A(ilin,Ramo.de(ilin))=+1;
    A(ilin,Ramo.para(ilin))=-1;
end
return;
end