function [fobj,num_crit]=  eval_solution(H_temp,H,w_k,nm)
format long
nm_2= length(H_temp(:,1))-length(H(:,1));%numero de medidas atualizado
G1=H_temp'*H_temp;
E=eye(nm_2+nm,nm_2+nm)-((H_temp)*(G1^-1)*H_temp');
dlmwrite('Covariancia.txt', E, 'delimiter', ' ', 'precision', '%.15f');
dlmwrite('nmed.txt',nm+nm_2);
dos('Crit_find_CPU.exe');
num_crit=(dlmread('Criticalidades.txt'))';

fobj=sum((w_k.*num_crit));

return

end

