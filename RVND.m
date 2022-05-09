function [sol_VND,fobj_VND,crit_VND] =RVND(actual_solution,actual_fobj,actual_crit,NB)
%global M
global start_test
global time_max
kmax= 4;% numero maximo de vizinhancas
index=0;
rand_numbers=(1:4);
k=zeros(1:4);
s_tupla= actual_solution;
s_fobj= actual_fobj;
s_crit = actual_crit;
options=4;

% %% Inicia Tabela Hash
% M = containers.Map('KeyType','char','ValueType','any'); 

%% Sorteio a primeira vizinhanca
sorteio=randi(options);
viz=rand_numbers(sorteio);
rand_numbers(sorteio)=0;
rand_numbers=(nonzeros(rand_numbers))';
index=index+1;
k(index)=viz;
options=options-1;

while index<=kmax && (toc(start_test)/60)<time_max
    
    if (k(index)==1)
       %[s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz1(s_tupla,s_fobj,s_crit);
       %[s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz1_altern(s_tupla);
       [s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz1_altern_hash_t(s_tupla);
    end
    
    if (k(index)==2)
       %[s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz2(s_tupla,s_fobj,s_crit);
       %[s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz2_altern(s_tupla);
       [s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz2_altern_hash_t(s_tupla);
        
    end
    
    if (k(index)==3)
       %[s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz3(s_tupla,s_fobj,s_crit); 
       %[s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz3_altern(s_tupla);
       [s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz3_altern_hash_t(s_tupla);
    end
    
    if (k(index)==4)
       %[s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz4(s_tupla,s_fobj,s_crit);
       %[s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz4_altern(s_tupla); 
       [s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz4_altern_hash_t(s_tupla);
    end
    
%     
%     if (k(index)==5)
%        [s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz5(s_tupla,s_fobj,s_crit); 
%     end
    
    %% Avalia resultado da vizinhanca atual
    if s_viz_fobj < s_fobj
        s_tupla= s_viz_sol;
        s_fobj= s_viz_fobj;
        s_crit = s_viz_crit;  
        
        index=1; % retorna à vizinhanca inicial
    else
        
        if options>0
           sorteio=randi(options);
           viz=rand_numbers(sorteio);
           rand_numbers(sorteio)=0;
           rand_numbers=(nonzeros(rand_numbers))';
           index=index+1;
           k(index)=viz;
           options=options-1;
        else
           index=index+1;
        end
         
    end
    
end

%% Retorna Solucao
    sol_VND=s_tupla;
    fobj_VND=s_fobj;
    crit_VND=s_crit;

end
