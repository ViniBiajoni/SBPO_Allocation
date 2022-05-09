function [sol_VND,fobj_VND,crit_VND] =VND(actual_solution,actual_fobj,actual_crit,NB)
%global M
global start_test
global time_max
kmax= 4;% numero maximo de vizinhancas
k=1;
s_tupla= actual_solution;
s_fobj= actual_fobj;
s_crit = actual_crit;


while k<=kmax && (toc(start_test)/60)<time_max
    
    if (k==1)
       %[s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz1(s_tupla,s_fobj,s_crit);
       %% Testes
       %[s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz1_altern(s_tupla);
       [s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz1_altern_hash_t(s_tupla);
    end
    
    if (k==2)
       %[s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz2(s_tupla,s_fobj,s_crit);
       %% Testes
       %[s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz2_altern(s_tupla);
       [s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz2_altern_hash_t(s_tupla);
    end
    
    if (k==3)
       %[s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz3(s_tupla,s_fobj,s_crit); 
       %% Testes
       %[s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz3_altern(s_tupla);
       [s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz3_altern_hash_t(s_tupla);
    end
    
    if (k==4)
       %[s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz4(s_tupla,s_fobj,s_crit);
       %% Testes
       %[s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz4_altern(s_tupla);
       [s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz4_altern_hash_t(s_tupla);
    end
    

%     if (k==5)
%        [s_viz_sol,s_viz_fobj,s_viz_crit]= est_viz5(s_tupla,s_fobj,s_crit); 
%     end
    
    %% Avalia resultado da vizinhanca atual
    if s_viz_fobj < s_fobj
        s_tupla= s_viz_sol;
        s_fobj= s_viz_fobj;
        s_crit = s_viz_crit;  
        
        k=1; % retorna à vizinhanca inicial
    else
        
        k=k+1;
    end
    
end

%% Retorna Solucao
    sol_VND=s_tupla;
    fobj_VND=s_fobj;
    crit_VND=s_crit;

end

