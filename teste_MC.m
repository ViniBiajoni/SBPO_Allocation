
%rng('default') % For reproducibility
p_fail=zeros(1,5)+0.1;
sum_fvalue=0;
N=1000000;

for num_casos=1:N
    x=zeros(1,5);
    for i=1:5
     r = rand;
     if (abs(r)>p_fail(i))
         x(i)=1;
     end
    end
    
    sum_fvalue= sum_fvalue + (1 - (1-x(1)*x(4))*(1-x(2)*x(5))*(1-x(1)*x(3)*x(5))*(1-x(2)*x(3)*x(4)));

end% teste Monte Carlo Basico



evaluation= sum_fvalue/N;
disp(evaluation);