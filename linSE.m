function [output_args]=linSE(input_args)
% linSE Linear State Estiomation algorithm
% measurement model: z=Hx+e
% x'=inv(HtH)Htz
% allow bypassing of reading statements

%clc

tol=1e-10;

if nargin < 1

    %sis_file=input('Arquivo da rede','s');
    %med_file=input('Arquivo de medidores','s');
    disp('Nenhum argumento de entrada');

    disp(nargin);

    %med_file='.\medidas\6BusCase01.med';

    %sis_file='.\sistemas\6Bus.cdf';

    %med_file='.\medidas\ieee14_london_teste.med';

    %med_file='.\medidas\ieee14_london_teste.med';

    %sis_file='.\sistemas\ieee14.cdf';

    %% IEEE 24 bus TEST SYSTEM
    sis_file='.\sistemas\ieee24.cdf';

    med_file='.\medidas\ieee24_gen_meet_2013.med';

    [Caso, Barra, Ramo]=ler_sistema(sis_file);

    [medidores]=ler_medidores(med_file);

    idz=logical( medidores.ok==0 & ( medidores.tipo==1 | medidores.tipo==2 | medidores.tipo==3) ); % medidores ativos de fluxo e injecao

    medidores.num = medidores.num(idz);

    medidores.de = medidores.de(idz);

    medidores.para = medidores.para(idz);

    medidores.tipo = medidores.tipo(idz);

    medidores.PMU_num = medidores.PMU_num(idz);

    medidores.leitura = medidores.leitura(idz);

    [H , order] = LinJ(medidores,Ramo,Barra,Caso);

    medidores.num = medidores.num(order);

    medidores.de = medidores.de(order);

    medidores.para = medidores.para(order);

    medidores.tipo = medidores.tipo(order);

    medidores.PMU_num = medidores.PMU_num(order);

    medidores.leitura = medidores.leitura(order);

    z=medidores.leitura;

    output_args{7}=order;


else

    disp('Mais de um argumento de entrada');

    disp(nargin);

    H=input_args{2};

    medidores=input_args{1};

    z=medidores.leitura;

end

m=size(H,1);


if (m == size(z,2))
    z=z';
end


I=eye(m);

G=(H'*H);

detG=det(G);

fprintf('\n O determinante de G eh: %5.3f \n',detG);

%pause

invG=inv(G);

x=invG*H'* z;

E=I-H*invG*H';

r=E*z;

% get critical measurements
cmeas_idx=logical(abs(diag(E)) < tol & abs(r) < tol);

if ~isempty(find(cmeas_idx,1))

    Ccrit=find(cmeas_idx)';

    cm=size(Ccrit,2);

    cmeas_idx_c=~cmeas_idx;

    rn=zeros(size(r));

    %get non critical residual sensitivity matrix
    E_noncrit=E(cmeas_idx_c,cmeas_idx_c);

    %get residuals of non-critical measurements
    rn(cmeas_idx_c)=abs(r(cmeas_idx_c)./sqrt(diag(E_noncrit)));
    %get non critical measurements
    %z_noncrit=z(cmeas_idx_c);

    %Critical measurements identification

    Cmeas=[];
    strCmeas=[];
    crit_size=size(Ccrit,2);

    for i=1:crit_size

        num=Ccrit(i);

        Cmeas{i}=medidores.num(num);

        switch medidores.tipo(num)

            case 1 % active power flow measurements

                strCmeas{i}={strcat('P',int2str(medidores.de(num)),'-',int2str(medidores.para(num)))};

            case 2 % active power injections measurements

                strCmeas{i}={strcat('P',int2str(medidores.para(num)))};

            case 3 % bus angle measurements

                strCmeas{i}={strcat('A',int2str(medidores.para(num)))};

        end

    end

else

    rn=abs(r./sqrt(diag(E)));

    cm=0;

    Cmeas=[];
    strCmeas=[];

end

%sort non critical measurements
[List,order]=sort(rn,'descend'); % zero measurements correspondig to critical measurements are placed in the end of the sorted list

save State r x E

Cset=[];

Cset{1}(1)=order(1);

i=2;

j=1;

k=1;

% get candidate Critical measurement sets

% Assumption: List has the form [ a a  b c d e e e f g h i j ]

while 1

    comp=round(List(i)-List(i-1));

    pij=abs(E(order(i),order(i-1)))/sqrt(E(order(i),order(i))*E(order(i-1),order(i-1)));

    if ( abs(pij - 1) < tol && comp < tol)

        j=j+1;

        Cset{k}(j)=order(i);

    else

        k=k+1;

        j=1;

        Cset{k}(j)=order(i);

    end

    if i < (m-cm)

        i=i+1;

    else

        break;

    end

end

Csize=k; %# of Critical Sets
% remove Cset with only one element

i=1;

while i < Csize

    if (size(Cset{i},2) < 2)

        for j=i+1:Csize

            Cset{j-1}=Cset{j};

        end

        Cset(j)=[];

        Csize=Csize-1;

    else

        i=i+1;

    end

end

if (size(Cset{i},2) < 2)

    Cset(i)=[];

    Csize=Csize-1;
end

Csize;

%Critical sets identification

Csets=[];
strCsets=[];

for k=1:Csize

    setsize=size(Cset{k},2);

    for i=1:setsize

        num=Cset{k}(i);

        Csets{k}(i)=medidores.num(num);

        switch medidores.tipo(num)

            case 1 % active power flow measurements

                strCsets{k}(i)={strcat('P',int2str(medidores.de(num)),'-',int2str(medidores.para(num)))};

            case 2 % active power injections measurements

                strCsets{k}(i)={strcat('P',int2str(medidores.para(num)))};

            case 3 % bus angle measurements

                strCsets{k}(i)={strcat('A',int2str(medidores.para(num)))};

        end

    end

end

output_args{1}=x;

output_args{2}=r;

output_args{3}=rn;

output_args{4}=E;

output_args{5}=medidores;

output_args{6}=H;

output_args{7}=Cmeas;

output_args{8}=Csets;

output_args{9}=strCmeas;

output_args{10}=strCsets;

% inserir strCset and strCset - OK!

end

% Convençoes:
% - Medidas de fluxo de potencia ativa:     tipo 1
% - Medidas de injecao de potencia ativa:   tipo 2
% - Medidas de angulo:                      tipo 3
% - Medidas de fluxo de potencia reativa:   tipo 4
% - Medidas de injecao de potencia reativa: tipo 5
% - Medidas de tensao:                      tipo 6
% - Medidas de corrente de ramo real        tipo 7
% - Medidas de corrente de ramo imag        tipo 8
% - Medidas de injecao de corrente real     tipo 9
% - Medidas de injecao de corrente imag     tipo 10
% - Medidas de fluxo de 1 => 2   ---> de=1 e para=2
% - Medidas de tensao na barra 3 ---> de=0 e para=3