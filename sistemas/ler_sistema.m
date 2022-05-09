function [Caso, Barra, Ramo]=ler_sistema(sis_file)
% Le arquivo de dados no formato IEEE Common Data Format.
%
% ler_sistema              -> Solicita o nome do arquivo de dados
% ler_sistema(cdf_file)    -> Utiliza o arquivo especificado
% 
% Define a estrutura dos dados das barra
Barra.num = [];
Barra.nome = {};
Barra.area = [];
Barra.zona = [];
Barra.tipo = [];
Barra.tensao = [];
Barra.angulo = [];
Barra.carga_MW = [];
Barra.carga_MVAR = [];
Barra.ger_MW = [];
Barra.ger_MVAR = [];
Barra.kVBase = [];
Barra.tensao_inic = [];
Barra.lim_max = [];
Barra.lim_min = [];
Barra.g = [];
Barra.b = [];
Barra.ctrl_rem = [];
Barra.P = [];
Barra.Q = [];
Barra.inj_P = [];
Barra.inj_Q = [];
Barra.estado_mod = [];
Barra.estado_ang = [];

% Define a estrutura dos dados dos ramos
Ramo.num = [];
Ramo.de = [];
Ramo.para = [];
Ramo.area = [];
Ramo.zona = [];
Ramo.circuito = [];
Ramo.tipo = [];
Ramo.r = [];
Ramo.x = [];
Ramo.b = [];
Ramo.rat_1 = [];
Ramo.rat_2 = [];
Ramo.rat_3 = [];
Ramo.barra_ctrl = [];
Ramo.lado = [];
Ramo.rel_tr = [];
Ramo.ang_tr = [];
Ramo.tap_min = [];
Ramo.tap_max = [];
Ramo.step = [];
Ramo.lim_min = [];
Ramo.lim_max = [];
Ramo.P_para = [];
Ramo.Q_para = [];
Ramo.P_de = [];
Ramo.Q_de = [];
Ramo.status = [];

% Define a estrutura do caso
if nargin < 1
    sis_file = input('Entre o nome do arquivo no formato IEEE CDF: ', 's');
end
% Abre o arquivo para leitura
fid = fopen(sis_file, 'r');
% Le e apresentar titulo do caso
while 1
    % Le a primeira linha do arquivo
    tline = fgetl(fid);
    % Se encontrar 'TAPE' le mais uma linha
    if findstr(tline,'TAPE');tline = fgetl(fid);end
    % Verifica se esta e a linha de titulo
    CASE_DATE=tline(2:9);
    if length(findstr(CASE_DATE,'/'))==2
        if length(tline) >= 30
            CASE_SENDER=tline(11:30);
        end
        if length(tline) >= 37
            MVA_BASE=str2double(tline(32:37));
        end
        if length(tline) >= 42
            YEAR=str2double(tline(39:42));
        end
        if length(tline) >= 44
            SEASON=tline(44);
        end
        if length(tline) >= 46
            CASE_IDENT=tline(46:length(tline));
        end
        break;
    end
end
% Procura a frase 'BUS DATA FOLLOWS'
while 1 
    tline = fgetl(fid);
    if  strcmp(tline(1:16), 'BUS DATA FOLLOWS'), break, end
end
n_ref = 0;
ibus = 0;
while 1
    tline = fgetl(fid);
    % Conclui se encontra o fim do bloco
    if  strcmp(tline(1:4), '-999'), break, end
    % Le dados das barras
    ibus = ibus + 1;
    Barra.num(ibus,1) = str2double(tline(1:4));
    Barra.nome(ibus,1) = cellstr(tline(6:17));
    Barra.area(ibus,1) = str2double(tline(19:20));
    Barra.zona(ibus,1) = str2double(tline(21:23));
    Barra.tipo(ibus,1) = str2double(tline(26));    % 25 - 26
    if Barra.tipo(ibus,1) == 3
        n_ref = n_ref + 1;
        if n_ref == 1
%             fprintf('Barra de Referencia: %d\n', Barra.num(ibus,1));
            Bref = Barra.num(ibus,1);
        else
%             fprintf('Barra de Referencia: %d\n', Barra.num(ibus,1));
            disp('AVISO: Mais de uma Barra de Referencia - Usando a primeira!')
        end
    end
    Barra.tensao(ibus,1) = str2double(tline(28:33));
    Barra.angulo(ibus,1) = pi/180*(str2double(tline(34:40)));
    if Barra.tipo(ibus,1) == 3
        if Barra.angulo(ibus,1) ~= 0
            fprintf('AVISO: Angulo da Barra de Referencia diferente de zero: %0.4f\n', Barra.angulo(ibus,1));
        end
    end
    Barra.carga_MW(ibus,1) = str2double(tline(41:49));
    Barra.carga_MVAR(ibus,1) = str2double(tline(50:58));
    Barra.ger_MW(ibus,1) = str2double(tline(59:67));
    Barra.ger_MVAR(ibus,1) = str2double(tline(68:75));
    Barra.kVBase(ibus,1) = str2double(tline(77:83));
    Barra.tensao_inic(ibus,1) = str2double(tline(85:90));
    Barra.lim_max(ibus,1) = str2double(tline(91:98));
    Barra.lim_min(ibus,1) = str2double(tline(99:106));
    Barra.g(ibus,1) = str2double(tline(107:114));
    Barra.b(ibus,1) = str2double(tline(115:122));
    Barra.ctrl_rem(ibus,1) = str2double(tline(124:127));
    Barra.estado_mod(ibus,1) = 1;
    Barra.estado_ang(ibus,1) = 0;
end
if n_ref == 0
    disp('AVISO: Nao ha Barra de Referencia! - Usando Barra 1')
    Bref = Barra.num(1,1);
end
% Procura a frase 'BRANCH DATA FOLLOWS'
while 1 
    tline = fgetl(fid);
    if  strcmp(tline(1:19), 'BRANCH DATA FOLLOWS'), break, end
end
% Ler dados dos ramos
ilin = 0;
while 1
    tline = fgetl(fid);
    % Conclui se encontra o fim do bloco
    if  strcmp(tline(1:4), '-999'), break, end
    ilin = ilin + 1;
    clin = length(tline);
    Ramo.num(ilin,1)=ilin;
    Ramo.de(ilin,1)=str2double(tline(1:4));
    Ramo.para(ilin,1)=str2double(tline(6:9));
    Ramo.area(ilin,1)=str2double(tline(11:12));
    Ramo.zona(ilin,1)=str2double(tline(13:15));    % 14
    Ramo.circuito(ilin,1)=str2double(tline(17));
    Ramo.tipo(ilin,1)=str2double(tline(19));
    Ramo.r(ilin,1)=str2double(tline(20:29));
    Ramo.x(ilin,1)=str2double(tline(30:39));            % 40
    Ramo.b(ilin,1)=str2double(tline(40:49));            % 41 - 50
    Ramo.rat_1(ilin,1)=str2double(tline(51:55));
    Ramo.rat_2(ilin,1)=str2double(tline(57:61));
    Ramo.rat_3(ilin,1)=str2double(tline(63:67));
    Ramo.barra_ctrl(ilin,1)=str2double(tline(69:72));
    Ramo.lado(ilin,1)=str2double(tline(74));
    Ramo.rel_tr(ilin,1)=str2double(tline(77:82));
    Ramo.ang_tr(ilin,1)=str2double(tline(84:90));
    Ramo.tap_min(ilin,1)=str2double(tline(91:97));
    Ramo.tap_max(ilin,1)=str2double(tline(98:104));
    Ramo.step(ilin,1)=str2double(tline(106:111));
    if clin <= 119
        Ramo.lim_min(ilin,1)=str2double(tline(113:clin));
    else
        Ramo.lim_min(ilin,1)=str2double(tline(113:119));
    end
    if clin <= 126
        Ramo.lim_max(ilin,1)=str2double(tline(120:clin));
    else
        Ramo.lim_max(ilin,1)=str2double(tline(120:126));
    end
    Ramo.status(ilin,1) = 1;
end
% Retorna os dados do caso
Caso.data = CASE_DATE;
Caso.origem = CASE_SENDER;
Caso.MVA = MVA_BASE;
Caso.ano = YEAR;
Caso.epoca = SEASON;
Caso.ident = CASE_IDENT;
Caso.NB = ibus;
Caso.NR = ilin;
Caso.Bref = Bref; 
return;

% DOCUMENT
% http://www.ee.washington.edu/research/pstca/formats/cdf.txt
% 
% Partial Description of the IEEE Common Data Format for the
% Exchange of Solved Load Flow Data
% The complete description can be found in the paper "Common Data
% Format for the Exchange of Solved Load Flow Data", Working Group on a
% Common Format for the Exchange of Solved Load Flow Data, _IEEE
% Transactions on Power Apparatus and Systems_, Vol. PAS-92, No. 6,
% November/December 1973, pp. 1916-1925.
% The data file has lines of up to 128 characters. The lines are grouped
% into sections with section headers. Data items are entered in specific
% columns. No blank items are allowed, enter zeros instead. Floating point
% items should have explicit decimal point. No implicit decimal points
% are used.
% Data type codes:
%       A - Alphanumeric (no special characters)
%       I - Integer
%       F - Floating point
%       * - Mandatory item
% Title Data
% ==========
% First card in file.
% CASE_DATE     Columns 2- 9    Date, in format DD/MM/YY with leading zeros.
%                               If no date provided, use 0b/0b/0b where b is
%                               blank.
% CASE_SENDER   Columns 11-30   Originator's name (A)
% MVA_BASE      Columns 32-37   MVA Base (F*)
% YEAR          Columns 39-42   Year (I)
% SEASON        Column  44      Season (S - Summer, W - Winter)
% CASE_IDENT    Column  46-73   Case identification (A)
% 
% Bus Data *
% ==========
% Section start card *:
% ---------------------
% Columns 1-16      BUS DATA FOLLOWS (not clear that any more than BUS in
%                   1-3 is significant) *
% Columns ?- ?      NNNNN ITEMS (column not clear, I would not count on this)
% 
% Bus data cards *:
% -----------------
% BUS_NO        Columns 1- 4    Bus number (I) *
% BUS_NAME      Columns 7-17    Name (A) (left justify) *
% BUS_AREA      Columns 19-20   Load flow area number (I) Don't use zero! *
% BUS_ZONE      Columns 21-23   Loss zone number (I)
% BUS_TYPE      Columns 25-26   Type (I) *
%                       0 - Unregulated (load, PQ)
%                       1 - Hold MVAR generation within voltage limits, (PQ)
%                       2 - Hold voltage within VAR limits (gen, PV)
%                       3 - Hold voltage and angle (swing, V-Theta) (must always
%                           have one)
% VOLT          Columns 28-33   Final voltage, p.u. (F) *
% ANGLE         Columns 34-40   Final angle, degrees (F) *
% L_MW          Columns 41-49   Load MW (F) *
% L_MVAR        Columns 50-59   Load MVAR (F) *
% G_MW          Columns 60-67   Generation MW (F) *
% G_MVAR        Columns 68-75   Generation MVAR (F) *
% KV_BASE       Columns 77-83   Base KV (F)
% D_VOLT        Columns 85-90   Desired volts (pu) (F) (This is desired remote
%                               voltage if this bus is controlling another bus.
% V_MAX         Columns 91-98   Maximum MVAR or voltage limit (F)
% V_MIN         Columns 99-106  Minimum MVAR or voltage limit (F)
% GPU           Columns 107-114 Shunt conductance G (per unit) (F) *
% BPU           Columns 115-122 Shunt susceptance B (per unit) (F) *
% R_CTRL        Columns 124-127 Remote controlled bus number
%
% Section end card:
% -----------------
% Columns 1- 4      -999
% 
% Branch Data *
% =============
% Section start card *:
% ---------------------
% Columns 1-16      BRANCH DATA FOLLOWS (not clear that any more than BRANCH
%                   is significant) *
% Columns 40?- ?    NNNNN ITEMS (column not clear, I would not count on this)
% 
% Branch data cards *:
% --------------------
% BRANCH_FROM   Columns 1- 4    Tap bus number (I) *
%                               For transformers or phase shifters, the side
%                               of the model the non-unity tap is on
% BRANCH_TO     Columns 6- 9    Z bus number (I) *
%                               For transformers and phase shifters, the side
%                               of the model the device impedance is on.
% BRANCH_AREA   Columns 11-12   Load flow area (I)
% BRANCH_ZONE   Columns 13-14   Loss zone (I)
% CIRCUIT_NO    Column 17       Circuit (I) * (Use 1 for single lines)
% BRANCH_TYPE   Column 19       Type (I) *
%                               0 - Transmission line
%                               1 - Fixed tap
%                               2 - Variable tap for voltage control (TCUL, LTC)
%                               3 - Variable tap (turns ratio) for MVAR control
%                               4 - Variable phase angle for MW control (phase shifter)
% RPU           Columns 20-29   Branch resistance R, per unit (F) *
% XPU           Columns 30-40   Branch reactance X, per unit (F) * No zero impedance lines
% BPU           Columns 41-50   Line charging B, per unit (F) * (total line charging, +B)
% RATING_1      Columns 51-55     Line MVA rating No 1 (I) Left justify!
% RATING_2      Columns 57-61     Line MVA rating No 2 (I) Left justify!
% RATING_3      Columns 63-67     Line MVA rating No 3 (I) Left justify!
% CTRL_BUS      Columns 69-72     Control bus number
% SIDE          Column 74         Side (I)
%                       0 - Controlled bus is one of the terminals
%                       1 - Controlled bus is near the tap side
%                       2 - Controlled bus is near the impedance side (Z bus)
% XTR_RATIO     Columns 77-82     Transformer final turns ratio (F)
% XTR_ANGLE     Columns 84-90     Transformer (phase shifter) final angle (F)
% TAP_MIN       Columns 91-97     Minimum tap or phase shift (F)
% TAP_MAX       Columns 98-104    Maximum tap or phase shift (F)
% STEP          Columns 106-111   Step size (F)
% LIMIT_MIN     Columns 113-119   Minimum voltage, MVAR or MW limit (F)
% LIMIT_MAX     Columns 120-126   Maximum voltage, MVAR or MW limit (F)
%
% Section end card:
% -----------------
% Columns 1- 4      -999
%
% Loss Zone Data
% ==============
% Section start card
% ------------------
% Columns 1-16      LOSS ZONES FOLLOWS (not clear that any more than LOSS
%                   is significant)
% Columns 40?- ?    NNNNN ITEMS (column not clear, I would not count on this)
%
% Loss Zone Cards:
% ----------------
% Columns 1- 3      Loss zone number (I)
% Columns 5-16      Loss zone name (A)
%
% Section end card:
% -----------------
% Columns 1- 3      -99
% 
% Interchange Data *
% ==================
% Section start card
% ------------------
% Columns 1-16      INTERCHANGE DATA FOLLOWS (not clear that any more than
%                   first word is significant).
% Columns 40?- ?    NNNNN ITEMS (column not clear, I would not count on this)
%
% Interchange Data Cards *:
% -------------------------
% Columns 1- 2      Area number (I) no zeros! *
% Columns 4- 7      Interchange slack bus number (I) *
% Columns 9-20      Alternate swing bus name (A)
% Columns 21-28     Area interchange export, MW (F) (+ = out) *
% Columns 30-35     Area interchange tolerance, MW (F) *
% Columns 38-43     Area code (abbreviated name) (A) *
% Columns 46-75     Area name (A)
%
% Section end card:
% -----------------
% Columns 1- 2      -9
%
% Tie Line Data
% =============
% Section start card
% ------------------
% Columns 1-16      TIE LINES FOLLOW (not clear that any more than TIE
%                   is significant)
% Columns 40?- ?    NNNNN ITEMS (column not clear, I would not count on this)
%
% Tie Line Cards:
% ---------------
% Columns 1- 4      Metered bus number (I)
% Columns 7-8       Metered area number (I)
% Columns 11-14     Non-metered bus number (I)
% Columns 17-18     Non-metered area number (I)
% Column 21         Circuit number
%
% Section end card:
% -----------------
% Columns 1- 3      -999

