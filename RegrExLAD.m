function out = RegrExLAD(GEM,D,Options)
       
%************************RegrExLAD Implementation*****************************

%Depends on:
%
% Gurobi, which can be found in www.gurobi.com (free academic licenses 
% available)
% contextmodel2COBRA function (provided in supplementary material)
% FVA function (provided in supplementary material)
%
%Arguments:
%
%     Required
%
% GEM: Metabolic model in COBRA structure. Blocked reactions should be 
% removed first, for instance using the reduceModel function of the COBRA 
% toolbox
% 
% D: Experimental data vector (values already mapped to reactions) must 
% be of same length as the number of reactions. It can be provided by the 
% genetorxn function
% 
%     Options:
%
% Options must be a structure with the following fields (not all required)
%
%     Options.Dstd: a vector containing the standard deviation of data, it is 
%     intended to be used as weighting factor for reactions with associated  
%     data, in the form (1/CV(D))*||V-D||_1, where CV=std(D)/mean(D)
% 
%     Options.essmet: A vector indicating the indexes of the metabolites that should be 
%     present in the final model
% 
%     Options.tasks: A vector indicating the indexes of the reactions that should be
%     active in the final model (i.e. flux value above epsilon)
% 
%     Options.blocked: A vector indicating the indexes of the reactions that should not
%     be active in the final model (i.e. zero flux value)
% 
%     Options.L_min: Minimum lambda value in the lambda sequence (default is set to 0)
% 
%     Options.L_max: Maximum lambda value in the lambda sequence (default is set to .2)
% 
%     Options.L_int: Step size used to generate the lambda sequence (default is set to
%     0.1)
%
%     Options.v_max: the maximum flux capacity of the reactions, default
%     set to 1 (the data vector must be scaled accordingly)
% 
%     Options.rxns_lb: a vector containing the lower bounds for each reaction (default
%     is set to 0, min value is 0, reversible reactions are splitted)
% 
%     Options.rxns_ub: a vector containing the upper bounds for each reaction (default 
%     is set to v_max)
% 
%     Options.time_limit: scalar indicating the time limit in seconds for the MIQP
%     (default is set to 60)
% 
%     Options.out_flag: selects if the gurobi solution's flag should be
%     included in the console, must be 0 or 1 (default 0)
% 
%     Options.v_eps: scalar indicating the threshold value to consider a reaction active,
%     i.e. a reaction i with flux value vi is active if vi>=epsilon
% 
%     Options.E_max: vector of length equal to the number of reactions indicating the 
%     maximum allowed error value between data and flux prediction, i.e. E=d-v,
%     (default is set to 1000)
% 
%     Options.contextname: character string indicating the name of the context to
%     display in the COBRA structure (default value is empty)
% 
%     Value:
%
%    Out: A structure with the fields,
% 
%         Flux: flux distribution
%         Cor: correlation value between fluxes and data
%         ZE: total error value (residual)
%         ActRxn: set of active (i.e. |Vi|>epsilon) reactions
%         ActRxnData: set of active reactions with associated data
%         Card: Number of active reactions
%         CardData: Number of active reactions with associated data
%         Lambda: used lambda-sequence
%
%    ContextCOBRA: A structure in COBRA format (see contextModel2COBRA
%    function for details
%
%**************************************************************************
%         Semidán (robaina@mpimp-golm.mpg.de), May, 2016
%**************************************************************************


%Argument evaluation
if ~exist('GEM','var'),
    error('The Genome-Scale model is missing...')
end

if ~exist('D','var'),
    error('The data vector is missing...')
end

if ~exist('Options','var'),
    Options = struct;
end

if ~isfield(GEM,'rev'),
    warning('rev (reversible reactions) field missing from GEM, constructing one from flux bounds')
    GEM.rev = zeros(length(GEM.rxns),1);
    GEM.rev(GEM.lb < 0) = 1;
end

if ~isfield(Options,'Dstd'),
    Dstd = 1;
else
    Dstd = Options.Dstd;
end

if ~isfield(Options,'essmet'), 
    essmet = [];
else
    essmet = Options.essmet;
end

if ~isfield(Options,'tasks'), 
    tasks = [];
else
    tasks = Options.tasks;
end

if ~isfield(Options,'blocked'), 
    blocked = [];
else
    blocked = Options.blocked;
end

if ~isfield(Options,'L_min'), 
    L_min = 0;
else
    L_min = Options.L_min;
end

if ~isfield(Options,'L_max'), 
    L_max = 0.1;
else
    L_max = Options.L_max;
end

if ~isfield(Options,'L_int'), 
    L_int = 0.01;
else
    L_int = Options.L_int;
end

if ~isfield(Options,'v_max'), 
    v_max = 1;
else
    v_max = Options.v_max;
end

if ~isfield(Options,'rxns_lb'), 
    rxns_lb = 0;
else
    rxns_lb = Options.rxns_lb;
end

if ~isfield(Options,'rxns_ub'), 
    rxns_ub = v_max;
else
    rxns_ub = Options.rxns_ub;
end

if ~isfield(Options,'MIPGap'), 
    MIPGap = 0.005;
else
    MIPGap = Options.MIPGap;
end

if ~isfield(Options,'OutFlag'), 
    OutFlag = 0;
else
    OutFlag = Options.OutFlag;
end

if ~isfield(Options,'v_eps'), 
    v_eps = v_max*1e-6;
else
    v_eps = Options.v_eps;
end

if ~isfield(Options,'E_max'), 
    E_max = v_max*1e3;
else
    E_max = Options.E_max;
end

if ~isfield(Options,'contextname'), 
    contextname = [];
else
    contextname = Options.contextname;
end

%Parse Data
S = GEM.S;
RevRxns = find(GEM.rev == 1);
IrrRxns = (setdiff(1:size(S,2),RevRxns))';

for i = 1:length(D),
    if length(Dstd) > 1,
        if isnan(Dstd(i)) || isinf(Dstd(i)),
           Dstd(i) = 1;
           warning('Error: Data std point Dstd(%d) is NaN or Inf and has been converted to 1',i);
        end
    end
end
if max(D) > v_max,
   D = D / max(D);
   Dstd = Dstd / max(D);
end
D = v_max*D;
D_or = D;

%Reaction partition: Irr_D, Irr_nD, Bfordat, Bfornodat, Rev_D, Rev_nD
Irr_D = IrrRxns(~isnan(D(IrrRxns)));
Irr_nD = IrrRxns(isnan(D(IrrRxns)));
Rev_D = RevRxns(~isnan(D(RevRxns)));
Rev_nD =  RevRxns(isnan(D(RevRxns)));
NIrr_D = length(Irr_D);
NRev_D = length(Rev_D);
NRev_nD = length(Rev_nD);
Rxns_Or = [Irr_D;Irr_nD;Rev_D;Rev_nD;Rev_D;Rev_nD];

%Stoichiometric Matrix reorganization
Sam = [S(:,Irr_D),S(:,Irr_nD),S(:,Rev_D),S(:,Rev_nD),-S(:,Rev_D),-S(:,Rev_nD)];  
N_Rxns = size(Sam,2);N_Mets = size(Sam,1);
N_Irr = length(IrrRxns);N_Rev = length(RevRxns);
D_Irr = D(Irr_D);
D_Rev = D(Rev_D);
D_For = D_Rev;

if length(Dstd) > 1,
    DstdIrr = Dstd(Irr_D);
    DstdFor = Dstd(Rev_D);
    DstdRev = DstdFor;
end

%Construction of L1 sequence
if L_min == L_max,
    L = L_min;
else
   L = L_min:L_int:L_max;
end

V = zeros(size(S,2),length(L));
D_cor = zeros(length(L),1);
D_dis = D_cor;
Card = D_cor;
CardData = D_cor;
iter_time = D_cor;

%Construction of bounding vectors
ub_rxns = v_max*ones(N_Rxns,1);
lb_rxns = zeros(N_Rxns,1);

if length(rxns_lb) ~= 1,
  for h = 1:size(S,2),
      lb_rxns(ismember(Rxns_Or,h)) = rxns_lb(h);
  end
end

if length(rxns_ub) ~= 1,
    for h = 1:size(S,2),
        ub_rxns(ismember(Rxns_Or,h)) = rxns_ub(h);
    end
end

if ~isempty(tasks),
   lb_rxns((ismember(Rxns_Or,tasks))) = v_eps;
end

if ~isempty(blocked),
   ub_rxns((ismember(Rxns_Or,blocked))) = v_eps;
end

Rev_lb = lb_rxns(RevRxns);
Rev_ub = ub_rxns(RevRxns);
lb = [zeros(2*NIrr_D + 4*NRev_D,1);lb_rxns;zeros(N_Rev,1)];
ub = [repmat(E_max,(2*NIrr_D + 4*NRev_D),1);ub_rxns;ones(N_Rev,1)];

%Include constraints on essential metabolites
if ~isempty(essmet),
    essmetmat = zeros(length(essmet),N_Rxns);
    for k = 1:length(essmet),
        essmetmat(k,Sam(essmet(k),:) > 0) = 1;
    end
    vecsense = [repmat('=',N_Mets,1);repmat('=',(NIrr_D+NRev_D),1);repmat('<',N_Rev,1);repmat('=',NRev_D,1);repmat('<',N_Rev,1);repmat('>',N_Rev,1);repmat('>',N_Rev,1);repmat('>',length(essmet),1)];
    b = [zeros(N_Mets,1);D_Irr;D_For;Rev_ub;zeros(NRev_D,1);zeros(N_Rev,1);Rev_lb;zeros(N_Rev,1);v_eps*ones(length(essmet),1)];
else
    vecsense = [repmat('=',N_Mets,1);repmat('=',(NIrr_D+NRev_D),1);repmat('<',N_Rev,1);repmat('=',NRev_D,1);repmat('<',N_Rev,1);repmat('>',N_Rev,1);repmat('>',N_Rev,1)];
    b = [zeros(N_Mets,1);D_Irr;D_For;Rev_ub;zeros(NRev_D,1);zeros(N_Rev,1);Rev_lb;zeros(N_Rev,1)];
end

%Construction of A matrix: E+irrdat, E-irrdat, E+for, E-for, E+rev, E-rev,Irr_D, Irr_nD, For_D,
%For_nD, Rev_D, Rev_nD, Y
A0 = [zeros(N_Mets,2*NIrr_D + 4*NRev_D),Sam,zeros(N_Mets,N_Rev)]; %SV=0
A1 = [diag(ones(NIrr_D,1)),-diag(ones(NIrr_D,1)),zeros(NIrr_D,4*NRev_D),diag(ones(NIrr_D,1)),zeros(NIrr_D,(N_Irr - NIrr_D)),zeros(NIrr_D,3*N_Rev)]; %V+E=D, IrrRxns
A2 = [zeros(NRev_D,2*NIrr_D),diag(ones(NRev_D,1)),-diag(ones(NRev_D,1)),zeros(NRev_D,2*NRev_D),zeros(NRev_D,N_Irr),diag(ones(NRev_D,1)),zeros(NRev_D,NRev_nD),zeros(NRev_D,N_Rev),diag(D_For),zeros(NRev_D,NRev_nD)]; %V+E+y*D=D, VforRxns
A3 = [zeros(N_Rev,2*NIrr_D + 4*NRev_D + N_Irr),diag(ones(N_Rev,1)),zeros(N_Rev,N_Rev),diag(Rev_ub)];%Vfor+y*Vformax<=Vformax
A4 = [zeros(NRev_D,2*NIrr_D),zeros(NRev_D,2*NRev_D),diag(ones(NRev_D,1)),-diag(ones(NRev_D,1)),zeros(NRev_D,N_Irr),zeros(NRev_D,N_Rev),diag(ones(NRev_D,1)),zeros(NRev_D,NRev_nD),diag(-D_Rev),zeros(NRev_D,NRev_nD)]; %V+E-y*D=0, VrevRxns
A5 = [zeros(N_Rev,2*NIrr_D + 4*NRev_D + N_Irr + N_Rev),diag(ones(N_Rev,1)),-diag(Rev_ub)]; %Vrev-y*Vrevmax<=0
A6 = [zeros(N_Rev,2*NIrr_D + 4*NRev_D + N_Irr),diag(ones(N_Rev,1)),zeros(N_Rev,N_Rev),diag(Rev_lb)];%Vfor+y*Vformin>=Vformin
A7 = [zeros(N_Rev,2*NIrr_D + 4*NRev_D + N_Irr + N_Rev),diag(ones(N_Rev,1)),-diag(Rev_lb)]; %Vrev-y*Vrevmin>=0

if ~isempty(essmet),
    A14 = zeros(length(essmet),size(A0,2));
    for k = 1:length(essmet),
        A14(k,:) = [zeros(1,2*NIrr_D+4*NRev_D),essmetmat(k,:),zeros(1,N_Rev)];
    end
Amat = [A0;A1;A2;A3;A4;A5;A6;A7];
elseif isempty(essmet), 
   Amat = [A0;A1;A2;A3;A4;A5;A6;A7];
end

%Integrates weighting based on standard deviation of each data point
if length(Dstd) > 1,
    W = 1./[DstdIrr;DstdIrr;DstdFor;DstdFor;DstdRev;DstdRev];
    W(W == Inf) = max(W);
else
    W = Dstd;
end

%Prepare gurobi structure
m.A = sparse(Amat);
m.rhs = b;
m.sense = vecsense;
m.vtype = [repmat('C',2*NIrr_D + 4*NRev_D + N_Rxns,1);repmat('B',N_Rev,1)];
m.lb = lb;
m.ub = ub;
params.OutputFlag = OutFlag;
params.Presolve = 2;
params.MIPGap = MIPGap;

%Iteration to determine optimum L
for i = 1:length(L),
    tic
    Lvec = repmat(L(i),N_Rxns,1);
    cvec = [W.*ones(2*NIrr_D + 4*NRev_D,1);Lvec;zeros(N_Rev,1)];
    
    %Solve MILP
    m.obj = cvec;
    gur = gurobi(m,params);
    
    %Reconstruct flux vector from gurobi solution
    X = gur.x;
    Virr_D = X((2*NIrr_D + 4*NRev_D + 1):(2*NIrr_D + 4*NRev_D + NIrr_D));
    Virr_nD = X((2*NIrr_D + 4*NRev_D + NIrr_D + 1):(2*NIrr_D + 4*NRev_D + N_Irr));
    Vfor_D = X((2*NIrr_D + 4*NRev_D + N_Irr + 1):(2*NIrr_D + 4*NRev_D + N_Irr + NRev_D));
    Vfor_nD = X((2*NIrr_D + 4*NRev_D + N_Irr + NRev_D + 1):(2*NIrr_D + 4*NRev_D + N_Irr + N_Rev));
    Vrev_D = X((2*NIrr_D + 4*NRev_D + N_Irr + N_Rev + 1):(2*NIrr_D + 4*NRev_D + N_Irr + N_Rev + NRev_D));
    Vrev_nD = X((2*NIrr_D + 4*NRev_D + N_Irr + N_Rev + NRev_D + 1):(2*NIrr_D + 4*NRev_D + N_Irr + 2*N_Rev));
    V(:,i) = zeros(size(S,2),1);V(Irr_D,i) = Virr_D;V(Irr_nD,i) = Virr_nD;
    V(Rev_D,i) = Vfor_D-Vrev_D;V(Rev_nD,i) = Vfor_nD - Vrev_nD;
    V((abs(V(:,i)) < v_eps),i) = 0;
    
    %Compute correlation between flux and data values
    if max(abs(V(:,i))) > 0,
        Cor = corrcoef(abs(V(~isnan(D_or),i)),D_or(~isnan(D_or)));
        D_cor(i) = Cor(1,2);
    end
    D_dis(i) = sum(sqrt((abs(V(~isnan(D_or),i)) - D_or(~isnan(D_or))).^2))/(NIrr_D + NRev_D);
    Card(i) = length(find(abs(V(:,i)) > 0));
    CardData(i) = length(intersect(find(abs(V(:,i)) > 0),find(~isnan(D_or))));
    iter_time(i) = toc; 
end

L_opt = find(D_cor == max(D_cor));
if length(Dstd) > 1,
    Weight = zeros(size(V,1),1);
    Weight(Irr_D) = W(1:NIrr_D);Weight(Rev_D) = W(2*NIrr_D + 1:2*NIrr_D + NRev_D);
    out.W = Weight;
    ZE = sum(W.*gur.x(1:2*NIrr_D + 4*NRev_D));
elseif length(Dstd) == 1,
    ZE = sum(gur.x(1:2*NIrr_D + 4*NRev_D));
end

%Build solution structure
out.Flux = V(:,L_opt);
out.Cor = D_cor(L_opt);
out.Res = D_dis(L_opt);
out.ActRxn = find(abs(V(:,L_opt)) > 0);
out.ActRxnData = intersect(find(abs(V(:,L_opt)) > 0),find(~isnan(D_or)));
out.Card = Card(L_opt);
out.CardData = CardData(L_opt);
out.Lambda = L(L_opt);
out.ZE = ZE;
out.ZV = sum(abs(out.Flux));
out.Time = iter_time;
if ~isempty(contextname),
   out.ContextCOBRA = contextmodel2COBRA(V(:,L_opt),GEM,contextname,0,0);
end
end
