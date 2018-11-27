function Sol = CorEx(GEM,C,Options)

%***************************The CorEx method*******************************
%**************************************************************************
%
%This function implements the CorEx method to obtain a context-specific
%model from a GEM and a core set of reactions.
%
%Required arguments are:
%  GEM: COBRA-like structure with the GEM
%  
%  C: an array with the indexes in the GEM of the reactions in the core set
%  
%Optional arguments are fields of the structure "Options" :
%  Z: the number of added non-core reactions, if added it must be equal or
%  greater than the optimal value found by a previous unconstrained CorEx
%  optimization.
%  
%  v_max: an upper bound for the maximum flux value of any reaction (in
%  reversibles the bound would be -Vmax <= v <= Vmax), default 1e3
%  
%  v_eps: the numerical threshold to consider a reaction active (default
%  1e-3)
%  
%  time_limit: a time limit (in seconds) for an early termination of the
%  MILP (default 60s
%
%  MIPGap: the relative gap between the lower and upper bound of the MIP
%  optimum. It serves as an early termination mechanism. Can be set from
%  1e-4 to inf.
%
%  out_flag: set to 1 for verbose
%
%Output (a MATLAB structure with the fields):
%   QualityCheck: counts the number of core (C) and non-core(P) reactions
%   in the CorEx model. Their value have to much among alternative optimal
%   networks
%
%   Z: the number of non-core reactions added by CorEx (i.e., its optimum)
%
%   Abinary: a binary array of the same size as the number of reactions in
%   the original GEM, 1 represent reactions that are included in the final,
%   context-specific model
%
%   AddedRxns: an array containing the indexes of the non-core reactions
%   added by CorEx
%
%   Xopt: the vector of optimal binary variable values returned by Gurobi.
%   It will be used by AltNet to facilite a warm start of the MILP
%  
%**************************************************************************
%           Semidán (robaina@mpimp-golm.mpg.de), May, 2016
%**************************************************************************

%Evaluate arguments
if ~exist('Options','var'),
    Options = struct;
end
if ~isfield(Options,'Z'),
    Z = [];
else
    Z = Options.Z;
end
if ~isfield(Options,'v_max'),
    v_max = 1e3;
else
    v_max = Options.v_max;
end
if ~isfield(Options,'v_eps'),
    v_eps = 1e-3;
else
    v_eps = Options.v_eps;
end
if ~isfield(Options,'time_limit'),
    time_limit = 60;
else
    time_limit = Options.time_limit;
end
if ~isfield(Options,'out_flag'),
    out_flag = 1;
else
    out_flag = Options.out_flag;
end
if ~isfield(Options,'MIPGap'),
    MIPGap = [];
else
    MIPGap = Options.MIPGap;
end

RevRxns = find(GEM.rev == 1);
IrrRxns = find(GEM.rev == 0);
N_rev = length(RevRxns);
N_irr = length(IrrRxns);
P = setdiff(1:length(GEM.rxns),C);

%Stoichiometric Matrix reorganization
Irr_C = intersect(C,IrrRxns);
Irr_P = intersect(P,IrrRxns);
Rev_C = intersect(C,RevRxns);
Rev_P = intersect(P,RevRxns);
NIrr_C = length(Irr_C);
NIrr_P = length(Irr_P);
NRev_C = length(Rev_C);
NRev_P = length(Rev_P);
S = GEM.S;
Sam = [S(:,Irr_C),S(:,Irr_P),S(:,Rev_C),S(:,Rev_P),-S(:,Rev_C),-S(:,Rev_P)];  
Rxns = size(Sam,2);
Mets = size(Sam,1);

%Construction of Amat (constraints) matrix (virr_C,virr_P,vfor_C,vfor_P,vrev_C,vrev_P,xirr(P),xRev(P),y)
A1 = [Sam,zeros(Mets,NIrr_P + NRev_P + N_rev)]; %SV=0
A2 = [eye(NIrr_C),zeros(NIrr_C,NIrr_P + 2*N_rev + NIrr_P + NRev_P + N_rev)]; %virr_C>=v_eps
A3 = [zeros(NRev_C,N_irr),eye(NRev_C),zeros(NRev_C,NRev_P),eye(NRev_C),zeros(NRev_C,NRev_P + NIrr_P + NRev_P + N_rev)]; %vfor_C+vrev_C>=v_eps

A4a = [zeros(NIrr_P,NIrr_C),eye(NIrr_P),zeros(NIrr_P,2*N_rev),-v_max*eye(NIrr_P),zeros(NIrr_P,NRev_P + N_rev)]; %virr_p-xirr*v_max<=0
A4b = [zeros(NIrr_P,NIrr_C),eye(NIrr_P),zeros(NIrr_P,2*N_rev),-v_eps*eye(NIrr_P),zeros(NIrr_P,NRev_P + N_rev)]; %virr_p-xirr*v_eps>=0

A5a = [zeros(NRev_P,N_irr + NRev_C),eye(NRev_P),zeros(NRev_P,NRev_C),eye(NRev_P),zeros(NRev_P,NIrr_P),-v_max*eye(NRev_P),zeros(NRev_P,N_rev)]; %(vfor_p + vback_p)-xRev*v_max<=0
A5b = [zeros(NRev_P,N_irr + NRev_C),eye(NRev_P),zeros(NRev_P,NRev_C),eye(NRev_P),zeros(NRev_P,NIrr_P),-v_eps*eye(NRev_P),zeros(NRev_P,N_rev)]; %(vfor_p + vback_p)-xRev*v_eps>=0

A6 = [zeros(N_rev,N_irr),eye(N_rev),zeros(N_rev,N_rev + NIrr_P + NRev_P),v_max*eye(N_rev)]; %vfor+y*v_max<=v_max
A7 = [zeros(N_rev,N_irr+N_rev),eye(N_rev),zeros(N_rev,NIrr_P + NRev_P),-v_max*eye(N_rev)]; %vrev-y*Vmax<=0
Amat = [A1;A2;A3;A4a;A4b;A5a;A5b;A6;A7];

%Construction of c, lb, ub, sense and b vectors
b = [zeros(Mets,1);v_eps*ones(NIrr_C + NRev_C,1);zeros(2*NIrr_P + 2*NRev_P,1);v_max*ones(N_rev,1);zeros(N_rev,1)];
vsense = [repmat('=',Mets,1);repmat('>',NIrr_C+NRev_C,1);repmat('<',NIrr_P,1);repmat('>',NIrr_P,1);repmat('<',NRev_P,1);repmat('>',NRev_P,1);repmat('<',2*N_rev,1)];
if ~isempty(Z),
   A9a = [zeros(1,Rxns),ones(1,NIrr_P+NRev_P),zeros(1,N_rev)]; %||x||_1<Z+v_eps
   A9b = [zeros(1,Rxns),ones(1,NIrr_P+NRev_P),zeros(1,N_rev)]; %||x||_1>Z-v_eps
   Amat = [Amat;A9a;A9b];
   b = [b;Z + 0.5;Z - 0.5];
   vsense = [vsense;'<';'>'];
end
c = [zeros(Rxns,1);ones(NIrr_P + NRev_P,1);zeros(N_rev,1)];
vtype = [repmat('C',Rxns,1);repmat('B',NIrr_P + NRev_P + N_rev,1)];
lb = zeros(size(Amat,2),1);
ub = [v_max*ones(Rxns,1);ones(NIrr_P + NRev_P + N_rev,1)];

%Build gurobi structure and solve MILP
model.A = sparse(Amat);
model.modelsense = 'min';
model.obj = c;
model.sense = vsense;
model.rhs = b;
model.lb = lb;
model.ub = ub;
model.vtype = vtype;
params.OutputFlag = out_flag;
params.FeasibilityTol = 1e-9;
params.IntFeasTol = 1e-9;
if ~isempty(time_limit),
   params.TimeLimit = time_limit;
end
if ~isempty(MIPGap),
   params.MIPGap = MIPGap;
end
gur = gurobi(model,params);

V = zeros(length(GEM.rxns),1);
V(Irr_C) = gur.x(1:NIrr_C);
V(Irr_P) = gur.x(NIrr_C + 1:N_irr);
V(Rev_C) = gur.x(N_irr + 1:N_irr + NRev_C) - gur.x(N_irr + N_rev + 1:N_irr + N_rev + NRev_C);
V(Rev_P) = gur.x(N_irr + NRev_C + 1:N_irr + NRev_C + NRev_P) - gur.x(N_irr + N_rev + NRev_C + 1:N_irr + N_rev + NRev_C + NRev_P);
V(abs(V) < (v_eps - 0.01*v_eps)) = 0;
A = zeros(length(V),1);A(abs(V) > 0) = 1;

QualityCheck(1,1) = length(find(abs(V(C)) > 0));
QualityCheck(1,2) = length(find(abs(V(P)) > 0));
QualityCheck(1,3) = length(find(abs(V) > 0));

%Construct solution structure
Sol.QualityCheck = [{'C Active Rxns','P Active Rxns','Total Active Rxns'};num2cell(QualityCheck)];
Sol.AddedRxns = find(abs(V(P)) > 0);
Sol.Z = length(find(abs(V) > 0)) - length(find(abs(V(C)) > 0));
Sol.Xopt = gur.x(Rxns + 1:Rxns + NIrr_P + NRev_P);
Sol.Abinary = A;
Sol.C = C;

end