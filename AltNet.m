function Sol = AltNet(GEM,SolCorEx,Options)

%***************************The AltNet procedure***************************
%**************************************************************************
% 
%This function implements the AltNet procedure to generate alternative
%optimal networks to CorEx and FastCORE.
% 
%Required arguments are:
%   GEM: COBRA-like structure with the GEM
%   
%   SolCorEx: the MATLAB structure returned by CorEx 
%   
% Optional values are fields of the structure "Options":
%   v_max: an upper bound for the maximum flux value of any reaction (in
%   reversibles the bound would be -Vmax <= v <= Vmax), default 1e3
%   
%   v_eps: the numerical threshold to consider a reaction active (default
%   1e-3)
%   
%   Nsamples: number of generated alternative optimal models
%   
%   time_limit: a time limit (in seconds) for an early termination of the
%   MILP 
%
%   MIPGap: the relative gap between the lower and upper bound of the MIP
%   optimum. It serves as an early termination mechanism. Can be set from
%   1e-4 to inf.
%   
%   iter_limit: an upper bound on the number of iterations in the main loop.
%   The loop continues until Nmodels is satisfied or the number of
%   iterations reaches IterLimit
%
%   out_flag: set to 1 for verbose, 0 otherwise
%
%   warm_start: if set to 1 then gurobi is fed with a binary vector used as
%   a warm start for the MILP. This vector corresponds to a previous
%   optimal vector of CorEx.
%   
%Output (a MATLAB structure with the fields):
%   QualityCheck: counts the number of core (C) and non-core(P) reactions
%   in the CorEx model. Their value have to much among alternative optimal
%   networks
%
%   Modmatrix: a binary matrix. Sampled models are coded in its columns. A
%   reaction is present in the final alternative optimal model if the
%   corresponding entry is 1
%
%   maxDiff: an array containing the Hamming distance between the sampled
%   alternative optimal models (columns in Modmatrix). Due to the nature of
%   the AltNet algorithm, the distance is computed between columns n + 1, n
%
%   TotalTime: total computational time in seconds
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
if ~isfield(Options,'Nsamples'),
    Nsamples = 10;
else
    Nsamples = Options.Nsamples;
end
if ~isfield(Options,'time_limit'),
    time_limit = [];
else
    time_limit = Options.time_limit;
end
if ~isfield(Options,'iter_limit'),
    iter_limit = 100;
else
    iter_limit = Options.iter_limit;
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
if ~isfield(Options,'warm_start'),
    warm_start = 0;
else
    warm_start = Options.warm_start;
end

Z = SolCorEx.Z;
Abinary = SolCorEx.Abinary;
C = SolCorEx.C;

%Categorize reactions into reversible and irreversible and core and
%non-core
RevRxns = find(GEM.rev == 1);
IrrRxns = find(GEM.rev == 0);
N_rev = length(RevRxns);
N_irr = length(IrrRxns);
P = setdiff(1:length(GEM.rxns),C);
S = GEM.S;

Irr_C = intersect(C,IrrRxns);
Irr_P = intersect(P,IrrRxns);
Rev_C = intersect(C,RevRxns);
Rev_P = intersect(P,RevRxns);
NIrr_C = length(Irr_C);
NIrr_P = length(Irr_P);
NRev_C = length(Rev_C);
NRev_P = length(Rev_P);

%Stoichiometric Matrix reorganization
Sam = [S(:,Irr_C),S(:,Irr_P),S(:,Rev_C),S(:,Rev_P),-S(:,Rev_C),-S(:,Rev_P)];  
N_rxns = size(Sam,2);
N_mets = size(Sam,1);

%Construction of Amat (constraints) matrix (virr_C,virr_P,vfor_C,vfor_P,vrev_C,vrev_P,xirr(P),xRev(P),y,delta+,delta-)
A1 = [Sam,zeros(N_mets,NIrr_P+NRev_P+N_rev+2*(NIrr_P+NRev_P))]; %SV=0
A2 = [eye(NIrr_C),zeros(NIrr_C,NIrr_P+2*N_rev+NIrr_P+NRev_P+N_rev+2*(NIrr_P+NRev_P))]; %virr_C>=epsilon 
A3 = [zeros(NRev_C,N_irr),eye(NRev_C),zeros(NRev_C,NRev_P),eye(NRev_C),zeros(NRev_C,NRev_P+NIrr_P+NRev_P+N_rev+2*(NIrr_P+NRev_P))]; %vfor_C+vrev_C>=epsilon

A4a = [zeros(NIrr_P,NIrr_C),eye(NIrr_P),zeros(NIrr_P,2*N_rev),-v_max*eye(NIrr_P),zeros(NIrr_P,NRev_P+N_rev+2*(NIrr_P+NRev_P))]; %virr_p-xirr*Vmax<=0
A4b = [zeros(NIrr_P,NIrr_C),eye(NIrr_P),zeros(NIrr_P,2*N_rev),-v_eps*eye(NIrr_P),zeros(NIrr_P,NRev_P+N_rev+2*(NIrr_P+NRev_P))]; %virr_p-xirr*epsilon>=0

A5a = [zeros(NRev_P,N_irr+NRev_C),eye(NRev_P),zeros(NRev_P,NRev_C),eye(NRev_P),zeros(NRev_P,NIrr_P),-v_max*eye(NRev_P),zeros(NRev_P,N_rev+2*(NIrr_P+NRev_P))]; %(vfor_p + vback_p)-xRev*Vmax<=0
A5b = [zeros(NRev_P,N_irr+NRev_C),eye(NRev_P),zeros(NRev_P,NRev_C),eye(NRev_P),zeros(NRev_P,NIrr_P),-v_eps*eye(NRev_P),zeros(NRev_P,N_rev+2*(NIrr_P+NRev_P))]; %(vfor_p + vback_p)-xRev*epsilon>=0

A7 = [zeros(N_rev,N_irr),eye(N_rev),zeros(N_rev,N_rev+NIrr_P+NRev_P),v_max*eye(N_rev),zeros(N_rev,+2*(NIrr_P+NRev_P))]; %vfor+y*Vmax<=Vmax
A8 = [zeros(N_rev,N_irr+N_rev),eye(N_rev),zeros(N_rev,NIrr_P+NRev_P),-v_max*eye(N_rev),zeros(N_rev,+2*(NIrr_P+NRev_P))]; %vrev-y*Vmax<=0
A9a = [zeros(1,N_rxns),ones(1,NIrr_P+NRev_P),zeros(1,N_rev+2*(NIrr_P+NRev_P))]; %||x||_1<Z+epsilon
A9b = [zeros(1,N_rxns),ones(1,NIrr_P+NRev_P),zeros(1,N_rev+2*(NIrr_P+NRev_P))]; %||x||_1>Z-epsilon

A10 = [zeros(NIrr_P + NRev_P,N_rxns),eye(NIrr_P + NRev_P),zeros(NIrr_P + NRev_P,N_rev),eye(NIrr_P + NRev_P),-eye(NIrr_P + NRev_P)]; %delta+ - delta- + x = Xopt
A11 = [zeros(NIrr_P + NRev_P,N_rxns + NIrr_P + NRev_P + N_rev),eye(NIrr_P + NRev_P),eye(NIrr_P + NRev_P)]; %delta+ + delta- <= 1

Amat = [A1;A2;A3;A4a;A4b;A5a;A5b;A7;A8;A9a;A9b;A10;A11];

%Construction of c, lb, ub, sense and b vectors
c = [zeros(N_rxns + NIrr_P + NRev_P + N_rev,1);ones(2*(NIrr_P + NRev_P),1)];
vsense = [repmat('=',N_mets,1);repmat('>',NIrr_C + NRev_C,1);repmat('<',NIrr_P,1);repmat('>',NIrr_P,1);repmat('<',NRev_P,1);repmat('>',NRev_P,1);repmat('<',2*N_rev,1);'<';'>';repmat('=',NIrr_P+NRev_P,1);repmat('<',NIrr_P+NRev_P,1)];
vtype = [repmat('C',N_rxns,1);repmat('B',NIrr_P+NRev_P + N_rev,1);repmat('B',2*(NIrr_P + NRev_P),1)];
lb = zeros(size(Amat,2),1);
ub = [v_max*ones(N_rxns,1);ones(NIrr_P + NRev_P + N_rev + 2*(NIrr_P + NRev_P),1)];

%Prepare gurobi structure
model.A = sparse(Amat);
model.modelsense = 'max';
model.obj = c;
model.sense = vsense;
model.lb = lb;
model.ub = ub;
model.vtype = vtype;
params.OutputFlag = out_flag;
if ~isempty(time_limit),
   params.TimeLimit = time_limit;
end
if ~isempty(MIPGap),
   params.MIPGap = MIPGap;
end
params.Threads = 1; %parallelization does not seem to improve the time at all!
params.FeasibilityTol = 1e-9;
params.IntFeasTol = 1e-9;

Vmatrix = zeros(length(GEM.rxns),Nsamples);
Modmatrix = zeros(length(GEM.rxns),Nsamples+1);
QualityCheck = zeros(Nsamples,2);
MaxDifferences = zeros(Nsamples,1);
wbar = waitbar(0,'Getting Alternative Networks...');
netcounter = 1;itercounter = 1;
Modmatrix(:,1) = Abinary;

%warm start for the binary variables
Xopt2 = SolCorEx.Xopt;
Xoptmatrix = zeros(length(Xopt2),Nsamples+1);
Xoptmatrix(:,1) = Xopt2;

tic
%Start Iteration
while netcounter < Nsamples && itercounter < iter_limit,
    waitbar(netcounter/Nsamples)
    %Solve MILP  
    model.rhs = [zeros(N_mets,1);v_eps*ones(NIrr_C+NRev_C,1);zeros(2*NIrr_P + 2*NRev_P,1);v_max*ones(N_rev,1);zeros(N_rev,1);Z + 0.5;Z - 0.5;Xopt2;ones(NIrr_P + NRev_P,1)];
   
    %Include warm start?
    if warm_start == 1,
       model.start = [NaN*ones(N_rxns,1);Xopt2;NaN*ones(N_rev + 2*(NIrr_P+NRev_P),1)];
    end

    try
      gur = gurobi(model,params);
      %Reconstruct flux vector
      Vmatrix(Irr_C,netcounter) = gur.x(1:NIrr_C);
      Vmatrix(Irr_P,netcounter) = gur.x(NIrr_C + 1:N_irr);
      Vmatrix(Rev_C,netcounter) = gur.x(N_irr + 1:N_irr+NRev_C) - gur.x(N_irr + N_rev + 1:N_irr + N_rev+NRev_C);
      Vmatrix(Rev_P,netcounter) = gur.x(N_irr + NRev_C + 1:N_irr + NRev_C + NRev_P) - gur.x(N_irr + N_rev + NRev_C + 1:N_irr + N_rev + NRev_C + NRev_P);
      Vmatrix(abs(Vmatrix) < (v_eps - 0.01*v_eps)) = 0;
      Modmatrix((abs(Vmatrix(:,netcounter)) > 0),netcounter + 1) = 1;
      %Update binary variables of the previous network for next iteration
      Xoptmatrix(:,netcounter + 1) = Xopt2;
      Xopt2 = gur.x(N_rxns+1:N_rxns+NIrr_P+NRev_P); 

    catch
        Modmatrix(:,netcounter + 1) = Modmatrix(:,netcounter);
    end
        MaxDifferences(netcounter) = sum(abs(Modmatrix(:,netcounter + 1) - Modmatrix(:,netcounter)));
        netcounter = netcounter + 1;

end

%Assesst quality of the networks: all networks must include the core set
%and be consistent (have non-zero flux only in all selected reactions)
for i = 1:Nsamples,
    QualityCheck(i,1) = length(find(abs(Vmatrix(C,i)) > 0));
    QualityCheck(i,2) = length(find(abs(Vmatrix(:,i)) > 0));
end
Modmatrix(:,sum(Modmatrix) == 0) = [];
QualityCheck = [{'Core Active Rxns','Total Active Rxns'};num2cell(QualityCheck)];
close(wbar)

%Construct solution structure
Sol.Modmatrix = unique(Modmatrix','rows')';
Sol.maxDiff = MaxDifferences;
Sol.QualityCheck = QualityCheck;
Sol.TotalTime = toc;
% Sol.Abinary=Abinary;
% Sol.Xoptmatrix=Xoptmatrix;
Sol.Vmatrix=Vmatrix;

end