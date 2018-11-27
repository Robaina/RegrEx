function Val=RegrExFVA(GEM,D,RegrExSol,Options)
      
%************************RegrEx FVA Implementation*************************

%Depends on Gurobi, which can be found in www.gurobi.com (free academic
%licenses available) FVA function (provided in supplementary material)
%
%Arguments:
%
% Required
%
% GEM: Metabolic model in COBRA structure. Blocked reactions should be
% removed first, for instance using the reduceModel function of the COBRA
% toolbox D: Experimental data vector (values already mapped to
% reactions)must be of same length as the number of reactions. It can be
% provided by the genetorxn function Vopt: flux vector obtained by RegrEx
%
% Optional
%
% Options must be a structure with the following fields (not all required)
%
%     Options.essmet: A vector indicating the indexes of the metabolites
%     that should be present in the final model 
%     
%     Options.tasks: A vector indicating the indexes of the reactions that 
%     should be active in the final model (i.e. flux value above eps) 
%     
%     Options.blocked: A vector indicating the indexes of the reactions that 
%     should not be active in the final model (i.e. zero flux value) 
%     
%     Options.Dstd: a vector containing the standard deviation of data, it
%     is intended to be used as weighting factor for reactions with
%     associated data, in the form (1/std(D))*||V-D||^2. Options.essmet: A
%     vector indicating the indexes of the metabolites that should be
%     present in the final model (This has been not widely tested)
%     
%     Options.rxns_lb: a vector containing the lower bounds for each reaction 
%     (default is set to 0,min value is 0, reversible reactions are splitted) 
%     
%     Options.rxns_ub: a vector containing the upper bounds for each reaction
%     (default is set to 1) 
%     
%     Options.time_limit: scalar indicating the time limit in seconds
%     for the MIQP (default is set to 60) 
%     
%     Options.eps: scalar indicating the threshold value to consider a 
%     reaction active, i.e. a reaction i with flux value vi is active if 
%     vi>=eps Options.E_max: vector of length equal to the number of 
%     reactions indicating the maximum allowed error value between data and 
%     flux prediction, i.e. E=d-v,(default is set to 1000)
%
%     Options.reaction: natural number indicating the index of a reaction,
%     RegrEx FVA will be applied just on this reaction and the 2 flux
%     distributions associated with the min and max will be returned.
%
%    
%Value:
%  
%  Vmin: vector containing the minimum flux values for each reaction under
%  optimal conditions, i.e. the minimum value in the alternate optima set
%  Vmax: vector containing the maximum flux values for each reaction under
%  optimal conditions, i.e. the maximum value in the alternate optima set
%
%**************************************************************************
%           Semidán (robaina@mpimp-golm.mpg.de), November, 2014
%**************************************************************************

%Argument evaluation

if ~exist('GEM','var'),
    error('The Genome-Scale model is missing...')
end

if ~exist('D','var'),
    error('The data vector is missing...')
end

if ~exist('RegrExSol','var'),
    error('The RegrEx solution structure is missing...')
end

if ~exist('Options','var'),
    Options=struct;
end

if ~isfield(Options,'tasks'), 
    tasks=[];
else
    tasks=Options.tasks;
end

if ~isfield(Options,'blocked'), 
    blocked=[];
else
    blocked=Options.blocked;
end

if ~isfield(Options,'rxns_lb'), 
    rxns_lb=0;
else
    rxns_lb=Options.rxns_lb;
end

if ~isfield(Options,'rxns_ub'), 
    rxns_ub=1;
else
    rxns_ub=Options.rxns_ub;
end

if ~isfield(Options,'time_limit'), 
    time_limit=60;
else
    time_limit=Options.time_limit;
end

if ~isfield(Options,'MaxCapacity'),
    MaxCapacity=1;
else
    MaxCapacity=Options.MaxCapacity;
end

if ~isfield(Options,'eps'), 
    eps=MaxCapacity*1e-6;
else
    eps=Options.eps;
end

if ~isfield(Options,'FeasTol'),
    FeasTol=1e-9;
else
    FeasTol=Options.FeasTol;
end

if ~isfield(Options,'E_max'), 
    E_max=MaxCapacity*1e3;
else
    E_max=Options.E_max;
end

if ~isfield(Options,'OutFlag'), 
    OutFlag=0;
else
    OutFlag=Options.OutFlag;
end

if ~isfield(Options,'deltaZE'),
    deltaZE=1e-4;
elseif isfield(Options,'deltaZE'),
    deltaZE=Options.deltaZE;
end

    S=GEM.S;
    Rev=find(GEM.rev==1);
    Irr=(setdiff(1:size(S,2),Rev))';
    Vmi=[];Vma=[];
    
    %Reaction partition: Irr_D, Irr_nD, Bfordat, Bfornodat, Rev_D, Rev_nD
    
    Irr_D=Irr((D(Irr)~=0));
    Irr_nD=Irr((D(Irr)==0));
    Rev_D=Rev((D(Rev)~=0));
    Rev_nD=Rev((D(Rev)==0));
    NIrr_D=length(Irr_D);
    NIrr_nD=length(Irr_nD);
    NRev_D=length(Rev_D);
    NRev_nD=length(Rev_nD);
    
    %Parse Data
   
    for i=1:length(D),
       if isnan(D(i)) || isinf(D(i)),
            D(i)=0;
            sprintf('Error: Data point D(%d) is NaN or Inf and has been converted to 0',i);
       end
    end
    if max(D)>MaxCapacity,
       D=D/max(D);
    end
    D=MaxCapacity*D;
    D=[D;D(Rev)];
    DIrr=D(Irr_D);
    DRev=D(Rev_D);
    DFor=DRev;
    ZV=RegrExSol.ZV;
    ZE=RegrExSol.ZE;
  
    %Stoichiometric Matrix reorganization
    
    Sam=[S(:,Irr_D),S(:,Irr_nD),S(:,Rev_D),S(:,Rev_nD),-S(:,Rev_D),-S(:,Rev_nD)];  
    Rxns=size(Sam,2);Mets=size(Sam,1);
    NIrr=length(Irr);NRev=length(Rev);
    Vmin=zeros(size(S,2),1);
    Vmax=Vmin;
    Rxns_Or=[Irr_D;Irr_nD;Rev_D;Rev_nD;Rev_D;Rev_nD];

    %Construction of sense, c and b vectors, definition of tasks and blocked reactions
    
    ubrxns=MaxCapacity*ones(Rxns,1);
    lbrxns=zeros(Rxns,1);
    if length(rxns_lb)~=1,
      for wbar=1:size(S,2),
          lbrxns(ismember(Rxns_Or,wbar))=rxns_lb(wbar);
      end
    end
    if length(rxns_ub)~=1,
      for wbar=1:size(S,2),
          ubrxns(ismember(Rxns_Or,wbar))=rxns_ub(wbar);
      end
    end
    if ~isempty(tasks),
        for i=1:length(tasks),
            if tasks(i)==454,
                lbrxns((ismember(Rxns_Or,tasks(i))))=eps;
            elseif tasks(i)~=454,
                lbrxns((ismember(Rxns_Or,tasks(i))))=0.001;
            end
        end
    end
    if ~isempty(blocked),
       ubrxns((ismember(Rxns_Or,blocked)))=eps;
    end
    
    Revub=ubrxns(Rev);
    lb=[zeros((2*NIrr_D+4*NRev_D),1);lbrxns;zeros(NRev,1)];
    ub=[E_max*ones((2*NIrr_D+4*NRev_D),1);ubrxns;ones(NRev,1)];
    vecsense=[repmat('=',Mets,1);repmat('=',(NIrr_D+NRev_D),1);repmat('<',NRev,1);repmat('=',NRev_D,1);repmat('<',NRev,1);'<';'>';repmat('=',1,1)];     
    b=[zeros(Mets,1);DIrr;DFor;Revub+eps;zeros(NRev_D,1);zeros(NRev,1)+eps;ZE+deltaZE*ZE;ZE-deltaZE*ZE;ZV];
   
    %Integrates weighting based on standard deviation of each data point
   
    if isfield(RegrExSol,'W'),
        W=RegrExSol.W;
        W=[W(Irr_D);W(Irr_D);W(Rev_D);W(Rev_D);W(Rev_D);W(Rev_D)]';
    elseif ~isfield(RegrExSol,'W'),
        W=1;
    end
 
    %Construction of A matrix: E+irr_D,E-irr_D,E+for_D,E-for_D,E+rev_D,E-rev_D,Virr_D,Virr_nD,Vfor_D,Vfor_nD,Vrev_D,Vrev_nD,x
    
    A0=[zeros(Mets,2*NIrr_D+4*NRev_D),Sam,zeros(Mets,NRev)]; %SV=0
    A1=[diag(ones(NIrr_D,1)),-diag(ones(NIrr_D,1)),zeros(NIrr_D,4*NRev_D),diag(ones(NIrr_D,1)),zeros(NIrr_D,NIrr_nD),zeros(NIrr_D,3*NRev)];%V+E=D for Irr^Data
    A2=[zeros(NRev_D,2*NIrr_D),diag(ones(NRev_D,1)),-diag(ones(NRev_D,1)),zeros(NRev_D,2*NRev_D),zeros(NRev_D,NIrr),diag(ones(NRev_D,1)),zeros(NRev_D,NRev_nD),zeros(NRev_D,NRev),diag(DFor),zeros(NRev_D,NRev_nD)];%V+E=D Vfor^Data
    A3=[zeros(NRev,2*NIrr_D+4*NRev_D+NIrr),diag(ones(NRev,1)),zeros(NRev,NRev),diag(Revub)]; %Vfor direction constrain ub
    A4=[zeros(NRev_D,2*NIrr_D),zeros(NRev_D,2*NRev_D),diag(ones(NRev_D,1)),-diag(ones(NRev_D,1)),zeros(NRev_D,NIrr),zeros(NRev_D,NRev),diag(ones(NRev_D,1)),zeros(NRev_D,NRev_nD),diag(-DRev),zeros(NRev_D,NRev_nD)];%V+E=D Vrev^Data
    A5=[zeros(NRev,2*NIrr_D+4*NRev_D+NIrr+NRev),diag(ones(NRev,1)),-diag(Revub)]; %Vrev direction constrain ub
    A6a=[W.*ones(1,2*NIrr_D+4*NRev_D),zeros(1,Rxns+NRev)];%(E+)+(E-) < ZE + eps
    A6b=[W.*ones(1,2*NIrr_D+4*NRev_D),zeros(1,Rxns+NRev)];%(E+)+(E-) > ZE - eps
    A7=[zeros(1,2*NIrr_D+4*NRev_D),ones(1,NIrr+NRev),ones(1,NRev),zeros(1,NRev)]; %||V||_1=ZV
    
    %Include carboxylation to oxygenation ratio: 1.5 <= V6/V85 <= 4 and
    %constraint biomass production
    if isfield(Options,'Carb2Oxy') && Options.Carb2Oxy=='T',   
     Rat1=zeros(1,Rxns);Rat1(ismember(Rxns_Or,6))=1;Rat1(ismember(Rxns_Or,85))=-4;
     Rat2=zeros(1,Rxns);Rat2(ismember(Rxns_Or,6))=1;Rat2(ismember(Rxns_Or,85))=-1.5;
     A8=[zeros(1,2*NIrr_D+4*NRev_D),Rat1,zeros(1,NRev)];
     A9=[zeros(1,2*NIrr_D+4*NRev_D),Rat2,zeros(1,NRev)];
     vecsense=[repmat('=',Mets,1);repmat('=',(NIrr_D+NRev_D),1);repmat('<',NRev,1);repmat('=',NRev_D,1);repmat('<',NRev,1);'<';'>';repmat('=',1,1);'<';'>'];     
     b=[zeros(Mets,1);DIrr;DFor;Revub+eps;zeros(NRev_D,1);zeros(NRev,1)+eps;ZE+deltaZE*ZE;ZE-deltaZE*ZE;ZV;0;0];
     Amat=[A0;A1;A2;A3;A4;A5;A6a;A6b;A7;A8;A9];
    elseif ~isfield(Options,'Carb2Oxy') || Options.Carb2Oxy=='F',
         Amat=[A0;A1;A2;A3;A4;A5;A6a;A6b;A7];
    end
    %Construction of Q matrix:

    Q=zeros(size(Amat,2),size(Amat,2));
    
    %***************************Main Loop Starts***************************
    %**********************************************************************
    
    if ~isfield(Options,'reaction'),
       wbar = waitbar(0,'Calculating Alternative Optima range...');
       RxnSeq=1:Rxns;
    elseif isfield(Options,'reaction'),
       RxnSeq=find(Rxns_Or==Options.reaction);
    end
    for i=1:length(RxnSeq),
        if ~isfield(Options,'reaction'),
          waitbar(i/length(RxnSeq))
        end
        
        %Select reaction to maximize
        cvec=zeros(size(Amat,2),1);cvec(2*NIrr_D+4*NRev_D+RxnSeq(i),1)=1;
        vtype=[repmat('C',2*NIrr_D+4*NRev_D+Rxns,1);repmat('B',NRev,1)];
        
        %Solve MILPs             
        
        m.Q=sparse(Q);
        m.obj=cvec;
        m.A=sparse(Amat);
        m.rhs=b;
        m.sense=vecsense;
        m.modelsense='max';
        m.vtype=vtype;
        m.lb=lb;
        m.ub=ub;
        params.Presolve=2;
        params.OutputFlag=OutFlag;
        params.FeasibilityTol=FeasTol;
        params.TimeLimit=time_limit;
        params.Threads=8;
        gurmax=gurobi(m,params);
        
        m.Q=sparse(Q);
        m.obj=cvec;
        m.A=sparse(Amat);
        m.rhs=b;
        m.sense=vecsense;
        m.modelsense='min';
        m.vtype=vtype;
        m.lb=lb;
        m.ub=ub;
        gurmin=gurobi(m,params);
       
        if ~isfield(Options,'reaction'),
            try
             Vma(i)=gurmax.objval;
            catch
                Vma(i)=nan;
            end
            try
             Vmi(i)=gurmin.objval;
            catch
                Vmi(i)=nan;
            end 
        elseif isfield(Options,'reaction'),
            if GEM.rev(Options.reaction),
                try
                  Vma(i)=gurmax.objval;
                  Vmi(i)=gurmin.objval;
                catch
                    Vma(i)=nan;
                    Vmi(i)=nan;
                end
            end
        end             
    end
    
    if ~isfield(Options,'reaction'),
        %reconstruct flux vector with non-split reversible reactions
        VIrrDma=Vma(1:NIrr_D);
        VIrrnDma=Vma((NIrr_D+1):NIrr);
        VforDma=Vma((NIrr+1):(NIrr+NRev_D));
        VfornDma=Vma((NIrr+NRev_D+1):(NIrr+NRev));
        VrevDma=Vma((NIrr+NRev+1):(NIrr+NRev+NRev_D));
        VrevnDma=Vma((NIrr+NRev+NRev_D+1):(NIrr+2*NRev));
        
        VIrrDmi=Vmi(1:NIrr_D);
        VIrrnDmi=Vmi((NIrr_D+1):NIrr);
        VforDmi=Vmi((NIrr+1):(NIrr+NRev_D));
        VfornDmi=Vmi((NIrr+NRev_D+1):(NIrr+NRev));
        VrevDmi=Vmi((NIrr+NRev+1):(NIrr+NRev+NRev_D));
        VrevnDmi=Vmi((NIrr+NRev+NRev_D+1):(NIrr+2*NRev));
        
        Vmin(Irr_D)=VIrrDmi;Vmin(Irr_nD)=VIrrnDmi;
        Vmax(Irr_D)=VIrrDma;Vmax(Irr_nD)=VIrrnDma;
        
        VminSplit=Vmin;VmaxSplit=Vmax;
        VminSplit(Rev_D)=VforDmi;VminSplit(Rev_nD)=VfornDmi;
        VmaxSplit(Rev_D)=VforDma;VmaxSplit(Rev_nD)=VfornDma;
        VmiS=nan*ones(size(GEM.S,2),1);VmaS=VmiS;
        VmiS(Rev_D)=VrevDmi;VmiS(Rev_nD)=VrevnDmi;
        VmaS(Rev_D)=VrevDma;VmaS(Rev_nD)=VrevnDma;
        VminSplit=[VminSplit;VmiS(~isnan(VmiS))];
        VmaxSplit=[VmaxSplit;VmaS(~isnan(VmaS))];
        
        for idRevD=1:NRev_D,
            if VrevDma(idRevD)>eps,
               Vmin(Rev_D(idRevD))=-VrevDma(idRevD);
            elseif VrevDma(idRevD)<=eps,
                Vmin(Rev_D(idRevD))=VforDmi(idRevD);
            end
        end
        for idRevnD=1:NRev_nD,
            if VrevnDma(idRevnD)>eps,
               Vmin(Rev_nD(idRevnD))=-VrevnDma(idRevnD);
            elseif VrevnDma(idRevnD)<=eps,
                Vmin(Rev_nD(idRevnD))=VfornDmi(idRevnD);
            end
        end
        for idRevD=1:NRev_D,
            if VforDma(idRevD)>eps,
               Vmax(Rev_D(idRevD))=VforDma(idRevD);
            elseif VforDma(idRevD)<=eps,
                Vmax(Rev_D(idRevD))=-VrevDmi(idRevD);
            end
        end
        for idRevnD=1:NRev_nD,
            if VfornDma(idRevnD)>eps,
               Vmax(Rev_nD(idRevnD))=VfornDma(idRevnD);
            elseif VfornDma(idRevnD)<=eps,
                Vmax(Rev_nD(idRevnD))=-VrevnDmi(idRevnD);
            end
        end
        
        %reconstruct flux vector with split reversible reactions
        
        
    elseif isfield(Options,'reaction'),
        if ~GEM.rev(Options.reaction),
            try
              Vmin=gurmin.objval;
            catch
                Vmin=nan;
            end
            try
               Vmax=gurmax.objval;
            catch
                Vmax=nan;
            end
        elseif GEM.rev(Options.reaction),
            if Vma(2)<=eps,
                Vmin=Vmi(1);
            elseif Vma(2)>eps,
                Vmin=-Vma(2);
            end
            if Vma(1)<=eps,
                Vmax=-Vmi(2);
            elseif Vma(1)>eps,
                Vmax=Vma(1);
            end    
        end
    end    
    
    if ~isfield(Options,'reaction'),
        close(wbar)
        Range(:,1)=Vmin;
        Range(:,2)=Vmax;
%         Range(abs(Val.Range)<eps)=0;
        Val.Range=sort(Range,2);
        Val.RangeSplit=[VminSplit,VmaxSplit];
    elseif isfield(Options,'reaction'),
        Range=[Vmin,Vmax];
        Val.Range=sort(Range,2);
%         Range(abs(Range)<eps)=0;
    end 

end
