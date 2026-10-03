% Estimate the constrained second-moment tensor from station observations.
% The design equation is G*m=d, where d contains temporal second moments and
% m=[mtt,mxt,myt,mxx,mxy,myy]. The 3-by-3 moment matrix is constrained to be
% positive semidefinite, and the residual norm is minimized with the LMI
% Control Toolbox functions FEASP and MINCX.
function [m2,Lc,Wc,vx,vy,tauc,info]=Step3_Inverse(Gmat,obs)
if nargin<2
    error('Step3_Inverse:NotEnoughInputs','Two inputs are required: G and d.');
end
Gmat=double(Gmat);
obs=double(obs(:));
if isempty(Gmat) || isempty(obs)
    error('Step3_Inverse:EmptyInput','G and d must not be empty.');
end
if size(Gmat,2)~=6
    error('Step3_Inverse:DesignMatrixSize','G must have exactly six columns.');
end
if size(Gmat,1)~=numel(obs)
    error('Step3_Inverse:SizeMismatch','G and d have incompatible dimensions.');
end
Energy=sum(obs.^2);

setlmis([]); % Initialize target objective
[Momvar(1),~,~]=lmivar(1,[1 1]);
[Momvar(2),~,~]=lmivar(1,[1 1]);
[Momvar(3),~,~]=lmivar(1,[1 1]);
[Momvar(4),~,~]=lmivar(1,[1 1]);
[Momvar(5),~,~]=lmivar(1,[1 1]);
[Momvar(6),~,~]=lmivar(1,[1 1]);

% Positive-semidefinite moment matrix:
% [mtt mxt myt; mxt mxx mxy; myt mxy myy] >= 0.
lmiterm([-1 1 1 Momvar(1)],1,1);
lmiterm([-1 1 2 Momvar(2)],1,1);
lmiterm([-1 1 3 Momvar(3)],1,1);
lmiterm([-1 2 2 Momvar(4)],1,1);
lmiterm([-1 2 3 Momvar(5)],1,1);
lmiterm([-1 3 3 Momvar(6)],1,1);

% Schur-complement form of ||G*m-d||^2 <= residualBound.
[residualVariable,~,~]=lmivar(1,[1 1]);
lmiterm([-2 1 1 residualVariable],1,1);
for parameterIndex=1:6
    lmiterm([-2 2 1 Momvar(parameterIndex)],Gmat(:,parameterIndex),1);
end
lmiterm([-2 2 1 0],-obs);
lmiterm([-2 2 2 0],eye(size(Gmat,1)));
lmiterm([-4 1 1 residualVariable],1,1);

% Keep mtt within the range of the observed second moments.
maxObservedMoment=max(obs);
% Constraint:
%     mtt <= maxObservedMoment
% Equivalent LMI:
%     mtt < maxObservedMoment
lmiterm([3 1 1 Momvar(1)],1,1);
lmiterm([-3 1 1 0],maxObservedMoment);


lmiSystem=getlmis;

% Find a feasible positive-semidefinite moment matrix.
[feasibilityValue,feasiblePoint]=feasp(lmiSystem);
if isempty(feasiblePoint) || feasibilityValue>1e-6
    error('Step3_Inverse:Infeasible','The LMI system is infeasible.');
end

% The seventh decision variable is the residual upper bound.
objective=zeros(7,1);
objective(7)=1;
[objectiveValue,solution]=mincx(...
    lmiSystem,objective,[1e-2,40,-100,10,0],feasiblePoint);
if isempty(solution)
    error('Step3_Inverse:OptimizationFailed','LMI optimization failed.');
end

m2=solution(1:6);
m2=m2(:);

if m2(1)<=0
    error('Step3_Inverse:InvalidDuration','The inferred mtt must be positive.');
end

% Characteristic spatial scales come from the 2-D spatial block.
spatialMoment=[m2(4),m2(5);m2(5),m2(6)];
spatialMoment=0.5*(spatialMoment+spatialMoment');
spatialEigenvalue=eig(spatialMoment);
spatialEigenvalue=max(spatialEigenvalue,0);
spatialEigenvalue=sort(spatialEigenvalue,'descend');
Lc=2*sqrt(spatialEigenvalue(1));
Wc=2*sqrt(spatialEigenvalue(2));

vx=m2(2)/m2(1);
vy=m2(3)/m2(1);
tauc=2*sqrt(m2(1));

pred=Gmat*m2;
residual=pred-obs;
info.MomentMatrix=[...
    m2(1),m2(2),m2(3); ...
    m2(2),m2(4),m2(5); ...
    m2(3),m2(5),m2(6)];
info.sMoment=spatialMoment;
info.eigen=spatialEigenvalue;
info.object=objectiveValue;
info.pred=pred;
info.res=residual;
info.Misfit=sum(residual.^2)/Energy;
info.rmsMisfit=sqrt(mean(residual.^2));
end
