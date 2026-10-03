function result=calDerived(m2,G,d)

m2=double(m2(:));
d=double(d(:));

if numel(m2)>=7
    m2=m2(1:6);
end

if numel(m2)~=6
    error('m2 must contain six parameters.');
end

if size(G,2)~=6
    error('G must have six columns.');
end

if size(G,1)~=numel(d)
    error('G and d have incompatible dimensions.');
end

mtt=m2(1);
mxt=m2(2);
myt=m2(3);
mxx=m2(4);
mxy=m2(5);
myy=m2(6);

if ~isfinite(mtt) || mtt<=0
    error('mtt must be positive.');
end

spatialMatrix=[mxx,mxy;mxy,myy];
spatialMatrix=0.5*(spatialMatrix+spatialMatrix.');

spatialEigenvalues=eig(spatialMatrix);
spatialEigenvalues=real(spatialEigenvalues);
spatialEigenvalues=max(spatialEigenvalues,0);
spatialEigenvalues=sort(spatialEigenvalues,'descend');

Lc=2*sqrt(spatialEigenvalues(1));
Wc=2*sqrt(spatialEigenvalues(2));

tauc=2*sqrt(mtt);

Vx=mxt/mtt;
Vy=myt/mtt;
mV0=hypot(Vx,Vy);

Vc=Lc/tauc;
L0=tauc*mV0;

if Lc>0
    ratio=L0/Lc;
else
    ratio=NaN;
end

predictedMoment=G*m2;
momentResidual=predictedMoment-d;

observedTauc=2*sqrt(max(d,0));
predictedTauc=2*sqrt(max(predictedMoment,0));
taucResidual=predictedTauc-observedTauc;

if sum(d.^2)>eps
    Misfit=sum(momentResidual.^2)/sum(d.^2);
else
    Misfit=NaN;
end

RMS_moment=sqrt(mean(momentResidual.^2));
RMS_tauc=sqrt(mean(taucResidual.^2));

result=struct();

result.m2=m2;
result.Lc=Lc;
result.Wc=Wc;
result.tauc=tauc;
result.Vx=Vx;
result.Vy=Vy;
result.V0=[Vx;Vy];
result.mV0=mV0;
result.Vc=Vc;
result.L0=L0;
result.ratio=ratio;

result.spatialMatrix=spatialMatrix;
result.spatialEigenvalues=spatialEigenvalues;

result.predictedMoment=predictedMoment;
result.momentResidual=momentResidual;
result.observedTauc=observedTauc;
result.predictedTauc=predictedTauc;
result.taucResidual=taucResidual;

result.Misfit=Misfit;
result.RMS_moment=RMS_moment;
result.RMS_tauc=RMS_tauc;

end
