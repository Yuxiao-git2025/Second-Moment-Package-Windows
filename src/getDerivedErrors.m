function errors=getDerivedErrors(m2,covariance)

m2=double(m2(:));
covariance=double(covariance);

if numel(m2)>=7
    m2=m2(1:6);
end

if numel(m2)~=6
    error('m2 must contain six parameters.');
end

if ~isequal(size(covariance),[6,6])
    error('Covariance matrix must be 6 by 6.');
end

% Symmetrize covariance.
covariance=0.5*(covariance+covariance.');

% Remove tiny asymmetric numerical errors.
covariance(abs(covariance)<1e-14*max(1,max(abs(covariance(:)))))=0;

% Derived quantities.
base=calDerived(m2,zeros(1,6),0);

% -------------------------------------------------------------------------
% Jacobian
% -------------------------------------------------------------------------

J=zeros(7,6);

% Numerical derivatives are safer here because Lc and Wc depend
% on eigenvalues and eigenvectors.
for kk=1:6

    step=1e-6*max(abs(m2(kk)),1);

    mp=m2;
    mm=m2;

    mp(kk)=mp(kk)+step;
    mm(kk)=mm(kk)-step;

    % Keep mtt positive.
    if kk==1 && mm(kk)<=0
        mm(kk)=0.5*m2(kk);
    end
    fp=calDerived(mp,zeros(1,6),0);
    fm=calDerived(mm,zeros(1,6),0);

    J(:,kk)=...
        ([fp.tauc;fp.Lc;fp.Wc;fp.Vx;fp.Vy;fp.mV0;fp.ratio]-...
         [fm.tauc;fm.Lc;fm.Wc;fm.Vx;fm.Vy;fm.mV0;fm.ratio])...
        /(mp(kk)-mm(kk));
end

derivedCovariance=J*covariance*J.';
derivedVariance=real(diag(derivedCovariance));
derivedVariance=max(derivedVariance,0);

stds=sqrt(derivedVariance);

errors=struct();

errors.sigmaTauc=stds(1);
errors.sigmaLc=stds(2);
errors.sigmaWc=stds(3);
errors.sigmaVx=stds(4);
errors.sigmaVy=stds(5);
errors.sigmaMv0=stds(6);
errors.sigmaRatio=stds(7);

errors.jacobian=J;
errors.covariance=derivedCovariance;

end
