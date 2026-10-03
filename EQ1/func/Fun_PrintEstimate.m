% =========================================================================
% Print inversion results
fprintf('============================================================\n');
fprintf('        SECOND-MOMENT INVERSION RESULTS\n');

fprintf('\nModel parameters:\n');
fprintf('  mtt = % .3e\n',m2(1));
fprintf('  mxt = % .3e\n',m2(2));
fprintf('  myt = % .3e\n',m2(3));
fprintf('  mxx = % .3e\n',m2(4));
fprintf('  mxy = % .3e\n',m2(5));
fprintf('  myy = % .3e\n',m2(6));
fprintf('  [% .3e  % .3e  % .3e]\n',m2(1),m2(2),m2(3));
fprintf('  [% .3e  % .3e  % .3e]\n',m2(2),m2(4),m2(5));
fprintf('  [% .3e  % .3e  % .3e]\n',m2(3),m2(5),m2(6));

fprintf('\nCharacteristic dimensions:\n');
fprintf('  Lc   = % .3f\n',Lc);
fprintf('  Wc   = % .3f\n',Wc);
fprintf('  tauc = % .3f\n',tauc);

fprintf('\nCentroid rupture velocity:\n');
fprintf('  Vx   = % .3f\n',V0(1));
fprintf('  Vy   = % .3f\n',V0(2));
fprintf('  |V0| = % .3f\n',mV0);

fprintf('\nDirectivity-related quantities:\n');
fprintf('  Vc    = % .3f\n',Vc);
fprintf('  L0    = % .3f\n',L0);
fprintf('  ratio = % .3f\n',ratio);

fprintf('\nMisfit statistics:\n');
fprintf('  Moment rel.misfit    = %.3e\n',Misfit);
fprintf('  Tauc rel.misfit      = %.3e\n',Misfit_tauc1);
