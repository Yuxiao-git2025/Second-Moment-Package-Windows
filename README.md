# Second-Moment-Package-Windows
A second-moment method is used to rapidly estimate rupture information for small and moderate earthquakes, employing deconvolution of EGFs to obtain ASTF.
---
```markdown
# Second-Moment Source Analysis

A MATLAB toolbox for measuring apparent source time functions, estimating station-based temporal second moments, performing constrained second-moment inversion, and evaluating uncertainty using bootstrap and azimuthal jackknife resampling.

The workflow is designed for earthquake source analysis using mainshock and empirical Green's function (EGF) waveform pairs. The final inversion estimates a second-moment tensor-like model and derives characteristic rupture quantities, including:

- Characteristic source duration
- Characteristic source length
- Characteristic source width
- Apparent centroid velocity
- Mean centroid velocity
- Directivity ratio

---

## 1. Overview

The complete processing sequence is:

```text
Input seismograms or example data
            |
            v
EGF-based apparent source time function measurement
            |
            v
Station-based temporal second moments
            |
            v
Ray tracing and Green's-function construction
            |
            v
Semidefinite constrained inversion
            |
            v
Second-moment model parameters
            |
            v
Derived source quantities
            |
            v
Bootstrap and azimuthal jackknife uncertainty analysis
```

The inversion is based on the linear relationship

\[
d = Gm_2,
\]

where:

- \(d\) is the vector of measured temporal second moments;
- \(G\) is the inverse operator or partial-derivative matrix;
- \(m_2\) is the six-component second-moment model vector.

The model vector is written as

\[
m_2 =
\begin{bmatrix}
m_{tt} &
m_{xt} &
m_{yt} &
m_{xx} &
m_{xy} &
m_{yy}
\end{bmatrix}^{T}.
\]

The implementation may internally use an additional dummy variable during the constrained optimization. The final physical model contains the first six elements.

---

# 2. Requirements

## 2.1 MATLAB

The workflow requires MATLAB with support for:

- Matrix operations
- Singular-value decomposition
- Numerical optimization routines used by `Step3_Inverse`
- Plotting functions
- Optional `tiledlayout` and `nexttile` plotting functions

The code should be run from the project root directory or from a MATLAB path that contains all source folders.

## 2.2 External or project-specific functions

The following functions must be available:

```text
Step1_Loading
Step2_Measure
Step3_Initial
Step3_Inverse
Step4_Bootstrap
Step4_Jackknife
calDerived
distaz
Fun_PlotASTFMap
Fun_PlotEstimate
Fun_PrintEstimate
```

The ray-tracing part also requires the TOPP executable and the associated velocity-model input files.

## 2.3 Input files

The default workflow expects files such as:

```text
Input/Data.mat
Input/MEASUREMENTS_Exmaple.mat
Input/Velocal.mat
```

The filename `MEASUREMENTS_Exmaple.mat` follows the current project naming convention. If the intended spelling is `MEASUREMENTS_Example.mat`, the filename and the corresponding `load` command should be changed consistently.

---

# 3. Main Control Parameters

The main MATLAB script begins with a set of control parameters.

```matlab
isData=1;
isMeasure=0;
isInt=0;
isplot=1;
isInverse=1;
niter=80;
```

## 3.1 Data input mode

```matlab
isData=1;
```

When `isData=1`, the script loads a prepared example data file:

```matlab
load('Input\Data.mat');
```

When `isData=0`, the script executes:

```matlab
Step1_Loading;
```

This mode is intended for reading and preparing Miniseed or other waveform data.

## 3.2 ASTF measurement mode

```matlab
isMeasure=0;
```

When `isMeasure=1`, the script calculates EGF-based apparent source time functions station by station:

```matlab
[t2,DONE,STF,GFsv,dhatsv,datasv,Tsv,T1sv,...
    epsv,epldsv,tpldsv,t0,t1,PhaseSv] = ...
    Step2_Measure( ...
        velEGF,velMS,niter,npMS,npEGF,dtsv,...
        stasm,compm,isInt);
```

When `isMeasure=0`, previously prepared measurements are loaded from:

```matlab
load('Input\MEASUREMENTS_Exmaple.mat');
```

The measurement stage produces the station-dependent temporal second moment vector `t2`.

## 3.3 Interactive integration mode

```matlab
isInt=0;
```

This parameter controls whether the integration limits of the apparent source time function are selected interactively.

- `isInt=0`: integrate over the default or automatically selected ASTF interval;
- `isInt=1`: perform an additional interactive picking procedure.

## 3.4 Plotting

```matlab
isplot=1;
```

When enabled, the script generates a station-summary map:

```matlab
Fun_PlotASTFMap(IJ,slon,slat,t2,STF,dtsv,t0,lone,late);
```

where `IJ` contains the stations for which `DONE==1`.

## 3.5 Inversion

```matlab
isInverse=1;
```

When enabled, the script performs:

1. Station-data assembly;
2. Velocity-model loading;
3. Ray tracing;
4. Construction of the inverse matrix;
5. Constrained inversion;
6. Calculation of derived quantities.

The parameter

```matlab
niter=80;
```

controls the number of iterations used during the ASTF measurement stage. Values between approximately 20 and 100 are normally appropriate, depending on the data quality and waveform complexity.

---

# 4. Data Preparation and Measurement

## 4.1 Effective station selection

After the measurement stage, only completed measurements are selected:

```matlab
IJ=find(DONE==1);
```

The following station-related arrays are then assembled:

```matlab
mlats=slat(IJ);
mlons=slon(IJ);
melevs=zeros(Nj,1);
d=t2(IJ);
phas=upper(char(PhaseSv(IJ)));
```

The variables are:

| Variable | Description |
|---|---|
| `mlats` | Station latitudes |
| `mlons` | Station longitudes |
| `melevs` | Station elevations, currently set to zero |
| `d` | Measured temporal second moments |
| `phas` | Selected seismic phase, for example `P` or `S` |
| `Nj` | Number of effective stations |

The data vector is converted to a row vector before calling the main inversion:

```matlab
d=d(:)';
```

The sampling-analysis functions subsequently standardize it back to a column vector.

## 4.2 Apparent source time functions

The measurement stage estimates the apparent source time function for each station using the mainshock and EGF waveforms. The important output is:

```matlab
t2
```

This quantity is interpreted as the measured temporal second moment, commonly written as:

\[
\mu^{(0,2)}.
\]

The source-time-function measurement stage may also return:

- STF waveforms;
- Truncated and aligned EGF waveforms;
- Fitted waveforms;
- Misfit curves;
- Duration estimates;
- Phase information;
- Integration limits.

---

# 5. Velocity Model and Ray Tracing

The inversion requires a velocity model. The current script loads:

```matlab
load('Input\Velocal.mat','Vp','topl');
Vs=Vp/1.73;
```

where:

- `Vp` is the P-wave velocity model;
- `Vs` is estimated using a fixed \(V_p/V_s\) ratio;
- `topl` contains the layer-top or model-depth information.

The source and fault geometry are specified by:

```matlab
strike=strike1;
dip=dip1;
```

The ray-tracing and inverse-operator construction are performed by:

```matlab
[G,takeoffs,taketimes]=Step3_Initial( ...
    mlats,mlons,melevs,late,lone,depe,...
    Vp,Vs,topl,phas,strike,dip);
```

The outputs are:

| Output | Description |
|---|---|
| `G` | Inverse operator or partial-derivative matrix |
| `takeoffs` | Ray take-off angles |
| `taketimes` | Ray travel times |

The matrix `G` must have six columns corresponding to the six model parameters in `m2`.

---

# 6. Semidefinite Constrained Inversion

The inversion is performed using:

```matlab
[m2,~,~,Vx,Vy,~,info]=Step3_Inverse(G,d');
```

The first output is the model vector. Depending on the implementation, additional outputs may include:

- Apparent centroid velocity components;
- Optimization status;
- Constraint information;
- Diagnostic quantities.

The final physical model is extracted using:

```matlab
m2=m2(1:6);
```

The inversion is based on a constrained least-squares problem. The spatial part of the model is represented by:

```matlab
Xvar=[m2(4),m2(5);
      m2(5),m2(6)];
```

The spatial matrix is expected to satisfy the physical constraints required by the semidefinite formulation.

---

# 7. Derived Source Quantities

The preferred way to calculate derived parameters is through:

```matlab
best=calDerived(m2,G,d);
```

The main derived quantities are defined as follows.

## 7.1 Characteristic duration

\[
\tau_c = 2\sqrt{m_{tt}}.
\]

In MATLAB:

```matlab
tauc=2*sqrt(m2(1));
```

## 7.2 Characteristic length and width

The spatial second-moment matrix is:

\[
X =
\begin{bmatrix}
m_{xx} & m_{xy}\\
m_{xy} & m_{yy}
\end{bmatrix}.
\]

Its singular values are calculated using:

```matlab
[U,S,V]=svd(Xvar);
```

The characteristic length and width are:

\[
L_c = 2\sqrt{\lambda_1},
\]

\[
W_c = 2\sqrt{\lambda_2},
\]

where \(\lambda_1\) and \(\lambda_2\) are the ordered singular values of the spatial matrix.

## 7.3 Apparent centroid velocity

The velocity components are calculated from:

\[
V_x = \frac{m_{xt}}{m_{tt}},
\]

\[
V_y = \frac{m_{yt}}{m_{tt}}.
\]

In MATLAB:

```matlab
V0=m2(2:3)/m2(1);
Vx=V0(1);
Vy=V0(2);
```

## 7.4 Mean centroid velocity

The magnitude of the apparent centroid velocity is:

\[
|V_0| = \sqrt{V_x^2+V_y^2}.
\]

In MATLAB:

```matlab
mV0=sqrt(sum(V0.^2));
```

## 7.5 Characteristic rupture velocity

The characteristic velocity associated with the spatial and temporal dimensions is:

\[
V_c = \frac{L_c}{\tau_c}.
\]

## 7.6 Centroid displacement scale

The centroid displacement scale is:

\[
L_0 = \tau_c |V_0|.
\]

## 7.7 Directivity ratio

The directivity ratio is:

\[
\frac{L_0}{L_c}
=
\frac{\tau_c |V_0|}{L_c}
=
\frac{|V_0|}{V_c}.
\]

It is interpreted as a normalized measure of unilateral or directed rupture behavior. Values close to zero indicate weak apparent directivity, while larger values indicate stronger directed centroid motion.

---

# 8. Sampling Uncertainty Analysis

The uncertainty analysis is controlled by:

```matlab
rng(2026);

isJack=true;
azband=20;

isBoot=true;
Nsample=500;
bconf=0.95;
```

The random seed allows the bootstrap calculation to be reproduced.

Before sampling analysis, the inputs are standardized:

```matlab
G=double(G);
d=double(d(:));
m2=double(m2(:));

if numel(m2)>=7
    m2=m2(1:6);
end

mlons=double(mlons(:));
mlats=double(mlats(:));
```

The best-fit derived values are calculated once:

```matlab
best=calDerived(m2,G,d);
```

---

# 9. Bootstrap Analysis

The bootstrap calculation is called using:

```matlab
if isBoot
    [mv0u,mv0l,bound2u,bound2l,...
        Lcu,Lcl,taucu,taucl,boot]=...
        Step4_Bootstrap( ...
            G,d,bconf,Nsample,...
            'MakeFigure',true,...
            'FigureNumber',5);
end
```

## 9.1 Bootstrap procedure

For each realization:

1. Randomly resample the observations with replacement;
2. Construct a resampled matrix `Gi`;
3. Construct a resampled data vector `di`;
4. Recalculate the constrained inversion;
5. Recalculate all derived quantities;
6. Store the result;
7. Remove failed realizations;
8. Calculate percentile confidence intervals.

The resampling is performed at the observation level. Therefore, the method preserves the original association between each row of `G` and the corresponding entry of `d`.

## 9.2 Bootstrap parameters

| Parameter | Description |
|---|---|
| `Nsample` | Number of bootstrap realizations |
| `bconf` | Confidence level, normally `0.95` |
| `MakeFigure` | Whether to display bootstrap histograms |
| `FigureNumber` | MATLAB figure number |

For example:

```matlab
Nsample=500;
bconf=0.95;
```

A larger value of `Nsample` produces smoother percentile estimates but increases computation time.

## 9.3 Bootstrap output

The structure `boot` contains the bootstrap realizations and confidence intervals.

Typical realization fields include:

```matlab
boot.m2
boot.tauc
boot.Lc
boot.Wc
boot.Vx
boot.Vy
boot.mV0
boot.Vc
boot.L0
boot.ratio
boot.bound2
boot.minvr
boot.Misfit
```

The field:

```matlab
boot.Nvalid
```

contains the number of successful bootstrap inversions.

Typical confidence-interval fields include:

```matlab
boot.taucCI
boot.LcCI
boot.WcCI
boot.VxCI
boot.VyCI
boot.mV0CI
boot.VcCI
boot.L0CI
boot.ratioCI
boot.bound2CI
boot.minvrCI
```

The confidence intervals are percentile intervals. They should always be displayed in the order:

```text
lower bound, upper bound
```

For example:

```text
Vy : [-1.61, 1.07]
```

If a confidence interval is calculated using quantiles, the implementation should enforce:

```matlab
limits=sort(limits(:));
```

This is particularly important for signed quantities such as `Vy`.

---

# 10. Azimuthal Jackknife Analysis

The azimuthal jackknife is called using:

```matlab
if isJack
    jack=Step4_Jackknife( ...
        G,d,m2,...
        late,lone,mlats,mlons,azband,...
        'MakeFigure',false,...
        'FigureNumber',6);
end
```

## 10.1 Jackknife procedure

The station azimuths are calculated from the source to each station. The full azimuth range is divided into groups of width `azband`.

For each azimuth group:

1. Delete all stations within that azimuth interval;
2. Retain the remaining stations;
3. Recalculate the constrained inversion;
4. Calculate derived source quantities;
5. Store the result;
6. Mark failed inversions as invalid;
7. Calculate jackknife uncertainties from the valid realizations.

For example:

```matlab
azband=20;
```

divides the 360-degree azimuth range into approximately 18 groups.

## 10.2 Jackknife parameters

| Parameter | Description |
|---|---|
| `azband` | Width of each deleted azimuth interval in degrees |
| `late` | Event latitude |
| `lone` | Event longitude |
| `mlats` | Station latitudes |
| `mlons` | Station longitudes |
| `MinObservations` | Minimum number of remaining observations |

The jackknife function should use `NaN` for failed or skipped realizations. Preallocating failed realizations with zeros would incorrectly classify them as valid solutions.

## 10.3 Jackknife output

The structure `jack` contains:

```matlab
jack.m2
jack.m2All
jack.tauc
jack.Lc
jack.Wc
jack.Vx
jack.Vy
jack.mV0
jack.ratio
jack.usedCount
jack.deletedCount
jack.azimuthLower
jack.azimuthUpper
jack.azimuthCenter
jack.valid
jack.Nvalid
jack.meanModel
jack.covarianceSample
jack.covarianceJK
jack.errors
```

The fields `jack.tauc`, `jack.Lc`, `jack.Vy`, and the other derived quantities contain the valid jackknife realizations.

The direct jackknife standard errors are available through:

```matlab
jack.errors.sigmaTauc
jack.errors.sigmaLc
jack.errors.sigmaWc
jack.errors.sigmaVx
jack.errors.sigmaVy
jack.errors.sigmaMv0
jack.errors.sigmaRatio
```

The classical delete-group jackknife covariance matrix is:

\[
C_{\mathrm{JK}}
=
\frac{K-1}{K}
\sum_{i=1}^{K}
(\theta_i-\bar{\theta})
(\theta_i-\bar{\theta})^T,
\]

where \(K\) is the number of valid jackknife realizations.

---

# 11. Interpreting Large Uncertainties

Large uncertainties in `Vy` and `L0/Lc` do not necessarily indicate a programming error.

The parameter

\[
V_y=\frac{m_{yt}}{m_{tt}}
\]

is sensitive to changes in both \(m_{yt}\) and \(m_{tt}\). It becomes particularly unstable when:

- The station azimuth distribution is uneven;
- One or more azimuth sectors strongly constrain the solution;
- The inverse matrix is poorly conditioned;
- The temporal second-moment measurements contain substantial scatter;
- The spatial matrix has a small second eigenvalue;
- Removing a small station group significantly changes the inversion.

The ratio

\[
\frac{L_0}{L_c}
=
\frac{\tau_c |V_0|}{L_c}
\]

is a nonlinear ratio and can become unstable when \(L_c\) varies strongly between jackknife realizations.

Therefore, large jackknife uncertainties should be interpreted as evidence of weak or uneven data constraints rather than automatically treated as a coding failure.

Recommended diagnostic quantities include:

```matlab
min(jack.Vy)
max(jack.Vy)

min(jack.Lc)
max(jack.Lc)

min(jack.ratio)
max(jack.ratio)

jack.Nvalid
```

A low valid-realization fraction should trigger further inspection of station distribution, data quality, matrix rank, and inversion constraints.

---

# 12. Recommended Execution Order

A complete run can be organized as follows:

```matlab
% 1. Load or prepare waveform data
Step1_Loading;

% 2. Measure ASTFs and temporal second moments
Step2_Measure;

% 3. Construct G and perform inversion
Step3_Initial;
Step3_Inverse;

% 4. Calculate best-fit derived quantities
best=calDerived(m2,G,d);

% 5. Bootstrap uncertainty
boot=Step4_Bootstrap(...);

% 6. Azimuthal jackknife uncertainty
jack=Step4_Jackknife(...);

% 7. Export or save results
save('Output/SecondMomentResults.mat',...
    'm2','best','boot','jack');
```

The current project uses the following compact sampling-analysis block:

```matlab
rng(2026);

isJack=true;
azband=20;

isBoot=true;
Nsample=500;
bconf=0.95;

G=double(G);
d=double(d(:));
m2=double(m2(:));

if numel(m2)>=7
    m2=m2(1:6);
end

mlons=double(mlons(:));
mlats=double(mlats(:));

best=calDerived(m2,G,d);

if isBoot
    [mv0u,mv0l,bound2u,bound2l,...
        Lcu,Lcl,taucu,taucl,boot]=...
        Step4_Bootstrap( ...
            G,d,bconf,Nsample,...
            'MakeFigure',true,...
            'FigureNumber',5);
end

if isJack
    jack=Step4_Jackknife( ...
        G,d,m2,late,lone,mlats,mlons,azband,...
        'MakeFigure',false,...
        'FigureNumber',6);
end
```

---

# 13. Suggested Output Reporting

A final report should include both best-fit values and uncertainty estimates.

Example:

```text
Parameter       Best estimate       Bootstrap 95% CI       Jackknife SE
-----------------------------------------------------------------------
tau_c           0.306               [0.232, 0.367]         0.031
L_c             0.399               [0.329, 1.525]         0.578
W_c             0.250               [0.046, 0.523]         0.450
V_x             0.823               [0.606, 2.233]         0.168
V_y             0.468               [-1.613, 1.066]        1.489
|V_0|           0.946               [0.729, 2.592]         0.782
L_0/L_c         0.727               [0.324, 0.985]         1.546
```

The bootstrap interval and jackknife standard error measure different aspects of uncertainty:

- Bootstrap intervals describe the distribution obtained by resampling observations with replacement.
- Jackknife errors describe the sensitivity of the solution to the deletion of azimuthal station groups.

Both are useful and should not be expected to be identical.

---

# 14. Reproducibility

The bootstrap analysis uses random resampling. To reproduce a result, set the random-number seed before calling the sampling functions:

```matlab
rng(2026);
```

If a different random seed is used, the bootstrap confidence intervals may change slightly, especially when `Nsample` is small.

For stable confidence intervals, use a sufficiently large number of bootstrap realizations, for example:

```matlab
Nsample=500;
```

or:

```matlab
Nsample=1000;
```

The computation may take longer because each bootstrap realization requires a new constrained inversion.

---

# 15. Saving Results

The main results can be saved using:

```matlab
if ~exist('Output','dir')
    mkdir('Output');
end

save('Output/SecondMomentResults.mat',...
    'm2','best','boot','jack');
```

Recommended saved variables are:

| Variable | Description |
|---|---|
| `m2` | Best-fit six-component model |
| `best` | Best-fit derived quantities |
| `boot` | Bootstrap realizations and confidence intervals |
| `jack` | Azimuthal jackknife realizations and uncertainties |
| `G` | Inverse operator |
| `d` | Station-based temporal second moments |
| `takeoffs` | Ray take-off angles |
| `taketimes` | Ray travel times |

For long-term reproducibility, also save:

- Velocity-model files;
- Fault geometry;
- Station coordinates;
- Phase assignments;
- Random-number seed;
- Main script version or Git commit hash.

---

# 16. Known Limitations

1. The velocity model currently estimates \(V_s\) using a fixed ratio:

   ```matlab
   Vs=Vp/1.73;
   ```

   This ratio may not be appropriate for all geological settings.

2. Station elevations are currently set to zero:

   ```matlab
   melevs=zeros(Nj,1);
   ```

   Topographic and elevation effects should be included when required.

3. Bootstrap uncertainty does not automatically include systematic errors in the velocity model, fault geometry, or ray tracing.

4. Jackknife results depend on the chosen azimuth-bin width.

5. Very large uncertainties in `Vy` or `L0/Lc` indicate that the corresponding quantities may not be strongly constrained by the available station distribution.

6. The `L0/Lc` ratio is nonlinear and should preferably be evaluated directly for every resampled or jackknife model rather than estimated only through a first-order covariance approximation.

---

# 17. License and Citation

If this code is used in a publication, please cite the associated methodological paper and acknowledge the data providers, velocity model, ray-tracing implementation, and empirical Green's function dataset.

A recommended citation section may be added here after the project paper and software license are finalized.

---
