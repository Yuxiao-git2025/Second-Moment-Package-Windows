# Second-Moment-Package-Windows

A MATLAB package for rapidly estimating earthquake rupture information from mainshock and empirical Green’s function (EGF) waveform pairs. The package uses EGF deconvolution to obtain apparent source time functions (ASTFs), calculates station-based temporal second moments, performs constrained second-moment inversion, and evaluates model uncertainty using bootstrap and azimuthal jackknife resampling.

---

## 1. Overview

The complete workflow is:

```text
Input waveform data
        |
        v
EGF deconvolution and ASTF measurement
        |
        v
Temporal second moments at individual stations
        |
        v
Ray tracing and construction of the inverse matrix
        |
        v
Constrained second-moment inversion
        |
        v
Derived rupture parameters
        |
        v
Bootstrap and azimuthal jackknife uncertainty analysis
```

The inversion solves the linear relation

$$
d = Gm_2,
$$

where:

- \(d\) is the vector of measured temporal second moments;
- \(G\) is the inverse operator;
- \(m_2\) is the second-moment model vector.

The physical model vector is

$$
m_2 =
\begin{bmatrix}
m_{tt} &
m_{xt} &
m_{yt} &
m_{xx} &
m_{xy} &
m_{yy}
\end{bmatrix}^{T}.
$$

The optimization routine may internally use an additional dummy variable. After inversion, only the first six elements are retained as the physical model.

---

## 2. Requirements and Input Data

### 2.1 Software requirements

The package requires:

- MATLAB;
- MATLAB Optimization Toolbox or the optimization functions used by `Step3_Inverse`;
- TOPP ray-tracing executable;
- Required MATLAB functions included in the project.

Important project functions include:

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

The project is intended for Windows and requires the TOPP executable and its runtime libraries to be available.

### 2.2 Input files

The default example workflow uses:

```text
Input/Data.mat
Input/MEASUREMENTS_Example.mat
Input/Velocal.mat
```

The input data should contain, or generate, the following variables:

| Variable | Description |
|---|---|
| `velMS` | Mainshock velocity seismograms |
| `velEGF` | EGF velocity seismograms |
| `stasm` | Station list |
| `compm` | Component list |
| `slat` | Station latitudes |
| `slon` | Station longitudes |
| `late` | Earthquake latitude |
| `lone` | Earthquake longitude |
| `depe` | Earthquake depth |
| `strike1` | Fault strike |
| `dip1` | Fault dip |
| `Vp` | P-wave velocity model |
| `topl` | Velocity-layer boundaries |

---

## 3. Main Control Parameters

The main script uses the following control parameters:

```matlab
isData=1;
isMeasure=0;
isInt=0;
isplot=1;
isInverse=1;
niter=80;
```

### Data input

```matlab
isData=1;
```

- `1`: load prepared example data from `Input/Data.mat`;
- `0`: execute `Step1_Loading` and read waveform files.

### ASTF measurement

```matlab
isMeasure=0;
```

- `1`: perform ASTF measurement using `Step2_Measure`;
- `0`: load previously measured results from `Input/MEASUREMENTS_Example.mat`.

### Integration limits

```matlab
isInt=0;
```

- `0`: use the default ASTF integration limits;
- `1`: perform an additional interactive selection of the integration limits.

### Plotting and inversion

```matlab
isplot=1;
isInverse=1;
```

- `isplot=1` enables station-summary and inversion plots;
- `isInverse=1` performs the ray tracing and constrained inversion.

The parameter

```matlab
niter=80;
```

controls the number of iterations used during ASTF measurement. Values between approximately 20 and 100 are generally suitable.

---

## 4. Processing Workflow

### 4.1 ASTF measurement

When `isMeasure=1`, the package calls:

```matlab
[t2,DONE,STF,GFsv,dhatsv,datasv,Tsv,T1sv,...
    epsv,epldsv,tpldsv,t0,t1,PhaseSv] = ...
    Step2_Measure( ...
        velEGF,velMS,niter,npMS,npEGF,dtsv,...
        stasm,compm,isInt);
```

The most important output is:

```matlab
t2
```

which represents the temporal second moment of the measured ASTF:

$$
\mu^{(0,2)}.
$$

Only successful measurements are selected:

```matlab
IJ=find(DONE==1);
```

The selected station data are then assembled as:

```matlab
mlats=slat(IJ);
mlons=slon(IJ);
d=t2(IJ);
phas=upper(char(PhaseSv(IJ)));
```

### 4.2 Ray tracing and matrix construction

The local velocity model is loaded using:

```matlab
load('Input\Velocal.mat','Vp','topl');
Vs=Vp/1.73;
```

The inverse matrix is constructed by:

```matlab
[G,takeoffs,taketimes]=Step3_Initial( ...
    mlats,mlons,melevs,late,lone,depe,...
    Vp,Vs,topl,phas,strike,dip);
```

The outputs are:

- `G`: inverse operator;
- `takeoffs`: ray take-off angles;
- `taketimes`: ray travel times.

The matrix `G` must have six physical model columns.

### 4.3 Constrained inversion

The inversion is performed using:

```matlab
[m2,~,~,Vx,Vy,~,info]=Step3_Inverse(G,d');
```

The physical model is retained as:

```matlab
m2=m2(1:6);
```

The spatial component of the model is:

```matlab
Xvar=[m2(4),m2(5);
      m2(5),m2(6)];
```

---

## 5. Derived Rupture Quantities

Derived quantities are calculated by:

```matlab
best=calDerived(m2,G,d);
```

### Characteristic duration

The characteristic source duration is

$$
\tau_c = 2\sqrt{m_{tt}}.
$$

In the code, this quantity is stored as:

```matlab
tauc
```

### Characteristic length and width

The spatial second-moment matrix is

$$
X =
\begin{bmatrix}
m_{xx} & m_{xy} \\
m_{xy} & m_{yy}
\end{bmatrix}.
$$

Its singular values are obtained from:

```matlab
[U,S,V]=svd(X);
```

The characteristic length and width are calculated as

$$
L_c = 2\sqrt{\lambda_1},
$$

and

$$
W_c = 2\sqrt{\lambda_2},
$$

where \(\lambda_1\) and \(\lambda_2\) are the ordered singular values of \(X\).

The corresponding variables are:

```matlab
Lc
Wc
```

### Apparent centroid velocity

The two components of the apparent centroid velocity are

$$
V_x = \frac{m_{xt}}{m_{tt}},
$$

and

$$
V_y = \frac{m_{yt}}{m_{tt}}.
$$

They are stored as:

```matlab
Vx
Vy
```

The velocity vector is

$$
\mathbf{V}_0 =
\begin{bmatrix}
V_x \\
V_y
\end{bmatrix}.
$$

### Mean centroid velocity

The magnitude of the centroid velocity is

$$
|V_0| = \sqrt{V_x^2+V_y^2}.
$$

This quantity is stored as:

```matlab
mV0
```

### Characteristic velocity

The characteristic rupture velocity is

$$
V_c = \frac{L_c}{\tau_c}.
$$

It is stored as:

```matlab
Vc
```

### Centroid displacement scale

The centroid displacement scale is

$$
L_0 = \tau_c |V_0|.
$$

It is stored as:

```matlab
L0
```

### Directivity ratio

The normalized directivity ratio is

$$
\frac{L_0}{L_c}
=
\frac{\tau_c |V_0|}{L_c}
=
\frac{|V_0|}{V_c}.
$$

It is stored as:

```matlab
ratio
```

A small value generally indicates weak apparent directivity, whereas a larger value indicates stronger directed centroid motion. The ratio should be interpreted together with its uncertainty.

---

## 6. Bootstrap and Azimuthal Jackknife

The uncertainty analysis is controlled by:

```matlab
rng(2026);

isJack=true;
azband=20;

isBoot=true;
Nsample=500;
bconf=0.95;
```

The random seed ensures that bootstrap results can be reproduced.

### 6.1 Bootstrap analysis

The bootstrap function is called as:

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

For each bootstrap realization, the procedure:

1. Resamples the rows of `G` and the corresponding elements of `d`;
2. Performs a new constrained inversion;
3. Calculates all derived quantities;
4. Rejects failed realizations;
5. Calculates percentile confidence intervals.

Important input parameters are:

| Parameter | Description |
|---|---|
| `Nsample` | Number of bootstrap realizations |
| `bconf` | Confidence level, for example `0.95` |
| `MakeFigure` | Whether to create bootstrap histograms |
| `FigureNumber` | MATLAB figure number |

The output structure `boot` contains realization arrays such as:

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
boot.Misfit
```

The corresponding confidence intervals are stored in:

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
```

Confidence intervals should always be reported as:

```text
[lower bound, upper bound]
```

For example:

```text
Vy : [-1.61, 1.07]
```

The implementation should use:

```matlab
limits=sort(limits(:));
```

before printing each interval.

### 6.2 Azimuthal jackknife

The azimuthal jackknife is called as:

```matlab
if isJack
    jack=Step4_Jackknife( ...
        G,d,m2,late,lone,mlats,mlons,azband,...
        'MakeFigure',false,...
        'FigureNumber',6);
end
```

For each azimuthal interval, the function:

1. Calculates the station azimuths;
2. Deletes stations within one azimuth bin;
3. Repeats the constrained inversion;
4. Calculates the derived quantities;
5. Stores the valid realization;
6. Computes jackknife uncertainties.

For example:

```matlab
azband=20;
```

uses azimuthal deletion bins with a width of 20 degrees.

The output structure contains:

```matlab
jack.m2
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
jack.covarianceSample
jack.covarianceJK
jack.errors
```

The direct jackknife standard errors are:

```matlab
jack.errors.sigmaTauc
jack.errors.sigmaLc
jack.errors.sigmaWc
jack.errors.sigmaVx
jack.errors.sigmaVy
jack.errors.sigmaMv0
jack.errors.sigmaRatio
```

Failed realizations must be stored as `NaN`, rather than zeros, so that they are not incorrectly treated as valid inversions.

---

## 7. Running the Uncertainty Analysis

The complete sampling-analysis block is:

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

A complete result can be saved using:

```matlab
if ~exist('Output','dir')
    mkdir('Output');
end

save('Output/SecondMomentResults.mat',...
    'm2','best','boot','jack','G','d');
```

The best-fit results are stored in `best`, while the uncertainty results are stored in `boot` and `jack`.

Large uncertainties, particularly for `Vy` or `L0/Lc`, may indicate:

- Uneven station-azimuth coverage;
- Strong dependence on individual stations;
- Poorly constrained inversion parameters;
- Instability of the spatial second-moment eigenvalues;
- Nonlinear amplification caused by the ratio \(L_0/L_c\).

These uncertainties should be interpreted as information about data resolution and model stability rather than automatically as programming errors.

---

## 8. Output, Limitations, and Reproducibility

A typical result table may contain:

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

Bootstrap confidence intervals and jackknife standard errors describe different sources of uncertainty:

- Bootstrap uncertainty is estimated by resampling observations with replacement.
- Jackknife uncertainty measures sensitivity to the deletion of azimuthal station groups.

The main limitations are:

1. The current implementation estimates

   ```matlab
   Vs=Vp/1.73;
   ```

   using a fixed \(V_p/V_s\) ratio.

2. Station elevations are currently set to zero:

   ```matlab
   melevs=zeros(Nj,1);
   ```

3. Bootstrap and jackknife calculations do not automatically include systematic errors from the velocity model, fault geometry, or ray tracing.

4. Jackknife results depend on the selected `azband`.

5. The nonlinear directivity ratio \(L_0/L_c\) may have substantially larger uncertainty than the individual parameters used to calculate it.

For reproducibility, save the random seed, input files, velocity model, fault geometry, station coordinates, and the Git commit used for the calculation.

---
