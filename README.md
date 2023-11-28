# MERMAID_waveform

Software to conduct the analysis and make the figures
in the paper **Waveform modeling of hydroacoustic teleseismic earthquake records from autonomous MERMAID floats**, by Sirawich Pipatprathanporn and Frederik J Simons

### Cited as

To be announced

Author: Sirawich Pipatprathanporn

Email:  sirawich@princeton.edu

## How to install the package

[0] Clone the repository

`git clone git@github.com:sirpipat/MERMAID_waveform.git`

[1] Install the following required dependent packages:

- [slepian_alpha](https://github.com/csdms-contrib/slepian_alpha)
- [slepian_oscar](https://github.com/csdms-contrib/slepian_oscar)
- [irisFetch](https://ds.iris.edu/ds/nodes/dmc/software/downloads/irisfetch.m/)
- [mermaid_buffer](https://github.com/sirpipat/MERMAID_buffer)
- [seizmo](https://github.com/sirpipat/seizmo)
- [TauP](https://www.seis.sc.edu/taup/) (Version 2.3.0 or later)
- [Instaseis](https://instaseis.net)
- [SPECFEM2D](https://github.com/SPECFEM/specfem2d)

[2] The following environmental variable must be set in the shell:

```
export MERMAID2=/where-you-clone-THIS-repoitory/
export MERMAID=/where-you-clone-MERMAID_Buffer-repoitory/
export ONEYEAR=/where-you-keep-the-buffer-files/
export SAC=/where-you-store-MERMAID-reports/
export IFILES=/where-you-keep-the-information-files/
export SFILES=/where-you-want-the-output-SAC-files-to-be/
export NCFILES=/where-you-keep-the-NetCDF-files/
export EPS=/where-you-want-plots-to-be-saved/
export SLEPIANS=/where-you-put-slepian_alpha-and-slepian_oscar/
export IRISFETCH=/where-you-keep-irisFetch/
export REMOTE2D=/where-you-keep-SPECFEM2D-simulations/
export SPECFEM2D=/where-you-keep-SPECFEM2D-software/
export SEIZMO=/where-you-clone-SEIZMO-repository/
export TAUP=/where-you-keep-TauP-package/
```

[3] Add the following paths `startup.m`, so that MATLAB recognizes the installed packages

```
addpath(genpath(getenv('SLEPIANS')))
addpath(genpath(getenv('MERMAID')))
addpath(genpath(getenv('MERMAID2')))
addpath(genpath(getenv('SEIZMO')))
addpath(genpath(getenv('IRISFETCH')))
javaaddpath(fullfile(getenv('TAUP'), 'lib', 'TauP-x.x.x.jar'))  % x.x.x is the TauP version
```

[4] Some of the packages may have their own `startup.m` files. To prevent MATLAB from using other startup files by accident, create another environmental variable to locate *your* startup file

```
export STARTUP=/where-you-put-your-start-up-file/
```

Then, add this line to *your* startup file

```
% ensure that running this does not cause MATLAB to use other startup file
addpath(genpath(getenv('STARTUP')))
```
