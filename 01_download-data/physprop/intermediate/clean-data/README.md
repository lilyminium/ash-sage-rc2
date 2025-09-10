
# Cleaning bad ThermoML data

This directory contains cases discovered of data mis-entered into the ThermoML data.

The *.py files read XML files from the ThermoML 2020 database, and saves out incorrect records into human-readable CSVs.

**Table of contents**
<!-- vscode-markdown-toc -->
* 1. [Files](#Files)
	* 1.1. [bad-br-data.csv](#bad-br-data.csv)
	* 1.2. [bad-dimethyl-carbonate-data.csv](#bad-dimethyl-carbonate-data.csv)
	* 1.3. [bad-p-data.csv](#bad-p-data.csv)
	* 1.4. [bad-ring-data.csv](#bad-ring-data.csv)
	* 1.5. [bad-ester-data.csv](#bad-ester-data.csv)
* 2. [Other problematic data points](#Otherproblematicdatapoints)
	* 2.1. [Incorrect temperature and pressure](#Incorrecttemperatureandpressure)
	* 2.2. [Baseline enthalpies of mixing](#Baselineenthalpiesofmixing)

<!-- vscode-markdown-toc-config
	numbering=true
	autoSave=true
	/vscode-markdown-toc-config -->
<!-- /vscode-markdown-toc -->


##  1. <a name='Files'></a>Files

The data points identified below were verified by manually comparing at last one data point in ThermoML, vs. the corresponding table in the paper. In the cases below we assume the paper is correct and that data was mis-entered into ThermoML.

Data points were identified by computing the property in simulation with an OpenFF force fields and investigating outlier values. No attempt was made to track down other publications in a series. As such, this is likely an incomplete list of incorrect data points in ThermoML at present.

###  1.1. <a name='bad-br-data.csv'></a>bad-br-data.csv

**DOIs**
- 10.1016/j.jct.2009.05.006
- 10.1016/j.jct.2006.05.004
- 10.1016/j.jct.2005.07.009

Each of these papers contains binary mixtures of bromoalkanes and esters. In each of these papers, the mole fraction $chi_1$ is specified of the ester, and the bromoalkane is (1-$chi_1$). However, some data points have been erroneously entered into ThermoML where $chi_1$ refers to the bromoalkane, effectively swapping the mole fraction. These data points are saved to `bad-br-data.csv`. Please note that only the densities have been manually checked against the paper; no manual validation of the excess molar volume has been done.

###  1.2. <a name='bad-dimethyl-carbonate-data.csv'></a>bad-dimethyl-carbonate-data.csv


**DOIs**
- 10.1021/je050052+

This paper calculates enthalpies of mixing of dimethyl carbonate and various alcohols. Here, the mole fraction $chi_1$ should refer to the proportion of dimethyl carbonate in the mixture. The methanol series in ThermoML has incorrectly used $chi_1$ to refer to methanol, effectively swapping the mole fraction.

###  1.3. <a name='bad-p-data.csv'></a>bad-p-data.csv

**DOIs**
- 10.1016/j.jct.2006.08.001

This paper calculates enthalpies of mixing of tributyl phosphate with alcohols. Here, the mole fraction $chi_1$ should refer to the proportion of tributyl phosphate with alcohols. The butan-1-ol series in ThermoML incorrectly uses $chi_1$ to refer to butan-1-ol, effectively swapping the mole fraction.

###  1.4. <a name='bad-ring-data.csv'></a>bad-ring-data.csv

**DOIs**
- 10.1021/je060026r

This paper calculates enthalpies of mixing of two alcohols (1-Methoxy-2-propanol or 2-Butoxy Ethanol) with organic solvents. $chi_1$ should refer to the alcohol. The 2-butoxyethan-1-ol series in ThermoML incorrectly uses $chi_1$ to refer to benzene or cyclohexane, effectively swapping the mole fractions.

###  1.5. <a name='bad-ester-data.csv'></a>bad-ester-data.csv

**DOIs**
- 10.1016/j.jct.2006.10.008

This paper calculates enthalpies of mixing of mixtures of methyl esters with chloroalkanes. In the paper the methyl esters comprise:
- HCOOCH3 (methyl methanoate)
- CH3COOCH3 (methyl ethanoate)
- CH3CH2COOCH3  (methyl propanoate)
- CH3(CH2)2COOCH3 (methyl butanoate)
- CH3(CH2)3COOCH3 (methyl pentanoate)

In ThermoML these have been incorrectly entered as ethyl methanoate, propyl methanoate, butyl methanoate, and pentyl methanoate. Methyl methanoate is fine.

##  2. <a name='Otherproblematicdatapoints'></a>Other problematic data points

Below are listed some other problematic data points that I have not written out to file, but account for in `../clean-data.py`.

###  2.1. <a name='Incorrecttemperatureandpressure'></a>Incorrect temperature and pressure

**DOIs**
- 10.1016/j.jct.2007.12.002

The paper above measures properties above 1 MPa, but somehow this is registering as atmospheric pressure.

###  2.2. <a name='Baselineenthalpiesofmixing'></a>Baseline enthalpies of mixing

There are several data points that *look like* enthalpies of mixing of binary mixtures, but are actually baselines of ternary or quarternary mixtures. These are easily identifiable as they have enthalpies of mixing of 0 kJ/mol.
