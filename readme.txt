
The date at which the reproducibility package was assembled
- October 19th 2024


The author(s) of the reproducibility package and their contact information: 
- Audrone Virbickaite, audrone.virbickaite@cunef.edu; audrone.virbickaite@gmail.com


The structure of the repository provided (e.g., if separating code, input data and output data):
- Codes must be run in consecutive order (00_, 01_, 02_ ...). 
- Data is loaded from 'data' folder
- Estimation output is saved in 'temp' folder then in later codes called upon to be read
- Figures and Tables are saved in the 'tables_and_figures' folder


The computing environment, language(s), licence(s) and package(s) necessary to run the reproducibility check (as well as their version); If additional information is needed to emulate the necessary environment (e.g., with conda), it should also be provided.
- R version 4.4.1 (2024-06-14 ucrt) -- "Race for Your Life"
- Packages: moments, xtable, coda, Rfast, mvnfast, zoo, truncnorm, mvtnorm, tidyquant, matrixcalc,
future.apply, plyr, xts, NMOF, Rglpk, MCS


The data being used and its format. Any relevant information regarding access to the data, origin, pre-processing, usage restrictions, etc. is to be provided.
- Cleaned and processed data is saved in .Rdata format. The original 1-minute data has been downloaded from the http://www.histdata.com/ website in csv format. 
Original files are not included as the total size is 5.67 GB. Data pre-processing details are available in the Online Appendix; after cleaning and pre-processing 
realized covariances have been calculated by using 10-minute sampling with 2-minute subsampling. Cleaned and processed data is saved in realized_measures.Rdata.
DGS10 file contains daily risk-free interest rate downloaded from FRED website: https://fred.stlouisfed.org/series/DGS10


Which code to run to produce specific tables and figures shown in the paper.
- 07_reproduce_tables_and_figures


The type of computer that was used for running the experiment and the expected runtime.
- Smaller functions were run on a personal computer 13th Gen Intel(R) Core(TM) i7-1360P 2.20 GHz, 32.0 GB (31.6 GB usable)
- 02_estimate_all_individual_copula_models and 05_portfolio codes were run on a server:   
Number of CPU(s) 0-63, GenuineIntel; Intel(R) Xeon(R) Gold 6326 CPU @ 2.90GHz; 131145156 kB memory
- Large codes took around a day (24 hours) to run on a Server


Any special setup that is required to run the reproducibility check (e.g., GPUs, parallel computing).
- By default some codes are calling for 4 core computing. 


For each sharable dataset, mention whether it is directly included in the replication kit or available elsewhere (repository, website).
- All datasets are included in the replication kit.


Some scripts create intermediary datasets from the raw data. As a default, include them in the replication kit. Please indicate which data files are intermediary and which part of the code generates them. If one cannot reproduce intermediary data files from the raw data because of bugs, time constraints, or insufficient CPU, the rest of the verification can still be carried out.
- Everything in the temp folder is intermediary. 
