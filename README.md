# Twophase-IC
Efficient estimation of semiparametric transformation model with interval-censored data in two-phase cohort studies
Author: Fei Gao
Date: 2023.07.17

Main program for simulations: TwoPhase_simu.cpp
  Compilation: g++ -I /folder_for_eigen/ -std=c++11 -O3  TwoPhase_simu.cpp -o TwoPhase_simu
  Job Submission: bsub -o out.txt TwoPhase_simu --nsub 1000 --nrep 1000 --hn 5 --seed 123 --r 0 --out n1000r0

Main program for data analysis: data.cpp
  Compilation: g++ -I /folder_for_eigen/ -std=c++11 -O3 -O3 data.cpp -o data
  Job submission: bsub ./data1 --hn 10 --r 1.0 --ngrid 20 --out datah10r1

Supporting programs
  TwoPhase.h: functions for simulation studies.
  TwoPhase_base.h: base functions for all programs.
  TwoPhase_discZ_int.h: functions for a special case with discrete Z and a model with an interaction of X and A. Used in HVTN504 data analysis.



