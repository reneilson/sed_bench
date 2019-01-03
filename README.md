# sed_bench
Benchmark for Simulated Electrical Disturbance  

In the literature, many papers about problems in Energy Power Quality (EPQ) use simulated data since is hard to find benchmarks with real disturbances. For this reason, some authors work with disturbance emulations using electrical devices (e.g. motors, solenoids), that is more close to the real anomalies. Moreover, to make a fair comparison with different approaches, for example in classification process, could be interesting have a benchmark (complete and standard) which could be used in all works, however, this is not available yet. In this context, this paper develops a benchmark of simulated electrical disturbances and provides it for the research academical community and developers. We named it as SED-Bench, and it contains nine (9) simple disturbances (interruption, sag, swell, oscillation, impulse,  flicker, notching, dc-level and harmonics), and eight (8) double disturbances, four (4) triple disturbances and the pure sine. Every disturbance is available in different datasets with fundamental frequency of 50Hz or 60Hz and SNR of 10dB, 20dB, 30dB, 50dB or without noise.

The data is divided as:
-bench_test
--60Hz
---Simple
---Multiple
--50Hz
---Simple
---Multiple
-bench_train
--60Hz
---Simple
---Multiple
--50Hz
---Simple
---Multiple

In each directory with data in .npy type of file there is data with 10dB with termination bc10, 20dB with termination bc20, 30dB with bc30, 50dB with bc50 and without noise with bc0. 

Each bc in the test folder has 200 data of each disturbance type and in the train folder has 1000 data of each disturbance type.

To summarize, the dataset is composed by:
• Simple Disturbances Without Noise (60Hz)
• Simple Disturbances With Noise (60Hz) (10dB, 20dB, 30dB, 50dB)
• Simple Disturbances Without Noise (50Hz)
• Simple Disturbances With Noise (50Hz) (10dB, 20dB, 30dB, 50dB)
• Multiple Disturbances Without Noise (60Hz)
• Multiple Disturbances With Noise (60Hz) (10dB, 20dB, 30dB, 50dB)
• Multiple Disturbances Without Noise (50Hz)
• Multiple Disturbances With Noise (50Hz) (10dB, 20dB, 30dB, 50dB)

Simple Disturbances:
-Sine
-DC Level
-Oscillation
-Impulsive Transient
-Flicker
-Sag
-Swell
-Harmonics
-Interruption
-Notching

Multiple Disturbances:
(Double)
-Sag + Harm
-Swell + Harm
-Sag + Notching
-Sag + Notching
-Sag + Osc
-Swell + Osc
-Sag + Imp
-Swell + Imp
(Triple)
-Sag + Harm + Osc
-Swell + Harm + Osc
-Sag + Harm + Notch
-Swell + Harm + Notch

