# CTA AI Design-Space Bounds

## Fixed Parameters

- `S` (`s1_deg`): `45.07529007428979` deg
- `B1` (`b1_fixed_m`): `8.041` m
- `C1` conditioned by fixed `S` and straight `TE(C0->C1)`
- no public `C2`: the inboard `TE(C1->C3)` blend is built with a hidden helper point
- `C3` is an active chord variable and `C4` is driven by the transition taper ratio `C4/C3`
- `TE(C3->C4)` is smooth and not constrained to remain straight
- `twist` is constant from Section 0 through Section 3
- `te_exact_segments`: `(0, 4)`

## Active Variables

| Parameter | Display Name | Symbol | Units | Lower | CTA | Upper |
| --- | --- | --- | --- | --- | --- | --- |
| span | Wing span | B2+B3 | m | 28 | 31.4585 | 35 |
| c1_root_chord | Body chord | C0 | m | 39 | 41.1795 | 43 |
| c2_c1_ratio | Transition-wing chord | C3 | m | 13 | 13.927 | 16 |
| c4_c3_ratio | Transition taper ratio | C4/C3 | - | 0.45 | 0.5578 | 0.6 |
| b2_span_ratio | Transition wing fraction B2 | B2/(B2+B3) | - | 0.13 | 0.142 | 0.21 |
| c4_c1_ratio | Wing tip chord | C5 | m | 0.8 | 0.8 | 1.8 |
| s2_deg | Transition-wing 50% chord sweep | S1 | deg | 15 | 34.6 | 45 |
| s3_deg | Outer-wing 25% chord sweep | S2 | deg | 22 | 24.7 | 33 |
| med_3_te_sweep_deg | Transition-wing trailing-edge sweep | med_3_TEswp | - | -10 | 0 | 25 |
| twist_c1_deg | Twist at C0/C3 | twist_C0/C3 | deg | 0.2 | 1 | 2 |
| twist_c3_deg | Twist at C4 | twist_C4 | deg | 0.1 | 0.8 | 1.5 |
| twist_c4_deg | Twist at C5 | twist_C5 | deg | 0.1 | 0.6 | 1 |
| c1_upper_cst_0 | C0 UPPER CST coefficient 0 | A_upper,0^C0 | - | 0.12182 | 0.18182 | 0.24182 |
| c1_upper_cst_1 | C0 UPPER CST coefficient 1 | A_upper,1^C0 | - | 0.110994 | 0.170994 | 0.230994 |
| c1_upper_cst_2 | C0 UPPER CST coefficient 2 | A_upper,2^C0 | - | 0.0955489 | 0.145549 | 0.195549 |
| c1_upper_cst_3 | C0 UPPER CST coefficient 3 | A_upper,3^C0 | - | 0.0946557 | 0.144656 | 0.194656 |
| c1_upper_cst_4 | C0 UPPER CST coefficient 4 | A_upper,4^C0 | - | 0.145836 | 0.185836 | 0.225836 |
| c1_upper_cst_5 | C0 UPPER CST coefficient 5 | A_upper,5^C0 | - | 0.216883 | 0.246883 | 0.276883 |
| c1_lower_cst_0 | C0 LOWER CST coefficient 0 | A_lower,0^C0 | - | 0.0302235 | 0.0902235 | 0.150224 |
| c1_lower_cst_1 | C0 LOWER CST coefficient 1 | A_lower,1^C0 | - | 0.0638262 | 0.123826 | 0.183826 |
| c1_lower_cst_2 | C0 LOWER CST coefficient 2 | A_lower,2^C0 | - | 0.105916 | 0.155916 | 0.205916 |
| c1_lower_cst_3 | C0 LOWER CST coefficient 3 | A_lower,3^C0 | - | 0.167536 | 0.217536 | 0.267536 |
| c1_lower_cst_4 | C0 LOWER CST coefficient 4 | A_lower,4^C0 | - | 0.223505 | 0.263505 | 0.303505 |
| c1_lower_cst_5 | C0 LOWER CST coefficient 5 | A_lower,5^C0 | - | 0.236405 | 0.266405 | 0.296405 |
| c2_upper_cst_0 | C3 UPPER CST coefficient 0 | A_upper,0^C3 | - | 0.0925355 | 0.152535 | 0.212535 |
| c2_upper_cst_1 | C3 UPPER CST coefficient 1 | A_upper,1^C3 | - | 0.178737 | 0.238737 | 0.298737 |
| c2_upper_cst_2 | C3 UPPER CST coefficient 2 | A_upper,2^C3 | - | 0.198247 | 0.248247 | 0.298247 |
| c2_upper_cst_3 | C3 UPPER CST coefficient 3 | A_upper,3^C3 | - | 0.156737 | 0.206737 | 0.256737 |
| c2_upper_cst_4 | C3 UPPER CST coefficient 4 | A_upper,4^C3 | - | 0.120324 | 0.160324 | 0.200324 |
| c2_upper_cst_5 | C3 UPPER CST coefficient 5 | A_upper,5^C3 | - | 0.088422 | 0.118422 | 0.148422 |
| c2_lower_cst_0 | C3 LOWER CST coefficient 0 | A_lower,0^C3 | - | 0.0514626 | 0.111463 | 0.171463 |
| c2_lower_cst_1 | C3 LOWER CST coefficient 1 | A_lower,1^C3 | - | 0.0815296 | 0.14153 | 0.20153 |
| c2_lower_cst_2 | C3 LOWER CST coefficient 2 | A_lower,2^C3 | - | 0.112304 | 0.162304 | 0.212304 |
| c2_lower_cst_3 | C3 LOWER CST coefficient 3 | A_lower,3^C3 | - | 0.130106 | 0.180106 | 0.230106 |
| c2_lower_cst_4 | C3 LOWER CST coefficient 4 | A_lower,4^C3 | - | 0.131876 | 0.171876 | 0.211876 |
| c2_lower_cst_5 | C3 LOWER CST coefficient 5 | A_lower,5^C3 | - | 0.11133 | 0.14133 | 0.17133 |
| c3_upper_cst_0 | C4 UPPER CST coefficient 0 | A_upper,0^C4 | - | 0.0319011 | 0.0919011 | 0.151901 |
| c3_upper_cst_1 | C4 UPPER CST coefficient 1 | A_upper,1^C4 | - | 0.0583346 | 0.118335 | 0.178335 |
| c3_upper_cst_2 | C4 UPPER CST coefficient 2 | A_upper,2^C4 | - | 0.106522 | 0.156522 | 0.206522 |
| c3_upper_cst_3 | C4 UPPER CST coefficient 3 | A_upper,3^C4 | - | 0.125268 | 0.175268 | 0.225268 |
| c3_upper_cst_4 | C4 UPPER CST coefficient 4 | A_upper,4^C4 | - | 0.116635 | 0.156635 | 0.196635 |
| c3_upper_cst_5 | C4 UPPER CST coefficient 5 | A_upper,5^C4 | - | 0.0860065 | 0.116006 | 0.146006 |
| c3_lower_cst_0 | C4 LOWER CST coefficient 0 | A_lower,0^C4 | - | 0 | 0.0505488 | 0.110549 |
| c3_lower_cst_1 | C4 LOWER CST coefficient 1 | A_lower,1^C4 | - | 0.0184803 | 0.0784803 | 0.13848 |
| c3_lower_cst_2 | C4 LOWER CST coefficient 2 | A_lower,2^C4 | - | 0.0583884 | 0.108388 | 0.158388 |
| c3_lower_cst_3 | C4 LOWER CST coefficient 3 | A_lower,3^C4 | - | 0.0791256 | 0.129126 | 0.179126 |
| c3_lower_cst_4 | C4 LOWER CST coefficient 4 | A_lower,4^C4 | - | 0.0811981 | 0.121198 | 0.161198 |
| c3_lower_cst_5 | C4 LOWER CST coefficient 5 | A_lower,5^C4 | - | 0.061361 | 0.091361 | 0.121361 |
| c4_upper_cst_0 | C5 UPPER CST coefficient 0 | A_upper,0^C5 | - | 0 | 0.0380992 | 0.0980992 |
| c4_upper_cst_1 | C5 UPPER CST coefficient 1 | A_upper,1^C5 | - | 0.0178194 | 0.0778194 | 0.137819 |
| c4_upper_cst_2 | C5 UPPER CST coefficient 2 | A_upper,2^C5 | - | 0.0647012 | 0.114701 | 0.164701 |
| c4_upper_cst_3 | C5 UPPER CST coefficient 3 | A_upper,3^C5 | - | 0.0845848 | 0.134585 | 0.184585 |
| c4_upper_cst_4 | C5 UPPER CST coefficient 4 | A_upper,4^C5 | - | 0.0809784 | 0.120978 | 0.160978 |
| c4_upper_cst_5 | C5 UPPER CST coefficient 5 | A_upper,5^C5 | - | 0.0546797 | 0.0846797 | 0.11468 |
| c4_lower_cst_0 | C5 LOWER CST coefficient 0 | A_lower,0^C5 | - | 0.0351757 | 0.0951757 | 0.155176 |
| c4_lower_cst_1 | C5 LOWER CST coefficient 1 | A_lower,1^C5 | - | 0.044252 | 0.104252 | 0.164252 |
| c4_lower_cst_2 | C5 LOWER CST coefficient 2 | A_lower,2^C5 | - | 0.0815437 | 0.131544 | 0.181544 |
| c4_lower_cst_3 | C5 LOWER CST coefficient 3 | A_lower,3^C5 | - | 0.0986002 | 0.1486 | 0.1986 |
| c4_lower_cst_4 | C5 LOWER CST coefficient 4 | A_lower,4^C5 | - | 0.09199 | 0.13199 | 0.17199 |
| c4_lower_cst_5 | C5 LOWER CST coefficient 5 | A_lower,5^C5 | - | 0.0622349 | 0.0922349 | 0.122235 |
