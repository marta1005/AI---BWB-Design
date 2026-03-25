# CTA AI Design-Space Bounds

## Fixed Parameters

- `S` (`s1_deg`): `66.87` deg
- `B1` (`b1_fixed_m`): `8.041` m
- `C1` conditioned by fixed `S` and straight `TE(C0->C1)`
- no public `C2`: the inboard `TE(C1->C3)` blend is built with a hidden helper point
- `C3` and `C4` are active chord variables
- `TE(C3->C4)` is smooth and not constrained to remain straight
- `twist` is constant from Section 0 through Section 3
- `te_exact_segments`: `(0, 4)`

## Active Variables

| Parameter | Display Name | Symbol | Units | Lower | Reference | Upper |
| --- | --- | --- | --- | --- | --- | --- |
| span | Wing span | B2+B3 | m | 30 | 31.459 | 35 |
| c1_root_chord | Body chord | C0 | m | 37 | 41.203 | 45 |
| c2_c1_ratio | Transition-wing chord | C3 | m | 13 | 13.927 | 16 |
| c3_c1_ratio | Outer-wing chord | C4 | m | 6.8 | 7.768 | 9.8 |
| b2_span_ratio | Transition wing fraction B2 | B2/(B2+B3) | - | 0.14 | 0.142026 | 0.23 |
| c4_c1_ratio | Wing tip chord | C5 | m | 0.8 | 0.8 | 1.8 |
| s2_deg | Sweep S1 | S1 | deg | 45 | 54.059 | 66 |
| s3_deg | Sweep S2 | S2 | deg | 27 | 27.71 | 40 |
| twist_c1_deg | Twist at C0/C3 | twist_C0/C3 | deg | 0.2 | 1 | 2 |
| twist_c3_deg | Twist at C4 | twist_C4 | deg | 0.1 | 0.8 | 1.5 |
| twist_c4_deg | Twist at C5 | twist_C5 | deg | 0.1 | 0.6 | 1 |
| c1_upper_cst_0 | C0 UPPER CST coefficient 0 | A_upper,0^C0 | - | 0.233948 | 0.233948 | 0.395411 |
| c1_upper_cst_1 | C0 UPPER CST coefficient 1 | A_upper,1^C0 | - | 0.280339 | 0.313004 | 0.460156 |
| c1_upper_cst_2 | C0 UPPER CST coefficient 2 | A_upper,2^C0 | - | 0.124486 | 0.124486 | 0.426503 |
| c1_upper_cst_3 | C0 UPPER CST coefficient 3 | A_upper,3^C0 | - | 0.139138 | 0.45129 | 0.45129 |
| c1_upper_cst_4 | C0 UPPER CST coefficient 4 | A_upper,4^C0 | - | 0.0461416 | 0.075379 | 0.127172 |
| c1_upper_cst_5 | C0 UPPER CST coefficient 5 | A_upper,5^C0 | - | 0 | 0.359262 | 0.359262 |
| c1_lower_cst_0 | C0 LOWER CST coefficient 0 | A_lower,0^C0 | - | 0.221796 | 0.221796 | 0.317829 |
| c1_lower_cst_1 | C0 LOWER CST coefficient 1 | A_lower,1^C0 | - | 0.164675 | 0.164675 | 0.369716 |
| c1_lower_cst_2 | C0 LOWER CST coefficient 2 | A_lower,2^C0 | - | 0.060559 | 0.060559 | 0.345116 |
| c1_lower_cst_3 | C0 LOWER CST coefficient 3 | A_lower,3^C0 | - | 0.106385 | 0.310561 | 0.310561 |
| c1_lower_cst_4 | C0 LOWER CST coefficient 4 | A_lower,4^C0 | - | -0.039576 | -0.039576 | 0.125291 |
| c1_lower_cst_5 | C0 LOWER CST coefficient 5 | A_lower,5^C0 | - | 0 | 0.224411 | 0.224411 |
| c2_upper_cst_0 | C3 UPPER CST coefficient 0 | A_upper,0^C3 | - | 0.205464 | 0.205464 | 0.378952 |
| c2_upper_cst_1 | C3 UPPER CST coefficient 1 | A_upper,1^C3 | - | 0.261981 | 0.283149 | 0.441798 |
| c2_upper_cst_2 | C3 UPPER CST coefficient 2 | A_upper,2^C3 | - | 0.112921 | 0.112921 | 0.414287 |
| c2_upper_cst_3 | C3 UPPER CST coefficient 3 | A_upper,3^C3 | - | 0.137933 | 0.403674 | 0.403674 |
| c2_upper_cst_4 | C3 UPPER CST coefficient 4 | A_upper,4^C3 | - | 0.0481062 | 0.073142 | 0.129136 |
| c2_upper_cst_5 | C3 UPPER CST coefficient 5 | A_upper,5^C3 | - | 0 | 0.322782 | 0.322782 |
| c2_lower_cst_0 | C3 LOWER CST coefficient 0 | A_lower,0^C3 | - | 0.193312 | 0.193312 | 0.304469 |
| c2_lower_cst_1 | C3 LOWER CST coefficient 1 | A_lower,1^C3 | - | 0.13482 | 0.13482 | 0.355538 |
| c2_lower_cst_2 | C3 LOWER CST coefficient 2 | A_lower,2^C3 | - | 0.048994 | 0.048994 | 0.331766 |
| c2_lower_cst_3 | C3 LOWER CST coefficient 3 | A_lower,3^C3 | - | 0.0991903 | 0.262945 | 0.262945 |
| c2_lower_cst_4 | C3 LOWER CST coefficient 4 | A_lower,4^C3 | - | -0.041814 | -0.041814 | 0.123438 |
| c2_lower_cst_5 | C3 LOWER CST coefficient 5 | A_lower,5^C3 | - | 0 | 0.187931 | 0.187931 |
| c3_upper_cst_0 | C4 UPPER CST coefficient 0 | A_upper,0^C4 | - | 0.17698 | 0.17698 | 0.365508 |
| c3_upper_cst_1 | C4 UPPER CST coefficient 1 | A_upper,1^C4 | - | 0.253294 | 0.253294 | 0.425563 |
| c3_upper_cst_2 | C4 UPPER CST coefficient 2 | A_upper,2^C4 | - | 0.101356 | 0.101356 | 0.407149 |
| c3_upper_cst_3 | C4 UPPER CST coefficient 3 | A_upper,3^C4 | - | 0.142786 | 0.356058 | 0.356058 |
| c3_upper_cst_4 | C4 UPPER CST coefficient 4 | A_upper,4^C4 | - | 0.0520489 | 0.070904 | 0.13239 |
| c3_upper_cst_5 | C4 UPPER CST coefficient 5 | A_upper,5^C4 | - | 0 | 0.286303 | 0.286303 |
| c3_lower_cst_0 | C4 LOWER CST coefficient 0 | A_lower,0^C4 | - | 0.164828 | 0.164828 | 0.291286 |
| c3_lower_cst_1 | C4 LOWER CST coefficient 1 | A_lower,1^C4 | - | 0.104965 | 0.104965 | 0.341163 |
| c3_lower_cst_2 | C4 LOWER CST coefficient 2 | A_lower,2^C4 | - | 0.037429 | 0.037429 | 0.323907 |
| c3_lower_cst_3 | C4 LOWER CST coefficient 3 | A_lower,3^C4 | - | 0.0985226 | 0.21533 | 0.255882 |
| c3_lower_cst_4 | C4 LOWER CST coefficient 4 | A_lower,4^C4 | - | -0.044051 | -0.044051 | 0.121916 |
| c3_lower_cst_5 | C4 LOWER CST coefficient 5 | A_lower,5^C4 | - | 0 | 0.151452 | 0.151452 |
| c4_upper_cst_0 | C5 UPPER CST coefficient 0 | A_upper,0^C5 | - | 0.17698 | 0.17698 | 0.365508 |
| c4_upper_cst_1 | C5 UPPER CST coefficient 1 | A_upper,1^C5 | - | 0.253294 | 0.253294 | 0.425563 |
| c4_upper_cst_2 | C5 UPPER CST coefficient 2 | A_upper,2^C5 | - | 0.101356 | 0.101356 | 0.407149 |
| c4_upper_cst_3 | C5 UPPER CST coefficient 3 | A_upper,3^C5 | - | 0.142786 | 0.356058 | 0.356058 |
| c4_upper_cst_4 | C5 UPPER CST coefficient 4 | A_upper,4^C5 | - | 0.0520489 | 0.070904 | 0.13239 |
| c4_upper_cst_5 | C5 UPPER CST coefficient 5 | A_upper,5^C5 | - | 0 | 0.286303 | 0.286303 |
| c4_lower_cst_0 | C5 LOWER CST coefficient 0 | A_lower,0^C5 | - | 0.164828 | 0.164828 | 0.291286 |
| c4_lower_cst_1 | C5 LOWER CST coefficient 1 | A_lower,1^C5 | - | 0.104965 | 0.104965 | 0.341163 |
| c4_lower_cst_2 | C5 LOWER CST coefficient 2 | A_lower,2^C5 | - | 0.037429 | 0.037429 | 0.323907 |
| c4_lower_cst_3 | C5 LOWER CST coefficient 3 | A_lower,3^C5 | - | 0.0985226 | 0.21533 | 0.255882 |
| c4_lower_cst_4 | C5 LOWER CST coefficient 4 | A_lower,4^C5 | - | -0.044051 | -0.044051 | 0.121916 |
| c4_lower_cst_5 | C5 LOWER CST coefficient 5 | A_lower,5^C5 | - | 0 | 0.151452 | 0.151452 |
