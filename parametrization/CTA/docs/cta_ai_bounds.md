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
| c1_upper_cst_0 | C0 UPPER CST coefficient 0 | A_upper,0^C0 | - | 0.261481 | 0.321481 | 0.381481 |
| c1_upper_cst_1 | C0 UPPER CST coefficient 1 | A_upper,1^C0 | - | 0.075418 | 0.135418 | 0.195418 |
| c1_upper_cst_2 | C0 UPPER CST coefficient 2 | A_upper,2^C0 | - | 0.150324 | 0.200324 | 0.250324 |
| c1_upper_cst_3 | C0 UPPER CST coefficient 3 | A_upper,3^C0 | - | 0.0538325 | 0.103833 | 0.153833 |
| c1_upper_cst_4 | C0 UPPER CST coefficient 4 | A_upper,4^C0 | - | 0.176301 | 0.216301 | 0.256301 |
| c1_upper_cst_5 | C0 UPPER CST coefficient 5 | A_upper,5^C0 | - | 0.132833 | 0.162833 | 0.192833 |
| c1_lower_cst_0 | C0 LOWER CST coefficient 0 | A_lower,0^C0 | - | 0.0579597 | 0.11796 | 0.17796 |
| c1_lower_cst_1 | C0 LOWER CST coefficient 1 | A_lower,1^C0 | - | 0.178236 | 0.238236 | 0.298236 |
| c1_lower_cst_2 | C0 LOWER CST coefficient 2 | A_lower,2^C0 | - | 0 | 0.0384797 | 0.0884797 |
| c1_lower_cst_3 | C0 LOWER CST coefficient 3 | A_lower,3^C0 | - | 0.2485 | 0.2985 | 0.3485 |
| c1_lower_cst_4 | C0 LOWER CST coefficient 4 | A_lower,4^C0 | - | 0.254061 | 0.294061 | 0.334061 |
| c1_lower_cst_5 | C0 LOWER CST coefficient 5 | A_lower,5^C0 | - | 0.0524508 | 0.0824508 | 0.112451 |
| c2_upper_cst_0 | C3 UPPER CST coefficient 0 | A_upper,0^C3 | - | 0.205994 | 0.265994 | 0.325994 |
| c2_upper_cst_1 | C3 UPPER CST coefficient 1 | A_upper,1^C3 | - | 0.203918 | 0.263918 | 0.323918 |
| c2_upper_cst_2 | C3 UPPER CST coefficient 2 | A_upper,2^C3 | - | 0.232616 | 0.282616 | 0.332616 |
| c2_upper_cst_3 | C3 UPPER CST coefficient 3 | A_upper,3^C3 | - | 0.119049 | 0.169049 | 0.219049 |
| c2_upper_cst_4 | C3 UPPER CST coefficient 4 | A_upper,4^C3 | - | 0.131835 | 0.171835 | 0.211835 |
| c2_upper_cst_5 | C3 UPPER CST coefficient 5 | A_upper,5^C3 | - | 0.0376101 | 0.0676101 | 0.0976101 |
| c2_lower_cst_0 | C3 LOWER CST coefficient 0 | A_lower,0^C3 | - | 0.105435 | 0.165435 | 0.225435 |
| c2_lower_cst_1 | C3 LOWER CST coefficient 1 | A_lower,1^C3 | - | 0.156822 | 0.216822 | 0.276822 |
| c2_lower_cst_2 | C3 LOWER CST coefficient 2 | A_lower,2^C3 | - | 0.0474823 | 0.0974823 | 0.147482 |
| c2_lower_cst_3 | C3 LOWER CST coefficient 3 | A_lower,3^C3 | - | 0.198314 | 0.248314 | 0.298314 |
| c2_lower_cst_4 | C3 LOWER CST coefficient 4 | A_lower,4^C3 | - | 0.0808849 | 0.120885 | 0.160885 |
| c2_lower_cst_5 | C3 LOWER CST coefficient 5 | A_lower,5^C3 | - | 0.0699319 | 0.0999319 | 0.129932 |
| c3_upper_cst_0 | C4 UPPER CST coefficient 0 | A_upper,0^C4 | - | 0.132497 | 0.192497 | 0.252497 |
| c3_upper_cst_1 | C4 UPPER CST coefficient 1 | A_upper,1^C4 | - | 0 | 0.0546772 | 0.114677 |
| c3_upper_cst_2 | C4 UPPER CST coefficient 2 | A_upper,2^C4 | - | 0.196107 | 0.246107 | 0.296107 |
| c3_upper_cst_3 | C4 UPPER CST coefficient 3 | A_upper,3^C4 | - | 0.0920357 | 0.142036 | 0.192036 |
| c3_upper_cst_4 | C4 UPPER CST coefficient 4 | A_upper,4^C4 | - | 0.116398 | 0.156398 | 0.196398 |
| c3_upper_cst_5 | C4 UPPER CST coefficient 5 | A_upper,5^C4 | - | 0.0286404 | 0.0586404 | 0.0886404 |
| c3_lower_cst_0 | C4 LOWER CST coefficient 0 | A_lower,0^C4 | - | 0.0107784 | 0.0707784 | 0.130778 |
| c3_lower_cst_1 | C4 LOWER CST coefficient 1 | A_lower,1^C4 | - | 0.0709301 | 0.13093 | 0.19093 |
| c3_lower_cst_2 | C4 LOWER CST coefficient 2 | A_lower,2^C4 | - | 0.0249666 | 0.0749666 | 0.124967 |
| c3_lower_cst_3 | C4 LOWER CST coefficient 3 | A_lower,3^C4 | - | 0.093058 | 0.143058 | 0.193058 |
| c3_lower_cst_4 | C4 LOWER CST coefficient 4 | A_lower,4^C4 | - | 0.0961378 | 0.136138 | 0.176138 |
| c3_lower_cst_5 | C4 LOWER CST coefficient 5 | A_lower,5^C4 | - | 0 | 0.0233196 | 0.0533196 |
| c4_upper_cst_0 | C5 UPPER CST coefficient 0 | A_upper,0^C5 | - | 0 | 0.0526879 | 0.112688 |
| c4_upper_cst_1 | C5 UPPER CST coefficient 1 | A_upper,1^C5 | - | 0.0686357 | 0.128636 | 0.188636 |
| c4_upper_cst_2 | C5 UPPER CST coefficient 2 | A_upper,2^C5 | - | 0.0377204 | 0.0877204 | 0.13772 |
| c4_upper_cst_3 | C5 UPPER CST coefficient 3 | A_upper,3^C5 | - | 0.0981183 | 0.148118 | 0.198118 |
| c4_upper_cst_4 | C5 UPPER CST coefficient 4 | A_upper,4^C5 | - | 0.0879116 | 0.127912 | 0.167912 |
| c4_upper_cst_5 | C5 UPPER CST coefficient 5 | A_upper,5^C5 | - | 0 | 0.0249661 | 0.0549661 |
| c4_lower_cst_0 | C5 LOWER CST coefficient 0 | A_lower,0^C5 | - | 0.135429 | 0.195429 | 0.255429 |
| c4_lower_cst_1 | C5 LOWER CST coefficient 1 | A_lower,1^C5 | - | 0 | 0.0390611 | 0.0990611 |
| c4_lower_cst_2 | C5 LOWER CST coefficient 2 | A_lower,2^C5 | - | 0.162915 | 0.212915 | 0.262915 |
| c4_lower_cst_3 | C5 LOWER CST coefficient 3 | A_lower,3^C5 | - | 0.0647407 | 0.114741 | 0.164741 |
| c4_lower_cst_4 | C5 LOWER CST coefficient 4 | A_lower,4^C5 | - | 0.10529 | 0.14529 | 0.18529 |
| c4_lower_cst_5 | C5 LOWER CST coefficient 5 | A_lower,5^C5 | - | 0 | 0.0295286 | 0.0595286 |
