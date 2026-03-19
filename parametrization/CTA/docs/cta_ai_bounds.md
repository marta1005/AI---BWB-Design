# CTA AI Design-Space Bounds

## Fixed Parameters

- `S` (`s1_deg`): `55.864059922087726` deg
- `C1` (`c1_root_chord`): `38.0` m
- `C2/C1` (`c2_c1_ratio`): `0.4394736842105263`
- `te_exact_segments`: `(0, 3)`

## Active Variables

| Parameter | Display Name | Symbol | Units | Lower | Reference | Upper |
| --- | --- | --- | --- | --- | --- | --- |
| span | Semi-span | b/2 | m | 35 | 39.5 | 45 |
| b1_span_ratio | Span segment B1 | B1/(b/2) | - | 0.18 | 0.202532 | 0.24 |
| b2_span_ratio | Span segment B2 | B2/(b/2) | - | 0.08 | 0.101266 | 0.13 |
| b3_span_ratio | Span segment B3 | B3/(b/2) | - | 0.64 | 0.696203 | 0.74 |
| s2_deg | Sweep S1 | S1 | deg | 52 | 57.1715 | 62 |
| s3_deg | Sweep S2 | S2 | deg | 28 | 30 | 34 |
| c3_c1_ratio | Section C3 chord ratio | C3/C1 | - | 0.24 | 0.276316 | 0.32 |
| c4_c1_ratio | Section C4 chord ratio | C4/C1 | - | 0.02 | 0.0295491 | 0.06 |
| twist_c1_deg | Twist at C1 | twist_C1 | deg | 0.2 | 1 | 2 |
| twist_c2_deg | Twist at C2 | twist_C2 | deg | 0.2 | 1 | 2 |
| twist_c3_deg | Twist at C3 | twist_C3 | deg | 0.1 | 0.8 | 1.8 |
| twist_c4_deg | Twist at C4 | twist_C4 | deg | 0.1 | 0.6 | 1.5 |
| c1_upper_cst_0 | C1 UPPER CST coefficient 0 | A_upper,0^C1 | - | 0.356118 | 0.375765 | 0.395411 |
| c1_upper_cst_1 | C1 UPPER CST coefficient 1 | A_upper,1^C1 | - | 0.280339 | 0.370247 | 0.460156 |
| c1_upper_cst_2 | C1 UPPER CST coefficient 2 | A_upper,2^C1 | - | 0.27605 | 0.351277 | 0.426503 |
| c1_upper_cst_3 | C1 UPPER CST coefficient 3 | A_upper,3^C1 | - | 0.139138 | 0.220184 | 0.301231 |
| c1_upper_cst_4 | C1 UPPER CST coefficient 4 | A_upper,4^C1 | - | 0.0461416 | 0.0866566 | 0.127172 |
| c1_upper_cst_5 | C1 UPPER CST coefficient 5 | A_upper,5^C1 | - | 0 | 0.032 | 0.0759049 |
| c1_lower_cst_0 | C1 LOWER CST coefficient 0 | A_lower,0^C1 | - | 0.278537 | 0.298183 | 0.317829 |
| c1_lower_cst_1 | C1 LOWER CST coefficient 1 | A_lower,1^C1 | - | 0.189899 | 0.279808 | 0.369716 |
| c1_lower_cst_2 | C1 LOWER CST coefficient 2 | A_lower,2^C1 | - | 0.194662 | 0.269889 | 0.345116 |
| c1_lower_cst_3 | C1 LOWER CST coefficient 3 | A_lower,3^C1 | - | 0.106385 | 0.187432 | 0.268478 |
| c1_lower_cst_4 | C1 LOWER CST coefficient 4 | A_lower,4^C1 | - | 0.0442607 | 0.0847757 | 0.125291 |
| c1_lower_cst_5 | C1 LOWER CST coefficient 5 | A_lower,5^C1 | - | 0 | 0.026 | 0.0699049 |
| c2_upper_cst_0 | C2 UPPER CST coefficient 0 | A_upper,0^C2 | - | 0.339659 | 0.359306 | 0.378952 |
| c2_upper_cst_1 | C2 UPPER CST coefficient 1 | A_upper,1^C2 | - | 0.261981 | 0.35189 | 0.441798 |
| c2_upper_cst_2 | C2 UPPER CST coefficient 2 | A_upper,2^C2 | - | 0.263833 | 0.33906 | 0.414287 |
| c2_upper_cst_3 | C2 UPPER CST coefficient 3 | A_upper,3^C2 | - | 0.137933 | 0.21898 | 0.300027 |
| c2_upper_cst_4 | C2 UPPER CST coefficient 4 | A_upper,4^C2 | - | 0.0481062 | 0.0886212 | 0.129136 |
| c2_upper_cst_5 | C2 UPPER CST coefficient 5 | A_upper,5^C2 | - | 0 | 0.03 | 0.0739049 |
| c2_lower_cst_0 | C2 LOWER CST coefficient 0 | A_lower,0^C2 | - | 0.265176 | 0.284822 | 0.304469 |
| c2_lower_cst_1 | C2 LOWER CST coefficient 1 | A_lower,1^C2 | - | 0.175721 | 0.265629 | 0.355538 |
| c2_lower_cst_2 | C2 LOWER CST coefficient 2 | A_lower,2^C2 | - | 0.181313 | 0.256539 | 0.331766 |
| c2_lower_cst_3 | C2 LOWER CST coefficient 3 | A_lower,3^C2 | - | 0.0991903 | 0.180237 | 0.261284 |
| c2_lower_cst_4 | C2 LOWER CST coefficient 4 | A_lower,4^C2 | - | 0.0424083 | 0.0829233 | 0.123438 |
| c2_lower_cst_5 | C2 LOWER CST coefficient 5 | A_lower,5^C2 | - | 0 | 0.024 | 0.0679049 |
| c3_upper_cst_0 | C3 UPPER CST coefficient 0 | A_upper,0^C3 | - | 0.333882 | 0.349695 | 0.365508 |
| c3_upper_cst_1 | C3 UPPER CST coefficient 1 | A_upper,1^C3 | - | 0.255449 | 0.340506 | 0.425563 |
| c3_upper_cst_2 | C3 UPPER CST coefficient 2 | A_upper,2^C3 | - | 0.258868 | 0.333009 | 0.407149 |
| c3_upper_cst_3 | C3 UPPER CST coefficient 3 | A_upper,3^C3 | - | 0.142786 | 0.221465 | 0.300145 |
| c3_upper_cst_4 | C3 UPPER CST coefficient 4 | A_upper,4^C3 | - | 0.0520489 | 0.0922196 | 0.13239 |
| c3_upper_cst_5 | C3 UPPER CST coefficient 5 | A_upper,5^C3 | - | 0 | 0.028 | 0.0684087 |
| c3_lower_cst_0 | C3 LOWER CST coefficient 0 | A_lower,0^C3 | - | 0.25966 | 0.275473 | 0.291286 |
| c3_lower_cst_1 | C3 LOWER CST coefficient 1 | A_lower,1^C3 | - | 0.17105 | 0.256106 | 0.341163 |
| c3_lower_cst_2 | C3 LOWER CST coefficient 2 | A_lower,2^C3 | - | 0.175626 | 0.249767 | 0.323907 |
| c3_lower_cst_3 | C3 LOWER CST coefficient 3 | A_lower,3^C3 | - | 0.0985226 | 0.177202 | 0.255882 |
| c3_lower_cst_4 | C3 LOWER CST coefficient 4 | A_lower,4^C3 | - | 0.0415742 | 0.0817449 | 0.121916 |
| c3_lower_cst_5 | C3 LOWER CST coefficient 5 | A_lower,5^C3 | - | 0 | 0.022 | 0.0624087 |
| c4_upper_cst_0 | C4 UPPER CST coefficient 0 | A_upper,0^C4 | - | 0.333882 | 0.349695 | 0.365508 |
| c4_upper_cst_1 | C4 UPPER CST coefficient 1 | A_upper,1^C4 | - | 0.255449 | 0.340506 | 0.425563 |
| c4_upper_cst_2 | C4 UPPER CST coefficient 2 | A_upper,2^C4 | - | 0.258868 | 0.333009 | 0.407149 |
| c4_upper_cst_3 | C4 UPPER CST coefficient 3 | A_upper,3^C4 | - | 0.142786 | 0.221465 | 0.300145 |
| c4_upper_cst_4 | C4 UPPER CST coefficient 4 | A_upper,4^C4 | - | 0.0520489 | 0.0922196 | 0.13239 |
| c4_upper_cst_5 | C4 UPPER CST coefficient 5 | A_upper,5^C4 | - | 0 | 0.028 | 0.0684087 |
| c4_lower_cst_0 | C4 LOWER CST coefficient 0 | A_lower,0^C4 | - | 0.25966 | 0.275473 | 0.291286 |
| c4_lower_cst_1 | C4 LOWER CST coefficient 1 | A_lower,1^C4 | - | 0.17105 | 0.256106 | 0.341163 |
| c4_lower_cst_2 | C4 LOWER CST coefficient 2 | A_lower,2^C4 | - | 0.175626 | 0.249767 | 0.323907 |
| c4_lower_cst_3 | C4 LOWER CST coefficient 3 | A_lower,3^C4 | - | 0.0985226 | 0.177202 | 0.255882 |
| c4_lower_cst_4 | C4 LOWER CST coefficient 4 | A_lower,4^C4 | - | 0.0415742 | 0.0817449 | 0.121916 |
| c4_lower_cst_5 | C4 LOWER CST coefficient 5 | A_lower,5^C4 | - | 0 | 0.022 | 0.0624087 |
