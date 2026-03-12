| Group | Parameter | Display Name | Symbol | Units | Normalization | Lower Bound | Reference | Upper Bound | Description |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| topology | span | Semi-span | b/2 | m | absolute | 20 | 39.5 | 60 | Half-span of the BWB planform. |
| topology | b1_span_ratio | Span segment B1 | B1/(b/2) | - | fraction of semi-span | 0.1 | 0.18 | 0.25 | Inboard span partition ratio. |
| topology | b2_span_ratio | Span segment B2 | B2/(b/2) | - | fraction of semi-span | 0.05 | 0.12 | 0.25 | Transition span partition ratio. |
| topology | b3_span_ratio | Span segment B3 | B3/(b/2) | - | fraction of semi-span | 0.45 | 0.7 | 0.8 | Outboard span partition ratio. |
| planform | le_root_x | Root LE x-position | x_LE,root | m | absolute | -5 | 0 | 8 | Leading-edge x-location at the root section. |
| planform | c1_root_chord | Root chord | C1 | m | absolute | 20 | 40 | 60 | Chord at the root control section. |
| planform | c2_c1_ratio | Section C2 chord ratio | C2/C1 | - | ratio to C1 | 0.55 | 0.7 | 0.85 | Chord ratio at section C2 relative to the root chord. |
| planform | c3_c1_ratio | Section C3 chord ratio | C3/C1 | - | ratio to C1 | 0.18 | 0.23 | 0.28 | Chord ratio at section C3 relative to the root chord. |
| planform | c4_c1_ratio | Section C4 chord ratio | C4/C1 | - | ratio to C1 | 0.06 | 0.08 | 0.09 | Chord ratio at section C4 relative to the root chord. |
| planform | s1_deg | Sweep segment S1 | S1 | deg | angle | 40 | 58 | 60 | Leading-edge sweep of segment S1. |
| planform | s2_deg | Sweep segment S2 | S2 | deg | angle | 40 | 50 | 60 | Leading-edge sweep of segment S2. |
| planform | s3_deg | Sweep segment S3 | S3 | deg | angle | 24 | 32 | 40 | Leading-edge sweep of segment S3. |
| planform | nose_blend_y | Root nose blend length | y_blend,nose | m | absolute | 0.5 | 2.5 | 6 | Spanwise blending length used to round the root nose. |
| class_function | cst_n1 | CST class exponent N1 | N1 | - | class-function exponent | 0.45 | 0.5 | 0.55 | Leading-edge exponent of the Kulfan class function. |
| class_function | cst_n2 | CST class exponent N2 | N2 | - | class-function exponent | 0.95 | 1 | 1.1 | Trailing-edge exponent of the Kulfan class function. |
| attitude | dihedral_deg | Global dihedral | Gamma | deg | angle | 0 | 0 | 12 | Global dihedral angle of the wing. |
| twist | twist_c1_deg | Twist at C1 | twist_C1 | deg | angle | -10 | 0 | 10 | Geometric twist at section C1. |
| twist | twist_c2_deg | Twist at C2 | twist_C2 | deg | angle | -10 | 0 | 10 | Geometric twist at section C2. |
| twist | twist_c3_deg | Twist at C3 | twist_C3 | deg | angle | -12 | -1 | 8 | Geometric twist at section C3. |
| twist | twist_c4_deg | Twist at C4 | twist_C4 | deg | angle | -15 | -3 | 10 | Geometric twist at section C4. |
| camber_mode | camber_c1 | Camber delta at C1 | dcamber_C1 | - | mode amplitude | -0.03 | 0 | 0.03 | Additional camber mode amplitude at section C1. |
| camber_mode | camber_c2 | Camber delta at C2 | dcamber_C2 | - | mode amplitude | -0.03 | 0 | 0.03 | Additional camber mode amplitude at section C2. |
| camber_mode | camber_c3 | Camber delta at C3 | dcamber_C3 | - | mode amplitude | -0.03 | 0 | 0.03 | Additional camber mode amplitude at section C3. |
| camber_mode | camber_c4 | Camber delta at C4 | dcamber_C4 | - | mode amplitude | -0.03 | 0 | 0.03 | Additional camber mode amplitude at section C4. |
| thickness_targets | c1_tc_max | Maximum thickness at C1 | (t/c)_max,C1 | - | t/c | 0.14 | 0.22 | 0.28 | Maximum thickness ratio target at section C1. |
| thickness_targets | c2_tc_max | Maximum thickness at C2 | (t/c)_max,C2 | - | t/c | 0.12 | 0.18 | 0.24 | Maximum thickness ratio target at section C2. |
| thickness_targets | c3_tc_max | Maximum thickness at C3 | (t/c)_max,C3 | - | t/c | 0.08 | 0.11 | 0.18 | Maximum thickness ratio target at section C3. |
| thickness_targets | c4_tc_max | Maximum thickness at C4 | (t/c)_max,C4 | - | t/c | 0.05 | 0.08 | 0.12 | Maximum thickness ratio target at section C4. |
| thickness_targets | c1_x_tmax | Thickness peak location at C1 | x_tmax,C1/c | - | x/c | 0.22 | 0.33 | 0.48 | Chordwise position of the maximum thickness at section C1. |
| thickness_targets | c2_x_tmax | Thickness peak location at C2 | x_tmax,C2/c | - | x/c | 0.2 | 0.31 | 0.45 | Chordwise position of the maximum thickness at section C2. |
| thickness_targets | c3_x_tmax | Thickness peak location at C3 | x_tmax,C3/c | - | x/c | 0.18 | 0.28 | 0.42 | Chordwise position of the maximum thickness at section C3. |
| thickness_targets | c4_x_tmax | Thickness peak location at C4 | x_tmax,C4/c | - | x/c | 0.14 | 0.24 | 0.35 | Chordwise position of the maximum thickness at section C4. |
| thickness_targets | c1_te_thickness | Trailing-edge thickness at C1 | t_TE,C1/c | - | t_TE/c_local | 0 | 0.002 | 0.01 | Trailing-edge thickness ratio at section C1. |
| thickness_targets | c2_te_thickness | Trailing-edge thickness at C2 | t_TE,C2/c | - | t_TE/c_local | 0 | 0.002 | 0.01 | Trailing-edge thickness ratio at section C2. |
| thickness_targets | c3_te_thickness | Trailing-edge thickness at C3 | t_TE,C3/c | - | t_TE/c_local | 0 | 0.0015 | 0.008 | Trailing-edge thickness ratio at section C3. |
| thickness_targets | c4_te_thickness | Trailing-edge thickness at C4 | t_TE,C4/c | - | t_TE/c_local | 0 | 0.001 | 0.006 | Trailing-edge thickness ratio at section C4. |
| cst_c1 | c1_upper_cst_0 | C1 UPPER CST coefficient 0 | A_upper,0^C1 | - | Bernstein coefficient | 0 | 0.35057 | 0.876426 | Kulfan CST Bernstein coefficient 0 on the upper surface of section C1. |
| cst_c1 | c1_upper_cst_1 | C1 UPPER CST coefficient 1 | A_upper,1^C1 | - | Bernstein coefficient | 0 | 0.36338 | 0.908451 | Kulfan CST Bernstein coefficient 1 on the upper surface of section C1. |
| cst_c1 | c1_upper_cst_2 | C1 UPPER CST coefficient 2 | A_upper,2^C1 | - | Bernstein coefficient | 0 | 0.333057 | 0.832643 | Kulfan CST Bernstein coefficient 2 on the upper surface of section C1. |
| cst_c1 | c1_upper_cst_3 | C1 UPPER CST coefficient 3 | A_upper,3^C1 | - | Bernstein coefficient | 0 | 0.20379 | 0.509475 | Kulfan CST Bernstein coefficient 3 on the upper surface of section C1. |
| cst_c1 | c1_upper_cst_4 | C1 UPPER CST coefficient 4 | A_upper,4^C1 | - | Bernstein coefficient | 0 | 0.0809508 | 0.45 | Kulfan CST Bernstein coefficient 4 on the upper surface of section C1. |
| cst_c1 | c1_upper_cst_5 | C1 UPPER CST coefficient 5 | A_upper,5^C1 | - | Bernstein coefficient | 0 | 0.05 | 0.45 | Kulfan CST Bernstein coefficient 5 on the upper surface of section C1. |
| cst_c1 | c1_lower_cst_0 | C1 LOWER CST coefficient 0 | A_lower,0^C1 | - | Bernstein coefficient | 0 | 0.140598 | 0.45 | Kulfan CST Bernstein coefficient 0 on the lower surface of section C1. |
| cst_c1 | c1_lower_cst_1 | C1 LOWER CST coefficient 1 | A_lower,1^C1 | - | Bernstein coefficient | 0 | 0.10518 | 0.45 | Kulfan CST Bernstein coefficient 1 on the lower surface of section C1. |
| cst_c1 | c1_lower_cst_2 | C1 LOWER CST coefficient 2 | A_lower,2^C1 | - | Bernstein coefficient | 0 | 0.0966736 | 0.45 | Kulfan CST Bernstein coefficient 2 on the lower surface of section C1. |
| cst_c1 | c1_lower_cst_3 | C1 LOWER CST coefficient 3 | A_lower,3^C1 | - | Bernstein coefficient | 0 | 0.0836959 | 0.45 | Kulfan CST Bernstein coefficient 3 on the lower surface of section C1. |
| cst_c1 | c1_lower_cst_4 | C1 LOWER CST coefficient 4 | A_lower,4^C1 | - | Bernstein coefficient | 0 | 0.0545472 | 0.45 | Kulfan CST Bernstein coefficient 4 on the lower surface of section C1. |
| cst_c1 | c1_lower_cst_5 | C1 LOWER CST coefficient 5 | A_lower,5^C1 | - | Bernstein coefficient | 0 | 0.01 | 0.45 | Kulfan CST Bernstein coefficient 5 on the lower surface of section C1. |
| cst_c2 | c2_upper_cst_0 | C2 UPPER CST coefficient 0 | A_upper,0^C2 | - | Bernstein coefficient | 0 | 0.274463 | 0.686158 | Kulfan CST Bernstein coefficient 0 on the upper surface of section C2. |
| cst_c2 | c2_upper_cst_1 | C2 UPPER CST coefficient 1 | A_upper,1^C2 | - | Bernstein coefficient | 0 | 0.2865 | 0.71625 | Kulfan CST Bernstein coefficient 1 on the upper surface of section C2. |
| cst_c2 | c2_upper_cst_2 | C2 UPPER CST coefficient 2 | A_upper,2^C2 | - | Bernstein coefficient | 0 | 0.273903 | 0.684758 | Kulfan CST Bernstein coefficient 2 on the upper surface of section C2. |
| cst_c2 | c2_upper_cst_3 | C2 UPPER CST coefficient 3 | A_upper,3^C2 | - | Bernstein coefficient | 0 | 0.170874 | 0.45 | Kulfan CST Bernstein coefficient 3 on the upper surface of section C2. |
| cst_c2 | c2_upper_cst_4 | C2 UPPER CST coefficient 4 | A_upper,4^C2 | - | Bernstein coefficient | 0 | 0.0667862 | 0.45 | Kulfan CST Bernstein coefficient 4 on the upper surface of section C2. |
| cst_c2 | c2_upper_cst_5 | C2 UPPER CST coefficient 5 | A_upper,5^C2 | - | Bernstein coefficient | 0 | 0.04 | 0.45 | Kulfan CST Bernstein coefficient 5 on the upper surface of section C2. |
| cst_c2 | c2_lower_cst_0 | C2 LOWER CST coefficient 0 | A_lower,0^C2 | - | Bernstein coefficient | 0 | 0.117786 | 0.45 | Kulfan CST Bernstein coefficient 0 on the lower surface of section C2. |
| cst_c2 | c2_lower_cst_1 | C2 LOWER CST coefficient 1 | A_lower,1^C2 | - | Bernstein coefficient | 0 | 0.087386 | 0.45 | Kulfan CST Bernstein coefficient 1 on the lower surface of section C2. |
| cst_c2 | c2_lower_cst_2 | C2 LOWER CST coefficient 2 | A_lower,2^C2 | - | Bernstein coefficient | 0 | 0.0807407 | 0.45 | Kulfan CST Bernstein coefficient 2 on the lower surface of section C2. |
| cst_c2 | c2_lower_cst_3 | C2 LOWER CST coefficient 3 | A_lower,3^C2 | - | Bernstein coefficient | 0 | 0.0699906 | 0.45 | Kulfan CST Bernstein coefficient 3 on the lower surface of section C2. |
| cst_c2 | c2_lower_cst_4 | C2 LOWER CST coefficient 4 | A_lower,4^C2 | - | Bernstein coefficient | 0 | 0.0446773 | 0.45 | Kulfan CST Bernstein coefficient 4 on the lower surface of section C2. |
| cst_c2 | c2_lower_cst_5 | C2 LOWER CST coefficient 5 | A_lower,5^C2 | - | Bernstein coefficient | 0 | 0.008 | 0.45 | Kulfan CST Bernstein coefficient 5 on the lower surface of section C2. |
| cst_c3 | c3_upper_cst_0 | C3 UPPER CST coefficient 0 | A_upper,0^C3 | - | Bernstein coefficient | 0 | 0.177117 | 0.45 | Kulfan CST Bernstein coefficient 0 on the upper surface of section C3. |
| cst_c3 | c3_upper_cst_1 | C3 UPPER CST coefficient 1 | A_upper,1^C3 | - | Bernstein coefficient | 0 | 0.178902 | 0.45 | Kulfan CST Bernstein coefficient 1 on the upper surface of section C3. |
| cst_c3 | c3_upper_cst_2 | C3 UPPER CST coefficient 2 | A_upper,2^C3 | - | Bernstein coefficient | 0 | 0.176151 | 0.45 | Kulfan CST Bernstein coefficient 2 on the upper surface of section C3. |
| cst_c3 | c3_upper_cst_3 | C3 UPPER CST coefficient 3 | A_upper,3^C3 | - | Bernstein coefficient | 0 | 0.112291 | 0.45 | Kulfan CST Bernstein coefficient 3 on the upper surface of section C3. |
| cst_c3 | c3_upper_cst_4 | C3 UPPER CST coefficient 4 | A_upper,4^C3 | - | Bernstein coefficient | 0 | 0.0436269 | 0.45 | Kulfan CST Bernstein coefficient 4 on the upper surface of section C3. |
| cst_c3 | c3_upper_cst_5 | C3 UPPER CST coefficient 5 | A_upper,5^C3 | - | Bernstein coefficient | 0 | 0.025 | 0.45 | Kulfan CST Bernstein coefficient 5 on the upper surface of section C3. |
| cst_c3 | c3_lower_cst_0 | C3 LOWER CST coefficient 0 | A_lower,0^C3 | - | Bernstein coefficient | 0 | 0.0779496 | 0.45 | Kulfan CST Bernstein coefficient 0 on the lower surface of section C3. |
| cst_c3 | c3_lower_cst_1 | C3 LOWER CST coefficient 1 | A_lower,1^C3 | - | Bernstein coefficient | 0 | 0.0547743 | 0.45 | Kulfan CST Bernstein coefficient 1 on the lower surface of section C3. |
| cst_c3 | c3_lower_cst_2 | C3 LOWER CST coefficient 2 | A_lower,2^C3 | - | Bernstein coefficient | 0 | 0.0497903 | 0.45 | Kulfan CST Bernstein coefficient 2 on the lower surface of section C3. |
| cst_c3 | c3_lower_cst_3 | C3 LOWER CST coefficient 3 | A_lower,3^C3 | - | Bernstein coefficient | 0 | 0.0440508 | 0.45 | Kulfan CST Bernstein coefficient 3 on the lower surface of section C3. |
| cst_c3 | c3_lower_cst_4 | C3 LOWER CST coefficient 4 | A_lower,4^C3 | - | Bernstein coefficient | 0 | 0.027903 | 0.45 | Kulfan CST Bernstein coefficient 4 on the lower surface of section C3. |
| cst_c3 | c3_lower_cst_5 | C3 LOWER CST coefficient 5 | A_lower,5^C3 | - | Bernstein coefficient | 0 | 0.005 | 0.45 | Kulfan CST Bernstein coefficient 5 on the lower surface of section C3. |
| cst_c4 | c4_upper_cst_0 | C4 UPPER CST coefficient 0 | A_upper,0^C4 | - | Bernstein coefficient | 0 | 0.132041 | 0.45 | Kulfan CST Bernstein coefficient 0 on the upper surface of section C4. |
| cst_c4 | c4_upper_cst_1 | C4 UPPER CST coefficient 1 | A_upper,1^C4 | - | Bernstein coefficient | 0 | 0.131983 | 0.45 | Kulfan CST Bernstein coefficient 1 on the upper surface of section C4. |
| cst_c4 | c4_upper_cst_2 | C4 UPPER CST coefficient 2 | A_upper,2^C4 | - | Bernstein coefficient | 0 | 0.136616 | 0.45 | Kulfan CST Bernstein coefficient 2 on the upper surface of section C4. |
| cst_c4 | c4_upper_cst_3 | C4 UPPER CST coefficient 3 | A_upper,3^C4 | - | Bernstein coefficient | 0 | 0.0896436 | 0.45 | Kulfan CST Bernstein coefficient 3 on the upper surface of section C4. |
| cst_c4 | c4_upper_cst_4 | C4 UPPER CST coefficient 4 | A_upper,4^C4 | - | Bernstein coefficient | 0 | 0.0342289 | 0.45 | Kulfan CST Bernstein coefficient 4 on the upper surface of section C4. |
| cst_c4 | c4_upper_cst_5 | C4 UPPER CST coefficient 5 | A_upper,5^C4 | - | Bernstein coefficient | 0 | 0.018 | 0.45 | Kulfan CST Bernstein coefficient 5 on the upper surface of section C4. |
| cst_c4 | c4_lower_cst_0 | C4 LOWER CST coefficient 0 | A_lower,0^C4 | - | Bernstein coefficient | 0 | 0.0560104 | 0.45 | Kulfan CST Bernstein coefficient 0 on the lower surface of section C4. |
| cst_c4 | c4_lower_cst_1 | C4 LOWER CST coefficient 1 | A_lower,1^C4 | - | Bernstein coefficient | 0 | 0.0367369 | 0.45 | Kulfan CST Bernstein coefficient 1 on the lower surface of section C4. |
| cst_c4 | c4_lower_cst_2 | C4 LOWER CST coefficient 2 | A_lower,2^C4 | - | Bernstein coefficient | 0 | 0.0313786 | 0.45 | Kulfan CST Bernstein coefficient 2 on the lower surface of section C4. |
| cst_c4 | c4_lower_cst_3 | C4 LOWER CST coefficient 3 | A_lower,3^C4 | - | Bernstein coefficient | 0 | 0.0285627 | 0.45 | Kulfan CST Bernstein coefficient 3 on the lower surface of section C4. |
| cst_c4 | c4_lower_cst_4 | C4 LOWER CST coefficient 4 | A_lower,4^C4 | - | Bernstein coefficient | 0 | 0.0191649 | 0.45 | Kulfan CST Bernstein coefficient 4 on the lower surface of section C4. |
| cst_c4 | c4_lower_cst_5 | C4 LOWER CST coefficient 5 | A_lower,5^C4 | - | Bernstein coefficient | 0 | 0.003 | 0.45 | Kulfan CST Bernstein coefficient 5 on the lower surface of section C4. |
