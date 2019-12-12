
            set terminal png size 700,500 truecolor
            set output "bwa.samtools.stats.plot/quals2.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set multiplot
             set rmargin 0; set lmargin 0; set tmargin 0; set bmargin 0; set origin 0.1,0.1; set size 0.4,0.8
            set yrange [0:50]
            set ylabel "Quality"
            set xlabel "Cycle (fwd reads)"
            plot '-' using 1:2:3 with filledcurve lt 1 lc rgb "#cccccc" t '25-75th percentile' , '-' using 1:2 with lines lc rgb "#000000" t 'Median', '-' using 1:2 with lines lt 1 t 'Mean'
        1	36	36
2	36	36
3	36	36
4	36	36
5	36	36
6	36	36
7	36	36
8	36	36
9	36	36
10	36	36
11	36	36
12	36	36
13	36	36
14	36	36
15	36	36
16	36	36
17	36	36
18	36	36
19	36	36
20	36	36
21	36	36
22	36	36
23	36	36
24	32	36
25	32	36
26	32	36
27	32	36
28	32	36
29	32	36
30	32	36
31	27	36
32	21	36
33	21	36
34	21	36
35	21	36
36	21	36
37	14	36
38	14	36
39	14	36
40	14	36
41	14	36
42	14	36
43	14	36
44	14	36
45	14	36
46	14	36
47	14	36
48	14	36
49	14	36
50	14	36
51	14	36
52	14	36
53	14	36
54	14	36
55	14	36
56	14	36
57	14	36
58	14	36
59	14	36
60	14	36
61	14	36
62	14	36
63	14	36
64	14	36
65	14	36
66	14	36
67	14	36
68	14	36
69	14	36
70	14	36
71	14	36
72	14	36
73	14	36
74	14	36
75	14	36
76	14	36
77	14	36
78	14	36
79	14	36
80	14	36
81	14	36
82	14	36
83	14	36
84	14	36
85	14	36
86	14	36
87	14	36
88	14	36
89	14	36
90	14	36
91	14	36
92	14	36
93	14	36
94	14	36
95	14	36
96	14	36
97	14	36
98	14	36
99	14	36
100	14	36
101	14	36
102	14	36
103	14	36
104	14	36
105	14	36
106	14	32
107	14	32
108	14	32
109	14	32
110	14	32
111	14	32
112	14	32
113	14	32
114	14	32
115	14	32
116	14	32
117	14	32
118	14	32
119	14	32
120	14	32
121	14	32
122	14	32
123	14	32
124	14	32
125	14	32
126	14	32
127	14	32
128	14	32
end
1	36
2	36
3	36
4	36
5	36
6	36
7	36
8	36
9	36
10	36
11	36
12	36
13	36
14	36
15	36
16	36
17	36
18	36
19	36
20	36
21	36
22	36
23	36
24	36
25	36
26	36
27	36
28	36
29	36
30	36
31	36
32	36
33	36
34	36
35	36
36	36
37	36
38	36
39	36
40	36
41	36
42	36
43	36
44	36
45	36
46	36
47	36
48	36
49	36
50	36
51	32
52	32
53	32
54	32
55	32
56	32
57	32
58	32
59	32
60	32
61	32
62	32
63	32
64	32
65	32
66	32
67	32
68	32
69	32
70	32
71	32
72	27
73	27
74	27
75	27
76	27
77	27
78	27
79	27
80	27
81	27
82	27
83	27
84	27
85	27
86	27
87	27
88	27
89	27
90	27
91	27
92	27
93	27
94	21
95	21
96	21
97	21
98	21
99	21
100	21
101	21
102	21
103	21
104	21
105	21
106	21
107	21
108	21
109	21
110	21
111	21
112	21
113	21
114	21
115	27
116	27
117	27
118	27
119	27
120	27
121	27
122	27
123	27
124	27
125	27
126	27
127	27
128	21
end
1	35.57
2	35.57
3	35.47
4	35.39
5	35.29
6	35.19
7	35.06
8	34.97
9	34.86
10	34.76
11	34.64
12	34.52
13	34.39
14	34.27
15	34.13
16	33.98
17	33.82
18	33.66
19	33.49
20	33.33
21	33.17
22	33.01
23	32.84
24	32.66
25	32.48
26	32.33
27	32.16
28	31.99
29	31.83
30	31.67
31	31.50
32	31.33
33	31.16
34	31.01
35	30.84
36	30.66
37	30.49
38	30.35
39	30.18
40	30.03
41	29.85
42	29.72
43	29.55
44	29.42
45	29.25
46	29.12
47	28.96
48	28.83
49	28.69
50	28.58
51	28.43
52	28.29
53	28.17
54	28.04
55	27.95
56	27.83
57	27.69
58	27.57
59	27.46
60	27.35
61	27.23
62	27.15
63	27.03
64	26.93
65	26.82
66	26.73
67	26.62
68	26.53
69	26.40
70	26.33
71	26.21
72	26.12
73	26.02
74	25.95
75	25.84
76	25.73
77	25.61
78	25.54
79	25.42
80	25.36
81	25.27
82	25.20
83	25.13
84	25.09
85	25.03
86	24.97
87	24.94
88	24.91
89	24.87
90	24.85
91	24.84
92	24.82
93	24.81
94	24.81
95	24.81
96	24.81
97	24.82
98	24.81
99	24.83
100	24.85
101	24.85
102	24.87
103	24.91
104	24.92
105	24.94
106	24.86
107	24.81
108	24.72
109	24.65
110	24.63
111	24.61
112	24.65
113	24.67
114	24.66
115	24.70
116	24.75
117	24.81
118	24.84
119	24.91
120	24.95
121	25.07
122	25.10
123	25.22
124	25.29
125	25.40
126	25.45
127	25.60
128	23.94
end

                set origin 0.55,0.1
                set size 0.4,0.8
                unset ytics
                set y2tics mirror
                set yrange [0:50]
                unset ylabel
                set xlabel "Cycle (rev reads)"
                set label "bwa.samtools.stats" at screen 0.5,0.95 center noenhanced
                plot '-' using 1:2:3 with filledcurve lt 1 lc rgb "#cccccc" t '25-75th percentile' , '-' using 1:2 with lines lc rgb "#000000" t 'Median', '-' using 1:2 with lines lt 2 t 'Mean'
            1	32	32
2	32	32
3	32	32
4	32	32
5	32	32
6	36	36
7	36	36
8	36	36
9	36	36
10	36	36
11	36	36
12	36	36
13	36	36
14	36	36
15	36	36
16	36	36
17	36	36
18	36	36
19	36	36
20	36	36
21	32	36
22	32	36
23	32	36
24	32	36
25	32	36
26	32	36
27	32	36
28	32	36
29	32	36
30	27	36
31	27	36
32	21	36
33	21	36
34	21	36
35	14	36
36	14	36
37	14	36
38	14	36
39	14	36
40	14	36
41	14	36
42	14	36
43	14	36
44	14	36
45	14	36
46	14	36
47	14	36
48	14	36
49	14	36
50	14	36
51	14	36
52	14	36
53	14	36
54	14	36
55	14	36
56	14	36
57	14	36
58	14	36
59	14	36
60	14	36
61	14	36
62	14	36
63	14	36
64	14	36
65	14	36
66	14	36
67	14	36
68	14	36
69	14	36
70	14	36
71	14	36
72	14	36
73	14	36
74	14	36
75	14	36
76	14	36
77	14	36
78	14	36
79	14	36
80	14	36
81	14	36
82	14	36
83	14	36
84	14	36
85	14	36
86	14	36
87	14	36
88	14	36
89	14	36
90	14	36
91	14	36
92	14	36
93	14	36
94	14	36
95	14	36
96	14	36
97	14	36
98	14	36
99	14	36
100	14	36
101	14	36
102	14	36
103	14	36
104	14	32
105	14	32
106	14	32
107	14	32
108	14	32
109	14	32
110	14	32
111	14	32
112	14	32
113	14	32
114	14	32
115	14	32
116	14	32
117	14	32
118	14	32
119	14	32
120	14	32
121	14	32
122	14	32
123	14	32
124	14	32
125	14	32
126	14	32
127	14	32
128	14	32
129	14	32
130	14	32
131	14	32
132	14	32
133	14	32
134	14	32
135	14	32
136	14	32
137	14	32
138	14	32
139	14	32
140	14	32
141	14	32
142	14	32
143	14	32
144	14	32
145	14	32
146	14	32
147	14	32
148	14	32
149	14	32
150	14	32
151	14	32
end
1	32
2	32
3	32
4	32
5	32
6	36
7	36
8	36
9	36
10	36
11	36
12	36
13	36
14	36
15	36
16	36
17	36
18	36
19	36
20	36
21	36
22	36
23	36
24	36
25	36
26	36
27	36
28	36
29	36
30	36
31	36
32	36
33	36
34	36
35	36
36	36
37	36
38	36
39	36
40	36
41	36
42	36
43	36
44	36
45	36
46	36
47	32
48	32
49	32
50	32
51	32
52	32
53	32
54	32
55	32
56	32
57	32
58	32
59	32
60	32
61	32
62	32
63	32
64	32
65	32
66	32
67	32
68	32
69	27
70	27
71	27
72	27
73	27
74	27
75	27
76	27
77	27
78	27
79	27
80	27
81	21
82	21
83	21
84	21
85	21
86	21
87	21
88	21
89	21
90	21
91	21
92	21
93	21
94	21
95	21
96	21
97	21
98	21
99	21
100	21
101	21
102	21
103	21
104	21
105	21
106	21
107	21
108	21
109	21
110	21
111	14
112	14
113	14
114	14
115	14
116	14
117	14
118	14
119	14
120	14
121	14
122	14
123	14
124	14
125	14
126	14
127	21
128	21
129	21
130	21
131	21
132	21
133	21
134	21
135	21
136	21
137	21
138	21
139	21
140	21
141	21
142	21
143	21
144	21
145	21
146	21
147	21
148	21
149	21
150	27
151	21
end
1	30.88
2	31.03
3	31.31
4	31.46
5	31.41
6	34.77
7	34.72
8	34.68
9	34.63
10	34.61
11	34.54
12	34.49
13	34.39
14	34.30
15	34.18
16	34.06
17	33.89
18	33.76
19	33.56
20	33.41
21	33.22
22	33.04
23	32.82
24	32.66
25	32.44
26	32.25
27	32.02
28	31.86
29	31.65
30	31.47
31	31.26
32	31.06
33	30.88
34	30.68
35	30.49
36	30.31
37	30.14
38	29.94
39	29.79
40	29.61
41	29.44
42	29.26
43	29.13
44	28.93
45	28.78
46	28.62
47	28.49
48	28.34
49	28.23
50	28.03
51	27.89
52	27.77
53	27.65
54	27.51
55	27.44
56	27.30
57	27.20
58	27.07
59	26.96
60	26.90
61	26.80
62	26.67
63	26.56
64	26.48
65	26.37
66	26.32
67	26.25
68	26.17
69	26.03
70	25.99
71	25.87
72	25.87
73	25.80
74	25.74
75	25.58
76	25.60
77	25.44
78	25.54
79	25.43
80	25.40
81	25.25
82	25.26
83	25.14
84	25.24
85	25.10
86	25.09
87	24.95
88	24.99
89	24.88
90	24.98
91	24.87
92	24.86
93	24.71
94	24.76
95	24.70
96	24.79
97	24.75
98	24.71
99	24.60
100	24.58
101	24.50
102	24.50
103	24.46
104	24.33
105	24.26
106	24.23
107	24.16
108	24.05
109	24.02
110	23.89
111	23.82
112	23.66
113	23.63
114	23.45
115	23.46
116	23.26
117	23.17
118	23.03
119	22.95
120	22.88
121	22.88
122	22.79
123	22.87
124	22.87
125	23.02
126	23.03
127	23.19
128	23.20
129	23.38
130	23.43
131	23.59
132	23.63
133	23.79
134	23.82
135	23.97
136	24.00
137	24.11
138	24.13
139	24.25
140	24.24
141	24.32
142	24.31
143	24.36
144	24.36
145	24.47
146	24.45
147	24.54
148	24.57
149	24.69
150	24.75
151	23.34
end
