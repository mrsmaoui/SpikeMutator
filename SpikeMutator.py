#!/usr/bin/env python

############# MODIFY PROGRAM PATHS ACCORDINGLY #############
SPIKE_PATH="/home/ubuntu/SpikeMutator"
SCWRL4_PATH="/home/ubuntu/SpikeMutator/extra_tools/SCWRL4.0"
GROMACE_BIN_PATH="/usr/local/gromacs/bin/gmx"
PDB2PQR_PATH="/home/ubuntu/SpikeMutator/extra_tools/pdb2pqr-1.7"
AQUASOL_PATH="/home/ubuntu/SpikeMutator/extra_tools/AquaSol_Complexes/bin"
############# YOU DONT NEED TO MODIFY PAST THIS POINT #############

from subprocess import call
import string
import sys 
import AquasolCPLX_routine
import time
import os

def main(argv):
	error = "ERROR: Please consult the README.txt FILE in your SpikeMutator directory to properly use this program"

	#verifying input variables are correct
	for i in range(len(argv)):
		if( i % 2 == 1):	#check odd values
			#exit if command doesn't exist
			if not (argv[i] == '-p' or argv[i] == '-m' or argv[i] == '-c' or argv[i] == '-a'):
				print error
				sys.exit(0)

	#get input parameters
	filename = 'spike'
	position = ''
	mutation = ''
	sequence = ''
	offset = 0
	conformation = ''
	action = 'mutation'

	fasta = open("spike.fasta", "r")
	sequence = fasta.readline().strip()

	try:
		for i in range(len(argv)):
			if(argv[i] == '-f'): filename = argv[i+1]
			elif (argv[i] == '-p'): position = argv[i+1]
			elif (argv[i] == '-m'): mutation = argv[i+1]
			elif (argv[i] == '-c'): conformation = argv[i+1]
			elif (argv[i] == '-a'): action = argv[i+1]

	except Exception:
		print error
		sys.exit(0)

	#check that input file is of pdb format

	if(sequence == ""):
		print error
		print "fasta file did not load properly"
		sys.exit(1)

	if(mutation != "" and mutation != "all" and ',' not in mutation):
		mutation = mutation.upper()
		#check that the mutation is a valid character
		mutationList = list("ARNDCQEGHILKMFPSTWYV")
		if(mutation not in mutationList):
			print "ERROR: Incorrect mutation choice"
			sys.exit(0)	

	if "down" not in conformation and "up" not in conformation:
		print error
		print "Kindly specify a conformation to use for the spike by using the -c flag: -c down, or -c up"
		sys.exit(1)

	if "energy" not in action and "mutation" not in action:
		print error
		print "You can only use the value 'energy' with the -a flag"
		sys.exit(1)

	#Determine which USAGE the user intends to perform
	if( (position != '') and (position != 'all') and ',' not in position and int(position) > 0 and len(mutation) == 1):
		params = {"filename":filename, "position":position, "mutation":mutation, 
			"sequence":sequence, "offset":offset, "conformation":conformation, "action":action}  	
                return mutate(params)

	else:
		print error
		sys.exit(0)


def calculateE(usage_path, filename):
	#Calculate Coulomb and LJ energies
	call("cd %s; %s/mdp/zeroMD.sh %s %s %s > %s/trash/result_%s.tmp 2>&1" % (usage_path, SPIKE_PATH, GROMACE_BIN_PATH, filename, SPIKE_PATH, usage_path, filename), shell=True)

	call("tail -10 %s/trash/result_%s.tmp > %s/trash/resultshort_%s.tmp" % (usage_path, filename, usage_path, filename), shell=True)

	fileGrid = open('%s/trash/resultshort_%s.tmp' % (usage_path, filename), 'r')
	fileLines = fileGrid.readlines()

	LJ = 0.0
	Coul = 0.0
	for one_line in fileLines:
		ones = string.split(one_line)
		if(len(ones) > 1 and ones[0] == 'LJ-14'):
			LJ = float(ones[1])
		if(len(ones) > 1 and ones[0] == 'Coulomb'):
			Coul = float(ones[2])

	#Calculate solvation energy
	temp = 300
	energy = AquasolCPLX_routine.applyAquaSol(AQUASOL_PATH, usage_path, PDB2PQR_PATH, SPIKE_PATH, temp, filename+'.pdb', filename)
	(Fw, Fv, Nw) = energy.split( )

	call("rm -f %s/trash/resultshort_%s.tmp %s/trash/result_%s.tmp" % (usage_path, filename, usage_path, filename), shell=True)

	return(Coul*0.238845897, LJ*0.238845897, Fw, Fv, Nw)


def mutateSequence(x, m, sequence, offset):
	originalSeq = list(sequence)
	originalSeq[ (x-offset) -1] = m
	return "".join(originalSeq)	
	

def mutate(params):

	filename = params['filename']
	mutation = params['mutation']
	position = int(params['position'])
	sequence = params['sequence']
	offset = params['offset']
	action = params['action']
	conformation = params['conformation']
	
        # Check if position is available in the 6VXX structure
	model_pos = mapRealToModel[position]
	if model_pos == 0:
		print("The selected position (%d) cannot be tested because it does not exist in the 6VXX protein structure" % model_pos)
		exit()

	#Perform Mutation on original sequence
	mutant_seq = mutateSequence(model_pos, mutation, sequence, offset)
	
	print(simulate(mutation, mutant_seq, position, filename, action))

	os.chdir(SPIKE_PATH)
	


def simulate(mutation, sequence, position, filename, action):
	os.chdir(SPIKE_PATH)	
	timestr = time.strftime("%Y%m%d-%H%M%S")
	usage_path = SPIKE_PATH + '/tmp/' + timestr

	#Create a workspace directory
	call('mkdir %s' % usage_path, shell=True)
	call('mkdir %s/trash' % usage_path, shell=True)

	#Generate .fasta file 
	call('echo %s > %s/%s_%d-%s.fasta' % (sequence, usage_path, filename,position,mutation), shell=True)
	
	#Perform mutation and generate structure file
	
	call("%s/Scwrl4 -i %s/spike_A.pdb -o %s/scwrl_spike_A_%d-%s.pdb -s %s/%s_%d-%s.fasta >/dev/null" % (SCWRL4_PATH, SPIKE_PATH, usage_path, position, mutation, usage_path, filename, position, mutation), shell=True)

	call("%s/Scwrl4 -i %s/spike_B.pdb -o %s/scwrl_spike_B_%d-%s.pdb -s %s/%s_%d-%s.fasta >/dev/null" % (SCWRL4_PATH, SPIKE_PATH, usage_path, position, mutation, usage_path, filename, position, mutation), shell=True)

	call("%s/Scwrl4 -i %s/spike_C.pdb -o %s/scwrl_spike_C_%d-%s.pdb -s %s/%s_%d-%s.fasta >/dev/null" % (SCWRL4_PATH, SPIKE_PATH, usage_path, position, mutation, usage_path, filename, position, mutation), shell=True)

	call("cat %s/scwrl_spike_A_%d-%s.pdb %s/scwrl_spike_B_%d-%s.pdb %s/scwrl_spike_C_%d-%s.pdb > %s/spike_%d-%s.pdb" % (usage_path, position, mutation, usage_path, position, mutation, usage_path, position, mutation, usage_path, position, mutation), shell=True)


	# Copy structure file into runs/ directory
	call("cp %s/spike_%d-%s.pdb %s/" % (usage_path, position, mutation, (SPIKE_PATH + '/structs/')), shell=True)
	
	resp = ""
	if action == 'energy':
                print "Position,Mutation,Coul,LJ,Solvation,Sequence"

		# Calculate Thermodynamics values
		(Coul, LJ, Fw, Fv, Nw) = calculateE(usage_path, "spike_%d-%s" % (position, mutation) )
		solvation = float(Fw) - float(Fv) - (float(Nw) * -0.0352669760297406)
		resp = str(position) + "," + mutation + "," + str(Coul) + "," + str(LJ) + "," + str(solvation) + "," + sequence

		#remove directory
		call('rm -rf %s/tmp/%s' % (SPIKE_PATH, timestr), shell=True)
	
		return resp

	return "Mutated Structure created at %s/spike_%d-%s.pdb" % ( (SPIKE_PATH + '/structs'), position, mutation)






mapRealToModel = dict()
mapRealToModel[1	] = 0
mapRealToModel[2	] = 0
mapRealToModel[3	] = 0
mapRealToModel[4	] = 0
mapRealToModel[5	] = 0
mapRealToModel[6	] = 0
mapRealToModel[7	] = 0
mapRealToModel[8	] = 0
mapRealToModel[9	] = 0
mapRealToModel[10	] = 0
mapRealToModel[11	] = 0
mapRealToModel[12	] = 0
mapRealToModel[13	] = 0
mapRealToModel[14	] = 0
mapRealToModel[15	] = 0
mapRealToModel[16	] = 0
mapRealToModel[17	] = 0
mapRealToModel[18	] = 0
mapRealToModel[19	] = 0
mapRealToModel[20	] = 0
mapRealToModel[21	] = 0
mapRealToModel[22	] = 0
mapRealToModel[23	] = 0
mapRealToModel[24	] = 0
mapRealToModel[25	] = 0
mapRealToModel[26	] = 0
mapRealToModel[27	] = 	1
mapRealToModel[28	] = 	2
mapRealToModel[29	] = 	3
mapRealToModel[30	] = 	4
mapRealToModel[31	] = 	5
mapRealToModel[32	] = 	6
mapRealToModel[33	] = 	7
mapRealToModel[34	] = 	8
mapRealToModel[35	] = 	9
mapRealToModel[36	] = 	10
mapRealToModel[37	] = 	11
mapRealToModel[38	] = 	12
mapRealToModel[39	] = 	13
mapRealToModel[40	] = 	14
mapRealToModel[41	] = 	15
mapRealToModel[42	] = 	16
mapRealToModel[43	] = 	17
mapRealToModel[44	] = 	18
mapRealToModel[45	] = 	19
mapRealToModel[46	] = 	20
mapRealToModel[47	] = 	21
mapRealToModel[48	] = 	22
mapRealToModel[49	] = 	23
mapRealToModel[50	] = 	24
mapRealToModel[51	] = 	25
mapRealToModel[52	] = 	26
mapRealToModel[53	] = 	27
mapRealToModel[54	] = 	28
mapRealToModel[55	] = 	29
mapRealToModel[56	] = 	30
mapRealToModel[57	] = 	31
mapRealToModel[58	] = 	32
mapRealToModel[59	] = 	33
mapRealToModel[60	] = 	34
mapRealToModel[61	] = 	35
mapRealToModel[62	] = 	36
mapRealToModel[63	] = 	37
mapRealToModel[64	] = 	38
mapRealToModel[65	] = 	39
mapRealToModel[66	] = 	40
mapRealToModel[67	] = 	41
mapRealToModel[68	] = 	42
mapRealToModel[69	] = 	43
mapRealToModel[70	] = 0
mapRealToModel[71	] = 0
mapRealToModel[72	] = 0
mapRealToModel[73	] = 0
mapRealToModel[74	] = 0
mapRealToModel[75	] = 0
mapRealToModel[76	] = 0
mapRealToModel[77	] = 0
mapRealToModel[78	] = 0
mapRealToModel[79	] = 0
mapRealToModel[80	] = 	44
mapRealToModel[81	] = 	45
mapRealToModel[82	] = 	46
mapRealToModel[83	] = 	47
mapRealToModel[84	] = 	48
mapRealToModel[85	] = 	49
mapRealToModel[86	] = 	50
mapRealToModel[87	] = 	51
mapRealToModel[88	] = 	52
mapRealToModel[89	] = 	53
mapRealToModel[90	] = 	54
mapRealToModel[91	] = 	55
mapRealToModel[92	] = 	56
mapRealToModel[93	] = 	57
mapRealToModel[94	] = 	58
mapRealToModel[95	] = 	59
mapRealToModel[96	] = 	60
mapRealToModel[97	] = 	61
mapRealToModel[98	] = 	62
mapRealToModel[99	] = 	63
mapRealToModel[100	] = 	64
mapRealToModel[101	] = 	65
mapRealToModel[102	] = 	66
mapRealToModel[103	] = 	67
mapRealToModel[104	] = 	68
mapRealToModel[105	] = 	69
mapRealToModel[106	] = 	70
mapRealToModel[107	] = 	71
mapRealToModel[108	] = 	72
mapRealToModel[109	] = 	73
mapRealToModel[110	] = 	74
mapRealToModel[111	] = 	75
mapRealToModel[112	] = 	76
mapRealToModel[113	] = 	77
mapRealToModel[114	] = 	78
mapRealToModel[115	] = 	79
mapRealToModel[116	] = 	80
mapRealToModel[117	] = 	81
mapRealToModel[118	] = 	82
mapRealToModel[119	] = 	83
mapRealToModel[120	] = 	84
mapRealToModel[121	] = 	85
mapRealToModel[122	] = 	86
mapRealToModel[123	] = 	87
mapRealToModel[124	] = 	88
mapRealToModel[125	] = 	89
mapRealToModel[126	] = 	90
mapRealToModel[127	] = 	91
mapRealToModel[128	] = 	92
mapRealToModel[129	] = 	93
mapRealToModel[130	] = 	94
mapRealToModel[131	] = 	95
mapRealToModel[132	] = 	96
mapRealToModel[133	] = 	97
mapRealToModel[134	] = 	98
mapRealToModel[135	] = 	99
mapRealToModel[136	] = 	100
mapRealToModel[137	] = 	101
mapRealToModel[138	] = 	102
mapRealToModel[139	] = 	103
mapRealToModel[140	] = 	104
mapRealToModel[141	] = 	105
mapRealToModel[142	] = 	106
mapRealToModel[143	] = 	107
mapRealToModel[144	] = 0
mapRealToModel[145	] = 0
mapRealToModel[146	] = 0
mapRealToModel[147	] = 0
mapRealToModel[148	] = 0
mapRealToModel[149	] = 0
mapRealToModel[150	] = 0
mapRealToModel[151	] = 0
mapRealToModel[152	] = 0
mapRealToModel[153	] = 0
mapRealToModel[154	] = 0
mapRealToModel[155	] = 0
mapRealToModel[156	] = 0
mapRealToModel[157	] = 0
mapRealToModel[158	] = 0
mapRealToModel[159	] = 0
mapRealToModel[160	] = 0
mapRealToModel[161	] = 0
mapRealToModel[162	] = 0
mapRealToModel[163	] = 0
mapRealToModel[164	] = 0
mapRealToModel[165	] = 	108
mapRealToModel[166	] = 	109
mapRealToModel[167	] = 	110
mapRealToModel[168	] = 	111
mapRealToModel[169	] = 	112
mapRealToModel[170	] = 	113
mapRealToModel[171	] = 	114
mapRealToModel[172	] = 	115
mapRealToModel[173	] = 0
mapRealToModel[174	] = 0
mapRealToModel[175	] = 0
mapRealToModel[176	] = 0
mapRealToModel[177	] = 0
mapRealToModel[178	] = 0
mapRealToModel[179	] = 0
mapRealToModel[180	] = 0
mapRealToModel[181	] = 0
mapRealToModel[182	] = 0
mapRealToModel[183	] = 0
mapRealToModel[184	] = 0
mapRealToModel[185	] = 0
mapRealToModel[186	] = 	116
mapRealToModel[187	] = 	117
mapRealToModel[188	] = 	118
mapRealToModel[189	] = 	119
mapRealToModel[190	] = 	120
mapRealToModel[191	] = 	121
mapRealToModel[192	] = 	122
mapRealToModel[193	] = 	123
mapRealToModel[194	] = 	124
mapRealToModel[195	] = 	125
mapRealToModel[196	] = 	126
mapRealToModel[197	] = 	127
mapRealToModel[198	] = 	128
mapRealToModel[199	] = 	129
mapRealToModel[200	] = 	130
mapRealToModel[201	] = 	131
mapRealToModel[202	] = 	132
mapRealToModel[203	] = 	133
mapRealToModel[204	] = 	134
mapRealToModel[205	] = 	135
mapRealToModel[206	] = 	136
mapRealToModel[207	] = 	137
mapRealToModel[208	] = 	138
mapRealToModel[209	] = 	139
mapRealToModel[210	] = 	140
mapRealToModel[211	] = 	141
mapRealToModel[212	] = 	142
mapRealToModel[213	] = 	143
mapRealToModel[214	] = 	144
mapRealToModel[215	] = 	145
mapRealToModel[216	] = 	146
mapRealToModel[217	] = 	147
mapRealToModel[218	] = 	148
mapRealToModel[219	] = 	149
mapRealToModel[220	] = 	150
mapRealToModel[221	] = 	151
mapRealToModel[222	] = 	152
mapRealToModel[223	] = 	153
mapRealToModel[224	] = 	154
mapRealToModel[225	] = 	155
mapRealToModel[226	] = 	156
mapRealToModel[227	] = 	157
mapRealToModel[228	] = 	158
mapRealToModel[229	] = 	159
mapRealToModel[230	] = 	160
mapRealToModel[231	] = 	161
mapRealToModel[232	] = 	162
mapRealToModel[233	] = 	163
mapRealToModel[234	] = 	164
mapRealToModel[235	] = 	165
mapRealToModel[236	] = 	166
mapRealToModel[237	] = 	167
mapRealToModel[238	] = 	168
mapRealToModel[239	] = 	169
mapRealToModel[240	] = 	170
mapRealToModel[241	] = 	171
mapRealToModel[242	] = 	172
mapRealToModel[243	] = 	173
mapRealToModel[244	] = 	174
mapRealToModel[245	] = 	175
mapRealToModel[246	] = 0
mapRealToModel[247	] = 0
mapRealToModel[248	] = 0
mapRealToModel[249	] = 0
mapRealToModel[250	] = 0
mapRealToModel[251	] = 0
mapRealToModel[252	] = 0
mapRealToModel[253	] = 0
mapRealToModel[254	] = 0
mapRealToModel[255	] = 0
mapRealToModel[256	] = 0
mapRealToModel[257	] = 0
mapRealToModel[258	] = 0
mapRealToModel[259	] = 0
mapRealToModel[260	] = 0
mapRealToModel[261	] = 0
mapRealToModel[262	] = 0
mapRealToModel[263	] = 	176
mapRealToModel[264	] = 	177
mapRealToModel[265	] = 	178
mapRealToModel[266	] = 	179
mapRealToModel[267	] = 	180
mapRealToModel[268	] = 	181
mapRealToModel[269	] = 	182
mapRealToModel[270	] = 	183
mapRealToModel[271	] = 	184
mapRealToModel[272	] = 	185
mapRealToModel[273	] = 	186
mapRealToModel[274	] = 	187
mapRealToModel[275	] = 	188
mapRealToModel[276	] = 	189
mapRealToModel[277	] = 	190
mapRealToModel[278	] = 	191
mapRealToModel[279	] = 	192
mapRealToModel[280	] = 	193
mapRealToModel[281	] = 	194
mapRealToModel[282	] = 	195
mapRealToModel[283	] = 	196
mapRealToModel[284	] = 	197
mapRealToModel[285	] = 	198
mapRealToModel[286	] = 	199
mapRealToModel[287	] = 	200
mapRealToModel[288	] = 	201
mapRealToModel[289	] = 	202
mapRealToModel[290	] = 	203
mapRealToModel[291	] = 	204
mapRealToModel[292	] = 	205
mapRealToModel[293	] = 	206
mapRealToModel[294	] = 	207
mapRealToModel[295	] = 	208
mapRealToModel[296	] = 	209
mapRealToModel[297	] = 	210
mapRealToModel[298	] = 	211
mapRealToModel[299	] = 	212
mapRealToModel[300	] = 	213
mapRealToModel[301	] = 	214
mapRealToModel[302	] = 	215
mapRealToModel[303	] = 	216
mapRealToModel[304	] = 	217
mapRealToModel[305	] = 	218
mapRealToModel[306	] = 	219
mapRealToModel[307	] = 	220
mapRealToModel[308	] = 	221
mapRealToModel[309	] = 	222
mapRealToModel[310	] = 	223
mapRealToModel[311	] = 	224
mapRealToModel[312	] = 	225
mapRealToModel[313	] = 	226
mapRealToModel[314	] = 	227
mapRealToModel[315	] = 	228
mapRealToModel[316	] = 	229
mapRealToModel[317	] = 	230
mapRealToModel[318	] = 	231
mapRealToModel[319	] = 	232
mapRealToModel[320	] = 	233
mapRealToModel[321	] = 	234
mapRealToModel[322	] = 	235
mapRealToModel[323	] = 	236
mapRealToModel[324	] = 	237
mapRealToModel[325	] = 	238
mapRealToModel[326	] = 	239
mapRealToModel[327	] = 	240
mapRealToModel[328	] = 	241
mapRealToModel[329	] = 	242
mapRealToModel[330	] = 	243
mapRealToModel[331	] = 	244
mapRealToModel[332	] = 	245
mapRealToModel[333	] = 	246
mapRealToModel[334	] = 	247
mapRealToModel[335	] = 	248
mapRealToModel[336	] = 	249
mapRealToModel[337	] = 	250
mapRealToModel[338	] = 	251
mapRealToModel[339	] = 	252
mapRealToModel[340	] = 	253
mapRealToModel[341	] = 	254
mapRealToModel[342	] = 	255
mapRealToModel[343	] = 	256
mapRealToModel[344	] = 	257
mapRealToModel[345	] = 	258
mapRealToModel[346	] = 	259
mapRealToModel[347	] = 	260
mapRealToModel[348	] = 	261
mapRealToModel[349	] = 	262
mapRealToModel[350	] = 	263
mapRealToModel[351	] = 	264
mapRealToModel[352	] = 	265
mapRealToModel[353	] = 	266
mapRealToModel[354	] = 	267
mapRealToModel[355	] = 	268
mapRealToModel[356	] = 	269
mapRealToModel[357	] = 	270
mapRealToModel[358	] = 	271
mapRealToModel[359	] = 	272
mapRealToModel[360	] = 	273
mapRealToModel[361	] = 	274
mapRealToModel[362	] = 	275
mapRealToModel[363	] = 	276
mapRealToModel[364	] = 	277
mapRealToModel[365	] = 	278
mapRealToModel[366	] = 	279
mapRealToModel[367	] = 	280
mapRealToModel[368	] = 	281
mapRealToModel[369	] = 	282
mapRealToModel[370	] = 	283
mapRealToModel[371	] = 	284
mapRealToModel[372	] = 	285
mapRealToModel[373	] = 	286
mapRealToModel[374	] = 	287
mapRealToModel[375	] = 	288
mapRealToModel[376	] = 	289
mapRealToModel[377	] = 	290
mapRealToModel[378	] = 	291
mapRealToModel[379	] = 	292
mapRealToModel[380	] = 	293
mapRealToModel[381	] = 	294
mapRealToModel[382	] = 	295
mapRealToModel[383	] = 	296
mapRealToModel[384	] = 	297
mapRealToModel[385	] = 	298
mapRealToModel[386	] = 	299
mapRealToModel[387	] = 	300
mapRealToModel[388	] = 	301
mapRealToModel[389	] = 	302
mapRealToModel[390	] = 	303
mapRealToModel[391	] = 	304
mapRealToModel[392	] = 	305
mapRealToModel[393	] = 	306
mapRealToModel[394	] = 	307
mapRealToModel[395	] = 	308
mapRealToModel[396	] = 	309
mapRealToModel[397	] = 	310
mapRealToModel[398	] = 	311
mapRealToModel[399	] = 	312
mapRealToModel[400	] = 	313
mapRealToModel[401	] = 	314
mapRealToModel[402	] = 	315
mapRealToModel[403	] = 	316
mapRealToModel[404	] = 	317
mapRealToModel[405	] = 	318
mapRealToModel[406	] = 	319
mapRealToModel[407	] = 	320
mapRealToModel[408	] = 	321
mapRealToModel[409	] = 	322
mapRealToModel[410	] = 	323
mapRealToModel[411	] = 	324
mapRealToModel[412	] = 	325
mapRealToModel[413	] = 	326
mapRealToModel[414	] = 	327
mapRealToModel[415	] = 	328
mapRealToModel[416	] = 	329
mapRealToModel[417	] = 	330
mapRealToModel[418	] = 	331
mapRealToModel[419	] = 	332
mapRealToModel[420	] = 	333
mapRealToModel[421	] = 	334
mapRealToModel[422	] = 	335
mapRealToModel[423	] = 	336
mapRealToModel[424	] = 	337
mapRealToModel[425	] = 	338
mapRealToModel[426	] = 	339
mapRealToModel[427	] = 	340
mapRealToModel[428	] = 	341
mapRealToModel[429	] = 	342
mapRealToModel[430	] = 	343
mapRealToModel[431	] = 	344
mapRealToModel[432	] = 	345
mapRealToModel[433	] = 	346
mapRealToModel[434	] = 	347
mapRealToModel[435	] = 	348
mapRealToModel[436	] = 	349
mapRealToModel[437	] = 	350
mapRealToModel[438	] = 	351
mapRealToModel[439	] = 	352
mapRealToModel[440	] = 	353
mapRealToModel[441	] = 	354
mapRealToModel[442	] = 	355
mapRealToModel[443	] = 	356
mapRealToModel[444	] = 	357
mapRealToModel[445	] = 0
mapRealToModel[446	] = 0
mapRealToModel[447	] = 	358
mapRealToModel[448	] = 	359
mapRealToModel[449	] = 	360
mapRealToModel[450	] = 	361
mapRealToModel[451	] = 	362
mapRealToModel[452	] = 	363
mapRealToModel[453	] = 	364
mapRealToModel[454	] = 	365
mapRealToModel[455	] = 0
mapRealToModel[456	] = 0
mapRealToModel[457	] = 0
mapRealToModel[458	] = 0
mapRealToModel[459	] = 0
mapRealToModel[460	] = 0
mapRealToModel[461	] = 0
mapRealToModel[462	] = 	366
mapRealToModel[463	] = 	367
mapRealToModel[464	] = 	368
mapRealToModel[465	] = 	369
mapRealToModel[466	] = 	370
mapRealToModel[467	] = 	371
mapRealToModel[468	] = 	372
mapRealToModel[469	] = 0
mapRealToModel[470	] = 0
mapRealToModel[471	] = 0
mapRealToModel[472	] = 0
mapRealToModel[473	] = 0
mapRealToModel[474	] = 0
mapRealToModel[475	] = 0
mapRealToModel[476	] = 0
mapRealToModel[477	] = 0
mapRealToModel[478	] = 0
mapRealToModel[479	] = 0
mapRealToModel[480	] = 0
mapRealToModel[481	] = 0
mapRealToModel[482	] = 0
mapRealToModel[483	] = 0
mapRealToModel[484	] = 0
mapRealToModel[485	] = 0
mapRealToModel[486	] = 0
mapRealToModel[487	] = 0
mapRealToModel[488	] = 0
mapRealToModel[489	] = 	373
mapRealToModel[490	] = 	374
mapRealToModel[491	] = 	375
mapRealToModel[492	] = 	376
mapRealToModel[493	] = 	377
mapRealToModel[494	] = 	378
mapRealToModel[495	] = 	379
mapRealToModel[496	] = 	380
mapRealToModel[497	] = 	381
mapRealToModel[498	] = 	382
mapRealToModel[499	] = 	383
mapRealToModel[500	] = 	384
mapRealToModel[501	] = 	385
mapRealToModel[502	] = 0
mapRealToModel[503	] = 	386
mapRealToModel[504	] = 	387
mapRealToModel[505	] = 	388
mapRealToModel[506	] = 	389
mapRealToModel[507	] = 	390
mapRealToModel[508	] = 	391
mapRealToModel[509	] = 	392
mapRealToModel[510	] = 	393
mapRealToModel[511	] = 	394
mapRealToModel[512	] = 	395
mapRealToModel[513	] = 	396
mapRealToModel[514	] = 	397
mapRealToModel[515	] = 	398
mapRealToModel[516	] = 	399
mapRealToModel[517	] = 	400
mapRealToModel[518	] = 	401
mapRealToModel[519	] = 	402
mapRealToModel[520	] = 	403
mapRealToModel[521	] = 	404
mapRealToModel[522	] = 	405
mapRealToModel[523	] = 	406
mapRealToModel[524	] = 	407
mapRealToModel[525	] = 	408
mapRealToModel[526	] = 	409
mapRealToModel[527	] = 	410
mapRealToModel[528	] = 	411
mapRealToModel[529	] = 	412
mapRealToModel[530	] = 	413
mapRealToModel[531	] = 	414
mapRealToModel[532	] = 	415
mapRealToModel[533	] = 	416
mapRealToModel[534	] = 	417
mapRealToModel[535	] = 	418
mapRealToModel[536	] = 	419
mapRealToModel[537	] = 	420
mapRealToModel[538	] = 	421
mapRealToModel[539	] = 	422
mapRealToModel[540	] = 	423
mapRealToModel[541	] = 	424
mapRealToModel[542	] = 	425
mapRealToModel[543	] = 	426
mapRealToModel[544	] = 	427
mapRealToModel[545	] = 	428
mapRealToModel[546	] = 	429
mapRealToModel[547	] = 	430
mapRealToModel[548	] = 	431
mapRealToModel[549	] = 	432
mapRealToModel[550	] = 	433
mapRealToModel[551	] = 	434
mapRealToModel[552	] = 	435
mapRealToModel[553	] = 	436
mapRealToModel[554	] = 	437
mapRealToModel[555	] = 	438
mapRealToModel[556	] = 	439
mapRealToModel[557	] = 	440
mapRealToModel[558	] = 	441
mapRealToModel[559	] = 	442
mapRealToModel[560	] = 	443
mapRealToModel[561	] = 	444
mapRealToModel[562	] = 	445
mapRealToModel[563	] = 	446
mapRealToModel[564	] = 	447
mapRealToModel[565	] = 	448
mapRealToModel[566	] = 	449
mapRealToModel[567	] = 	450
mapRealToModel[568	] = 	451
mapRealToModel[569	] = 	452
mapRealToModel[570	] = 	453
mapRealToModel[571	] = 	454
mapRealToModel[572	] = 	455
mapRealToModel[573	] = 	456
mapRealToModel[574	] = 	457
mapRealToModel[575	] = 	458
mapRealToModel[576	] = 	459
mapRealToModel[577	] = 	460
mapRealToModel[578	] = 	461
mapRealToModel[579	] = 	462
mapRealToModel[580	] = 	463
mapRealToModel[581	] = 	464
mapRealToModel[582	] = 	465
mapRealToModel[583	] = 	466
mapRealToModel[584	] = 	467
mapRealToModel[585	] = 	468
mapRealToModel[586	] = 	469
mapRealToModel[587	] = 	470
mapRealToModel[588	] = 	471
mapRealToModel[589	] = 	472
mapRealToModel[590	] = 	473
mapRealToModel[591	] = 	474
mapRealToModel[592	] = 	475
mapRealToModel[593	] = 	476
mapRealToModel[594	] = 	477
mapRealToModel[595	] = 	478
mapRealToModel[596	] = 	479
mapRealToModel[597	] = 	480
mapRealToModel[598	] = 	481
mapRealToModel[599	] = 	482
mapRealToModel[600	] = 	483
mapRealToModel[601	] = 	484
mapRealToModel[602	] = 	485
mapRealToModel[603	] = 	486
mapRealToModel[604	] = 	487
mapRealToModel[605	] = 	488
mapRealToModel[606	] = 	489
mapRealToModel[607	] = 	490
mapRealToModel[608	] = 	491
mapRealToModel[609	] = 	492
mapRealToModel[610	] = 	493
mapRealToModel[611	] = 	494
mapRealToModel[612	] = 	495
mapRealToModel[613	] = 	496
mapRealToModel[614	] = 	497
mapRealToModel[615	] = 	498
mapRealToModel[616	] = 	499
mapRealToModel[617	] = 	500
mapRealToModel[618	] = 	501
mapRealToModel[619	] = 	502
mapRealToModel[620	] = 	503
mapRealToModel[621	] = 0
mapRealToModel[622	] = 0
mapRealToModel[623	] = 0
mapRealToModel[624	] = 0
mapRealToModel[625	] = 0
mapRealToModel[626	] = 0
mapRealToModel[627	] = 0
mapRealToModel[628	] = 0
mapRealToModel[629	] = 0
mapRealToModel[630	] = 0
mapRealToModel[631	] = 0
mapRealToModel[632	] = 0
mapRealToModel[633	] = 0
mapRealToModel[634	] = 0
mapRealToModel[635	] = 0
mapRealToModel[636	] = 0
mapRealToModel[637	] = 0
mapRealToModel[638	] = 0
mapRealToModel[639	] = 0
mapRealToModel[640	] = 0
mapRealToModel[641	] = 	504
mapRealToModel[642	] = 	505
mapRealToModel[643	] = 	506
mapRealToModel[644	] = 	507
mapRealToModel[645	] = 	508
mapRealToModel[646	] = 	509
mapRealToModel[647	] = 	510
mapRealToModel[648	] = 	511
mapRealToModel[649	] = 	512
mapRealToModel[650	] = 	513
mapRealToModel[651	] = 	514
mapRealToModel[652	] = 	515
mapRealToModel[653	] = 	516
mapRealToModel[654	] = 	517
mapRealToModel[655	] = 	518
mapRealToModel[656	] = 	519
mapRealToModel[657	] = 	520
mapRealToModel[658	] = 	521
mapRealToModel[659	] = 	522
mapRealToModel[660	] = 	523
mapRealToModel[661	] = 	524
mapRealToModel[662	] = 	525
mapRealToModel[663	] = 	526
mapRealToModel[664	] = 	527
mapRealToModel[665	] = 	528
mapRealToModel[666	] = 	529
mapRealToModel[667	] = 	530
mapRealToModel[668	] = 	531
mapRealToModel[669	] = 	532
mapRealToModel[670	] = 	533
mapRealToModel[671	] = 	534
mapRealToModel[672	] = 	535
mapRealToModel[673	] = 	536
mapRealToModel[674	] = 	537
mapRealToModel[675	] = 	538
mapRealToModel[676	] = 	539
mapRealToModel[677	] = 0
mapRealToModel[678	] = 0
mapRealToModel[679	] = 0
mapRealToModel[680	] = 0
mapRealToModel[681	] = 0
mapRealToModel[682	] = 0
mapRealToModel[683	] = 0
mapRealToModel[684	] = 0
mapRealToModel[685	] = 0
mapRealToModel[686	] = 0
mapRealToModel[687	] = 0
mapRealToModel[688	] = 0
mapRealToModel[689	] = 	540
mapRealToModel[690	] = 	541
mapRealToModel[691	] = 	542
mapRealToModel[692	] = 	543
mapRealToModel[693	] = 	544
mapRealToModel[694	] = 	545
mapRealToModel[695	] = 	546
mapRealToModel[696	] = 	547
mapRealToModel[697	] = 	548
mapRealToModel[698	] = 	549
mapRealToModel[699	] = 	550
mapRealToModel[700	] = 	551
mapRealToModel[701	] = 	552
mapRealToModel[702	] = 	553
mapRealToModel[703	] = 	554
mapRealToModel[704	] = 	555
mapRealToModel[705	] = 	556
mapRealToModel[706	] = 	557
mapRealToModel[707	] = 	558
mapRealToModel[708	] = 	559
mapRealToModel[709	] = 	560
mapRealToModel[710	] = 	561
mapRealToModel[711	] = 	562
mapRealToModel[712	] = 	563
mapRealToModel[713	] = 	564
mapRealToModel[714	] = 	565
mapRealToModel[715	] = 	566
mapRealToModel[716	] = 	567
mapRealToModel[717	] = 	568
mapRealToModel[718	] = 	569
mapRealToModel[719	] = 	570
mapRealToModel[720	] = 	571
mapRealToModel[721	] = 	572
mapRealToModel[722	] = 	573
mapRealToModel[723	] = 	574
mapRealToModel[724	] = 	575
mapRealToModel[725	] = 	576
mapRealToModel[726	] = 	577
mapRealToModel[727	] = 	578
mapRealToModel[728	] = 	579
mapRealToModel[729	] = 	580
mapRealToModel[730	] = 	581
mapRealToModel[731	] = 	582
mapRealToModel[732	] = 	583
mapRealToModel[733	] = 	584
mapRealToModel[734	] = 	585
mapRealToModel[735	] = 	586
mapRealToModel[736	] = 	587
mapRealToModel[737	] = 	588
mapRealToModel[738	] = 	589
mapRealToModel[739	] = 	590
mapRealToModel[740	] = 	591
mapRealToModel[741	] = 	592
mapRealToModel[742	] = 	593
mapRealToModel[743	] = 	594
mapRealToModel[744	] = 	595
mapRealToModel[745	] = 	596
mapRealToModel[746	] = 	597
mapRealToModel[747	] = 	598
mapRealToModel[748	] = 	599
mapRealToModel[749	] = 	600
mapRealToModel[750	] = 	601
mapRealToModel[751	] = 	602
mapRealToModel[752	] = 	603
mapRealToModel[753	] = 	604
mapRealToModel[754	] = 	605
mapRealToModel[755	] = 	606
mapRealToModel[756	] = 	607
mapRealToModel[757	] = 	608
mapRealToModel[758	] = 	609
mapRealToModel[759	] = 	610
mapRealToModel[760	] = 	611
mapRealToModel[761	] = 	612
mapRealToModel[762	] = 	613
mapRealToModel[763	] = 	614
mapRealToModel[764	] = 	615
mapRealToModel[765	] = 	616
mapRealToModel[766	] = 	617
mapRealToModel[767	] = 	618
mapRealToModel[768	] = 	619
mapRealToModel[769	] = 	620
mapRealToModel[770	] = 	621
mapRealToModel[771	] = 	622
mapRealToModel[772	] = 	623
mapRealToModel[773	] = 	624
mapRealToModel[774	] = 	625
mapRealToModel[775	] = 	626
mapRealToModel[776	] = 	627
mapRealToModel[777	] = 	628
mapRealToModel[778	] = 	629
mapRealToModel[779	] = 	630
mapRealToModel[780	] = 	631
mapRealToModel[781	] = 	632
mapRealToModel[782	] = 	633
mapRealToModel[783	] = 	634
mapRealToModel[784	] = 	635
mapRealToModel[785	] = 	636
mapRealToModel[786	] = 	637
mapRealToModel[787	] = 	638
mapRealToModel[788	] = 	639
mapRealToModel[789	] = 	640
mapRealToModel[790	] = 	641
mapRealToModel[791	] = 	642
mapRealToModel[792	] = 	643
mapRealToModel[793	] = 	644
mapRealToModel[794	] = 	645
mapRealToModel[795	] = 	646
mapRealToModel[796	] = 	647
mapRealToModel[797	] = 	648
mapRealToModel[798	] = 	649
mapRealToModel[799	] = 	650
mapRealToModel[800	] = 	651
mapRealToModel[801	] = 	652
mapRealToModel[802	] = 	653
mapRealToModel[803	] = 	654
mapRealToModel[804	] = 	655
mapRealToModel[805	] = 	656
mapRealToModel[806	] = 	657
mapRealToModel[807	] = 	658
mapRealToModel[808	] = 	659
mapRealToModel[809	] = 	660
mapRealToModel[810	] = 	661
mapRealToModel[811	] = 	662
mapRealToModel[812	] = 	663
mapRealToModel[813	] = 	664
mapRealToModel[814	] = 	665
mapRealToModel[815	] = 	666
mapRealToModel[816	] = 	667
mapRealToModel[817	] = 	668
mapRealToModel[818	] = 	669
mapRealToModel[819	] = 	670
mapRealToModel[820	] = 	671
mapRealToModel[821	] = 	672
mapRealToModel[822	] = 	673
mapRealToModel[823	] = 	674
mapRealToModel[824	] = 	675
mapRealToModel[825	] = 	676
mapRealToModel[826	] = 	677
mapRealToModel[827	] = 	678
mapRealToModel[828	] = 0
mapRealToModel[829	] = 0
mapRealToModel[830	] = 0
mapRealToModel[831	] = 0
mapRealToModel[832	] = 0
mapRealToModel[833	] = 0
mapRealToModel[834	] = 0
mapRealToModel[835	] = 0
mapRealToModel[836	] = 0
mapRealToModel[837	] = 0
mapRealToModel[838	] = 0
mapRealToModel[839	] = 0
mapRealToModel[840	] = 0
mapRealToModel[841	] = 0
mapRealToModel[842	] = 0
mapRealToModel[843	] = 0
mapRealToModel[844	] = 0
mapRealToModel[845	] = 0
mapRealToModel[846	] = 0
mapRealToModel[847	] = 0
mapRealToModel[848	] = 0
mapRealToModel[849	] = 0
mapRealToModel[850	] = 0
mapRealToModel[851	] = 0
mapRealToModel[852	] = 0
mapRealToModel[853	] = 0
mapRealToModel[854	] = 	679
mapRealToModel[855	] = 	680
mapRealToModel[856	] = 	681
mapRealToModel[857	] = 	682
mapRealToModel[858	] = 	683
mapRealToModel[859	] = 	684
mapRealToModel[860	] = 	685
mapRealToModel[861	] = 	686
mapRealToModel[862	] = 	687
mapRealToModel[863	] = 	688
mapRealToModel[864	] = 	689
mapRealToModel[865	] = 	690
mapRealToModel[866	] = 	691
mapRealToModel[867	] = 	692
mapRealToModel[868	] = 	693
mapRealToModel[869	] = 	694
mapRealToModel[870	] = 	695
mapRealToModel[871	] = 	696
mapRealToModel[872	] = 	697
mapRealToModel[873	] = 	698
mapRealToModel[874	] = 	699
mapRealToModel[875	] = 	700
mapRealToModel[876	] = 	701
mapRealToModel[877	] = 	702
mapRealToModel[878	] = 	703
mapRealToModel[879	] = 	704
mapRealToModel[880	] = 	705
mapRealToModel[881	] = 	706
mapRealToModel[882	] = 	707
mapRealToModel[883	] = 	708
mapRealToModel[884	] = 	709
mapRealToModel[885	] = 	710
mapRealToModel[886	] = 	711
mapRealToModel[887	] = 	712
mapRealToModel[888	] = 	713
mapRealToModel[889	] = 	714
mapRealToModel[890	] = 	715
mapRealToModel[891	] = 	716
mapRealToModel[892	] = 	717
mapRealToModel[893	] = 	718
mapRealToModel[894	] = 	719
mapRealToModel[895	] = 	720
mapRealToModel[896	] = 	721
mapRealToModel[897	] = 	722
mapRealToModel[898	] = 	723
mapRealToModel[899	] = 	724
mapRealToModel[900	] = 	725
mapRealToModel[901	] = 	726
mapRealToModel[902	] = 	727
mapRealToModel[903	] = 	728
mapRealToModel[904	] = 	729
mapRealToModel[905	] = 	730
mapRealToModel[906	] = 	731
mapRealToModel[907	] = 	732
mapRealToModel[908	] = 	733
mapRealToModel[909	] = 	734
mapRealToModel[910	] = 	735
mapRealToModel[911	] = 	736
mapRealToModel[912	] = 	737
mapRealToModel[913	] = 	738
mapRealToModel[914	] = 	739
mapRealToModel[915	] = 	740
mapRealToModel[916	] = 	741
mapRealToModel[917	] = 	742
mapRealToModel[918	] = 	743
mapRealToModel[919	] = 	744
mapRealToModel[920	] = 	745
mapRealToModel[921	] = 	746
mapRealToModel[922	] = 	747
mapRealToModel[923	] = 	748
mapRealToModel[924	] = 	749
mapRealToModel[925	] = 	750
mapRealToModel[926	] = 	751
mapRealToModel[927	] = 	752
mapRealToModel[928	] = 	753
mapRealToModel[929	] = 	754
mapRealToModel[930	] = 	755
mapRealToModel[931	] = 	756
mapRealToModel[932	] = 	757
mapRealToModel[933	] = 	758
mapRealToModel[934	] = 	759
mapRealToModel[935	] = 	760
mapRealToModel[936	] = 	761
mapRealToModel[937	] = 	762
mapRealToModel[938	] = 	763
mapRealToModel[939	] = 	764
mapRealToModel[940	] = 	765
mapRealToModel[941	] = 	766
mapRealToModel[942	] = 	767
mapRealToModel[943	] = 	768
mapRealToModel[944	] = 	769
mapRealToModel[945	] = 	770
mapRealToModel[946	] = 	771
mapRealToModel[947	] = 	772
mapRealToModel[948	] = 	773
mapRealToModel[949	] = 	774
mapRealToModel[950	] = 	775
mapRealToModel[951	] = 	776
mapRealToModel[952	] = 	777
mapRealToModel[953	] = 	778
mapRealToModel[954	] = 	779
mapRealToModel[955	] = 	780
mapRealToModel[956	] = 	781
mapRealToModel[957	] = 	782
mapRealToModel[958	] = 	783
mapRealToModel[959	] = 	784
mapRealToModel[960	] = 	785
mapRealToModel[961	] = 	786
mapRealToModel[962	] = 	787
mapRealToModel[963	] = 	788
mapRealToModel[964	] = 	789
mapRealToModel[965	] = 	790
mapRealToModel[966	] = 	791
mapRealToModel[967	] = 	792
mapRealToModel[968	] = 	793
mapRealToModel[969	] = 	794
mapRealToModel[970	] = 	795
mapRealToModel[971	] = 	796
mapRealToModel[972	] = 	797
mapRealToModel[973	] = 	798
mapRealToModel[974	] = 	799
mapRealToModel[975	] = 	800
mapRealToModel[976	] = 	801
mapRealToModel[977	] = 	802
mapRealToModel[978	] = 	803
mapRealToModel[979	] = 	804
mapRealToModel[980	] = 	805
mapRealToModel[981	] = 	806
mapRealToModel[982	] = 	807
mapRealToModel[983	] = 	808
mapRealToModel[984	] = 	809
mapRealToModel[985	] = 	810
mapRealToModel[986	] = 	811
mapRealToModel[987	] = 	812
mapRealToModel[988	] = 	813
mapRealToModel[989	] = 	814
mapRealToModel[990	] = 	815
mapRealToModel[991	] = 	816
mapRealToModel[992	] = 	817
mapRealToModel[993	] = 	818
mapRealToModel[994	] = 	819
mapRealToModel[995	] = 	820
mapRealToModel[996	] = 	821
mapRealToModel[997	] = 	822
mapRealToModel[998	] = 	823
mapRealToModel[999	] = 	824
mapRealToModel[1000	] = 	825
mapRealToModel[1001	] = 	826
mapRealToModel[1002	] = 	827
mapRealToModel[1003	] = 	828
mapRealToModel[1004	] = 	829
mapRealToModel[1005	] = 	830
mapRealToModel[1006	] = 	831
mapRealToModel[1007	] = 	832
mapRealToModel[1008	] = 	833
mapRealToModel[1009	] = 	834
mapRealToModel[1010	] = 	835
mapRealToModel[1011	] = 	836
mapRealToModel[1012	] = 	837
mapRealToModel[1013	] = 	838
mapRealToModel[1014	] = 	839
mapRealToModel[1015	] = 	840
mapRealToModel[1016	] = 	841
mapRealToModel[1017	] = 	842
mapRealToModel[1018	] = 	843
mapRealToModel[1019	] = 	844
mapRealToModel[1020	] = 	845
mapRealToModel[1021	] = 	846
mapRealToModel[1022	] = 	847
mapRealToModel[1023	] = 	848
mapRealToModel[1024	] = 	849
mapRealToModel[1025	] = 	850
mapRealToModel[1026	] = 	851
mapRealToModel[1027	] = 	852
mapRealToModel[1028	] = 	853
mapRealToModel[1029	] = 	854
mapRealToModel[1030	] = 	855
mapRealToModel[1031	] = 	856
mapRealToModel[1032	] = 	857
mapRealToModel[1033	] = 	858
mapRealToModel[1034	] = 	859
mapRealToModel[1035	] = 	860
mapRealToModel[1036	] = 	861
mapRealToModel[1037	] = 	862
mapRealToModel[1038	] = 	863
mapRealToModel[1039	] = 	864
mapRealToModel[1040	] = 	865
mapRealToModel[1041	] = 	866
mapRealToModel[1042	] = 	867
mapRealToModel[1043	] = 	868
mapRealToModel[1044	] = 	869
mapRealToModel[1045	] = 	870
mapRealToModel[1046	] = 	871
mapRealToModel[1047	] = 	872
mapRealToModel[1048	] = 	873
mapRealToModel[1049	] = 	874
mapRealToModel[1050	] = 	875
mapRealToModel[1051	] = 	876
mapRealToModel[1052	] = 	877
mapRealToModel[1053	] = 	878
mapRealToModel[1054	] = 	879
mapRealToModel[1055	] = 	880
mapRealToModel[1056	] = 	881
mapRealToModel[1057	] = 	882
mapRealToModel[1058	] = 	883
mapRealToModel[1059	] = 	884
mapRealToModel[1060	] = 	885
mapRealToModel[1061	] = 	886
mapRealToModel[1062	] = 	887
mapRealToModel[1063	] = 	888
mapRealToModel[1064	] = 	889
mapRealToModel[1065	] = 	890
mapRealToModel[1066	] = 	891
mapRealToModel[1067	] = 	892
mapRealToModel[1068	] = 	893
mapRealToModel[1069	] = 	894
mapRealToModel[1070	] = 	895
mapRealToModel[1071	] = 	896
mapRealToModel[1072	] = 	897
mapRealToModel[1073	] = 	898
mapRealToModel[1074	] = 	899
mapRealToModel[1075	] = 	900
mapRealToModel[1076	] = 	901
mapRealToModel[1077	] = 	902
mapRealToModel[1078	] = 	903
mapRealToModel[1079	] = 	904
mapRealToModel[1080	] = 	905
mapRealToModel[1081	] = 	906
mapRealToModel[1082	] = 	907
mapRealToModel[1083	] = 	908
mapRealToModel[1084	] = 	909
mapRealToModel[1085	] = 	910
mapRealToModel[1086	] = 	911
mapRealToModel[1087	] = 	912
mapRealToModel[1088	] = 	913
mapRealToModel[1089	] = 	914
mapRealToModel[1090	] = 	915
mapRealToModel[1091	] = 	916
mapRealToModel[1092	] = 	917
mapRealToModel[1093	] = 	918
mapRealToModel[1094	] = 	919
mapRealToModel[1095	] = 	920
mapRealToModel[1096	] = 	921
mapRealToModel[1097	] = 	922
mapRealToModel[1098	] = 	923
mapRealToModel[1099	] = 	924
mapRealToModel[1100	] = 	925
mapRealToModel[1101	] = 	926
mapRealToModel[1102	] = 	927
mapRealToModel[1103	] = 	928
mapRealToModel[1104	] = 	929
mapRealToModel[1105	] = 	930
mapRealToModel[1106	] = 	931
mapRealToModel[1107	] = 	932
mapRealToModel[1108	] = 	933
mapRealToModel[1109	] = 	934
mapRealToModel[1110	] = 	935
mapRealToModel[1111	] = 	936
mapRealToModel[1112	] = 	937
mapRealToModel[1113	] = 	938
mapRealToModel[1114	] = 	939
mapRealToModel[1115	] = 	940
mapRealToModel[1116	] = 	941
mapRealToModel[1117	] = 	942
mapRealToModel[1118	] = 	943
mapRealToModel[1119	] = 	944
mapRealToModel[1120	] = 	945
mapRealToModel[1121	] = 	946
mapRealToModel[1122	] = 	947
mapRealToModel[1123	] = 	948
mapRealToModel[1124	] = 	949
mapRealToModel[1125	] = 	950
mapRealToModel[1126	] = 	951
mapRealToModel[1127	] = 	952
mapRealToModel[1128	] = 	953
mapRealToModel[1129	] = 	954
mapRealToModel[1130	] = 	955
mapRealToModel[1131	] = 	956
mapRealToModel[1132	] = 	957
mapRealToModel[1133	] = 	958
mapRealToModel[1134	] = 	959
mapRealToModel[1135	] = 	960
mapRealToModel[1136	] = 	961
mapRealToModel[1137	] = 	962
mapRealToModel[1138	] = 	963
mapRealToModel[1139	] = 	964
mapRealToModel[1140	] = 	965
mapRealToModel[1141	] = 	966
mapRealToModel[1142	] = 	967
mapRealToModel[1143	] = 	968
mapRealToModel[1144	] = 	969
mapRealToModel[1145	] = 	970
mapRealToModel[1146	] = 	971
mapRealToModel[1147	] = 	972

main(sys.argv)

