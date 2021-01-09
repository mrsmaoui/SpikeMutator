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
from spike_closed import mapRealToModel
from spike_open import mapRealToModelOpenA
from spike_open import mapRealToModelOpenB
from spike_open import mapRealToModelOpenC
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

	# Closed structure
	fasta = open("spikes/closed/spike.fasta", "r")
	sequence = fasta.readline().strip()

	# Open structure
	fasta = open("spikes/open/spike_A_open.fasta", "r")
	sequence_A_open = fasta.readline().strip()

	fasta = open("spikes/open/spike_B_open.fasta", "r")
	sequence_B_open = fasta.readline().strip()

	fasta = open("spikes/open/spike_C_open.fasta", "r")
	sequence_C_open = fasta.readline().strip()

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

	if(sequence == "" or sequence_A_open == "" or sequence_B_open == "" or sequence_C_open == ""):
		print error
		print "Fasta file(s) did not load properly. Make sure they exist in the spikes/ directory"
		sys.exit(1)

	if(mutation != "" and mutation != "all" and ',' not in mutation):
		mutation = mutation.upper()
		#check that the mutation is a valid character
		mutationList = list("ARNDCQEGHILKMFPSTWYV")
		if(mutation not in mutationList):
			print "ERROR: Incorrect mutation choice"
			sys.exit(0)	

	if "closed" not in conformation and "open" not in conformation:
		print error
		print "Kindly specify a conformation to use for the spike by using the -c flag: -c closed, or -c open"
		sys.exit(1)

	if "energy" not in action and "mutation" not in action:
		print error
		print "You can only use the value 'energy' with the -a flag"
		sys.exit(1)

	#Determine which USAGE the user intends to perform
	if( (position != '') and (position != 'all') and ',' not in position and int(position) > 0 and len(mutation) == 1):
		
		if conformation == "open":
			params = {"filename":filename, "position":position, "mutation":mutation, 
				"sequence_A":sequence_A_open, "sequence_B":sequence_B_open, "sequence_C":sequence_C_open, "offset":offset, "action":action}  	
			return mutateOpen(params)

		if conformation == "closed":
			params = {"filename":filename, "position":position, "mutation":mutation, "sequence":sequence, "offset":offset, "action":action}  	
			return mutateClosed(params)

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
	

def mutateClosed(params):

	filename = params['filename']
	mutation = params['mutation']
	position = int(params['position'])
	sequence = params['sequence']
	offset = params['offset']
	action = params['action']
	
	# Check if position is available in the 6VXX structure
	model_pos = mapRealToModel[position]
	if model_pos == 0:
		print("The selected position (%d) cannot be tested because it does not exist in the 6VXX protein structure" % model_pos)
		exit()

	#Perform Mutation on original sequence
	mutant_seq = mutateSequence(model_pos, mutation, sequence, offset)
	
	print(simulateClosed(mutation, mutant_seq, position, filename, action))

	os.chdir(SPIKE_PATH)
	

def mutateOpen(params):

	filename = params['filename']
	mutation = params['mutation']
	position = int(params['position'])
	sequence_A_open = params['sequence_A']
	sequence_B_open = params['sequence_B']
	sequence_C_open = params['sequence_C']
	offset = params['offset']
	action = params['action']
	
	# Check if position is available in the 6VYB structure
	model_posA = mapRealToModelOpenA[position]
	model_posB = mapRealToModelOpenB[position]
	model_posC = mapRealToModelOpenC[position]

	if model_posA == 0 or model_posB == 0 or model_posC == 0:
		print("The selected position (%d) cannot be tested because it does not exist in one or more of the the 6VYB protein chain structures" % model_pos)
		exit()

	#Perform Mutation on original sequence
	mutant_seq_A = mutateSequence(model_posA, mutation, sequence_A_open, offset)
	mutant_seq_B = mutateSequence(model_posB, mutation, sequence_B_open, offset)
	mutant_seq_C = mutateSequence(model_posC, mutation, sequence_C_open, offset)
	
	print(simulateOpen(mutation, mutant_seq_A, mutant_seq_B, mutant_seq_C, position, filename, action))

	os.chdir(SPIKE_PATH)


def simulateClosed(mutation, sequence, position, filename, action):
	os.chdir(SPIKE_PATH)	
	timestr = time.strftime("%Y%m%d-%H%M%S")
	usage_path = SPIKE_PATH + '/tmp/' + timestr

	#Create a workspace directory
	call('mkdir %s' % usage_path, shell=True)
	call('mkdir %s/trash' % usage_path, shell=True)

	#Generate .fasta file 
	call('echo %s > %s/%s_%d-%s.fasta' % (sequence, usage_path, filename,position,mutation), shell=True)
	
	#Perform mutation and generate structure file
	
	call("%s/Scwrl4 -i %s/spikes/closed/spike_A.pdb -o %s/scwrl_spike_A_%d-%s.pdb -s %s/%s_%d-%s.fasta >/dev/null" % (SCWRL4_PATH, SPIKE_PATH, usage_path, position, mutation, usage_path, filename, position, mutation), shell=True)

	call("%s/Scwrl4 -i %s/spikes/closed/spike_B.pdb -o %s/scwrl_spike_B_%d-%s.pdb -s %s/%s_%d-%s.fasta >/dev/null" % (SCWRL4_PATH, SPIKE_PATH, usage_path, position, mutation, usage_path, filename, position, mutation), shell=True)

	call("%s/Scwrl4 -i %s/spikes/closed/spike_C.pdb -o %s/scwrl_spike_C_%d-%s.pdb -s %s/%s_%d-%s.fasta >/dev/null" % (SCWRL4_PATH, SPIKE_PATH, usage_path, position, mutation, usage_path, filename, position, mutation), shell=True)

	call("cat %s/scwrl_spike_A_%d-%s.pdb %s/scwrl_spike_B_%d-%s.pdb %s/scwrl_spike_C_%d-%s.pdb > %s/spike_%d-%s.pdb" % (usage_path, position, mutation, usage_path, position, mutation, usage_path, position, mutation, usage_path, position, mutation), shell=True)


	# Copy structure file into runs/ directory
	call("cp %s/spike_%d-%s.pdb %s/" % (usage_path, position, mutation, (SPIKE_PATH + '/structs/spike_closed_%d-%s.pdb'), position, mutation), shell=True)
	
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

	return "Mutated Structure created at %s/spike_closed_%d-%s.pdb" % ( (SPIKE_PATH + '/structs'), position, mutation)


def simulateOpen(mutation, sequence_A, sequence_B, sequence_C, position, filename, action):
	os.chdir(SPIKE_PATH)	
	timestr = time.strftime("%Y%m%d-%H%M%S")
	usage_path = SPIKE_PATH + '/tmp/' + timestr

	#Create a workspace directory
	call('mkdir %s' % usage_path, shell=True)
	call('mkdir %s/trash' % usage_path, shell=True)

	#Generate .fasta file 
	call('echo %s > %s/%s_A_%d-%s.fasta' % (sequence_A, usage_path, filename,position,mutation), shell=True)
	call('echo %s > %s/%s_B_%d-%s.fasta' % (sequence_B, usage_path, filename,position,mutation), shell=True)
	call('echo %s > %s/%s_C_%d-%s.fasta' % (sequence_C, usage_path, filename,position,mutation), shell=True)
	
	#Perform mutation and generate structure file
	
	call("%s/Scwrl4 -i %s/spikes/open/spike_A_open.pdb -o %s/scwrl_spike_A_%d-%s.pdb -s %s/%s_A_%d-%s.fasta >/dev/null" % (SCWRL4_PATH, SPIKE_PATH, usage_path, position, mutation, usage_path, filename, position, mutation), shell=True)

	call("%s/Scwrl4 -i %s/spikes/open/spike_B_open.pdb -o %s/scwrl_spike_B_%d-%s.pdb -s %s/%s_B_%d-%s.fasta >/dev/null" % (SCWRL4_PATH, SPIKE_PATH, usage_path, position, mutation, usage_path, filename, position, mutation), shell=True)

	call("%s/Scwrl4 -i %s/spikes/open/spike_C_open.pdb -o %s/scwrl_spike_C_%d-%s.pdb -s %s/%s_C_%d-%s.fasta >/dev/null" % (SCWRL4_PATH, SPIKE_PATH, usage_path, position, mutation, usage_path, filename, position, mutation), shell=True)

	call("cat %s/scwrl_spike_A_%d-%s.pdb %s/scwrl_spike_B_%d-%s.pdb %s/scwrl_spike_C_%d-%s.pdb > %s/spike_%d-%s.pdb" % (usage_path, position, mutation, usage_path, position, mutation, usage_path, position, mutation, usage_path, position, mutation), shell=True)


	# Copy structure file into runs/ directory
	call("cp %s/spike_%d-%s.pdb %s/" % (usage_path, position, mutation, (SPIKE_PATH + '/structs/spike_open_%d-%s.pdb'), position, mutation), shell=True)
	
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

	return "Mutated Structure created at %s/spike_open_%d-%s.pdb" % ( (SPIKE_PATH + '/structs'), position, mutation)


main(sys.argv)

