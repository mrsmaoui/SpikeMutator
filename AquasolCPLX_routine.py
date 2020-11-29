#!/usr/bin/env python
import os
import math
import string
import fileinput
import tempfile
import sys
from subprocess import call
import os

def energyInWater(AQUASOL_PATH, MAPOR_PATH, temp, pdb_File, pqr_arg):
	#print '\t--> CREATING AQUASOL %s_pbl.in FILE IN WATER' % pqr_arg
        template = open(MAPOR_PATH + '/mdp/PBL_template.in', 'r')
        template_lines = template.readlines()
        template.close()
        template_lines[16] = "File Input Mol     :   %s.pqr                   ! Input file for charges, PDB or PQR format\n" % pqr_arg
        template_lines[21] = "File Output LOG    :   %s_pbl.log               ! Output log file\n" % pqr_arg
        template_lines[36] = "Temperature        :    %s.00                    ! Temperature for the run\n" % temp

        file_pbl = open(pqr_arg+'_pbl.in', 'w')
        for line in template_lines:
                file_pbl.write(line)
        file_pbl.close()
        #print '\t--> AQUASOL %s_pbl.in FILE CREATED' % pqr_arg


        #print '\t--> RUNNING AQUASOL ON %s AND OBTAINING GRID ENERGY' % pqr_arg
        #call("/data/pkgs/AquaSol_Complexes/bin/AquaSolCplx.gfort " + ('%s_pbl.in ' % (pqr_arg)), shell=True)

        call(AQUASOL_PATH + "/AquaSolCplx.gfort " + ('%s_pbl.in > %s_pbl.res' % (pqr_arg, pqr_arg)), shell=True)

        fileGrid = open(pqr_arg+'_pbl.res', 'r')
        fileLines = fileGrid.readlines()
        gridE = 0.0
        units = ''
        waterM = 0.0
        for one_line in fileLines:
                ones = string.split(one_line)
                if(len(ones) ==  10):
                        if (ones[0] + ' ' + ones[1] + ' ' + ones[2]) == 'Free energy (total)':
                                gridE = float(ones[4])
                                units = ones[5]

                if(len(ones) == 5):
                        if(ones[0] + ' ' + ones[1] + ' ' + ones[2] + ' ' + ones[3]) == 'Number of water molecule:':
                                waterM = float(ones[4])
                                break

        fileGrid.close()
        #print '\t--> GRID FREE ENERGY CALCULATED = %s %s %s' % (gridE, units, waterM)

        # Cleaning files
        call("rm -f " + '%s_pbl.in %s_pbl.res %s_pbl.log' % (pqr_arg, pqr_arg, pqr_arg), shell=True)
        rounded_energy = "%0.2f %0.2f" % (gridE, waterM)
	return rounded_energy

def energyInVacuum(AQUASOL_PATH, MAPOR_PATH, temp, pdb_File, pqr_arg):
	#print '\t--> CREATING AQUASOL %s_pbl.in IN VACUUM FILE' % pqr_arg
        template = open(MAPOR_PATH + '/mdp/PBL_template_vacuum.in', 'r')
        template_lines = template.readlines()
        template.close()
        template_lines[16] = "File Input Mol     :   %s.pqr                   ! Input file for charges, PDB or PQR format\n" % pqr_arg
        template_lines[21] = "File Output LOG    :   %s_pbl.log               ! Output log file\n" % pqr_arg
        template_lines[36] = "Temperature        :    %s.00                    ! Temperature for the run\n" % temp

        file_pbl = open(pqr_arg+'_pbl.in', 'w')
        for line in template_lines:
                file_pbl.write(line)
        file_pbl.close()
        #print '\t--> AQUASOL %s_pbl.in FILE CREATED' % pqr_arg


        #print '\t--> RUNNING AQUASOL ON %s AND OBTAINING GRID ENERGY' % pqr_arg
        #call("/data/pkgs/AquaSol_Complexes/bin/AquaSolCplx.gfort " + ('%s_pbl.in ' % (pqr_arg)), shell=True)

        call(AQUASOL_PATH + "/AquaSolCplx.gfort " + ('%s_pbl.in > %s_pbl.res' % (pqr_arg, pqr_arg)), shell=True)

        fileGrid = open(pqr_arg+'_pbl.res', 'r')
        fileLines = fileGrid.readlines()
        gridE = 0.0
        units = ''
        for one_line in fileLines:
                ones = string.split(one_line)
                if(len(ones) ==  10):
                        if (ones[0] + ' ' + ones[1] + ' ' + ones[2]) == 'Free energy (total)':
                                gridE = float(ones[4])
                                units = ones[5]
				break	

        fileGrid.close()
        #print '\t--> GRID FREE ENERGY CALCULATED = %s %s' % (gridE, units)

        # Cleaning files
        call("rm -f " + '%s_pbl.in %s_pbl.res %s_pbl.log %s.pqr' % (pqr_arg, pqr_arg, pqr_arg, pqr_arg), shell=True)
        rounded_energy = "%0.2f" % (gridE)
	return rounded_energy



def applyAquaSol(AQUASOL_PATH, usage_path, pdb2pqr_path, MAPOR_PATH, temp, pdb_File, pqr_arg):
	os.chdir(usage_path)
        call(pdb2pqr_path + "/pdb2pqr.py --ff=AMBER " + pdb_File + ' ' + pqr_arg + '.pqr > test.test', shell=True)
        #call(pdb2pqr_path + "/pdb2pqr.py --ff=AMBER " + pdb_File + ' ' + pqr_arg + '.pqr > test.test >/dev/null 1>/dev/null 2>/dev/null', shell=True)

	e1 = string.split(energyInWater(AQUASOL_PATH, MAPOR_PATH, temp, pdb_File, pqr_arg))
	e2 = energyInVacuum(AQUASOL_PATH, MAPOR_PATH, temp, pdb_File, pqr_arg)

	return e1[0] + " " + e2 + " " + e1[1]
