#!/opt/local/bin//python
from numpy import *
import gzip
import scripts
import scripts1
import scripts2
import scripts3
import constants
import subprocess
import sys

def history_maker(ids,positions,file_string,path,binary):
	"""creates the full history of interesting stars
	ids: array containing the ids of interesting stars
	binintfile: binint file
	collfile: collision log file
	mergefile: merger file
	filestring: the filestring from cmc
	binary: whether there was any primordial binaries or not -> 0: no binary, 1: yes binary"""
	binintfile=path+file_string+'.binint.log'
	collfile=path+file_string+'.collision.log'
	mergefile=path+file_string+'.semergedisrupt.log'
	history_interactions={}
	#make the collision dictionary by reading and sorting the collision log file
	coll=scripts1.collision(collfile)
	#make the pure stellar evolution dictionary
	se=scripts2.bse_int(mergefile)

	if binary:
		#make the binint dictionary by reading from the binint log file
		binint=scripts3.read_binint(binintfile) 


	#now go through the ids of interest and populate the history file with every kind of interactions
	for i in range(len(ids)):
		if binary:
			#read the binint history of the id in this dictionary
			binint_id=scripts3.binary_interaction(binint,ids[i])
			bininteract={'binint': binint_id,
				'binint_coll': {}
				}
			#cross_correlate with collisions: if binint resulted in a merger this gives if any of the stars had any collisions before this
			if len(binint_id.keys())>0:
				for j in binint_id.keys():
					#print j
					if binint_id[j]['merge']==1:
						for k in range(len(binint_id[j]['mergeids'])):
							cross_id=int(binint_id[j]['mergeids'][k])
							#print cross_id
							binint_coll=scripts1.call_collision_tree(coll,cross_id)
							bininteract['binint_coll'][k]=binint_coll
		else:
			bininteract = {}
			
		collision=scripts1.call_collision_tree(coll,ids[i])

		if se.has_key(ids[i]):
			se_id=se[ids[i]]
		else:
			se_id={}

		history_interactions[ids[i]]={'binint': bininteract,
				'coll': collision,
				'se': se_id,
				'position': positions[i]
				}
	
	return history_interactions


def print_files_of_history_for_collisions(star_ids, positions, filestring, binary, prefixstring='coll_history'):
	"""Takes one star id and finds its dynamical history.  Once the history is obtained, it writes this history in a verbose 
	format in files named as the prefixstring+star id.  """
	import numpy as np
	units = scripts.read_units(filestring)
	h = history_maker(star_ids,positions,filestring,binary)
	for i in xrange(len(star_ids)):
		if len(h[star_ids[i]]['coll'])>0:  #some collision did actually take place
			pos_end = positions[i] * units[0]['l_pc']
			writefilename = filestring+'.'+prefixstring+'.'+str(star_ids[i])+'.dat'
			writefile = open(writefilename, 'w')
			writefile.write("#1.time(Myr) 2.final_mass(Msun) 3.final_id 4.interaction 5.position_of_interaction(pc) 6.position_at_end_of_sim(pc) 7.parent_masses1(Msun) 8.parent_ids1 9.parent_masses2(Msun) 10.parent_ids2 11.etc....  \n")
			t_myr, pos_interaction, interaction, par_masses, par_ids, final_mass,final_id = [], [], [], [], [], [], [], 
			for j in h[star_ids[i]]['coll'].keys():
				t_myr.append( h[star_ids[i]]['coll'][j]['coll_params']['time'] * units[0]['t_myr'] )
				print t_myr
				pos_interaction.append( h[star_ids[i]]['coll'][j]['coll_params']['position'] * units[0]['l_pc'] )
				interaction.append( h[star_ids[i]]['coll'][j]['coll_params']['interaction'] )
				par_masses.append( h[star_ids[i]]['coll'][j]['coll_params']['parents']['masses'] )
				par_ids.append( h[star_ids[i]]['coll'][j]['coll_params']['parents']['IDs'] )
				final_mass.append( h[star_ids[i]]['coll'][j]['coll_params']['mass'] )
				final_id.append(j)
			sorted = np.argsort(t_myr)
			for k in xrange(len(sorted)):
				writefile.write("%f %f %ld %s %f %f " %( t_myr[sorted[k]], final_mass[sorted[k]], final_id[sorted[k]], interaction[sorted[k]], pos_interaction[sorted[k]], pos_end, ))
				for l in xrange(len(par_masses[sorted[k]])):
					writefile.write("%f %ld " %( par_masses[sorted[k]][l], par_ids[sorted[k]][l], ))
				writefile.write("\n")
			writefile.close()

					 	
def print_files_of_history_for_collisions_all_sims(basename, nstring='c'):
	import create_table
	import glob
	dirs = dirs = create_table.get_the_dirs(basename, nstring)
	for i in xrange(len(dirs)):
		filestring = dirs[i]+'/initial'
		print filestring
		snapstring = filestring+'.snap*.dat.gz'
		snaps = glob.glob(snapstring)
		sorted_snaps = sort(snaps)
		print sorted_snaps
		if len(sorted_snaps)>0:
			snapno = sorted_snaps[-1].split('/')[-1].split('snap')[1].split('.dat')[0]
			print 'reading snapno', snapno
			units = scripts.read_units(filestring)
			temp_t_myr = scripts.find_t_myr(filestring, snapno)
			if temp_t_myr >= 9000.:
				binary = 0
				binstring = filestring+'.bin.dat'
				binfiles = glob.glob(binstring)
				if len(binfiles)>0:
					binary = 1
					print 'there are binaries: setting binary=', binary
				else:
					binary=0
					print 'all singles: setting binary=', binary
				print 't_myr', temp_t_myr
				temp_bss = scripts.find_BSS(filestring,snapno, 0.001, 1.1)
				ids, pos = temp_bss[3], temp_bss[4]
				print 'printing collision history files for BSs in', dirs[i]
				print_files_of_history_for_collisions(ids, pos, filestring, binary, 'coll_history_test')
	
				




		


#def classifying_BSS(h, star_id, binary, m_cut, t_now, filestring, z):
#	"""
#	h: history dictionary made by history maker
#	star_id: id of the BSS that is of interest
#	binary: were their primordial binaries in the simulation?
#	m_cut: the mass cut above which we define MS stars to be BSSs for t=t_now
#	t_now: last snapshot when BSSs are obtained, t_now is in code units
#	divide all of the BSSs in the following categories depending on what made it a BSS:
#		All collisions
#		SS collision
#		SE driven merger or disruption
#		binary evolution mass transfer
#		note time when it became a BSS and what was the interaction that made it
#		to do that:
#		find the m_cut at t=t_now
#		go back in time and find after what interaction m_star >= m_cut: this is the interaction that 
#			made the BSS
#			get the time and position for the interaction and also the type of the interaction
#		When the star is not found to have a binary interaction, collision, or merger and disruption then 
#			the only other way it can grow in mass is if it had a clean undisturbed mass transfer in a 
#			binary
#			in that case, I need to get these masses and evolve them using BSE to find the onset time 
#			of the mass transfer.  This will be the time when the BSS started to be created.
#	"""
#	print 'star_id', star_id, 't_now', t_now
#	units = scripts.read_units(filestring)
#	bss_coll, bss_ss_coll, bss_bs_coll, bss_bb_coll = 0, 0, 0, 0
#	bss_se, bss_se_merger, bss_se_disruption = 0, 0, 0
#	bss_had_binint = 0
#	bss_se_binint, bss_mtb_binint = 0, 0
#	bss_mtb, bss_mtb_pure = 0, 0
#	SS_coll = 0
#	BB_coll = 0
#	BS_coll = 0
#	t_dict = {}
#	m_dict = {}
#	interaction_dict = {}
#	position_dict = {}
#	channel_dict = {}
#	count = 0
#	#first look for collisions
#	for i in h[star_id]['coll'].keys():
#		m_coll = h[star_id]['coll'][i]['coll_params']['mass']
#		t_coll = h[star_id]['coll'][i]['coll_params']['time']
#		if m_coll >= m_cut:
#			count += 1
#			m_dict[count] = m_coll
#			t_dict[count] = t_coll
#			interaction_dict[count] = h[star_id]['coll'][i]['coll_params']['interaction']
#			position_dict[count] = h[star_id]['coll'][i]['coll_params']['position']
#			channel_dict[count] = 'collision'
#			print channel_dict[count], t_dict[count], star_id
#
#	#then for SE merger and disruption
#	if len(h[star_id]['se'].keys())>0:
#		m_se = h[star_id]['se']['mass']
#		t_se = h[star_id]['se']['time']
#		if m_se >= m_cut:
#			count += 1
#			m_dict[count] = m_se
#			t_dict[count] = t_se
#			interaction_dict[count] = h[star_id]['se']['description']
#			position_dict[count] = h[star_id]['se']['position']
#			channel_dict[count] = 'se'
#			print channel_dict[count], t_dict[count], star_id
#
#	#Now figure out which of these interactions actually made the BSS
#	print interaction_dict
#	t_earliest = t_now+1. #initialization: the plus 1 is to avoid the rare event where on the last timestep 
#				#the BSS is formed.  In that case due to round off, t_now will seem like a later 
#				#time compared to the BSS formation time
#	if len(interaction_dict.keys())>0:
#		print 'len', len(interaction_dict.keys()), 't_now', t_now
#		for i in interaction_dict.keys():
#			if t_dict[i] < t_earliest:
#				#finding which interaction first made it a BSS if there were multiple interactions
#				actual = i
#				t_earliest = t_dict[i]
#	else:
#		actual = 0   #this should happen for any MTB: no coll/se interaction is found
#	print 'actual', actual
#
#
#	if actual>0:
#		actual_interaction = interaction_dict[actual]
#		actual_t = t_dict[actual]
#		actual_channel = channel_dict[actual]
#		actual_position = position_dict[actual]
#	else:
#		actual_interaction = 'mtb'
#		actual_t = t_now
#		actual_channel = 'mtb'
#		actual_position = h[star_id]['position']
#
#
#	if actual_channel == 'collision':
#		bss_coll = 1
#		if actual_interaction == 'single-single':
#			bss_ss_coll = 1
#		elif actual_interaction == 'binary-single':
#			bss_bs_coll = 1
#		elif actual_interaction == 'binary-binary':
#			bss_bb_coll = 1
#	if actual_channel == 'se':
#		bss_se = 1
#		if actual_interaction == 'mergers with one star gone':
#			bss_se_merger = 1
#		elif actual_interaction == 'both stars intact':
#			bss_se_disruption = 1
#
#	#check if this star ever had a binary interaction before: may be it is not important
#	if binary==1:
#		if (h[star_id]['binint']['binint'].keys()) > 0:
#			for i in h[star_id]['binint']['binint'].keys():
#				t_binint = h[star_id]['binint']['binint'][i]['interaction']['type']['time'] 
#				if t_binint < t_earliest:
#					bss_had_binint = 1
#
#	#check whether it may have been mass transfer in a binary
#	if bss_coll == 0 and bss_se == 0:
#		bss_mtb = 1
#		############################################################################
#		#now need to find the time when this actually happened if they are pristine
#		############################################################################
#		if bss_had_binint == 0:
#			bss_mtb_pure = 1
#			t_mtb = find_mtb_bss_beginning_time(star_id, z, filestring)
#			t_mtb = t_mtb/units[0]['t_myr']
#			if t_mtb < actual_t:
#				actual_t = t_mtb
#
#	
#	#mixed with binary interaction before SE or MTB: is it important?
#	if bss_se == 1 and bss_had_binint == 1:
#		bss_se_binint = 1
#	if bss_mtb == 1 and bss_had_binint == 1:
#		bss_mtb_binint = 1
#
#	return bss_coll, bss_ss_coll, bss_bs_coll, bss_bb_coll, bss_se, bss_se_merger, bss_se_disruption, bss_had_binint, bss_mtb, bss_mtb_pure, bss_se_binint, bss_mtb_binint, actual_t, actual_position, star_id
			


def classifying_BSS(h, star_id, binary, m_cut, t_now, filestring, z):
	"""
	h: history dictionary made by history maker
	star_id: id of the BSS that is of interest
	binary: were their primordial binaries in the simulation?
	m_cut: the mass cut above which we define MS stars to be BSSs for t=t_now
	t_now: last snapshot when BSSs are obtained, t_now is in code units
	divide all of the BSSs in the following categories depending on what made it a BSS:
		All collisions
		SS collision
		SE driven merger or disruption
		binary evolution mass transfer
		note time when it became a BSS and what was the interaction that made it
		to do that:
		find the m_cut at t=t_now
		go back in time and find after what interaction m_star >= m_cut: this is the interaction that made the BSS
		get the time and position for the interaction and also the type of the interaction
		When the star is not found to have a binary interaction, collision, or merger and disruption then 
			the only other way it can grow in mass is if it had a clean undisturbed mass transfer in a binary
			in that case, I need to get these masses and evolve them using BSE to find the onset time of the mass transfer.  This will be the time when the BSS started to be created.
	"""
	#print 'star_id', star_id, 't_now', t_now
	units = scripts.read_units(filestring)
	bss_coll, bss_ss_coll, bss_bs_coll, bss_bb_coll = 0, 0, 0, 0
	bss_se, bss_se_merger, bss_se_disruption = 0, 0, 0
	bss_had_binint = 0
	bss_se_binint, bss_mtb_binint = 0, 0
	bss_mtb, bss_mtb_pure = 0, 0
	SS_coll = 0
	BB_coll = 0
	BS_coll = 0
	t_dict = {}
	m_dict = {}
	interaction_dict = {}
	position_dict = {}
	channel_dict = {}
	count = 0
	interacted, t_of_interaction1, t_of_interaction2, t_of_interaction3, t_of_interaction = 0, 0., 0., 0., 0.
	#first look for collisions
	print h[star_id]['coll'].keys()
	for i in h[star_id]['coll'].keys():
		print i
		interacted = 1
		m_coll = h[star_id]['coll'][i]['coll_params']['mass']
		t_coll = h[star_id]['coll'][i]['coll_params']['time']
		t_of_interaction1 = t_coll
		if m_coll >= m_cut:
			count += 1
			m_dict[count] = m_coll
			t_dict[count] = t_coll
			interaction_dict[count] = h[star_id]['coll'][i]['coll_params']['interaction']
			position_dict[count] = h[star_id]['coll'][i]['coll_params']['position']
			channel_dict[count] = 'collision'
			print channel_dict[count], t_dict[count], star_id

	#then for SE merger and disruption
	print h[star_id]['se'].keys()
	if len(h[star_id]['se'].keys())>0:
		interacted = 1
		t_of_interaction2 = h[star_id]['se']['time']
		m_se = h[star_id]['se']['mass']
		t_se = h[star_id]['se']['time']
		if m_se >= m_cut:
			count += 1
			m_dict[count] = m_se
			t_dict[count] = t_se
			interaction_dict[count] = h[star_id]['se']['description']
			position_dict[count] = h[star_id]['se']['position']
			channel_dict[count] = 'se'
			print channel_dict[count], t_dict[count], star_id

	#Now figure out which of these interactions actually made the BSS
	#print interaction_dict
	t_earliest = t_now+1. #initialization: the plus 1 is to avoid the rare event where on the last timestep the BSS is formed.  In that case due to round off, t_now will seem like a later 
	#time compared to the BSS formation time
	if len(interaction_dict.keys())>0:
		print 'len', len(interaction_dict.keys()), 't_now', t_now, interaction_dict.keys()
		for i in interaction_dict.keys():
			if t_dict[i] < t_earliest:
				#finding which interaction first made it a BSS if there were multiple interactions
				actual = i
				t_earliest = t_dict[i]
				print actual
	else:
		actual = 0   #this should happen for any MTB: no coll/se interaction is found
	print 'actual', actual


	if actual>0:
		actual_interaction = interaction_dict[actual]
		actual_t = t_dict[actual]
		actual_channel = channel_dict[actual]
		actual_position = position_dict[actual]
	else:
		actual_interaction = 'mtb'
		actual_t = t_now
		actual_channel = 'mtb'
		actual_position = h[star_id]['position']


	if actual_channel == 'collision':
		bss_coll = 1
		if actual_interaction == 'single-single':
			bss_ss_coll = 1
		elif actual_interaction == 'binary-single':
			bss_bs_coll = 1
		elif actual_interaction == 'binary-binary':
			bss_bb_coll = 1
	if actual_channel == 'se':
		bss_se = 1
		if actual_interaction == 'mergers with one star gone':
			bss_se_merger = 1
		elif actual_interaction == 'both stars intact':
			bss_se_disruption = 1

	#check if this star ever had a binary interaction before: may be it is not important
	if binary==1:
		if len(h[star_id]['binint']['binint'].keys()) > 0:
			for i in h[star_id]['binint']['binint'].keys():
				t_binint = h[star_id]['binint']['binint'][i]['interaction']['type']['time']
				t_of_interaction3 = t_binint
				if t_binint < t_earliest:
					bss_had_binint = 1

	#check whether it may have been mass transfer in a binary
	t_of_interaction = max([t_of_interaction1, t_of_interaction2, t_of_interaction3])
	print 'bss_coll', bss_coll, 'bss_se', bss_se, 'bss_had_binint', bss_had_binint, 'interacted', interacted, 't_of_int', t_of_interaction

	try:
		if bss_coll == 0 and bss_se == 0:
			bss_mtb = 1
			actual_t=-100   ##Shi: Because I don't know how to implement bse here
			############################################################################
			#now need to find the time when this actually happened if they are pristine
			############################################################################
			if bss_had_binint == 0 and interacted == 0:
				bss_mtb_pure = 1
			#	t_mtb = find_mtb_bss_beginning_time(star_id, z, filestring, t_of_interaction)
			#	t_mtb = t_mtb/units[0]['t_myr']
			#	if t_mtb < actual_t:
			#		actual_t = t_mtb
			#else:
			#	t_mtb = find_mtb_bss_beginning_time(star_id, z, filestring, t_of_interaction)
			#	t_mtb = t_mtb/units[0]['t_myr']
			#	if t_mtb < actual_t:
			#		actual_t = t_mtb
	except (UnboundLocalError, IndexError):
		print sys.exc_info()[0]
		print sys.exc_info()[1]
		actual_t = t_now + 1.  #these are cases where t_formation is not found so pile them up outside the integration time


	
	#mixed with binary interaction before SE or MTB: is it important?
	if bss_se == 1 and bss_had_binint == 1:
		bss_se_binint = 1
	if bss_mtb == 1 and bss_had_binint == 1:
		bss_mtb_binint = 1

	return bss_coll, bss_ss_coll, bss_bs_coll, bss_bb_coll, bss_se, bss_se_merger, bss_se_disruption, bss_had_binint, bss_mtb, bss_mtb_pure, bss_se_binint, bss_mtb_binint, actual_t, actual_position, star_id



#def find_mtb_bss_beginning_time(star_id, z, filestring):
#	"""finds the binary properties"""
#	snapfile = filestring+'.snap0000.dat.gz'
#	print snapfile
#	snapatzero = gzip.open(snapfile,'rb')
#	star_id_string = ' '+str(int(star_id))+' '
#	print star_id_string
#	try:
#		for line in snapatzero:
#			if line.rfind(star_id_string) > -1:
#				values = line.split()
#				print values
#				#make sure things make sense
#				if float(values[7]) == 0.:
#					print 'cannot be: supposed to be a binary'
#					raise StopIteration()
#				else:
#					m1, m2 = float(values[8]), float(values[9])
#					a = float(values[12])
#					e = float(values[13])
#					R1, R2 = float(values[21]), float(values[22])
#					kstar1, kstar2 = float(values[17]), float(values[18])
#					print m1
#					raise StopIteration()
#	except StopIteration:
#		print 'obtained binary properties'
#		print m1, m2, a, e, R1, R2, kstar1, kstar2
#	
#	#now create the binary.in file for this one
#	tphysf = 12000.
#	tb = scripts.get_period(a*constants.AU, m1*constants.Msun,m2*constants.Msun) 
#	f=open('binary.in' , 'w')
#	f.write("%f %f %f %f %d %d %f %f\n"  %(m1, m2, 12000., tb, kstar1, kstar2, z, e))
#	#f.write("0.5 0.0 1.0 3.0 0.5 0\n0 1 0 1 0 0 1.8 29769\n0.05 0.01 0.02\n190.0 0.125 1.0 1.5 0.001 10.0 -1.0\n1113.89139868144343\n117.3629102921285 4.94763450979127661 0.0583570371042554659\n-282.92588170378781 9.13749733018044452 1.08589687349379727\n# alternatives for top line ...\n7.816 4.387 15000. 1964.18453 1 1 0.02 0.0\n2.9 0.9 15000. 8.0 1 1 0.02 0.7\n")
#	f.write("0.5 0.0 1.0 3.0 0.5 0\n0 1 0 1 0 0 1.8 29769\n0.05 0.01 0.02\n190.0 0.125 1.0 1.5 0.001 10.0 -1.0\n")
#	f.close()
#	dataout,dataerr=subprocess.Popen([r"./bse"],stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
#	a = dataout.split('\n')
#	for i in range(len(a)):
#		if a[i].rfind('BEG RCHE')>-1:
#			print a[i]
#			t_mtb = float(a[i].split()[0])
#			print t_mtb
#	return t_mtb




def find_mtb_bss_beginning_time(star_id, z, filestring, t, BSEEXEC='/Users/shiye/Documents/bse/bse'):
	"""finds the binary properties"""
	from glob import glob
	#finding star after it came into being. For non-collisional or non-merger stars it is at t=0
	#Otherwise it must be at some later time
	#t here is given as the time of collision/merger event that created the star and changed its id
	snapstring = filestring+'.snap*.dat.gz'
	snaps = glob(snapstring)
	print snaps
	snapcounter, found_snap = 0, 0
	while found_snap==0:
		snapfile = snaps[snapcounter]
		snap = gzip.open(snapfile,'rb')
		line = snap.readline()
		t_snap = float(line.split('=')[1].split()[0])
		if t_snap>=t:
			actual_snapfile = snapfile
			found_snap = 1
		snapcounter = snapcounter+1

	#snapfile = filestring+'.snap0000.dat.gz'
	print actual_snapfile
	snapatzero = gzip.open(actual_snapfile,'rb')
	star_id_string = ' '+str(int(star_id))+' '
	print star_id_string
	try:
		for line in snapatzero:
			if line.rfind(star_id_string) > -1:
				values = line.split()
				print values
				#make sure things make sense
				if float(values[7]) == 0.:
					print 'cannot be: supposed to be a binary'
					raise StopIteration()
				else:
					m1, m2 = float(values[8]), float(values[9])
					a = float(values[12])
					e = float(values[13])
					R1, R2 = float(values[21]), float(values[22])
					kstar1, kstar2 = float(values[17]), float(values[18])
					print m1
					raise StopIteration()
	except StopIteration:
		print 'obtained binary properties' 
		print m1, m2, a, e, R1, R2, kstar1, kstar2
	
	#now create the binary.in file for this one
	tphysf = 12000.
	tb = scripts.get_period(a*constants.AU, m1*constants.Msun,m2*constants.Msun) 
	f=open('/Users/shiye/Documents/ClusterGroup/BSSproject/bss_binary.in' , 'w')
	f.write("%f %f %f %f %d %d %f %f\n"  %(m1, m2, 12000., tb, kstar1, kstar2, z, e))
	f.write("0.5 0.0 1.0 3.0 0.5\n0 1 0 1 0 0 1.8 29769\n0.05 0.01 0.02\n190.0 0.125 1.0 1.5 0.001 10.0 -1.0\n")
	f.close()
	dataout,dataerr=subprocess.Popen([BSEEXEC],stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
	print dataout, dataerr
	a = dataout.split('\n')
	for i in range(len(a)):
		if a[i].rfind('BEG RCHE')>-1:
			print a[i]
			t_mtb = float(a[i].split()[0])
			print t_mtb
	return t_mtb



#def classify_all_bss_stars(filestring, snapno, z, mcut):
#	#first get time now from the snap
#	units = scripts.read_units(filestring)
#	snap = filestring+'.snap'+snapno+'.dat.gz'
#	f_snap = gzip.open(snap,'rb')
#	line_snap = f_snap.readline()
#	t_snap = float(line_snap.split()[1].split('=')[1]) 
#	print 't_snap', t_snap
#
#	#find rc at t_snap for core population determination
#	dynfile = filestring+'.dyn.dat'
#	print dynfile
#	dyndata = loadtxt(dynfile)
#	print dyndata[-1]
#	try:
#		for i in range(len(dyndata)-1, -1, -1):   #start from the last line, it will be quicker
#			if dyndata[i,0]==t_snap:
#				rc = dyndata[i,7]
#				raise StopIteration()
#	except StopIteration:
#		print 'found rc', rc, 'at time', t_snap
#
#	#were there primordial binaries or any binary interactions:
#	from glob import glob
#	binintstring=filestring+'.binint.log'
#	#binint=glob('*.binint.log')
#	binint=glob(binintstring)
#	binary=1
#	if len(binint)==0:
#		binary=0
#	print 'binary', binary	
#	BSSs = scripts.find_BSS(filestring,snapno, z, mcut)
#	
#	bss_ids, bss_positions, m_to = BSSs[3], BSSs[4], BSSs[5]
#	h = history_maker(bss_ids, bss_positions, filestring, binary)
#	all_bss = len(bss_ids)
#	#now classify them one by one 
#	#how many are in the core
#	core_count = 0
#	for i in range(len(bss_positions)):
#		if bss_positions[i] < rc:
#			core_count += 1
#	print 'core_count', core_count
#	#initialize the arrays
#	classified = zeros(( len(bss_ids), 15 ))
#	classified_core = zeros(( core_count, 15 ))
#	core_count = 0
#	m_cutoff = m_to*mcut
#	for i in range(len(bss_ids)):
#		classified[i] = classifying_BSS(h, bss_ids[i], binary, m_cutoff, t_snap, filestring, z)
#		if bss_positions[i]<=rc:
#			classified_core[core_count] = classified[i]
#			core_count += 1
#	print classified_core
#
#	file_all = filestring+'.all_bss_formation_properties_one_sim.dat'
#	file_all = open(file_all, 'w')
#	file_all.write("#1.bss_coll 2.bss_ss_coll 3.bss_bs_coll 4.bss_bb_coll 5.bss_se 6.bss_se_merger 7.bss_se_disruption 8.bss_had_binint 9.bss_mtb 10.bss_mtb_pure 11.bss_se_binint 12.bss_mtb_binint 13.actual_t(Myr) 14.actual_position 15.star_id\n")
#	for i in range(len(classified)):
#		for j in range(12):
#			file_all.write("%d " %(classified[i,j],))
#		age = (t_snap-classified[i,12])*units[0]['t_myr']
#		file_all.write("%f " %( age ))
#		file_all.write("%f " %(classified[i,13],))
#		file_all.write("%d\n" %(classified[i,14],))
#	file_all.close()
#
#	file_core = filestring+'.core_bss_formation_properties_one_sim.dat'
#	file_core = open(file_core, 'w')
#	file_core.write("#1.bss_coll 2.bss_ss_coll 3.bss_bs_coll 4.bss_bb_coll 5.bss_se 6.bss_se_merger 7.bss_se_disruption 8.bss_had_binint 9.bss_mtb 10.bss_mtb_pure 11.bss_se_binint 12.bss_mtb_binint 13.actual_t(Myr) 14.actual_position 15.star_id\n")
#	for i in range(len(classified_core)):
#		for j in range(12):
#			file_core.write("%d " %(classified_core[i,j],))
#		file_core.write("%f " %((t_snap-classified_core[i,12])*units[0]['t_myr'] ))
#		file_core.write("%f " %(classified_core[i,13],))
#		file_core.write("%d\n" %(classified_core[i,14],))
#	file_core.close()
#
#	
#
#	total_coll = sum(classified[:,0])
#	total_ss_coll = sum(classified[:,1])
#	total_bs_coll = sum(classified[:,2])
#	total_bb_coll = sum(classified[:,3])
#	total_se = sum(classified[:,4])
#	total_se_merger = sum(classified[:,5])
#	total_se_disruption = sum(classified[:,6])
#	total_had_binint = sum(classified[:,7])
#	total_mtb = sum(classified[:,8])
#	total_mtb_pure = sum(classified[:,9])
#	total_se_binint = sum(classified[:,10])
#	total_mtb_binint = sum(classified[:,11])
#	largest_dynamic_t = t_snap - min(classified[:,12])
#	median_dynamic_age = t_snap - median(classified[:,12])
#	print classified[:,12]
#
#	total_core_coll = sum(classified_core[:,0])
#	total_core_ss_coll = sum(classified_core[:,1])
#	total_core_bs_coll = sum(classified_core[:,2])
#	total_core_bb_coll = sum(classified_core[:,3])
#	total_core_se = sum(classified_core[:,4])
#	total_core_se_merger = sum(classified_core[:,5])
#	total_core_se_disruption = sum(classified_core[:,6])
#	total_core_had_binint = sum(classified_core[:,7])
#	total_core_mtb = sum(classified_core[:,8])
#	total_core_mtb_pure = sum(classified_core[:,9])
#	total_core_se_binint = sum(classified[:,10])
#	total_core_mtb_binint = sum(classified[:,11])
#	if total_core_coll > 0:
#		largest_dynamic_t_c = t_snap - min(classified_core[:,12])
#		median_dynamic_age_c = t_snap - median(classified_core[:,12])
#	else:
#		largest_dynamic_t_c = -1.
#		median_dynamic_age_c = -1.
#	print classified_core[:,12]
#
#	total = [all_bss, total_coll, total_ss_coll, total_bs_coll, total_bb_coll, total_se, total_se_merger, total_se_disruption, total_had_binint, total_mtb, total_mtb_pure, total_se_binint, total_mtb_binint, largest_dynamic_t, median_dynamic_age]
#
#	total_core = [core_count, total_core_coll, total_core_ss_coll, total_core_bs_coll, total_core_bb_coll, total_core_se, total_core_se_merger, total_core_se_disruption, total_core_had_binint, total_core_mtb, total_core_mtb_pure, total_core_se_binint, total_core_mtb_binint, largest_dynamic_t_c, median_dynamic_age_c]
#
#	return total, total_core

	


def classify_all_bss_stars(filestring, snapno, z, mcut):
	#first get time now from the snap
	units = scripts.read_units(filestring)
	snap = filestring+'.snap'+snapno+'.dat.gz'
	f_snap = gzip.open(snap,'rb')
	line_snap = f_snap.readline()
	t_snap = float(line_snap.split()[1].split('=')[1]) 
	print 't_snap', t_snap
	
	#find rc at t_snap for core population determination
	dynfile = filestring+'.dyn.dat'
	print dynfile
	dyndata = loadtxt(dynfile)
	print dyndata[-1]
	try:
		for i in range(len(dyndata)-1, -1, -1):   #start from the last line, it will be quicker
			if dyndata[i,0]==t_snap:
				rc = dyndata[i,7]
				raise StopIteration()
	except StopIteration:
		print 'found rc', rc, 'at time', t_snap

	#were there primordial binaries or any binary interactions:
	from glob import glob
	binintstring=filestring+'.binint.log'
	#binint=glob('*.binint.log')
	binint=glob(binintstring)
	binary=1
	if len(binint)==0:
		binary=0
	print 'binary', binary
	print 'calling scripts.find_BSS'
	BSSs = scripts.find_BSS(filestring,snapno, z, mcut)
	
	bss_ids, bss_positions, m_to = BSSs[3], BSSs[4], BSSs[5]
	print 'calling history_maker'
	h = history_maker(bss_ids, bss_positions, filestring, binary)
	all_bss = len(bss_ids)
	#now classify them one by one 
	#how many are in the core
	core_count = 0
	for i in range(len(bss_positions)):
		if bss_positions[i] < rc:
			core_count += 1
	print 'core_count', core_count
	#initialize the arrays
	classified = zeros(( len(bss_ids), 15 ))
	classified_core = zeros(( core_count, 15 ))
	core_count = 0
	m_cutoff = m_to*mcut
	for i in range(len(bss_ids)):
		classified[i] = classifying_BSS(h, bss_ids[i], binary, m_cutoff, t_snap, filestring, z)
		if bss_positions[i]<=rc:
			classified_core[core_count] = classified[i]
			core_count += 1
	print classified_core

	file_all = filestring+'.all_bss_formation_properties_one_sim.dat'
	file_all = open(file_all, 'w')
	file_all.write("#1.bss_coll 2.bss_ss_coll 3.bss_bs_coll 4.bss_bb_coll 5.bss_se 6.bss_se_merger 7.bss_se_disruption 8.bss_had_binint 9.bss_mtb 10.bss_mtb_pure 11.bss_se_binint 12.bss_mtb_binint 13.actual_t(Myr) 14.actual_position 15.star_id\n")
	for i in range(len(classified)):
		for j in range(12):
			file_all.write("%d " %(classified[i,j],))
		age = (t_snap-classified[i,12])*units[0]['t_myr']
		file_all.write("%f " %( age ))
		file_all.write("%f " %(classified[i,13],))
		file_all.write("%d\n" %(classified[i,14],))
	file_all.close()

	file_core = filestring+'.core_bss_formation_properties_one_sim.dat'
	file_core = open(file_core, 'w')
	file_core.write("#1.bss_coll 2.bss_ss_coll 3.bss_bs_coll 4.bss_bb_coll 5.bss_se 6.bss_se_merger 7.bss_se_disruption 8.bss_had_binint 9.bss_mtb 10.bss_mtb_pure 11.bss_se_binint 12.bss_mtb_binint 13.actual_t(Myr) 14.actual_position 15.star_id\n")
	for i in range(len(classified_core)):
		for j in range(12):
			file_core.write("%d " %(classified_core[i,j],))
		file_core.write("%f " %((t_snap-classified_core[i,12])*units[0]['t_myr'] ))
		file_core.write("%f " %(classified_core[i,13],))
		file_core.write("%d\n" %(classified_core[i,14],))
	file_core.close()

	

	total_coll = sum(classified[:,0])
	total_ss_coll = sum(classified[:,1])
	total_bs_coll = sum(classified[:,2])
	total_bb_coll = sum(classified[:,3])
	total_se = sum(classified[:,4])
	total_se_merger = sum(classified[:,5])
	total_se_disruption = sum(classified[:,6])
	total_had_binint = sum(classified[:,7])
	total_mtb = sum(classified[:,8])
	total_mtb_pure = sum(classified[:,9])
	total_se_binint = sum(classified[:,10])
	total_mtb_binint = sum(classified[:,11])
	largest_dynamic_t = t_snap - min(classified[:,12])
	median_dynamic_age = t_snap - median(classified[:,12])
	print classified[:,12]

	total_core_coll = sum(classified_core[:,0])
	total_core_ss_coll = sum(classified_core[:,1])
	total_core_bs_coll = sum(classified_core[:,2])
	total_core_bb_coll = sum(classified_core[:,3])
	total_core_se = sum(classified_core[:,4])
	total_core_se_merger = sum(classified_core[:,5])
	total_core_se_disruption = sum(classified_core[:,6])
	total_core_had_binint = sum(classified_core[:,7])
	total_core_mtb = sum(classified_core[:,8])
	total_core_mtb_pure = sum(classified_core[:,9])
	total_core_se_binint = sum(classified[:,10])
	total_core_mtb_binint = sum(classified[:,11])
	if total_core_coll > 0:
		largest_dynamic_t_c = t_snap - min(classified_core[:,12])
		median_dynamic_age_c = t_snap - median(classified_core[:,12])
	else:
		largest_dynamic_t_c = -1.
		median_dynamic_age_c = -1.
	print classified_core[:,12]

	total = [all_bss, total_coll, total_ss_coll, total_bs_coll, total_bb_coll, total_se, total_se_merger, total_se_disruption, total_had_binint, total_mtb, total_mtb_pure, total_se_binint, total_mtb_binint, largest_dynamic_t, median_dynamic_age]

	total_core = [core_count, total_core_coll, total_core_ss_coll, total_core_bs_coll, total_core_bb_coll, total_core_se, total_core_se_merger, total_core_se_disruption, total_core_had_binint, total_core_mtb, total_core_mtb_pure, total_core_se_binint, total_core_mtb_binint, largest_dynamic_t_c, median_dynamic_age_c]

	return total, total_core




		

def classify_all_bss_stars_all_sims(basename, z, mcut, writefile, problemfile):
	import create_table 
	import glob 
	dirs = create_table.get_the_dirs(basename, 'c4-')
	writefile = open(writefile, 'a')
	problemfile = open(problemfile, 'a')
	writefile.write("#1.all_bss 2.total_coll 3.total_ss_coll 4.total_bs_coll 5.total_bb_coll 6.total_se 7.total_se_merger 8.total_se_disruption 9.total_had_binint 10.total_mtb 11.total_mtb_pure 12.total_se_binint 13.total_mtb_binint 14.age_oldest 15.age_median 16.core_count 17.total_core_coll 18.total_core_ss_coll 19.total_core_bs_coll 20.total_core_bb_coll 21.total_core_se 22.total_core_se_merger 23.total_core_se_disruption 24.total_core_had_binint 25.total_core_mtb 26.total_core_mtb_pure 27.total_core_se_binint 28.total_core_mtb_binint 29.age_oldest_core 30.age_median_core 31.run\n")
	
	for i in range(len(dirs)):
		string = dirs[i]+'/initial'
		print 'string'
		snapstring = string+'.snap*.dat.gz'
		snaps = glob.glob(snapstring)
		sorted_snaps = sort(snaps)
		print sorted_snaps
		
		if len(snaps)>0:
			snapno = sorted_snaps[-1].split('/')[-1].split('snap')[1].split('.dat')[0]
			print snapno
			
			try:	
				units = scripts.read_units(string)
				temp_t_myr = scripts.find_t_myr(string, snapno)
				if temp_t_myr >= 9000.:
					print 't_myr', temp_t_myr
					total, total_core = classify_all_bss_stars(string, snapno, z, mcut)
	
					all_bss = total[0]
					total_coll, total_ss_coll, total_bs_coll, total_bb_coll = total[1], total[2], total[3], total[4]
					total_se, total_se_merger, total_se_disruption = total[5], total[6], total[7]
					total_had_binint = total[8]
					total_mtb, total_mtb_pure = total[9], total[10]
					total_se_binint, total_mtb_binint = total[11], total[12]
					age_oldest, age_median = total[13]*units[0]['t_myr'], total[14]*units[0]['t_myr']
	
					writefile.write( "%d " %(all_bss,) )
					writefile.write( "%d %d %d %d " %(total_coll, total_ss_coll, total_bs_coll, total_bb_coll,) )
					writefile.write( "%d %d %d " %(total_se, total_se_merger, total_se_disruption,) )
					writefile.write( "%d " %(total_had_binint,) )
					writefile.write( "%d %d " %(total_mtb, total_mtb_pure,) )
					writefile.write( "%d %d " %(total_se_binint, total_mtb_binint,) )
					writefile.write( "%f %f " %(age_oldest, age_median,) )
	
					core_count = total_core[0]
					total_core_coll, total_core_ss_coll, total_core_bs_coll, total_core_bb_coll = total_core[1], total_core[2], total_core[3], total_core[4]
					total_core_se, total_core_se_merger, total_core_se_disruption = total_core[5], total_core[6], total_core[7]
					total_core_had_binint = total_core[8]
					total_core_mtb, total_core_mtb_pure = total_core[9], total_core[10]
					total_core_se_binint, total_core_mtb_binint = total_core[11], total_core[12]
					age_oldest_core, age_median_core = total_core[13]*units[0]['t_myr'], total_core[14]*units[0]['t_myr']
	
					writefile.write( "%d " %(core_count,))
					writefile.write( "%d %d %d %d " %(total_core_coll, total_core_ss_coll, total_core_bs_coll, total_core_bb_coll,))
					writefile.write( "%d %d %d " %(total_core_se, total_core_se_merger, total_core_se_disruption,))
					writefile.write( "%d " %(total_core_had_binint,))
					writefile.write( "%d %d " %(total_core_mtb, total_core_mtb_pure,) )
					writefile.write( "%d %d " %(total_core_se_binint, total_core_mtb_binint,) )
					writefile.write( "%f %f " %(age_oldest_core, age_median_core,) )
	
					writefile.write( "%s\n" %(dirs[i],) )
				else:
					problemfile.write("%s: ran only upto %f Myr. Discarded.snapno %s\n" %(dirs[i], temp_t_myr, snapno))


			except (ValueError, IOError):
				problemfile.write("%s\n" %(dirs[i],))

	writefile.close()
	problemfile.close()


def just_print(total, total_core, directorystring, writefilename):
	"""If I want to run only one directory to extract the BS info and then want to print them in the same BS files that gets created after running 
	classify_all_bss_stars_all_sims, then use this nifty routine.  
	1. run classify_all_bss_stars first to obtain total, and total_core
	2. run just_print after supplying the file this line needs to be printed to and the directory specifying the model.  
	"""
	convstring = directorystring+'/initial'
	units = scripts.read_units(convstring)
	writefile = open(writefilename, 'a')  #Note that this is appending to the file.  Be careful
	writefile.write("#1.all_bss 2.total_coll 3.total_ss_coll 4.total_bs_coll 5.total_bb_coll 6.total_se 7.total_se_merger 8.total_se_disruption 9.total_had_binint 10.total_mtb 11.total_mtb_pure 12.total_se_binint 13.total_mtb_binint 14.age_oldest 15.age_median 16.core_count 17.total_core_coll 18.total_core_ss_coll 19.total_core_bs_coll 20.total_core_bb_coll 21.total_core_se 22.total_core_se_merger 23.total_core_se_disruption 24.total_core_had_binint 25.total_core_mtb 26.total_core_mtb_pure 27.total_core_se_binint 28.total_core_mtb_binint 29.age_oldest_core 30.age_median_core 31.run\n")
	all_bss = total[0]
	total_coll, total_ss_coll, total_bs_coll, total_bb_coll = total[1], total[2], total[3], total[4]
	total_se, total_se_merger, total_se_disruption = total[5], total[6], total[7]
	total_had_binint = total[8]
	total_mtb, total_mtb_pure = total[9], total[10]
	total_se_binint, total_mtb_binint = total[11], total[12]
	age_oldest, age_median = total[13]*units[0]['t_myr'], total[14]*units[0]['t_myr']

	writefile.write( "%d " %(all_bss,) )
	writefile.write( "%d %d %d %d " %(total_coll, total_ss_coll, total_bs_coll, total_bb_coll,) )
	writefile.write( "%d %d %d " %(total_se, total_se_merger, total_se_disruption,) )
	writefile.write( "%d " %(total_had_binint,) )
	writefile.write( "%d %d " %(total_mtb, total_mtb_pure,) )
	writefile.write( "%d %d " %(total_se_binint, total_mtb_binint,) )
	writefile.write( "%f %f " %(age_oldest, age_median,) )

	core_count = total_core[0]
	total_core_coll, total_core_ss_coll, total_core_bs_coll, total_core_bb_coll = total_core[1], total_core[2], total_core[3], total_core[4]
	total_core_se, total_core_se_merger, total_core_se_disruption = total_core[5], total_core[6], total_core[7]
	total_core_had_binint = total_core[8]
	total_core_mtb, total_core_mtb_pure = total_core[9], total_core[10]
	total_core_se_binint, total_core_mtb_binint = total_core[11], total_core[12]
	age_oldest_core, age_median_core = total_core[13]*units[0]['t_myr'], total_core[14]*units[0]['t_myr']

	writefile.write( "%d " %(core_count,))
	writefile.write( "%d %d %d %d " %(total_core_coll, total_core_ss_coll, total_core_bs_coll, total_core_bb_coll,))
	writefile.write( "%d %d %d " %(total_core_se, total_core_se_merger, total_core_se_disruption,))
	writefile.write( "%d " %(total_core_had_binint,))
	writefile.write( "%d %d " %(total_core_mtb, total_core_mtb_pure,) )
	writefile.write( "%d %d " %(total_core_se_binint, total_core_mtb_binint,) )
	writefile.write( "%f %f " %(age_oldest_core, age_median_core,) )

	writefile.write( "%s\n" %(directorystring,) )

#def extract_cluster_properties_for_classified_bss(bssfilename, writefilename):
#	"""Once classify_all_bss_stars_all_sims has finished running, each simulation directory should 
#	have a file where the branching ratios are already calculated.  Also the full paths of the simulation 
#	files are there.  This is just a conenient routine to extract some cluster global properties and write 
#	another file with that data also.  Now only interested in rhoc, rc, vcrms to calculate Gamma.  """
#	from glob import glob
#	import create_table as table
#	writefile = open(writefilename, 'w')
#	writefile.write("#1.all_bss 2.total_coll 3.total_ss_coll 4.total_bs_coll 5.total_bb_coll 6.total_se 7.total_se_merger 8.total_se_disruption 9.total_had_binint 10.total_mtb 11.total_mtb_pure 12.total_se_binint 13.total_mtb_binint 14.age_oldest(Myr) 15.age_median(Myr) 16.core_count 17.total_core_coll 18.total_core_ss_coll 19.total_core_bs_coll 20.total_core_bb_coll 21.total_core_se 22.total_core_se_merger 23.total_core_se_disruption 24.total_core_had_binint 25.total_core_mtb 26.total_core_mtb_pure 27.total_core_se_binint 28.total_core_mtb_binint 29.age_oldest_core(Myr) 30.age_median_core(Myr) 31.rc(pc) 32.rhoc(Msun/pc^3) 33.vcrms(km/s) 34.Gamma 35.fbf 36.fbcf 37.tf(Myr) 37.run\n")
#	bssfile = open(bssfilename, 'r')
#	bssfile.seek(0)
#	bssfile.readline()
#	for line in bssfile:
#		a=line.split()
#		datafilestring=a[-1]
#		dynfile=datafilestring+'/initial.dyn.dat'
#		data = loadtxt(dynfile)
#		convstring=datafilestring+'/initial'
#		units=scripts.read_units(convstring)
#		print units
#		rc = mean(data[-10:-1,7])*units[0]['l_pc']
#		rhoc = mean(data[-10:-1,21])*units[0]['m_msun']/(units[0]['l_pc'])**3.
#		vcrms = mean(data[-10:-1,23])*units[0]['l_cgs']/units[0]['nbt_cgs']/10.**5.
#		Gamma = rhoc**2. * rc**3. / vcrms
#		print 'rc, rhoc, vcrms, Gamma', rc, rhoc, vcrms, Gamma
#		fbi, fbf = table.get_fb_sim(datafilestring, 'initial')
#		fbci, fbcf = table.get_fbc_sim(datafilestring, 'initial')
#		tf = table.get_tf_sim(datafilestring, 'initial')
#
#		for j in range(30):
#			writefile.write("%s " %(a[j],))
#		writefile.write("%f %f %f %f %f %f %f %s\n" %(rc, rhoc, vcrms, Gamma, fbf, fbcf, tf, datafilestring,))
#		
#		#bssfile=datafilestring+'/initial.all_bss_formation_properties_one_sim.dat'
#		#allbssdata=loadtxt(allbssfile)
#		#print allbssdata, len(allbssdata)
#		#try:
#		#	oldest_all = max(allbssdata[:,12])
#		#except IndexError:
#		#	print 'only one BSS'
#		#	oldest_all = allbssdata[12]
#
#
#		#corebssfile=datafilestring+'/initial.core_bss_formation_properties_one_sim.dat'
#		#corebssdata=loadtxt(corebssfile)
#		#try:
#		#	oldest_core = max(corebssdata[:,12])
#		#except IndexError:
#		#	print 'only one BSSc'
#		#	oldest_core = corebssdata[12]
#		
#		#for i in range(13):
#		#	writefile.write("%s " %(a[i],))
#		#writefile.write("%f " %(oldest_all,))
#		#for i in range(13):
#		#	writefile.write("%s " %(a[i+13]))
#		#writefile.write("%f " %(oldest_core))
#		#writefile.write("%f %f %f %f " %(rc, rhoc, vcrms, Gamma,))
#		#writefile.write("%s\n" %(datafilestring))
#
#	writefile.close()



def extract_cluster_properties_for_classified_bss(bssfilename, writefilename):
	"""Once classify_all_bss_stars_all_sims has finished running, each simulation directory should 
	have a file where the branching ratios are already calculated.  Also the full paths of the simulation 
	files are there.  This is just a conenient routine to extract some cluster global properties and write 
	another file with that data also.  Now only interested in rhoc, rc, vcrms to calculate Gamma.  """
	from glob import glob
	import create_table as table
	writefile = open(writefilename, 'w')
	writefile.write("#1.all_bss 2.total_coll 3.total_ss_coll 4.total_bs_coll 5.total_bb_coll 6.total_se 7.total_se_merger 8.total_se_disruption 9.total_had_binint 10.total_mtb 11.total_mtb_pure 12.total_se_binint 13.total_mtb_binint 14.age_oldest(Myr) 15.age_median(Myr) 16.core_count 17.total_core_coll 18.total_core_ss_coll 19.total_core_bs_coll 20.total_core_bb_coll 21.total_core_se 22.total_core_se_merger 23.total_core_se_disruption 24.total_core_had_binint 25.total_core_mtb 26.total_core_mtb_pure 27.total_core_se_binint 28.total_core_mtb_binint 29.age_oldest_core(Myr) 30.age_median_core(Myr) 31.rc(pc) 32.rhoc(Msun/pc^3) 33.vcrms(km/s) 34.Gamma 35.fbi 36.fbf 37.fbci 38.fbcf 39.tf(Myr) 40.HB_count 41.HB_count_core 42.run\n")
	bssfile = open(bssfilename, 'r')
	bssfile.seek(0)
	bssfile.readline()
	for line in bssfile:
		if line.rfind('#')==-1:
			print 'line:', line
			a=line.split()
			datafilestring=a[-1]
			dynfile=datafilestring+'/initial.dyn.dat'
			data = loadtxt(dynfile)
			convstring=datafilestring+'/initial'
			units=scripts.read_units(convstring)
			print units
			rc = mean(data[-10:-1,7])*units[0]['l_pc']
			rhoc = mean(data[-10:-1,21])*units[0]['m_msun']/(units[0]['l_pc'])**3.
			vcrms = mean(data[-10:-1,23])*units[0]['l_cgs']/units[0]['nbt_cgs']/10.**5.
			Gamma = rhoc**2. * rc**3. / vcrms
			print 'rc, rhoc, vcrms, Gamma', rc, rhoc, vcrms, Gamma
			fbi, fbf = table.get_fb_sim(datafilestring, 'initial')
			fbci, fbcf = table.get_fbc_sim(datafilestring, 'initial')
			tf = table.get_tf_sim(datafilestring, 'initial')
	
			snapstring = datafilestring+'/initial.snap*.dat.gz' 
			snaps = glob(snapstring)
			snapno = snaps[-1].split('/')[-1].split('snap')[1].split('.')[0]
			startype_sing_count, startype_bin_count,t_myr = scripts.find_startype(convstring,snapno, 4 , 'HB')
			startype_count = startype_sing_count + startype_bin_count
			startype_sing_count_core, startype_bin_count_core,t_myr = scripts.find_startype_core(convstring,snapno, rc/units[0]['l_pc'],4 ,'HB')
			startype_count_core = startype_sing_count_core + startype_bin_count_core
		

			for j in range(30):
				writefile.write("%s " %(a[j],))
			writefile.write("%f %f %f %f %f %f %f %f %f %d %d %s\n" %(rc, rhoc, vcrms, Gamma, fbi, fbf, fbci, fbcf, tf, startype_count, startype_count_core, datafilestring))
		
		#bssfile=datafilestring+'/initial.all_bss_formation_properties_one_sim.dat'
		#allbssdata=loadtxt(allbssfile)
		#print allbssdata, len(allbssdata)
		#try:
		#	oldest_all = max(allbssdata[:,12])
		#except IndexError:
		#	print 'only one BSS'
		#	oldest_all = allbssdata[12]


		#corebssfile=datafilestring+'/initial.core_bss_formation_properties_one_sim.dat'
		#corebssdata=loadtxt(corebssfile)
		#try:
		#	oldest_core = max(corebssdata[:,12])
		#except IndexError:
		#	print 'only one BSSc'
		#	oldest_core = corebssdata[12]
		
		#for i in range(13):
		#	writefile.write("%s " %(a[i],))
		#writefile.write("%f " %(oldest_all,))
		#for i in range(13):
		#	writefile.write("%s " %(a[i+13]))
		#writefile.write("%f " %(oldest_core))
		#writefile.write("%f %f %f %f " %(rc, rhoc, vcrms, Gamma,))
		#writefile.write("%s\n" %(datafilestring))

	writefile.close()






def branching_ratio(history, binary):
	all=[]
	binint=[]
	pure_binint=[]
	coll=[]
	pure_coll=[]
	merger=[]
	pure_merger=[]
	binint_coll=[]
	binint_merger=[]
	pure_mtb=[]
	
	for i in history.keys():
		all+=[history[i]['position']]
		if binary==1 and len(history[i]['binint']['binint'].keys())>0:
		#if len(history[i]['binint']['binint'].keys())>0:
			binint+=[history[i]['position']]
		
		if binary==1 and len(history[i]['binint']['binint'].keys())>0 and len(history[i]['coll'].keys())==0 and len(history[i]['se'].keys())==0:
			pure_binint+=[history[i]['position']]

		if len(history[i]['coll'].keys())>0:
			coll+=[history[i]['position']]

		if binary==1 and len(history[i]['binint']['binint'].keys())==0 and len(history[i]['coll'].keys())>0 and len(history[i]['se'].keys())==0:
			pure_coll+=[history[i]['position']]
			
		if binary==1 and len(history[i]['se'].keys())>0:
			merger+=[history[i]['position']]

		if binary==1 and len(history[i]['binint']['binint'].keys())==0 and len(history[i]['coll'].keys())==0 and len(history[i]['se'].keys())>0:
			pure_merger+=[history[i]['position']]

		if binary==1 and len(history[i]['coll'].keys())>0 and len(history[i]['binint']['binint'].keys())>0:
			binint_coll+=[history[i]['position']]

		if binary==1 and len(history[i]['se'].keys())>0 and len(history[i]['binint']['binint'].keys())>0:
			binint_merger+=[history[i]['position']]
		if binary==1 and len(history[i]['se'].keys())==0 and len(history[i]['coll'].keys())==0 and len(history[i]['binint']['binint'].keys())==0:
			pure_mtb+=[history[i]['position']]

	if binary==0:
		pure_coll = coll

	branching_ratio={'binint': float(len(binint))/float(len(history.keys())),
			'pure_binint': float(len(pure_binint))/float(len(history.keys())),
			'coll': float(len(coll))/float(len(history.keys())),
			'pure_coll': float(len(pure_coll))/float(len(history.keys())),
			'merger': float(len(merger))/float(len(history.keys())),
			'pure_merger': float(len(pure_merger))/float(len(history.keys())),
			'pure_mtb': float(len(pure_mtb))/float(len(history.keys()))
			}

	return branching_ratio



def branching_ratio_plot(history,r_reference):
	all=[]
	binint=[]
	pure_binint=[]
	coll=[]
	pure_coll=[]
	merger=[]
	pure_merger=[]
	binint_coll=[]
	binint_merger=[]
	
	for i in history.keys():
		all+=[history[i]['position']]
		if len(history[i]['binint']['binint'].keys())>0:
			binint+=[history[i]['position']]
		
		if len(history[i]['binint']['binint'].keys())>0 and len(history[i]['coll'].keys())==0 and len(history[i]['se'].keys())==0:
			pure_binint+=[history[i]['position']]

		if len(history[i]['coll'].keys())>0:
			coll+=[history[i]['position']]

		if len(history[i]['binint']['binint'].keys())==0 and len(history[i]['coll'].keys())>0 and len(history[i]['se'].keys())==0:
			pure_coll+=[history[i]['position']]
			
		if len(history[i]['se'].keys())>0:
			merger+=[history[i]['position']]

		if len(history[i]['binint']['binint'].keys())==0 and len(history[i]['coll'].keys())==0 and len(history[i]['se'].keys())>0:
			pure_merger+=[history[i]['position']]

		if len(history[i]['coll'].keys())>0 and len(history[i]['binint']['binint'].keys())>0:
			binint_coll+=[history[i]['position']]

		if len(history[i]['se'].keys())>0 and len(history[i]['binint']['binint'].keys())>0:
			binint_merger+=[history[i]['position']]

	branching_ratio={'binint': float(len(binint))/float(len(history.keys())),
			'pure_binint': float(len(pure_binint))/float(len(history.keys())),
			'coll': float(len(coll))/float(len(history.keys())),
			'pure_coll': float(len(pure_coll))/float(len(history.keys())),
			'merger': float(len(merger))/float(len(history.keys())),
			'pure_merger': float(len(pure_merger))/float(len(history.keys()))
			}
	branching={'all': all,
		'binint': binint,
		'pure_binint': pure_binint,
		'coll': coll,
		'pure_coll': pure_coll,
		'merger': merger,
		'pure_merger': pure_merger,
		'branching_ratio': branching_ratio,
		'r_dist': r_reference
		}


	#now make plot	
	nobins=20
	min=0
	max=5

	r_dist_hist, lower_edges = histogram(branching['r_dist'],bins=nobins,range=(min,max),normed=False)
	binwidth = lower_edges[-1]/len(lower_edges)
	lower_edges=lower_edges[:(len(lower_edges)-1)]
	r_dist_cumhist = r_dist_hist.cumsum()

	all_hist, lower_edges = histogram(branching['all'],bins=nobins,range=(min,max),normed=False)
	lower_edges=lower_edges[:(len(lower_edges)-1)]

	norm_all_hist = zeros(nobins)
	for i in range(nobins):
		norm_all_hist[i] = float(all_hist[i]) / float(r_dist_hist[i])

	all_cumhist = all_hist.cumsum()
	norm_all_cum = zeros(nobins)
	for i in range(nobins):
		norm_all_cum[i] = float(all_cumhist[i]) / float(r_dist_cumhist[i])
	
	binint_hist, lower_edges = histogram(branching['binint'],bins=nobins,range=(min,max),normed=False)
	lower_edges=lower_edges[:(len(lower_edges)-1)]

	norm_binint_hist = zeros(nobins)
	for i in range(nobins):
		norm_binint_hist[i] = float(binint_hist[i]) / float(r_dist_hist[i])

	binint_cumhist = binint_hist.cumsum()
	norm_binint_cum = zeros(nobins)
	for i in range(nobins):
		norm_binint_cum[i] = float(binint_cumhist[i]) / float(r_dist_cumhist[i])
	
	pure_binint_hist, lower_edges = histogram(branching['pure_binint'],bins=nobins,range=(min,max),normed=False)
	lower_edges=lower_edges[:(len(lower_edges)-1)]
	
	norm_pure_binint_hist=zeros(nobins)
	for i in range(nobins):
		norm_pure_binint_hist[i] = float(pure_binint_hist[i]) / float(r_dist_hist[i])

	pure_binint_cumhist = pure_binint_hist.cumsum()
	norm_pure_binint_cum=zeros(nobins)
	for i in range(nobins):
		norm_pure_binint_cum[i] = float(pure_binint_cumhist[i]) / float(r_dist_cumhist[i])

	coll_hist, lower_edges = histogram(branching['coll'],bins=nobins,range=(min,max),normed=False)
	lower_edges=lower_edges[:(len(lower_edges)-1)]

	norm_coll_hist = zeros(nobins)
	for i in range(nobins):
		norm_coll_hist[i] = float(coll_hist[i]) / float(r_dist_hist[i])

	coll_cumhist = coll_hist.cumsum()
	norm_coll_cum = zeros(nobins)
	for i in range(nobins):
		norm_coll_cum[i] = float(coll_cumhist[i]) / float(r_dist_cumhist[i])
	
	pure_coll_hist, lower_edges = histogram(branching['pure_coll'],bins=nobins,range=(min,max),normed=False)
	lower_edges=lower_edges[:(len(lower_edges)-1)]

	norm_pure_coll_hist = zeros(nobins)
	for i in range(nobins):
		norm_pure_coll_hist[i] = float(pure_coll_hist[i]) / float(r_dist_hist[i])

	pure_coll_cumhist = pure_coll_hist.cumsum()
	norm_pure_coll_cum = zeros(nobins)
	for i in range(nobins):
		norm_pure_coll_cum[i] = float(pure_coll_cumhist[i]) / float(r_dist_cumhist[i])
	
	merger_hist, lower_edges = histogram(branching['merger'],bins=nobins,range=(min,max),normed=False)
	lower_edges=lower_edges[:(len(lower_edges)-1)]

	norm_merger_hist = zeros(nobins)
	for i in range(nobins):
		norm_merger_hist[i] = float(merger_hist[i]) / float(r_dist_hist[i])

	merger_cumhist = merger_hist.cumsum()
	norm_merger_cum = zeros(nobins)
	for i in range(nobins):
		norm_merger_cum[i] = float(merger_cumhist[i]) / float(r_dist_cumhist[i])
	
	pure_merger_hist, lower_edges = histogram(branching['pure_merger'],bins=nobins,range=(min,max),normed=False)
	lower_edges=lower_edges[:(len(lower_edges)-1)]

	norm_pure_merger_hist = zeros(nobins)
	for i in range(nobins):
		norm_pure_merger_hist[i] = float(pure_merger_hist[i]) / float(r_dist_hist[i])

	pure_merger_cumhist = pure_merger_hist.cumsum()
	norm_pure_merger_cum = zeros(nobins)
	for i in range(nobins):
		norm_pure_merger_cum[i] = float(pure_merger_cumhist[i]) / float(r_dist_cumhist[i])


	import gracePlot
	gpl=gracePlot.gracePlot()
	gpl.hold()

	gpl.plot(lower_edges+0.5*binwidth, norm_all_hist, symbols=1, styles=1)
	gpl.plot(lower_edges+0.5*binwidth, norm_binint_hist, symbols=1, styles=1)
	gpl.plot(lower_edges+0.5*binwidth, norm_pure_binint_hist, symbols=1, styles=1)
	gpl.plot(lower_edges+0.5*binwidth, norm_coll_hist, symbols=1, styles=1)
	gpl.plot(lower_edges+0.5*binwidth, norm_pure_coll_hist, symbols=1, styles=1)
	gpl.plot(lower_edges+0.5*binwidth, norm_merger_hist, symbols=1, styles=1)
	gpl.plot(lower_edges+0.5*binwidth, norm_pure_merger_hist, symbols=1, styles=1)
	
	return branching['branching_ratio']


def find_giants(snapfile):
	snapdata=loadtxt(snapfile)
	r_giants=[]
	for i in range(len(snapdata)):
		#single RGBs
		if snapdata[i,7]==0 and snapdata[i,14]==4:
			r_giants.append(snapdata[i,2])
			print len(r_giants)
		elif snapdata[i,7]==1:
			if snapdata[i,17]==4 or snapdata[i,18]==4:
				r_giants.append(snapdata[i,2])
				print len(r_giants)
	
	return r_giants

def find_r_dist(snapfile):
	snapdata=loadtxt(snapfile)
	r=[]
	for i in range(len(snapdata)):
		r.append(snapdata[i,2])

	return r




			
def run_them_all(run):
	classify_all_bss_stars_all_sims('/Volumes/MyBook_raid2/gridruns/RG_Solar_fixed_rvir_more/runs', 0.001, 1.1, 'RG_Solar_fixed_rvir_more_bss_branching_new.dat', 'RG_Solar_fixed_rvir_more_bss_branching_problem_new.dat')
	classify_all_bss_stars_all_sims('/Volumes/MyBook_raid2/gridruns/RG_Solar_fixed_rvir_more1/runs', 0.001, 1.1, 'RG_Solar_fixed_rvir_more1_bss_branching_new.dat', 'RG_Solar_fixed_rvir_more1_bss_branching_problem_new.dat')

		
