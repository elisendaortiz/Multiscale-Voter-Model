#!/usr/bin/python
# -*- coding: utf-8 -*-
import time
from operator import itemgetter
import numpy as np
import random as rnd
import sys, traceback
import os

def load_network(net):
	global G,N,neighbors,states,node_labels,t,results,consensus1,consensus0,consensus_dict,coord,numgroups,groups,degree,states_groups,node_group,cumm_degree_group
	## Read network
	neighbors = {}
	for line in np.loadtxt(path_network_files+net+'.edge'):
		i = int(line[0])
		j = int(line[1])		
		if i not in neighbors:
			neighbors[i] = []
		if j not in neighbors:
			neighbors[j] = []
		if j not in neighbors[i]:
			neighbors[i].append(j)
		else:
			continue
#			print "Error!", i, j #This errors rise when the edgelist has double connections eg. 1 4 and afterwards it appears 4 1.
		if i not in neighbors[j]:
			neighbors[j].append(i)
		else:
			continue
#			print "Error!", i, j #This errors rise when the edgelist has double connections eg. 1 4 and afterwards it appears 4 1.
	degree = {}
	for i in sorted(neighbors):
		degree[i] = len(neighbors[i])
	N=len(degree)
	node_labels=list(degree.keys())
	# Read embedding coordinates (r,theta) 
	coord = {}
   	with open(path_network_files+net+'.inf_coord') as f:  
		for line in f.readlines():
		    li = line.lstrip()
		    if not li.startswith("#"):
			key= int(line.split()[0])
			coord[key]=[float(line.split()[3]),float(line.split()[2])] #{nodelabel:[r,theta]}


def make_groups():
	global G,N,neighbors,states,node_labels,t,results,consensus1,consensus0,consensus_dict,coord,numgroups,groups,degree,states_groups,node_group,cumm_degree_group
	#To make the groups of adjoining nodes in S1 we sort nodes by their theta's so to ensure they are adjoining, then and allocate these into groups of certain size=num_elements_per_group
	#Notice we are "rounding" down num_elements_per_group, so the last group will always contain the extra nodes whenever the division isn't exact
	#We create two objects. 
	# 1st) groups={group_1:[node,..node_whatever], ..., group_last:[nodes]}, to use when calculating the overall state of the group based on the nodes it contains
	# 2nd) node_group ={nodelabel: group_1, ...etc} to quick access the group to which each node belongs when calculating the probability of copying states
	coord_list=[]
	for key in sorted(coord):
		coord_list.append((key,coord[key][0],coord[key][1]))
	coord_list=sorted(coord_list,key=itemgetter(2))
	groups={} 
	for i in range(1,numgroups+1):
		groups[i]=[]

	num_elements_per_group=int(N/numgroups) 
	igroup=1
	node_group=dict.fromkeys(node_labels) 
	for index in range(len(coord_list)):
		node=coord_list[index][0]
		groups[igroup].append(node)
		node_group[node]=igroup	
		if index == int((igroup*num_elements_per_group)-1):
			if igroup >= numgroups:
				igroup = numgroups
			else:
				igroup= igroup+1


def groups_states():
	global G,N,neighbors,states,node_labels,t,results,consensus1,consensus0,consensus_dict,coord,numgroups,groups,degree,states_groups,node_group,cumm_degree_group
	states_groups={} # {group_1: state_of_group,..., etc }
	cumm_degree_group={}  # {group_1: cummulative_degree_of_group,..., etc }
	for igroup in sorted(groups):
		cum_num=0.0
		cum_den=0.0
		for node in groups[igroup]:
			cum_num=cum_num+(states[node]*degree[node])
			cum_den=cum_den+degree[node]
		degree_weigthed_state=cum_num/cum_den
		states_groups[igroup]=degree_weigthed_state
		cumm_degree_group[igroup]=cum_den#/float(len(groups[igroup]))#num nodes in group we do it this way coz the last group may vary
		

def init(params):
	global G,N,neighbors,states,node_labels,t,results,consensus1,consensus0,consensus_dict,coord,numgroups,groups,degree,states_groups,node_group,cumm_degree_group
	#Create a dictionary of states of neighbors and initialize a fraction of them in the ON state (=1), while leaving the rest OFF (=0)
	num_nodes_state1=int(round(float(params[1])*N))
	if num_nodes_state1 ==0:
		print 'In order to initiate the dynamics at least 1 node needs to be in state ON. Please satisfy:\n ini_frac_state1 >=', 1.0/float(N),'\n'
		sys.exit()	
	states=dict.fromkeys(node_labels, int(0))
	rnd.shuffle(node_labels)
	for i in range(1,num_nodes_state1+1):
		states[node_labels[i-1]]=int(1)
	#Initialize time:
	t=0.0
	

def voter_step():
	global G,N,neighbors,states,node_labels,t,results,consensus1,consensus0,consensus_dict,coord,numgroups,groups,degree,states_groups,node_group,cumm_degree_group
	#Select a random node 'i' in the network and a random neighbor 'j' of 'i'
	irand= int(round(rnd.random()*(len(node_labels)-1))) # random_scaled_value = min + (value * (max - min)) where range=[0,len(nodelabels)-1]. El -1 es perq python compta des del 0.
	node_i=node_labels[irand]
	node_i_neighbors= neighbors[node_i]
	jrand= int(round(rnd.random()*(len(node_i_neighbors)-1)))
	node_j=node_i_neighbors[jrand]

	#Update state of node 'i' to that of node 'j' with probability: Pi-->j = 1- | state_j - state_of_j's_group |
	j_group= node_group[node_j]	
	Pj=1.0-abs(float(states[node_j])-float(states_groups[j_group]))
	if Pj <0.0:
		print 'Warning: Negative probability of copying neighbor Pj'
		print node_j, states[node_j],states_groups[j_group]
		
		
	if rnd.random() <= Pj: #Condition to copy state of node 'j'
		if states[node_i] != states[node_j]: # IF 'i' CHANGES STATE by COPYING 'j', then we update states of node 'i' and its group, otherwise they stay the same eventhough we've copied 'j'
		#	#Update state of node 'i'
			states[node_i]=states[node_j]
		#	#Update state of group of node 
			i_group=node_group[node_i]
			if int(states[node_i]) == 1: #When node 'i' new state is 1.
				states_groups[i_group]=states_groups[i_group] + (float(degree[node_i])/float(cumm_degree_group[i_group]))
				if states_groups[i_group] >1.000 and (states_groups[i_group] -1.00) <= 0.000000001 :
					states_groups[i_group]=1.00						
			elif int(states[node_i]) == 0:#When node 'i' new state is 0.
				states_groups[i_group]=states_groups[i_group] - (float(degree[node_i])/float(cumm_degree_group[i_group]))
				if states_groups[i_group] <= 0.000000001:
					states_groups[i_group] =0.00

	#Increase time by 1/N so that N iterations are equiv to 1 time step and for each time step all nodes on average have had the opportunity to switch states.
	t=t+float(1.0/float(N))

def observe():
	global G,N,neighbors,states,node_labels,t,results,consensus1,consensus0,consensus_dict,coord,numgroups,groups,degree,states_groups,node_group,cumm_degree_group
	Nstate1=sum(value == 1 for value in states.values())
	Nstate0=N-Nstate1	
	results.append((t,Nstate1,Nstate0))
	if Nstate1 == N: 
		consensus1= True
	else:
		consensus1 = False
	if Nstate0 ==N:
		consensus0= True
	else:
		consensus0=False	


def stats_1run(irun,itera):
	global G,N,neighbors,states,node_labels,t,results,consensus1,consensus0,consensus_dict,coord,numgroups,groups,degree,states_groups,node_group,cumm_degree_group
	#Gives dictionary with {irun: (type_of_consensus, maxtime,maxiter)} 
	#where: type_of_consensus =1 if all nodes ON (and 0 if all OFF), maxtime = time to reach consensus, maxiter= nº of iterations used
	last=len(results)-1
	if consensus1 == True:
		consensus_dict[irun]=(1,results[last][0],itera)
	else:
		consensus_dict[irun]=(0,results[last][0],itera)


def stats_allruns(params):
	global G,N,neighbors,states,node_labels,t,results,consensus1,consensus0,consensus_dict,coord,numgroups,groups,degree,states_groups,node_group,cumm_degree_group
	Ncons1=0.0
	Ncons0=0.0
	T1=[]
	T0=[]
	for item in consensus_dict.values():		# iterate over list of tuples (type_of_consensus, maxtime, maxiter). Each tuple is of one run=key.
		if int(item[2]) !=int(N*params[3]): 	#Only when consensus has been really reached so (niter NEQ maxiter)
			if item[0] == 1:
				Ncons1=Ncons1+1
				T1.append(float(item[1]))
			if item[0] == 0:
				Ncons0=Ncons0+1
				T0.append(float(item[1]))
	cons1_prob=float(Ncons1)/float(params[2])
	cons0_prob=float(Ncons0)/float(params[2])
	if Ncons1 >0:
		av_time1=np.mean(np.array(T1))
		stdev_time1=np.std(np.array(T1))
	else:
		#print 'No case of consensus of state= ON'
		av_time1=str('Unknown') #In fact 'Unknown' means the av.time would be infinity since it was never reached this sort of consensus
		stdev_time1=str('Unknown')
	if Ncons0 >0:
		av_time0=np.mean(np.array(T0))
		stdev_time0=np.std(np.array(T0))
	else:
		#print 'No case of consensus of state= OFF'
		av_time0=str('Unknown')
		stdev_time0=str('Unknown')

	if save_output_files==True:
			#       0            1        2          3           4        5         6	7
		print>>outf1,groupsize, params[1],av_time1,stdev_time1,av_time0,stdev_time0,cons1_prob,cons0_prob


#~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
# MAIN CODE
#~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
start_time = time.time()

net = sys.argv[1]
groupsize = int(sys.argv[2]) #You can provide this or numgroups and modify the 2 lines after print '-----',net, '-----' to run over numgroups.

#PARAMETERS AND VARIABLES:++++++++++++++++++++++++++++++++
ini_frac_state1=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
# We control how many times "all nodes on average" have changed state, that is how many times we do N iterations, which means maxiter fixes the number of time steps of the dynamics.
# This is so because time advances as delta_time=1/N, so every N iterations we have 1 time_step. At every completed time_step, all nodes on average, have had the opportunity to change state
# because there has been N iterations in this time_step. This means total_num_iterations = N*maxiter
maxiter=10.0**4 	#10.0**4 
maxruns=1000 		#Per nar be hauria de ser in (10**3 - 10**5) | This is how many realizations we do for every initial condition for statistical purposes.

#INPUTS:
path_network_files='./'

#OUTPUTS:
save_output_files=True 
outf1 = open('./'+net+'_output_r='+str(int(groupsize))+'.txt', "w")
print>>outf1, '# groupsize (in nodes)| rho_ini_on| <T>on | std_<T>on | <T>off | std_<T>off | Prob_consensus_on| Prob_consensus_off | '
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#-----	CODE	--------------
load_network(net)
print '-----',net, '-----'
numgroups=int(N/groupsize)
#Divide network into groups
if numgroups > N:
	print 'WARNING: Number of groups cannot exceed number of nodes in the network, choose a smaller value of numgroups parameter.'
	sys.exit()
make_groups()

#Loop over initial conditions--------------------------
icond=1
for condini in ini_frac_state1:
	print 'r=',groupsize,'ini_frac=',condini
	params=[net, condini,maxruns,maxiter]
	rnd.seed(icond)
	# Loop over runs/realizations (1run goes from time=0 to time=tmax) ------------------
	irun=1
	consensus_dict={}
	consensus_times=[]

	while irun <=params[2]:
		# 1 run: 
		init(params)
		groups_states()
		results=[]
		observe()
		itera=1
		#Loop over time iterations (iteration =time-steps until reaching consensus) ---------------
		while consensus1== consensus0 ==False and itera <= float(N*params[3]):
			voter_step()
			if int(itera) % int(N) == 0:
				#We get the values to store in 'results' only at each completed time step (every N iterations)
				observe()
			if consensus1 ==True or consensus0 ==True:
				#print 'iterations needed for consensus=',itera
				break
			elif int(itera)==int(N*params[3]):
				#print 'maxiter reached without consensus'
				break
			itera=itera+1
			
		#End Loop over time iterations -------------------------------------------------------------
		if int(irun) % int(100) == 0:
			print 'realizations=',irun,'| time=',(time.time() - start_time),'sec'

		stats_1run(irun,itera)
		del results,states,states_groups		
		irun=irun+1

	#End Loop over runs------------------------------------------------------------------
	stats_allruns(params)
	del consensus_dict,params 
	icond=icond+1
#End Loop over initial conditions-------------------------	

#-----END CODE--------------	

	


