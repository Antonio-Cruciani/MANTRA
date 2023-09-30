import os
import subprocess
from subprocess import Popen, PIPE, STDOUT
import json




	
def exact_meausures(inputList,topt):
	for path_opt in topt:
		for dataset in inputList:
			dd = {'ds':dataset,'topt':path_opt}
			with open('exp.json', 'w+') as f:
				json.dump(dd, f)

			cmd =  ["julia", "experiment_exact.jl"]
			p = Popen(cmd, stdout=PIPE, stderr=STDOUT)
			for line in p.stdout:
				print(line)
			print("#----------------------------#")
				

	


def estimations_sampling(inputList,epsilon,sample_list,topt,sampler):
	for path_opt in topt:
		j = 0
		for eps in epsilon:
			for dataset in inputList:
				for samp in sampler:
					dd = {'ds':dataset,'topt':path_opt,"eps":eps,"ss":sample_list[j],'algo':"ob",'sampler':samp}
					with open('exp.json', 'w+') as f:
						json.dump(dd, f)
					cmd =  ["julia", "experiment_bernstein.jl"]
					for i in range(10):
						p = Popen(cmd, stdout=PIPE, stderr=STDOUT)
						for line in p.stdout:
							print(line)
						print("#----------------------------#")
			j+=1


	


'''
inputList = [
	"16_brain_100206_90.txt",
    "17_brain_100206_70.txt",
    "01_hypertext.txt",
    "02_highschool.txt",
    "03_hospital_ward.txt",
    "04_college_msg.txt",
    "05_wiki_elections.txt",
    "06_highschool.txt",
    "07_digg_reply.txt",
    "08_infectious.txt",
    "09_primary_school.txt",
    "10_facebook_wall.txt",
    "11_slashdot_reply.txt",
    "12_highschool.txt",
    "13_topology.txt",
    "14_SMS.txt",
    "21_mathoverflow.txt",
    "20_askubuntu.txt",
    "22_superuser.txt"

]

inputList = [
    "20_askubuntu.txt",
    "22_superuser.txt"

]
'''
inputList = [
    "21_mathoverflow.txt"


]

epsilon = [0.1,0.07,0.05,0.01]
sample_list = [100,350,750,1000]
#topt = ["pfm","sh","sfm"]
topt = ["sh","sfm"]
sampler = ["bernstein","wub","cmcera"]

#estimations_sampling(inputList,epsilon,sample_list,topt,sampler)
exact_meausures(inputList,topt)