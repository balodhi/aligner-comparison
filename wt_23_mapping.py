import time as time
import numpy as np
from time import sleep
import os
import json
import logging
import logging.config



def runAlgorithm(commands,logger):
	numOfRepetition = commands['datainfo']['repetition']
	fastq_file_prefix_list = commands['datainfo']['fastq_file_prefix']
	for i in range(len(fastq_file_prefix_list)):
		fastq_file_prefix = fastq_file_prefix_list[i]
		timedurationoutputfile = fastq_file_prefix +"_time.log"
		logger.info('creating log file: '+timedurationoutputfile)
		elapsed_time = []
		print("Starting prefix:",fastq_file_prefix)
		logger.info('starting prefix: '+fastq_file_prefix)

		'''
		commands = {}
		commands['hisat'] = {}
		commands['hisat']['run_aligner'] = True
		commands['hisat']['build_index'] = False
		commands['hisat']['index_command'] = False
		commands['hisat']['aligner'] = "hisat -p 18 -x hisat1 -1 ./Fastq_files/"+ fastq_file_prefix +"_R1.fastq -2 ./Fastq_files/"+ fastq_file_prefix +"_R2.fastq -S "+ fastq_file_prefix +"_hisat.sam"
		
		commands['hisat2'] = {}
		commands['hisat2']['run_aligner'] = True
		commands['hisat2']['build_index'] = False
		commands['hisat2']['index_command'] = False
		commands['hisat2']['aligner'] = "hisat2 -p 18 -x hisat2 -1 ./Fastq_files/"+ fastq_file_prefix +"_R1.fastq -2 ./Fastq_files/"+ fastq_file_prefix +"_R2.fastq -S "+ fastq_file_prefix +"_hisat2.sam"

		commands['tophat'] = {}
		commands['tophat']['run_aligner'] = True
		commands['tophat']['build_index'] = False
		commands['tophat']['index_command'] = False
		commands['tophat']['aligner'] = "tophat --output-dir ./tophat1_"+ fastq_file_prefix +" -p 18 -G genes.gtf Ath_reference ./Fastq_files/"+ fastq_file_prefix +"_R1.fastq ./Fastq_files/"+ fastq_file_prefix +"_R2.fastq"
		
		commands['tophat2'] = {}
		commands['tophat2']['run_aligner'] = True
		commands['tophat2']['build_index'] = False
		commands['tophat2']['index_command'] = False
		commands['tophat2']['aligner'] = "tophat2 --output-dir ./tophat2_"+ fastq_file_prefix +" -p 18 -G genes.gtf Ath_reference ./Fastq_files/"+ fastq_file_prefix +"_R1.fastq ./Fastq_files/"+ fastq_file_prefix +"_R2.fastq"
		
		commands['star'] = {}
		commands['star']['run_aligner'] = True
		commands['star']['build_index'] = False
		commands['star']['index_command'] = False
		commands['star']['aligner'] = "STAR --runThreadN 18 --genomeDir ./star_index/ --sjdbGTFfile ./genes.gtf --readFilesIn ./Fastq_files/"+ fastq_file_prefix +"_R1.fastq ./Fastq_files/"+ fastq_file_prefix +"_R2.fastq --outFileNamePrefix ./"+ fastq_file_prefix +"_STAR_"

		commands['crac'] = {}
		commands['crac']['run_aligner'] = True
		commands['crac']['build_index'] = False
		commands['crac']['index_command'] = False
		commands['crac']['aligner'] = "crac -i crac -k 22 -r ./Fastq_files/"+ fastq_file_prefix +"_R1.fastq ./Fastq_files/"+ fastq_file_prefix +"_R2.fastq --sam "+ fastq_file_prefix +"_crac.sam --summary "+ fastq_file_prefix +"_crac_summary --nb-threads 18"

		commands['gsnap'] = {}
		commands['gsnap']['run_aligner'] = True
		commands['gsnap']['build_index'] = False
		commands['gsnap']['index_command'] = False
		commands['gsnap']['aligner'] = "gsnap -D ./gsnap_index -d gmap --novelsplicing 1 -t 18 ./Fastq_files/"+ fastq_file_prefix +"_R1.fastq ./Fastq_files/"+ fastq_file_prefix +"_R2.fastq > "+ fastq_file_prefix +"_gsnap.sam"

		commands['soapsplice'] = {}
		commands['soapsplice']['run_aligner'] = True
		commands['soapsplice']['build_index'] = False
		commands['soapsplice']['index_command'] = False
		commands['soapsplice']['aligner'] = "/usr/local/bin/SOAPsplice-v1.10/bin/soapsplice -d ./soapsplice_index/soapsplice.fa.index -1 ./Fastq_files/"+ fastq_file_prefix +"_R1.fastq -2 ./Fastq_files/"+ fastq_file_prefix +"_R2.fastq -I 101 -o ./"+ fastq_file_prefix +"_soapsplice.sam -p 18 -f 2"

		commands['rapmap'] = {}
		commands['rapmap']['run_aligner'] = True
		commands['rapmap']['build_index'] = False
		commands['rapmap']['index_command'] = False
		commands['rapmap']['aligner'] = "rapmap quasimap -i ./rapmap_index/rapmap -1 ./Fastq_files/"+ fastq_file_prefix +"_R1.fastq -2 ./Fastq_files/"+ fastq_file_prefix +"_R2.fastq -t 18 -o "+ fastq_file_prefix +"_rapmap.sam"

		commands['rum'] = {}
		commands['rum']['run_aligner'] = True
		commands['rum']['build_index'] = False
		commands['rum']['index_command'] = False
		commands['rum']['aligner'] = "rum_runner align --index-dir ./rum_index/Arabidopsis/ --name "+ fastq_file_prefix +" -o ./RUM --chunks 18 ./Fastq_files/"+ fastq_file_prefix +"_R1.fastq ./Fastq_files/"+ fastq_file_prefix +"_R2.fastq"

		commands['subread'] = {}
		commands['subread']['run_aligner'] = True
		commands['subread']['build_index'] = False
		commands['subread']['index_command'] = False
		commands['subread']['aligner'] = "subjunc -i ./subread_index/subread -r ./Fastq_files/"+ fastq_file_prefix +"_R1.fastq -R ./Fastq_files/"+ fastq_file_prefix +"_R2.fastq -T 18 --allJunctions --SAMoutput -o "+ fastq_file_prefix +"_subread.sam"

		'''

		fid = open(timedurationoutputfile,'w')
		for key in commands:
			
			if commands[key]['command_type'] == 'aligner' and commands[key]['run_aligner'] == True:
				if commands[key]['build_index']:
					logger.info('building index for '+str(key))
					logger.info('command: '+str(commands[key]['index_command']))
					print('building index...')
					build_index(commands[key],str(key))
				print("Start command :",str(key))
				logger.info('formatting the aligner command. replace all variables in command')
				command = formatCommand(commands[key]['aligner_command'],commands[key]['aligner_variables'],commands['datainfo'],fastq_file_prefix)
				##index check
				for i in range(numOfRepetition):
					print(" "*10,"command (",str(key) ,") and repeatition (",str(i),")"," "*10,'\n')
					print(command)
					print("*"*30)
					if commands[key] == 'rum':
						logger.info('deleting the RUM output folder if present')
						os.system(commands['misc']['delRum'])
					logger.info('starting aligner: '+str(key))
					logger.info('command: '+str(command))
					start_time = time.time()
					os.system(command)
					end_time = time.time()
					elapsed_time.append(end_time-start_time)		
					print("Time Elapsed for the command:",elapsed_time[i])
					print("*"*30)
					logger.info('aligner repetition '+str(i) + ' is completed')
					sleep(0.05)
				print("*"*30)
				stdev = np.std(elapsed_time)
				meanT = np.mean(elapsed_time)
				print ("Standard Deviation of execuation:", stdev)
				print ("Mean of execuation:", meanT)
				print("*"*30)
				elapsed_time = []
				fid.write("*"*30+"\n")
				fid.write(command+"\n")
				fid.write(str(stdev)+"\n")
				fid.write(str(meanT)+"\n")
				fid.write("*"*30+"\n")
				logger.info('time file written. aligner: '+str(key))
			#break
		fid.close()

def build_index(aligner,alignerName):

	if not os.path.exists(alignerName+'_index'):
		os.makedirs(alignerName+'_index')
	os.system(aligner['index_command'])
	print('build complete')



def main():
	setup_logging()
	logger = logging.getLogger(__name__)
	logger.info('Starting the program')
	jsonfile = './commands.json'
	organismName = str(raw_input('Provide the organism name: ') or 'Arabidopsis')
	logger.info('organism Name'+organismName+' was provided')
	if not os.path.exists(os.path.join('.',organismName)):
		print('There is no data of ',organismName,'. Exiting the program !!!')
		logger.info('no information about organism. Program Exit')
		exit()

	print('Reading Commands....')
	logger.info('Started reading the json configuration file')
	commands = readCommands(jsonfile)
	
	default_fastq_file_prefix_list = commands['datainfo']['fastq_file_prefix']
	fastq_file_prefix_list = str(raw_input("Enter a list of fastq files prefix,separated by a comma. Hit enter when you finished:\n"))
	
	if  fastq_file_prefix_list:
		fastq_file_prefix_list = fastq_file_prefix_list.replace(" "," ").split(",")
		fastq_file_prefix_list = filter(None,fastq_file_prefix_list)
		logger.info('list of prefixes are given'+str(fastq_file_prefix_list))
	else:
		print('No list of fastq prefix is provided. Selecting the default list. \n', default_fastq_file_prefix_list)
		logger.info('no list was given loading the default list: '+str(default_fastq_file_prefix_list))
	
	
	os.chdir(os.path.join('.',organismName))
	logger.info('change directory')

	runAlgorithm(commands,logger)

def readCommands(jsonfile):
	with open(jsonfile) as file:
		data = json.loads(file.read())
	return data

def formatCommand(command,commandVariables,datainfo,fastq_file_prefix):
	for i in range(0,len(commandVariables)):
		#print(commandVariables[i],'@@@@')
		if commandVariables[i] == 'fastq_file_prefix':
			command = command.replace('fastq_file_prefix',fastq_file_prefix)
		else:
			command = command.replace(commandVariables[i],datainfo[commandVariables[i]])
		#print(command)
	
	return command


def setup_logging(default_path='logging.json',default_level=logging.INFO,env_key='LOG_CFG'):
    """Setup logging configuration

    """
    path = default_path
    value = os.getenv(env_key, None)
    if value:
        path = value
    if os.path.exists(path):
        with open(path, 'rt') as f:
            config = json.load(f)
        logging.config.dictConfig(config)
    else:
        logging.basicConfig(level=default_level)

if __name__ == "__main__":
	main()




