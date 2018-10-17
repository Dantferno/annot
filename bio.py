#!/bin/python
import subprocess 
import time
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
start_time = time.time()

seqname = subprocess.Popen('grep -oP \'(?<=>)\w+\' seq.fasta', shell=True, stdout=subprocess.PIPE)
output, err = seqname.communicate()
seqname = str(output.decode('ascii').strip('\n'))
subprocess.call('mkdir \'{}\''.format(seqname), shell=True)


	
def getORF(seqname):
	#EMBOSS getorf, minsize 180, find 0 (entre codon stop)
	subprocess.call('getorf seq.fasta ./{0}/ORF_{0}.fasta -minsize 180'.format(seqname), shell=True)
	

def coordonneORF(seqname):
	global bp
	sizes = subprocess.Popen('cut -d\'[\' -f 2 ./{0}/ORF_{0}.fasta | cut -d\']\' -f 1 -s'.format(seqname), shell=True, stdout=subprocess.PIPE)
	output, err = sizes.communicate()
	sizes = str(output.decode('UTF-8').rstrip())
	b = sizes.splitlines()
	bp = pd.DataFrame(b, columns=['Coordonne'])
	bp.index += 1


def sizeORF(seqname, bp):
	orflist = []
	orfaa= []
	with open('./{0}/ORF_{0}.fasta'.format(seqname),'r') as f :
		orfs = SeqIO.parse(f, 'fasta')
		for orf in orfs:
			orflist.append(len(orf.seq))
	bp['Nbr AA'] = pd.Series(orflist, index=bp.index)
	for i in orflist:
		orfaa.append(i*3)
	bp['Nbr nucl'] = pd.Series(orfaa, index=bp.index)

	
def blastORF(seqname):
	print('***Blasting ORF traduction contre NR***')
	#BLASTp tout les ORFs contre nr en ligne avec une evalue de 1e-10, max target seq = 5000 
	subprocess.call('blastp -db nr -remote -query ./{0}/ORF_{0}.fasta -out ./{0}/blastNR_{0} -evalue 1e-10 -outfmt \"6 qseqid sacc evalue bitscore length pident\" -max_target_seqs 5000'.format(seqname), shell=True) 
	
	
def orfblastcount(seqname):
	global orfcount
	orfcount = subprocess.Popen('grep {0} ./{0}/ORF_{0}.fasta | wc -l '.format(seqname), shell=True, stdout=subprocess.PIPE)
	output, err = orfcount.communicate()
	orfcount = int(output)
	print('Nombre d\'ORF : {}'.format(orfcount))
	wc = [0]
	for i in range(1, orfcount+1):
		cmd = subprocess.Popen('grep {0}_{1} ./{0}/blastNR_{0} | wc -l'.format(seqname, i), shell=True, stdout=subprocess.PIPE)
		output, err = cmd.communicate()
		wc.append(int(output))
		print('ORF{0}: {1} hit'.format(i, wc[i]))
	bp['Nbr blastp evalue=1e-10'] = pd.Series(wc[1:], index=bp.index)
	print(bp.to_latex())

	
def ORFchoice():
	global orfchoice
	orfchoice = input('Choisi l\'orf étudié : ')
	
	
def parse_blast(seqname, orfchoice):
	#réparti les hits blast dans des fichiers différents selon l'ORF 
	subprocess.call('mkdir ./{0}/ORF{1}'.format(seqname,orfchoice), shell=True) 
	subprocess.call('grep ^{0}_{1} ./{0}/blastNR_{0} > ./{0}/ORF{1}/blastNR_{0}_{1}'.format(seqname, orfchoice), shell=True) 
	
def plotblast(seqname, orfchoice):
	blastresultNR = pd.read_csv('./{0}/ORF{1}/blastNR_{0}_{1}'.format(seqname, orfchoice), sep='\t', names=['qseqid', 'qass', 'evalue', 'bitscore', 'length', 'pident'])
	evalue = blastresultNR['evalue']
	bitscore = blastresultNR['bitscore']
	length = blastresultNR['length']
	pident = blastresultNR['pident']
	fig = plt.figure(figsize=[10, 10])
	ax1 = fig.add_subplot(411)
	ax1.plot(evalue, label='Evalue',color='y')
	ax1.set(ylabel='Evalue', title='Résultats blastp ORF{0} de {1}'.format(orfchoice, seqname))
	plt.yscale('log')
	ax2 = fig.add_subplot(412)
	ax2.plot(bitscore)
	ax2.set(ylabel='bitscore')
	ax3 = fig.add_subplot(413)
	ax3.plot(length)
	ax3.set(ylabel='longueur alignement')
	ax4 = fig.add_subplot(414)
	ax4.plot(pident)
	ax4.set(xlabel='{} résultats blastp contre NR ordonnés par e-value croissante'.format(len(evalue)), ylabel='pourcentage d\'identité')
	plt.savefig('./{0}/ORF{1}/blastNRplot.pdf'.format(seqname, orfchoice), bbox_inches='tight')
	
def parse_orf(seqname, orfchoice):
	with open('./{0}/ORF_{0}.fasta'.format(seqname), 'r') as f, open('./{0}/ORF{1}/ORF{1}.fasta'.format(seqname, orfchoice), 'w+') as g :
		for ORF in SeqIO.parse(f, 'fasta') :
			nameorf = '{0}_{1}'.format(seqname, orfchoice)
			if ORF.id == nameorf:
				SeqIO.write(ORF, g, 'fasta')
				
def blastSW(seqname, orfchoice):
	print('***Blasting ORF traduction contre NR***')
	#BLASTp tout les ORFs contre SW en ligne avec une evalue de 1e-10, max target seq = 5000 
	subprocess.call('blastp -db swissprot -remote -query ./{0}/ORF{1}/ORF{1}.fasta -out ./{0}/ORF{1}/blastNR_{0} -evalue 1e-10 -outfmt \"6 qseqid sseqid evalue bitscore length pident\" -max_target_seqs 5000'.format(seqname, orfchoice), shell=True) 



# getORF(seqname)
coordonneORF(seqname)
sizeORF(seqname, bp)
# blastORF(seqname)
orfblastcount(seqname)
# ORFchoice()
# parse_blast(seqname, orfchoice)
# plotblast(seqname, orfchoice)
# swissprotblast(seqname, orfchoice)
print("--- %s seconds ---" % (time.time() - start_time))
#inteproscan
