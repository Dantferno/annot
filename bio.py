#!/bin/python
import subprocess 
import time
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from latex import build_pdf
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
	bp.index.name = 'ORF'


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
	#orfcount =  le nombre d'ORF
	orfcount = subprocess.Popen('grep {0} ./{0}/ORF_{0}.fasta | wc -l '.format(seqname), shell=True, stdout=subprocess.PIPE)
	output, err = orfcount.communicate()
	orfcount = int(output)
	print('Nombre d\'ORF : {}'.format(orfcount))
	#Conserve le nombre de blast pourchaque ORF dans wc[]
	wc = [0]
	for i in range(1, orfcount+1):
		cmd = subprocess.Popen('grep {0}_{1} ./{0}/blastNR_{0} | wc -l'.format(seqname, i), shell=True, stdout=subprocess.PIPE)
		output, err = cmd.communicate()
		wc.append(int(output))
		print('ORF{0}: {1} hit'.format(i, wc[i]))
	#Nbre de blast ==>Dataframe bp 
	bp['Nbr blastp evalue=1e-10'] = pd.Series(wc[1:], index=bp.index)


	
def ORFchoice():
	global orfchoice
	orfchoice = input('Choisi l\'orf étudié : ')
	
	
def parse_blast(seqname, orfchoice):
	#réparti les hits blast dans des fichiers différents selon l'ORF 
	subprocess.call('mkdir ./{0}/ORF{1}'.format(seqname,orfchoice), shell=True) 
	subprocess.call('grep ^{0}_{1} ./{0}/blastNR_{0} > ./{0}/ORF{1}/blastNR_{0}_{1}'.format(seqname, orfchoice), shell=True) 
	
def plotblast(seqname, orfchoice, bp):
##if seqname.contnain('_') 
	blastresultNR = pd.read_csv('./{0}/ORF{1}/blastNR_{0}_{1}'.format(seqname, orfchoice), sep='\t', names=['qseqid', 'qass', 'evalue', 'bitscore', 'length', 'pident'])
	evalue = blastresultNR['evalue']
	bitscore = blastresultNR['bitscore']
	length = blastresultNR['length']
	pident = blastresultNR['pident']
	bp.at[int(orfchoice), 'best evalue'] = evalue[0]
	bp.at[int(orfchoice), 'bitscore'] = bitscore[0]
	bp.at[int(orfchoice), 'length'] = length[0]
	bp.at[int(orfchoice), 'pident'] = pident[0]
	fig = plt.figure(figsize=[10, 10])
	
	markersize = 1
	x = [i for i in range(len(bitscore))]
	
	ax1 = fig.add_subplot(411)
	ax1.plot(evalue)
	ax1.set(ylabel='Evalue', title='Résultats blastp ORF{0} de {1}'.format(orfchoice, seqname))
	plt.yscale('log')
	ax2 = fig.add_subplot(412)
	ax2.scatter(x, bitscore, s=markersize)
	ax2.set(ylabel='bitscore')
	ax3 = fig.add_subplot(413)
	ax3.scatter(x, length, s=markersize)
	ax3.set(ylabel='longueur alignement')
	ax4 = fig.add_subplot(414)
	ax4.scatter(x, pident, s=markersize)
	ax4.set(xlabel='{} résultats blastp contre NR ordonnés par e-value croissante'.format(len(evalue)), ylabel='pourcentage d\'identité')
	plt.savefig('./{0}/ORF{1}/blastNRplot.pdf'.format(seqname, orfchoice), bbox_inches='tight')
	
def parse_orf(seqname, orfchoice):
	with open('./{0}/ORF_{0}.fasta'.format(seqname), 'r') as f, open('./{0}/ORF{1}/ORF{1}.fasta'.format(seqname, orfchoice), 'w+') as g :
		for ORF in SeqIO.parse(f, 'fasta') :
			nameorf = '{0}_{1}'.format(seqname, orfchoice)
			if ORF.id == nameorf:
				SeqIO.write(ORF, g, 'fasta')
				
def blastSW(seqname, orfchoice):
	print('***Blasting ORF traduction contre SW***')
	#BLASTp tout les ORFs contre SW en ligne avec une evalue de 1e-10, max target seq = 5000 
	subprocess.call('blastp -db swissprot -remote -query ./{0}/ORF{1}/ORF{1}.fasta -out ./{0}/ORF{1}/blastSW_{0} -evalue 1e-10 -outfmt \"6 qseqid sseqid evalue bitscore length pident\" -max_target_seqs 5000'.format(seqname, orfchoice), shell=True) 

	
def blastSWcount(seqname, orfchoice, bp):
	cmd = subprocess.Popen('grep {0} ./{0}/ORF{1}/blastSW_{0} | wc -l'.format(seqname, orfchoice), shell=True, stdout=subprocess.PIPE)		
	output, err = cmd.communicate()
	blastswcount = int(output)
	#Nbre de blast ==>Dataframe bp 
	bp.at[int(orfchoice), 'Nbre hit blastp SW'] = blastswcount

	
def tableau(seqname, bp, orfchoice):
	with open('./{0}/ORF{1}/tableau.csv'.format(seqname, orfchoice), 'w+') as f:
		print(bp)
		bp.to_csv(f)
	
	
	
def tableauLATEX(seqname, bp, orfchoice):
	with open('./{0}/ORF{1}/mytable.tex'.format(seqname, orfchoice), 'w+') as f:
		f.write('\\documentclass{slides} \n')
		f.write('\\usepackage{booktabs} \n')
		f.write('\\begin{document} \n')
		f.write(bp.to_latex())
		f.write('\n\\end{document}')
		build_pdf(f, texinputs=[])

    	
getORF(seqname)
coordonneORF(seqname)
sizeORF(seqname, bp)
blastORF(seqname)
orfblastcount(seqname)
ORFchoice()
parse_blast(seqname, orfchoice)
parse_orf(seqname, orfchoice)
plotblast(seqname, orfchoice, bp)
blastSW(seqname, orfchoice)
blastSWcount(seqname, orfchoice, bp)
tableau(seqname, bp, orfchoice)
print("--- %s seconds ---" % (time.time() - start_time))
