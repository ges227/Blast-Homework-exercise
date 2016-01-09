#import sys
import os.path
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
blast_prog= "/Users/GaelleSop/Desktop/Bioinformatics_Computing_BCHB_524/Bioinformatics_Computing/ncbi-blast-2.2.31+/bin/blastp"
blast_query= "./ex20blastdb/drosoph-ribosome.fasta"
blast_db= "./ex20blastdb/yeast-ribosome.fasta"
if not os.path.exists("ex20results.xml"):
	cmdline= NcbiblastpCommandline(cmd=blast_prog, query=blast_query, db=blast_db,outfmt=5,out="ex20results.xml")
	stdout,stderr=cmdline()


result_handle= open("ex20results.xml")
blastresultdict={}
for blast_result in NCBIXML.parse(result_handle):
	print "\n\n","DROSOPHILA QUERY:",blast_result.query
	querynum_start=blast_result.query.find("gi|")
	querynum_fin= blast_result.query.find("|ref")
	querynum=blast_result.query[querynum_start+3:querynum_fin]
	print querynum
	num1hsp=True
	hsp=blast_result.alignments[0].hsps[0]
	if hsp.expect<1e-5:
		blastresultdict[querynum]=(hsp.expect,hsp.score,blast_result.alignments[0].title)
		print "YEAST BEST MATCH:",blast_result.alignments[0].title, 
		print "Alignment best e-value:",hsp.expect
		print "Alignment best score", hsp.score
		print "Alignment best positives", hsp.positives
		print "Alignment best gaps", hsp.gaps
	if querynum not in blastresultdict:
		print "No good yeast match found for this protein query."
result_handle.close()


#sortedblastresults_eval= sorted(blastresultdict.items(), key = lambda p: p[1][0])
blastresults_sortedbyscore= sorted(blastresultdict.items(), key = lambda p: p[1][1], reverse=True)
#for item in blastresults_sortedbyscore:
#	print item
print "\n\n***********\nGI ID of DROSOPHILA protein that is most highly conserved in both databases:", blastresults_sortedbyscore[0][0]
print "Corresponding Yeast Match:\n", blastresultdict[blastresults_sortedbyscore[0][0]][2]
#print "Corresponding evalue is",blastresults_sortedbyscore[0][1][0] 
#print "Corresponding Score is",blastresults_sortedbyscore[0][1][1] 
