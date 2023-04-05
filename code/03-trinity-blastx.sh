# Scripts/Code for Zostera/Laby Transcriptome Assembly and Annotation Transcriptome Annotation
/usr/local/bin/Trinity-v2.6.5/util/TrinityStats.pl Trinity.fasta

#Output of TrinityStats:
#Total ‘genes’ = 767,295
#Total transcripts: 1,065,804
#Percent GC: 47.84%
#N50: 1075
#Median Contig length: 349
#Average Contig length: 654.71

#download all Zostera marina protein sequences from GenBank (as an AA Fasta file)
#named Zostera_marina_prot_ncbi.fasta
#Make the blastx database using this AA fasta file
makeblastdb –in Zostera_marina_prot_ncbi.fasta –out Zostera_marina_prot_db –parse_seqids –dbtype prot

#Start the blastx of the Trinity output fasta file against this newly made blastdb
nohup blastx -query Trinity.fasta -db Zostera_marina_prot_db -evalue 1e-10 -outfmt 6 -max_target_seqs 1 > Zm_transcriptome_annotation &
#There are errors in the blastx output files: “Error: (1431.1) FASTA-Reader: Warning: FASTA-Reader: Title is #very long: 1509 characters (max is 1000)”. This has to do with too much information in the fasta headers.
cut -d ' ' -f1 Zostera_marina_prot_ncbi.fasta > Relabeled_Zostera_marina_prot_ncbi.fasta
cut -d ' ' -f1 Trinity.fasta > Relabeled_Trinity.fasta
#These scripts removed the extra information in the fasta headers, leaving only the IDs.
#remake the blastdb and restart the blastx
makeblastdb –in Zostera_marina_prot_ncbi.fasta –out Zostera_marina_prot_db –parse_seqids –dbtype prot
#Start the blastx of the Trinity output fasta file against this newly made blastdb
nohup blastx -query Relabeled_Trinity.fasta -db Relabeled_Zostera_marina_prot_db -evalue 1e-10 -outfmt 6 -max_target_seqs 1 > ZmBlastx_transcriptome_output &

#filter Zoster matched contigs into their own fasta file
filter_fasta.py –f Trinity.fasta –o Zostera_contigs.fasta –s ZmBlastx_transcriptome_output .

#filter the non-Zostera matched contigs into their own fasta file
filter_fasta.py –f Trinity.fasta –o nonZostera_contigs.fasta –s ZmBlastx_transcriptome_output –n
#The Trinity.fasta file has 1,065,804 contigs in it. 281,600 contigs are “Zostera” contigs based on significant #identity to the Zostera marina genome. The nonZostera_contigs.fasta file has 784,204 contigs in it. all contigs #are accounted for

#blastx of Zostera and nonZostera contigs against the SwissProt database
diamond makedb --in /data1/blastdb/uniprot_sprot.fasta –d swissprot
nice nohup diamond blastx –d swissprot –q nonZostera_contigs.fasta –o nonZostera_SwissProt_output --sensitive –f 6 –e 1e-10 –k 1 &
#NONZOSTERA contigs have 87,550 hits on the swissprot database (out of 784,204). Next, will repeat this #diamond blastx using the Zostera contigs
-nice nohup diamond blastx –d swissprot –q Zostera_contigs.fasta –o Zostera_SwissProt_output --sensitive –f 6 –e 1e-5 –k 1 &
#ZOSTERA contigs have 270061 hits on the swissprot database (out of 281,600)

#blast non-zostera contigs to Hondaea genome
nice nohup blastx -query nonZostera_contigs.fasta -db Hondaea_prot_db -evalue 1e-10 -outfmt 6 -max_target_seqs 1 > nonZostera_Hondaea_genome_BLAST_annotation
