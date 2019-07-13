# `data` Subdirectory Structure

Raw data used for [analyses](https://github.com/eimd-2019/project-EWD-transcriptomics/tree/master/analyses). Repository contents and relevant metadata below.

`trinity` output:
- [genes.counts.matrix.txt](https://github.com/eimd-2019/project-EWD-transcriptomics/blob/master/data/genes.counts.matrix.txt): Gene counts. Row names are gene IDs, and columns are gene counts for different samples.

*Z. marina*:
- [Zostera_contigs.fasta](https://gannet.fish.washington.edu/spartina/2019-07-05-eimd/project-EWD-transcriptomics/data/Zostera_contigs.fasta): FASTA file output from `trinity`
- [Zostera_SwissProt_e5_output](https://github.com/eimd-2019/project-EWD-transcriptomics/blob/master/data/Zostera_SwissProt_e5_output): Tab-delimited `blastx` output
- [Zostera-blast-sep.tab](https://github.com/eimd-2019/project-EWD-transcriptomics/blob/master/data/Zostera-blast-sep.tab): Intermediate file with Uniprot information separated by `|` now separated by tabs
- [Zostera-blast-sort.tab](https://github.com/eimd-2019/project-EWD-transcriptomics/blob/master/data/Zostera-blast-sort.tab): Sorted version of above file.
- [Zostera-blast-annot.tab](https://gannet.fish.washington.edu/spartina/2019-07-05-eimd/project-EWD-transcriptomics/data/Zostera-blast-annot.tab):
- [Zostera-blast-annot-geneIDOnly.tab](https://github.com/eimd-2019/project-EWD-transcriptomics/blob/master/data/Zostera-blast-annot-geneIDOnly.tab): Gene IDs isolated from nnotated `blastx` output
- [Zostera-blast-annot-withGeneID.tab](https://gannet.fish.washington.edu/spartina/2019-07-05-eimd/project-EWD-transcriptomics/data/Zostera-blast-annot-withGeneID.tab):
- [Zostera-blast-annot-withGeneID-noIsoforms.tab](https://github.com/eimd-2019/project-EWD-transcriptomics/blob/master/data/Zostera-blast-annot-withGeneID-noIsoforms.tab): Annotated `blastx` output with gene ID information and no isoforms
- [Zostera-blast-annot-withGeneID-noIsoforms-geneCounts.tab](https://github.com/eimd-2019/project-EWD-transcriptomics/blob/master/data/Zostera-blast-annot-withGeneID-noIsoforms-geneCounts.tab): Annotated `blastx` output with gene IDs, no isoforms, and gene count information

non-*Zostera spp*:
- [nonZostera_contigs.fasta](https://gannet.fish.washington.edu/spartina/2019-07-05-eimd/project-EWD-transcriptomics/data/nonZostera_contigs.fasta): FASTA file output from `trinity`
- [nonZostera_SwissProt_e5_outputBOX.txt](https://github.com/eimd-2019/project-EWD-transcriptomics/blob/master/data/nonZostera_SwissProt_e5_outputBOX.txt): Tab-delimited `blastx` output
- [nonZostera-blast-sep.tab](https://github.com/eimd-2019/project-EWD-transcriptomics/blob/master/data/nonZostera-blast-sep.tab): Intermediate file with Uniprot information separated by `|` now separated by tabs
- [nonZostera-blast-sort.tab](https://github.com/eimd-2019/project-EWD-transcriptomics/blob/master/data/nonZostera-blast-sort.tab): Sorted version of above file.
- [nonZostera-blast-annot.tab](https://gannet.fish.washington.edu/spartina/2019-07-05-eimd/project-EWD-transcriptomics/data/nonZostera-blast-annot.tab):
- [nonZostera-blast-annot-geneIDOnly.tab](https://github.com/eimd-2019/project-EWD-transcriptomics/blob/master/data/nonZostera-blast-annot-geneIDOnly.tab): Gene IDs isolated from nnotated `blastx` output
- [nonZostera-blast-annot-withGeneID.tab](https://github.com/eimd-2019/project-EWD-transcriptomics/blob/master/data/nonZostera-blast-annot-withGeneID.tab):
- [nonZostera-blast-annot-withGeneID-noIsoforms.tab](https://github.com/eimd-2019/project-EWD-transcriptomics/blob/master/data/nonZostera-blast-annot-withGeneID-noIsoforms.tab): Annotated `blastx` output with gene ID information and no isoforms
- [nonZostera-blast-annot-withGeneID-noIsoforms-geneCounts.tab](https://github.com/eimd-2019/project-EWD-transcriptomics/blob/master/data/nonZostera-blast-annot-withGeneID-noIsoforms-geneCounts.tab): Annotated `blastx` output with gene IDs, no isoforms, and gene count information
