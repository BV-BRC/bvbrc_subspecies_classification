Release Notes
----------------
Small changes in the code to tolerate Blast crashes when searching matches for very short amino acid sequences.

VIGOR3 Dependencies
-------------------
VIGOR3 requires the following Perl modules not included in the current standard core set:
Bio::SeqIO (BioPerl)
DBI
DBD::SQLite
LWP::Simple (used only by adhoc/getNucGB.pl)
LWP::UserAgent (used only by adhoc/getNucGB.pl)

VIGOR3 is also dependent upon the following programs being installed on the same computer:

SQLite database manager
NCBI Blast Package v.2.2.26 (VIGOR3 is currently incompatible with newer Blast versions that use different command-line options.)
CD-HIT (http://weizhongli-lab.org/cd-hit/)
ClustalW2 (http://www.clustal.org/clustal2/)
Muscle (http://www.drive5.com/muscle/)

Deploying VIGOR3
----------------

1. un-tar VIGOR3.tgz.  This creates a directory named VIGOR3 containing the VIGOR
   software and reference databases.

	$ tar xzvf VIGOR3.tgz -C /mypath

2. define a scratch space for vigor

	$ # path used here is an example, any directory will do 
	$ mkdir /mypath/VIGOR3/tempspace
	$ chmod 777 /mypath/VIGOR3/tempspace

3. define a symbolic link for the scratch space (You could also edit createSymLinks.csh 
   to suit your environment and create the symbolic links with this script - save the 
   edited version somewhere else, so you can re-use it directly on any following install).

	$ # symbolic link requires FULL path
	$ cd /mypath/VIGOR3
	$ chmod 775 prod3
	$ ln -s /mypath/VIGOR3/tempspace prod3/vigorscratch

4. define symbolic links for external programs

	$ cd /mypath/VIGOR3
	$ ln -s /usr/local/bin/perl     prod3/perl
	$ ln -s /usr/local/bin/blastall prod3/blastall
	$ ln -s /usr/local/bin/bl2seq   prod3/bl2seq
	$ ln -s /usr/local/bin/formatdb prod3/formatdb
	$ ln -s /usr/local/bin/fastacmd prod3/fastacmd
	$ ln -s /usr/local/bin/clustalw prod3/clustalw2
	$ ln -s /usr/local/bin/muscle   prod3/muscle
	$ chmod 555 prod3
	
	notes:

	1. muscle is only used by the programs in "dbutils", it is not required
	   by VIGOR.
	2. the dbutils directory under prod3 contains utility programs used to support
	   the creation of reference databases for VIGOR
	3. the adhoc directory under prod3 contains a handful of adhoc programs
	   created during the project, these programs use many of VIGOR's library
	   functions but are not part of VIGOR
	4. three additional pipeline programs are contained in the prod3 directory
		a. rna_finder - used by the JCVI pipeline to annotate non-coding genes
		b. tblUTR - used by the JCVI pipeline to extend gene boundsries to
		   include the UTRs
		c. hmm3Evidence - used by the JCVI pipeline to suppply HMM3 evidence
		   supporting the functional annotation of the gene

 
Running VIGOR3
--------------

Example:

	$ /mypath/prod3/VIGOR3.pl -d yfv -i samples/westnile.fasta -o test/westnile
	(sample fasta and output files can be found in the samples directory)

Usage:
  -- allow VIGOR to choose the reference database
  $./VIGOR3.pl -i inputfasta -o outputprefix

  -- tell VIGOR which reference database to use
  $./VIGOR3.pl -d refdb -i inputfasta -o outputprefix

Command Line Options
  -a auto-select the reference database, equivalent to "-d any", default behavior unless
      overridden by -d or -G, (-A is a synonym for this option)
  -d <ref db>, specify the reference database to be used, (-D is a synonym for this option)
  -e <evalue>, override the default evalue used to identify potential genes, the default
     is usually 1E-5, but varies by reference database
  -c <pct ref> minimum coverage of reference product (0-100) required to report a gene, by
     default coverage is ignored
  -C complete (linear) genome (do not treat edges as gaps)
  -0 (zero) complete circular genome (allows gene to span origin)
  -f <0, 1, or 2>, frameshift sensitivity, 0=ignore frameshifts, 1=normal (default), 2=sensitive
  -G <genbank file>, use a genbank file as the reference database, caution: VIGOR genbank
     parsing is fairly rudimentary and many genbank files are unparseable.  Partial genes will
     be ignored. Note: genbank files do not record enough information to handle RNA editing
  -i <input fasta>, path to fasta file of genomic sequences to be annotated, (-I is a synonym
      for this option)
  -l do NOT use locus_tags in TBL file output (incompatible with -L)
  -L USE locus_tags in TBL file output (incompatible with -l)
  -o <output prefix>, prefix for outputfile files, e.g. if the ouput prefix is /mydir/anno
     VIGOR will create output files /mydir/anno.tbl, /mydir/anno.stats, etc., (-O is a synonym
     for this option)
  -P <parameter=value~~...~~paramaeter=value>, override default values of VIGOR parameters
  -j turn off JCVI rules, JCVI rules treat gaps and ambiguity codes conservatively, use
     this option to relax these constraints and produce a more speculative annotation
  -m ignore reference match requirements (coverage/identity/similarity), sometimes useful
     when running VIGOR to evaluate raw contigs and rough draft sequences
  -s <gene size> minimum size (aa) of product required to report a gene, by default size is
     ignored

Outputs:
  outputprefix.rpt - summary of program results
  outputprefix.stats - run statistics (per genome sequence) in tab-delimited format
  outputprefix.cds - fasta file of predicted CDSs
  outputprefix.pep - fasta file of predicted proteins
  outputprefix.tbl - predicted features in GenBank tbl format
  outputprefix.aln - alignment of predicted protein to reference, and reference protein to genome
  outputprefix.fs - subset of aln report for those genes with potential sequencing issues
  outputprefix.at - potential sequencing issues in tab-delimited format

Reference Databases:
  Name        Description                                   (Synonyms)
  any         any virus                                     (vda)
  cov_abcdx   Alpha/Beta/Gamma/Delta/Unclassified Cov*      
  veev        Alphaviruses (VEEV/EEEV)                      (alpha, eeev)
  bunya       Bunyaviridae                                  
  hanta       Bunyaviridae Hantavirus                       (hantavirus)
  obunya      Bunyaviridae Orthobunyavirus                  
  bunya_misc  Bunyaviridae miscellaneous                    
  cdv         Canine Distemper Virus                        
  gcv         Coronavirus                                   (cov)
  gcv_g1a     Coronavirus Group 1A                          (cov_g1a)
  gcv_g1b     Coronavirus Group 1B                          (cov_g1b)
  gcv_g2a     Coronavirus Group 2A                          (cov_g2a)
  gcv_g2b     Coronavirus Group 2B (SARS)                   (cov_g2b, cov_sars, gcv_sars, sars)
  gcv_g2cd    Coronavirus Group 2C & 2D                     (cov_g2c, cov_g2cd, cov_g2d, gcv_g2c, gcv_g2d)
  gcv_g3      Coronavirus Group 3                           (cov_g3)
  filo        Filoviridae (Ebola/Marburg)                   (ebola, marburg)
  giv         Flu                                           (flu, flumb, fluutr, giv2, giv3, piv, swiv)
  giv_a       Flu A                                         (flu_a)
  giv_b       Flu B                                         (flu_b)
  giv_c       Flu C                                         (flu_c)
  hrv         Human Rhinovirus/Enterovirus                  (entero, hrv2, rhino)
  hadv        Human adenovirus                              
  hadv_a      Human adenovirus A                            
  hadv_b      Human adenovirus B                            
  hadv_c      Human adenovirus C                            
  hadv_d      Human adenovirus D                            
  hadv_e      Human adenovirus E                            
  hadv_f      Human adenovirus F                            
  hadv_g      Human adenovirus G                            
  hhv         Human herpesvirus                             (hsv)
  hhv1        Human herpesvirus 1                           (hsv1)
  hhv2        Human herpesvirus 2                           
  hhv3        Human herpesvirus 3 (Varicellovirus)          (var)
  hhv4        Human herpesvirus 4                           
  hhv5        Human herpesvirus 5                           
  msl         Measles / Morbillivirus                       (measles)
  mpv         Metapneumovirus (MPV)                         
  mmp         Mumps / Rubulavirus                           (mumps)
  norv        Norovirus                                     (noro)
  norv_1      Norovirus I                                   (noro1, noro_1, norv1)
  norv_2      Norovirus II                                  (noro2, noro_2, norv2)
  norv_misc   Norovirus miscellaneous                       
  norv_mur    Norovirus murine                              
  rabies      Rabies                                        
  rsv         Respiratory syntactical virus (RSV)           
  respiro     Respirovirus                                  (resp)
  hpiv_1      Respirovirus HPIV-1                           (hpiv1)
  hpiv_3      Respirovirus HPIV-3                           (hpiv3)
  sendai      Respirovirus Sendai                           
  rtv         Rotavirus                                     (rota)
  rtv_a       Rotavirus A                                   (rota_a)
  rtv_b       Rotavirus B                                   (rota_b)
  rtv_c       Rotavirus C                                   (rota_c)
  rtv_f       Rotavirus F                                   (rota_f)
  rtv_g       Rotavirus G                                   (rota_g)
  rbl         Rubella                                       (rubella)
  sapo        Sapovirus                                     
  wnv         West Nile Virus - Both Lineage I and II       
  wnvI        West Nile Virus - Lineage I                   
  wnvII       West Nile Virus - Lineage II                  
  yfv         Yellow Fever / Japanese encephalitis (JEV)    (jev)

* non-standard grouping, must be invoked directly, not included in "any virus" via -A or as a
  subset of other -D specifications
