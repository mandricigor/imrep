# ImReP


## One sentence summary

ImReP is a computational method for rapid and accurate profiling of the adaptive immune repertoire from regular RNA-Seq data.

## Brief description

ImReP is able to quantify individual immune response based on a recombination landscape of genes encoding B and T cell receptors (BCR and TCR).  ImReP is able to efficiently extract TCR- and BCR- derived reads from the RNA-Seq data and accurately assemble clonotypes (defined as clones with identical CDR3 amino acid sequences)  and detect corresponding V(D)J recombinations.Using CAST clustering technique, ImReP is able to correct assembled clonotypes for PCR and sequencing errors.




# Tutorial 
ImReP tutorial is available at https://github.com/mandricigor/imrep/wiki





## Contact

This software was developed by [Igor Mandric](https://github.com/mandricigor) and [Serghei Mangul](https://sergheimangul.wordpress.com/). Please do not hesitate to contact us (imandric1@student.gsu.edu, smangul@ucla.edu) if you have any comments, suggestions, or clarification requests regarding the tutorial or if you would like to contribute to this resource.


## Publication

Mangul, Serghei and Mandric Igor, et al., “Profiling adaptive immune repertoires across multiple human tissues by RNA Sequencing“, bioRxiv (2016)

## Presentations

Slides from the talk at PSB 2017 workshop on Open Data for Discovery Science  are available here : https://sergheimangul.files.wordpress.com/2017/01/psb2017_public1.pdf 


## The Atlas of Immunoglobilin Repertoires (TheAIR)

Using ImreP we have created the Atlas of Immunoglobilin Repertoires (TheAIR). The AIR is a collection of CDR3 from Immunoglobulin (Ig) receptor repetoires. CDR3s are assembled from 8,555 human RNA-seq samples across 544 individuals from 53 tissues from the Genotype-Tissue Expression (GTEx v6) project using ImReP.The AIR has one of the largest collection of CDR3 sequences (n=3.6 million) and tissue types (n=53). TheAIR is freely available at : https://smangul1.github.io/TheAIR/


# Releases

## ImReP 0.3 release 01/04/2017

ImReP 0.3 is available for download here.  This is a new release with the following fixes and changes:

- the -o/ –overlapLen allows to set up the minimal overlap length for reads to be assembled into CDR3 clonotype (second stage of ImReP)
- the -extendedOutput options allows producing two additional files: full_cdr3.txt and partial_cdr3.txt containing detailed information for each read. Details about the format are available at https://github.com/mandricigor/imrep/wiki/ImReP-output
- toy example composed of simulated receptor-derived reads is distributed with the tool. Fasta file with simulated reads is located under example directory

## ImReP 0.2 release 11/30/2016

ImReP 0.2 is available for download here.  This is a new release with the following fixes and changes:

- We have added detailed tutorial available at https://github.com/mandricigor/imrep/wiki
- the -t option has been added to control the stringency of CDR3 clustering.  The value can be from 0.0 to 1.0. Note that thresholds near 1.0 are more liberal and result in more CDR3 to be reported. The default value is 0.2.
- the -c option has been added. Chain types should be separated by comma. The default value is IGH,IGK,IGL,TRA,TRB,TRD,TRG.
- the –noOverlapStep option. This allows to skip the second stage of the ImReP assembly.

## ImReP 0.1 release 11/21/2016

ImReP 0.1 is available for download here. The first public release of ImReP. Because this is the first release, the manual is very limited. Only the basic options have been described, but we plan to update it frequently. If you have any questions about how ImReP works, please contact  Igor Mandric (mandric.igor@gmail.com) and Serghei Mangul (smangul@ucla.edu)


