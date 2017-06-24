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

Using ImReP we have created the Atlas of Immunoglobilin Repertoires (TheAIR). The AIR is a collection of CDR3 from Immunoglobulin (Ig) receptor repetoires. CDR3s are assembled from 8,555 human RNA-seq samples across 544 individuals from 53 tissues from the Genotype-Tissue Expression (GTEx v6) project using ImReP.The AIR has one of the largest collection of CDR3 sequences (n=3.6 million) and tissue types (n=53). TheAIR is freely available at : https://smangul1.github.io/TheAIR/


# Releases

To get the latest version, use 
```
git clone https://github.com/mandricigor/imrep.git
```
Information about releases is availabe here 
https://github.com/mandricigor/imrep/releases





