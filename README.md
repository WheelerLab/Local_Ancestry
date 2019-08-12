# Local Ancestry

My final project for Advanced Bioinformatics (BIOI 500) was a comparison of three local ancestry estimation softwares: LAMP-LD, RFMix, and ELAI. This repository builds upon the work in the original class project and runs six populations of different ancestral backgrounds and proportions with these softwares to compare accuracy and resource usage. We work to convert the genotypic format requirements of these softwares and streamline their usage. We will use these local ancestry softwares in real genotypic and transcriptomic data from the Multi-Ethnic Study of Atherosclerosis (MESA) and Modeling the Epidemiologic Transition Study to improve gene expression prediction in diverse populations. Analyses performed include:

* [Generating ancestral proportions from 1000G populations](https://github.com/aandaleon/Local_Ancestry/blob/master/01_get_1000G_proportions.R)
* [Simulating genotypes from 1000G CEU and YRI](https://github.com/aandaleon/Local_Ancestry/blob/master/02_sample_refs.R)
* [Also see adsim scripts and notes](https://github.com/WheelerLab/Local_Ancestry/tree/master/adsim)
 
<center>

|      | NAT | CEU | YRI |
|------|-----|-----|-----|
| MXL  | 60% | 40% | 0%  |
| PUR  | 16% | 84% | 0%  |
| ACB  | 3%  | 10% | 87% |
| ASW  | 1%  | 26% | 73% |
| EVEN | 33% | 33% | 33% |
| [BRYC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4289685/) | 18% | 65% | 6%  |

</center>

* Converting to software format and run
  * [VCF to LAMP-LD](https://github.com/aandaleon/Local_Ancestry/blob/master/03_LAMP_LD.sh)
  * [VCF to RFMix](https://github.com/aandaleon/Local_Ancestry/blob/master/04_RFMix.py)
  * [VCF to ELAI](https://github.com/aandaleon/Local_Ancestry/blob/master/05_ELAI.sh)
  * [VCF to MOSAIC](https://github.com/WheelerLab/Local_Ancestry/blob/master/06_MOSAIC.sh)
  * [VCF to Loter](https://github.com/WheelerLab/Local_Ancestry/blob/master/07_Loter.sh)
* Measure [accuracy](https://github.com/WheelerLab/Local_Ancestry/tree/master/accuracy_estimation) and resource use of softwares

## Links
* BIOI 500 summaries
  * [Paper](https://docs.google.com/document/d/1YVKdQc5WmQK9lJFTgaxiBgNVxHsKWZ4x_5nJwLcCf-Q/edit?usp=sharing)
  * [Presentation](https://docs.google.com/presentation/d/1LJXkkStQNTmNCr_zg6MTAqe_HDryur0adjWiqwoM06g/edit?usp=sharing)
* LAMP-LD 
  * [Paper](https://academic.oup.com/bioinformatics/article/28/10/1359/212139)
  * [Software](http://lamp.icsi.berkeley.edu/lamp/lampld/)
* RFMix 
  * [Paper](https://www.sciencedirect.com/science/article/pii/S0002929713002899?via%3Dihub)
  * [Software](https://sites.google.com/site/rfmixlocalancestryinference/)
* ELAI
  * [Paper](http://www.genetics.org/content/196/3/625.long)
  * [Software](http://www.haplotype.org/elai.html)
* MOSAIC
  * [Paper](https://www.genetics.org/content/212/3/869)
  * [Software](https://maths.ucd.ie/~mst/MOSAIC/)
* Loter
  * [Paper](https://academic.oup.com/mbe/article/35/9/2318/5040668)
  * [Software](https://github.com/bcm-uga/Loter)

