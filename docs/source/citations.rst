===================
Citing and citation
===================

If you use ``cazomevolve`` in your work _please_ cite our work (including the provided DOI), as well as 
the specfic version of the tool you used. This is not only helpful to us as the developers to get our 
work out into the world, but it is also essential for the reproducibility and integrity of scientific research.

**Citation:**

   Hobbs, Emma. E. M., Gloster, Tracey, M., Pritchard, Leighton (2023) cazomevolve, _GitHub_. DOI: 10.5281/zenodo.6614827

``cazomevovle`` depends on a number of tools. To recognise the contributions that the 
authors and developers have made, please also cite the following:

**CAZy:**

``cazomevolve`` uses the CAZy family classifications establised and curated by the CAZy database.

   Elodie Drula and others, The carbohydrate-active enzyme database: functions and literature, Nucleic Acids Research, Volume 50, Issue D1, 7 January 2022, Pages D571–D577, https://doi.org/10.1093/nar/gkab1045

**cazy_webscraper:**

   Hobbs, E. E. M., Gloster, T. M., and Pritchard, L. (2022) 'cazy_webscraper: local compilation and interrogation of comprehensive CAZyme datasets', bioRxiv, https://doi.org/10.1101/2022.12.02.518825

For additional citations when using ``cazy_webscraper``, see the ``cazy_webscraper`` `documentation <https://cazy-webscraper.readthedocs.io/en/latest/citation.html>`_.

**ncbi-genome-download:**

   Blin et al. (2017) ncbi-genome-download, https://github.com/kblin/ncbi-genome-download

**dbCAN:**

dbCAN version 2:

   Zhang H, Yohe T, Huang L, Entwistle S, Wu P, Yang Z, Busk PK, Xu Y, Yin Y. dbCAN2: a meta server for automated carbohydrate-active enzyme annotation. Nucleic Acids Res. 2018 Jul 2;46(W1):W95-W101. doi: 10.1093/nar/gky418. PMID: 29771380; PMCID: PMC6031026.

dbCAN version 3:

If using dbCAN version 3, cite the publication for version 2 as well as eCAMI:

   Xu J, Zhang H, Zheng J, Dovoedo P, Yin Y. eCAMI: simultaneous classification and motif identification for enzyme annotation. Bioinformatics. 2020 Apr 1;36(7):2068-2075. doi: 10.1093/bioinformatics/btz908. PMID: 31794006.

dbCAN version 4:

   Zheng J, Ge Q, Yan Y, Zhang X, Huang L, Yin Y. dbCAN3: automated carbohydrate-active enzyme and substrate annotation. Nucleic Acids Res. 2023 May 1:gkad328. doi: 10.1093/nar/gkad328. Epub ahead of print. PMID: 37125649.

**BLAST Score Ratio:**

    Rasko DA, Myers GS, Ravel J. Visualization of comparative genomic analyses by BLAST score ratio. BMC Bioinformatics. 2005 Jan 5;6:2. doi: 10.1186/1471-2105-6-2. PMID: 15634352; PMCID: PMC545078.

**DIAMOND**

   Buchfink B, Xie C, Huson DH. Fast and sensitive protein alignment using DIAMOND. Nat Methods. 2015 Jan;12(1):59-60. doi: 10.1038/nmeth.3176. Epub 2014 Nov 17. PMID: 25402007.

**BLAST**

   Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. Basic local alignment search tool. Journal of molecular biology. 1990;215(3):403–10.

**Explore CAZy family sequence diversity and CAZomes:**

``cazomevolve`` uses several packages to visualise and interrogate the dataset.

   @article{Waskom2021,
      doi={10.21105/joss.03021},
      url={https://doi.org/10.21105/joss.03021},
      year={2021},
      publisher={The Open Journal},
      volume={6},
      number={60},
      pages={3021},
      author={Michael L. Waskom},
      title={seaborn: statistical data visualization},
      journal={Journal of Open Source Software}
   }

   @article{Virtanen2020,
      author ={Virtanen, Pauli and Gommers, Ralf and Oliphant, Travis E. and
                  Haberland, Matt and Reddy, Tyler and Cournapeau, David and
                  Burovski, Evgeni and Peterson, Pearu and Weckesser, Warren and
                  Bright, Jonathan and {van der Walt}, St{\'e}fan J. and
                  Brett, Matthew and Wilson, Joshua and Millman, K. Jarrod and
                  Mayorov, Nikolay and Nelson, Andrew R. J. and Jones, Eric and
                  Kern, Robert and Larson, Eric and Carey, C J and
                  Polat, {\.I}lhan and Feng, Yu and Moore, Eric W. and
                  {VanderPlas}, Jake and Laxalde, Denis and Perktold, Josef and
                  Cimrman, Robert and Henriksen, Ian and Quintero, E. A. and
                  Harris, Charles R. and Archibald, Anne M. and
                  Ribeiro, Ant{\^o}nio H. and Pedregosa, Fabian and
                  {van Mulbregt}, Paul and {SciPy 1.0 Contributors}},
      title  ={{{SciPy} 1.0: Fundamental Algorithms for Scientific
                  Computing in Python}},
      journal={Nature Methods},
      year   ={2020},
      volume ={17},
      pages  ={261--272},
   }

   @article{scikit-learn,
      title={Scikit-learn: Machine Learning in {P}ython},
      author={Pedregosa, F. and Varoquaux, G. and Gramfort, A. and Michel, V.
               and Thirion, B. and Grisel, O. and Blondel, M. and Prettenhofer, P.
               and Weiss, R. and Dubourg, V. and Vanderplas, J. and Passos, A. and
               Cournapeau, D. and Brucher, M. and Perrot, M. and Duchesnay, E.},
      journal={Journal of Machine Learning Research},
      volume={12},
      pages={2825--2830},
      year={2011}
   }