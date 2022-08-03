# Overestimated-stochasticity
Stochastic process in microbial community assembly is overestimated due to random sampling.

- The seed communities containing 10<sup>4</sup> taxa (i.e. OTUs) and 10<sup>8</sup> organisms (i.e. sequences) based on taxonomic composition were constructed in this study. The ecological processes of dispersal and drift were simulated by the R function `cons_seedcomm.R`, generating a series of seed communities with different β-diversity.

- Samples with different numbers (5000 to 200,000) of individuals were randomly picked from the seed communities to simulate real communities. β-diversity and stochasticity of mock communities influenced by the different sequencing depths. The R function `cons_mockcomm.R` was used to construct mock communities.

- **Null models** have been developed to quantitatively disentangle the relative importance of deterministic vs. stochastic processes in structuring the compositional variations of biological communities. In this study, two contrasting null models are employed. Therefore, the parameter *abundance* = "shuffle" is added to `NST::tNST` and renamed as `tNSTmod.R`. These two null models represent the heterogeneity and homogeneity of null communities, respectively, corresponding to the options "shuffle" and "region" for *abundance*.

- In addition to the stochastic ratio approach, the **RC<SUB>bray</SUB>** metric was also employed to quantify the contribution of different ecological processes to the compositional variations of microbial communities. Because it was technically almost impossible to simulate the phylogenetic relationships representing the community assembly process of mock communities, only RC<SUB>bray</SUB> was used. The R function `Raup_Crick_Abundance.r` provided by [Stegen et al.](https://github.com/stegen/Stegen_etal_ISME_2013) for RC<SUB>bray</SUB> metric analysis.

The primary data used in the study are available at [https://github.com/KaiMa-endeavour/Overestimated-stochasticity/releases/tag/data](https://github.com/KaiMa-endeavour/Overestimated-stochasticity/releases/tag/data).

## References

Guo, X., Feng, J., Shi, Z., Zhou, X., Yuan, M., Tao, X., Hale, L., Yuan, T., Wang, J., Qin, Y., Zhou, A., Fu, Y., Wu, L., He, Z., Van Nostrand, J.D., Ning, D., Liu, X., Luo, Y., Tiedje, J.M., Yang, Y., Zhou, J., 2018. Climate warming leads to divergent succession of grassland microbial communities. Nature Clim Change 8, 813–818. https://doi.org/10.1038/s41558-018-0254-2

May, F., Gerstner, K., McGlinn, D.J., Xiao, X., Chase, J.M., 2018. mobsim: An r package for the simulation and measurement of biodiversity across spatial scales. Methods in Ecology and Evolution 9, 1401–1408. https://doi.org/10.1111/2041-210x.12986

Ning, D., Deng, Y., Tiedje, J.M., Zhou, J., 2019. A general framework for quantitatively assessing ecological stochasticity. Proc. Natl. Acad. Sci. U.S.A. 116, 16892–16898. https://doi.org/10.1073/pnas.1904623116

Raup, D.M., Crick, R.E., 1979. Measurement of Faunal Similarity in Paleontology. Journal of Paleontology 53, 1213–1227.

Stegen, J.C., Lin, X., Fredrickson, J.K., Chen, X., Kennedy, D.W., Murray, C.J., Rockhold, M.L., Konopka, A., 2013. Quantifying community assembly processes and identifying features that impose them. The ISME Journal 7, 2069–2079. https://doi.org/10.1038/ismej.2013.93

Stegen, J.C., Lin, X., Fredrickson, J.K., Konopka, A.E., 2015. Estimating and mapping ecological processes influencing microbial community assembly. Frontiers in Microbiology 6. https://doi.org/10.3389/fmicb.2015.00370

Zhou, J., Deng, Y., Zhang, P., Xue, K., Liang, Y., Van Nostrand, J.D., Yang, Y., He, Z., Wu, L., Stahl, D.A., Hazen, T.C., Tiedje, J.M., Arkin, A.P., 2014. Stochasticity, succession, and environmental perturbations in a fluidic ecosystem. Proceedings of the National Academy of Sciences of the United States of America 111, E836. https://doi.org/10.1073/pnas.1324044111

Zhou, J., Ning, D., 2017. Stochastic Community Assembly: Does It Matter in Microbial Ecology? Microbiology and Molecular Biology Reviews 81. https://doi.org/10.1128/mmbr.00002-17


