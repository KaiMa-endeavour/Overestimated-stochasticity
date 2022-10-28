# Overestimated-stochasticity
Stochastic process in microbial community assembly is overestimated due to random sampling.

Ma, K., & Tu, Q. (2022). Random sampling associated with microbial profiling leads to overestimated stochasticity inference in community assembly. Frontiers in Microbiology, 13, 1011269. https://doi.org/10.3389/fmicb.2022.1011269

## Methods
- The seed communities containing 10<sup>4</sup> taxa and 10<sup>8</sup> organisms based on taxonomic composition were constructed in this study. The ecological processes, such as dispersal and drift, were simulated by the R function `cons_seedcomm`, generating a series of seed communities with different β-diversity.

- Individuals with different numbers (5000 to 200,000) were randomly picked from the seed community to simulate the real communities. β-diversity and stochasticity of mock communities influenced by the different sequencing depths. The R function `cons_mockcomm` was used to construct mock communities.

- **Null models** have been developed to quantitatively disentangle the relative importance of deterministic vs. stochastic processes in structuring the compositional variations of biological communities. In this study, two contrasting null models are employed. Therefore, the parameter *abundance* = "shuffle" is added to `NST::tNST` and renamed as `tNSTmod`. The two null models generates heterogeneous and homogeneous null communities, corresponding to the *abundance* options "shuffle" and "region", respectively.

- In addition to the stochastic ratio approach, the **RC<SUB>bray</SUB>** metric was also employed to quantify the contribution of different ecological processes to the compositional variations of microbial communities. Because it was technically almost impossible to simulate the phylogenetic relationships representing the community assembly process of mock communities, only RC<SUB>bray</SUB> was used. The R function `Raup_Crick_Abundance` provided by [Stegen et al.](https://github.com/stegen/Stegen_etal_ISME_2013) for RC<SUB>bray</SUB> metric analysis.

The primary data used in the study are available at [https://github.com/KaiMa-endeavour/Overestimated-stochasticity/releases](https://github.com/KaiMa-endeavour/Overestimated-stochasticity/releases).

## References

Guo, X., Feng, J., Shi, Z., Zhou, X., Yuan, M., Tao, X., Hale, L., Yuan, T., Wang, J., Qin, Y., Zhou, A., Fu, Y., Wu, L., He, Z., Van Nostrand, J. D., Ning, D., Liu, X., Luo, Y., Tiedje, J. M., … Zhou, J. (2018). Climate warming leads to divergent succession of grassland microbial communities. Nature Climate Change, 8(9), Article 9. https://doi.org/10.1038/s41558-018-0254-2

May, F., Gerstner, K., McGlinn, D. J., Xiao, X., & Chase, J. M. (2018). mobsim: An r package for the simulation and measurement of biodiversity across spatial scales. Methods in Ecology and Evolution, 9(6), 1401–1408. https://doi.org/10.1111/2041-210x.12986

Ning, D., Deng, Y., Tiedje, J. M., & Zhou, J. (2019). A general framework for quantitatively assessing ecological stochasticity. Proceedings of the National Academy of Sciences of the United States of America, 116(34), 16892–16898. https://doi.org/10.1073/pnas.1904623116

Raup, D. M., & Crick, R. E. (1979). Measurement of Faunal Similarity in Paleontology. Journal of Paleontology, 53(5), 1213–1227.

Stegen, J. C., Lin, X., Fredrickson, J. K., Chen, X., Kennedy, D. W., Murray, C. J., Rockhold, M. L., & Konopka, A. (2013). Quantifying community assembly processes and identifying features that impose them. The ISME Journal, 7(11), 2069–2079. https://doi.org/10.1038/ismej.2013.93

Stegen, J. C., Lin, X., Fredrickson, J. K., & Konopka, A. E. (2015). Estimating and mapping ecological processes influencing microbial community assembly. Frontiers in Microbiology, 6, 370. https://doi.org/10.3389/fmicb.2015.00370

Zhou, J., Deng, Y., Zhang, P., Xue, K., Liang, Y., Van Nostrand, J. D., Yang, Y., He, Z., Wu, L., Stahl, D. A., Hazen, T. C., Tiedje, J. M., & Arkin, A. P. (2014). Stochasticity, succession, and environmental perturbations in a fluidic ecosystem. Proceedings of the National Academy of Sciences of the United States of America, 111(9), E836. https://doi.org/10.1073/pnas.1324044111

Zhou, J., & Ning, D. (2017). Stochastic Community Assembly: Does It Matter in Microbial Ecology? Microbiology and Molecular Biology Reviews, 81(4), e00002-17. https://doi.org/10.1128/mmbr.00002-17

