---
name: Psychometrics
topic: Psychometric Models and Methods
maintainer: Patrick Mair, Yves Rosseel, Kathrin Gruber
email: mair@fas.harvard.edu
version: 2023-12-15
source: https://github.com/cran-task-views/Psychometrics/
---


Psychometrics is concerned with theory and techniques of psychological
measurement. Psychometricians have also worked collaboratively with
those in the field of statistics and quantitative methods to develop
improved ways to organize, analyze, and scale corresponding data. Since
much functionality is already contained in base R and there is
considerable overlap between tools for psychometry and tools described
in other views, we only give a brief overview of packages that are closely related to
psychometric methodology.

Contributions are always welcome and encouraged, either via e-mail to the
maintainer or by submitting an issue or pull request in the GitHub
repository linked above.


### Item Response Theory (IRT):

-   The `r pkg("eRm", priority = "core")` package fits
    extended Rasch models, i.e. the ordinary Rasch model for dichotomous
    data (RM), the linear logistic test model (LLTM), the rating scale
    model (RSM) and its linear extension (LRSM), the partial credit
    model (PCM) and its linear extension (LPCM) using conditional ML
    estimation. Missing values are allowed.
-   The package `r pkg("ltm", priority = "core")` also fits
    the simple RM. Additionally, functions for estimating Birnbaum's 2-
    and 3-parameter models based on a marginal ML approach are
    implemented as well as the graded response model for polytomous
    data, and the linear multidimensional logistic model.
-   The `r pkg("mirt", priority = "core")` estimates
    dichotomous and polytomous response data using unidimensional and
    multidimensional latent trait models under the IRT paradigm.
    Exploratory and confirmatory models can be estimated with quadrature
    (EM) or stochastic (MHRM) methods. Confirmatory bi-factor and
    two-tier analyses are available for modeling item testlets. Multiple
    group analysis and mixed effects designs also are available for
    detecting differential item functioning and modeling item and person
    covariates.
-   `r pkg("TAM", priority = "core")` fits unidimensional and
    multidimensional item response models and also includes multifaceted
    models, latent regression models and options for drawing plausible
    values.
-   `r pkg("Dire")` fits weighted latent variable linear models, estimating score
    distributions for groups of people, in an IRT framework,
    following Cohen and Jiang (1999). It can then draw plausible values.
-   `r pkg("PLmixed")` fits (generalized) linear mixed
    models (GLMM) with factor structures.
-   `r pkg("MLCIRTwithin")` provides a flexible framework
    for the estimation of discrete two-tier IRT models for the analysis
    of dichotomous and ordinal polytomous item responses.
-   `r pkg("IRTShiny")` provides an interactive shiny
    application for IRT analysis.
-   Some additional uni- and multidimensional item response models
    (especially for locally dependent item responses) and some
    exploratory methods (DETECT, LSDM, model-based reliability) are
    included in `r pkg("sirt")`.
-   The `r pkg("pcIRT")` estimates the multidimensional
    polytomous Rasch model and the Mueller's continuous rating scale
    model.
-   An implementation of the partial credit model with response styles
    is given in the `r pkg("PCMRS")`.
-   `r pkg("MultiLCIRT")` estimates IRT models under (1)
    multidimensionality assumption, (2) discreteness of latent
    traits, (3) binary and ordinal polytomous items.
-   Conditional maximum likelihood estimation via the EM algorithm and
    information-criterion-based model selection in binary mixed Rasch
    models are implemented in the `r pkg("psychomix")`
    package.
-   The `r pkg("PP")` package includes estimation of (MLE,
    WLE, MAP, EAP, ROBUST) person parameters for the 1,2,3,4-PL model
    and the GPCM (generalized partial credit model). The parameters are
    estimated under the assumption that the item parameters are known
    and fixed. The package is useful e.g. in the case that items from an
    item pool/item bank with known item parameters are administered to a
    new population of test-takers and an ability estimation for every
    test-taker is needed.
-   The `r pkg("equateIRT")` package computes direct, chain
    and average (bisector) equating coefficients with standard errors
    using Item Response Theory (IRT) methods for dichotomous items.
    `r pkg("equateMultiple")` can be used for equating of
    multiple forms using IRT methods.
-   `r pkg("kequate")` implements the kernel method of test
    equating using the CB, EG, SG, NEAT CE/PSE and NEC designs,
    supporting gaussian, logistic and uniform kernels and unsmoothed and
    pre-smoothed input data.
-   The `r pkg("EstCRM")` package calibrates the parameters
    for Samejima's Continuous IRT Model via EM algorithm and Maximum
    Likelihood. It allows to compute item fit residual statistics, to
    draw empirical 3D item category response curves, to draw theoretical
    3D item category response curves, and to generate data under the CRM
    for simulation studies.
-   The `r pkg("difR")` package contains several traditional
    methods to detect DIF in dichotomously scored items. Both uniform
    and non-uniform DIF effects can be detected, with methods relying
    upon item response models or not. Some methods deal with more than
    one focal group.
-   The package `r pkg("lordif")` provides a logistic
    regression framework for detecting various types of DIF.
-   `r pkg("DIFplus")` allows users to implement extensions
    of the Mantel-Haenszel DIF detection procedures in the presence of
    multilevel data.
-   `r pkg("DIFlasso")` implements a penalty approach to
    differential item functioning in Rasch models. It can handle
    settings with multiple (metric) covariates.
-   `r pkg("GPCMlasso")` provides a function to detect DIF
    in generalized partial credit models (GPCM).
-   `r pkg("DIFtree")` performs recursive partitioning for
    simultaneous selection of items and variables that induce DIF in
    dichotomous or polytomous items.
-   `r pkg("DIFboost")` can be used for DIF detection in
    Rasch models by boosting techniques.
-   A set of functions to perform differential item and item functioning
    analyses is implemented in the `r pkg("DFIT")` package.
    It includes functions to use the Monte Carlo item parameter
    replication (IPR) approach for obtaining the associated statistical
    significance tests cut-off points.
-   The `r pkg("difNLR")` package uses nonlinear regression
    to estimate DIF.
-   The `r pkg("catR")` package allows for computarized
    adaptive testing using IRT methods.
-   The `r pkg("mirtCAT")` package provides tools to
    generate an HTML interface for creating adaptive and non-adaptive
    educational and psychological tests using the shiny package.
    Suitable for applying unidimensional and multidimensional
    computerized adaptive tests using IRT methodology and for creating
    simple questionnaires forms to collect response data directly in R.
-   `r pkg ("D3mirt")` for identifying, estimating, and plotting descriptive multidimensional item response theory models, restricted to 3D and dichotomous or polytomous data that fit the two-parameter logistic model or the graded response model.
-   `r pkg("xxIRT")` is implementation of related to IRT and
    computer-based testing.
-   `r pkg("Rirt")` estimates the 3-parameter-logistic model, generalized partial credit model, and graded response model.
-   Explicit calculation (not estimation) of Rasch item parameters
    (dichotomous and polytomous) by means of a pairwise comparison
    approach can be done using the `r pkg("pairwise")`
    package.
-   Multilevel Rasch models can be estimated using
    `r pkg("lme4")`, `r pkg("nlme")`, and
    `r pkg("MCMCglmm")` with crossed or partially crossed
    random effects. `r pkg("GLMMRR")` adds some flexibility
    in terms of link functions, whereas `r pkg("ordinal")`
    can be used for polytomous models. An infrastructure for estimating
    tree-structured item response models of the GLMM family using
    `r pkg("lme4")` is provided in
    `r pkg("irtrees")`.
-   Nonparametric IRT analysis can be computed by means if the
    `r pkg("mokken", priority = "core")` package. It includes
    an automated item selection algorithm, and various checks of model
    assumptions.
-   Nonparametric IRT for nonmonotonic IRFs of proximity data can be
    fitted using the `r pkg("mudfold")` package.
-   `r pkg("RaschSampler")` allows the construction of exact
    Rasch model tests by generating random zero-one matrices with given
    marginals.
-   Statistical power simulation for testing the Rasch model based on a
    three-way ANOVA design with mixed classification can be carried out
    using `r pkg("pwrRasch")`.
-   Tools to assess model fit and identify misfitting items for Rasch
    models and PCMS are implemented in `r pkg("iarm")`. It
    includes item fit statistics, ICCs, item-restscore association,
    conditional likelihood ratio tests, assessment of measurement error,
    estimates of the reliability and test targeting.
-   `r pkg("cacIRT")` computes classification accuracy and
    consistency under Item Response Theory. Implements total score and
    latent trait IRT methods as well as total score kernel-smoothed
    methods.
-   The package `r pkg("irtoys")` provides a simple common
    interface to the estimation of item parameters in IRT models for
    binary responses with three different programs (ICL, BILOG-MG, and
    ltm, and a variety of functions useful with IRT models.
-   The `r pkg("CDM")` estimates several cognitive diagnosis
    models (DINA, DINO, GDINA, RRUM, LCDM, pGDINA, mcDINA), the general
    diagnostic model (GDM) and structured latent class analysis (SLCA).
-   Gaussian ordination, related to logistic IRT and also approximated
    as maximum likelihood estimation through canonical correspondence
    analysis is implemented in various forms in the package
    `r pkg("VGAM")`.
-   `r pkg("immer")` implements some item response models
    for multiple ratings, including the hierarchical rater model and a
    wrapper function to the commercial FACETS program.
-   The `r pkg("latdiag")` package produces commands to
    drive the dot program from graphviz to produce a graph useful in
    deciding whether a set of binary items might have a latent scale
    with non-crossing ICCs.
-   The purpose of the `r pkg("rpf")` package is to factor
    out logic and math common to IRT fitting, diagnostics, and analysis.
    It is envisioned as core support code suitable for more specialized
    IRT packages to build upon.
-   `r pkg("WrightMap")` provides graphical tools for
    plotting item-person maps.
-   `r pkg("irtDemo")` includes a collection of shiny
    applications to demonstrate or to explore fundamental IRT concepts.
    `r pkg("ifaTools")` is a shiny interface to IRT with
    `r pkg("OpenMx", priority = "core")`.
-   IRT utility functions described in the Baker/Kim book are included
    in `r pkg("birtr")`.
-   Convenience functions to use and automate IRT modeling for judgement
    data are implemented in `r pkg("jrt")`.
-   The `r pkg("conquestr")` package allows users to call
    ACER ConQuest from within R.
-   `r pkg("TestDesign")` implements optimal test design
    approaches for fixed and adaptive test construction.
-   `r pkg("PROsetta")` provides functions for performing scale
    linking between instruments, based on the PROsetta Stone methodology.
-   `r pkg("maat")` performs simulations for multiple administrations adaptive testing.

### Correspondence Analysis (CA), Optimal Scaling:

-   The package `r pkg("ca", priority = "core")` comprises
    two parts, one for simple correspondence analysis and one for
    multiple and joint correspondence analysis.
-   Simple and canonical CA are provided by the package
    `r pkg("anacor", priority = "core")`, including
    confidence ellipsoids. It allows for different scaling methods such
    as standard scaling, Benzecri scaling, centroid scaling, and Goodman
    scaling.
-   Homogeneity analysis aka multiple CA and various Gifi extensions can
    be computed by means of the `r pkg("Gifi", priority = "core")` package,
    which replaces `r pkg("homals")`. This
    package includes various other optimal scaling methods such as
    Morals (monotone regression), Princals (nonlinear PCA), Homals (multiple correspondence analysis), etc.
-   Simple and multiple correspondence analysis can be performed using
    `corresp()` and `mca()` in package `r pkg("MASS")`.
-   The package `r pkg("ade4", priority = "core")` contains
    an extensive set of functions covering, e.g., principal components,
    simple and multiple, fuzzy, non symmetric, and decentered
    correspondence analysis. Additional functionality is provided at
    [Bioconductor](http://www.bioconductor.org/) in the package `made4`
    (see also [here](http://bioinf.ucd.ie/people/aedin/R/) ).
-   The package `r pkg("cocorresp")` fits predictive and
    symmetric co-correspondence analysis (CoCA) models to relate one
    data matrix to another data matrix.
-   Apart from several factor analytic methods
    `r pkg("FactoMineR")` performs CA including
    supplementary row and/or column points and multiple correspondence
    analysis (MCA) with supplementary individuals, supplementary
    quantitative variables and supplementary qualitative variables.
-   Package `r pkg("vegan", priority = "core")` supports all
    basic ordination methods, including non-metric multidimensional
    scaling. The constrained ordination methods include constrained
    analysis of proximities, redundancy analysis, and constrained
    (canonical) and partially constrained correspondence analysis.
-   `r pkg("cabootcrs")` computes bootstrap confidence
    regions for CA.
-   `r pkg("cncaGUI")` implements a GUI with which users can
    construct and interact with canonical (non-symmetrical) CA.
-   SVD based multivariate exploratory methods such as PCA, CA, MCA (as
    well as a Hellinger form of CA), generalized PCA are implemented in
    `r pkg("ExPosition")`. The package also allows for
    supplementary data projection.
-   `r pkg("cds")` can be used for constrained dual scaling
    for detecting response styles.
-   `r pkg("CAvariants")` provides six variants of two-way
    CA: simple, singly ordered, doubly ordered, non-symmetrical, singly
    ordered non-symmetrical ca, and doubly ordered non-symmetrical.
-   `r pkg("MCAvariants")` provides MCA and ordered MCA via
    orthogonal polynomials.
-   Specific and class specific MCA on survey-like data can be fitted
    using `r pkg("soc.ca")`.
-   `r pkg("optiscale")` provides tools for performing an
    optimal scaling transformation on a data vector.
-   A general framework of optimal scaling methods is implemented in the
    `r pkg("aspect")`.
-   `r pkg("candisc")`: Visualizing generalized canonical discriminant and canonical correlation analysis.

### Factor Analysis (FA), Principal Component Analysis (PCA):

-   Exploratory FA is the package stats as function `factanal()` and
    `fa()` and `fa.poly()` (ordinal data) in
    `r pkg("psych", priority = "core")`.
-   `r pkg("BayesFM")` for Bayesian EFA.
-   `r pkg("esaBcv")` estimates the number of latent factors
    and factor matrix.
-   `r pkg("SparseFactorAnalysis")` scales count and binary
    data with sparse FA.
-   `r pkg("EFAutilities")` computes robust standard errors
    and factor correlations under a variety of conditions.
-   `r pkg("faoutlier")` implements influential case
    detection methods for FA and SEM.
-   The package `r pkg("psych")` includes functions such as
    `fa.parallel()` and `VSS()` for estimating the appropriate number of
    factors/components as well as `ICLUST()` for item clustering.
-   PCA can be fitted with `prcomp()` (based on `svd()`, preferred) as
    well as `princomp()` (based on `eigen()` for compatibility with
    S-PLUS). Additional rotation methods for FA based on gradient
    projection algorithms can be found in the package
    `r pkg("GPArotation")`. The package
    `r pkg("nFactors")` produces a non-graphical solution to
    the Cattell scree test. Some graphical PCA representations can be
    found in the `r pkg("psy", priority = "core")` package.
    `r pkg("paran")` implements Horn's test of principal
    components/factors.
-   FA and PCA with supplementary individuals and supplementary
    quantitative/qualitative variables can be performed using the
    `r pkg("FactoMineR")` package whereas
    `r pkg("MCMCpack")` has some options for sampling from
    the posterior for ordinal and mixed factor models.
-   The `r pkg("Gifi")` package implements Princals, a PCA
    version for mixed-scale level input data.
-   `r pkg("nsprcomp")` and
    `r pkg("elasticnet")` fit sparse PCA.
-   Threeway PCA models (Tucker, Parafac/Candecomp) can be fitted using
    `r pkg("PTAk")`, `r pkg("ThreeWay")`, and
    `r pkg("multiway")`.
-   Independent component analysis (ICA) can be computed using
    `r pkg("fastICA")`, `r pkg("ica")`, and `r pkg("eegkit")` (designed for EEG data).
-   A desired number of robust principal components can be computed with
    the `r pkg("pcaPP")` package.
-   `r pkg("bpca")` implements 2D and 3D biplots of
    multivariate data based on PCA and diagnostic tools of the quality
    of the reduction.
-   `r pkg("missMDA")` provides imputation of incomplete
    continuous or categorical datasets in principal component analysis
    (PCA), multiple correspondence analysis (MCA) model, or multiple
    factor analysis (MFA) model.

### Structural Equation Models (SEM):

-   The package `r pkg("lavaan", priority = "core")` can be
    used to estimate a large variety of multivariate statistical models,
    including path analysis, confirmatory factor analysis, structural
    equation modeling and growth curve models. It includes the lavaan
    model syntax which allows users to express their models in a compact
    way and allows for ML, GLS, WLS, robust ML using Satorra-Bentler
    corrections, and FIML for data with missing values. It fully
    supports for meanstructures and multiple groups and reports
    standardized solutions, fit measures, modification indices and more
    as output.
-   The `r pkg("OpenMx", priority = "core")` package allows for the estimation
    of a wide variety of advanced multivariate statistical models. It
    consists of a library of functions and optimizers that allow you to
    quickly and flexibly define an SEM model and estimate parameters
    given observed data.
-   The `r pkg("sem")` package fits
    general (i.e., latent-variable) SEMs by FIML, and structural
    equations in observed-variable models by 2SLS. Categorical variables
    in SEMs can be accommodated via the `r pkg("polycor")`
    package.
-   `r pkg("tidySEM")` provides a tidy workflow for generating, estimating, reporting, and plotting structural equation models using lavaan, OpenMx, or Mplus.
-   `r pkg("SEMsens")` performs sensitivity analysis for omitted confounders in structural equation models using meta-heuristic optimization methods. 
-   `r pkg("lslx")` fits semi-confirmatory SEM via penalized
    likelihood with elastic net or minimax concave penalty.
-   The `r pkg("lavaan.survey")` package allows for complex
    survey structural equation modeling (SEM). It fits structural
    equation models (SEM) including factor analysis, multivariate
    regression models with latent variables and many other latent
    variable models while correcting estimates, standard errors, and
    chi-square-derived fit measures for a complex sampling design. It
    incorporates clustering, stratification, sampling weights, and
    finite population corrections into a SEM analysis.
-   The `r pkg("nlsem")` package fits nonlinear structural
    equation mixture models using the EM algorithm. Three different
    approaches are implemented: LMS (Latent Moderated Structural
    Equations), SEMM (Structural Equation Mixture Models), and NSEMM
    (Nonlinear Structural Equations Mixture Models).
-   A collection of functions for conducting meta-analysis using a
    structural equation modeling (SEM) approach via OpenMx is provided
    by the `r pkg("metaSEM")` package.
-   A general implementation of a computational framework for latent
    variable models (including structural equation models) is given in
    `r pkg("lava")`.
-   The `r pkg("pls")` package can be used for partial
    least-squares estimation. The package `r pkg("cSEM")`
    fits structural equation models using composite based approaches (e.g., PLS).
-   `r pkg("simsem")` is a package designed to aid in Monte
    Carlo simulations using SEM (for methodological investigations,
    power analyses and much more).
-   `r pkg("Sim.DiffProc")` provides a framework for
    parallelized Monte Carlo simulation-estimation in multidimensional
    continuous-time models, which have been implemented as SEM.
-   `r pkg("semTools")` is a package of add on functions
    that can aid in fitting SEMs in R (for example one function
    automates imputing missing data, running imputed datasets and
    combining the results from these datasets).
-   `r pkg("semPlot")` produces path diagrams and visual
    analysis for outputs of various SEM packages.
-   `r pkg("plotSEMM")` for graphing nonlinear relations
    among latent variables from structural equation mixture models.
-   `r pkg("influence.SEM")` implements outlier, leverage
    diagnostics, and case influence for SEM.
-   `r pkg("piecewiseSEM")` fits piecewise SEM.
-   `r pkg("rsem")` implements robust SEM with missing data
    and auxiliary variables.
-   `r pkg("regsem")` performs Regularization on SEM.
-   Recursive partitioning (SEM trees, SEM forests) is implemented in
    `r pkg("semtree")`.
-   `r pkg("lsl")` conducts SEM via penalized likelihood
    (latent structure learning).
-   `r pkg("MIIVsem")` contains functions for estimating
    structural equation models using instrumental variables.
-   The `r pkg("systemfit")` package implements a wider
    variety of estimators for observed-variables models, including
    nonlinear simultaneous-equations models.
-   `r pkg("STARTS")` contains functions for estimating the
    STARTS model.
-   Interfaces between R and other SEM software:
    `r pkg("REQS")`, `r pkg("MplusAutomation")`,
    and `r pkg("lisrelToR")`.

### Multidimensional Scaling (MDS):

-   The `r pkg("smacof", priority = "core")` package provides
    many approaches to metric and nonmetric MDS, including extensions
    for MDS with external constraints, spherical MDS, asymmetric MDS,
    three-way MDS (INDSCAL/IDIOSCAL), Bentler-Weeks model,
    unidimensional scaling, Procrustes, inverse MDS.
-   `r pkg("smacofx")` for flexible MDS analyses including MULTISCALE, Sammon mapping, ALSCAL, local MDS, elastic scaling, Box-Cox MDS, POST-MDS, curvilinear component and distance analysis, etc.
-   `r pkg("cops")` for cluster optimized prozimity scaling pronouncing the clustered appearance of the configuration.
-   `r pkg("stops")` provides a collection of methods that fit nonlinear distance transformations in multidimensional scaling (MDS) and trade-off the fit with structure considerations to find optimal parameters also known as structure optimized proximity scaling.
-   `r pkg("MASS")` and stats provide functionalities for
    computing classical MDS using the `cmdscale()` function. Sammon
    mapping `sammon()` and non-metric MDS `isoMDS()` are other relevant
    functions.
-   Nonmetric MDS can also be computed with `metaMDS()` in
    `r pkg("vegan")`. Furthermore,
    `r pkg("labdsv")` and `r pkg("ecodist")`
    provide the function `nmds()`. Also, the
    `r pkg("ExPosition")` implements a function for metric
    MDS.
-   Principal coordinate analysis can be computed with `capscale()` in
    `r pkg("vegan")`; in `r pkg("labdsv")` and
    `r pkg("ecodist")` using `pco()` and with `dudi.pco()`
    in `r pkg("ade4")`.
-   INDSCAL is also implemented in the `r pkg("SensoMineR")`
    package.
-   The package `r pkg("MLDS")` allows for the computation
    of maximum likelihood difference scaling (MLDS).
-   `r pkg("DistatisR")` implements the DiSTATIS/CovSTATIS
    3-way metric MDS approach.
-   `r pkg("munfold")` provides functions for metric
    unfolding.
-   The `r pkg("asymmetry")` package implements the
    slide-vector model for asymmetric MDS.
-   `r pkg("semds")` fits asymmetric and three-way MDS
    within an SEM framework.
-   `r pkg("cops")` performs cluster optimized proximity
    scaling which refers to MDS methods that aim at pronouncing the
    clustered appearance of the configuration.

### Classical Test Theory (CTT):

-   The `r pkg("CTT", priority = "core")` package can be used
    to perform a variety of tasks and analyses associated with classical
    test theory: score multiple-choice responses, perform reliability
    analyses, conduct item analyses, and transform scores onto different
    scales.
-   For multilevel model ICC for slope heterogeneity see
    `r pkg("iccbeta")`.
-   An interactive shiny application for CTT is provided by
    `r pkg("CTTShiny")`.
-   The `r pkg("cocron")` package provides functions to
    statistically compare two or more alpha coefficients based on either
    dependent or independent groups of individuals.
-   The `r pkg("betafunctions")` package includes an
    implementation of the so-called "Livingston and Lewis" approach to
    classification accuracy and consistency.
-   Cronbach alpha, kappa coefficients, and intra-class correlation
    coefficients (ICC) can be found in the `r pkg("psy")`
    package. Functions for ICC computation can be also found in the
    packages `r pkg("psych")`, and `r pkg("ICC")`.
-   A number of routines for scale construction and reliability analysis
    useful for personality and experimental psychology are contained in
    the package `r pkg("psych")`.
-   `r pkg("subscore")` can be used for computing subscores
    in CTT and IRT.
-   The quantifying construct validity procedure is implemented in
    `r pkg("qcv")`.

### Knowledge Structure Analysis:

-   `r pkg("DAKS")` provides functions and example datasets
    for the psychometric theory of knowledge spaces. This package
    implements data analysis methods and procedures for simulating data
    and transforming different formulations in knowledge space theory.
-   The `r pkg("kst")` package contains basic functionality
    to generate, handle, and manipulate deterministic knowledge
    structures based on sets and relations. Functions for fitting
    probabilistic knowledge structures are included in the
    `r pkg("pks")` package.

### Latent Class and Profile Analysis:

-   LCA with random effects can be performed with the package
    `r pkg("randomLCA")`. In addition, the package
    `r pkg("e1071")` provides the function `lca()`. Another
    package is `r pkg("poLCA")` for polytomous variable
    latent class analysis. LCA can also be fitted using
    `r pkg("flexmix")` which optionally allows for the
    inclusion of concomitant variables and latent class regression.
-   `r pkg("LCAvarsel")` implements variable selection for
    LCA.
-   `r pkg("tidyLPA")` is a user-friendly implementation of
    latent profile analysis.
-   `r pkg("ClustVarLV")` clusters variables around latent
    variables.
-   `r pkg("multilevLCA")` for single-level and multilevel latent class analysis with covariates.

### Paired Comparisons, Rankings, Ratings:

-   Bradley-Terry models for paired comparisons are implemented in the
    package `r pkg("BradleyTerry2")` and in
    `r pkg("eba")`. The latter allows for the computation of
    elimination-by-aspects models.
-   Recursive partitioning trees for Bradley-Terry models are
    implemented in `r pkg("psychotree", priority = "core")`.
-   `r pkg("BTLLasso")` allows one to include
    subject-specific and object-specific covariates into paired
    comparison models shrinks the effects using Lasso.
-   `r pkg("prefmod", priority = "core")` fits loglinear
    Bradley-Terry models (LLBT) and pattern models for paired
    comparisons, rankings, and ratings.
-   `r github("jvparidon/lmerMultiMember")` can fit the loglinear
    Bradley-Terry model as a mixed-effects model (GLMM) using `lme4`.
-   `r pkg("pcFactorStan")` provides convenience functions
    and pre-programmed Stan models related to the pairwise comparison
    factor model.
-   `r pkg("PLMIX")` fits finite mixtures of Plackett-Luce
    models for partial top rankings/orderings within the Bayesian
    framework.
-   A variety of unfolding techniques for rankings and ratings are
    implemented in `r pkg("smacof")`.
-   Thurstonian IRT models for forced-choice items can be fitted with
    `r pkg("thurstonianIRT")`.

### Network Psychometrics:

-   Estimation of a sparse inverse covariance matrix using a lasso
    penalty (graphical lasso) can be achieved using
    `r pkg("glasso")`.
-   `r pkg("networktools")` includes assorted tools for
    network analysis (bridge centrality, impact, and goldbricker).
-   Bootstrap methods to assess accuracy and stability of estimated
    network structures and centrality indices are implemented in
    `r pkg("bootnet")`.
-   Permutation tests for network comparisons are implemented in
    `r pkg("NetworkComparisonTest")`.
-   Model-based recursive partitioning for networks:
    `r pkg("networktree")`.
-   Network structures for multilevel and graphical vector
    autoregression models can be obtain using
    `r pkg("mlVAR")` and
    `r pkg("graphicalVAR")`.
-   `r pkg("mgm")` estimates time-varying k-order mixed
    graphical models.
-   `r pkg("EstimateGroupNetwork")` can be used to
    simultaneously estimate networks from different groups or classes
    via joint graphical lasso.
-   Various implementations for Ising models:
    `r pkg("IsingSampler")`,
    `r pkg("IsingFit")`.
-   `r pkg("lvnet")` simultaneously estimates factor and
    network models.
-   Approaches for SEM and Confirmatory Network Analysis are implemented
    in `r pkg("psychonetrics")`. This includes multi-group
    (dynamic) SEM in combination with confirmatory network models from
    cross-sectional, time-series and panel data.
-   Network models for longitudinal data estimated within an SEM
    framework: `r pkg("gimme")`.
-   The `r pkg("qgraph")` package can be used to visualize
    data as networks.
-   `r pkg("NetworkToolbox")` implements network analysis
    and graph theory measures used in neuroscience, cognitive science,
    and psychology. Methods include various filtering methods and
    approaches such as threshold, dependency, information filtering
    networks, and efficiency-cost optimization.
-   `r pkg("NetworkToolbox")` implements methods and measures for brain, cognitive, and psychometric network analysis.
-   `r pkg("EGAnet")` implements the Exploratory Graph Analysis (EGA) framework for dimensionality and psychometric assessment.

### Bayesian Psychometrics:

-   `r pkg("blavaan", priority = "core")` fits a variety of
    Bayesian latent variable models, including confirmatory factor
    analysis, structural equation models, and latent growth curve
    models.
-   An analytical framework for latent variables with different Bayesian
    learning methods, including the partially confirmatory factor
    analysis and partially confirmatory IRT is implemented in
    `r pkg("LAWBL")`.
-   Bayesian approaches for estimating item and person parameters by
    means of Gibbs-Sampling are included in
    `r pkg("MCMCpack")`. In addition, the
    `r pkg("pscl")` package allows for Bayesian IRT and roll
    call analysis.
-   `r pkg("LNIRT")` is a package for log-normal response
    time IRT modeling for responses and response times, estimated
    with MCMC.
-   `r pkg("edstan")` provides convenience functions and
    preprogrammed Stan models related to IRT.
-   `r pkg("fourPNO")` can be used for Bayesian 4-PL IRT
    estimation.
-   Gibbs sampling for Bayesian estimation of (Exploratory) Reduced
    Reparameterized Unified Models are implemented in
    `r pkg("rrum")` and `r pkg("errum")`.
-   For Bayesian estimation of the (exploratory) DINA (deterministic
    input, noisy and gate) see `r pkg("dina")` and
    `r pkg("edina")`.
-   `r pkg("slcm")` provides an implementation of the exploratory
    Sparse Latent Class Model (SLCM) for Binary Data.
-   `r pkg("ohoegdm")` performs a Bayesian estimation of the
    ordinal exploratory Higher-order General Diagnostic Model (OHOEGDM) 
    for Polytomous Data.
-   Data package containing coded item and q matrices used in various
    psychometric publications: `r pkg("edmdata")`.
-   `r pkg("BayesLCA")` implements Bayesian LCA.

### Other Related Packages:

-   The `r pkg("psychotools")` provides an infrastructure
    for psychometric modeling such as data classes (e.g., for paired
    comparisons) and basic model fitting functions (e.g., for Rasch and
    Bradley-Terry models).
-   `r pkg("quickpsy")` is a package developed to quickly
    fit and plot psychometric functions for multiple conditions.
-   `r pkg("cNORM")` provides methods for generating
    regression based continuous norms. The approach does not rely on
    prior distribution assumptions and is thus non-parametric, but it
    can be combined with Box-Cox power transformations for
    semi-parametrically modelling the data as well.
-   A system for the management, assessment, and psychometric analysis
    of data from educational and psychological tests is implemented in
    `r pkg("dexter")`, with multi-stage test calibration in
    `r pkg("dexterMST")`, and a GUI via `r pkg("dextergui")`. 
-   Psychometric mixture models based on flexmix infrastructure are
    provided by means of the `r pkg("psychomix")` package
    (at the moment Rasch mixture models and Bradley-Terry mixture
    models).
-   The `r pkg("equate")` package contains functions for
    non-IRT equating under both random groups and nonequivalent groups
    with anchor test designs. Mean, linear, equipercentile and
    circle-arc equating are supported, as are methods for univariate and
    bivariate presmoothing of score distributions. Specific equating
    methods currently supported include Tucker, Levine observed score,
    Levine true score, Braun/Holland, frequency estimation, and chained
    equating.
-   Interactive shiny application for analysis of educational tests and
    their items are provided by the
    `r pkg("ShinyItemAnalysis")` package.
-   Coefficients for interrater reliability and agreements can be
    computed with the `r pkg("irr")`.
-   Statistical tools for the analysis of psychophysical data are implemented in `r pkg("psyphy")` and `r pkg("MixedPsy")`. 
-   Functions and example datasets for Fechnerian scaling of discrete
    object sets are provided by `r pkg("fechner")`. It
    computes Fechnerian distances among objects representing subjective
    dissimilarities, and other related information.
-   The `r pkg("mediation")` allows both parametric and
    nonparametric causal mediation analysis. It also allows researchers
    to conduct sensitivity analysis for certain parametric models.
-   Causal mediation analysis using natural effect models can be performed using `r pkg("medflex")`.
-   The package `r pkg("multiplex")` is especially designed
    for social networks with relations at different levels. In this
    sense, the program has effective ways to treat multiple networks
    data sets with routines that combine algebraic structures like the
    partially ordered semigroup with the existing relational bundles
    found in multiple networks. An algebraic approach for two-mode
    networks is made through Galois derivations between families of the
    pair of subsets.
-   Social Relations Analyses for round robin designs are implemented in
    the `r pkg("TripleR")` package. It implements all
    functionality of the SOREMO software, and provides new functions
    like the handling of missing values, significance tests for single
    groups, or the calculation of the self enhancement index.
-   Fitting and testing multinomial processing tree models, a class of
    statistical models for categorical data with latent parameters, can
    be performed using the `r pkg("mpt")` package. These
    parameters are the link probabilities of a tree-like graph and
    represent the cognitive processing steps executed to arrive at
    observable response categories.The `r pkg("MPTinR")`
    package provides a user-friendly way for analysis of multinomial
    processing tree (MPT) models. The `r pkg("TreeBUGS")`
    package provides user-friendly methods to fit Bayesian hierarchical
    MPT models (beta-MPT and latent-trait MPT) and implements
    posterior-predictive checks, summary plots, correlations and
    regressions for person-level MPT parameters.
-   Beta regression for modeling beta-distributed dependent variables,
    e.g., rates and proportions, is available in
    `r pkg("betareg")`.
-   The `r pkg("cocor")` package provides functions to
    compare two correlations based on either dependent or independent
    groups.
-   The `r pkg("profileR")` package provides a set of tools
    that implement profile analysis and cross-validation techniques.
-   The `r pkg("TestScorer")` package provides a GUI for
    entering test items and obtaining raw and transformed scores. The
    results are shown on the console and can be saved to a tabular text
    file for further statistical analysis. The user can define his own
    tests and scoring procedures through a GUI.
-   `r pkg("wCorr")` calculates Pearson, Spearman,
    tetrachoric polychoric, and polyserial correlation coefficients, in
    weighted or unweighted form.
-   The `r pkg("gtheory")` package fits univariate and
    multivariate generalizability theory (G-theory) models.
-   The `r pkg("GDINA")` package estimates various cognitive
    diagnosis models (CDMs) within the generalized deterministic inputs,
    noisy and gate (G-DINA) model and the sequential G-DINA model
    framework. It can also be used to conduct Q-matrix validation, item
    and model fit statistics, model comparison at the test and item
    level and differential item functioning. A graphical user interface
    is also provided.
-   Simulation routines for cognitive diagnostic model DINA and rRUM are
    implemented in `r pkg("simcdm")`.
-   `r pkg("TestDataImputation")` for missing item responses
    imputation for test and assessment data.
-   `r pkg("LAM")` includes some procedures for latent
    variable modeling with a particular focus on multilevel data.
-   `r pkg("psychTools")` contains tools to accompany the
    `r pkg("psych")` package.
-   `r pkg("ata")` provides a collection of psychometric
    methods to process item metadataand use target assessment and
    measurement blueprint constraints to assemble a test form.
-   `r pkg("heplots")`: Visualizing hypothesis tests in multivariate linear models with hypothesis error plots.
-   `r pkg("simlandr")` provides tools to estimate generalized
    potential landscapes from formal psychological models.
    `r pkg("fitlandr")` provides nonparametric methods to estimate 
    bivariate vector fields and generalized potential landscapes from 
    intensive longitudinal data. `r pkg("Isinglandr")` can be used to
    estimate generalized potential landscapes from Ising networks 
    generated from the `r pkg("IsingFit")` package.


### Links
-   [Journal of Statistical Software Special Volume: Psychometrics in R](http://www.jstatsoft.org/v20/)
-   [Journal of Statistical Software Special Volume: Psychoco - Psychometric Computing in R](http://www.jstatsoft.org/v48/)
-   [Journal of Statistical Software Special Volume: Festschrift for Jan de Leeuw](http://www.jstatsoft.org/v73/)
-   [Tools for computational cognitive science](https://github.com/djnavarro)
