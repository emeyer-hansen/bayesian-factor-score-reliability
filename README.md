[_metadata_:author]:- "Emil Niclas Meyer-Hansen"
[_metadata_:date]:- "2026-11-26"
[_metadata_:tags]:- "markdown metadata"
# Bayesian Factor Score Reliability: Conceptualizing a Bayesian Measure of Factor Score Indeterminacy

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24552942.svg)](https://doi.org/10.5281/zenodo.24552942)
[![Release](https://img.shields.io/github/v/release/emeyer-hansen/bayesian-factor-score-reliability)](https://github.com/emeyer-hansen/bayesian-factor-score-reliability/releases/latest)

*Emil Niclas Meyer-Hansen*

2026-11-26

A formal Bayesian (re)conceptualization of the factor score reliability (FSR) index is developed by the author for the social sciences in the research paper *Bayesian Factor Score Reliability: Conceptualizing a Bayesian Measure of Factor Score Indeterminacy* (see the Abstract below). The paper is currently a work-in-progress. It subscribes to the [*open science standard*](https://www.cos.io/open-science), is partly licensed under the CC BY 4.0 and partly under the GPL-3.0 (or later) and GPL-2.0 (see the License below), and it is made freely available in PDF-format.


## Table of Contents
- [Abstract](#abstract)
- [FAQ](#faq)
- [Changelog](#changelog)
- [Copyright & License](#license)
- [Citation](#citation)
- [References](#references)


<a id="abstract"></a>
## Abstract

The indeterminacy of scores is a fundamental problem in factor analysis. Since a high indeterminacy makes results based on factor scores unreliable, it is necessary for researchers to demonstrate that the degree of indeterminacy is low. Within the frequentist statistical framework, this can be assessed with the factor score reliability (FSR) index, but for the increasingly popular Bayesian statistical framework, no corresponding measure currently exists. The lack of a Bayesian version prevents researchers from using the comparative advantages of that framework, such as being able to incorporate prior information into the FSR index and formulate coherent probabilistic statements from its posterior distribution. This paper solves that limitation by providing a novel Bayesian (re)conceptualization of this measure of factor score indeterminacy and free software for its implementation. The mathematical properties of this Bayesian FSR (BSFR) index is compared to its frequentist counterpart, and its ease of application is demonstrated in relation to the honesty-{\allowbreak}humility personality trait using empirical data.[©](#license)

**Keywords**: *Factor score reliability; Factor score (in)determinacy; Factor scores; Factor analysis; Structural equation modeling; Covariance structure analysis; Multivariate analysis; Incomplete data analysis; Bayesian inference; Frequentist inference; Honesty-humility; IPIP-HEXACO*


<a id="faq"></a>
## FAQ
- Where can I download the latest version of the paper? The newest version of the paper is always provided as the file named 'Article.pdf'.
- How can I reproduce the findings of this paper? The findings of the research paper can be reproduced in R using the files provided in the 'Archive'-folder. Identify the version of the paper that you want to replicate, download the files, specify your working directory and begin executing the code provided in the 'analysis.R'-file. Due to Github-limitations on file size, there are restrictions on what data and plots can be provided here.
- Can I download these materials elsewhere? No, all materials are currently only provided here on Github.


<a id="changelog"></a>
## Changelog
- **2026-11-26-08 8:12 CEST**
    - [Version 2026-11-26-08] - Working Paper (Initial Release).


<a id="license"></a>
## Copyright & License (Addendum)

\textbf{Copyright}: 2026 The author. \textbf{License}: Except where otherwise indicated, all contents of this document and associated files are licensed under the *Creative Commons Attribution 4.0 International License* ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)). All data and software, including but *not* necessarily limited to, source code, executable code, code snippets, code chunks, algorithms, and/or scripts, attributable to this document and/or its associated files are expressly excluded from the foregoing license, and unless otherwise indicated, are instead licensed under the *GNU General Public License, version 3* ([GPL-3.0](https://www.gnu.org/licenses/gpl-3.0.html)) but permission is granted to re-license these materials under later versions. The processed data sampled from the Open-Source Psychometrics Project (2014; provided in Petras & Dantlgraber, 2022) is excluded from the foregoing licenses and is instead licensed under the *GNU General Public License, version 2* ([GPL-2.0](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)). By engaging with this document and/or any associated files, which include, but are *not* necessarily limited to, downloading, using, viewing, and/or distributing any of them, in parts of whole, you agree to comply with the applicable license terms for the respective content types.


<a id="citation"></a>
## Citation

Building on Bayesian inference and similar conceptualizations made for other scientific disciplines, this formal conceptualization of Bayesian factor score reliability (BFSR) index is an original (re)conceptualization by Emil Niclas Meyer-Hansen, conceived as part of the research paper associated with this project. For correspondence, contact the author via email: [emil098meyerhansen@gmail.com](mailto:emil098meyerhansen@gmail.com)

Please, if you use, refer to, modify, and/or continue the development of this BFSR index, provide proper reference and citation to its founding author. An example of proper citation is provided below:
```
Meyer-Hansen, E. N. (2026): 'Bayesian Factor Score Reliability: Conceptualizing a Bayesian Measure of Factor Score Indeterminacy', Github, Working paper (v2026-11-26-08).
```

For LaTeX users, a BibTeX entry is provided below:
```
@unpublished{,
  title = {Bayesian Factor Score Reliability: Conceptualizing a Bayesian Measure of Factor Score Indeterminacy},
  author = {Emil Niclas Meyer-Hansen},
  publisher = {Github},
  year = {2026},
  doi = {},
  pubstate = {Working Paper (v2026-11-26-08)}
}
```


<a id="references"></a>
## References
- Open-Source Psychometrics Project. (2014). Answers to the IPIP HEXACO equivalent scales [Dataset]. (Accessed 2026/10/11)
- Petras, N., Dantlgraber, M., & Reips, U.-D. (2024). Illustrating psychometric tests, scales, and constructs: An R package for item pool visualization. Behavior Research Methods, 56(2), 639-650. doi: 10.3758/s13428-022-02052-7
