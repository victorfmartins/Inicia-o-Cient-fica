# Scientific Initiation
This repository was made to share my ongoing works in neural signal processing analysis

### What is this project?
It is part of my graduation requirements under the [Molecular Science Course](https://www.cecm.usp.br/#) at [University of São Paulo](https://www5.usp.br/).
I have chosen to work in the analysis of neuro time series from electrode recording in the mouse brain. I work with a team from the Federal University of Rio Grande do Norte in electro-neurophysiology analysis. My work is to study the effect of sample size in several standard metrics of the literature. I read papers and reverse engineer its analysis looking for sample size bias. The goal is to introduce new metric that indicates the extent of sample size bias in each one. [Lucas Tavares](https://github.com/lucaase) and [Adriano Tort](https://scholar.google.com.br/citations?hl=en&user=Z7lq_2gAAAAJ&view_op=list_works&sortby=pubdate) supervises my studies.

### Contents
There are five scripts name coded as <number_metricName>. Each one presents a toy model and its analysis over sample size bias of a given metric. Some of them are followed by analysis of real data recorded by the [Computational Neurophysiology Laboratory](https://tortlab.github.io/) team or taken from public online repositories such as [CRCNS](https://crcns.org/).

SimulationResults folder contains saved data from length computations to prevent time consuming recalculations. Each subfolder received the name of the script that generated the data.

DATA.mat contains the real data recorded by the Labs and used in the five main scripts.

MVGC, CircStat, boudedline, randraw.m, and eegfilt.m are toolboxes and scripts from third parties gathere from [MathWorks](https://www.mathworks.com/) website

### 📝 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

Made by Victor Martins [See my linkedin](https://www.linkedin.com/in/victor-franco-martins-1503a417b)