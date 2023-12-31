# woolmixmod: Untangling the Zooarchaeological Evidence for Intensive Wool Production in Northern Mesopotamia, c. 4500-1500 cal. BC

Contact information: Max D. Price (max.d.price@durham.ac.uk) and Jesse Wolfhagen (jwolfhag@purdue.edu)

This GitHub repository maintains the code and data files necessary to reproduce the analyses of the manuscript "Wool They, Won’t They: Untangling the Zooarchaeological Evidence for Intensive Wool Production in Northern Mesopotamia c. 4500-1500 cal. BC”, authored by Max D. Price and Jesse Wolfhagen, submitted to the *Journal of Anthropological Archaeology*. The methods rely on [Stan](https://mc-stan.org/) for Bayesian analysis and are written for use in the R statistical computing framework.

Bayesian modeling results are saved as three files for the metric data and two files for the aging data. These files contain the posterior distributions of key model parameters at different scales of inference (the site-level and period-level). Additionally, the metric analysis contains a "specimen-level" result that displays the posterior mean membership probability for each specimen (i.e., the calculated probability that a specimen is immature, female, or male).

References for the data sources are available in the associated bibtex file, as referenced in the two data files (wool_mandible_data.csv and woolmixmod_demographic_data.csv).

**NOTE**: Currently, the code and data *do not* include measurement data from three sites that are included in the paper's analysis (Jebel Aruda, Tell es-Sweyhat, and El Qitar).

# Abstract	

Wool was of central importance in ancient Mesopotamia, a fact amply demonstrated by texts documenting the large-scale production, processing, and exchange of wool. It has long been thought that, at some point in the Chalcolithic or Bronze Age, production of caprine husbandry was reorganized to support the burgeoning wool industry. Zooarchaeologists have devised a number of methodological techniques to examine this hypothetical process. Here, we offer a critical examination of zooarchaeological kill-off patterns and biometrical data from sites dating to the Chalcolithic through Middle Bronze Age (c. 4500-1500 BC) in northern Mesopotamia/eastern Anatolia. We synthesize the existing data and use Bayesian statistical modeling to estimate three indicators of intensive wool production: the retention of older animals, the ratio of males to females, and the size of sheep. The data provide no definitive answer. While we confirm previous assessments that show an increase in sheep size in the 4th millennium BC, we find no clear pattern in the ratio of adult males to adult females. Kill-off data show remarkably subtle shifts in the slaughter of older caprines at a regional level rather than a dramatic shift toward intensive wool production. However, an in-depth analysis of an urban network in the Karababa Basin of SE Turkey in the 3rd millennium BC is consistent with an economic focus on wool production.
