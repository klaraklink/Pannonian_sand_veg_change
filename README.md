# Klinkovská et al. (2025) Preslia: Vegetation remains, specialists fade: changes in Pannonian sand grasslands

## Supplementary code and data to the article:

Klinkovská, K., Axmanová, I., Borovyk, D., Ćuk, M., Kopčan, B., Valachovič, M., Clark, A. T. & Chytrý, M. (2025) Vegetation remains, specialists fade: changes in Pannonian sand grasslands. Preslia.

### Data

The data contains vascular plant species composition data from resurveyed vegetation plots of sand vegetation in the Vienna Basin (Czech Republic, Slovakia, Austria, 48°12'–48°58'N, 16°39'–17°19'E). The dataset consists of 86 vegetation plots first surveyed between 1952 and 2010 (with the most plots sampled in 1994–1997 and 2001–2006) and resurveyed in 2023–2024. The historical plots were relocated using the geographical coordinates where available, original locality descriptions, and information on slope, aspect and elevation. We measured the geographical coordinates of each plot using GPS with a location uncertainty of 2–5 m. In 2024, most of the plots were also located using a differential GPS with an accuracy of ca. 5 cm in each corner of the plot. Most plots were squares of 16 m^2^ in the Czech Republic, 25 m^2^ in Slovakia and 100 m^2^ in Austria. During the repeated sampling, we always used the same plot size as in the first sampling.

The percentage cover of each vascular plant species was estimated. The methodology of cover estimation differed across different sources. The cover was either estimated using a cover scale, usually the seven-grade or the nine-grade Braun-Blanquet scale (Westhoff & van der Maarel 1978), or in percentages. Prior to the analysis, different cover scales were transformed into the seven-grade Braun-Blanquet scale, which uses the lowest number of cover classes. Then, we converted the cover classes to percentages as the mean value of the given cover class interval.

The header data structure follows that of the ReSurveyEurope Database (<http://euroveg.org/eva-database-re-survey-europe>).

The data on species composition, environmental variables and species characteristics are provided as csv files with columns separated by commas or excel files:

-   `pannonian_sands_resurvey_spe.csv` contains the percentage covers of vascular plant species in the resurveyed plots. Taxonomic concepts and nomenclature follow Kaplan et al. (2019), or Euro+Med PlantBase (Euro+Med 2025) in the case of taxa not occurring in the Czech Republic.
-   `pannonian_sands_resurvey_head.csv` contains header data for the vegetation plots.
-   `annuals.csv` contains the list of annual species which were excluded from the analysis due to the possible phenological differences during the first and the second surveys.
-   `diag.csv` contains the list of species diagnostic for sand vegetation, dry grasslands and ruderal vegetation according to the national vegetation classifications (Mucina et al. 1993, Valachovič et al. 1995, Chytrý 2007) and corresponding habitat types according to Chytrý et al. (2010) and Šuvada (2023).
-   `red_list.csv` contains the list of species included in the Red Lists of vascular plants of the Czech Republic (Grulich 2017), Slovakia (Eliáš 2015) and Austria (Schratt-Ehrendorfer et al. 2022).
-   `aliens.xlsx` contains the list of alien species according to the national checklists of alien plants (Medvecká et al. 2012, Pyšek et al. 2020, Glaser et al. In preparation).
-   `woody.cscv` contains the list of trees and shrubs according to Dřevojan (2020).
-   `Ellenberg_disturbance.xlsx` contains the data on species indicator values used in the analysis (Tichý et al. 2023, Midolo et al. 2023) and `taxonomy_eu_convert.csv` is a conversion table used to match our nomenclature with that of the species indicator values data.

The data on resurveyed vegetation plots are also stored in the ReSurveyEurope database (CZ_0039, Knollová et al. 2024; <http://euroveg.org/eva-database-re-survey-europe>).

\### Scripts

-   `1_classification.R`: Classification of the historical vegetation plots into three vegetation types (pioneer sand vegetation, closed acidophilous and ruderalized sand grasslands and closed basiphilous sand grasslands)
-   `2_ordinations.R`: Changes in species composition using PCoA and db-RDA.
-   `3_univariate_models`: Changes in species richness, proportions of sand and dry grassland specialists, ruderal, threatened, alien and woody species per plot.
-   `4_eivs_permutation_test.R`: Changes in Ellenberg-type indicator values and disturbance indicator values.
