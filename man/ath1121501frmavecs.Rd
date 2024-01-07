\name{ath1121501frmavecs}
\docType{data}
\alias{ath1121501frmavecs}
\title{Vectors used by fRMA for ath1121501}
\description{Vectors allowing to apply the fRMA procedure on data obtained with Affymetrix Arabidopsis ATH1 Genome Array and annotated with the ATH1-121501 CDF file (GPL198 platform). This data object was automatically created by the package frmaTools version 1.52.0.}
\usage{data(ath1121501frmavecs)}
\format{A list with 6 elements.
  \tabular{ll}{
    normVec \tab normalization vector \cr
    probeVec \tab probe effect vector \cr
    probeVarWithin \tab within batch probe variance \cr
    probeVarBetween \tab between batch probe variance \cr
    probesetSD \tab within probeset standard deviation \cr
    medianSE \tab median standard errors \cr
}}
\examples{
 data(ath1121501frmavecs) 
 str(ath1121501frmavecs)
}
\details{
The vectors were computed based on 100 triplicates, orginating from 100 different data series available on the GPL198 platform on Gene Expression omnibus (GEO), as recommended by the authors of the frozen Robust Multiarray Analysis (McCall et al., 2010 - DOI: 10.1093/biostatistics/kxp059; McCall et Irizarry, 2011 - DOI: 10.1186/1471-2105-12-369). Those 100 triplicates were related to a large variety of study subjects, encompassing biotic treatments, abiotic stresses, control conditions (in seedlings), development (different tissues, organs, in control condition) and chemical treatments (hormones, growth regulators, ...). The following table provides the GEO sample ID corresponding to those triplicates, alongside with the GEO sample title, GEO series ID, the class (biotic, abiotic, control, development, chemical), genotype, biomaterial (tissues, organs, cells,..), age and treatment. All this information was manually curated.
  \tabular{llllllll}{
    sample \tab title \tab series \tab class \tab genotype \tab biomaterial \tab age \tab treatment \cr
    GSM508439 \tab pathogen infection: TuMV inoculated - RNA fraction: total RNA - rep1 \tab GSE20278 \tab biotic \tab Col-0 \tab rosette leaf \tab 3.5-week-old \tab Turnip mosaic virus \cr
    GSM508440 \tab pathogen infection: TuMV inoculated - RNA fraction: total RNA - rep2 \tab GSE20278 \tab biotic \tab Col-0 \tab rosette leaf \tab 3.5-week-old \tab Turnip mosaic virus \cr
    GSM508441 \tab pathogen infection: TuMV inoculated - RNA fraction: total RNA - rep3 \tab GSE20278 \tab biotic \tab Col-0 \tab rosette leaf \tab 3.5-week-old \tab Turnip mosaic virus \cr
    GSM39203 \tab Col_Chitin1 \tab GSE2169 \tab biotic \tab Col-0 \tab seedling \tab 2-week-old \tab chitin \cr
    GSM39204 \tab Col_Chitin2 \tab GSE2169 \tab biotic \tab Col-0 \tab seedling \tab 2-week-old \tab chitin \cr
    GSM39205 \tab Col_Chitin3 \tab GSE2169 \tab biotic \tab Col-0 \tab seedling \tab 2-week-old \tab chitin \cr
    GSM48125 \tab Col_CSC1 \tab GSE2538 \tab biotic \tab Col-0 \tab seedling \tab 2-week-old \tab crab shell chitin \cr
    GSM48126 \tab Col_CSC2 \tab GSE2538 \tab biotic \tab Col-0 \tab seedling \tab 2-week-old \tab crab shell chitin \cr
    GSM48127 \tab Col_CSC3 \tab GSE2538 \tab biotic \tab Col-0 \tab seedling \tab 2-week-old \tab crab shell chitin \cr
    GSM133079 \tab JD AT+EO COL WT 12H INFECTED \tab GSE5686 \tab biotic \tab Col-0 \tab rosette leaf \tab 29-day-old \tab Erysiphe orontii \cr
    GSM133095 \tab JD AT+EO COL WT EXP2 12H INFECTED \tab GSE5686 \tab biotic \tab Col-0 \tab rosette leaf \tab 29-day-old \tab Erysiphe orontii \cr
    GSM133101 \tab JD AT+EO TIME EXP3 EO INF 12H \tab GSE5686 \tab biotic \tab Col-0 \tab rosette leaf \tab 29-day-old \tab Erysiphe orontii \cr
    GSM1111986 \tab AB Infected Col-0 leaf , biological replicate 1 \tab GSE45690 \tab biotic \tab Col-0 \tab leaf \tab 4-week-old \tab Alternaria brassicicola \cr
    GSM1111987 \tab AB Infected Col-0 leaf , biological replicate 2 \tab GSE45690 \tab biotic \tab Col-0 \tab leaf \tab 4-week-old \tab Alternaria brassicicola \cr
    GSM1111988 \tab AB Infected Col-0 leaf , biological replicate 3 \tab GSE45690 \tab biotic \tab Col-0 \tab leaf \tab 4-week-old \tab Alternaria brassicicola \cr
    GSM1144681 \tab Plant inoculated with PsJN inactivated, biological rep 1 \tab GSE47092 \tab biotic \tab Col-0 \tab seedling \tab 13-day-old \tab Burkholderia phytofirmans inactivated|(biotic treatment) \cr
    GSM1144682 \tab Plant inoculated with PsJN inactivated, biological rep 2 \tab GSE47092 \tab biotic \tab Col-0 \tab seedling \tab 13-day-old \tab Burkholderia phytofirmans inactivated|(biotic treatment) \cr
    GSM1144683 \tab Plant inoculated with PsJN inactivated, biological rep 3 \tab GSE47092 \tab biotic \tab Col-0 \tab seedling \tab 13-day-old \tab Burkholderia phytofirmans inactivated|(biotic treatment) \cr
    GSM1521091 \tab Full AC2 at T1, biological rep1 \tab GSE62180 \tab biotic \tab Col-0 \tab rosette leaf \tab 4-week-old \tab Cabbage leaf curl virus full AC2 protein \cr
    GSM1521092 \tab Full AC2 at T1, biological rep2 \tab GSE62180 \tab biotic \tab Col-0 \tab rosette leaf \tab 4-week-old \tab Cabbage leaf curl virus full AC2 protein \cr
    GSM1521093 \tab Full AC2 at T1, biological rep3 \tab GSE62180 \tab biotic \tab Col-0 \tab rosette leaf \tab 4-week-old \tab Cabbage leaf curl virus full AC2 protein \cr
    GSM1532918 \tab Cocultivation at 8 hrs, biological rep1_leaf \tab GSE62749 \tab biotic \tab Col-0 \tab leaf \tab 3-week-old \tab Agrobacterium tumefaciens C58 \cr
    GSM1532919 \tab Cocultivation at 8 hrs, biological rep2_leaf \tab GSE62749 \tab biotic \tab Col-0 \tab leaf \tab 3-week-old \tab Agrobacterium tumefaciens C58 \cr
    GSM1532922 \tab Cocultivation at 8 hrs, biological rep1_root \tab GSE62749 \tab biotic \tab Col-0 \tab leaf \tab 3-week-old \tab Agrobacterium tumefaciens C58 \cr
    GSM1970346 \tab Col-0 roots, Fusarium treatment, replicate 2 \tab GSE75928 \tab biotic \tab Col-0 \tab root \tab 4-week-old \tab Fusarium oxysporum \cr
    GSM1970347 \tab Col-0 roots, Fusarium treatment, replicate 3 \tab GSE75928 \tab biotic \tab Col-0 \tab root \tab 4-week-old \tab Fusarium oxysporum \cr
    GSM1970348 \tab Col-0 roots, Fusarium treatment, replicate 4 \tab GSE75928 \tab biotic \tab Col-0 \tab root \tab 4 week old \tab Fusarium oxysporum \cr
    GSM2325661 \tab seedling roots, 3h OGs, bio rep 1 \tab GSE87216 \tab biotic \tab Col-0 \tab root \tab 1-week-old \tab pectin-derived oligogalacturonides \cr
    GSM2325662 \tab seedling roots, 3h OGs, bio rep 2 \tab GSE87216 \tab biotic \tab Col-0 \tab root \tab 1-week-old \tab pectin-derived oligogalacturonides \cr
    GSM2325663 \tab seedling roots, 3h OGs, bio rep 3 \tab GSE87216 \tab biotic \tab Col-0 \tab root \tab 1-week-old \tab pectin-derived oligogalacturonides \cr
    GSM151700 \tab Col-0 48 hpi, biological replicate 1 \tab GSE6556 \tab biotic \tab Col-0 \tab leaf \tab 5-week-old \tab Pseudomonas syringae pv. tomato DC3000 (avrRpt2) \cr
    GSM151701 \tab Col-0 48 hpi, biological replicate 2 \tab GSE6556 \tab biotic \tab Col-0 \tab leaf \tab 5-week-old \tab Pseudomonas syringae pv. tomato DC3000 (avrRpt2) \cr
    GSM151702 \tab Col-0 48 hpi, biological replicate 3 \tab GSE6556 \tab biotic \tab Col-0 \tab leaf \tab 5-week-old \tab Pseudomonas syringae pv. tomato DC3000 (avrRpt2) \cr
    GSM828864 \tab Leafminer damaged A.thaliana seedlings rep1 \tab GSE33505 \tab biotic \tab Col-0 \tab leaf \tab 4-week-old \tab Liriomyza huidobrensis \cr
    GSM828865 \tab Leafminer damaged A.thaliana seedlings rep2 \tab GSE33505 \tab biotic \tab Col-0 \tab leaf \tab 4-week-old \tab Liriomyza huidobrensis \cr
    GSM828866 \tab Leafminer damaged A.thaliana seedlings rep3 \tab GSE33505 \tab biotic \tab Col-0 \tab leaf \tab 4-week-old \tab Liriomyza huidobrensis \cr
    GSM298415 \tab A.thaliana Col-0, Flg22 4h rep1 \tab GSE11807 \tab biotic \tab Col-0 \tab leaf \tab 6-week-old \tab flagelin22 \cr
    GSM298418 \tab A.thaliana Col-0, Flg22 4h rep3 \tab GSE11807 \tab biotic \tab Col-0 \tab leaf \tab 6-week-old \tab flagelin22 \cr
    GSM298420 \tab A.thaliana Col-0, Flg22 4h rep2 \tab GSE11807 \tab biotic \tab Col-0 \tab leaf \tab 6-week-old \tab flagelin22 \cr
    GSM3449511 \tab Col-0/Psm-1, biological rep1 \tab GSE121886 \tab biotic \tab Col-0 \tab leaf \tab 4-week-old \tab Pseudomonas syringae pv. maculicola ES4326 \cr
    GSM3449512 \tab Col-0/Psm-2, biological rep2 \tab GSE121886 \tab biotic \tab Col-0 \tab leaf \tab 4-week-old \tab Pseudomonas syringae pv. maculicola ES4326 \cr
    GSM3449513 \tab Col-0/Psm-3, biological rep3 \tab GSE121886 \tab biotic \tab Col-0 \tab leaf \tab 4-week-old \tab Pseudomonas syringae pv. maculicola ES4326 \cr
    GSM436317 \tab Col PSTDC3000 24h 1st rep \tab GSE17500 \tab biotic \tab Col-0 \tab leaf \tab 5-week-old \tab Pseudomonas syringae pv. tomato DC3000 \cr
    GSM436318 \tab Col PSTDC3000 24h 2nd rep \tab GSE17500 \tab biotic \tab Col-0 \tab leaf \tab 5-week-old \tab Pseudomonas syringae pv. tomato DC3000 \cr
    GSM436319 \tab Col PSTDC3000 24h 3rd rep \tab GSE17500 \tab biotic \tab Col-0 \tab leaf \tab 5-week-old \tab Pseudomonas syringae pv. tomato DC3000 \cr
    GSM469414 \tab DeVos_1-2_Treatment_Rep1_ATH1 \tab GSE18960 \tab biotic \tab Col-0 \tab leaf \tab 3-week-old \tab Myzus persicae saliva \cr
    GSM469416 \tab DeVos_1-4_Treatment_Rep2_ATH1 \tab GSE18960 \tab biotic \tab Col-0 \tab leaf \tab 3-week-old \tab Myzus persicae saliva \cr
    GSM469418 \tab DeVos_1-6_Treatment_Rep3_ATH1 \tab GSE18960 \tab biotic \tab Col-0 \tab leaf \tab 3-week-old \tab Myzus persicae saliva \cr
    GSM469763 \tab Mitra_2-7_Col0-PsmES4326-24hpi_Rep1_ATH1 \tab GSE18978 \tab biotic \tab Col-0 \tab leaf \tab 4-week-old \tab Pseudomonas syringae pv. maculicola ES4326 \cr
    GSM469764 \tab Mitra_2-8_Col0-PsmES4326-24hpi_Rep2_ATH1 \tab GSE18978 \tab biotic \tab Col-0 \tab leaf \tab 4-week-old \tab Pseudomonas syringae pv. maculicola ES4326 \cr
    GSM469765 \tab Mitra_2-9_Col0-PsmES4326-24hpi_Rep3_ATH1 \tab GSE18978 \tab biotic \tab Col-0 \tab leaf \tab 4-week-old \tab Pseudomonas syringae pv. maculicola ES4326 \cr
    GSM545384 \tab colptobis_2 with array type ATH1 from Affymetrix \tab GSE21920 \tab biotic \tab Col-0 \tab rosette leaf \tab 4-week-old \tab Pseudomonas syringae pv. tomato  \cr
    GSM545385 \tab colptobis_1 with array type ATH1 from Affymetrix \tab GSE21920 \tab biotic \tab Col-0 \tab rosette leaf \tab 4-week-old \tab Pseudomonas syringae pv. tomato  \cr
    GSM545386 \tab colptobis_3 with array type ATH1 from Affymetrix \tab GSE21920 \tab biotic \tab Col-0 \tab rosette leaf \tab 4-week-old \tab Pseudomonas syringae pv. tomato  \cr
    GSM605343 \tab Pseudomonas-treatment 1 \tab GSE24552 \tab biotic \tab Col-0 \tab seedling \tab 18-day-old \tab Pseudomonas sp. G62 \cr
    GSM605344 \tab Pseudomonas-treatment 2 \tab GSE24552 \tab biotic \tab Col-0 \tab seedling \tab 18-day-old \tab Pseudomonas sp. G62 \cr
    GSM605345 \tab Pseudomonas-treatment 3 \tab GSE24552 \tab biotic \tab Col-0 \tab seedling \tab 18-day-old \tab Pseudomonas sp. G62 \cr
    GSM921521 \tab Arab_GiantCell_Root_14d_BioRep1 \tab GSE37553 \tab biotic \tab Col-0 \tab giant cell \tab 3-week-old \tab Meloidogyne incognita \cr
    GSM921522 \tab Arab_GiantCell_Root_14d_BioRep2 \tab GSE37553 \tab biotic \tab Col-0 \tab giant cell \tab 3-week-old \tab Meloidogyne incognita \cr
    GSM921523 \tab Arab_GiantCell_Root_14d_BioRep3 \tab GSE37553 \tab biotic \tab Col-0 \tab giant cell \tab 3-week-old \tab Meloidogyne incognita \cr
    GSM5597522 \tab Col 6h wound A \tab GSE184793 \tab abiotic \tab Col-0 \tab leaf \tab  \tab wounding \cr
    GSM5597523 \tab Col 6h wound B \tab GSE184793 \tab abiotic \tab Col-0 \tab leaf \tab  \tab wounding \cr
    GSM5597524 \tab Col 6h wound C \tab GSE184793 \tab abiotic \tab Col-0 \tab leaf \tab  \tab wounding \cr
    GSM2769932 \tab Wild-type seedlings, primed and triggered, biological rep2 \tab GSE103398 \tab abiotic \tab Col-0 \tab seedling \tab 1-week-old \tab heat primed and triggered (37°C 1 hour + normal temperature 3 days + 37°C 1 hour) \cr
    GSM2769933 \tab Wild-type seedlings, primed and triggered, biological rep3 \tab GSE103398 \tab abiotic \tab Col-0 \tab seedling \tab 1-week-old \tab heat primed and triggered (37°C 1 hour + normal temperature 3 days + 37°C 1 hour) \cr
    GSM2769931 \tab Wild-type seedlings, primed and triggered, biological rep1 \tab GSE103398 \tab abiotic \tab Col-0 \tab seedling \tab 1-week-old \tab heat primed and triggered (37°C 1 hour + normal temperature 3 days + 37°C 1 hour) \cr
    GSM265411 \tab Arabidopsis, whole roots, -Fe, replicate 1 \tab GSE10496 \tab abiotic \tab Col-0 \tab root \tab 6-day-old \tab iron deficiency \cr
    GSM265412 \tab Arabidopsis, whole roots, -Fe, replicate 2 \tab GSE10496 \tab abiotic \tab Col-0 \tab root \tab 6-day-old \tab iron deficiency \cr
    GSM265413 \tab Arabidopsis, whole roots, -Fe, replicate 3 \tab GSE10496 \tab abiotic \tab Col-0 \tab root \tab 6-day-old \tab iron deficiency \cr
    GSM2340100 \tab Low blue light, rep1 \tab GSE87770 \tab abiotic \tab Col-0 \tab petiole \tab 29-day-old \tab low blue light \cr
    GSM2340101 \tab Low blue light, rep2 \tab GSE87770 \tab abiotic \tab Col-0 \tab petiole \tab 29-day-old \tab low blue light \cr
    GSM2340102 \tab Low blue light, rep3 \tab GSE87770 \tab abiotic \tab Col-0 \tab petiole \tab 29-day-old \tab low blue light \cr
    GSM1444203 \tab Col upon hypoxia under light submergence for 48 h, biological rep1 \tab GSE59719 \tab abiotic \tab Col-0 \tab leaf \tab 4-week-old \tab submergence (hypoxia) \cr
    GSM1444204 \tab Col upon hypoxia under light submergence for 48 h, biological rep2 \tab GSE59719 \tab abiotic \tab Col-0 \tab leaf \tab 4-week-old \tab submergence (hypoxia) \cr
    GSM1444205 \tab Col upon hypoxia under light submergence for 48 h, biological rep3 \tab GSE59719 \tab abiotic \tab Col-0 \tab leaf \tab 4-week-old \tab submergence (hypoxia) \cr
    GSM290628 \tab rosettes_drought_rep1 \tab GSE11538 \tab abiotic \tab Col-0 \tab rosette \tab 4-week-old \tab drought \cr
    GSM290629 \tab rosettes_drought_rep2 \tab GSE11538 \tab abiotic \tab Col-0 \tab rosette \tab 4-week-old \tab drought \cr
    GSM290630 \tab rosettes_drought_rep3 \tab GSE11538 \tab abiotic \tab Col-0 \tab rosette \tab 4-week-old \tab drought \cr
    GSM451832 \tab WT_3h_-P_rep1 \tab GSE18071 \tab abiotic \tab Col-0 \tab shoot \tab 2-week-old \tab phosphorus deficiency \cr
    GSM451833 \tab WT_3h_-P_rep2 \tab GSE18071 \tab abiotic \tab Col-0 \tab shoot \tab 2-week-old \tab phosphorus deficiency \cr
    GSM451834 \tab WT_3h_-P_rep3 \tab GSE18071 \tab abiotic \tab Col-0 \tab shoot \tab 2-week-old \tab phosphorus deficiency \cr
    GSM604644 \tab seedling 24h light, (+Lin), rep1 \tab GSE24517 \tab abiotic \tab Col-0 \tab seedling \tab 1-week-old \tab lincomycin \cr
    GSM604645 \tab seedling 24h light, (+Lin), rep2 \tab GSE24517 \tab abiotic \tab Col-0 \tab seedling \tab 1-week-old \tab lincomycin \cr
    GSM604646 \tab seedling 24h light, (+Lin), rep3 \tab GSE24517 \tab abiotic \tab Col-0 \tab seedling \tab 1-week-old \tab lincomycin \cr
    GSM852939 \tab leaf_Col_ozone_2d_biological rep1 \tab GSE34667 \tab abiotic \tab Col-0 \tab leaf \tab 5-week-old \tab ozone \cr
    GSM852940 \tab leaf_Col_ozone_2d_biological rep2 \tab GSE34667 \tab abiotic \tab Col-0 \tab leaf \tab 5-week-old \tab ozone \cr
    GSM852941 \tab leaf_Col_ozone_2d_biological rep3 \tab GSE34667 \tab abiotic \tab Col-0 \tab leaf \tab 5-week-old \tab ozone \cr
    GSM989196 \tab wild type (Col) glucose treatment, biological rep1 \tab GSE40245 \tab abiotic \tab Col-0 \tab seedling \tab 3-day-old \tab glucose \cr
    GSM989197 \tab wild type (Col) glucose treatment, biological rep2 \tab GSE40245 \tab abiotic \tab Col-0 \tab seedling \tab 3-day-old \tab glucose \cr
    GSM989198 \tab wild type (Col) glucose treatment, biological rep3 \tab GSE40245 \tab abiotic \tab Col-0 \tab seedling \tab 3-day-old \tab glucose \cr
    GSM1415488 \tab light-hs, rep1 \tab GSE58616 \tab abiotic \tab Col-0 \tab seedling \tab 13-day-old \tab heat shock \cr
    GSM1415489 \tab light-hs, rep2 \tab GSE58616 \tab abiotic \tab Col-0 \tab seedling \tab 13-day-old \tab heat shock \cr
    GSM1415490 \tab light-hs, rep3 \tab GSE58616 \tab abiotic \tab Col-0 \tab seedling \tab 13-day-old \tab heat shock \cr
    GSM1446768 \tab TiO2 treated germinant, biological rep1 \tab GSE59809 \tab abiotic \tab Col-0 \tab seedling \tab 12-day-old \tab titanium dioxide \cr
    GSM1446769 \tab TiO2 treated germinant, biological rep2 \tab GSE59809 \tab abiotic \tab Col-0 \tab seedling \tab 12-day-old \tab titanium dioxide \cr
    GSM1446770 \tab TiO2 treated germinant, biological rep3 \tab GSE59809 \tab abiotic \tab Col-0 \tab seedling \tab 12-day-old \tab titanium dioxide \cr
    GSM347118 \tab Arabidopsis rosette_Light Wounded_Rep3 \tab GSE13803 \tab abiotic \tab Col-0 \tab rosette \tab 4-week-old \tab wounding \cr
    GSM347112 \tab Arabidopsis rosette_Light Wounded_Rep1 \tab GSE13803 \tab abiotic \tab Col-0 \tab rosette \tab 4-week-old \tab wounding \cr
    GSM347114 \tab Arabidopsis rosette_Light Wounded_Rep2 \tab GSE13803 \tab abiotic \tab Col-0 \tab rosette \tab 4-week-old \tab wounding \cr
    GSM4125057 \tab wt Cd100 rep 1 \tab GSE138943 \tab abiotic \tab Col-0 \tab seedling \tab 2-week-old \tab cadmium \cr
    GSM4125058 \tab wt Cd100 rep 2 \tab GSE138943 \tab abiotic \tab Col-0 \tab seedling \tab 2-week-old \tab cadmium \cr
    GSM4125059 \tab wt Cd100 rep 3 \tab GSE138943 \tab abiotic \tab Col-0 \tab seedling \tab 2-week-old \tab cadmium \cr
    GSM392180 \tab WT-16ºC-RepA \tab GSE15689 \tab abiotic \tab Col-0 \tab seedling \tab 10-day-old \tab cold \cr
    GSM392263 \tab WT-16ºC-RepB \tab GSE15689 \tab abiotic \tab Col-0 \tab seedling \tab 10-day-old \tab cold \cr
    GSM392273 \tab WT-16ºC-RepC \tab GSE15689 \tab abiotic \tab Col-0 \tab seedling \tab 10-day-old \tab cold \cr
    GSM5584760 \tab col submerged for 48h, biological rep1 \tab GSE184340 \tab abiotic \tab Col-0 \tab leaf \tab 4-week-old \tab submergence (hypoxia) \cr
    GSM5584761 \tab col submerged for 48h, biological rep2 \tab GSE184340 \tab abiotic \tab Col-0 \tab leaf \tab 4-week-old \tab submergence (hypoxia) \cr
    GSM5584762 \tab col submerged for 48h, biological rep3 \tab GSE184340 \tab abiotic \tab Col-0 \tab leaf \tab 4-week-old \tab submergence (hypoxia) \cr
    GSM539316 \tab WT-sib, treatment N-L+, 1 biological rep \tab GSE21601 \tab abiotic \tab Col-0 \tab seedling \tab 2-week-old \tab nitrogen deficiency \cr
    GSM539317 \tab WT-sib, treatment N-L+, 2 biological rep \tab GSE21601 \tab abiotic \tab Col-0 \tab seedling \tab 2-week-old \tab nitrogen deficiency \cr
    GSM539318 \tab WT-sib, treatment N-L+, 3 biological rep \tab GSE21601 \tab abiotic \tab Col-0 \tab seedling \tab 2-week-old \tab nitrogen deficiency \cr
    GSM566862 \tab Split KNO3 roots - 2hours - repeat 2 \tab GSE22966 \tab abiotic \tab Col-0 \tab root \tab 2-week-old \tab root splitting (contrast in nitrogen supply) \cr
    GSM566863 \tab Split KNO3 roots - 2hours - repeat 3 \tab GSE22966 \tab abiotic \tab Col-0 \tab root \tab 2-week-old \tab root splitting (contrast in nitrogen supply) \cr
    GSM566861 \tab Split KNO3 roots - 2hours - repeat 1 \tab GSE22966 \tab abiotic \tab Col-0 \tab root \tab 2-week-old \tab root splitting (contrast in nitrogen supply) \cr
    GSM643107 \tab leaf 30\%Inhibition rep1 [Arabidopsis] \tab GSE26197 \tab abiotic \tab Col-0 \tab first fully expanded leaf \tab  \tab photosynthesis at 30\% \cr
    GSM643108 \tab leaf 30\%Inhibition rep2 [Arabidopsis] \tab GSE26197 \tab abiotic \tab Col-0 \tab first fully expanded leaf \tab  \tab photosynthesis at 30\% \cr
    GSM643110 \tab leaf 30\%Inhibition rep4 [Arabidopsis] \tab GSE26197 \tab abiotic \tab Col-0 \tab first fully expanded leaf \tab  \tab photosynthesis at 30\% \cr
    GSM901075 \tab osmotic_10day_leaf_rep1 \tab GSE36789 \tab abiotic \tab Col-0 \tab leaf \tab 10-day-old \tab mannitol \cr
    GSM901076 \tab osmotic_10day_leaf_rep2 \tab GSE36789 \tab abiotic \tab Col-0 \tab leaf \tab 10-day-old \tab mannitol \cr
    GSM901077 \tab osmotic_10day_leaf_rep3 \tab GSE36789 \tab abiotic \tab Col-0 \tab leaf \tab 10-day-old \tab mannitol \cr
    GSM2896650 \tab Columbia, biological rep1 \tab GSE108376 \tab control \tab Col-0 \tab seedling \tab 12-day-old \tab Murashige and Skoog medium; 16h photoperiod; 75 microE light intensity \cr
    GSM2896651 \tab Columbia, biological rep2 \tab GSE108376 \tab control \tab Col-0 \tab seedling \tab 12-day-old \tab Murashige and Skoog medium; 16h photoperiod; 75 microE light intensity \cr
    GSM2896652 \tab Columbia, biological rep3 \tab GSE108376 \tab control \tab Col-0 \tab seedling \tab 12-day-old \tab Murashige and Skoog medium; 16h photoperiod; 75 microE light intensity \cr
    GSM3639524 \tab seedings under mock treatment, biological rep1 \tab GSE127830 \tab control \tab Col-0 \tab seedings growing on filter paper \tab 10-day-old \tab filter paper medium; 12h photoperiod; 120 microE light intensity; 60 \% relative humidity; 25 °C; mock treatment (DMSO) \cr
    GSM3639525 \tab seedings under mock treatment, biological rep2 \tab GSE127830 \tab control \tab Col-0 \tab seedings growing on filter paper \tab 10-day-old \tab filter paper medium; 12h photoperiod; 120 microE light intensity; 60 \% relative humidity; 25 °C; mock treatment (DMSO) \cr
    GSM3639526 \tab seedings under mock treatment, biological rep3 \tab GSE127830 \tab control \tab Col-0 \tab seedings growing on filter paper \tab 10-day-old \tab filter paper medium; 12h photoperiod; 120 microE light intensity; 60 \% relative humidity; 25 °C; mock treatment (DMSO) \cr
    GSM538647 \tab Seedlings untreated, biological rep1 \tab GSE21556 \tab control \tab Col-0 \tab seedling \tab 8-day-old \tab 1:3 vermiculite:soil; 16h photoperiod; 80-100 microE light intensity; 50-60 \% relative humidity; 22 °C \cr
    GSM538648 \tab Seedlings untreated, biological rep2 \tab GSE21556 \tab control \tab Col-0 \tab seedling \tab 8-day-old \tab 1:3 vermiculite:soil; 16h photoperiod; 80-100 microE light intensity; 50-60 \% relative humidity; 22 °C \cr
    GSM538649 \tab Seedlings untreated, biological rep3 \tab GSE21556 \tab control \tab Col-0 \tab seedling \tab 8-day-old \tab 1:3 vermiculite:soil; 16h photoperiod; 80-100 microE light intensity; 50-60 \% relative humidity; 22 °C \cr
    GSM542673 \tab control_4h-1 \tab GSE21786 \tab control \tab Col-0 \tab seedling \tab 5-day-old \tab Murashige and Skoog agar plate with sucrose; 22 °C \cr
    GSM542674 \tab control_4h-2 \tab GSE21786 \tab control \tab Col-0 \tab seedling \tab 5-day-old \tab Murashige and Skoog agar plate with sucrose; 22 °C \cr
    GSM542675 \tab control_4h-3 \tab GSE21786 \tab control \tab Col-0 \tab seedling \tab 5-day-old \tab Murashige and Skoog agar plate with sucrose; 22 °C \cr
    GSM1028173 \tab Col-0 biological rep1 \tab GSE41958 \tab control \tab Col-0 \tab seedling \tab 15-day-old \tab Murashige and Skoog medium with sucrose; 16h photoperiod; 5 hours zeitgeber time \cr
    GSM1028174 \tab Col-0 biological rep2 \tab GSE41958 \tab control \tab Col-0 \tab seedling \tab 15-day-old \tab Murashige and Skoog medium with sucrose; 16h photoperiod; 5 hours zeitgeber time \cr
    GSM1028175 \tab Col-0 biological rep3 \tab GSE41958 \tab control \tab Col-0 \tab seedling \tab 15-day-old \tab Murashige and Skoog medium with sucrose; 16h photoperiod; 5 hours zeitgeber time \cr
    GSM250982 \tab Col-0 without treatment, biological rep1 \tab GSE9957 \tab control \tab Col-0 \tab shoot \tab 3-week-old \tab 14h photoperiod; 22 °C \cr
    GSM250983 \tab Col-0 without treatment, biological rep2 \tab GSE9957 \tab control \tab Col-0 \tab shoot \tab 3-week-old \tab 14h photoperiod; 22 °C \cr
    GSM250984 \tab Col-0 without treatment, biological rep3 \tab GSE9957 \tab control \tab Col-0 \tab shoot \tab 3-week-old \tab 14h photoperiod; 22 °C \cr
    GSM1202382 \tab Col-0, WT, at HL 0h, biological rep1 \tab GSE49596 \tab control \tab Col-0 \tab shoot \tab 12-day-old \tab half Murashige and Skoog medium with sucrose; 24h photoperiod; 10 microE light intensity \cr
    GSM1202383 \tab Col-0, WT, at HL 0h, biological rep2 \tab GSE49596 \tab control \tab Col-0 \tab shoot \tab 12-day-old \tab half Murashige and Skoog medium with sucrose; 24h photoperiod; 10 microE light intensity \cr
    GSM1202384 \tab Col-0, WT, at HL 0h, biological rep3 \tab GSE49596 \tab control \tab Col-0 \tab shoot \tab 12-day-old \tab half Murashige and Skoog medium with sucrose; 24h photoperiod; 10 microE light intensity \cr
    GSM1257966 \tab WT complete, biological rep1 \tab GSE52046 \tab control \tab Col-0 \tab seedling \tab 1-week-old \tab complete medium \cr
    GSM1257967 \tab WT complete, biological rep2 \tab GSE52046 \tab control \tab Col-0 \tab seedling \tab 1-week-old \tab complete medium \cr
    GSM1257968 \tab WT complete, biological rep3 \tab GSE52046 \tab control \tab Col-0 \tab seedling \tab 1-week-old \tab complete medium \cr
    GSM128757 \tab Mittler_2-1_wildtype_Rep1_ATH1 \tab GSE5530 \tab control \tab Col-0 \tab seedling \tab 5-day old \tab Murashige and Skoog liquid medium \cr
    GSM128758 \tab Mittler_2-2_wildtype_Rep2_ATH1 \tab GSE5530 \tab control \tab Col-0 \tab seedling \tab 5-day old \tab Murashige and Skoog liquid medium \cr
    GSM128759 \tab Mittler_2-3_wildtype_Rep3_ATH1 \tab GSE5530 \tab control \tab Col-0 \tab seedling \tab 5-day old \tab Murashige and Skoog liquid medium \cr
    GSM1875295 \tab WT_rep1 \tab GSE72954 \tab control \tab Col-0 \tab seedling \tab 9-day-old \tab germination medium with glucose; 16h photoperiod;  23 °C day/ 21 °C night \cr
    GSM1875296 \tab WT_rep2 \tab GSE72954 \tab control \tab Col-0 \tab seedling \tab 9 day old \tab germination medium with glucose; 16h photoperiod;  23 °C day/ 21 °C night \cr
    GSM1875297 \tab WT_rep3 \tab GSE72954 \tab control \tab Col-0 \tab seedling \tab 9 day old \tab germination medium with glucose; 16h photoperiod;  23 °C day/ 21 °C night \cr
    GSM9595 \tab Col_control_I_1 \tab GSE629 \tab control \tab Col-0 \tab seedling \tab 1-week-old \tab Murashige and Skoog liquid medium with sucrose; 16h photoperiod; 24 °C \cr
    GSM9596 \tab Col_control_I_2 \tab GSE629 \tab control \tab Col-0 \tab seedling \tab 1-week-old \tab Murashige and Skoog liquid medium with sucrose; 16h photoperiod; 24 °C \cr
    GSM9597 \tab Col_control_I_3 \tab GSE629 \tab control \tab Col-0 \tab seedling \tab 1-week-old \tab Murashige and Skoog liquid medium with sucrose; 16h photoperiod; 24 °C \cr
    GSM2144831 \tab Col-0, biological rep1 \tab GSE81218 \tab control \tab Col-0 \tab seedling \tab 10-day-old \tab Murashige and Skoog agar plate; 16h photoperiod; 120 microE light intensity; 21 °C; mock treatment (water) \cr
    GSM2144832 \tab Col-0, biological rep2 \tab GSE81218 \tab control \tab Col-0 \tab seedling \tab 10-day-old \tab Murashige and Skoog agar plate; 16h photoperiod; 120 microE light intensity; 21 °C; mock treatment (water) \cr
    GSM2144833 \tab Col-0, biological rep3 \tab GSE81218 \tab control \tab Col-0 \tab seedling \tab 10-day-old \tab Murashige and Skoog agar plate; 16h photoperiod; 120 microE light intensity; 21 °C; mock treatment (water) \cr
    GSM591834 \tab WT mock treated_biological rep 1 \tab GSE24052 \tab control \tab Col-0 -800GSTF8::LUC \tab seedling \tab 4-day-old \tab Murashige and Skoog agar plate; 16h photoperiod; 21 °C \cr
    GSM591835 \tab WT mock treated_biological rep 2 \tab GSE24052 \tab control \tab Col-0 -800GSTF8::LUC \tab seedling \tab 4-day-old \tab Murashige and Skoog agar plate; 16h photoperiod; 21 °C \cr
    GSM591836 \tab WT mock treated_biological rep 3 \tab GSE24052 \tab control \tab Col-0 -800GSTF8::LUC \tab seedling \tab 4-day-old \tab Murashige and Skoog agar plate; 16h photoperiod; 21 °C \cr
    GSM617578 \tab Col, rep1 \tab GSE25134 \tab control \tab Col-0 \tab seedling \tab 4.5-day-old \tab half Murashige and Skoog medium with sucrose; 0h photoperiod \cr
    GSM617579 \tab Col, rep2 \tab GSE25134 \tab control \tab Col-0 \tab seedling \tab 4.5-day-old \tab half Murashige and Skoog medium with sucrose; 0h photoperiod \cr
    GSM617580 \tab Col, rep3 \tab GSE25134 \tab control \tab Col-0 \tab seedling \tab 4.5-day-old \tab half Murashige and Skoog medium with sucrose; 0h photoperiod \cr
    GSM656490 \tab Col-0, 0 hrs, bio rep 2 \tab GSE26679 \tab control \tab Col-0 \tab shoot \tab 4-week-old \tab Metromix soil; 9h photoperiod; 23 °C \cr
    GSM656491 \tab Col-0, 0 hrs, bio rep 3 \tab GSE26679 \tab control \tab Col-0 \tab shoot \tab 4-week-old \tab Metromix soil; 9h photoperiod; 23 °C \cr
    GSM656492 \tab Col-0, 0 hrs, bio rep 4 \tab GSE26679 \tab control \tab Col-0 \tab shoot \tab 4-week-old \tab Metromix soil; 9h photoperiod; 23 °C \cr
    GSM738872 \tab Col-0 WT seedling at 5 dpg, biological rep2 \tab GSE29814 \tab control \tab Col-0 \tab seedling \tab 5-day-old \tab Murashige and Skoog agar plate; 18h photoperiod; 21 °C \cr
    GSM738873 \tab Col-0 WT seedling at 5 dpg, biological rep3 \tab GSE29814 \tab control \tab Col-0 \tab seedling \tab 5-day-old \tab Murashige and Skoog agar plate; 18h photoperiod; 21 °C \cr
    GSM738874 \tab Col-0 WT seedling at 5 dpg, biological rep4 \tab GSE29814 \tab control \tab Col-0 \tab seedling \tab 5-day-old \tab Murashige and Skoog agar plate; 18h photoperiod; 21 °C \cr
    GSM761617 \tab Columbia_0min_biorep1 \tab GSE30702 \tab control \tab Col-0 \tab seedling \tab 1-week-old \tab Murashige and Skoog medium with sucrose; 22 °C \cr
    GSM761620 \tab Columbia_0min_biorep4 \tab GSE30702 \tab control \tab Col-0 \tab seedling \tab 1-week-old \tab Murashige and Skoog medium with sucrose; 22 °C \cr
    GSM761621 \tab Columbia_0min_biorep5 \tab GSE30702 \tab control \tab Col-0 \tab seedling \tab 1-week-old \tab Murashige and Skoog medium with sucrose; 22 °C \cr
    GSM856033 \tab seedling at 14 day_mock treatment_biological replica 1 \tab GSE34837 \tab control \tab Col-0 \tab seedling \tab 2-week-old \tab Murashige and Skoog agar plate; 16h photoperiod; 22 °C; mock treatment (ethanol) \cr
    GSM856034 \tab seedling at 14 day_mock treatment_biological replica 2 \tab GSE34837 \tab control \tab Col-0 \tab seedling \tab 2-week-old \tab Murashige and Skoog agar plate; 16h photoperiod; 22 °C; mock treatment (ethanol) \cr
    GSM856035 \tab seedling at 14 day_mock treatment_biological replica 3 \tab GSE34837 \tab control \tab Col-0 \tab seedling \tab 2-week-old \tab Murashige and Skoog agar plate; 16h photoperiod; 22 °C; mock treatment (ethanol) \cr
    GSM686101 \tab wild type Col-0, rep1 \tab GSE27704 \tab control \tab Col-0 \tab seedling \tab 3.25-day-old \tab Murashige and Skoog medium; 72 hours darkness + 6 hours 100 microE light intensity; 23°C \cr
    GSM686102 \tab wild type Col-0, rep2 \tab GSE27704 \tab control \tab Col-0 \tab seedling \tab 3.25-day-old \tab Murashige and Skoog medium; 72 hours darkness + 6 hours 100 microE light intensity; 23°C \cr
    GSM686103 \tab wild type Col-0, rep3 \tab GSE27704 \tab control \tab Col-0 \tab seedling \tab 3.25-day-old \tab Murashige and Skoog medium; 72 hours darkness + 6 hours 100 microE light intensity; 23°C \cr
    GSM1304062 \tab Leaves mock treated 3 h, biological rep1 \tab GSE53957 \tab control \tab Col-0 \tab rosette leaf \tab 3-week-old \tab soil; 11h photoperiod; 100 microE light intensity; 70\\\% humidity; 21 °C; mock treatment (methanol) \cr
    GSM1304063 \tab Leaves mock treated 3 h, biological rep2 \tab GSE53957 \tab control \tab Col-0 \tab rosette leaf \tab 3-week-old \tab soil; 11h photoperiod; 100 microE light intensity; 70\% humidity; 21 °C; mock treatment (methanol) \cr
    GSM1304064 \tab Leaves mock treated 3 h, biological rep3 \tab GSE53957 \tab control \tab Col-0 \tab rosette leaf \tab 3-week-old \tab soil; 11h photoperiod; 100 microE light intensity; 70\% humidity; 21 °C; mock treatment (methanol) \cr
    GSM2482338 \tab Col-0 developing seeds, 10-11 days after flower opening, biological replicate 1 [re-analysis] \tab GSE94763 \tab developmental \tab Col-0 \tab silique \tab 10-11 days after flower opening \tab control \cr
    GSM2482339 \tab Col-0 developing seeds, 10-11 days after flower opening, biological replicate 2 [re-analysis] \tab GSE94763 \tab developmental \tab Col-0 \tab silique \tab 10-11 days after flower opening \tab control \cr
    GSM2482340 \tab Col-0 developing seeds, 10-11 days after flower opening, biological replicate 3 [re-analysis] \tab GSE94763 \tab developmental \tab Col-0 \tab silique \tab 10-11 days after flower opening \tab control \cr
    GSM142750 \tab MJ001_ATH1_A1-jones-WT1 \tab GSE6165 \tab developmental \tab Col-0 \tab small excised segments of primary root tissue corresponding to root-hair elongation zone \tab 12-day-old \tab control \cr
    GSM142751 \tab MJ001_ATH1_A2-jones-WT2 \tab GSE6165 \tab developmental \tab Col-0 \tab small excised segments of primary root tissue corresponding to root-hair elongation zone \tab 12-day-old \tab control \cr
    GSM142754 \tab MJ001_ATH1_A5-jones-WT-Rep3 \tab GSE6165 \tab developmental \tab Col-0 \tab small excised segments of primary root tissue corresponding to root-hair elongation zone \tab 12-day-old \tab control \cr
    GSM184901 \tab Arabidopsis, root cells, stele, standard conditions, replicate 1 \tab GSE7641 \tab developmental \tab Col-0 WOL::GFP \tab root stele \tab 6-day-old \tab control \cr
    GSM184902 \tab Arabidopsis, root cells, stele, standard conditions, replicate 2 \tab GSE7641 \tab developmental \tab Col-0 WOL::GFP \tab root stele \tab 6-day-old \tab control \cr
    GSM184903 \tab Arabidopsis, root cells, stele, standard conditions, replicate 3 \tab GSE7641 \tab developmental \tab Col-0 WOL::GFP \tab root stele \tab 6-day-old \tab control \cr
    GSM131831 \tab Quick_A27_0-0hr_Rep2_ATH1 \tab GSE5639 \tab developmental \tab Col-0 \tab developing leaf insertions 19-21 \tab 4-week old \tab control \cr
    GSM131833 \tab Quick_A53_0-0hr_Rep3_ATH1 \tab GSE5639 \tab developmental \tab Col-0 \tab developing leaf insertions 19-21 \tab 4-week old \tab control \cr
    GSM131835 \tab Quick_A79_0-0hr_Rep4_ATH1 \tab GSE5639 \tab developmental \tab Col-0 \tab developing leaf insertions 19-21 \tab 4-week old \tab control \cr
    GSM1511283 \tab wild-type (Col-0) 8 DAF Seed 1 \tab GSE61684 \tab developmental \tab Col-0 \tab developing seeds dissected from sliques \tab 8 days after floweing \tab control \cr
    GSM1511284 \tab wild-type (Col-0) 8 DAF Seed 2 \tab GSE61684 \tab developmental \tab Col-0 \tab developing seeds dissected from sliques \tab 8 days after floweing \tab control \cr
    GSM1511285 \tab wild-type (Col-0) 8 DAF Seed 3 \tab GSE61684 \tab developmental \tab Col-0 \tab developing seeds dissected from sliques \tab 8 days after floweing \tab control \cr
    GSM133978 \tab Birnbaum_1-8_StageII-1_Rep1_ATH1 \tab GSE5749 \tab developmental \tab Col-0 \tab lateral root stage II \tab  \tab control \cr
    GSM133979 \tab Birnbaum_1-9_StageII-2_Rep2_ATH1 \tab GSE5749 \tab developmental \tab Col-0 \tab lateral root stage II \tab  \tab control \cr
    GSM133980 \tab Birnbaum_1-10_StageII-3_Rep3_ATH1 \tab GSE5749 \tab developmental \tab Col-0 \tab lateral root stage II \tab  \tab control \cr
    GSM133767 \tab Lindsey_1-19_torpedo-basal_Rep4_ATH1 \tab GSE5730 \tab developmental \tab Col-0 \tab torpedo embryo, basal tissue \tab  \tab control \cr
    GSM133768 \tab Lindsey_1-20_torpedo-basal_Rep5_ATH1 \tab GSE5730 \tab developmental \tab Col-0 \tab torpedo embryo, basal tissue \tab  \tab control \cr
    GSM133769 \tab Lindsey_1-21_torpedo-basal_Rep6_ATH1 \tab GSE5730 \tab developmental \tab Col-0 \tab torpedo embryo, basal tissue \tab  \tab control \cr
    GSM128782 \tab Somerville_1-5_flower-GC6_Rep2_ATH1 \tab GSE5533 \tab developmental \tab Col-0 \tab Developed flowers and unopened buds \tab 29-day-old \tab control \cr
    GSM128783 \tab Somerville_1-6_flower-GH5_Rep1_ATH1 \tab GSE5533 \tab developmental \tab Col-0 \tab Developed flowers and unopened buds \tab 29-day-old \tab control \cr
    GSM128784 \tab Somerville_1-7_flower-GH6_Rep2_ATH1 \tab GSE5533 \tab developmental \tab Col-0 \tab Developed flowers and unopened buds \tab 29-day-old \tab control \cr
    GSM1296372 \tab WT(Col-0), biological rep1 \tab GSE53580 \tab developmental \tab Col-0 \tab basal stem \tab 12-17 cm height plants \tab control \cr
    GSM1296373 \tab WT(Col-0), biological rep2 \tab GSE53580 \tab developmental \tab Col-0 \tab basal stem \tab 12-17 cm height plants \tab control \cr
    GSM1296374 \tab WT(Col-0), biological rep3 \tab GSE53580 \tab developmental \tab Col-0 \tab basal stem \tab 12-17 cm height plants \tab control \cr
    GSM1133319 \tab Day Time Control 1 \tab GSE46621 \tab developmental \tab Col-0 \tab shoot \tab 15-day-old \tab control \cr
    GSM1133320 \tab Day Time Control 2 \tab GSE46621 \tab developmental \tab Col-0 \tab shoot \tab 15-day-old \tab control \cr
    GSM1133321 \tab Day Time Control 3 \tab GSE46621 \tab developmental \tab Col-0 \tab shoot \tab 15-day-old \tab control \cr
    GSM1091681 \tab Columbia Clipped_Biol_Rep 1 \tab GSE44781 \tab developmental \tab Col-0 \tab secondary meristem \tab at bolting \tab clipped \cr
    GSM1091682 \tab Columbia Clipped_Biol_Rep 2 \tab GSE44781 \tab developmental \tab Col-0 \tab secondary meristem \tab at bolting \tab clipped \cr
    GSM1091683 \tab Columbia Clipped_Biol_Rep 3 \tab GSE44781 \tab developmental \tab Col-0 \tab secondary meristem \tab at bolting \tab clipped \cr
    GSM871248 \tab pericycle control 1 \tab GSE35580 \tab developmental \tab Col-0 \tab pericycle \tab 1-week-old \tab control \cr
    GSM871249 \tab pericycle control 2 \tab GSE35580 \tab developmental \tab Col-0 \tab pericycle \tab 1-week-old \tab control \cr
    GSM871250 \tab pericycle control 3 \tab GSE35580 \tab developmental \tab Col-0 \tab pericycle \tab 1-week-old \tab control \cr
    GSM131540 \tab ATGE_26_A \tab GSE5630 \tab developmental \tab Col-0 \tab cauline leaf \tab  \tab control \cr
    GSM131541 \tab ATGE_26_B \tab GSE5630 \tab developmental \tab Col-0 \tab cauline leaf \tab  \tab control \cr
    GSM131542 \tab ATGE_26_C \tab GSE5630 \tab developmental \tab Col-0 \tab cauline leaf \tab  \tab control \cr
    GSM390161 \tab lateral nectary, stage 14-15 rep1 \tab GSE15601 \tab developmental \tab Col-0 \tab lateral nectaries (from stage 14-15 flowers) \tab 30-35-day-old \tab control \cr
    GSM390162 \tab lateral nectary, stage 14-15 rep2 \tab GSE15601 \tab developmental \tab Col-0 \tab lateral nectaries (from stage 14-15 flowers) \tab 30-35-day-old \tab control \cr
    GSM390163 \tab lateral nectary, stage 14-15 rep3 \tab GSE15601 \tab developmental \tab Col-0 \tab lateral nectaries (from stage 14-15 flowers) \tab 30-35-day-old \tab control \cr
    GSM2044875 \tab Empty vector plants without DEX treatment at T9, biological rep1 \tab GSE77153 \tab developmental \tab Col-0 VP16-GR \tab hypocotyl \tab 6-day-old \tab darkness; mock treatment (DMSO) \cr
    GSM2044876 \tab Empty vector plants without DEX treatment at T9, biological rep2 \tab GSE77153 \tab developmental \tab Col-0 VP16-GR \tab hypocotyl \tab 6-day-old \tab darkness; mock treatment (DMSO) \cr
    GSM2044877 \tab Empty vector plants without DEX treatment at T9, biological rep3 \tab GSE77153 \tab developmental \tab Col-0 VP16-GR \tab hypocotyl \tab 6-day-old \tab darkness; mock treatment (DMSO) \cr
    GSM184503 \tab Pericycle root cells 2hr KCl control treated, biological rep1 \tab GSE7631 \tab developmental \tab Col-0 \tab pericycle root cells \tab 12-day-old \tab control \cr
    GSM184504 \tab Pericycle root cells 2hr KCl control treated, biological rep2 \tab GSE7631 \tab developmental \tab Col-0 \tab pericycle root cells \tab 12-day-old \tab control \cr
    GSM184505 \tab Pericycle root cells 2hr KCl control treated, biological rep3 \tab GSE7631 \tab developmental \tab Col-0 \tab pericycle root cells \tab 12-day-old \tab control \cr
    GSM1509558 \tab ScionFlowerBud1 \tab GSE61631 \tab developmental \tab Col-0 \tab scion (flower bud) \tab  \tab control \cr
    GSM1509559 \tab ScionFlowerBud2 \tab GSE61631 \tab developmental \tab Col-0 \tab scion (flower bud) \tab  \tab control \cr
    GSM1509560 \tab ScionFlowerBud3 \tab GSE61631 \tab developmental \tab Col-0 \tab scion (flower bud) \tab  \tab control \cr
    GSM1535437 \tab Cells whole embryo, biological rep 1 \tab GSE60242 \tab developmental \tab Col-0 \tab embryo \tab 8-32 cell stage \tab control \cr
    GSM1535438 \tab Cells whole embryo, biological rep 2 \tab GSE60242 \tab developmental \tab Col-0 \tab embryo \tab 8-32 cell stage \tab control \cr
    GSM1535439 \tab Cells whole embryo, biological rep 3 \tab GSE60242 \tab developmental \tab Col-0 \tab embryo \tab 8-32 cell stage \tab control \cr
    GSM1289220 \tab Col-0 plants under control conditions, rep1 \tab GSE53308 \tab developmental \tab Col-0 \tab rosette \tab 28-day-old \tab control \cr
    GSM1289221 \tab Col-0 plants under control conditions, rep2 \tab GSE53308 \tab developmental \tab Col-0 \tab rosette \tab 28-day-old \tab control \cr
    GSM1289222 \tab Col-0 plants under control conditions, rep3 \tab GSE53308 \tab developmental \tab Col-0 \tab rosette \tab 28-day-old \tab control \cr
    GSM433643 \tab 4hPT, biological rep2 \tab GSE17343 \tab developmental \tab Col-0 \tab pollen tube \tab grown for 4 hours \tab liquid pollen tube growth medium \cr
    GSM433644 \tab 4hPT, biological rep3 \tab GSE17343 \tab developmental \tab Col-0 \tab pollen tube \tab grown for 4 hours \tab liquid pollen tube growth medium \cr
    GSM433645 \tab 4hPT, biological rep4 \tab GSE17343 \tab developmental \tab Col-0 \tab pollen tube \tab grown for 4 hours \tab liquid pollen tube growth medium \cr
    GSM325126 \tab auxin-treated seedlings, biological rep1 \tab GSE12964 \tab chemical \tab Col-0 \tab seedling \tab 10-day-old \tab 2,4-dichlorophenoxyacetic acid \cr
    GSM325127 \tab auxin-treated seedlings, biological rep2 \tab GSE12964 \tab chemical \tab Col-0 \tab seedling \tab 10-day-old \tab 2,4-dichlorophenoxyacetic acid \cr
    GSM325128 \tab auxin-treated seedlings, biological rep3 \tab GSE12964 \tab chemical \tab Col-0 \tab seedling \tab 10-day-old \tab 2,4-dichlorophenoxyacetic acid \cr
    GSM373534 \tab SA treated, biological rep1 \tab GSE14961 \tab chemical \tab Col-0 \tab seedling \tab 10-day-old \tab salicylic acid \cr
    GSM373535 \tab SA treated, biological rep2 \tab GSE14961 \tab chemical \tab Col-0 \tab seedling \tab 10-day-old \tab salicylic acid \cr
    GSM373536 \tab SA treated, biological rep3 \tab GSE14961 \tab chemical \tab Col-0 \tab seedling \tab 10-day-old \tab salicylic acid \cr
    GSM469825 \tab Gronlund_1-21_BR3+GA-180mins_Rep3_ATH1 \tab GSE18985 \tab chemical \tab Col-0 \tab hypocotyl \tab 9-day-old \tab gibberellin 4 \cr
    GSM469832 \tab Gronlund_1-14_BR2+GA-180mins_Rep2_ATH1 \tab GSE18985 \tab chemical \tab Col-0 \tab hypocotyl \tab 9-day-old \tab gibberellin 4 \cr
    GSM469839 \tab Gronlund_1-7_BR1+GA-180mins_Rep1_ATH1 \tab GSE18985 \tab chemical \tab Col-0 \tab hypocotyl \tab 9-day-old \tab gibberellin 4 \cr
    GSM679546 \tab Mature 1-1 \tab GSE27508 \tab chemical \tab Col-0 \tab hypocotyl \tab 35-day-old \tab 2,4-dichlorophenoxyacetic acid \cr
    GSM679547 \tab Mature 1-2 \tab GSE27508 \tab chemical \tab Col-0 \tab hypocotyl \tab 35-day-old \tab 2,4-dichlorophenoxyacetic acid \cr
    GSM679548 \tab Mature 1-3 \tab GSE27508 \tab chemical \tab Col-0 \tab hypocotyl \tab 35-day-old \tab 2,4-dichlorophenoxyacetic acid \cr
    GSM713255 \tab ABA_3 \tab GSE28800 \tab chemical \tab Col-0 \tab seedling \tab 12-day-old \tab abscisic acid \cr
    GSM713256 \tab ABA_2 \tab GSE28800 \tab chemical \tab Col-0 \tab seedling \tab 12-day-old \tab abscisic acid \cr
    GSM713257 \tab ABA_1 \tab GSE28800 \tab chemical \tab Col-0 \tab seedling \tab 12-day-old \tab abscisic acid \cr
    GSM744705 \tab SAM_L-AOPP_Rep1 \tab GSE30092 \tab chemical \tab Col-0 \tab shoot apical meristem \tab 4-day-old \tab 2-aminoxy-3-phenylpropionic acid \cr
    GSM744706 \tab SAM_L-AOPP_Rep2 \tab GSE30092 \tab chemical \tab Col-0 \tab shoot apical meristem \tab 4-day-old \tab 2-aminoxy-3-phenylpropionic acid \cr
    GSM744707 \tab SAM_L-AOPP_Rep3 \tab GSE30092 \tab chemical \tab Col-0 \tab shoot apical meristem \tab 4-day-old \tab 2-aminoxy-3-phenylpropionic acid \cr
    GSM1053027 \tab Roots at T2 naxillin, rep1 \tab GSE42896 \tab chemical \tab Col-0 \tab root \tab 3-day-old \tab naxillin \cr
    GSM1053028 \tab Roots at T2 naxillin, rep2 \tab GSE42896 \tab chemical \tab Col-0 \tab root \tab 3-day-old \tab naxillin \cr
    GSM1053029 \tab Roots at T2 naxillin, rep3 \tab GSE42896 \tab chemical \tab Col-0 \tab root \tab 3-day-old \tab naxillin \cr
    GSM1399594 \tab PEO-IAA replicate1 \tab GSE58028 \tab chemical \tab Col-0 \tab seedling \tab 1-week-old \tab 2-(1H-Indol-3-yl)-4-oxo-4-phenyl-butyric acid \cr
    GSM1399595 \tab PEO-IAA replicate2 \tab GSE58028 \tab chemical \tab Col-0 \tab seedling \tab 1-week-old \tab 2-(1H-Indol-3-yl)-4-oxo-4-phenyl-butyric acid \cr
    GSM1399596 \tab PEO-IAA replicate3 \tab GSE58028 \tab chemical \tab Col-0 \tab seedling \tab 1-week-old \tab 2-(1H-Indol-3-yl)-4-oxo-4-phenyl-butyric acid \cr
    GSM1516365 \tab Col-0 leaves, MeJA treatment, replicate 2 \tab GSE61884 \tab chemical \tab Col-0 \tab leaf \tab 4-week-old \tab methyl-jasmonates \cr
    GSM1516366 \tab Col-0 leaves, MeJA treatment, replicate 3 \tab GSE61884 \tab chemical \tab Col-0 \tab leaf \tab 4-week-old \tab methyl-jasmonates \cr
    GSM1516367 \tab Col-0 leaves, MeJA treatment, replicate 4 \tab GSE61884 \tab chemical \tab Col-0 \tab leaf \tab 4-week-old \tab methyl-jasmonates \cr
    GSM269475 \tab Col-0, 24h after BTH, rep-A \tab GSE10646 \tab chemical \tab Col-0 \tab plant \tab 5-week-old \tab benzothiadiazole \cr
    GSM269477 \tab Col-0, 24h after BTH, rep-B \tab GSE10646 \tab chemical \tab Col-0 \tab plant \tab 5-week-old \tab benzothiadiazole \cr
    GSM269478 \tab Col-0, 24h after BTH, rep-C \tab GSE10646 \tab chemical \tab Col-0 \tab plant \tab 5-week-old \tab benzothiadiazole \cr
    GSM297833 \tab tunicamycin_rep2 \tab GSE11758 \tab chemical \tab Col-0 \tab mature leaf \tab 5-week-old \tab tunicamycin \cr
    GSM297839 \tab tunicamycin_rep3 \tab GSE11758 \tab chemical \tab Col-0 \tab mature leaf \tab 5-week-old \tab tunicamycin \cr
    GSM297843 \tab tunicamycin_rep5 \tab GSE11758 \tab chemical \tab Col-0 \tab mature leaf \tab 5-week-old \tab tunicamycin \cr
    GSM323078 \tab Col-0_norflurazon_1 \tab GSE12887 \tab chemical \tab Col-0 \tab seedling \tab 5-day-old \tab norflurazon \cr
    GSM323079 \tab Col-0_norflurazon_2 \tab GSE12887 \tab chemical \tab Col-0 \tab seedling \tab 5-day-old \tab norflurazon \cr
    GSM323080 \tab Col-0_norflurazon_3 \tab GSE12887 \tab chemical \tab Col-0 \tab seedling \tab 5-day-old \tab norflurazon \cr
    GSM347944 \tab Col-0_DCA_6d_rep1 \tab GSE13833 \tab chemical \tab Col-0 \tab shoot \tab 2-week-old \tab 3,5-dichloroanthranilic acid \cr
    GSM347945 \tab Col-0_DCA_6d_rep2 \tab GSE13833 \tab chemical \tab Col-0 \tab shoot \tab 2-week-old \tab 3,5-dichloroanthranilic acid \cr
    GSM347946 \tab Col-0_DCA_6d_rep3 \tab GSE13833 \tab chemical \tab Col-0 \tab shoot \tab 2-week-old \tab 3,5-dichloroanthranilic acid \cr
    GSM469812 \tab Rylott_1-3_TNT-treated-seedlings_Rep3_ATH1 \tab GSE18983 \tab chemical \tab Col-0 \tab seedling \tab 1.05-week-old \tab trinitrotoluen \cr
    GSM469813 \tab Rylott_1-2_TNT-treated-seedlings_Rep2_ATH1 \tab GSE18983 \tab chemical \tab Col-0 \tab seedling \tab 1.05-week-old \tab trinitrotoluen \cr
    GSM469814 \tab Rylott_1-1_TNT-treated-seedlings \tab GSE18983 \tab chemical \tab Col-0 \tab seedling \tab 1.05-week-old \tab trinitrotoluen \cr
    GSM702591 \tab Skipsey_1-16_CMP_24hr_Rep1_ATH1 \tab GSE28431 \tab chemical \tab Col-0 \tab root \tab 2-week-old \tab 4-chloro-6-methyl-2-phenylpyrimidine \cr
    GSM702592 \tab Skipsey_1-17_CMP_24hr_Rep2_ATH1 \tab GSE28431 \tab chemical \tab Col-0 \tab root culture \tab 2-week-old \tab 4-chloro-6-methyl-2-phenylpyrimidine \cr
    GSM702593 \tab Skipsey_1-18_CMP_24hr_Rep3_ATH1 \tab GSE28431 \tab chemical \tab Col-0 \tab root culture \tab 2-week-old \tab 4-chloro-6-methyl-2-phenylpyrimidine \cr
    GSM4693848 \tab Columbia0_Pakerine_1 \tab GSE155026 \tab chemical \tab Col-0 \tab seedling \tab 8-day-old \tab pakerine \cr
    GSM4693849 \tab Columbia0_Pakerine_2 \tab GSE155026 \tab chemical \tab Col-0 \tab seedling \tab 8-day-old \tab pakerine \cr
    GSM4693850 \tab Columbia0_Pakerine_3 \tab GSE155026 \tab chemical \tab Col-0 \tab seedling \tab 8-day-old \tab pakerine \cr
    GSM952979 \tab seedling culture at 4 h_EtOH_biol rep 1 \tab GSE38965 \tab chemical \tab Col-0 \tab seedling \tab 9-day-old \tab ethanol \cr
    GSM952980 \tab seedling culture at 4 h_EtOH_biol rep 2 \tab GSE38965 \tab chemical \tab Col-0 \tab seedling \tab 9-day-old \tab ethanol \cr
    GSM952981 \tab seedling culture at 4 h_EtOH_biol rep 3 \tab GSE38965 \tab chemical \tab Col-0 \tab seedling \tab 9-day-old \tab ethanol \cr
    GSM1375895 \tab NAA +AA 1 \tab GSE57140 \tab chemical \tab Col-0:LUC \tab seedling \tab 13-day-old \tab 1-naphthaleneacetic acid + antimycin \cr
    GSM1375896 \tab NAA +AA 2 \tab GSE57140 \tab chemical \tab Col-0:LUC \tab seedling \tab 13-day-old \tab 1-naphthaleneacetic acid + antimycin \cr
    GSM1375897 \tab NAA +AA 3 \tab GSE57140 \tab chemical \tab Col-0:LUC \tab seedling \tab 13-day-old \tab 1-naphthaleneacetic acid + antimycin \cr
    GSM1541887 \tab LTN006-2, aza-dc Rep1 \tab GSE63131 \tab chemical \tab Col-0 \tab seedling \tab 16-day-old \tab 5-aza-2'-deoxycytidine; cold \cr
    GSM1541888 \tab LTN007-2, aza-dc Rep2 \tab GSE63131 \tab chemical \tab Col-0 \tab seedling \tab 16-day-old \tab 5-aza-2'-deoxycytidine; cold \cr
    GSM1541889 \tab LTN008-2, aza-dc Rep3 \tab GSE63131 \tab chemical \tab Col-0 \tab seedling \tab 16-day-old \tab 5-aza-2'-deoxycytidine; cold \cr
    GSM157382 \tab Sakakibara_1-1_TZ-treatment-wild_Rep1_ATH1 \tab GSE6832 \tab chemical \tab Col-0 \tab shoot \tab 2-week-old \tab trans-zeatin \cr
    GSM157383 \tab Sakakibara_1-2_TZ-treatment-wild_Rep2_ATH1 \tab GSE6832 \tab chemical \tab Col-0 \tab shoot \tab 2-week-old \tab trans-zeatin \cr
    GSM157384 \tab Sakakibara_1-3_TZ-treatment-wild_Rep3_ATH1 \tab GSE6832 \tab chemical \tab Col-0 \tab shoot \tab 2-week-old \tab trans-zeatin \cr
  }
}
\keyword{datasets}