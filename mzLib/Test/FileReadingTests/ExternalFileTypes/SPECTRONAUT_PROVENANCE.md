# Spectronaut test fixtures: provenance

These files are cut from public MassIVE deposits. Every line is the deposited bytes, CRLF included
(`.gitattributes` marks them `-text`). The header line and a selection of whole data lines were kept; nothing
was edited. Two files were renamed so they end with their type's `GetFileExtension()`, and the source column
gives each original name.

The cutting script is `tools/cut_fixtures.py` in the readers project (`E:\CodeReview\readers`). Rerunning it
reproduces these bytes from the sources below.

Sources were fetched from `massive-ftp.ucsd.edu` (FTPS). For the multi-GB exports, only byte ranges were
fetched: the first 64 KB, or 64 KB at a stated offset. Their SHA-256 is of the fetched bytes.

| Fixture | Rows | Bytes | SHA-256 | Source | Source SHA-256 | Selection |
|---|---|---|---|---|---|---|
| `Spectronaut19_PG_Report.tsv` | 80 | 17225 | `bb0dd79dd0ba3780` | MSV000096226 quant/20241010_162601_220424_AM_MouseBrain_DR_WholeProteome_Report.tsv | `70ffec0a22f23cf4` | Every row of 8 protein groups, file order kept: A2A432 = plain single protein; O55074;Q7TN79 = collapsed gene list (2 members, 1 gene); O08580;Q61539 = collapsed existence list (2 members, 1 value); A2A5Y0;Q61765;Q9D646;Q9D646 = duplicate id in PG.ProteinGroups, deduplicated in PG.ProteinAccessions, gene list with empty entries; P01675;P01676;P01677;P01678;P01679;P04940;P04941;P04944 = 8 members, no genes; A2A5Z6 = NaN in every run (0 precursors quantified); P01631 = blank PG.Genes; P00761 = blank PG.Protein Existence (also blank genes) |
| `Spectronaut19_ConditionSetup.tsv` | 10 | 790 | `1f264f852500c74d` | MSV000096226 metadata/220424_AM_MouseBrain_DR_WholeProteome_ConditionSetup.tsv | `1f264f852500c74d` | The whole file. |
| `SpectronautPrecursor_pep_Report.tsv` | 60 | 10984 | `1f4d7488e2b1fe9b` | MSV000098115 search/20250610_150722_240320_SD374566_DUBactivity_Report_pep.tsv (first 64 KB) | `588b46fb10fe96bd` | The first 40 rows, plus every head row with a modification or NaN PG.Quantity, capped at 70. File order is kept. |
| `Spectronaut10_fragment_Report.csv` | 30 | 13774 | `2f3df5c8918265eb` | MSV000087085 quant/20181105_Aging_stoich_Spectronaut10_CR_reprocess_Report.csv (first 64 KB) | `f2b0f4868e82371b` | Every fragment row of 5 precursors (FG.Id): the first 4 and the first with a [+O] delta. _K[Acetyl (K)]K[Acetyl (K)]LLAK[Acetyl D3 (K)]DLR_.2; _K[Acetyl (K)]K[Acetyl D3 (K)]LLAK[Acetyl D3 (K)]DLR_.2; _K[Acetyl D3 (K)]K[Acetyl (K)]LLAK[Acetyl D3 (K)]DLR_.2; _K[Acetyl D3 (K)]K[Acetyl D3 (K)]LLAK[Acetyl D3 (K)]DLR_.2; _K[Acetyl D3 (K)]K[Acetyl (K)]DLK[Acetyl (K)]M[+O]SVLISK[Acetyl D3 (K)]LNDTE_.3 |
| `Spectronaut10_quoted_Report.csv` | 8 | 4251 | `f570e2723a354131` | MSV000087085 quant/20181105_AgingStudy_Spectronaut_reprocess2_WTCR5KOCR5_2of4.csv (first 64 KB) | `63cbc2b76e9d1add` | The first 8 data rows. |
| `Spectronaut_PG_Report.xls` | 18 | 2562 | `019678d4ac64a98a` | MSV000094212 quant/20220626_220231_220625_intestine_regeneration_distal_Mclass_finalLib_Report.xls (first 64 KB) | `494abe3a25117cc5` | The first 15 rows, plus the first 3 rows with NaN PG.Quantity. File order is kept. |
| `Spectronaut_TRICReport.tsv` | 63 | 20609 | `febf162b894ceb0f` | MSV000090777 quant/Results/120min/20220513_165438_2022_0513_SSL2_ND_InjvsNon_120min_Spec16_Library_v01/2022_0513_SSL2_ND_InjvsNon_120min_Spec16_Library_v01_TRICReport.tsv (head, plus byte ranges at size*3/16 and size*5/16) | `d311c4346ab53171` | Whole transition groups (every rank row of one transition_group_id). From the head (run 0): the first 5 groups and the first with Oxidation. From each slice (runs 1 and 2): the first 2 decoy groups and the first 3 target groups. Slice-edge groups are dropped. Sources: head, MSV000090777_TRICReport_k3_off111588926.slice64k, MSV000090777_TRICReport_k5_off185981543.slice64k |
| `Spectronaut_PG_Pivot_Report.tsv` | 19 | 3521 | `031fd752c5739417` | MSV000090777 quant/Results/120min/20220513_165438_2022_0513_SSL2_ND_InjvsNon_120min_Spec16_Library_v01/2022_0513_SSL2_ND_InjvsNon_120min_Spec16_Library_v01_Report_Protein_Quant_Pivot.xls (whole file; same bytes, renamed .tsv) | `60b4f16b1df55532` | The first 15 rows, plus the first 5 rows with a Filtered cell and the first 5 with a NaN cell. File order is kept. |
| `Spectronaut_custom_Pivot_Report.csv` | 10 | 6774 | `4d36f1922d595316` | MSV000090777 quant/Results/120min/20220513_165438_2022_0513_SSL2_ND_InjvsNon_120min_Spec16_Library_v01/2022_0513_SSL2_ND_InjvsNon_120min_Spec16_Library_v01_Report_Birgit.csv (first 64 KB) | `e99abfd6a27b15c0` | The first 10 data rows. |
| `Spectronaut_PTMSiteReport.tsv` | 23 | 7378 | `8c532ddde5ffd41e` | MSV000090777 quant/Results/120min/20220722_135200_22_0722_SSL2_ND_120min_Inj_Library_v01/22_0722_SSL2_ND_120min_Inj_Library_v01_-PTMSiteReport.tsv (first 64 KB) | `297dde1f1b1b5216` | The first 20 rows, plus the first 3 Oxidation (M) sites. File order is kept. |

Datasets: MSV000096226 (= PXD057243), MSV000098115, MSV000087085, MSV000094212 and MSV000090777, all public at MassIVE.
