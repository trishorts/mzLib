# mzTab test fixtures: where each one comes from

These are the fixtures for `MzTabFile` (`TestMzTabFile`, `TestSupportedFileExtensions`). Each is real output from a
search or quantification tool, as published. Each excerpt is a selection of whole original lines, in source order,
with the source's line endings: the MaxQuant file is CRLF and the others LF. `.gitattributes` marks them `-text`, so
git keeps those bytes.

An excerpt keeps every line before the PSM header except PRT rows past the first three (and any PEP rows), the PSH
header, and the chosen PSM rows. A chosen row brings with it the adjacent rows that share its PSM_ID, so no multi-row
PSM is cut in half.

**The one edit:** a username inside a local path (`c:/users/<name>/`, `/Users/<name>/`) is replaced with `user`. That
keeps a submitter's or developer's account name out of mzLib. Only the MaxQuant `ms_run[n]-location` lines, the
jmzTab `custom[2]` lines and the mzML's `<sourceFile>` line changed. `BSA1_F1.mzML` is an indexed mzML, whose index
gives each spectrum's byte offset, so its replacement is `username`, the same length as the name it replaces, and
every offset stays valid.

| Fixture | Bytes | SHA-256 | Source | Source SHA-256 | Licence |
|---|---|---|---|---|---|
| `OpenMS3.7_ProteomicsLFQ_1_out.mzTab` | 65,324 | `1d39a3457e82e27f5e5cdc90400a48e7f72ebe2b654bab6033a7748a93a92761` | OpenMS `src/tests/topp/ProteomicsLFQ_1_out.mzTab` at `fbb25191d7659e9aa0326456cc187b7f9119409e`, whole | identical | BSD-3-Clause (OpenMS) |
| `BSA1_F1.mzML` | 5,788,763 | `3e66d6b29e7dd923d9eb3a81e95700c7d2ed2805b357434c008ab9838cf11b08` | OpenMS `share/OpenMS/examples/FRACTIONS/BSA1_F1.mzML` at the same commit, whole but for the scrubbed path: the spectra of the mzTab's `ms_run[1]` | `c3eb9b13f423d89997d14ca68470951026a387afa3959b7dd1c59ceb9c16b2f4` | BSD-3-Clause (OpenMS) |
| `MaxQuant1.6.3.1_PXD011194_excerpt.mzTab` | 34,214 | `a03234f015ea7a2e3e8c6e04ced3f392c5e14cffaaebe82e02886dda3ef742f5` | PRIDE PXD011194 `mzTab.mztab.gz` (MaxQuant's own), `https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/09/PXD011194/` | `681765693b94f81e730bcaee7a86a87b2432afdf633cf7a38a980476a70902b9` (decompressed) | CC0 (PRIDE Archive) |
| `ProteomeDiscoverer2.5_MSV000096151_excerpt.mzTab` | 11,388 | `15be289a6ce2ca1979ce2d4cb164c69ee769b6bb95f971fe3a2698a03787c5a3` | MassIVE MSV000096151 `ccms_result/rs240609_1505_KY_10plex_01.mzTab` (Proteome Discoverer 2.5, through MassIVE's conversion) | `f12ec3a4470c7145df29a66fb621749180c7394ded0ed01623c1bde36d31a145` | CC0 1.0 (MassIVE) |
| `MassIVE_MSV000078886_excerpt.mzTab` | 5,951 | `585928fe64f83d9702d34fc93a1676dde3bcc8263b3cd434ded8212bfef66a58` | MassIVE MSV000078886 `result/RESULT/results.mzTab` (MassIVE's TSV-to-mzTab conversion) | `a4b5ba3f368d552c54c6af1bd6ab4e2699f0196bee92ab9dce8e9b21532c55cb` | CC0 1.0 (MassIVE) |
| `Mascot2.3_jmzTab_small_excerpt.mzTab` | 12,334 | `80e4d13b05e531e7594e5ca5e67e9673e1dd001e22e1c9c14cb5a965cb2caee9` | jmzTab `jmztab-parser/src/test/resources/testset/mzTabs/small.mztab` at `db11063f2026d6fe33502b8983aee9b083e3a7ba` (Mascot 2.3.01, converted from mzIdentML) | `bbd044476eea67f102b5988ab29deefde9818535d0c499e213e8938c1dea3f96` | Apache-2.0 (jmzTab, PRIDE-Utilities, EMBL-EBI) |
| `PrideConverter_jmzTab_PRIDE_Example.mzTab` | 63,092 | `340270a8deee42a299b30ffe0ef34f1c99567e62bded669098bc471e8c7c9a04` | jmzTab `jmztab-model/src/test/resources/testset/PRIDE_Example.mztab` at `6d544578ea420b4fab6fffebc8e00544445f2226`, whole but for the scrubbed line | `a23742162910244bbf52ae5ed68e00282095c8f004cce2e79669b3cc46cc1560` | Apache-2.0 (jmzTab, PRIDE-Utilities, EMBL-EBI) |
| `pyteomics_mzTabM_test_mztab2.mzTab` | 8,381 | `f2f7a95b05676c9b5c69531a7289e5f3adbbc28269710318613564d6d59b8d7f` | pyteomics `tests/test_mztab2.mztab` at `024ae2e75aba6776ca0d4be3a06334173b694fd7`, whole (mzTab-M 2.0.0-M) | identical | Apache-2.0 (pyteomics) |

## The PSM rows each excerpt keeps

Row numbers count PSM lines in the source, starting at 1.

- **MaxQuant (85 rows, 49 PSMs).** Rows 1-40, which hold PSM_ID 1 (`ms_run[16]:index=3508` on two proteins; opt scan
  995) and the first CHEMMODs. Then: 65, a Cys peptide with no modifications under the fixed Carbamidomethyl; 363, a
  Cys peptide with an Oxidation; 448, the N-terminal acetyl `0-CHEMMOD:42.0105646863`; 1611, a `, `-joined cell; 3655,
  a `:` in an accession; 6933 and 17291, `REV__` decoys; 37266, a `-` pre. Each is widened to its PSM_ID's rows.
- **Proteome Discoverer (30 rows, 17 PSMs).** Rows 1-12, where PSM_ID 4 is rows 3-5 on three proteins. Then the two
  later exact copies of PSM_ID 4 (source rows 468-470 and 888-890), with their neighbours (466-472, 886-892).
- **MassIVE (41 rows, 41 PSMs).** Rows 1-31: nominal `CHEMMOD:-n` and `CHEMMOD:+n`, `-` termini, pyro-Glu `1-UNIMOD:28`.
  Then 48 and 87, `XXX_` decoys in a file with no decoy column; 788-793, PSM_ID 790 written for two spectra on adjacent
  rows; 1083-1084, PSM_ID 1082 likewise.
- **Mascot (36 rows, 36 PSMs).** Rows 1-30: UNIMOD and neutral-loss params, bare-number spectra_refs. Then 86-88, a
  `-` pre; 409-411, an `M` pre.

A script cut the excerpts from the sources by these rules. Applying the same rules to the source files (the checksums
above) gives the same bytes.

## Licences

- OpenMS: BSD-3-Clause, https://github.com/OpenMS/OpenMS/blob/develop/LICENSE
- PRIDE PXD011194: CC0, as the PRIDE Archive v3 API reports for the project.
- MassIVE MSV000096151 and MSV000078886: CC0 1.0, as each dataset page states.
- jmzTab: Apache License 2.0, https://github.com/PRIDE-Utilities/jmztab/blob/master/LICENSE. The two files are
  redistributed unmodified apart from the excerpting and the scrubbed path described above.
- pyteomics: Apache License 2.0, https://github.com/levitsky/pyteomics/blob/master/LICENSE. Redistributed unmodified.
