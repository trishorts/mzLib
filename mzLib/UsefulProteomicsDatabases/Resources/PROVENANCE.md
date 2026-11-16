# Vendored ontology snapshots

These files are pinned third-party controlled vocabularies, embedded into
`UsefulProteomicsDatabases.dll` and read by `ControlledVocabulary`.

They are deliberately **pinned**, not downloaded: resolving a term against a moving target means
the same experiment annotated a year apart gets different accessions, and a corpus built that way
silently stops joining to itself. Refreshing one is a reviewable change to this repository.

| File | Source | Version | Fetched |
|---|---|---|---|
| `psi-ms.obo` | https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo | `data-version: 4.1.258` | 2026-08-07 |
| `pride_cv.obo` | https://raw.githubusercontent.com/PRIDE-Utilities/pride-ontology/master/pride_cv.obo | `data-version: releases/2026-06-19` | 2026-08-07 |

## Refreshing

```bash
curl -sL -o psi-ms.obo   https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo
curl -sL -o pride_cv.obo https://raw.githubusercontent.com/PRIDE-Utilities/pride-ontology/master/pride_cv.obo
```

Then update the versions above and **run the tests**. `TestControlledVocabulary` pins specific
accession-to-name mappings on purpose: if a refresh changes what `Q Exactive` or `HCD` resolves to,
that must surface as a failing test and a deliberate decision, not as a silent shift in every file
written afterwards. `AgreesWithTheMzmlReadersOwnAccessionTable` additionally checks the snapshot
still agrees with `Readers.Mzml.DissociationDictionary`.

## Licence

Both are licensed **CC BY 4.0** (https://creativecommons.org/licenses/by/4.0/), which requires
attribution to be conveyed to downstream recipients. Because mzLib ships them inside a NuGet
package, that attribution lives in `THIRD-PARTY-NOTICES.md` at the repository root, which is
included in the package.

# EBI's public Aspera key

`aspera_tokenauth_id_rsa` is the key EMBL-EBI publishes for **anonymous** Aspera downloads from its public
archives, PRIDE among them (EBI knowledge-base article KB0011565). It is embedded into
`UsefulProteomicsDatabases.dll` and written out, for the length of one download, by
`PrideArchiveClient.DownloadFileViaAsperaAsync`, because `ascp` reads its key from a file.

It is a private key in form only. EBI hands the same key to everyone, so it opens nothing that is not already
public, and it is **not a secret**: committing it leaks nothing, and a secret scanner that flags it is mistaken.
pridepy, EBI's own PRIDE client, redistributes it the same way.

| File | Source | Version | Fetched | Size | sha256 |
|---|---|---|---|---|---|
| `aspera_tokenauth_id_rsa` | https://raw.githubusercontent.com/PRIDE-Archive/pridepy/fbd30d43366c69a93a11600f09bb46fe33f96a93/pridepy/aspera/key/aspera_tokenauth_id_rsa | pridepy commit `fbd30d43` (2025-12-22, the key's last change) | 2026-11-16 | 1,679 B | `4f90eaab1230f9793d357ea7669ff9205edd30864a0bcf5ef3eb42a2720de1f1` |

The key is RSA. EBI's older DSA key (`asperaweb_id_dsa.openssh`, still bundled with Aspera Connect) and Connect's own
ascp 4.4.3 both fail to authenticate against EBI's server now; ascp 4.4.4 or newer with this key works on Windows and
Linux (checked 2026-10-11 and 2026-11-03).

## Refreshing

```bash
curl -fsSL -o aspera_tokenauth_id_rsa https://raw.githubusercontent.com/PRIDE-Archive/pridepy/master/pridepy/aspera/key/aspera_tokenauth_id_rsa
sha256sum aspera_tokenauth_id_rsa
```

Refresh only when EBI changes its key, which shows as every Aspera download failing with
`PrideAsperaFailure.LoginRefused`. That symptom alone is not proof: EBI refused this unchanged key for about three days
from 2026-10-30. Compare the SHA-256 with the row above before replacing anything.

## Licence

pridepy is licensed **Apache-2.0** (https://www.apache.org/licenses/LICENSE-2.0); its notice is in
`THIRD-PARTY-NOTICES.md` at the repository root.
