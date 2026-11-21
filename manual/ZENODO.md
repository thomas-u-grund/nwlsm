# Zenodo record for the nwlsm manual

Metadata ready to paste into a **new** Zenodo record (the first for this manual; it is not a new version of the nwsaom or nwergm record). Nothing has been uploaded.

**DOI.** None yet. On zenodo.org start a new upload, press *Reserve DOI*, write the reserved all-versions DOI into `\manualdoi` in `nwlsm_manual.tex` (it is empty now, and an empty DOI prints nothing), rebuild with `make -C manual pdf` and upload the rebuilt PDF to the reserved record. The citation text always uses the all-versions DOI (master RELEASING.md B.10).

**Labelling.** The title and the citation of the manual contain no version number. The Zenodo version label is the series "X.Y" (here "1.0"; internal to Zenodo, not part of the citation), and the description says that the manual describes nwlsm X.Y.x. The PDF title page has one line saying which series it describes (`\manualseries`); the exact package version of the examples is `\pkgversion`. The software version is cited separately, as printed by `nwlsm version`.

**Patch releases.** The manual of the series 1.0 holds for 1.0.x; a patch release changes no documented option or number. Title and citation stay unchanged; no new Zenodo version is needed for a patch. A new Zenodo version only for 1.4 or 2.0 with changed documented behaviour (then change `\manualseries`).

## Metadata

**Title.** nwlsm: Latent space models for networks in Stata. Manual

**Version.** 1.0 (the series label, internal to Zenodo: the manual describes nwlsm 1.0.x; not part of the citation)

**Publication date.** The day of the upload.

**Resource type.** Publication, sub-type Software documentation (as the nwsaom and nwergm records).

**Creators.** Grund, Thomas U. (RWTH Aachen University), ORCID 0000-0001-6751-5429 (as on the other records; remove it in the form if the author does not want it).

**License.** Creative Commons Attribution 4.0 International (CC-BY-4.0), as the other manual records. The code stays under GPL-3.0-only.

**Access.** Open.

**Description** (plain text, one paragraph):

> This is the user manual and reference of nwlsm, series 1.0 (examples produced with version 1.0.4 and the nwcommands core 3.33.0), a Stata package for latent space models of networks: models in which every actor has a position in a latent space of one to ten dimensions and the probability of a tie depends on the distance between the positions (Euclidean model) or on their inner product (bilinear model). nwlsm builds on the nwcommands core package and follows the R package latentnet (Krivitsky and Handcock) in its model, priors, minimum Kullback-Leibler positions, latent clusters and BIC, and the latent space model of Hoff, Raftery and Handcock; it adds the latent space model for social influence of Sweet and Adhikari. The posterior is sampled by MCMC in a compiled plugin, with parallel chains. The manual describes the installation, the model and its priors, the data and the syntax, every option, worked examples (the karate club, Sampson's monks, the Lazega law firm, and simulated data for the influence model) for covariates, random sender, receiver and sociality effects, the bilinear term, latent clusters and the choice of their number by BIC, valued networks (Poisson, binomial, normal) and social influence, post-estimation, diagnostics and goodness of fit, the plot options (name, saving, export), the comparison of models, the validation against latentnet and by simulation, and the limits of the package. nwlsm requires the nwcommands core package (version 3.29.0 or later) and Stata 17 or later and is meant for researchers who work in Stata. It is not a replacement for latentnet, which remains the reference implementation, and it is not affiliated with or endorsed by its authors. The software is free software under the GNU General Public License, version 3. The package and the manual were prepared with the help of an AI coding assistant working under the author's direction, as stated in the manual. Please cite this manual (all-versions DOI) and report the exact version printed by nwlsm version.

**Keywords.** Stata; nwcommands; nwlsm; latent space model; latent position model; social network analysis; Bayesian; MCMC; latent clusters; social influence; latentnet

**Language.** English

**Related identifiers.** (relation, identifier, type)

| Relation | Identifier | Type |
|---|---|---|
| Is supplement to | https://github.com/thomas-u-grund/nwlsm/releases/tag/v1.0.4 | URL |
| Is supplement to | https://github.com/thomas-u-grund/nwlsm | URL |
| Requires | https://github.com/thomas-u-grund/nwcommands | URL |
| Is related to (manual of the sibling package nwsaom) | 10.5281/zenodo.23623807 | DOI |
| Is related to (manual of the sibling package nwergm) | 10.5281/zenodo.23661331 | DOI |

The first identifier points to the GitHub release v1.0.4: confirm that this release is published before the upload (if the package has a later 1.0.x release by then the manual still holds; use the tag you want to point to). The two sibling links are optional.

**Other fields.** Funding: none. Communities: none.

## Steps for the owner

1. Read the PDF once (`manual/nwlsm_manual.pdf`, 60 pages); the examples were produced with Stata 17 on macOS with version 1.0.4 of the package and core 3.33.0.
2. On zenodo.org create the new record, reserve the DOI, put it into `\manualdoi`, rebuild (`make -C manual pdf`), upload the PDF and fill the fields above.
3. Commit the rebuilt PDF (`git add -f manual/nwlsm_manual.pdf`).
4. Afterwards: the DOI goes into the citation texts (README of the package, `CITATION.cff`, the website page *How to cite* of nwlsm and the manual link of its menu), as for nwergm; none of these is changed yet.
