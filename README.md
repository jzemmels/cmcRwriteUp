# cmcRwriteUp
Write-up for cmcR package

For the data used in this paper, see: [https://github.com/CSAFE-ISU/cartridgeCaseScans/](https://github.com/CSAFE-ISU/cartridgeCaseScans/)

Below is an explanation of non-standard files/folders included with the submission.

- data folder: contains data to demonstrate the cmcR package and reproduce various results discussed in the paper. A description of each file in this folder follows.
  - fadul1-1.x3p, fadul1-2.x3p, fadu2-1.x3p: Two same-source and one different-source cartridge case scans. They are used to demonstrate usage of the cmcR package.
  - cmcCountData.Rdata file: contains CMC count data from the Fadul et al. (2011) set of 40 cartridge case scans. These data are used to generate Figures 12 and 13.
  - cmcCountData_script.R file: a script to reproduce the cmcCountData.RData file
- images folder: contains manually-created images and photos used in the paper
- figures folder: contains automatically-generated images (while running zemmels-vanderplas-hofmann.R) used in the paper
