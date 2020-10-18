## ----localDataDir, include=FALSE----------------------------------------------
if(!dir.exists("data")){
  dir.create("data")
}
if(!file.exists("data/fadul1-1.x3p")){
  library(dplyr) # pipe not defined yet
  download.file("https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/DownloadMeasurement/2d9cc51f-6f66-40a0-973a-a9292dbee36d", destfile = "data/fadul1-1.x3p", mode = "wb")
}
if(!file.exists("data/fadul1-2.x3p")){
  download.file("https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/DownloadMeasurement/cb296c98-39f5-46eb-abff-320a2f5568e8", destfile = "data/fadul1-2.x3p", mode = "wb")
}
if(!file.exists("data/fadul2-1.x3p")){
  download.file("https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/DownloadMeasurement/8ae0b86d-210a-41fd-ad75-8212f9522f96", destfile = "data/fadul2-1.x3p", mode = "wb")
}


## ---- derivativeImagesDir,include=FALSE---------------------------------------
if(!dir.exists("derivatives")){
  dir.create("derivatives")
}


## ----setup,echo=FALSE,message=FALSE,warning=FALSE-----------------------------
knitr::opts_chunk$set(cache = T, dpi = 300, fig.width = 8, fig.height = 4, out.width = "\\textwidth", dpi = 300)
library(cmcR) # remotes::install_github("CSAFE-ISU/cmcR")
library(tidyverse)
library(x3ptools)
library(rgl)
set.seed(4132020)


## ----eval=FALSE,echo=TRUE-----------------------------------------------------
#> library(cmcR)
#> set.seed(4132020)
#> 
#> nrbtd_url <- "https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/"
#> 
#> fadul1.1_id <- "DownloadMeasurement/2d9cc51f-6f66-40a0-973a-a9292dbee36d"
#> fadul1.2_id <- "DownloadMeasurement/cb296c98-39f5-46eb-abff-320a2f5568e8"
#> fadul2.2_id <- "DownloadMeasurement/8ae0b86d-210a-41fd-ad75-8212f9522f96"
#> 
#> # Code to download breech face impressions:
#> nrbtd_url <- "https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/"
#> download.file(
#>   file.path(nrbtd_url, fadul1.1_id), destfile = "data/fadul1-1.x3p", mode = "wb")
#> download.file(
#>   file.path(nrbtd_url, fadul1.2_id), destfile = "data/fadul1-2.x3p", mode = "wb")
#> download.file(
#>   file.path(nrbtd_url, fadul2.1_id), destfile = "data/fadul2-1.x3p", mode = "wb")


## ----echo=FALSE,fig.cap='\\label{fig:ccPair_separated} A cartridge case pair with visible breech face impressions under a microscrope. Forensic practitioners determine by visual comparison how similar these impressions are to conclude whether the pair is a match. It can be especially challenging to determine whether a marking is an \\emph{individual characteristic} caused by contact with the breech face of a particular firearm or a \\emph{subclass characteristic} shared across cartridge cases manufactured by the same manufacturer within a short time frame \\citep{firearmTraining}.',fig.align='center',fig.pos='htbp',out.width="\\textwidth"----
knitr::include_graphics("images/cartridgeCasePair_separated.PNG")


## ----echo=FALSE,fig.cap='\\label{fig:ccPair_combined} A cartridge case pair with visible breech face impressions under a microscrope.  A thin line can be seen separating the two views. The degree to which the markings coincide is used to conclude whether the pair comes from the same source.',fig.pos='htbp',out.width="\\textwidth"----
knitr::include_graphics("images/cartridgeCasePair_comparison_with_line.PNG")


## ---- fadul1-1Screenshot,include=FALSE----------------------------------------
fadul1.1 <- x3ptools::read_x3p("data/fadul1-1.x3p")

#apply low-pass filter to reduce noise in scan:
surface1 <- fadul1.1 %>%
  cmcR::preProcess_gaussFilter(wavelength = 16,filtertype = "lp")

surface1 <- surface1$surface.matrix

params <- rgl::r3dDefaults

zoom <- .7
size <- c(300,300)

params$windowRect <- c(40, 125, 40 + size[1], 125 + size[2])
params$userMatrix <- diag(c(1, 1, 1, 1))
params$zoom <- zoom

#for some reason the first rgl device opened doesn't plot anything, but
#subsequent devices do...
open3d(params = params)
rgl.close()

#opens blank "canvas" upon which we can add lights, surfaces, etc.
open3d(params = params)

#removes any previously declared lights in scene
rgl.pop("lights")

#set-up two lights for scene -- a lot of experimentation possible here
light3d(x = -1,y = 1,z = 2,viewpoint.rel = TRUE,ambient = "white",diffuse = "white",specular = "white")
light3d(x = 0,y = 0,z = 10,ambient = "grey60",diffuse = "grey50",specular = "grey60",viewpoint.rel = TRUE)

#setup surface visualization
multiply <- 1 #x3ptools::image_x3p default to exaggerate relief
z <- multiply * surface1 # Exaggerate the relief
yidx <- ncol(z):1
y <- fadul1.1$header.info$incrementY * yidx
x <- fadul1.1$header.info$incrementX * (1:nrow(z))

# emission, specular, ambient affect how the surface interacts with lights --
# again, a lot of possible experimentation
surface3d(x, y, z, back = "filled",emission = "grey30",specular = "grey50",ambient = "grey10")

x3ptools::x3p_snapshot(file = "derivatives/fadul1-1.png")

rgl.close()


## ----fadul1-2Screenshot,include=FALSE-----------------------------------------
fadul1.2 <- x3ptools::read_x3p("data/fadul1-2.x3p")

surface2 <- fadul1.2 %>%
  cmcR::preProcess_gaussFilter(wavelength = 16,filtertype = "lp")
#opens blank "canvas" upon which we can add lights, surfaces, etc.
open3d(params = params)

surface2 <- surface2$surface.matrix

#removes any previously declared lights in scene
rgl.pop("lights")

#set-up two lights for scene -- a lot of experimentation possible here
light3d(x = -1,y = 1,z = 2,viewpoint.rel = TRUE,ambient = "white",diffuse = "white",specular = "white")
light3d(x = 0,y = 0,z = 10,ambient = "grey60",diffuse = "grey50",specular = "grey60",viewpoint.rel = TRUE)

#setup surface visualization
multiply <- 1 #x3ptools::image_x3p default to exaggerate relief
z <- multiply * surface2 # Exaggerate the relief
yidx <- ncol(z):1
y <- fadul1.2$header.info$incrementY * yidx
x <- fadul1.2$header.info$incrementX * (1:nrow(z))

# emission, specular, ambient affect how the surface interacts with lights --
# again, a lot of possible experimentation
surface3d(x, y, z, back = "filled",emission = "grey30",specular = "grey50",ambient = "grey10")

x3ptools::x3p_snapshot(file = "derivatives/fadul1-2.png")

rgl.close()


## ---- rawBFs,echo=FALSE,fig.cap='\\label{fig:cartridgeCasePair} Unprocessed surface matrices of the known-match Fadul 1-1 (left) and Fadul 1-2 (right) \\citep{fadul_empirical_nodate}. The observations in the corners of these surface matrices are artifacts of the staging area in which these scans were taken. The holes on the interior of the primer surfaces are caused by the firing pin striking the primer during the firing process. The region of the primer around this hole does not come into uniform contact with the breech face of the firearm.', fig.subcap=c('',''),fig.align='center',fig.pos='htbp',out.width='.49\\linewidth',out.height='.49\\linewidth'----
knitr::include_graphics(c("derivatives/fadul1-1.png","derivatives/fadul1-2.png"))


## ----load-data, include = F, cache = T----------------------------------------

fadul1.1 <- x3ptools::read_x3p("data/fadul1-1.x3p") %>%
  cmcR::preProcess_cropBFExterior(radiusOffset = -30,
                                  agg_function = median) %>%
  cmcR::preProcess_filterBFInterior(radiusOffset = 200) %>%
  cmcR::preProcess_removeBFTrend(statistic = "quantile",
                                 tau = .5,
                                 method = "fn") %>%
  cmcR::preProcess_gaussFilter() %>%
  x3ptools::sample_x3p()

fadul1.2 <- x3ptools::read_x3p("data/fadul1-2.x3p") %>%
  cmcR::preProcess_cropBFExterior(radiusOffset = -30,
                                  agg_function = median) %>%
  cmcR::preProcess_filterBFInterior(radiusOffset = 200) %>%
  cmcR::preProcess_removeBFTrend(statistic = "quantile",
                                 tau = .5,
                                 method = "fn") %>%
  cmcR::preProcess_gaussFilter() %>%
  x3ptools::sample_x3p()


## ----cmc-ccf, include = F, cache = T------------------------------------------
kmComparison <- cmcR::cellCCF_bothDirections(x3p1 = fadul1.1,
                                             x3p2 = fadul1.2,
                                             cellNumHoriz = 8,
                                             minObservedProp = .1,
                                             regionToCellProp = 9,
                                             ccfMethod = "fftThenPairwise")

kmCMC <- cmcR::cmcFilter_improved(kmComparison,
                                    ccf_thresh = .5,
                                    dx_thresh = 20,
                                    theta_thresh = 6,
                                    missingTheta_decision = "fail")


## ----cmc-filter, include = F, cache = T---------------------------------------
# if (!file.exists("data/kmcmc.Rdata")) {
#   save(kmCMC, file = "data/kmcmc.Rdata")
# } else {
#   load("data/kmcmc.Rdata")
# }


## ----include = F,cache = T----------------------------------------------------
# fadul1.1_downsampled <- x3ptools::x3p_sample(fadul1.1,m = 2)
# 
# fadul1.1_ransac <- fadul1.1_downsampled
# fadul1.1_ransac$surface.matrix <- fadul1.1_ransac$surface.matrix %>%
#   cmcR::preProcess_ransac(ransacInlierThresh = 1e-5) %>%
#   cmcR::preProcess_levelBF(useResiduals = TRUE)
# 
# fadul1.1_cropped <- fadul1.1_ransac
# fadul1.1_cropped$surface.matrix <- fadul1.1_cropped$surface.matrix %>%
#   cmcR::preProcess_cropWS(croppingThresh = 1)
# fadul1.1_cropped$header.info$sizeY <- ncol(fadul1.1_cropped$surface.matrix)
# fadul1.1_cropped$header.info$sizeX <- nrow(fadul1.1_cropped$surface.matrix)
# 
# fadul1.1_noFP <- fadul1.1_cropped
# fadul1.1_noFP$surface.matrix <- fadul1.1_noFP$surface.matrix %>%
#   cmcR::preProcess_removeFPCircle(aggregationFunction = function(x,na.rm){median(x,na.rm) - 11})
# 
# fadul1.1_filtered <- fadul1.1_noFP
# fadul1.1_filtered$surface.matrix <- fadul1.1_filtered$surface.matrix %>%
#   cmcR::preProcess_gaussFilter(res = fadul1.1_filtered$header.info$incrementY,
#                                wavelength = c(16,500))

fadul1.1_original <- x3ptools::read_x3p("data/fadul1-1.x3p")

fadul1.1_croppedExt <- cmcR::preProcess_cropBFExterior(x3p = fadul1.1_original,
                                            radiusOffset = -30,
                                            agg_function = median) 

fadul1.1_croppedInt <- cmcR::preProcess_filterBFInterior(fadul1.1_croppedExt,
                                                         radiusOffset = 200)

fadul1.1_medRemoved <-   cmcR::preProcess_removeBFTrend(fadul1.1_croppedInt,
                                                        statistic = "quantile",
                                                        tau = .5,
                                                        method = "fn")



fadul1.1_downsampled <- x3ptools::sample_x3p(fadul1.1_medRemoved,
                                             m = 2)

fadul1.1_bpFiltered<- cmcR::preProcess_gaussFilter(x3p = fadul1.1_downsampled,
                                                   wavelength = c(16,500),
                                                   filtertype = "bp")


## ---- echo = F,warning = F,message = F,cache = T,fig.cap='\\label{fig:processingPipeline} Illustration of \\pkg{cmcR} preprocessing pipeline designed to automate the manual cleaning in the CMC papers. \\svp{At each stage, the amount of variability in height across the scan decreases as extraneous sources of noise are removed.}',fig.align='center',fig.pos='htbp',out.width='\\textwidth', message = F, warning = F----

preProcessingPlot <- cmcR::x3pListPlot(list("(1) Original" = fadul1.1_original,
                                            "(2) Crop exterior/interior" = fadul1.1_croppedInt,
                                            "(3) Level surface" = fadul1.1_medRemoved,
                                            "(4) Band-pass filter" = fadul1.1_bpFiltered),
                                       type = "list",
                                       legend.quantiles = c(0,.5,1)) %>%
  map2(.x = .,
       .y = list(element_text(),element_blank(),element_blank(),element_blank()),
       .f = ~ .x + theme(legend.position = "bottom",
                         legend.title = .y) +
         ggplot2::guides(fill = ggplot2::guide_colourbar(barheight = grid::unit(.3,"in"),
                                                         barwidth = grid::unit(1.5,"in"),
                                                         label.theme = ggplot2::element_text(size = 7),
                                                         title.theme = ggplot2::element_text(size = 10),
                                                         title.position = "top",
                                                         frame.colour = "black",
                                                         ticks.colour = "black"),
                         colour = FALSE) +
         scale_fill_gradientn(colours = rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')),
                              values = scales::rescale(quantile(.x[[1]]$value,c(0,.01,.025,.1,.25,.5,.75,0.9,.975,.99,1),na.rm = TRUE)),
                              breaks = c(round(min(.x[[1]]$value*1e6,na.rm = TRUE),2),
                                         0,
                                         round(max(.x[[1]]$value*1e6,na.rm=TRUE),2)),
                              limits = c(1.01*min(.x[[1]]$value*1e6,na.rm = TRUE),
                                         1.01*max(.x[[1]]$value*1e6,na.rm=TRUE)),
                              na.value = "gray80") +
         ggplot2::labs(fill = expression("Height ["*mu*"m]")))

gridExtra::grid.arrange(preProcessingPlot$`(1) Original`,
                        preProcessingPlot$`(2) Crop exterior/interior`,
                        preProcessingPlot$`(3) Level surface`,
                        preProcessingPlot$`(4) Band-pass filter`,
                        widths = unit(c(1,1,1,1),units = "null"))


## ----include=FALSE,eval=FALSE-------------------------------------------------
#> fadul1.1 <- selectBFImpression_sample_x3p(x3p_path = paste0(nrbtd_url,fadul1.1_id),
#>                                           ransacInlierThresh = 1e-6, #.1 microns
#>                                           ransacIters = 300,
#>                                           m = 2, #sample_x3p downsample rate
#>                                           gaussFilterWavelength = c(16,500),
#>                                           gaussFilterType = "bp") #band-pass filter


## ----echo=FALSE,cache = T,fig.cap='\\label{fig:processedScans} Fadul 1-1 and Fadul 1-2 after preprocessing. Similar striated markings are now easier to visually identify on both surfaces.',fig.align='center',fig.pos='htbp',out.width='\\textwidth', message = F, warning = F----

cmcR::x3pListPlot(x3pList = list("Fadul 1-1" = fadul1.1,
                                 "Fadul 1-2" = fadul1.2),
                  # x3pList = list("Fadul 1-1" = fadul1.1$x3p,
                  #"Fadul 1-2" = fadul1.2$x3p),
                  type = "faceted",
                  rotate = 90,
                  legend.quantiles = c(0,.01,.2,.5,.8,.99,1)) +
  guides(fill = guide_colourbar(barheight = grid::unit(2.6,"inches"),
                                label.theme = element_text(size = 8),
                                title.theme = ggplot2::element_text(size = 12),
                                frame.colour = "black",
                                ticks.colour = "black")) +
  theme(legend.position = c(1.11,.551),plot.margin = ggplot2::margin(c(0,3,.2,0),unit = "cm"))


## ---- echo=FALSE,fig.cap='\\label{fig:cmc_illustration} Illustration of comparing a "cell" in the reference cartridge case scan (left) to a larger region in a questioned cartridge case scan (right). Every cell in the reference cartridge case is similarly paired with a region in the questioned. To determine the rotation at which the two cartridge cases align, the cell/region pairs are compared for various rotations of the questioned cartridge case.',fig.align='center',fig.pos='htbp',out.width='.75\\textwidth'----

knitr::include_graphics("images/cmc_illustration.PNG")


## ----echo=TRUE,eval=FALSE-----------------------------------------------------
#> kmComparison <- cellCCF_bothDirections(x3p1 = fadul1.1,
#>                                        x3p2 = fadul1.2,
#>                                        thetas = seq(-30,30,by = 3),
#>                                        cellNumHoriz = 8,
#>                                        cellNumVert = 8,
#>                                        minObservedProp = .1)


## ----include=FALSE,eval=FALSE-------------------------------------------------
#> kmComparison$comparison_1to2$ccfResults %>%
#>   topResultsPerCell() %>%
#>   head()

## ----echo=FALSE,warning=F,message=F,eval=TRUE,cache = T-----------------------
kmComparison$comparison_1to2$ccfResults %>%
  topResultsPerCell() %>%
  ungroup() %>%
  mutate(cellIndex = cmcR:::linear_to_matrix(index = (cellNum %% ceiling(sqrt(max(cellNum)))) +
                                               floor((ceiling(sqrt(max(cellNum)))^2 - cellNum)/ceiling(sqrt(max(cellNum))))*ceiling(sqrt(max(cellNum))) +
                                               ifelse(cellNum %% ceiling(sqrt(max(cellNum))) == 0,ceiling(sqrt(max(cellNum))),0),
                                             nrow = ceiling(sqrt(max(cellNum))),
                                             byrow = TRUE)) %>%
  mutate(`Cell index` = cellIndex,
         `Pairwise-complete corr.` = round(ccf,3),
         `FFT-based corr.` = round(fft.ccf,3),
         `Non-missing percent.` = round(100*nonMissingProportion,2)) %>%
  select(-c(cellNum,cellID)) %>%
  select(c(`Cell index`,`Pairwise-complete corr.`,`FFT-based corr.`,dx,dy,theta,`Non-missing percent.`)) %>%
  filter(theta == -24) %>%
  arrange(`Cell index`) %>%
  head() %>%
  knitr::kable(caption = "\\label{tab:cellCCF} Example of output from correlation cell comparison procedure between Fadul 1-1 and Fadul 1-2 rotated by -24 degrees. Due to the large proportion of missing values that are replaced to compute the FFT-based correlation, the pairwise-complete correlation is most often greater than the FFT-based correlation.",
               format = "latex",
               align = c("|c","c","c","r","r","r","c|"),
               col.names = c("Cell Index",
                             "Pairwise-comp. corr.",
                             "FFT-based corr.",
                             "$\\Delta$x",
                             "$\\Delta$y",
                             "$\\theta$",
                             "Non-NA percent."),
               escape = FALSE) %>%
  kableExtra::kable_styling(latex_options = "hold_position",
                            position = "center")


## ----echo=FALSE,cache = T-----------------------------------------------------
# fadul2.1 <- cmcR::selectBFImpression_sample_x3p("data/fadul2-1.x3p",
#                                                 gaussFilterWavelength = c(16,500))

fadul2.1 <- x3ptools::read_x3p("data/fadul2-1.x3p") %>%
  preProcess_cropBFExterior(radiusOffset = -30,
                 agg_function = median) %>%
  preProcess_filterBFInterior(radiusOffset = 200) %>%
  preProcess_removeBFTrend(statistic = "quantile",
                              tau = .5,
                              method = "fn") %>%
  cmcR::preProcess_gaussFilter() %>%
  x3ptools::sample_x3p()

knmComparison <- cmcR::cellCCF_bothDirections(x3p1 = fadul1.1,
                                              x3p2 = fadul2.1,
                                              #x3p1 = fadul1.2$x3p,
                                              #x3p2 = fadul2.1$x3p,
                                              cellNumHoriz = 8,
                                              minObservedProp = .1,
                                              regionToCellProp = 9,
                                              ccfMethod = "fftThenPairwise")

knmCMC <- cmcR::cmcFilter_improved(knmComparison,
                                   ccf_thresh = .5,
                                   dx_thresh = 20,
                                   theta_thresh = 6,
                                   missingTheta_decision = "fail")


## ----cache = T,include=FALSE--------------------------------------------------
knmCMCPlot <- cmcR::cmcPlot(fadul1.1,
              fadul2.1,
              cellCCF_bothDirections_output = knmComparison,
              cmcFilter_improved_output = knmCMC,
              type = "faceted",
              x3pNames = c("Fadul 1-1","Fadul 2-1"),
              legend.quantiles = c(0,.01,.2,.5,.8,.99,1),
              height.colors = colorspace::desaturate(c('#7f3b08','#b35806',
                                                       '#e08214','#fdb863',
                                                       '#fee0b6','#f7f7f7',
                                                       '#d8daeb','#b2abd2',
                                                       '#8073ac','#542788',
                                                       '#2d004b'),
                                                     amount = .75),
              cell.colors = c("#a60b00","#1b03a3"),
              cell.alpha = .15,
              na.value = "gray80")

#We want Fadul 1-1 to be the reference cartridge case for all 3 of these plots
#to enforce visual consistency. The original CMC method applied to Fadul 1-1 and
#Fadul 1-2 results in two CMC counts (one for each comparison direction), the
#minimum of which is associated with Fadul 1-2 being the reference cartridge
#case. Unless changed, the cmcPlot function (by design) will then create a Top
#Vote CMC plot with Fadul 1-2 as the reference cartridge case (which we don't
#want). Since we desire Fadul 1-1 to be the reference for this particular
#illustrative example, we will artificially inflate the number CMCs in the Fadul
#1-2 vs. Fadul 1-1 direction (i.e., with Fadul 1-2 as the reference) so as to
#force Fadul 1-1 to be plotted as the reference.

kmCMC_fake <- kmCMC

kmCMC_fake$originalMethodCMCs$comparison_2to1 <- bind_rows(kmCMC_fake$originalMethodCMCs$comparison_2to1,
                                                    kmCMC_fake$originalMethodCMCs$comparison_2to1)

kmCMCPlot <- cmcR::cmcPlot(fadul1.1,
                          fadul1.2,
                          # fadul1.2$x3p,
                          # fadul1.1$x3p,
                          cellCCF_bothDirections_output = kmComparison,
                          cmcFilter_improved_output = kmCMC_fake,
                          # cellCCF_bothDirections_output = kmComparison_swapped,
                          # cmcFilter_improved_output = kmCMC_swapped,
                          x3pNames = c("Fadul 1-1","Fadul 1-2"),
                          legend.quantiles = c(0,.01,.2,.5,.8,.99,1),
                          height.colors = colorspace::desaturate(c('#7f3b08',
                                                                   '#b35806',
                                                                   '#e08214',
                                                                   '#fdb863',
                                                                   '#fee0b6',
                                                                   '#f7f7f7',
                                                                   '#d8daeb',
                                                                   '#b2abd2',
                                                                   '#8073ac',
                                                                   '#542788',
                                                                   '#2d004b'),
                                                                 amount = .75),
        cell.colors = c("#a60b00","#1b03a3"),
        cell.alpha = .15,
        na.value = "gray80")

knmCMCPlot_list_comparison1to2 <- cmcR::cmcPlot(fadul1.1,
                                 fadul2.1,
                                 cellCCF_bothDirections_output = knmComparison,
                                 cmcFilter_improved_output = knmCMC,
                                 type = "list",
                                 x3pNames = c("Fadul 1-1","Fadul 2-1"),
                                 legend.quantiles = c(0,.01,.2,.5,.8,.99,1),
                                 height.colors = colorspace::desaturate(c('#7f3b08','#b35806',
                                                                          '#e08214','#fdb863',
                                                                          '#fee0b6','#f7f7f7',
                                                                          '#d8daeb','#b2abd2',
                                                                          '#8073ac','#542788',
                                                                          '#2d004b'),
                                                                        amount = .75),
                                 cell.colors = c("black","black"),
                                 cell.alpha = .15,
                                 cell.text.size = 7,
                                 na.value = "gray90")

kmCMCPlot_list_comparison1to2 <- cmcR::cmcPlot(fadul1.1,
                                fadul1.2,
                                # fadul1.2$x3p,
                                # fadul1.1$x3p,
                                cellCCF_bothDirections_output = kmComparison,
                                cmcFilter_improved_output = kmCMC_fake,
                                # cellCCF_bothDirections_output = kmComparison_swapped,
                                # cmcFilter_improved_output = kmCMC_swapped,
                                x3pNames = c("Fadul 1-1","Fadul 1-2"),
                                legend.quantiles = c(0,.01,.2,.5,.8,.99,1),
                                height.colors = colorspace::desaturate(c('#7f3b08',
                                                                         '#b35806',
                                                                         '#e08214',
                                                                         '#fdb863',
                                                                         '#fee0b6',
                                                                         '#f7f7f7',
                                                                         '#d8daeb',
                                                                         '#b2abd2',
                                                                         '#8073ac',
                                                                         '#542788',
                                                                         '#2d004b'),
                                                                       amount = .75),
                                cell.colors = c("black","black"),
                                cell.alpha = .15,
                                cell.text.size = 7,
                                na.value = "gray90",
                                type = "list")

# kmCMCPlot_list_comparison1to2$originalMethodCMCs$`Fadul 1-1`$layers[[2]]$mapping$fill <- NULL

# kmCMCPlot_list_comparison1to2$originalMethodCMCs$`Fadul 1-1`$layers[[2]]$mapping$colour <- "black"
# 
# kmCMCPlot_list_comparison1to2$originalMethodCMCs$`Fadul 1-1`$layers[[2]]$mapping$size <-


## ----echo=TRUE,eval=FALSE-----------------------------------------------------
#> kmCMC <- cmcR::cmcFilter_improved(kmComparison,
#>                                   ccf_thresh = .5,
#>                                   dx_thresh = 20,
#>                                   dy_thresh = 20,
#>                                   theta_thresh = 6)


## ----echo=FALSE,warning=FALSE,message=FALSE,cache = F,fig.cap='\\label{fig:topVoteCMCPlot} CMC results for the comparison between Fadul 1-1 and Fadul 1-2 using the original method of \\citet{song_proposed_2013}. Fadul 1-1 (left) acts as the reference cartridge case to which Fadul 1-2 (right) is compared. The 19 blue cells on the right show the phase registration of each cell identified as congruent. The 12 red cells on the right show where cells \\emph{not} identified as congruent achieve the maximum pairwise-complete correlation across all rotations of Fadul 1-2 considered.',fig.align='center',fig.pos='htbp',out.width='\\textwidth'----

kmCMCPlot$originalMethodCMCs


## ----echo=FALSE,warning=FALSE,message=FALSE,cache = F,fig.cap='\\label{fig:highCMCPlot} CMC results for the comparison between Fadul 1-1 and Fadul 1-2 using the High CMC method, which tends to assign a higher CMC count to cartridge case pairs that pass its more stringent criteria (discussed in section REF). This matching cartridge case pair passes the High CMC criterion, meaning 26 cells are identified as congruent while only 6 are not.',fig.align='center',fig.pos='htbp',out.width='\\textwidth'----

kmCMCPlot$highCMCs


## ----echo=FALSE,warning=FALSE,message=FALSE,cache = F,fig.cap='\\label{fig:knmCMCPlot} CMC results for the comparison between the non-match pair Fadul 1-1 and Fadul 2-1. This pair failed the High CMC criterion and also yielded 0 CMCs under the original method of \\citet{song_proposed_2013}. As such, all 32 cells considered are classified as non-congruent.', fig.align='center',fig.pos='htbp',out.width='\\textwidth'----

knmCMCPlot$originalMethodCMCs


## ----echo=FALSE,fig.cap='\\label{fig:topVoteCMC_sensitivity} CMC count relative frequencies under the original method of \\citet{song_proposed_2013} for various $T_{\\Delta x}, T_{\\Delta y}$, and $T_{\\text{CCF}}$ values and $T_{\\theta} = 6$ degrees. Thresholds become more stringent as the $T_{\\Delta x}, T_{\\Delta y}$ values decrease and $T_{\\text{CCF}}$ values increase. AUCs based on the minimum CMC count classification threshold and Between/Within group variance ratios are shown for each combination. AUC $= 1$ corresponds to perfect separation of the match and non-match CMC count distributions and a larger Between/Within group variance ratio value is preferred.', fig.align='left',fig.pos='htbp',out.width='\\textwidth'----

load("data/cmcCountData.RData")

cmcs_varRatios <-  cmcs_downSampled %>%
  pivot_longer(cols = c(originalMethodCMCs,highCMCs),names_to = "cmcType",values_to = "cmcCount") %>%
  group_by(theta_thresh,ccf_thresh,dx_thresh,missingTheta_decision,compareInitialAndHighThetas,cmcType,type) %>%
  mutate(group_var = var(cmcCount)) %>%
  ungroup(type) %>%
  mutate(withinGroup_var = mean(group_var),
         betweenGroup_var = var(cmcCount)) %>%
  mutate(varRatio = betweenGroup_var/withinGroup_var)

cmcs_aucs <- cmcs_downSampled %>%
  group_by(theta_thresh,ccf_thresh,dx_thresh,missingTheta_decision,compareInitialAndHighThetas) %>%
  dplyr::group_split()  %>%
  map_dfr(function(paramConditionedResults = .x){
    
    topVoteROC <- paramConditionedResults %>%
      pROC::roc(response = type,
                predictor = originalMethodCMCs,
                levels = c("non-match","match"),
                quiet = TRUE)
    
    highROC <- paramConditionedResults %>%
      pROC::roc(response = type,
                predictor = highCMCs,
                levels = c("non-match","match"),
                quiet = TRUE)
    
    paramConditionedResults %>%
      mutate(topVoteCMC_AUC = rep(topVoteROC$auc),
             highCMC_AUC = rep(highROC$auc))
  }) 

topVoteCMC_labels <- cmcs_aucs %>%
  ungroup() %>%
  filter(missingTheta_decision == "replace" & compareInitialAndHighThetas == TRUE & theta_thresh == 6) %>%
  filter((ccf_thresh %in% c(.4,.5) | (ccf_thresh > .39 & ccf_thresh < .41)) & dx_thresh %in% c(10,20,30)) %>%
  inner_join(cmcs_varRatios %>%
               ungroup() %>%
               filter(cmcType == "originalMethodCMCs") %>%
               select(ccf_thresh,dx_thresh,theta_thresh,varRatio) %>%
               distinct(),by = c("theta_thresh","ccf_thresh","dx_thresh")) %>%
  select(ccf_thresh,dx_thresh,topVoteCMC_AUC,varRatio) %>%
  distinct() %>%
  rename(`Trans. Thresh` = dx_thresh,
         `CCF Thresh` = ccf_thresh) %>%
  mutate(label = sprintf("AUC: %.4f\nVar. Ratio: %.2f", topVoteCMC_AUC, varRatio))
  

cmcs_aucs %>%
  ungroup() %>%
  filter(missingTheta_decision == "replace" & compareInitialAndHighThetas == TRUE & theta_thresh == 6) %>% 
  filter((ccf_thresh %in% c(.4,.5) | (ccf_thresh > .39 & ccf_thresh < .41)) & dx_thresh %in% c(10,20,30)) %>%
  select(-c(missingTheta_decision,compareInitialAndHighThetas)) %>%
  distinct() %>% 
  rename(`Trans. Thresh` = dx_thresh,
         `CCF Thresh` = ccf_thresh) %>%
  ggplot() +
  geom_bar(aes(x = originalMethodCMCs,
               y = ..prop..,
               fill = type),
           stat = "count",
           alpha = .7) +
  geom_label(data = topVoteCMC_labels,
             aes(x = 15,
                 y = .7,
                 label = label),
             size = 4) +
  facet_grid(rows = vars(`CCF Thresh`),
             cols = vars(`Trans. Thresh`),
             scales = "free_y",labeller = label_both) +
  scale_fill_manual(values = c("#40B0A6","#E1BE6A")) +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() + 
  theme(legend.position = "bottom",
        strip.text = element_text(size = 7)) +
  xlab("Original Method CMC Count") +
  ylab("Relative Frequency")


## ----echo=FALSE, fig.cap='\\label{fig:highCMC_sensitivity} High CMC count relative frequencies for various $T_{\\Delta x}, T_{\\Delta y}$, and $T_{\\text{CCF}}$ values and $T_{\\theta} = 6$ degrees, \\code{missingTheta\\_decision = "fail"}, and \\code{compareInitialAndHighThetas = TRUE}. Thresholds become more stringent as the $T_{\\Delta x}, T_{\\Delta y}$ values decrease and $T_{\\text{CCF}}$ values increase. AUCs based on the minimum CMC count classification threshold and Between/Within group variance ratios are shown for each combination.',fig.align='left',fig.pos='htbp',out.width='\\textwidth'----
highCMC_conservative_labels <- cmcs_aucs %>%
  ungroup() %>%
  filter(missingTheta_decision == "fail" & compareInitialAndHighThetas == TRUE & theta_thresh == 6) %>%
  filter((ccf_thresh %in% c(.4,.5) | (ccf_thresh > .39 & ccf_thresh < .41)) & dx_thresh %in% c(10,20,30)) %>%
  inner_join(cmcs_varRatios %>%
               ungroup() %>%
               filter(cmcType == "highCMCs" & missingTheta_decision == "fail" & compareInitialAndHighThetas == TRUE) %>%
               select(ccf_thresh,dx_thresh,theta_thresh,varRatio) %>%
               distinct(),
             by = c("theta_thresh","ccf_thresh","dx_thresh")) %>%
  select(ccf_thresh,dx_thresh,highCMC_AUC,varRatio) %>%
  distinct() %>%
  rename(`Trans. Thresh` = dx_thresh,
         `CCF Thresh` = ccf_thresh) %>%
  mutate(label = sprintf("AUC: %.4f\nVar. Ratio: %.2f", highCMC_AUC, varRatio))

cmcs_aucs %>%
  ungroup() %>%
  filter(missingTheta_decision == "fail" & compareInitialAndHighThetas == TRUE & theta_thresh == 6) %>% 
  filter((ccf_thresh %in% c(.4,.5) | (ccf_thresh > .39 & ccf_thresh < .41)) & dx_thresh %in% c(10,20,30)) %>%
  distinct() %>%
  rename(`Trans. Thresh` = dx_thresh,
         `CCF Thresh` = ccf_thresh) %>%
  ggplot() +
  geom_bar(aes(x = highCMCs,
               y = ..prop..,
               fill = type),
           stat = "count",
           alpha = .7) +
  geom_label(data = highCMC_conservative_labels,
             aes(x = 22.5,
                 y = .7,
                 label = label),
             size = 4) +
  facet_grid(rows = vars(`CCF Thresh`),
             cols = vars(`Trans. Thresh`),
             scales = "free_y",labeller = label_both) +
  scale_fill_manual(values = c("#40B0A6","#E1BE6A")) +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() + 
  theme(legend.position = "bottom",
        strip.text = element_text(size = 7)) +
  xlab("High CMC Count") +
  ylab("Relative Frequency")


## ----echo=FALSE, fig.cap='\\label{fig:highCMC_comparison_sensitivity} High CMC count relative frequencies for various interpretations on the High CMC method description given in \\citet{tong_improved_2015} with $T_{\\Delta x} = T_{\\Delta y} = 20$, $T_{\\text{CCF}} = .5$, and $T_{\\theta} = 6$ degrees. In particlar, compares CMC results for different the 4 combinations of \\code{missingTheta\\_decision = "fail"} or \\code{"replace"} and \\code{compareInitialAndHighThetas = TRUE} or \\code{FALSE} (see \\@sec:highMethod for details). Different interpretations of the High CMC method description lead to different results, particularly with respect to the behavior of the non-match CMC count distribution. Additionally, as evidenced by the stillness of the matching CMC count distribution, this illustrates how the High CMC criterion functions essentially as a preliminary classification rule for non-matches vs. matches.',fig.align='left',fig.pos='htbp',out.width='\\textwidth'----

highCMC_comparison_labels <- cmcs_aucs %>%
  ungroup() %>%
  filter(missingTheta_decision %in% c("replace","fail")  & theta_thresh == 6) %>%
  filter(ccf_thresh == .5 & dx_thresh == 20) %>%
  inner_join(cmcs_varRatios %>%
               ungroup() %>%
               filter(missingTheta_decision %in% c("replace","fail")  & theta_thresh == 6 & ccf_thresh == .5 & dx_thresh == 20 & cmcType == "highCMCs") %>%
               select(missingTheta_decision,compareInitialAndHighThetas,ccf_thresh,theta_thresh,dx_thresh,varRatio) %>%
               distinct(),
             by = c("missingTheta_decision","compareInitialAndHighThetas","ccf_thresh","theta_thresh","dx_thresh")) %>%  
  select(missingTheta_decision,compareInitialAndHighThetas,highCMC_AUC,varRatio)  %>%
  distinct() %>%
  rename(`Missing Theta` = missingTheta_decision,
         `Compare Initial/High Thetas` = compareInitialAndHighThetas) %>%
  mutate(label = sprintf("AUC: %.4f\nVar. Ratio: %.2f", highCMC_AUC, varRatio))

cmcs_aucs %>%
  ungroup() %>%
  filter(missingTheta_decision %in% c("replace","fail") & theta_thresh == 6) %>%
  filter(ccf_thresh == .5 & dx_thresh == 20)  %>%
  distinct() %>%
  rename(`Missing Theta` = missingTheta_decision,
         `Compare Initial/High Thetas` = compareInitialAndHighThetas) %>%
  ggplot() +
  geom_bar(aes(x = highCMCs,
               y = ..prop..,
               fill = type),
           stat = "count",
           alpha = .7) +
  geom_label(data = highCMC_comparison_labels,
             aes(x = 22.5,
                 y = .5,
                 label = label),
             size = 4) +
  facet_grid(rows = vars(`Missing Theta`),
             cols = vars(`Compare Initial/High Thetas`),
             scales = "free_y",labeller = label_both) +
  scale_fill_manual(values = c("#40B0A6","#E1BE6A")) +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() + 
  theme(legend.position = "bottom",
        strip.text = element_text(size = 7)) +
  xlab("High CMC Count") +
  ylab("Relative Frequency")


## ----echo=FALSE, fig.cap='\\label{fig:highCMC_noLeveling_sensitivity} High CMC count relative frequencies for various $T_{\\Delta x}, T_{\\Delta y}$, and $T_{\\text{CCF}}$ values and $T_{\\theta} = 6$ degrees, \\code{missingTheta\\_decision = "fail"}, and \\code{compareInitialAndHighThetas = TRUE} where the breech face impression surfaces were not leveled during preprocessing. Compared to \\autoref{fig:highCMC_sensitivity}, this demonstrates that High CMC method becomes less effective if this processing step is skipped. Modularization of the processing procedures into a "pipeline" is useful for uncovering such sensitivites.',fig.align='left',fig.pos='htbp',out.width='\\textwidth'----
cmcs_noLeveling_varRatios <-  cmcs_downSampled_noLeveling %>%
  pivot_longer(cols = c(originalMethodCMCs,highCMCs),names_to = "cmcType",values_to = "cmcCount") %>%
  group_by(theta_thresh,ccf_thresh,dx_thresh,missingTheta_decision,compareInitialAndHighThetas,cmcType,type) %>%
  mutate(group_var = var(cmcCount)) %>%
  ungroup(type) %>%
  mutate(withinGroup_var = mean(group_var),
         betweenGroup_var = var(cmcCount)) %>%
  mutate(varRatio = betweenGroup_var/withinGroup_var)

cmcs_noLeveling_aucs <- cmcs_downSampled_noLeveling %>%
  group_by(theta_thresh,ccf_thresh,dx_thresh,missingTheta_decision,compareInitialAndHighThetas) %>%
  dplyr::group_split()  %>%
  map_dfr(function(paramConditionedResults = .x){
    
    topVoteROC <- paramConditionedResults %>%
      pROC::roc(response = type,
                predictor = originalMethodCMCs,
                levels = c("non-match","match"),
                quiet = TRUE)
    
    highROC <- paramConditionedResults %>%
      pROC::roc(response = type,
                predictor = highCMCs,
                levels = c("non-match","match"),
                quiet = TRUE)
    
    paramConditionedResults %>%
      mutate(topVoteCMC_AUC = rep(topVoteROC$auc),
             highCMC_AUC = rep(highROC$auc))
  })

highCMC_noLeveling_labels <- cmcs_noLeveling_aucs %>%
  ungroup() %>%
  filter(missingTheta_decision == "fail" & compareInitialAndHighThetas == TRUE & theta_thresh == 6) %>%
  filter((ccf_thresh %in% c(.4,.5) | (ccf_thresh > .39 & ccf_thresh < .41)) & dx_thresh %in% c(10,20,30)) %>%
  inner_join(cmcs_noLeveling_varRatios %>%
               ungroup() %>%
               filter(cmcType == "highCMCs" & missingTheta_decision == "fail" & compareInitialAndHighThetas == TRUE) %>%
               select(ccf_thresh,dx_thresh,theta_thresh,varRatio) %>%
               distinct(),
             by = c("theta_thresh","ccf_thresh","dx_thresh")) %>%
  select(ccf_thresh,dx_thresh,highCMC_AUC,varRatio) %>%
  distinct() %>%
  rename(`Trans. Thresh` = dx_thresh,
         `CCF Thresh` = ccf_thresh) %>%
  mutate(label = sprintf("AUC: %.4f\nVar. Ratio: %.2f", highCMC_AUC, varRatio))

cmcs_noLeveling_aucs %>%
  ungroup() %>%
  filter(missingTheta_decision == "fail" & compareInitialAndHighThetas == TRUE & theta_thresh == 6) %>% 
  filter((ccf_thresh %in% c(.4,.5) | (ccf_thresh > .39 & ccf_thresh < .41)) & dx_thresh %in% c(10,20,30)) %>%
  distinct() %>%
  rename(`Trans. Thresh` = dx_thresh,
         `CCF Thresh` = ccf_thresh) %>%
  ggplot() +
  geom_bar(aes(x = highCMCs,
               y = ..prop..,
               fill = type),
           stat = "count",
           alpha = .7) +
  geom_label(data = highCMC_noLeveling_labels,
             aes(x = 22.5,
                 y = .7,
                 label = label),
             size = 4) +
  facet_grid(rows = vars(`CCF Thresh`),
             cols = vars(`Trans. Thresh`),
             scales = "free_y",labeller = label_both) +
  scale_fill_manual(values = c("#40B0A6","#E1BE6A")) +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() + 
  theme(legend.position = "bottom",
        strip.text = element_text(size = 7)) +
  xlab("High CMC Count") +
  ylab("Relative Frequency")

