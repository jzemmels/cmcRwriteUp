## ----localDataDir, include=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------
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


## ---- derivativeImagesDir,include=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------
if(!dir.exists("derivatives")){
  dir.create("derivatives")
}


## ----setup,echo=FALSE,message=FALSE,warning=FALSE----------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(cache = T, dpi = 300, fig.width = 8, fig.height = 4, out.width = "\\textwidth", dpi = 300)
library(cmcR) # remotes::install_github("CSAFE-ISU/cmcR")
library(tidyverse)
library(x3ptools) # remotes::install_github("heike/x3ptools")
library(rgl)
library(ggpcp)


## ----eval=FALSE,echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------
## # SVP comment: Should do this in a tidy way with less code if possible...
## library(cmcR)
## 
## fadul1.1_id <- "DownloadMeasurement/2d9cc51f-6f66-40a0-973a-a9292dbee36d"
## # Same source comparison
## fadul1.2_id <- "DownloadMeasurement/cb296c98-39f5-46eb-abff-320a2f5568e8"
## # Different source comparison
## fadul2.1_id <- "DownloadMeasurement/8ae0b86d-210a-41fd-ad75-8212f9522f96"
## 
## #Code to download breech face impressions:
## 
## # Aside: while the URL says "NRBTD", it's
## #actually the NIST Ballistics Toolmark Research Database (so their URL
## #is mistaken)
## 
## nbtrd_url <- "https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement"
## download.file(
##   file.path(nbtrd_url , fadul1.1_id), destfile = "data/fadul1-1.x3p", mode = "wb")
## download.file(
##   file.path(nbtrd_url , fadul1.2_id), destfile = "data/fadul1-2.x3p", mode = "wb")
## download.file(
##   file.path(nbtrd_url, fadul2.1_id), destfile = "data/fadul2-1.x3p", mode = "wb")


## ----eval=FALSE,echo=TRUE----------------------------------------------------------------------------------------------------------------------------------------------------------------------
## library(cmcR)
## 
## nbtrd_url <- "https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement"
## 
## x3p_ids <- c("DownloadMeasurement/2d9cc51f-6f66-40a0-973a-a9292dbee36d",
##              "DownloadMeasurement/cb296c98-39f5-46eb-abff-320a2f5568e8",
##              "DownloadMeasurement/8ae0b86d-210a-41fd-ad75-8212f9522f96")
## 
## file_names <- c("fadul1-1.x3p","fadul1-2.x3p","fadul2-1.x3p")
## 
## purrr::walk2(.x = x3p_ids,
##              .y = file_names,
##              .f = function(x3p_id,file_name){
##                download.file(url = file.path(nbtrd_url, x3p_id),
##                              destfile = paste0("data/",file_name),mode = "wb")
##              })


## ----echo=FALSE,fig.cap='\\label{fig:ccPair_combined} A cartridge case pair with visible breech face impressions under a microscrope.  A thin line can be seen separating the two views. The degree to which the markings coincide is used to conclude whether the pair comes from the same source.',fig.pos='htbp',out.width="\\textwidth"----
knitr::include_graphics("images/cartridgeCasePair_comparison_with_line.PNG")


## ---- fadul1-1Screenshot,include=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------
fadul1.1 <- x3ptools::x3p_read("data/fadul1-1.x3p")

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


## ----fadul1-2Screenshot,include=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------
fadul1.2 <- x3ptools::x3p_read("data/fadul1-2.x3p")

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


## ---- rawBFs,echo=FALSE,fig.cap='\\label{fig:cartridgeCasePair} Unprocessed surface matrices of the known-match Fadul 1-1 (left) and Fadul 1-2 (right) \\citep{fadul_empirical_2011}. The observations in the corners of these surface matrices are artifacts of the staging area in which these scans were taken. The holes on the interior of the primer surfaces are caused by the firing pin striking the primer during the firing process. The region of the primer around this hole does not come into uniform contact with the breech face of the firearm.', fig.subcap=c('',''),fig.align='center',fig.pos='htbp',out.width='.49\\linewidth',out.height='.49\\linewidth'----
knitr::include_graphics(c("derivatives/fadul1-1.png","derivatives/fadul1-2.png"))


## ----load-data, include = F, cache = T---------------------------------------------------------------------------------------------------------------------------------------------------------

fadul1.1 <- x3ptools::x3p_read("data/fadul1-1.x3p") %>%
  cmcR::preProcess_crop(region = "exterior",
                        radiusOffset = -30) %>%
  cmcR::preProcess_crop(region = "interior",
                        radiusOffset = 200) %>%
  cmcR::preProcess_removeTrend(statistic = "quantile",
                                 tau = .5,
                                 method = "fn") %>%
  cmcR::preProcess_gaussFilter() %>%
  x3ptools::sample_x3p()

fadul1.2 <- x3ptools::x3p_read("data/fadul1-2.x3p") %>%
  cmcR::preProcess_crop(region = "exterior",
                        radiusOffset = -30) %>%
  cmcR::preProcess_crop(region = "interior",
                        radiusOffset = 200) %>%
  cmcR::preProcess_removeTrend(statistic = "quantile",
                                 tau = .5,
                                 method = "fn") %>%
  cmcR::preProcess_gaussFilter() %>%
  x3ptools::sample_x3p()


## ----cmc-ccf, include = F, cache = T-----------------------------------------------------------------------------------------------------------------------------------------------------------
kmComparisonFeatures <- purrr::map_dfr(seq(-30,30,by = 3),
                                       ~ comparison_allTogether(reference = fadul1.1,
                                                                target = fadul1.2,
                                                                numCells = 64,
                                                                maxMissingProp = .85,
                                                                theta = .))

kmComparisonFeatures_rev <- purrr::map_dfr(seq(-30,30,by = 3),
                                           ~ comparison_allTogether(reference = fadul1.2,
                                                                    target = fadul1.1,
                                                                    numCells = 64,
                                                                    maxMissingProp = .85,
                                                                    theta = .))

kmComparison_cmcs <- kmComparisonFeatures %>%
  mutate(originalMethodClassif = decision_CMC(cellIndex = cellIndex,
                                              x = x,
                                              y = y,
                                              theta = theta,
                                              corr = pairwiseCompCor,
                                              xThresh = 20,
                                              thetaThresh = 6,
                                              corrThresh = .5),
         highCMCClassif = decision_CMC(cellIndex = cellIndex,
                                              x = x,
                                              y = y,
                                              theta = theta,
                                              corr = pairwiseCompCor,
                                              xThresh = 20,
                                              thetaThresh = 6,
                                              corrThresh = .5,
                                              tau = 1))

kmComparison_cmcs_rev <- kmComparisonFeatures_rev %>%
  mutate(originalMethodClassif = decision_CMC(cellIndex = cellIndex,
                                              x = x,
                                              y = y,
                                              theta = theta,
                                              corr = pairwiseCompCor,
                                              xThresh = 20,
                                              thetaThresh = 6,
                                              corrThresh = .5),
         highCMCClassif = decision_CMC(cellIndex = cellIndex,
                                              x = x,
                                              y = y,
                                              theta = theta,
                                              corr = pairwiseCompCor,
                                              xThresh = 20,
                                              thetaThresh = 6,
                                              corrThresh = .5,
                                              tau = 1))

bind_rows(kmComparison_cmcs,
          kmComparison_cmcs_rev) %>%
  filter(highCMCClassif == "CMC") %>%
  group_by(cellIndex) %>%
  filter(pairwiseCompCor == max(pairwiseCompCor))


## ----cache=FALSE, include=F--------------------------------------------------------------------------------------------------------------------------------------------------------------------
fadul1.1_original <- x3ptools::x3p_read("data/fadul1-1.x3p")

fadul1.1_croppedExt <- cmcR::preProcess_crop(fadul1.1_original,
                                             region = "exterior",
                                             radiusOffset = -30)

fadul1.1_croppedInt <- cmcR::preProcess_crop(fadul1.1_croppedExt,
                                             region = "interior",
                                             radiusOffset = 200)

fadul1.1_medRemoved <-   cmcR::preProcess_removeTrend(fadul1.1_croppedInt,
                                                        statistic = "quantile",
                                                        tau = .5,
                                                        method = "fn")



fadul1.1_downsampled <- x3ptools::sample_x3p(fadul1.1_medRemoved,
                                             m = 2)

fadul1.1_bpFiltered<- cmcR::preProcess_gaussFilter(x3p = fadul1.1_downsampled,
                                                   wavelength = c(16,500),
                                                   filtertype = "bp")


## ---- echo = F,warning = F,message = F,cache = T,fig.cap='\\label{fig:processingPipeline} Illustration of the  preprocessing pipeline implemented in \\CRANpkg{cmcR}.  At each stage, the variability in height across the scan decreases as extraneous sources of noise are removed.',fig.align='center',fig.pos='htbp',out.width='\\textwidth', message = F, warning = F----

preProcessingPlot <- cmcR::x3pListPlot(list(fadul1.1_original,
                                            fadul1.1_croppedInt,
                                            fadul1.1_medRemoved,
                                            fadul1.1_bpFiltered) %>%
                                         set_names(c("(1) Original \n x3p_read()",
                                                    "(2) Crop exterior/interior \n preProcess_crop()",
                                                    "(3) Level surface \n preProcess_removeTrend()",
                                                    "(4) Band-pass filter \n preProcess_gaussFilter()")),
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


## ---- echo=TRUE,eval=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------
## # Step (1)
## fadul1.1 <- x3ptools::x3p_read("data/fadul1-1.x3p")


## ---- echo=TRUE,eval=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------
## # Step (2)
## fadul1.1_cropped <- fadul1.1%>%
##   cmcR::preProcess_crop(region = "exterior") %>%
##   cmcR::preProcess_crop(region = "interior")


## ---- echo=TRUE,eval=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------
## # Step (3)
## fadul1.1_deTrended <- fadul1.1_cropped %>%
##   preProcess_removeTrend(statistic = "quantile",
##                          tau = .5,
##                          method = "fn")


## ---- echo=TRUE,eval=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------
## # Step (4)
## fadul1.1_processed <- fadul1.1_deTrended %>%
##   preProcess_gaussFilter(filtertype = "bp",
##                          wavelength = c(16,500)) %>%
##   x3ptools::x3p_sample(m = 2)


## ----echo=FALSE,cache = T,fig.cap='\\label{fig:processedScans} Fadul 1-1 and Fadul 1-2 after preprocessing. Similar striated markings are now easier to visually identify on both surfaces. It is now clearer that one of the scans needs to be rotated to align better with the other.',fig.align='center',fig.pos='htbp',out.width='\\textwidth', message = F, warning = F----

cmcR::x3pListPlot(x3pList = list("Fadul 1-1" = fadul1.1,
                                 "Fadul 1-2" = fadul1.2),
                  # x3pList = list("Fadul 1-1" = fadul1.1$x3p,
                  #"Fadul 1-2" = fadul1.2$x3p),
                  type = "faceted",
                  rotate = 90,
                  legend.quantiles = c(0,.01,.2,.5,.8,.99,1)) +
  guides(fill = guide_colourbar(barheight = grid::unit(2.6,"inches"),
                                label.theme = element_text(size = 7),
                                title.theme = ggplot2::element_text(size = 9),
                                frame.colour = "black",
                                ticks.colour = "black")) +
  theme(legend.position = c(1.11,.551),plot.margin = ggplot2::margin(c(0,3,.2,0),unit = "cm"))


## ---- echo=FALSE,fig.cap='\\label{fig:cmc_illustration} Illustration of comparing a cell in the reference cartridge case scan (left) to a larger region in a questioned cartridge case scan (right). Every one of the cells in the reference cartridge case is similarly paired with a region in the questioned cartridge case.  To determine the rotation at which the two cartridge cases align, the cell-region pairs are compared for various rotations of the questioned cartridge case.',fig.align='center',fig.pos='htbp',out.width='.75\\textwidth'----

knitr::include_graphics("images/cmc_illustration.PNG")


## ----echo=TRUE,eval=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------
## kmComparisonFeatures <- purrr::map_dfr(seq(-30,30,by = 3),
##                                        ~ comparison_allTogether(reference = fadul1.1,
##                                                                 target = fadul1.2,
##                                                                 numCells = 64,
##                                                                 maxMissingProp = .85,
##                                                                 theta = .))
## 
## kmComparisonFeatures_rev <- purrr::map_dfr(seq(-30,30,by = 3),
##                                            ~ comparison_allTogether(reference = fadul1.2,
##                                                                     target = fadul1.1,
##                                                                     numCells = 64,
##                                                                     maxMissingProp = .85,
##                                                                     theta = .))


## ----echo=FALSE,warning=F,message=F,eval=TRUE,cache = T----------------------------------------------------------------------------------------------------------------------------------------
kmComparisonFeatures %>%
  mutate(`Cell index` = cellIndex,
         `Pairwise-complete corr.` = round(pairwiseCompCor,3),
         `FFT-based corr.` = round(fft_ccf,3)) %>%
  select(c(`Cell index`,`Pairwise-complete corr.`,`FFT-based corr.`,x,y,theta)) %>%
  filter(theta == -24) %>%
  arrange(`Cell index`) %>%
  head(5) %>%
  knitr::kable(caption = "\\label{tab:cellCCF} Example of output from correlation cell comparison procedure between Fadul 1-1 and Fadul 1-2 rotated by -24 degrees. Due to the large proportion of missing values that are replaced to compute the FFT-based correlation, the pairwise-complete correlation is most often greater than the FFT-based correlation.",
               format = "latex",
               align = c("|c","c","c","r","r","r|"),
               col.names = c("Cell Index",
                             "Pairwise-comp. corr.",
                             "FFT-based corr.",
                             "$\\Delta$x",
                             "$\\Delta$y",
                             "$\\theta$"),
               escape = FALSE, 
               booktabs = TRUE) %>%
  kableExtra::kable_styling(latex_options = "striped", 
                            position = "center", 
                            stripe_color="lightgray")


## ----echo=FALSE,cache = T----------------------------------------------------------------------------------------------------------------------------------------------------------------------
fadul2.1 <- x3ptools::x3p_read("data/fadul2-1.x3p") %>%
  cmcR::preProcess_crop(region = "exterior",
                        radiusOffset = -30) %>%
  cmcR::preProcess_crop(region = "interior",
                        radiusOffset = 200) %>%
  preProcess_removeTrend(statistic = "quantile",
                              tau = .5,
                              method = "fn") %>%
  cmcR::preProcess_gaussFilter() %>%
  x3ptools::sample_x3p()

knmComparisonFeatures <- purrr::map_dfr(seq(-30,30,by = 3),
                                        ~ comparison_allTogether(reference = fadul1.1,
                                                                 target = fadul2.1,
                                                                 numCells = 64,
                                                                 maxMissingProp = .85,
                                                                 theta = .))

knmComparisonFeatures_rev <- purrr::map_dfr(seq(-30,30,by = 3),
                                            ~ comparison_allTogether(reference = fadul2.1,
                                                                     target = fadul1.1,
                                                                     numCells = 64,
                                                                     maxMissingProp = .85,
                                                                     theta = .))

knmComparison_cmcs <- knmComparisonFeatures %>%
  mutate(originalMethodClassif = decision_CMC(cellIndex = cellIndex,
                                              x = x,
                                              y = y,
                                              theta = theta,
                                              corr = pairwiseCompCor,
                                              xThresh = 20,
                                              thetaThresh = 6,
                                              corrThresh = .5),
         highCMCClassif = decision_CMC(cellIndex = cellIndex,
                                              x = x,
                                              y = y,
                                              theta = theta,
                                              corr = pairwiseCompCor,
                                              xThresh = 20,
                                              thetaThresh = 6,
                                              corrThresh = .5,
                                              tau = 1))

knmComparison_cmcs_rev <- knmComparisonFeatures_rev %>%
  mutate(originalMethodClassif = decision_CMC(cellIndex = cellIndex,
                                              x = x,
                                              y = y,
                                              theta = theta,
                                              corr = pairwiseCompCor,
                                              xThresh = 20,
                                              thetaThresh = 6,
                                              corrThresh = .5),
         highCMCClassif = decision_CMC(cellIndex = cellIndex,
                                              x = x,
                                              y = y,
                                              theta = theta,
                                              corr = pairwiseCompCor,
                                              xThresh = 20,
                                              thetaThresh = 6,
                                              corrThresh = .5,
                                              tau = 1))


## ----echo=TRUE,eval=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------
## kmComparison_cmcs <- kmComparisonFeatures %>%
##   mutate(originalMethodClassif = decision_CMC(cellIndex = cellIndex,
##                                               x = x,
##                                               y = y,
##                                               theta = theta,
##                                               corr = pairwiseCompCor,
##                                               xThresh = 20,
##                                               thetaThresh = 6,
##                                               corrThresh = .5),
##          highCMCClassif = decision_CMC(cellIndex = cellIndex,
##                                        x = x,
##                                        y = y,
##                                        theta = theta,
##                                        corr = pairwiseCompCor,
##                                        xThresh = 20,
##                                        thetaThresh = 6,
##                                        corrThresh = .5,
##                                        tau = 1))


## ----echo=FALSE,warning=FALSE,message=FALSE,cache = F,fig.align='center',fig.pos='htbp',out.width="\\textwidth",fig.cap='\\label{fig:topVoteCMCPlot} CMC results for the comparison between Fadul 1-1 and Fadul 1-2 using the original method. The two plots in the top row show the 17 CMCs when Fadul 1-1 is treated as the ``reference" cartridge case to which Fadul 1-2 (the ``target") is compared. The second row shows the 18 CMCs when the roles are reversed. Red cells indicate where cells \\emph{not} identified as congruent achieve the maximum pairwise-complete correlation across all rotations of the target scan. '----

library(patchwork)

kmCMCPlot <- cmcR::cmcPlot(reference = fadul1.1,
                            target = fadul1.2,
                            reference_v_target_CMCs = kmComparison_cmcs,
                            target_v_reference_CMCs = kmComparison_cmcs_rev,
                            type = "faceted",
                            x3pNames = c("Fadul 1-1","Fadul 1-2"),
                            legend.quantiles = c(0,.01,.2,.5,.8,.99,1),
                            cell.colors = c("#a60b00","#1b03a3"),
                            cell.alpha = .15,
                            na.value = "grey100")



plt1 <- kmCMCPlot$originalMethodCMCs_reference_v_target

plt1$layers[[4]]$aes_params$size <- 1.5

plt2 <- kmCMCPlot$originalMethodCMCs_target_v_reference +
  theme(legend.text = element_text(size = 4.5),
        legend.title = element_text(size = 4.5))

plt2$layers[[4]]$aes_params$size <- 1.5

plt3 <- cowplot::get_legend(plt2)

plt1 <- plt1  +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",
        strip.text = element_blank()) +
  geom_text(data = data.frame(x3p = c("Fadul 1-1","Fadul 1-2"),
                              type = c("\n(Reference)","\n(Target)"),
                              x = c(1750,1750),
                              y = c(1500,1500)),
            aes(x = x,y = y,label = paste0(x3p,type)),size = 1.75)

plt2 <- plt2  +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",
        strip.text = element_blank()) +
  geom_text(data = data.frame(x3p = c("Fadul 1-1","Fadul 1-2"),
                              type = c("\n(Target)","\n(Reference)"),
                              x = c(1750,1750),
                              y = c(1500,1500)),
            aes(x = x,y = y,label = paste0(x3p,type)),size = 1.75)


plt <- (plt1 / plt2 / plt3) +
  plot_layout(ncol = 1,
              heights =  c(2,1.8,.5),
              widths = c(1,.9,1))

ggsave(filename = "cmcR_files/figure-latex/kmOriginalMethod.pdf",plot = plt)

invisible(knitr::plot_crop("cmcR_files/figure-latex/kmOriginalMethod.pdf",quiet = TRUE))

knitr::include_graphics(path = "cmcR_files/figure-latex/kmOriginalMethod.pdf")


## ----echo=FALSE,warning=FALSE,message=FALSE,cache = F, fig.align='center', fig.pos='htbp', out.width="\\textwidth", fig.cap='\\label{fig:highCMCPlot} Applying the High CMC method to the comparison of Fadul 1-1 and Fadul 1-2 results in 19 CMCs when Fadul 1-1 is treated as the reference (top) and 18 CMCs when Fadul 1-2 is treated as the reference (bottom). Although the individual comparisons do not yield considerably more CMCs than under the original CMC method, \\citet{tong_improved_2015} indicate that the High CMCs from both comparisons are combined as the final High CMC count (each cell is counted at most once). Combining the results means that the High CMC method tends to produce higher CMC counts than the original CMC method. In this example, the combined High CMC count is 24 CMCs.'----

plt1 <- kmCMCPlot$highCMC_reference_v_target


plt2 <- kmCMCPlot$highCMC_target_v_reference +
  theme(legend.text = element_text(size = 4.5),
        legend.title = element_text(size = 4.5))

plt3 <- cowplot::get_legend(plt2)

plt1$layers[[4]]$aes_params$size <- 1.5
plt2$layers[[4]]$aes_params$size <- 1.5

plt1 <- plt1  +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",
        strip.text = element_blank()) +
  geom_text(data = data.frame(x3p = c("Fadul 1-1","Fadul 1-2"),
                              type = c("\n(Reference)","\n(Target)"),
                              x = c(1750,1750),
                              y = c(1500,1500)),
            aes(x = x,y = y,label = paste0(x3p,type)),size = 1.75)

plt2 <- plt2  +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",
        strip.text = element_blank()) +
  geom_text(data = data.frame(x3p = c("Fadul 1-1","Fadul 1-2"),
                              type = c("\n(Target)","\n(Reference)"),
                              x = c(1750,1750),
                              y = c(1500,1500)),
            aes(x = x,y = y,label = paste0(x3p,type)),size = 1.75)


plt <- (plt1 / plt2 / plt3) +
  plot_layout(ncol = 1,
              heights =  c(2,2,.5),
              widths = c(1,1,1))

ggsave(filename = "cmcR_files/figure-latex/kmHighCMC.pdf",plot = plt)

invisible(knitr::plot_crop("cmcR_files/figure-latex/kmHighCMC.pdf"))

knitr::include_graphics(path = "cmcR_files/figure-latex/kmHighCMC.pdf")


## ----echo=FALSE,warning=FALSE,message=FALSE,cache = F, fig.align='center',fig.pos='htbp',out.width="\\textwidth", fig.cap='\\label{fig:knmCMCPlot} Applying both decision rules to the comparison between the non-match pair Fadul 1-1 and Fadul 2-1 results in 1 CMC under the original method (shown above) and 0 CMCs under the High CMC method (not shown). The seemingly random behavior of the red cells exemplifies the assumption that cells in a non-match comparison do not exhibit an observable pattern. Random chance should be the prevailing factor in classifying non-match cells as CMCs.'----

knmCMCPlot <- cmcR::cmcPlot(reference = fadul1.1,
                            target = fadul2.1,
                            reference_v_target_CMCs = knmComparison_cmcs,
                            target_v_reference_CMCs = knmComparison_cmcs_rev,
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
                            na.value = "grey100")

plt1 <- knmCMCPlot$originalMethodCMCs_reference_v_target +
  theme(legend.text = element_text(size = 4.5),
        legend.title = element_text(size = 4.5))


plt2 <- knmCMCPlot$originalMethodCMCs_target_v_reference

plt3 <- cowplot::get_legend(plt1)

plt1$layers[[4]]$aes_params$size <- 1.5
plt2$layers[[4]]$aes_params$size <- 1.5

plt1 <- plt1  +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",
        strip.text = element_blank()) +
  geom_text(data = data.frame(x3p = c("Fadul 1-1","Fadul 2-1"),
                              type = c("\n(Reference)","\n(Target)"),
                              x = c(1750,1750),
                              y = c(1500,1500)),
            aes(x = x,y = y,label = paste0(x3p,type)),size = 1.75)

plt2 <- plt2  +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",
        strip.text = element_blank()) +
  geom_text(data = data.frame(x3p = c("Fadul 1-1","Fadul 2-1"),
                              type = c("\n(Target)","\n(Reference)"),
                              x = c(1750,1750),
                              y = c(1500,1500)),
            aes(x = x,y = y,label = paste0(x3p,type)),size = 1.75)


plt <- (plt1 / plt2 / plt3) +
  plot_layout(ncol = 1,
              heights =  c(2,2,.5),
              widths = c(1,1,1))

ggsave(filename = "cmcR_files/figure-latex/knmOriginalMethod.pdf",plot = plt)

invisible(knitr::plot_crop("cmcR_files/figure-latex/knmOriginalMethod.pdf",quiet = TRUE))

knitr::include_graphics(path = "cmcR_files/figure-latex/knmOriginalMethod.pdf")


## ----echo=FALSE,fig.cap='\\label{fig:decisionRuleSensitivity_comparison} CMC count relative frequencies under the original method and the High CMC method for $T_{\\Delta x} = 20 = T_{\\Delta y}$ pixels, $T_{\\text{CCF}} = .5$, and $T_{\\theta} = 6$ degrees. AUC $= 1.00$ corresponds to perfect separation of the match and non-match CMC count distributions. We can see that, for this set of processing parameters, the High CMC method yields higher CMC counts for known matches that the original method while known non-matches have the same distribution under both methods.', fig.align='left',fig.pos='htbp',out.width='\\textwidth'----

load("data/cmcCountData.RData")

cmcCountData %>%
  ungroup() %>%
  filter(thetaThresh == 6 & 
           corThresh == .5 &
           transThresh == 20 &
           trendRemoved == TRUE) %>%
  group_by(thetaThresh,corThresh,transThresh,type) %>%
  mutate(n = n/sum(n),
         decisionRule = factor(decisionRule,levels = c("originalMethodCMCs","highCMCs"))) %>%
  ungroup() %>%
  rename(`Trans. Thresh` = transThresh,
         `CCF Thresh` = corThresh) %>%
  mutate(label = sprintf("AUC: %.2f\nVar. Ratio: %.2f", AUC, varRatio)) %>%
  ggplot() +
  geom_bar(aes(x = cmcCount,
               y = n,
               fill = type),
           stat = "identity",
           alpha = .7) +
  geom_label(aes(x = 15,
                 y = .25,
                 label = label),
             size = 4) +
  facet_grid(rows = vars(decisionRule),
             labeller = labeller(decisionRule = c("High CMC","Original Method") %>% set_names(c("highCMCs","originalMethodCMCs")))) +
  scale_fill_manual(values = c("#40B0A6","#E1BE6A")) +
  guides(fill = guide_legend(title = "Type",
                             override.aes = list(alpha = 1))) +
  theme_bw() + 
  theme(legend.position = "bottom",
        strip.text = element_text(size = 7)) +
  xlab("CMC Count") +
  ylab("Relative Frequency")


## ----echo=FALSE, fig.cap='\\label{fig:cmc_sensitivityScatter} Variance ratios under are plotted for different parameter settings. High variance ratios are indicative of a a good separation between CMC counts for known matching pairs and known-non matching pairs. The High CMC method generally performs better than the original method. Removing the trend during preprocessing, even though not explicitly described as a preprocessing step in the CMC papers, has a major impact on the effectiveness of the CMC method. In this setting, translation thresholds $T_x, T_y \\in [15,20]$, a rotation threshold $T_\\theta = 6$, and a CCF threshold $T_{\\text{CCF}} \\in [.4,5]$ lead to a separation of results. ',fig.align='left',fig.pos='htbp',out.width='\\textwidth'----

cmcCountData %>%
  mutate(trendRemoved = factor(trendRemoved)) %>%
  ggplot(aes(x = transThresh,
             y = varRatio,
             colour = corThresh)) +
  geom_point() +
  scale_colour_gradient(low = "#a1d99b",
                        high = "#00441b",
                        breaks = seq(.35,.6,by = .05)) +
  facet_grid(thetaThresh ~ decisionRule + trendRemoved,
             labeller = labeller(decisionRule = c("High CMC","Original Method") %>% set_names(c("highCMCs","originalMethodCMCs")),
                                 thetaThresh = c("Theta Thresh.: 3","Theta Thresh.: 6") %>% setNames(c(3,6)),
                                 trendRemoved = c("Trend Removed: TRUE","Trend Removed: FALSE") %>% setNames(c(TRUE,FALSE)))) +
  xlab("Translation Threshold") +
  ylab("variance ratio") +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(colour = guide_colorbar(title = "CCF Threshold",
                                 barwidth =  8,
                                 title.hjust = -1,
                                 title.vjust = .825,
                                 frame.colour = "black",
                                 ticks.colour = "black"))

