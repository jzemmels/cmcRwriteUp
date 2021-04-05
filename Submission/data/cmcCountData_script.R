## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(cmcR)
library(x3ptools)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

#create list of all x3p object filenames
fileList <- map(1:10,~ system(paste0("ls data/Fadul_",.,"/cc/Fadul*.x3p"),
                              intern=TRUE))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#create 63 match pair names
kmCombinationNames <- fileList %>%
  map(~ str_extract_all(string = .,
                        pattern = "Fadul\\ .*\\.x3p")) %>% #extract file name from location string
  map(~ str_remove_all(string = .,
                       pattern = ".x3p| ")) %>% #remove .x3p extension and whitespace
  map(~ combn(.,2,simplify = FALSE) %>%
        map(~ paste0(.,collapse="_"))) %>% #create vector of each pair of km file names
  unlist()

#create 717 non-match pair names
knmCombinationNames <- map2(fileList,1:10,
                            function(fileNames = .x,firearmNumber = .y){
                              bfNames <- fileNames %>%
                                map(~ str_extract_all(string = .,
                                                      pattern = "Fadul\\ .*\\.x3p")) %>% #extract file name from location string
                                map(~ str_remove_all(string = .,
                                                     pattern = ".x3p| ")) %>% #remove .x3p extension and whitespace
                                unlist() %>%
                                paste0("firearm",firearmNumber,"_",.)
                            })  %>%
  #      })  %>%
  purrr::flatten() %>%
  combn(2,simplify=FALSE) %>%
  map(function(pair){
    pair[[2]]
    # # names(pair)[1]
    firearmNumber <- pair[[1]] %>%
      stringr::str_extract("firearm[0-9]{1,2}")
    
    #Firearm 10 is difficult to tell apart from firearm 1 in the regex
    #expressions above, so it's handled separately:
    if(any(str_detect(pair,pattern = "Fadul10\\-.*$")) &
       any(str_detect(pair,pattern = "Fadul1\\-.*$")) | 
       any(str_detect(pair,pattern = "Fadul[EL]")) & #Fadul1 E and Fadul L are both from firearm 10
       any(str_detect(pair,pattern = "Fadul1\\-.*$")) | 
       any(str_detect(pair,pattern = "Fadul[EL]")) & #Fadul1 E and Fadul L are both from firearm 10
       any(str_detect(pair,pattern = "FadulF")) | 
       any(str_detect(pair,pattern = "Fadul10\\-.*$")) & #Fadul1 E and Fadul L are both from firearm 10
       any(str_detect(pair,pattern = "FadulF"))){
      return(pair)
    }
    if(str_detect(pair[[2]],pattern = firearmNumber)){
      return(NA)
    }
    else{
      return(pair)
    }
  })   %>%
  discard(~ all(is.na(.))) %>% 
  map_chr(~ paste0(.[[1]],"_",.[[2]]))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Preprocess scans -- returns a list of lists, one per firearm. Each of these
#firearm lists contains processed x3ps. Comment-out the preProcess_removeTrend
#function call to skip the de-trending step.
selectedBreechFaces <- fileList %>%
  purrr::map(function(kmSet){
    breechFaceNames <- kmSet %>%
      map(~ str_extract_all(string = .,
                            pattern = "Fadul\\ .*\\.x3p")) %>% #extract file name from location string
      map(~ str_remove_all(string = .,
                           pattern = ".x3p| ")) %>% #remove .x3p extension and whitespace
      unlist()
    
    selectedBreechFace <- kmSet %>% map(function(filePath){
      
      x3p <- x3ptools::read_x3p(filePath) %>%
        cmcR::preProcess_crop(region = "exterior",
                              radiusOffset = -30) %>%
        cmcR::preProcess_crop(region = "interior",
                              radiusOffset = 200) %>%
        cmcR::preProcess_removeTrend(statistic = "quantile",
                                     tau = .5,
                                     method = "fn") %>%
        cmcR::preProcess_gaussFilter() %>%
        x3ptools::x3p_sample(m=2)
      
      return(x3p)
      
    }) %>%
      setNames(breechFaceNames)
    
    return(selectedBreechFace)
    
  }) %>%
  setNames(paste0("firearm",1:10))

#Repeat preprocessing, but skipping de-trending step
selectedBreechFaces_noTrendRemoval <- fileList %>%
  purrr::map(function(kmSet){
    breechFaceNames <- kmSet %>%
      map(~ str_extract_all(string = .,
                            pattern = "Fadul\\ .*\\.x3p")) %>% #extract file name from location string
      map(~ str_remove_all(string = .,
                           pattern = ".x3p| ")) %>% #remove .x3p extension and whitespace
      unlist()
    
    selectedBreechFace <- kmSet %>% map(function(filePath){
      
      x3p <- x3ptools::read_x3p(filePath) %>%
        cmcR::preProcess_crop(region = "exterior",
                              radiusOffset = -30) %>%
        cmcR::preProcess_crop(region = "interior",
                              radiusOffset = 200) %>%
        # cmcR::preProcess_removeTrend(statistic = "quantile",
        #                              tau = .5,
        #                              method = "fn") %>%
        cmcR::preProcess_gaussFilter() %>%
        x3ptools::x3p_sample(m=2)
      
      return(x3p)
      
    }) %>%
      setNames(breechFaceNames)
    
    return(selectedBreechFace)
    
  }) %>%
  setNames(paste0("firearm",1:10))



## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Perform comparison procedure. This is currently so to allow for simultaneous
#comparison of multiple versions of the preprocessed scans (i.e., if you want to
#compare across different preprocessing conditions). Just add these to the
#list(selectedBreechFaces) object below.

#BE WARNED: comparing a set of 780 cartridge case pairs takes about 7 hours (on
#our machine, at least) because it computes the comparison results from 780
#comparisons: each comparing across 21 rotations, with about 30 cell/region
#pairs compared per theta value, and all of this is done twice for each pair.
#This means that there are about 1,000,000 cell/region comparisons performed in
#total (each of whic takes about .025 seconds)
comparisonData <- map(list(
  selectedBreechFaces,
  selectedBreechFaces_noTrendRemoval
),
function(selectBFs){
  # browser()
  
  kmPairCombinations <- selectBFs %>%
    map(function(imSet){
      kmPairs <- combn(x = imSet,m = 2,simplify = FALSE)
    }) %>%
    flatten() %>%
    setNames(kmCombinationNames)
  
  knmPairCombinations <- selectBFs %>%
    map2(.x = .,
         .y = fileList,
         function(kmSet,fileNames){
           bfNames <- fileNames %>%
             map(~ str_extract_all(string = .,
                                   pattern = "Fadul\\ .*\\.x3p")) %>% #extract file name from location string
             map(~ str_remove_all(string = .,
                                  pattern = ".x3p| ")) %>% #remove .x3p extension and whitespace
             unlist()
           
           setNames(kmSet,bfNames)
         }) %>%
    map2(.x = .,
         .y = names(.),
         function(im,imName){
           setNames(im, nm = names(im) %>%
                      map(~ paste0(.,"_",imName)))
         })  %>%
    purrr::flatten() %>%
    combn(2,simplify=FALSE) %>%
    map(function(pair){
      # names(pair)[1]
      firearmNumber <- names(pair)[1] %>%
        stringr::str_extract("firearm[0-9]{1,2}")
      
      # str_detect(names(pair)[2],pattern = firearmNumber)
      if(any(str_detect(names(pair),pattern = "10")) & 
         any(str_detect(names(pair),pattern = "1$"))){
        return(pair)
      }
      if(str_detect(names(pair)[2],pattern = firearmNumber)){
        return(NA)
      }
      else{
        return(pair)
      }
    }) %>%
    discard(~ all(is.na(.))) %>%
    setNames(knmCombinationNames)
  
  kmCorrs_64Cells_fftThenPairwise <- kmPairCombinations %>%
    map2_dfr(.x = .,
             .y = kmCombinationNames,
             function(pair,pairName){
               
               kmComparisonFeatures <- purrr::map_dfr(seq(-30,30,by = 3),
                                                      ~ comparison_allTogether(reference = pair[[1]],
                                                                               target = pair[[2]],
                                                                               numCells = 64,
                                                                               maxMissingProp = .9,
                                                                               theta = .)) %>%
                 mutate(direction = "1to2",
                        type = "match",
                        pairName = pairName)
               
               kmComparisonFeatures_rev <- purrr::map_dfr(seq(-30,30,by = 3),
                                                          ~ comparison_allTogether(reference = pair[[2]],
                                                                                   target = pair[[1]],
                                                                                   numCells = 64,
                                                                                   maxMissingProp = .9,
                                                                                   theta = .))  %>%
                 mutate(direction = "2to1",
                        type = "match",
                        pairName = pairName)
               
               return(bind_rows(kmComparisonFeatures,kmComparisonFeatures_rev))
               
             })
  
  print("fftThenPairwise KM pairs done")
  
  knmCorrs_64Cells_fftThenPairwise <- knmPairCombinations %>%
    map2_dfr(.x = .,
             .y = knmCombinationNames,
             function(pair,pairName){
               
               knmComparisonFeatures <- purrr::map_dfr(seq(-30,30,by = 3),
                                                       ~ comparison_allTogether(reference = pair[[1]],
                                                                                target = pair[[2]],
                                                                                numCells = 64,
                                                                                maxMissingProp = .9,
                                                                                theta = .)) %>%
                 mutate(direction = "1to2",
                        type = "non-match",
                        pairName = pairName)
               
               knmComparisonFeatures_rev <- purrr::map_dfr(seq(-30,30,by = 3),
                                                           ~ comparison_allTogether(reference = pair[[2]],
                                                                                    target = pair[[1]],
                                                                                    numCells = 64,
                                                                                    maxMissingProp = .9,
                                                                                    theta = .))  %>%
                 mutate(direction = "2to1",
                        type = "non-match",
                        pairName = pairName)
               
               return(bind_rows(knmComparisonFeatures,knmComparisonFeatures_rev))
               
             })
  
  print("fftThenPairwise KNM pairs done")
  
  return(bind_rows(kmCorrs_64Cells_fftThenPairwise,
                   knmCorrs_64Cells_fftThenPairwise))
})


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Calculate the CMCs across different CCF, translation, and theta thresholds.
#This can be repeated for any versions of the preprocessed scans.

#BE WARNED: This will also take a while to run (although less than the
#comparison procedure). Also, the cmcResults data frame will be about 28 GB
#(about 140,000,000*2 rows) due to the large number of combinations of thresholds
#considered and the fact that we're considering 2 sets of preprocessed scans
cmcResults <- comparisonData %>%
  map2_dfr(.x = .,
           .y = c(TRUE,FALSE),
           function(comparisonDF,trendRemoved){
             comparisonDF%>%
               group_by(pairName,direction) %>%
               group_split() %>%
               map_dfr(function(pairCorrResults){
                 
                 map_dfr(purrr::cross(.l = list("ccf" = seq(.35,.6,by = .025),
                                                "translation" = c(5,10,15,20,25,30),
                                                "theta" = c(3,6))),
                         function(thresholds){
                           
                           pairCorrResults %>%
                             mutate(originalMethodClassif = decision_CMC(cellIndex,x,y,theta,pairwiseCompCor,
                                                                         xThresh = thresholds$translation,thetaThresh = thresholds$theta,corrThresh = thresholds$ccf),
                                    highCMCClassif = decision_CMC(cellIndex,x,y,theta,pairwiseCompCor,
                                                                  xThresh = thresholds$translation,thetaThresh = thresholds$theta,corrThresh = thresholds$ccf,
                                                                  tau = 1),
                                    transThresh = thresholds$translation,
                                    thetaThresh = thresholds$theta,
                                    corThresh = thresholds$ccf) %>%
                             return()
                         }) %>%
                   return()
               }) %>%
               mutate(trendRemoved = trendRemoved)
           })


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Combine the CMC results from both comparison directions for each set of
# translation thresholds and each cartridge case pair. This will result in final
# CMC counts under the original method Song (2013) and the High CMC method for
# each cartridge case pair. Also, calculate the variance ratios and AUCs
# associated with each set of processing settings
cmcResults_combinedDirections <- cmcResults %>%
  group_by(pairName,transThresh,corThresh,thetaThresh,trendRemoved) %>%
  group_split() %>%
  purrr::map_dfr(function(pairResults){
    
    cmcs <- cmcR::decision_combineDirections(pairResults %>%
                                               filter(direction == "1to2"),
                                             pairResults %>%
                                               filter(direction == "2to1"),
                                             thetaThresh = 6)
    
    originalMethodCount <- min(c(nrow(cmcs$originalMethodCMCs[[1]]),nrow(cmcs$originalMethodCMCs[[2]])))
    
    pairResults %>%
      select(c(pairName,transThresh,corThresh,thetaThresh)) %>%
      distinct() %>%
      mutate(originalMethodCMCs = originalMethodCount,
             highCMCs = nrow(cmcs$highCMCs))
  }) %>%
  mutate(type = ifelse(stringr::str_detect(pairName,"firearm"),"non-match","match")) %>%
  pivot_longer(cols = c("originalMethodCMCs","highCMCs"),
               names_to = "decisionRule",
               values_to = "cmcCount") %>%
  ungroup() %>%
  group_by(transThresh,corThresh,thetaThresh,trendRemoved,decisionRule) %>%
  group_split() %>%
  map_dfr(function(paramConditionedResults){
    
    decisionRuleAUC <- paramConditionedResults %>%
      pROC::roc(response = type,
                predictor = cmcCount,
                levels = c("non-match","match"),
                quiet = TRUE)
    
    paramConditionedResults %>%
      mutate(AUC = as.numeric(decisionRuleAUC$auc))
  }) %>%
  ungroup()


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Summarize the AUCs and variance ratios for each set of thresholds.
cmcCountData <-  cmcResults_combinedDirections %>%
  group_by(transThresh,corThresh,thetaThresh,trendRemoved,decisionRule,type,cmcCount) %>%
  summarise(n = n(),
            varRatio = unique(varRatio),
            AUC = unique(AUC)) %>%
  ungroup()

