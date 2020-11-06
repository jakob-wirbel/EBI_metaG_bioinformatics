---
title: "Comparative Metagenome Analysis - EBI Metagenomics Bioinformatics"
output:
  html_document:
    toc: yes
    toc_float: yes
    number_sections: yes
    keep_md: true
    df_print: paged
author: "Jakob Wirbel and Georg Zeller"
---

In this practical, we are going to explore statistical testing methods for 
metagenomic data, how to visualize the results, and how to train machine 
learning models using the `SIAMCAT` package.



# Setup

## Preparing the R environment 

In order to get started, we should first prepare our `R` environment and load
the packages we will need later on. Additionally, the data used in this 
practical are stored on the EMBL servers and we can set the base path for the
downloads.


```r
library("tidyverse") # for general data wrangling and plotting
library("SIAMCAT")   # for statistical and ML analyses

data.location <- 'https://embl.de/download/zeller/metaG_course/'
```

## Loading the Data

In this practical, we are going to have a look at the data from
[Zeller et al. _MSB_ 2014](https://doi.org/10.15252/msb.20145645). In this 
study, the authors recruited patients with **colorectal cancer (CRC)** 
and **healthy controls (CTR)** and performed shotgun metagenomic sequencing 
of fecal samples. The raw data have already been pre-processed and analyzed 
with the [mOTUs2](https://doi.org/10.1038/s41467-019-08844-4) taxonomic 
profiler.

### Features

First, we are going to load the taxonomic profiles and store them as a matrix.


```r
fn.feat.fr  <- paste0(data.location, '/motus_profiles/FR-CRC.motus')
tax.profiles <- read.table(fn.feat.fr, sep = '\t', quote = '', 
                           comment.char = '', skip = 61, 
                           stringsAsFactors = FALSE, check.names = FALSE, 
                           row.names = 1, header = TRUE)
tax.profiles <- as.matrix(tax.profiles)
```

The taxonomic profiles contain absolute counts and can easily be transformed
into relative abundances using the `prop.table` function:


```r
rel.tax.profiles <- prop.table(tax.profiles, 2)
```

### Metadata

Additionally, we also need the information which sample belongs to which group.
Therefore, we are loading the metadata table as well:


```r
fn.meta.fr  <- paste0(data.location, '/metadata/meta_FR-CRC.tsv')
df.meta <- read_tsv(fn.meta.fr)
df.meta
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Sample_ID"],"name":[1],"type":["chr"],"align":["left"]},{"label":["External_ID"],"name":[2],"type":["chr"],"align":["left"]},{"label":["Age"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Gender"],"name":[4],"type":["chr"],"align":["left"]},{"label":["BMI"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["Country"],"name":[6],"type":["chr"],"align":["left"]},{"label":["AJCC_stage"],"name":[7],"type":["chr"],"align":["left"]},{"label":["TNM_stage"],"name":[8],"type":["chr"],"align":["left"]},{"label":["Localization"],"name":[9],"type":["chr"],"align":["left"]},{"label":["Study"],"name":[10],"type":["chr"],"align":["left"]},{"label":["Sampling_rel_to_colonoscopy"],"name":[11],"type":["chr"],"align":["left"]},{"label":["FOBT"],"name":[12],"type":["chr"],"align":["left"]},{"label":["Group"],"name":[13],"type":["chr"],"align":["left"]},{"label":["Library_Size"],"name":[14],"type":["dbl"],"align":["right"]}],"data":[{"1":"CCIS00146684ST-4-0","2":"FR-726","3":"72","4":"F","5":"25","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"35443944"},{"1":"CCIS00281083ST-3-0","2":"FR-060","3":"53","4":"M","5":"32","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"19307896"},{"1":"CCIS02124300ST-4-0","2":"FR-568","3":"35","4":"M","5":"23","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"42141246"},{"1":"CCIS02379307ST-4-0","2":"FR-828","3":"67","4":"M","5":"28","6":"FRA","7":"I","8":"T1N0M0","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"4829533"},{"1":"CCIS02856720ST-4-0","2":"FR-027","3":"74","4":"M","5":"27","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CTR","14":"34294675"},{"1":"CCIS03473770ST-4-0","2":"FR-192","3":"29","4":"M","5":"24","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"20262319"},{"1":"CCIS03857607ST-4-0","2":"FR-349","3":"61","4":"M","5":"23","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"29917559"},{"1":"CCIS05314658ST-4-0","2":"FR-169","3":"65","4":"F","5":"24","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"28453362"},{"1":"CCIS06260551ST-3-0","2":"FR-200","3":"58","4":"M","5":"24","6":"FRA","7":"IV","8":"T3N1M1","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"15598843"},{"1":"CCIS07277498ST-4-0","2":"FR-276","3":"63","4":"M","5":"NA","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"NA","13":"CTR","14":"31510809"},{"1":"CCIS07539127ST-4-0","2":"FR-460","3":"77","4":"F","5":"22","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"32025452"},{"1":"CCIS07648107ST-4-0","2":"FR-053","3":"62","4":"F","5":"21","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"30986245"},{"1":"CCIS08668806ST-3-0","2":"FR-214","3":"63","4":"M","5":"26","6":"FRA","7":"NA","8":"NA","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"NAA","14":"8467141"},{"1":"CCIS09568613ST-4-0","2":"FR-400","3":"67","4":"M","5":"21","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"35354460"},{"1":"CCIS10706551ST-3-0","2":"FR-193","3":"25","4":"M","5":"24","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"28855308"},{"1":"CCIS10793554ST-4-0","2":"FR-606","3":"64","4":"M","5":"26","6":"FRA","7":"NA","8":"NA","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"NAA","14":"18941421"},{"1":"CCIS11015875ST-4-0","2":"FR-788","3":"62","4":"M","5":"26","6":"FRA","7":"III","8":"T3N1M0","9":"Rectum","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"26824567"},{"1":"CCIS11019776ST-4-0","2":"FR-548","3":"58","4":"M","5":"22","6":"FRA","7":"NA","8":"NA","9":"LC/RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"ADA","14":"22249555"},{"1":"CCIS11354283ST-4-0","2":"FR-212","3":"67","4":"M","5":"27","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"28417199"},{"1":"CCIS11362406ST-4-0","2":"FR-398","3":"44","4":"F","5":"26","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"46574413"},{"1":"CCIS11558985ST-4-0","2":"FR-902","3":"56","4":"F","5":"NA","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"17480615"},{"1":"CCIS12370844ST-4-0","2":"FR-823","3":"81","4":"F","5":"39","6":"FRA","7":"I","8":"T2N0M0","9":"Rectum","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"1981993"},{"1":"CCIS12656533ST-4-0","2":"FR-654","3":"51","4":"M","5":"30","6":"FRA","7":"IV","8":"T2N1M1","9":"Rectum","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"21293790"},{"1":"CCIS13047523ST-4-0","2":"FR-003","3":"70","4":"M","5":"22","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"29557894"},{"1":"CCIS14449628ST-4-0","2":"FR-404","3":"59","4":"F","5":"20","6":"FRA","7":"III","8":"T4N1M0","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"22367285"},{"1":"CCIS15704761ST-4-0","2":"FR-901","3":"56","4":"F","5":"NA","6":"FRA","7":"IV","8":"T4N1M1","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"19267185"},{"1":"CCIS15794887ST-4-0","2":"FR-024","3":"37","4":"F","5":"18","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"19330062"},{"1":"CCIS16326685ST-4-0","2":"FR-542","3":"46","4":"F","5":"29","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"24000063"},{"1":"CCIS16383318ST-4-0","2":"FR-139","3":"61","4":"F","5":"24","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"36238309"},{"1":"CCIS16561622ST-4-0","2":"FR-030","3":"54","4":"M","5":"26","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"41280291"},{"1":"CCIS17669415ST-4-0","2":"FR-125","3":"72","4":"F","5":"37","6":"FRA","7":"I","8":"T1N0M0","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"30538779"},{"1":"CCIS19142497ST-3-0","2":"FR-051","3":"68","4":"M","5":"23","6":"FRA","7":"NA","8":"NA","9":"LC/RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"ADA","14":"12933986"},{"1":"CCIS20795251ST-4-0","2":"FR-221","3":"71","4":"F","5":"23","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"25358270"},{"1":"CCIS21126322ST-4-0","2":"FR-121","3":"59","4":"M","5":"25","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"36157541"},{"1":"CCIS21278152ST-4-0","2":"FR-414","3":"63","4":"M","5":"25","6":"FRA","7":"II","8":"T3N0M0","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"23458416"},{"1":"CCIS22275061ST-4-0","2":"FR-191","3":"62","4":"M","5":"29","6":"FRA","7":"NA","8":"NA","9":"Rectum","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"ADA","14":"42358438"},{"1":"CCIS22416007ST-4-0","2":"FR-815","3":"82","4":"M","5":"24","6":"FRA","7":"IV","8":"T4N1M1","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"42178363"},{"1":"CCIS22906510ST-20-0","2":"FR-826","3":"77","4":"M","5":"23","6":"FRA","7":"NA","8":"NA","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"NAA","14":"40587980"},{"1":"CCIS22958137ST-20-0","2":"FR-506","3":"60","4":"M","5":"26","6":"FRA","7":"IV","8":"T4N1M1","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"46226004"},{"1":"CCIS23164343ST-4-0","2":"FR-316","3":"62","4":"M","5":"23","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"32281417"},{"1":"CCIS24254057ST-4-0","2":"FR-767","3":"69","4":"F","5":"24","6":"FRA","7":"IV","8":"T4N1M1","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"61932153"},{"1":"CCIS24898163ST-4-0","2":"FR-282","3":"59","4":"M","5":"31","6":"FRA","7":"NA","8":"NA","9":"Rectum","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"NAA","14":"5647333"},{"1":"CCIS25399172ST-4-0","2":"FR-419","3":"60","4":"M","5":"27","6":"FRA","7":"NA","8":"NA","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"NAA","14":"13361371"},{"1":"CCIS27304052ST-3-0","2":"FR-001","3":"52","4":"F","5":"20","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"19441026"},{"1":"CCIS27927933ST-4-0","2":"FR-115","3":"72","4":"M","5":"24","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"44183270"},{"1":"CCIS28384594ST-4-0","2":"FR-382","3":"89","4":"M","5":"24","6":"FRA","7":"NA","8":"NA","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"ADA","14":"25528379"},{"1":"CCIS29210128ST-4-0","2":"FR-166","3":"64","4":"M","5":"27","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"21601905"},{"1":"CCIS29688262ST-20-0","2":"FR-305","3":"50","4":"M","5":"25","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"38854418"},{"1":"CCIS31434951ST-20-0","2":"FR-643","3":"66","4":"M","5":"24","6":"FRA","7":"NA","8":"NA","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"NAA","14":"13840048"},{"1":"CCIS32105356ST-4-0","2":"FR-510","3":"53","4":"M","5":"30","6":"FRA","7":"NA","8":"NA","9":"LC/RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"NAA","14":"23774880"},{"1":"CCIS32452666ST-4-0","2":"FR-105","3":"68","4":"M","5":"26","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"25019747"},{"1":"CCIS33816588ST-4-0","2":"FR-652","3":"64","4":"M","5":"21","6":"FRA","7":"III","8":"T3N1M0","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"30493004"},{"1":"CCIS34055159ST-4-0","2":"FR-430","3":"80","4":"M","5":"29","6":"FRA","7":"I","8":"T2N0M0","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"21250090"},{"1":"CCIS34604008ST-4-0","2":"FR-129","3":"63","4":"F","5":"27","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"33724649"},{"1":"CCIS35092938ST-4-0","2":"FR-500","3":"62","4":"F","5":"21","6":"FRA","7":"NA","8":"NA","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"NAA","14":"20654081"},{"1":"CCIS35100175ST-4-0","2":"FR-812","3":"73","4":"M","5":"40","6":"FRA","7":"I","8":"T1N0M0","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"25119675"},{"1":"CCIS36699628ST-4-0","2":"FR-142","3":"68","4":"F","5":"25","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"31159348"},{"1":"CCIS36797902ST-4-0","2":"FR-198","3":"63","4":"M","5":"25","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"19681060"},{"1":"CCIS37250421ST-4-0","2":"FR-496","3":"59","4":"F","5":"21","6":"FRA","7":"NA","8":"NA","9":"Rectum","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"NAA","14":"22577160"},{"1":"CCIS38765456ST-20-0","2":"FR-723","3":"79","4":"F","5":"22","6":"FRA","7":"IV","8":"T4N1M1","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"38274855"},{"1":"CCIS40244499ST-3-0","2":"FR-202","3":"87","4":"F","5":"32","6":"FRA","7":"II","8":"T4N0M0","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"24505848"},{"1":"CCIS41222843ST-4-0","2":"FR-298","3":"73","4":"F","5":"24","6":"FRA","7":"III","8":"T3N1M0","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"42370435"},{"1":"CCIS41288781ST-4-0","2":"FR-817","3":"74","4":"F","5":"19","6":"FRA","7":"II","8":"T3N0M0","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"20807788"},{"1":"CCIS41548810ST-4-0","2":"FR-451","3":"78","4":"M","5":"19","6":"FRA","7":"IV","8":"T3N0M1","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"25909920"},{"1":"CCIS41692898ST-4-0","2":"FR-558","3":"67","4":"F","5":"28","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"35216845"},{"1":"CCIS41806458ST-4-0","2":"FR-664","3":"59","4":"F","5":"NA","6":"FRA","7":"NA","8":"NA","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"NAA","14":"27606591"},{"1":"CCIS44093303ST-4-0","2":"FR-721","3":"62","4":"F","5":"20","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CTR","14":"34539245"},{"1":"CCIS44676181ST-4-0","2":"FR-161","3":"48","4":"F","5":"21","6":"FRA","7":"NA","8":"NA","9":"Rectum","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"NAA","14":"19955843"},{"1":"CCIS44743950ST-4-0","2":"FR-596","3":"66","4":"M","5":"22","6":"FRA","7":"NA","8":"NA","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"NAA","14":"38044446"},{"1":"CCIS44757994ST-4-0","2":"FR-770","3":"75","4":"M","5":"37","6":"FRA","7":"III","8":"T3N1M0","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"22702296"},{"1":"CCIS45571137ST-3-0","2":"FR-116","3":"55","4":"F","5":"20","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"15000419"},{"1":"CCIS45793747ST-4-0","2":"FR-684","3":"67","4":"M","5":"25","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"39708809"},{"1":"CCIS46047672ST-4-0","2":"FR-825","3":"65","4":"M","5":"24","6":"FRA","7":"IV","8":"T4N0M1","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"18277769"},{"1":"CCIS46467422ST-4-0","2":"FR-830","3":"54","4":"F","5":"22","6":"FRA","7":"III","8":"T3N1M0","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"5534486"},{"1":"CCIS47284573ST-4-0","2":"FR-667","3":"67","4":"M","5":"25","6":"FRA","7":"III","8":"T4N1M0","9":"Rectum","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"18352816"},{"1":"CCIS47745468ST-4-0","2":"FR-054","3":"63","4":"M","5":"27","6":"FRA","7":"NA","8":"NA","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"ADA","14":"36205688"},{"1":"CCIS48174381ST-4-0","2":"FR-751","3":"66","4":"M","5":"25","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"24187553"},{"1":"CCIS48507077ST-4-0","2":"FR-507","3":"72","4":"F","5":"32","6":"FRA","7":"NA","8":"NA","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"NAA","14":"26359362"},{"1":"CCIS48579360ST-4-0","2":"FR-829","3":"52","4":"F","5":"23","6":"FRA","7":"NA","8":"NA","9":"Rectum","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"NAA","14":"10676683"},{"1":"CCIS48725289ST-4-0","2":"FR-170","3":"62","4":"F","5":"21","6":"FRA","7":"I","8":"T1N0M0","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"15677334"},{"1":"CCIS50003399ST-4-0","2":"FR-194","3":"66","4":"F","5":"28","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"29725483"},{"1":"CCIS50148151ST-4-0","2":"FR-503","3":"87","4":"F","5":"15","6":"FRA","7":"III","8":"T2N1M0","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"24166049"},{"1":"CCIS50369211ST-20-0","2":"FR-215","3":"53","4":"F","5":"33","6":"FRA","7":"NA","8":"NA","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"NAA","14":"9842615"},{"1":"CCIS50471204ST-4-0","2":"FR-213","3":"61","4":"F","5":"23","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"38401964"},{"1":"CCIS50561855ST-4-0","2":"FR-474","3":"60","4":"F","5":"26","6":"FRA","7":"NA","8":"NA","9":"LC/RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"ADA","14":"14986182"},{"1":"CCIS51595129ST-4-0","2":"FR-399","3":"63","4":"F","5":"22","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"27642810"},{"1":"CCIS51667829ST-4-0","2":"FR-113","3":"63","4":"M","5":"26","6":"FRA","7":"NA","8":"NA","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"NAA","14":"35840774"},{"1":"CCIS52234160ST-4-0","2":"FR-473","3":"64","4":"M","5":"28","6":"FRA","7":"NA","8":"NA","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"ADA","14":"22346644"},{"1":"CCIS52370277ST-4-0","2":"FR-827","3":"53","4":"M","5":"25","6":"FRA","7":"IV","8":"T4N1M1","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"26193276"},{"1":"CCIS53043478ST-4-0","2":"FR-551","3":"72","4":"M","5":"22","6":"FRA","7":"IV","8":"T3N1M1","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"17962548"},{"1":"CCIS53355328ST-4-0","2":"FR-293","3":"76","4":"M","5":"NA","6":"FRA","7":"I","8":"T2N0M0","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"17433179"},{"1":"CCIS53557295ST-4-0","2":"FR-666","3":"45","4":"M","5":"24","6":"FRA","7":"NA","8":"NA","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"NAA","14":"14671389"},{"1":"CCIS54027808ST-4-0","2":"FR-835","3":"52","4":"M","5":"33","6":"FRA","7":"NA","8":"NA","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"NAA","14":"22030805"},{"1":"CCIS55230578ST-4-0","2":"FR-734","3":"64","4":"M","5":"30","6":"FRA","7":"IV","8":"T3N1M1","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"27804303"},{"1":"CCIS55531770ST-4-0","2":"FR-790","3":"45","4":"F","5":"24","6":"FRA","7":"NA","8":"NA","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"NAA","14":"6745142"},{"1":"CCIS56503244ST-3-0","2":"FR-162","3":"68","4":"M","5":"21","6":"FRA","7":"NA","8":"NA","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"ADA","14":"13595352"},{"1":"CCIS58234805ST-4-0","2":"FR-268","3":"70","4":"F","5":"26","6":"FRA","7":"III","8":"T3N1M0","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"11012918"},{"1":"CCIS59132091ST-4-0","2":"FR-617","3":"80","4":"F","5":"23","6":"FRA","7":"I","8":"T1N0M0","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"47653407"},{"1":"CCIS59903910ST-4-0","2":"FR-820","3":"76","4":"F","5":"34","6":"FRA","7":"NA","8":"NA","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"NAA","14":"24858698"},{"1":"CCIS61287323ST-4-0","2":"FR-376","3":"76","4":"M","5":"25","6":"FRA","7":"II","8":"T3N0M0","9":"Rectum","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"35663674"},{"1":"CCIS62605362ST-3-0","2":"FR-504","3":"56","4":"F","5":"23","6":"FRA","7":"II","8":"T3N0M0","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"58613057"},{"1":"CCIS62794166ST-4-0","2":"FR-780","3":"68","4":"M","5":"23","6":"FRA","7":"IV","8":"T3N1M1","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"18137582"},{"1":"CCIS63448910ST-4-0","2":"FR-390","3":"61","4":"F","5":"34","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"29010067"},{"1":"CCIS63468405ST-4-0","2":"FR-132","3":"69","4":"F","5":"25","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"31230015"},{"1":"CCIS63910149ST-4-0","2":"FR-719","3":"49","4":"F","5":"23","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"7403037"},{"1":"CCIS64773582ST-4-0","2":"FR-195","3":"63","4":"F","5":"23","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"43574702"},{"1":"CCIS64785924ST-20-0","2":"FR-173","3":"59","4":"F","5":"25","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"5953194"},{"1":"CCIS65479369ST-4-0","2":"FR-772","3":"69","4":"F","5":"25","6":"FRA","7":"IV","8":"T3N1M1","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"7732781"},{"1":"CCIS70398272ST-4-0","2":"FR-792","3":"55","4":"M","5":"24","6":"FRA","7":"NA","8":"NA","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"NAA","14":"36295986"},{"1":"CCIS71301801ST-4-0","2":"FR-281","3":"50","4":"F","5":"24","6":"FRA","7":"I","8":"T1N0M0","9":"Rectum","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"45032701"},{"1":"CCIS71578391ST-4-0","2":"FR-187","3":"70","4":"M","5":"25","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"34251424"},{"1":"CCIS72607085ST-4-0","2":"FR-824","3":"74","4":"F","5":"26","6":"FRA","7":"IV","8":"T3NxM1","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"21729304"},{"1":"CCIS74239020ST-4-0","2":"FR-156","3":"62","4":"M","5":"24","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"28128807"},{"1":"CCIS74726977ST-3-0","2":"FR-026","3":"66","4":"M","5":"24","6":"FRA","7":"NA","8":"NA","9":"LC/RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"NAA","14":"62926247"},{"1":"CCIS76563044ST-4-0","2":"FR-672","3":"68","4":"M","5":"36","6":"FRA","7":"NA","8":"NA","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"NAA","14":"35200771"},{"1":"CCIS76845094ST-20-0","2":"FR-450","3":"44","4":"M","5":"20","6":"FRA","7":"IV","8":"T3N1M1","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"72543335"},{"1":"CCIS77100899ST-4-0","2":"FR-211","3":"73","4":"M","5":"26","6":"FRA","7":"NA","8":"NA","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"ADA","14":"32162657"},{"1":"CCIS77252613ST-4-0","2":"FR-783","3":"65","4":"M","5":"26","6":"FRA","7":"I","8":"T1N0M0","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"43032257"},{"1":"CCIS78100604ST-4-0","2":"FR-728","3":"64","4":"F","5":"30","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"36588707"},{"1":"CCIS78318719ST-4-0","2":"FR-768","3":"69","4":"F","5":"30","6":"FRA","7":"II","8":"T3N0M0","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"26554667"},{"1":"CCIS79210440ST-3-0","2":"FR-039","3":"65","4":"M","5":"30","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CTR","14":"17338108"},{"1":"CCIS80834637ST-4-0","2":"FR-716","3":"60","4":"F","5":"28","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"44830160"},{"1":"CCIS81139242ST-4-0","2":"FR-118","3":"64","4":"F","5":"20","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"29841558"},{"1":"CCIS81710917ST-20-0","2":"FR-539","3":"49","4":"M","5":"20","6":"FRA","7":"NA","8":"NA","9":"Rectum","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"NAA","14":"48035727"},{"1":"CCIS81735969ST-20-0","2":"FR-312","3":"68","4":"M","5":"26","6":"FRA","7":"NA","8":"NA","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"ADA","14":"34764065"},{"1":"CCIS81887263ST-4-0","2":"FR-302","3":"59","4":"M","5":"27","6":"FRA","7":"I","8":"T2N0M0","9":"Rectum","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"43180950"},{"1":"CCIS82146115ST-4-0","2":"FR-628","3":"51","4":"M","5":"24","6":"FRA","7":"IV","8":"T3N1M1","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"43960694"},{"1":"CCIS82507866ST-3-0","2":"FR-040","3":"57","4":"M","5":"24","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"18262334"},{"1":"CCIS82944710ST-20-0","2":"FR-730","3":"38","4":"F","5":"22","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"4215880"},{"1":"CCIS83445808ST-4-0","2":"FR-495","3":"66","4":"M","5":"29","6":"FRA","7":"NA","8":"NA","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"NAA","14":"19501458"},{"1":"CCIS83574003ST-4-0","2":"FR-241","3":"60","4":"M","5":"24","6":"FRA","7":"NA","8":"NA","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"ADA","14":"38304657"},{"1":"CCIS83870198ST-4-0","2":"FR-344","3":"48","4":"F","5":"24","6":"FRA","7":"I","8":"T1N0M0","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"33343498"},{"1":"CCIS84543192ST-4-0","2":"FR-722","3":"45","4":"F","5":"28","6":"FRA","7":"I","8":"T2N0M0","9":"Rectum","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"14973046"},{"1":"CCIS85214191ST-3-0","2":"FR-208","3":"63","4":"M","5":"22","6":"FRA","7":"IV","8":"T3N1M1","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"26950900"},{"1":"CCIS87116798ST-4-0","2":"FR-310","3":"85","4":"F","5":"20","6":"FRA","7":"IV","8":"T3N1M1","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"22643954"},{"1":"CCIS87167916ST-4-0","2":"FR-223","3":"79","4":"M","5":"30","6":"FRA","7":"IV","8":"T4N1M1","9":"Rectum","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"20100401"},{"1":"CCIS87252800ST-4-0","2":"FR-505","3":"73","4":"M","5":"17","6":"FRA","7":"I","8":"T2N0M0","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"11753596"},{"1":"CCIS87290971ST-4-0","2":"FR-196","3":"71","4":"M","5":"29","6":"FRA","7":"NA","8":"NA","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"ADA","14":"27257848"},{"1":"CCIS87605453ST-4-0","2":"FR-328","3":"73","4":"F","5":"26","6":"FRA","7":"II","8":"T3N0M0","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"39141969"},{"1":"CCIS88007743ST-4-0","2":"FR-393","3":"67","4":"F","5":"23","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"29101486"},{"1":"CCIS88317640ST-4-0","2":"FR-759","3":"74","4":"M","5":"27","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"33771666"},{"1":"CCIS90164298ST-4-0","2":"FR-294","3":"84","4":"M","5":"23","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"30094526"},{"1":"CCIS90166425ST-4-0","2":"FR-449","3":"78","4":"F","5":"23","6":"FRA","7":"NA","8":"NA","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"ADA","14":"27271452"},{"1":"CCIS90443472ST-4-0","2":"FR-249","3":"63","4":"M","5":"30","6":"FRA","7":"NA","8":"NA","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"ADA","14":"31948167"},{"1":"CCIS90903952ST-3-0","2":"FR-218","3":"71","4":"F","5":"22","6":"FRA","7":"NA","8":"NA","9":"LC/RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"ADA","14":"22788163"},{"1":"CCIS91228662ST-4-0","2":"FR-275","3":"63","4":"M","5":"31","6":"FRA","7":"I","8":"T1N0M0","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"44044852"},{"1":"CCIS93040568ST-20-0","2":"FR-682","3":"65","4":"M","5":"30","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"6436717"},{"1":"CCIS94417875ST-3-0","2":"FR-110","3":"59","4":"F","5":"25","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"18113211"},{"1":"CCIS94496512ST-4-0","2":"FR-229","3":"64","4":"M","5":"28","6":"FRA","7":"NA","8":"NA","9":"Rectum","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"NAA","14":"17750885"},{"1":"CCIS94603952ST-4-0","2":"FR-025","3":"80","4":"F","5":"21","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"32451093"},{"1":"CCIS95097901ST-4-0","2":"FR-696","3":"52","4":"M","5":"24","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"24039454"},{"1":"CCIS95409808ST-4-0","2":"FR-152","3":"63","4":"F","5":"21","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"31515564"},{"1":"CCIS96387239ST-4-0","2":"FR-626","3":"66","4":"M","5":"20","6":"FRA","7":"NA","8":"NA","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"NAA","14":"15901300"},{"1":"CCIS98482370ST-3-0","2":"FR-052","3":"53","4":"F","5":"30","6":"FRA","7":"NA","8":"NA","9":"NA","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CTR","14":"18276972"},{"1":"CCIS98512455ST-4-0","2":"FR-459","3":"63","4":"M","5":"29","6":"FRA","7":"IV","8":"T4N1M1","9":"RC","10":"FR-CRC","11":"BEFORE","12":"Negative","13":"CRC","14":"20597661"},{"1":"CCIS98832363ST-4-0","2":"FR-552","3":"55","4":"F","5":"25","6":"FRA","7":"III","8":"T3N1M0","9":"LC","10":"FR-CRC","11":"BEFORE","12":"Positive","13":"CRC","14":"30762491"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
table(df.meta$Group)
## 
## ADA CRC CTR NAA 
##  15  53  61  27
```

We are interested in the comparison between control samples (`CTR`) and 
colorectal cancer samples (`CRC`), so we first remove the other samples, 
which represent advanced adenoma (`ADA`) or non-advanced adenoma (`NAA`).
Also, we transform the metadata into a data.frame object (which is easier 
for some analyses later on).


```r
df.meta <- df.meta %>% 
  filter(Group %in% c('CRC', 'CTR')) %>% 
  as.data.frame()
rownames(df.meta) <- df.meta$Sample_ID

# restrict the taxonomic profiles to CRC and CTR samples
tax.profiles <- tax.profiles[,rownames(df.meta)]
rel.tax.profiles <- rel.tax.profiles[,rownames(df.meta)]
```

## Feature filtering

Currently, the matrix of taxonomic profiles contains `14213` different 
bacterial species. Of those, not all will be relevant for our question, since
some are present only in a handful of samples (low prevalence) or at extremely
low abundance. Therefore, it can make sense to filter your taxonomic profiles
before you begin the analysis. Here, we could for example use the maximum
species abundance as a filtering criterion. All species that have a relative
abundance of at least `1e-03` in at least one of the samples will be kept, 
the rest is filtered out. 



```r
species.max.value <- apply(rel.tax.profiles, 1, max)
f.idx <- which(species.max.value > 1e-03)
rel.tax.profiles.filt <- rel.tax.profiles[f.idx,]
```

Additionally, the mOTUs2 profiler can also estimate how much of the relative 
abundance cannot be classified. We also filter out this share of "Unknown".


```r
# unclassified are indicated by -1 in the mOTUs2 profiles
idx.um <- which(rownames(rel.tax.profiles.filt) == '-1')
rel.tax.profiles.filt <- rel.tax.profiles.filt[-idx.um,]
```


# Association Testing

Now that we have set up everyting, we can test all microbial species 
for statistically significant differences. In order to do so, we perform a 
Wilcoxon test on each individual bacterial species.


```r
p.vals <- rep_len(1, nrow(rel.tax.profiles.filt))
names(p.vals) <- rownames(rel.tax.profiles.filt)
stopifnot(all(rownames(df.meta) == colnames(rel.tax.profiles.filt)))
for (i in rownames(rel.tax.profiles.filt)){
  x <- rel.tax.profiles.filt[i,]
  y <- df.meta$Group
  t <- wilcox.test(x~y)
  p.vals[i] <- t$p.value
}
head(sort(p.vals))
## Fusobacterium nucleatum subsp. animalis [ref_mOTU_v25_01001] 
##                                                 1.072932e-07 
##                  Dialister pneumosintes [ref_mOTU_v25_03630] 
##                                                 5.047528e-06 
##                    Hungatella hathewayi [ref_mOTU_v25_03436] 
##                                                 3.079069e-05 
##           Olsenella sp. Marseille-P2300 [ref_mOTU_v25_10001] 
##                                                 2.283988e-04 
##                         Clostridium sp. [ref_mOTU_v25_03680] 
##                                                 2.756465e-04 
##                   [Eubacterium] eligens [ref_mOTU_v25_02375] 
##                                                 3.395879e-04
```

The species with the most significant effect seems to be _Fusobacterium
nucleatum_, so let us take a look at the distribution of this species:


```r
species <- 'Fusobacterium nucleatum subsp. animalis [ref_mOTU_v25_01001]'
df.plot <- tibble(fuso=rel.tax.profiles.filt[species,],
                  label=df.meta$Group)
df.plot %>% 
  ggplot(aes(x=label, y=fuso)) + 
    geom_boxplot() + 
    xlab('') + 
    ylab('F. nucleatum rel. ab.')
```

![](comparative_metagenome_analysis_files/figure-html/fuso-1.png)<!-- -->

Let us remember that log-scales are important when visualizing relative 
abundance data!


```r
df.plot %>% 
  ggplot(aes(x=label, y=log10(fuso + 1e-05))) + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.08) + 
    xlab('') + 
    ylab('log10(F. nucleatum rel. ab.)')
```

![](comparative_metagenome_analysis_files/figure-html/fuso_2-1.png)<!-- -->

# SIAMCAT Association Testing

We can also use the `SIAMCAT` R package to test for differential abundance and
produce standard visualizations.


```r
library("SIAMCAT")
```

Within `SIAMCAT`, the data are stored in the `SIAMCAT` object which contains 
the feature matrix, the metadata, and information about the groups you want to
compare.


```r
sc.obj <- siamcat(feat=rel.tax.profiles, meta=df.meta, 
                  label='Group', case='CRC')
## + starting create.label
## Label used as case:
##    CRC
## Label used as control:
##    CTR
## + finished create.label.from.metadata in 0.026 s
## + starting validate.data
## +++ checking overlap between labels and features
## + Keeping labels of 114 sample(s).
## +++ checking sample number per class
## +++ checking overlap between samples and metadata
## + finished validate.data in 0.444 s
```

We can use `SIAMCAT` for feature filtering as well:


```r
sc.obj <- filter.features(sc.obj, filter.method = 'abundance', cutoff = 1e-03)
## Features successfully filtered
sc.obj <- filter.features(sc.obj, filter.method = 'prevalence', 
                          cutoff = 0.05, feature.type = 'filtered')
## Features successfully filtered
sc.obj
## siamcat-class object
## label()                Label object:         61 CTR and 53 CRC samples
## filt_feat()            Filtered features:    839 features after abundance, prevalence filtering
## 
## contains phyloseq-class experiment-level object @phyloseq:
## phyloseq@otu_table()   OTU Table:            [ 14213 taxa and 114 samples ]
## phyloseq@sam_data()    Sample Data:          [ 114 samples by 13 sample variables ]
```

Now, we can test the filtered feature for differential abundance with `SIAMCAT`:


```r
sc.obj <- check.associations(sc.obj, detect.lim = 1e-05, 
                             fn.plot = './associations.pdf')
```
![](comparative_metagenome_analysis_files/figure-html/sc_assoc_testing_real-1.png)<!-- -->

# Exercises: Visualization 


* The associations metrics computed by `SIAMCAT` are stored in the `SIAMCAT` 
object and can be extracted by using `associations(sc.obj)`, if you want to 
have a closer look at the results for yourself. Plot a volcano plot of the 
associations between cancer and controls using the output from `SIAMCAT`.

* **Optional** Create a ordination plot for our data and colour the samples 
by group. How would you interpret the results? Try out different ecological 
distances. How does the choice of distance affect the group separation?
(**Tip**: make sure to check out the `vegdist` function in the **vegan** 
package and also the `pco` function in the **labdsv** package)

# Machine learning with SIAMCAT

## Machine learning workflow

The `SIAMCAT` machine learning workflow consists of several steps:

<!--html_preserve--><div id="htmlwidget-33e389e6993a9f62baa2" style="width:672px;height:480px;" class="grViz html-widget"></div>
<script type="application/json" data-for="htmlwidget-33e389e6993a9f62baa2">{"x":{"diagram":"\ndigraph siamcat_workflow {\n\n  # a \"graph\" statement\n  graph [overlap = true, fontsize = 10]\n\n  # several \"node\" statements\n  node [shape = box,\n  fontname = Helvetica,\n  label=\"filter.features\"]\n  C\n\n  node [shape = box,\n  fontname = Helvetica,\n  label=\"normalize.features\"]\n  F\n\n  node [shape = box,\n  fontname = Helvetica,\n  label=\"create.data.split\"]\n  G\n\n  node [shape = box,\n  fontname = Helvetica,\n  label=\"train.model\"]\n  H\n\n  node [shape = box,\n  fontname = Helvetica,\n  label=\"make.predictions\"]\n  I\n\n  node [shape = box,\n  fontname = Helvetica,\n  label=\"evaluate.prediction\"]\n  J\n\n  node [shape = box,\n  fontname = Helvetica,\n  label=\"model.evaluation.plot\"]\n  K\n\n  node [shape = box,\n  fontname = Helvetica,\n  label=\"model.interpretation.plot\"]\n  L\n\n  # several \"edge\" statements\n  C->F\n  F->G G->H H->I I->J J->K\n  J->L\n  }\n","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->

Since we already created a `SIAMCAT` object and filtered the raw data, we can
start directly with the next step.

## Normalization

`SIAMCAT` offers a few normalization approaches that can be useful for
subsequent statistical modeling in the sense that they transform features in
a way that can increase the accuracy of the resulting models. Importantly,
these normalization techniques do not make use of any label information
(patient status), and can thus be applied up front to the whole data set 
(and outside of the following cross validation).


```r
sc.obj <- normalize.features(sc.obj, norm.method = 'log.std',
                             norm.param = list(log.n0=1e-05, sd.min.q=0))
## Features normalized successfully.
sc.obj
## siamcat-class object
## label()                Label object:         61 CTR and 53 CRC samples
## filt_feat()            Filtered features:    839 features after abundance, prevalence filtering
## associations()         Associations:         Results from association testing
##                                              with 11 significant features at alpha 0.05
## norm_feat()            Normalized features:  839 features normalized using log.std
## 
## contains phyloseq-class experiment-level object @phyloseq:
## phyloseq@otu_table()   OTU Table:            [ 14213 taxa and 114 samples ]
## phyloseq@sam_data()    Sample Data:          [ 114 samples by 13 sample variables ]
```

## Cross validation split

Cross validation is a technique to assess how well an ML model would generalize 
to external data by partionining the dataset into training and test sets.
Here, we split the dataset into 10 parts and then train a model on 9 of these
parts and use the left-out part to test the model. The whole process is 
repeated 10 times.


```r
sc.obj <- create.data.split(sc.obj, num.folds = 10, num.resample = 10)
## Features splitted for cross-validation successfully.
```


## Model Training

Now, we can train a
[LASSO logistic regression classifier](https://www.jstor.org/stable/2346178)
in order to distinguish CRC cases and controls.


```r
sc.obj <- train.model(sc.obj, method='lasso')
## Trained lasso models successfully.
```

## Predictions

This function will automatically apply the models trained in cross validation 
to their respective test sets and aggregate the predictions across the whole 
data set.


```r
sc.obj <- make.predictions(sc.obj)
## Made predictions successfully.
```

## Model Evaluation 

Calling the `evaluate.predictions` funtion will result in an assessment of
precision and recall as well as in ROC analysis, both of which can be plotted
as a pdf file using the `model.evaluation.plot` funtion (the name of/path to
the pdf file is passed as an argument).


```r
sc.obj <- evaluate.predictions(sc.obj)
## Evaluated predictions successfully.
model.evaluation.plot(sc.obj, fn.plot = './eval_plot.pdf')
## Plotted evaluation of predictions successfully to: ./eval_plot.pdf
```

![](./eval_plot.png)

## Model interpretation

Finally, the `model.interpretation.plot` function will plot characteristics 
of the models (i.e. model coefficients or feature importance) alongside the 
input data aiding in understanding how / why the model works (or not).


```r
model.interpretation.plot(sc.obj, consens.thres = 0.7,
                          fn.plot = './interpretation_plot.pdf')
## Successfully plotted model interpretation plot to: ./interpretation_plot.pdf
```

![](./interpretation_plot.png)


# Exercise: Prediction on external data

On the EMBL cluster, there is another dataset from a colorectal cancer 
metagenomic study. The study population was recruited in Austria, so you can
find the data under:

```r
fn.feat.at <- paste0(data.location, '/motus_profiles/AT-CRC.motus')
fn.meta.at <- paste0(data.location, '/metadata/meta_AT-CRC.tsv')
```

* Apply the trained model on this dataset and check the model performance 
on the external dataset.

* Train a `SIAMCAT` model on the Austrian dataset and apply it to the French 
dataset. How does the model transfer on the external dataset compare between
the two datasets? Compare also the feature weights when training on the French
or Austrian dataset.  
**Note**: You can supply several SIAMCAT objects to the function 
`model.evaluation.plot` and compare two ROC curves in the same plot.

# Exercise: Taxonomic vs functional predictors

In addition to the taxonomic profiles, we also created functional profiles
for the French dataset. You can find it under:

```r
fn.feat.fr.kegg <- paste0(data.location, 
                          '/functional_profiles/KEGG_kos_FR-CRC.tsv')
```

* Explore the distribution of the functional data (abundance distribution, 
prevalence, etc.) and compare it to what you observe with the taxonomic 
profiles. Which filtering regime would make sense for functional data?

* Use `SIAMCAT` to train a model based on the functional KEGG profiles and 
compare it to the one trained on taxonomic profiles.  
**Note** Since the functional profiles will have many thousands features,
it makes sense to perform feature selection on your dataset. You can 
supply the parameters for the feature selection to the `train.model` function
in `SIAMCAT`.

# Further information

You can find more information about `SIAMCAT` on https://siamcat.embl.de 
or on Bioconductor under 
https://www.bioconductor.org/packages/release/bioc/html/SIAMCAT.html

There you can also find several vignettes which go into more detail about 
different applications for `SIAMCAT`.

# SessionInfo


```r
sessionInfo()
## R version 4.0.2 (2020-06-22)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Catalina 10.15.5
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] DiagrammeR_1.0.6.1 SIAMCAT_1.11.0     phyloseq_1.32.0    mlr_2.18.0        
##  [5] ParamHelpers_1.14  forcats_0.5.0      stringr_1.4.0      dplyr_1.0.2       
##  [9] purrr_0.3.4        readr_1.4.0        tidyr_1.1.2        tibble_3.0.4      
## [13] ggplot2_3.3.2      tidyverse_1.3.0   
## 
## loaded via a namespace (and not attached):
##   [1] colorspace_1.4-1    ellipsis_0.3.1      XVector_0.28.0     
##   [4] fs_1.5.0            rstudioapi_0.11     farver_2.0.3       
##   [7] fansi_0.4.1         lubridate_1.7.9     xml2_1.3.2         
##  [10] codetools_0.2-18    splines_4.0.2       PRROC_1.3.1        
##  [13] knitr_1.30          ade4_1.7-16         jsonlite_1.7.1     
##  [16] pROC_1.16.2         broom_0.7.2         gridBase_0.4-7     
##  [19] cluster_2.1.0       dbplyr_2.0.0        compiler_4.0.2     
##  [22] httr_1.4.2          backports_1.2.0     assertthat_0.2.1   
##  [25] Matrix_1.2-18       cli_2.1.0           visNetwork_2.0.9   
##  [28] htmltools_0.5.0     prettyunits_1.1.1   tools_4.0.2        
##  [31] igraph_1.2.6        gtable_0.3.0        glue_1.4.2         
##  [34] reshape2_1.4.4      LiblineaR_2.10-8    fastmatch_1.1-0    
##  [37] Rcpp_1.0.5          parallelMap_1.5.0   Biobase_2.48.0     
##  [40] cellranger_1.1.0    vctrs_0.3.4         Biostrings_2.56.0  
##  [43] multtest_2.44.0     ape_5.4-1           nlme_3.1-150       
##  [46] iterators_1.0.13    xfun_0.19           rvest_0.3.6        
##  [49] lifecycle_0.2.0     XML_3.99-0.5        beanplot_1.2       
##  [52] zlibbioc_1.34.0     MASS_7.3-53         scales_1.1.1       
##  [55] hms_0.5.3           parallel_4.0.2      biomformat_1.16.0  
##  [58] rhdf5_2.32.4        RColorBrewer_1.1-2  BBmisc_1.11        
##  [61] curl_4.3            yaml_2.2.1          gridExtra_2.3      
##  [64] stringi_1.5.3       S4Vectors_0.26.1    corrplot_0.84      
##  [67] foreach_1.5.1       checkmate_2.0.0     permute_0.9-5      
##  [70] BiocGenerics_0.34.0 shape_1.4.5         rlang_0.4.8        
##  [73] pkgconfig_2.0.3     matrixStats_0.57.0  evaluate_0.14      
##  [76] lattice_0.20-41     Rhdf5lib_1.10.1     htmlwidgets_1.5.2  
##  [79] labeling_0.4.2      tidyselect_1.1.0    plyr_1.8.6         
##  [82] magrittr_1.5        R6_2.5.0            IRanges_2.22.2     
##  [85] generics_0.1.0      DBI_1.1.0           pillar_1.4.6       
##  [88] haven_2.3.1         withr_2.3.0         mgcv_1.8-33        
##  [91] survival_3.2-7      modelr_0.1.8        crayon_1.3.4       
##  [94] rmarkdown_2.5       progress_1.2.2      grid_4.0.2         
##  [97] readxl_1.3.1        data.table_1.13.2   vegan_2.5-6        
## [100] infotheo_1.2.0      reprex_0.3.0        digest_0.6.27      
## [103] stats4_4.0.2        munsell_0.5.0       glmnet_4.0-2
```
