
#####PACKAGES##################################################################

library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library("tidyverse")
library("MLmetrics")
library("ROCR")
library("mlmRev")
library("lme4")
library("rstanarm")
library("DescTools")
library("cmdstanr")
library("posterior")
library("rethinking")
library("AMR")
library("caret")
library("data.table")
library("devtools")
library("amrabxlookup")
library("corrplot")
library("bayesplot")
library("dbarts")
library("glue")
library("ggridges")
library("tidymodels")
library("pak")
library("devtools")
library("truncnorm")
library("missMethods")
library("dichromat")
library("RColorBrewer")
setwd("/Users/alexhoward/Documents/GitHub/MIMIC-IV")

#####FUNCTIONS#################################################################

#Read-in and cleaning
micro_clean <- function(file_location,microbiologyevents_filename) {
  
  path_to_data <- file_location
  
  read_in <- function(file_name) {
    
    file_path <- file.path(path_to_data, file_name)
    df <- fread(file_path)
    
  }
  
  df <- read_in("microbiologyevents3.csv")
  
  df <- df %>% mutate(org_name = str_remove(org_name,"POSITIVE FOR"),#---------------------------------------------------Cleaning
                      org_name = str_remove(org_name,"PRESUMPTIVELY"),
                      org_name = str_remove(org_name,"PRESUMPTIVE"),
                      org_name = str_remove(org_name,"PROBABLE"),
                      org_name = str_remove(org_name,"IDENTIFICATION"),
                      org_name = str_remove(org_name,"RESEMBLING"),
                      org_name = str_remove(org_name,"SEEN"),
                      org_name = str_remove(org_name,"MODERATE"),
                      org_name = str_remove(org_name,"FEW"),
                      org_name = str_remove(org_name,"BETA"),
                      org_name = str_remove(org_name,"METHICILLIN RESISTANT"),
                      org_name = str_remove(org_name,"NUTRITIONALLY VARIANT"),
                      org_name = str_remove(org_name,"NOT C. PERFRINGENS OR C. SEPTICUM"),
                      org_name = str_remove(org_name,"-LACTAMASE POSITIVE"),
                      org_name = str_remove(org_name,"-LACTAMASE NEGATIVE"),
                      org_name = str_remove(org_name,"VIRAL ANTIGEN"),
                      org_name = str_remove(org_name,"CANDIDA INCONSPICUA"),
                      org_name = str_remove(org_name,"/POSADASII"),
                      org_name = str_remove(org_name,"NOT FUMIGATUS, FLAVUS OR NIGER"),
                      org_name = str_remove(org_name,"MRSA POSITIVE"),
                      org_name = str_remove(org_name,"MRSA NEGATIVE"),
                      org_name = str_remove(org_name,"HISTOLYTICA/DISPAR"),
                      org_name = case_when(grepl("NON-FERMENTER",org_name)~"PSEUDOMONADALES",                     
                                           TRUE~org_name),
                      org_name = case_when(grepl("ABIOTROPHIA/GRANULICATELLA",org_name)~"STREPTOCOCCUS",                     
                                           TRUE~org_name),
                      org_name = case_when(grepl("S. AUREUS POSITIVE",org_name)~"STAPHYLOCOCCUS AUREUS",                     
                                           TRUE~org_name),
                      org_name = case_when(grepl("ASPERGILLUS FUMIGATUS COMPLEX",org_name)~"ASPERGILLUS FUMIGATUS",                     
                                           TRUE~org_name),
                      org_name = case_when(grepl("(CRYPTOSPORIDIUM PARVUM OOCYSTS|CUNNINGHAMELLA BERTHOLLETIAE|EPIDERMOPHYTON FLOCCOSUM|EXOPHIALA JEANSELMEI COMPLEX|SCEDOSPORIUM|NEOASCOCHYTA DESMAZIERI|NEOSCYTALIDIUM DIMIDIATUM|LOMENTOSPORA|NEUROSPORA|PERONEUTYPA SCOPARIA|SPOROTHRIX SCHENCKII COMPLEX|ZYGOSACCHAROMYCES FERMENTATI)",org_name)~"UNKNOWN FUNGUS",                     
                                           TRUE~org_name)
  ) %>%#-------------------------------------------------------------------------------------------------------------------Removal of AMR package non-interpretable rows
    filter(!grepl("(CANCELLED|VIRUS|SIMPLEX|PARAINFLUENZA|INFLUENZA A|INFLUENZA B|TICK|AFB GROWN|GRAM VARIABLE RODS|HYMENOLEPIS)",org_name)) %>% 
    mutate(ab_name=AMR::as.ab(ab_name)) %>%#-------------------------------------------------------------------------------AMR package parsing of antimicrobial and organism names
    mutate(org_name=AMR::as.mo(org_name)) %>% 
    amrabxlookup::transpose_microbioevents(
      key_columns = c('subject_id','micro_specimen_id','isolate_num','org_name','ab_itemid','test_name','test_seq'),#------Transpose AST results to columns
      required_columns =c('subject_id','chartdate',"hadm_id","order_provider_id",
                          "charttime","micro_specimen_id","spec_itemid","spec_type_desc",
                          "storedate","storetime","test_itemid","test_name","org_itemid",
                          "isolate_num","org_name","comments",'test_seq'),
      transpose_key_column = 'ab_name',
      transpose_value_column = 'interpretation',
      fill = "N/A",
      non_empty_filter_column='subject_id') %>%
    add_column(AMX=NA, AMC=NA, TIC=NA,PME=NA, FOS=NA,TMP=NA,#---------------------------------------------------------------Add missing AMR package-recognised antimicrobial agent columns
               MFX=NA, NOR=NA,CPD = NA, FOX1=NA,TEC=NA,TLV=NA,ORI=NA,
               TGC=NA,AZM=NA,ATM=NA,CRB=NA,CTX=NA,CPT=NA,SPT=NA,TZD=NA,ERV=NA,OMC=NA,FDX=NA,
               CZT=NA,LEX=NA,CLR=NA,DAL=NA,CZA=NA,NOV=NA,ETP=NA,
               MTR=NA,QDA=NA,TEM=NA,COL=NA,CHL=NA,BPR=NA,CEC=NA) %>%
    mutate(org_fullname = AMR::mo_fullname(org_name),#----------------------------------------------------------------------Additional organism categorisation columns
           org_kingdom = AMR::mo_kingdom(org_name),
           org_phylum = AMR::mo_phylum(org_name),
           org_class = AMR::mo_class(org_name),
           org_order = AMR::mo_order(org_name),
           org_family = AMR::mo_family(org_name),
           org_genus = AMR::mo_genus(org_name),
           org_species = AMR::mo_species(org_name),
           org_gram = AMR::mo_gramstain(org_name),
           org_o2 = AMR::mo_oxygen_tolerance(org_name),
           org_path = AMR::mo_pathogenicity(org_name),
           UTI = case_when(grepl("URINE",spec_type_desc)~TRUE,
                           TRUE~FALSE)) %>%
    relocate(PEN,OXA,AMP,AMX,PIP,TIC,CRB,PME,SAM,AMC,TZP,TEM,#--------------------------------------------------------------AST column reorganisation
             ATM,
             LEX,CZO,CEC,CXM,FOX1,CTX,CRO,CAZ,CPD,FEP,CPT,BPR,CZA,CZT,
             ETP,MEM,IPM,
             LVX,MFX,CIP,NOR,
             GEN,TOB,SPT,
             TMP,SXT,
             COL,NIT,FOS,NOV,CHL,
             TGC,ERV,OMC,TCY,
             ERY,CLR,AZM,CLI,QDA,
             LNZ,TZD,TEC,VAN,DAL,TLV,ORI,DAP,RIF,
             FDX,MTR,
             .before = "AMK"
    ) %>% relocate(AMK,.after = "GEN") %>% 
    mutate_at(vars(PEN:MTR),as.sir)
  
  df %>% mutate(#----------------------------------------------------------------------------------------------------------Addition of breakpoint interpretation and UTI columns
    guideline=rep("Original CLSI",nrow(df)),
    urine_interp = case_when(spec_type_desc=="URINE" &
                               !is.na(org_name) &
                               (org_path=="Potentially pathogenic" |
                                  grepl("(discontinued|MIXED)",comments))~ "Possible UTI",
                             spec_type_desc=="URINE" &
                               !is.na(org_name) &
                               org_path=="Pathogenic" &
                               comments=="" ~ "Probable UTI",
                             TRUE ~ "Unlikely UTI"),
    AMPC=case_when(grepl("Citrobacter braakii",org_fullname) |#-----------------------------------------------------------Addition of chromosomal AmpC column
                     grepl("Citrobacter freundii",org_fullname) |
                     grepl("Citrobacter gillenii",org_fullname) |
                     grepl("Citrobacter murliniae",org_fullname) |
                     grepl("Citrobacter rodenticum",org_fullname) |
                     grepl("Citrobacter sedlakii",org_fullname) |
                     grepl("Citrobacter werkmanii",org_fullname) |
                     grepl("Citrobacter youngae",org_fullname) |
                     grepl("Enterobacter",org_fullname) |
                     grepl("Hafnia alvei",org_fullname) |
                     grepl("Klebsiella aerogenes",org_fullname) |
                     grepl("Morganella morganii",org_fullname) |
                     grepl("Providencia",org_fullname) |
                     grepl("Serratia",org_fullname) |
                     org_order=="Enterobacterales"& (CAZ=="R"|CAZ=="I")&FEP=="S"~"R",
                   (org_order=="Enterobacterales"& (CAZ=="R" & is.na(FEP))) |
                     (org_order=="Enterobacterales"& (is.na(CAZ) & FEP=="S" )) ~ NA,
                   TRUE~"S")
  ) %>% relocate(comments,.before="guideline")
  
}

#Intrinsic resistance population
intr_mic <- function(df) {
  
  x <- custom_eucast_rules(genus=="Enterococcus"~cephalosporins=="R",#----------------------------------------------Add custom rules
                           genus=="Enterococcus"~aminoglycosides=="R",
                           genus=="Enterococcus"~macrolides=="R",
                           genus=="Enterococcus"~lincosamides=="R",
                           fullname=="Enterococcus faecium"~carbapenems=="R",
                           genus=="Enterococcus"&AMP=="R"~SAM=="R",
                           genus=="Staphylococcus"&OXA=="S"~AMC=="S",
                           genus=="Staphylococcus"&OXA=="S"~SAM=="S",
                           genus=="Staphylococcus"&OXA=="S"~TZP=="S",
                           genus=="Staphylococcus"&OXA=="S"~AMC=="S",
                           genus=="Staphylococcus"&OXA=="S"~cephalosporins=="S",
                           genus=="Staphylococcus"&OXA=="S"~carbapenems=="S",
                           genus=="Staphylococcus"~CAZ=="R",
                           genus=="Staphylococcus"&OXA=="R"~AMC=="R",
                           genus=="Staphylococcus"&OXA=="R"~SAM=="R",
                           genus=="Staphylococcus"&OXA=="R"~TZP=="R",
                           genus=="Staphylococcus"&OXA=="R"~AMC=="R",
                           genus=="Staphylococcus"&OXA=="R"~cephalosporins_1st=="R",
                           genus=="Staphylococcus"&OXA=="R"~cephalosporins_2nd=="R",
                           genus=="Staphylococcus"&OXA=="R"~cephalosporins_3rd=="R",
                           genus=="Staphylococcus"&OXA=="R"~cephalosporins_4th=="R",
                           genus=="Staphylococcus"&OXA=="R"~carbapenems=="R",
                           genus=="Streptococcus"&PEN=="S"~aminopenicillins=="S",
                           genus=="Streptococcus"&PEN=="S"~ureidopenicillins=="S",
                           genus=="Streptococcus"&PEN=="S"~cephalosporins_except_caz=="S",
                           kingdom=="Fungi"~aminoglycosides=="R",
                           kingdom=="Fungi"~antimycobacterials=="R",
                           kingdom=="Fungi"~betalactams=="R",
                           kingdom=="Fungi"~quinolones=="R",
                           kingdom=="Fungi"~lincosamides=="R",
                           kingdom=="Fungi"~macrolides=="R",
                           kingdom=="Fungi"~oxazolidinones=="R",
                           kingdom=="Fungi"~polymyxins=="R",
                           kingdom=="Fungi"~streptogramins=="R",
                           kingdom=="Fungi"~tetracyclines=="R",
                           kingdom=="Fungi"~trimethoprims=="R",
                           kingdom=="Fungi"~glycopeptides=="R",
                           kingdom=="Fungi"~MTR=="R",
                           kingdom=="Fungi"~FDX=="R",
                           kingdom=="Fungi"~NIT=="R",
                           kingdom=="Fungi"~FOS=="R",
                           kingdom=="Fungi"~NOV=="R",
                           kingdom=="Fungi"~CHL=="R",
                           kingdom=="Fungi"~DAP=="R",
                           genus=="Enterobacter"~AMP=="R",
                           genus=="Enterobacter"~AMX=="R",
                           genus=="Enterobacter"~SAM=="R",
                           genus=="Enterobacter"~AMC=="R",
                           genus=="Enterobacter"~LEX=="R",
                           genus=="Enterobacter"~CZO=="R",
                           PIP=="S"~TZP=="S",
                           phylum=="Pseudomonadota"~DAP=="R",
                           org_o2=="aerobe"~MTR=="R")
  
  df %>% #---------------------------------------------------------------------------------------------------------Fill intrinsic resistance
    eucast_rules(col_mo = "org_name",
                 ampc_cephalosporin_resistance = "R",
                 rules="all",
                 custom_rules = x) %>% 
    mutate(MTR=case_when(org_o2!="anaerobe"~"R",
                         TRUE~MTR),
           cleaning = rep("(w/intr. R)",nrow(df)))
  
}

#Bayesian modelling and simulation
res_sim <- function(df,col,condition,col2,condition2,antibiotic,alpha_prior,beta_prior,antimicrobial_name,extra="") {
  
  antibiotic <- enquo(antibiotic)
  col <- enquo(col)
  col2 <- enquo(col2)
  
  df$isolate_id <- as.character(df$org_name)#------------------------------------------------------------------------Unique isolate id column
  df$isolate_id[!is.na(df$isolate_id)] <- 1:sum(!is.na(df$org_name))
  
  x <- nrow(df %>%#--------------------------------------------------------------------------------------------------Number of observed Rs
              dplyr::filter(grepl(condition,!!col) &
                              grepl(condition2,!!col2) &
                              !!antibiotic=="R"))
  
  N <- nrow(df %>%#--------------------------------------------------------------------------------------------------Total number of observed results
              dplyr::filter(grepl(condition,!!col) &
                              grepl(condition2,!!col2) &
                              !is.na(!!antibiotic)))
  
  if(N>1) {
    
    #Bayesian calculations
    p <- seq( from=0 , to=1 , length.out=1e4 )
    
    posterior_alpha <- alpha_prior + x
    
    posterior_beta <- beta_prior + N - x
    
    mean_prob <- posterior_alpha / (posterior_alpha + posterior_beta)
    mode_prob <- (posterior_alpha - 1) / (posterior_alpha + posterior_beta - 2)
    
    prior <- (p ^ (alpha_prior - 1)) * ((1 - p) ^ (beta_prior - 1))
    
    likelihood <- (p ^ x) * ( (1 - p) ^ (N - x)  )
    
    posterior <-  p ^ (posterior_alpha - 1) * (1 - p) ^ (posterior_beta - 1)
    
    if (mean(posterior) != 0 & mean(posterior) != 1) {
      
      #Sampling posterior distribution
      
      prior_samples <- sample( p , prob=prior , size=1e4 , replace=TRUE )
      prior_samples <- tibble(Probability = prior_samples,Distribution=rep("Prior",length(prior_samples)))
      
      likelihood_samples <- sample( p , prob=likelihood , size=1e4 , replace=TRUE )
      likelihood_samples <- tibble(Probability = likelihood_samples,Distribution=rep("Likelihood",length(likelihood_samples)))
      
      post_samples <- sample( p , prob=posterior , size=1e4 , replace=TRUE )
      post_samples <- tibble(Probability = post_samples,Distribution=rep("Posterior",length(post_samples)))
      
      #Prior, likelihood and posterior density plot
      
      prior_2 <- prior/max(prior)
      prior_plot <- tibble(Density = prior_2,Distribution=rep("Prior",length(prior_2)),Probability=p)
      
      likelihood_2 <- likelihood/max(likelihood)
      likelihood_plot <- tibble(Density = likelihood_2,Distribution=rep("Likelihood",length(likelihood_2)),Probability=p)
      
      posterior_2 <- posterior/max(posterior)
      post_plot <- tibble(Density = posterior_2,Distribution=rep("Posterior",length(posterior_2)),Probability=p)
      
      post_df <- rbind(prior_plot,likelihood_plot,post_plot)
      post_df$Distribution <- factor(post_df$Distribution, levels=c("Prior","Likelihood","Posterior"))
      
      print(ggplot(post_df,aes(y=Density,x=Probability,group=Distribution,fill=Distribution,color=Distribution)) +
        geom_line() +
        theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank()) +
        geom_line(size=.5) +
        geom_ribbon(data=subset(post_df,Probability>0 & Probability<1),aes(x=Probability,ymax=Density),ymin=0,alpha=0.3) +
        scale_fill_manual(name='', values=c("Prior" = "red", "Likelihood" = "green4","Posterior"="blue")) +
        guides(color = FALSE) +
        labs(title=glue("Probability: {antimicrobial_name} resistance in {condition}{extra}")))
      
      N_star <- nrow(df %>%
                       dplyr::filter(grepl(condition,!!col) &
                                       grepl(condition2,!!col2) &
                                       is.na(!!antibiotic)))
      
      #Posterior predictive bar plot
      
      BetaBinom <- Vectorize(function(x_star){
        log.val <- lchoose(N_star, x_star) + lbeta(posterior_alpha+x_star,posterior_beta+N_star-x_star) - lbeta(posterior_alpha,posterior_beta)
        return(exp(log.val))
      })
      
      
      post_predictive <- BetaBinom(1:N_star)
      plot(1:N_star,BetaBinom(1:N_star),type="h",col="darkblue",xlab="Estimated prevalence of resistance",ylab="Probability density",
           main = glue("Estimated prevalence of {antimicrobial_name} resistance in {N_star} {condition}{extra} isolates"),cex.axis= 1.5,cex.lab=1.5,lwd=4)
      
      samples <- sample( p , prob=posterior , size=1e4 , replace=TRUE )
      
    } else {
      
      samples = mean(posterior)
      
    }
    
    #Summary statistics
    
    simu <- rbetabinom(nrow(df %>%
                              dplyr::filter(grepl(condition,!!col) &
                                              grepl(condition2,!!col2) &
                                              is.na(!!antibiotic))),
                       size = 1 ,
                       shape1 = posterior_alpha,
                       shape2 = posterior_beta)
    
    n_likelihood <- nrow(df %>%
                           dplyr::filter(grepl(condition,!!col) &
                                           grepl(condition2,!!col2) &
                                           !is.na(!!antibiotic)))
    n_simulated <- nrow(df %>%
                          dplyr::filter(grepl(condition,!!col) &
                                          grepl(condition2,!!col2) &
                                          is.na(!!antibiotic)))
    
    
    if(mean(likelihood)!=0) {
      
      likelihood_samples <- sample( p , prob=likelihood , size=1e4 , replace=TRUE )
      
    } else {
      
      likelihood_samples <- 0
      
    }
    
    
    
    prior_amr_rate <- alpha_prior/(alpha_prior+beta_prior)
    mean_likelihood <- mean(rbinom(1e4,
                                   size = 1 ,
                                   prob = likelihood_samples))
    mean_posterior <- posterior_alpha / (posterior_alpha + posterior_beta)
    mode_posterior <- (posterior_alpha - 1) / (posterior_alpha + posterior_beta - 2)
    HPDI_posterior <- ifelse(mean(samples)==0,NA,data.frame(t(HPDI(samples,prob=0.95))))
    HPDI_samples <- data.frame(cbind(n_likelihood,n_simulated,prior_amr_rate,mean_likelihood,mean_posterior,mode_posterior,HPDI_posterior))
    rownames(HPDI_samples) <- glue("{condition}_{antimicrobial_name}{extra}")
    
    assign(glue("{condition}_{antimicrobial_name}{extra}"),HPDI_samples,envir = .GlobalEnv)
    
    #Result simulation
    
    simu <- ifelse(simu==1,"R","S")
    
    target <- df %>% dplyr::filter(grepl(condition,!!col) &
                                     grepl(condition2,!!col2)) %>% 
      select(isolate_id,!!antibiotic) %>% 
      mutate(!!antibiotic := as.character(!!antibiotic))
    
    target[is.na(target)] <- simu
    
    target <- target %>% distinct(isolate_id,.keep_all = T)
    
    df <- df %>% mutate(!!antibiotic := as.character(!!antibiotic)) %>% 
      rows_update(target,by=c("isolate_id"))
    
    sample_df <- tibble(cbind(tibble(samples=samples),tibble(org_name=rep(glue("{condition}"),length(samples)))))
    
    sample_df <-
      sample_df %>%
      group_by(org_name) %>%
      mutate(outlier = samples < quantile(samples, .25) - 1.5*IQR(samples) | samples > quantile(samples, .75) + 1.5*IQR(samples)) %>%
      ungroup
    
    assign(glue("{condition}_{antimicrobial_name}{extra}df"),sample_df,envir = .GlobalEnv)
    
    df %>% select(-isolate_id)
    
  } else {
    
    #Output for insufficient results to inform likelihood
    
    print(glue("Insufficient results to calculate {antimicrobial_name} resistance likelihood for {condition}"))
    
    missing <- data.frame(Antimicrobial=glue("{antimicrobial_name}"),Organism=glue("{condition}"))
    
    assign(glue("missing"),missing,envir = .GlobalEnv)
    
    missings <- rbind(missings,missing)
    
    assign(glue("missings"),missings,envir = .GlobalEnv)
    
    samples <- rep(1,1e4)
    
    HPDI_samples <- data.frame(matrix(nrow=1,ncol=7))
    rownames(HPDI_samples) <- glue("{condition}_{antimicrobial_name}{extra}")
    colnames(HPDI_samples) <- c("n_likelihood","n_simulated","prior_amr_rate",
                                "mean_likelihood","mean_posterior","mode_posterior",
                                "HPDI_posterior")
    
    assign(glue("{condition}_{antimicrobial_name}{extra}"),HPDI_samples,envir = .GlobalEnv)
    
    sample_df <- tibble(cbind(tibble(samples=rep(5,length(samples)))),tibble(org_name=rep(glue("{condition}"),length(samples))),
                        tibble(outlier=rep(FALSE,length(samples))))
    
    assign(glue("{condition}_{antimicrobial_name}{extra}_df"),sample_df,envir = .GlobalEnv)
    
    
    df %>% select(-isolate_id)
    
  }
  
}

#Observed prevalence simulation
norm_sim <- function(df,col,condition,col2,condition2,antibiotic,alpha_prior,beta_prior,antimicrobial_name,extra="") {
  
  antibiotic <- enquo(antibiotic)
  col <- enquo(col)
  col2 <- enquo(col2)
  
  df$isolate_id <- as.character(df$org_name)
  df$isolate_id[!is.na(df$isolate_id)] <- 1:sum(!is.na(df$org_name))
  
  x <- nrow(df %>%
              dplyr::filter(grepl(condition,!!col) &
                              grepl(condition2,!!col2) &
                              !!antibiotic=="R"))
  
  N <- nrow(df %>%
              dplyr::filter(grepl(condition,!!col) &
                              grepl(condition2,!!col2) &
                              !is.na(!!antibiotic)))
  
  if(N>1) {
    
    AMR_rate <- x/N
    
    
    simu <- rbinom(nrow(df %>%
                          dplyr::filter(grepl(condition,!!col) &
                                          grepl(condition2,!!col2) &
                                          is.na(!!antibiotic))),
                   size = 1,prob = AMR_rate)
    
    
    n_measured <- nrow(df %>%
                         dplyr::filter(grepl(condition,!!col) &
                                         grepl(condition2,!!col2) &
                                         !is.na(!!antibiotic)))
    n_simulated <- nrow(df %>%
                          dplyr::filter(grepl(condition,!!col) &
                                          grepl(condition2,!!col2) &
                                          is.na(!!antibiotic)))
    
    CIs <- prop.test(x,N)
    CIs <- CIs$conf.int
    CIs <- data.frame(t(CIs))
    colnames(CIs) <- c("lower_95_ci","upper_95_ci")
    
    simu <- ifelse(simu==1,"R","S")
    
    HPDI_AMR_rate <- data.frame(cbind(n_measured,n_simulated,AMR_rate,CIs))
    rownames(HPDI_AMR_rate) <- glue("{condition}_{antimicrobial_name}{extra}")
    
    assign(glue("{condition}_{antimicrobial_name}{extra}"),HPDI_AMR_rate,envir = .GlobalEnv)
    
    target <- df %>% dplyr::filter(grepl(condition,!!col) &
                                     grepl(condition2,!!col2)) %>% 
      select(isolate_id,!!antibiotic) %>% 
      mutate(!!antibiotic := as.character(!!antibiotic))
    
    target[is.na(target)] <- simu
    
    target <- target %>% distinct(isolate_id,.keep_all = T)
    
    df <- df %>% mutate(!!antibiotic := as.character(!!antibiotic)) %>% 
      rows_update(target,by=c("isolate_id"))
    
    sample_df <- tibble(cbind(tibble(AMR_rate=AMR_rate),tibble(org_name=rep(glue("{condition}"),length(AMR_rate)))),CIs)
    
    assign(glue("{condition}_{antimicrobial_name}{extra}df"),sample_df,envir = .GlobalEnv)
    
    df %>% select(-isolate_id)
    
  } else {
    
    print(glue("Insufficient results to calculate {antimicrobial_name} resistance likelihood for {condition}"))
    
    missing <- data.frame(Antimicrobial=glue("{antimicrobial_name}"),Organism=glue("{condition}"))
    
    assign(glue("missing"),missing,envir = .GlobalEnv)
    
    missings <- rbind(missings,missing)
    
    assign(glue("missings"),missings,envir = .GlobalEnv)
    
    AMR_rate <- 1
    
    HPDI_AMR_rate <- data.frame(matrix(nrow=1,ncol=5))
    rownames(HPDI_AMR_rate) <- glue("{condition}_{antimicrobial_name}{extra}")
    colnames(HPDI_AMR_rate) <- c("n_measured","n_simulated","AMR_rate","lower_95_ci","upper_95_ci")
    
    assign(glue("{condition}_{antimicrobial_name}{extra}_norms_"),AMR_rate,envir = .GlobalEnv)
    
    sample_df <- tibble(cbind(tibble(AMR_rate=rep(0,length(AMR_rate)))),tibble(org_name=rep(glue("{condition}"),length(AMR_rate))),
                        tibble(lower_95_ci=rep(0,length(AMR_rate))),tibble(upper_95_ci=rep(1,length(AMR_rate))))
    
    assign(glue("{condition}_{antimicrobial_name}{extra}_norms_df"),sample_df,envir = .GlobalEnv)
    
    df %>% select(-isolate_id)
    
  }
  
}

#Individual antimicrobial simulations
PEN_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Streptococcus pneumoniae",org_fullname,"",PEN,1,3,"Benzylpenicillin",) %>%
    res_sim(org_fullname,"Staphylococcus aureus",org_fullname,"",PEN,3,1,"Benzylpenicillin",)
  
  PEN_summary <- data.frame(rbind(
    `Streptococcus pneumoniae_Benzylpenicillin`,
    `Staphylococcus aureus_Benzylpenicillin`
  ))
  
  PENdf <- data.frame(rbind(
    `Streptococcus pneumoniae_Benzylpenicillindf`,
    `Staphylococcus aureus_Benzylpenicillindf`
  ))
  
  PENdf$org_name <- factor(PENdf$org_name, levels = PENdf %>% 
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  PEN_plot <- ggplot(PENdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = PENdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Benzylpenicillin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    scale_fill_brewer(palette = "Spectral") +
    theme_classic() +
    theme(legend.position = "none")
  
  print(PEN_plot)
  
  assign("PEN_summary",PEN_summary,envir = .GlobalEnv)
  
  df
  
}
PEN_normal <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    norm_sim(org_fullname,"Streptococcus pneumoniae",org_fullname,"",PEN,1,1,"Benzylpenicillin","_norms_") %>%
    norm_sim(org_fullname,"Staphylococcus aureus",org_fullname,"",PEN,1,1,"Benzylpenicillin","_norms_")
  
  PEN_norms_ <- data.frame(rbind(
    `Streptococcus pneumoniae_Benzylpenicillin_norms_`,
    `Staphylococcus aureus_Benzylpenicillin_norms_`))
  
  PEN_norms_df <- data.frame(rbind(
    `Streptococcus pneumoniae_Benzylpenicillin_norms_df`,
    `Staphylococcus aureus_Benzylpenicillin_norms_df`
  ))
  
  PEN_norms_df$org_name <- factor(PEN_norms_df$org_name, levels = PEN_norms_df %>% 
                                    distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  assign("PEN_norms_summary",PEN_norms_,envir = .GlobalEnv)
  
  
  df
  
}
AMP_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Streptococcus pneumoniae",org_fullname,"",AMP,1,3,"Ampicillin") %>% 
    res_sim(org_fullname,"^Enterococcus$",org_fullname,"",AMP,3,7,"Ampicillin") %>%
    res_sim(org_fullname,"Escherichia coli",org_fullname,"",AMP,3,3,"Ampicillin") %>%
    res_sim(org_fullname,"Proteus mirabilis",org_fullname,"",AMP,3,7,"Ampicillin")
  
  
  AMP_summary <- data.frame(rbind(
    `Streptococcus pneumoniae_Ampicillin`,
    `^Enterococcus$_Ampicillin`,
    `Escherichia coli_Ampicillin`,
    `Proteus mirabilis_Ampicillin`
  ))
  
  AMPdf <- data.frame(rbind(
    `Streptococcus pneumoniae_Ampicillindf`,
    `^Enterococcus$_Ampicillindf`,
    `Escherichia coli_Ampicillindf`,
    `Proteus mirabilis_Ampicillindf`
  ))
  
  AMPdf$org_name <- factor(AMPdf$org_name, levels = AMPdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  AMP_plot <- ggplot(AMPdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = AMPdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Ampicillin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(AMP_plot)
  
  assign("AMP_summary",AMP_summary,envir = .GlobalEnv)
  
  df
  
}
AMP_normal <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    norm_sim(org_fullname,"^Enterococcus$",org_fullname,"",AMP,1,1,"Ampicillin","_norms_") %>%
    norm_sim(org_fullname,"Escherichia coli",org_fullname,"",AMP,1,1,"Ampicillin","_norms_") %>%
    norm_sim(org_fullname,"Proteus mirabilis",org_fullname,"",AMP,1,1,"Ampicillin","_norms_") %>%
    norm_sim(org_fullname,"Streptococcus pneumoniae",org_fullname,"",AMP,1,1,"Ampicillin","_norms_")
  
  AMP_norms_ <- data.frame(rbind(
    `Streptococcus pneumoniae_Ampicillin_norms_`,
    `^Enterococcus$_Ampicillin_norms_`,
    `Escherichia coli_Ampicillin_norms_`,
    `Proteus mirabilis_Ampicillin_norms_`
  ))
  
  AMP_norms_df <- data.frame(rbind(
    `Streptococcus pneumoniae_Ampicillin_norms_df`,
    `^Enterococcus$_Ampicillin_norms_df`,
    `Escherichia coli_Ampicillin_norms_df`,
    `Proteus mirabilis_Ampicillin_norms_df`
  ))
  
  AMP_norms_df$org_name <- factor(AMP_norms_df$org_name, levels = AMP_norms_df %>%
                                    distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  assign("AMP_norms_summary",AMP_norms_,envir = .GlobalEnv)
  
  df
  
}
OXA_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  
  df <- df %>%
    res_sim(org_fullname,"Staphylococcus aureus",org_fullname,"",OXA,3,7,"Oxacillin",)
  
  OXA_summary <- data.frame(rbind(
    `Staphylococcus aureus_Oxacillin`
  ))
  
  OXAdf <- data.frame(rbind(
    `Staphylococcus aureus_Oxacillindf`
  ))
  
  OXAdf$org_name <- factor(OXAdf$org_name, levels = OXAdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  OXA_plot <- ggplot(OXAdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = OXAdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Oxacillin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(OXA_plot)
  
  assign("OXA_summary",OXA_summary,envir = .GlobalEnv)
  
  df
  
}
OXA_normal <- function(df) {
  
  par(mfrow = c(1,1))
  
  
  df <- df %>%
    norm_sim(org_fullname,"Staphylococcus aureus",org_fullname,"",OXA,1,1,"Oxacillin","_norms_")
  
  OXA_norms_summary <- data.frame(rbind(
    `Staphylococcus aureus_Oxacillin_norms_`
  ))
  
  OXA_norms_df <- data.frame(rbind(
    `Staphylococcus aureus_Oxacillin_norms_df`
  ))
  
  OXA_norms_df$org_name <- factor(OXA_norms_df$org_name, levels = OXA_norms_df %>%
                                    distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  assign("OXA_norms_summary",OXA_norms_summary,envir = .GlobalEnv)
  
  df
  
}
SAM_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Escherichia coli",org_fullname,"",SAM,3,7,"Ampicillin-sulbactam",) %>%
    res_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",SAM,3,7,"Ampicillin-sulbactam",) %>%
    res_sim(org_fullname,"Proteus mirabilis",org_fullname,"",SAM,3,7,"Ampicillin-sulbactam",)
  
  SAM_summary <- data.frame(rbind(
    `Escherichia coli_Ampicillin-sulbactam`,
    `Klebsiella pneumoniae_Ampicillin-sulbactam`,
    `Proteus mirabilis_Ampicillin-sulbactam`
  ))
  
  SAMdf <- data.frame(rbind(
    `Escherichia coli_Ampicillin-sulbactamdf`,
    `Klebsiella pneumoniae_Ampicillin-sulbactamdf`,
    `Proteus mirabilis_Ampicillin-sulbactamdf`
  ))
  
  SAMdf$org_name <- factor(SAMdf$org_name, levels = SAMdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  SAM_plot <- ggplot(SAMdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = SAMdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Ampicillin-sulbactam"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(SAM_plot)
  
  assign("SAM_summary",SAM_summary,envir = .GlobalEnv)
  
  df
  
}
SAM_normal <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    norm_sim(org_fullname,"Escherichia coli",org_fullname,"",SAM,1,1,"Ampicillin-sulbactam","_norms_") %>%
    norm_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",SAM,1,1,"Ampicillin-sulbactam","_norms_") %>%
    norm_sim(org_fullname,"Proteus mirabilis",org_fullname,"",SAM,1,1,"Ampicillin-sulbactam","_norms_")
  
  SAM_norms_summary <- data.frame(rbind(
    `Escherichia coli_Ampicillin-sulbactam_norms_`,
    `Klebsiella pneumoniae_Ampicillin-sulbactam_norms_`,
    `Proteus mirabilis_Ampicillin-sulbactam_norms_`
  ))
  
  SAM_norms_df <- data.frame(rbind(
    `Escherichia coli_Ampicillin-sulbactam_norms_df`,
    `Klebsiella pneumoniae_Ampicillin-sulbactam_norms_df`,
    `Proteus mirabilis_Ampicillin-sulbactam_norms_df`
  ))
  
  SAM_norms_df$org_name <- factor(SAM_norms_df$org_name, levels = SAM_norms_df %>%
                                    distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  assign("SAM_norms_summary",SAM_norms_summary,envir = .GlobalEnv)
  
  df
  
}
TZP_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Escherichia coli",org_fullname,"",TZP,1,3,"Piperacillin-tazobactam",) %>%
    res_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",TZP,1,3,"Piperacillin-tazobactam",) %>%
    res_sim(org_fullname,"Proteus mirabilis",org_fullname,"",TZP,1,3,"Piperacillin-tazobactam",) %>%
    res_sim(org_fullname,"Pseudomonas aeruginosa",org_fullname,"",TZP,1,3,"Piperacillin-tazobactam","")
  
  TZP_summary <- data.frame(rbind(
    `Escherichia coli_Piperacillin-tazobactam`,
    `Klebsiella pneumoniae_Piperacillin-tazobactam`,
    `Proteus mirabilis_Piperacillin-tazobactam`,
    `Pseudomonas aeruginosa_Piperacillin-tazobactam`
  ))
  
  TZPdf <- data.frame(rbind(
    `Escherichia coli_Piperacillin-tazobactamdf`,
    `Klebsiella pneumoniae_Piperacillin-tazobactamdf`,
    `Proteus mirabilis_Piperacillin-tazobactamdf`,
    `Pseudomonas aeruginosa_Piperacillin-tazobactamdf`
  ))
  
  TZPdf$org_name <- factor(TZPdf$org_name, levels = TZPdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  TZP_plot <- ggplot(TZPdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = TZPdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Piperacillin-tazobactam"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1)+
    theme_classic() +
    theme(legend.position = "none") 
  
  print(TZP_plot)
  
  assign("TZP_summary",TZP_summary,envir = .GlobalEnv)
  
  df
}
TZP_normal <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    norm_sim(org_fullname,"Escherichia coli",org_fullname,"",TZP,1,1,"Piperacillin-tazobactam","_norms_") %>%
    norm_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",TZP,1,1,"Piperacillin-tazobactam","_norms_") %>%
    norm_sim(org_fullname,"Proteus mirabilis",org_fullname,"",TZP,1,1,"Piperacillin-tazobactam","_norms_") %>%
    norm_sim(org_fullname,"Pseudomonas aeruginosa",org_fullname,"",TZP,1,1,"Piperacillin-tazobactam","_norms_")
  
  TZP_norms_summary <- data.frame(rbind(
    `Escherichia coli_Piperacillin-tazobactam_norms_`,
    `Klebsiella pneumoniae_Piperacillin-tazobactam_norms_`,
    `Proteus mirabilis_Piperacillin-tazobactam_norms_`,
    `Pseudomonas aeruginosa_Piperacillin-tazobactam_norms_`
  ))
  
  TZP_norms_df <- data.frame(rbind(
    `Escherichia coli_Piperacillin-tazobactam_norms_df`,
    `Klebsiella pneumoniae_Piperacillin-tazobactam_norms_df`,
    `Proteus mirabilis_Piperacillin-tazobactam_norms_df`,
    `Pseudomonas aeruginosa_Piperacillin-tazobactam_norms_df`
  ))
  
  TZP_norms_df$org_name <- factor(TZP_norms_df$org_name, levels = TZP_norms_df %>%
                                    distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  assign("TZP_norms_summary",TZP_norms_summary,envir = .GlobalEnv)
  
  df
}
CZO_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Escherichia coli",org_fullname,"",CZO,3,7,"Cefazolin") %>%
    res_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",CZO,3,7,"Cefazolin") %>%
    res_sim(org_fullname,"Proteus mirabilis",org_fullname,"",CZO,3,7,"Cefazolin")
  
  CZO_summary <- data.frame(rbind(
    `Escherichia coli_Cefazolin`,
    `Klebsiella pneumoniae_Cefazolin`,
    `Proteus mirabilis_Cefazolin`
  ))
  
  CZOdf <- data.frame(rbind(
    `Escherichia coli_Cefazolindf`,
    `Klebsiella pneumoniae_Cefazolindf`,
    `Proteus mirabilis_Cefazolindf`
  ))
  
  CZOdf$org_name <- factor(CZOdf$org_name, levels = CZOdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  CZO_plot <- ggplot(CZOdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = CZOdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Cefazolin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(CZO_plot)
  
  assign("CZO_summary",CZO_summary,envir = .GlobalEnv)
  
  df
}
CZO_normal <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    norm_sim(org_fullname,"Escherichia coli",org_fullname,"",CZO,1,1,"Cefazolin","_norms_") %>%
    norm_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",CZO,1,1,"Cefazolin","_norms_") %>%
    norm_sim(org_fullname,"Proteus mirabilis",org_fullname,"",CZO,1,1,"Cefazolin","_norms_")
  
  CZO_norms_summary <- data.frame(rbind(
    `Escherichia coli_Cefazolin_norms_`,
    `Klebsiella pneumoniae_Cefazolin_norms_`,
    `Proteus mirabilis_Cefazolin_norms_`
  ))
  
  CZO_norms_df <- data.frame(rbind(
    `Escherichia coli_Cefazolin_norms_df`,
    `Klebsiella pneumoniae_Cefazolin_norms_df`,
    `Proteus mirabilis_Cefazolin_norms_df`
  ))
  
  CZO_norms_df$org_name <- factor(CZO_norms_df$org_name, levels = CZO_norms_df %>%
                                    distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  assign("CZO_norms_summary",CZO_norms_summary,envir = .GlobalEnv)
  
  df
}#Resource dotplot

numbsims <- data.frame(cbind(numbsims,rep("BEAR",nrow(numbnorms))))
numbnorms <- data.frame(cbind(numbnorms,rep("EOP",nrow(numbnorms))))
numbsims[,2] <- as.numeric(numbsims[,2])
numbnorms[,2] <- as.numeric(numbnorms[,2])

colnames(numbsims) <- c("Antimicrobial","n","Estimator")
colnames(numbnorms) <- c("Antimicrobial","n","Estimator")
numbs <- data.frame(rbind(numbsims,numbnorms))

numbs$Antimicrobial <- factor(numbs$Antimicrobial, levels = numbs %>% 
                                arrange(desc(n)) %>%
                                distinct(Antimicrobial) %>% unlist()) %>% 
  fct_rev()

numbdiffs <- data.frame(numbsims[,2]-numbnorms[,2])
colnames(numbdiffs) <- c("n2")
numbdiffs$Antimicrobial <- numbsims$ab_name
numbdiffs <- numbdiffs %>% mutate(model = case_when(
  n2>0 ~ "BEAR",
  n2<0 ~ "EOP"
))

numbdiffs$n2 <- abs(numbdiffs$n2)
numbdiffs$Antimicrobial <- numbsims$Antimicrobial
numbdiffs$BEAR <- numbsims$n
numbdiffs$EOP <- numbnorms$n
numbdiffs <- numbdiffs %>% mutate(n = case_when(
  model=="BEAR"~BEAR+5,
  model=="EOP"~EOP+5
))

resource_plot <- ggplot(numbs, aes(x=Antimicrobial,y=n)) +
  geom_line(aes(group=Antimicrobial),alpha=0.5)+
  geom_point(aes(color=Estimator),size=4) +
  coord_flip() +
  ggtitle(glue("Minimum number of observed antimicrobial susceptibility results per pathogen required to avoid extreme estimation errors"))+
  xlab("Antimicrobial agent") +
  ylab("Number of observed results required") +
  ylim(0,110) +
  scale_color_manual(values=c("blue","green3"))+
  geom_text(data = numbdiffs, aes(color = model, 
                                  label = as.character(glue("+{n2}"))),
            size = 3,hjust=0.5)

print(resource_plot)
CXM_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Escherichia coli",org_fullname,"",CXM,1,3,"Cefuroxime") %>%
    res_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",CXM,1,3,"Cefuroxime") %>%
    res_sim(org_fullname,"Proteus mirabilis",org_fullname,"",CXM,1,3,"Cefuroxime")
  
  CXM_summary <- data.frame(rbind(
    `Escherichia coli_Cefuroxime`,
    `Klebsiella pneumoniae_Cefuroxime`,
    `Proteus mirabilis_Cefuroxime`
  ))
  
  CXMdf <- data.frame(rbind(
    `Escherichia coli_Cefuroximedf`,
    `Klebsiella pneumoniae_Cefuroximedf`,
    `Proteus mirabilis_Cefuroximedf`
  ))
  
  CXMdf$org_name <- factor(CXMdf$org_name, levels = CXMdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  CXM_plot <- ggplot(CXMdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = CXMdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Cefuroxime"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(CXM_plot)
  
  assign("CXM_summary",CXM_summary,envir = .GlobalEnv)
  
  df
}
CXM_normal <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    norm_sim(org_fullname,"Escherichia coli",org_fullname,"",CXM,1,1,"Cefuroxime","_norms_") %>%
    norm_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",CXM,1,1,"Cefuroxime","_norms_") %>%
    norm_sim(org_fullname,"Proteus mirabilis",org_fullname,"",CXM,1,1,"Cefuroxime","_norms_")
  
  CXM_norms_summary <- data.frame(rbind(
    `Escherichia coli_Cefuroxime_norms_`,
    `Klebsiella pneumoniae_Cefuroxime_norms_`,
    `Proteus mirabilis_Cefuroxime_norms_`
  ))
  
  CXM_norms_df <- data.frame(rbind(
    `Escherichia coli_Cefuroxime_norms_df`,
    `Klebsiella pneumoniae_Cefuroxime_norms_df`,
    `Proteus mirabilis_Cefuroxime_norms_df`
  ))
  
  CXM_norms_df$org_name <- factor(CXM_norms_df$org_name, levels = CXM_norms_df %>%
                                    distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  assign("CXM_norms_summary",CXM_norms_summary,envir = .GlobalEnv)
  
  df
}
CRO_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Escherichia coli",org_fullname,"",CRO,1,3,"Ceftriaxone") %>%
    res_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",CRO,1,3,"Ceftriaxone") %>% 
    res_sim(org_fullname,"Proteus mirabilis",org_fullname,"",CRO,1,3,"Ceftriaxone")
  
  CRO_summary <- data.frame(rbind(
    `Escherichia coli_Ceftriaxone`,
    `Klebsiella pneumoniae_Ceftriaxone`,
    `Proteus mirabilis_Ceftriaxone`
  ))
  
  CROdf <- data.frame(rbind(
    `Escherichia coli_Ceftriaxonedf`,
    `Klebsiella pneumoniae_Ceftriaxonedf`,
    `Proteus mirabilis_Ceftriaxonedf`
  ))
  
  CROdf$org_name <- factor(CROdf$org_name, levels = CROdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  CRO_plot <- ggplot(CROdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = CROdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Ceftriaxone"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(CRO_plot)
  
  assign("CRO_summary",CRO_summary,envir = .GlobalEnv)
  
  df
}
CRO_normal <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    norm_sim(org_fullname,"Escherichia coli",org_fullname,"",CRO,1,1,"Ceftriaxone","_norms_") %>%
    norm_sim(org_fullname,"Proteus mirabilis",org_fullname,"",CRO,1,1,"Ceftriaxone","_norms_") %>%
    norm_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",CRO,1,1,"Ceftriaxone","_norms_")
  
  CRO_norms_summary <- data.frame(rbind(
    `Escherichia coli_Ceftriaxone_norms_`,
    `Klebsiella pneumoniae_Ceftriaxone_norms_`,
    `Proteus mirabilis_Ceftriaxone_norms_`
  ))
  
  CRO_norms_df <- data.frame(rbind(
    `Escherichia coli_Ceftriaxone_norms_df`,
    `Klebsiella pneumoniae_Ceftriaxone_norms_df`,
    `Proteus mirabilis_Ceftriaxone_norms_df`
  ))
  
  CRO_norms_df$org_name <- factor(CRO_norms_df$org_name, levels = CRO_norms_df %>%
                                    distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  assign("CRO_norms_summary",CRO_norms_summary,envir = .GlobalEnv)
  
  df
}
CAZ_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Escherichia coli",org_fullname,"",CAZ,1,3,"Ceftazidime") %>%
    res_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",CAZ,1,3,"Ceftazidime") %>%
    res_sim(org_fullname,"Proteus mirabilis",org_fullname,"",CAZ,1,3,"Ceftazidime") %>%
    res_sim(org_fullname,"Pseudomonas aeruginosa",org_fullname,"",CAZ,1,3,"Ceftazidime")
  
  CAZ_summary <- data.frame(rbind(
    `Escherichia coli_Ceftazidime`,
    `Klebsiella pneumoniae_Ceftazidime`,
    `Proteus mirabilis_Ceftazidime`,
    `Pseudomonas aeruginosa_Ceftazidime`
  ))
  
  CAZdf <- data.frame(rbind(
    `Escherichia coli_Ceftazidimedf`,
    `Klebsiella pneumoniae_Ceftazidimedf`,
    `Proteus mirabilis_Ceftazidimedf`,
    `Pseudomonas aeruginosa_Ceftazidimedf`
  ))
  
  CAZdf$org_name <- factor(CAZdf$org_name, levels = CAZdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  CAZ_plot <- ggplot(CAZdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = CAZdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Ceftazidime"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(CAZ_plot)
  
  assign("CAZ_summary",CAZ_summary,envir = .GlobalEnv)
  
  df
}
CAZ_normal <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    norm_sim(org_fullname,"Escherichia coli",org_fullname,"",CAZ,1,1,"Ceftazidime","_norms_") %>%
    norm_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",CAZ,1,1,"Ceftazidime","_norms_") %>%
    norm_sim(org_fullname,"Proteus mirabilis",org_fullname,"",CAZ,1,1,"Ceftazidime","_norms_") %>%
    norm_sim(org_fullname,"Pseudomonas aeruginosa",org_fullname,"",CAZ,1,1,"Ceftazidime","_norms_")
  
  CAZ_norms_summary <- data.frame(rbind(
    `Escherichia coli_Ceftazidime_norms_`,
    `Klebsiella pneumoniae_Ceftazidime_norms_`,
    `Proteus mirabilis_Ceftazidime_norms_`,
    `Pseudomonas aeruginosa_Ceftazidime_norms_`
  ))
  
  CAZ_norms_df <- data.frame(rbind(
    `Escherichia coli_Ceftazidime_norms_df`,
    `Klebsiella pneumoniae_Ceftazidime_norms_df`,
    `Proteus mirabilis_Ceftazidime_norms_df`,
    `Pseudomonas aeruginosa_Ceftazidime_norms_df`
  ))
  
  CAZ_norms_df$org_name <- factor(CAZ_norms_df$org_name, levels = CAZ_norms_df %>%
                                    distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  assign("CAZ_norms_summary",CAZ_norms_summary,envir = .GlobalEnv)
  
  df
}
FEP_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Citrobacter freundii",org_fullname,"",FEP,1,3,"Cefepime") %>%
    res_sim(org_fullname,"Enterobacter cloacae",org_fullname,"",FEP,1,3,"Cefepime") %>%
    res_sim(org_fullname,"Escherichia coli",org_fullname,"",FEP,1,3,"Cefepime") %>%
    res_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",FEP,1,3,"Cefepime") %>%
    res_sim(org_fullname,"Morganella morganii",org_fullname,"",FEP,1,3,"Cefepime") %>%
    res_sim(org_fullname,"Proteus mirabilis",org_fullname,"",FEP,1,3,"Cefepime") %>%
    res_sim(org_fullname,"Pseudomonas aeruginosa",org_fullname,"",FEP,1,3,"Cefepime") %>%
    res_sim(org_fullname,"Serratia marcescens",org_fullname,"",FEP,1,3,"Cefepime")
  
  FEP_summary <- data.frame(rbind(
    `Citrobacter freundii_Cefepime`,
    `Enterobacter cloacae_Cefepime`,
    `Escherichia coli_Cefepime`,
    `Klebsiella pneumoniae_Cefepime`,
    `Morganella morganii_Cefepime`,
    `Proteus mirabilis_Cefepime`,
    `Pseudomonas aeruginosa_Cefepime`,
    `Serratia marcescens_Cefepime`
  ))
  
  FEPdf <- data.frame(rbind(
    `Citrobacter freundii_Cefepimedf`,
    `Enterobacter cloacae_Cefepimedf`,
    `Escherichia coli_Cefepimedf`,
    `Klebsiella pneumoniae_Cefepimedf`,
    `Morganella morganii_Cefepimedf`,
    `Proteus mirabilis_Cefepimedf`,
    `Pseudomonas aeruginosa_Cefepimedf`,
    `Serratia marcescens_Cefepimedf`
  ))
  
  FEPdf$org_name <- factor(FEPdf$org_name, levels = FEPdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  FEP_plot <- ggplot(FEPdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = FEPdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Cefepime"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(FEP_plot)
  
  assign("FEP_summary",FEP_summary,envir = .GlobalEnv)
  
  df
}
FEP_normal <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    norm_sim(org_fullname,"Citrobacter freundii",org_fullname,"",FEP,1,1,"Cefepime","_norms_") %>%
    norm_sim(org_fullname,"Enterobacter cloacae",org_fullname,"",FEP,1,1,"Cefepime","_norms_") %>%
    norm_sim(org_fullname,"Escherichia coli",org_fullname,"",FEP,1,1,"Cefepime","_norms_") %>%
    norm_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",FEP,1,1,"Cefepime","_norms_") %>%
    norm_sim(org_fullname,"Morganella morganii",org_fullname,"",FEP,1,1,"Cefepime","_norms_") %>%
    norm_sim(org_fullname,"Proteus mirabilis",org_fullname,"",FEP,1,1,"Cefepime","_norms_") %>%
    norm_sim(org_fullname,"Pseudomonas aeruginosa",org_fullname,"",FEP,1,1,"Cefepime","_norms_") %>%
    norm_sim(org_fullname,"Serratia marcescens",org_fullname,"",FEP,1,1,"Cefepime","_norms_")
  
  FEP_norms_summary <- data.frame(rbind(
    `Citrobacter freundii_Cefepime_norms_`,
    `Enterobacter cloacae_Cefepime_norms_`,
    `Escherichia coli_Cefepime_norms_`,
    `Klebsiella pneumoniae_Cefepime_norms_`,
    `Morganella morganii_Cefepime_norms_`,
    `Proteus mirabilis_Cefepime_norms_`,
    `Pseudomonas aeruginosa_Cefepime_norms_`,
    `Serratia marcescens_Cefepime_norms_`
  ))
  
  FEP_norms_df <- data.frame(rbind(
    `Citrobacter freundii_Cefepime_norms_df`,
    `Enterobacter cloacae_Cefepime_norms_df`,
    `Escherichia coli_Cefepime_norms_df`,
    `Klebsiella pneumoniae_Cefepime_norms_df`,
    `Morganella morganii_Cefepime_norms_df`,
    `Proteus mirabilis_Cefepime_norms_df`,
    `Pseudomonas aeruginosa_Cefepime_norms_df`,
    `Serratia marcescens_Cefepime_norms_df`
  ))
  
  FEP_norms_df$org_name <- factor(FEP_norms_df$org_name, levels = FEP_norms_df %>%
                                    distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  assign("FEP_norms_summary",FEP_norms_summary,envir = .GlobalEnv)
  
  df
}
MEM_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Citrobacter freundii",org_fullname,"",MEM,1,3,"Meropenem") %>%
    res_sim(org_fullname,"Enterobacter cloacae",org_fullname,"",MEM,1,3,"Meropenem") %>%
    res_sim(org_fullname,"Escherichia coli",org_fullname,"",MEM,1,3,"Meropenem") %>%
    res_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",MEM,1,3,"Meropenem") %>%
    res_sim(org_fullname,"Morganella morganii",org_fullname,"",MEM,1,3,"Meropenem") %>%
    res_sim(org_fullname,"Proteus mirabilis",org_fullname,"",MEM,1,3,"Meropenem") %>%
    res_sim(org_fullname,"Pseudomonas aeruginosa",org_fullname,"",MEM,1,3,"Meropenem") %>%
    res_sim(org_fullname,"Serratia marcescens",org_fullname,"",MEM,1,3,"Meropenem")
  
  MEM_summary <- data.frame(rbind(
    `Citrobacter freundii_Meropenem`,
    `Enterobacter cloacae_Meropenem`,
    `Escherichia coli_Meropenem`,
    `Klebsiella pneumoniae_Meropenem`,
    `Morganella morganii_Meropenem`,
    `Proteus mirabilis_Meropenem`,
    `Pseudomonas aeruginosa_Meropenem`,
    `Serratia marcescens_Meropenem`
  ))
  
  MEMdf <- data.frame(rbind(
    `Citrobacter freundii_Meropenemdf`,
    `Enterobacter cloacae_Meropenemdf`,
    `Escherichia coli_Meropenemdf`,
    `Klebsiella pneumoniae_Meropenemdf`,
    `Morganella morganii_Meropenemdf`,
    `Proteus mirabilis_Meropenemdf`,
    `Pseudomonas aeruginosa_Meropenemdf`,
    `Serratia marcescens_Meropenemdf`
  ))
  
  MEMdf$org_name <- factor(MEMdf$org_name, levels = MEMdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  MEM_plot <- ggplot(MEMdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = MEMdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Meropenem"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(MEM_plot)
  
  assign("MEM_summary",MEM_summary,envir = .GlobalEnv)
  
  df
}
MEM_normal <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    norm_sim(org_fullname,"Citrobacter freundii",org_fullname,"",MEM,1,1,"Meropenem","_norms_") %>%
    norm_sim(org_fullname,"Enterobacter cloacae",org_fullname,"",MEM,1,1,"Meropenem","_norms_") %>%
    norm_sim(org_fullname,"Escherichia coli",org_fullname,"",MEM,1,1,"Meropenem","_norms_") %>%
    norm_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",MEM,1,1,"Meropenem","_norms_") %>%
    norm_sim(org_fullname,"Morganella morganii",org_fullname,"",MEM,1,1,"Meropenem","_norms_") %>%
    norm_sim(org_fullname,"Proteus mirabilis",org_fullname,"",MEM,1,1,"Meropenem","_norms_") %>%
    norm_sim(org_fullname,"Pseudomonas aeruginosa",org_fullname,"",MEM,1,1,"Meropenem","_norms_") %>%
    norm_sim(org_fullname,"Serratia marcescens",org_fullname,"",MEM,1,1,"Meropenem","_norms_")
  
  MEM_norms_summary <- data.frame(rbind(
    `Citrobacter freundii_Meropenem_norms_`,
    `Enterobacter cloacae_Meropenem_norms_`,
    `Escherichia coli_Meropenem_norms_`,
    `Klebsiella pneumoniae_Meropenem_norms_`,
    `Morganella morganii_Meropenem_norms_`,
    `Proteus mirabilis_Meropenem_norms_`,
    `Pseudomonas aeruginosa_Meropenem_norms_`,
    `Serratia marcescens_Meropenem_norms_`
  ))
  
  MEM_norms_df <- data.frame(rbind(
    `Citrobacter freundii_Meropenem_norms_df`,
    `Enterobacter cloacae_Meropenem_norms_df`,
    `Escherichia coli_Meropenem_norms_df`,
    `Klebsiella pneumoniae_Meropenem_norms_df`,
    `Morganella morganii_Meropenem_norms_df`,
    `Proteus mirabilis_Meropenem_norms_df`,
    `Pseudomonas aeruginosa_Meropenem_norms_df`,
    `Serratia marcescens_Meropenem_norms_df`
  ))
  
  MEM_norms_df$org_name <- factor(MEM_norms_df$org_name, levels = MEM_norms_df %>%
                                    distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  assign("MEM_norms_summary",MEM_norms_summary,envir = .GlobalEnv)
  
  df
}
CIP_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Citrobacter freundii",org_fullname,"",CIP,3,7,"Ciprofloxacin") %>%
    res_sim(org_fullname,"Enterobacter cloacae",org_fullname,"",CIP,3,7,"Ciprofloxacin") %>%
    res_sim(org_fullname,"Escherichia coli",org_fullname,"",CIP,3,7,"Ciprofloxacin") %>%
    res_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",CIP,3,7,"Ciprofloxacin") %>%
    res_sim(org_fullname,"Morganella morganii",org_fullname,"",CIP,3,7,"Ciprofloxacin") %>%
    res_sim(org_fullname,"Proteus mirabilis",org_fullname,"",CIP,3,7,"Ciprofloxacin") %>%
    res_sim(org_fullname,"Pseudomonas aeruginosa",org_fullname,"",CIP,1,3,"Ciprofloxacin") %>%
    res_sim(org_fullname,"Serratia marcescens",org_fullname,"",CIP,3,7,"Ciprofloxacin")
  
  CIP_summary <- data.frame(rbind(
    `Citrobacter freundii_Ciprofloxacin`,
    `Enterobacter cloacae_Ciprofloxacin`,
    `Escherichia coli_Ciprofloxacin`,
    `Klebsiella pneumoniae_Ciprofloxacin`,
    `Morganella morganii_Ciprofloxacin`,
    `Proteus mirabilis_Ciprofloxacin`,
    `Pseudomonas aeruginosa_Ciprofloxacin`,
    `Serratia marcescens_Ciprofloxacin`
  ))
  
  CIPdf <- data.frame(rbind(
    `Citrobacter freundii_Ciprofloxacindf`,
    `Enterobacter cloacae_Ciprofloxacindf`,
    `Escherichia coli_Ciprofloxacindf`,
    `Klebsiella pneumoniae_Ciprofloxacindf`,
    `Morganella morganii_Ciprofloxacindf`,
    `Proteus mirabilis_Ciprofloxacindf`,
    `Pseudomonas aeruginosa_Ciprofloxacindf`,
    `Serratia marcescens_Ciprofloxacindf`
  ))
  
  CIPdf$org_name <- factor(CIPdf$org_name, levels = CIPdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  CIP_plot <- ggplot(CIPdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = CIPdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Ciprofloxacin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(CIP_plot)
  
  assign("CIP_summary",CIP_summary,envir = .GlobalEnv)
  
  df
}
CIP_normal <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    norm_sim(org_fullname,"Citrobacter freundii",org_fullname,"",CIP,1,1,"Ciprofloxacin","_norms_") %>%
    norm_sim(org_fullname,"Enterobacter cloacae",org_fullname,"",CIP,1,1,"Ciprofloxacin","_norms_") %>%
    norm_sim(org_fullname,"Escherichia coli",org_fullname,"",CIP,1,1,"Ciprofloxacin","_norms_") %>%
    norm_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",CIP,1,1,"Ciprofloxacin","_norms_") %>%
    norm_sim(org_fullname,"Morganella morganii",org_fullname,"",CIP,1,1,"Ciprofloxacin","_norms_") %>%
    norm_sim(org_fullname,"Proteus mirabilis",org_fullname,"",CIP,1,1,"Ciprofloxacin","_norms_") %>%
    norm_sim(org_fullname,"Pseudomonas aeruginosa",org_fullname,"",CIP,1,1,"Ciprofloxacin","_norms_") %>%
    norm_sim(org_fullname,"Serratia marcescens",org_fullname,"",CIP,1,1,"Ciprofloxacin","_norms_")
  
  CIP_norms_summary <- data.frame(rbind(
    `Citrobacter freundii_Ciprofloxacin_norms_`,
    `Enterobacter cloacae_Ciprofloxacin_norms_`,
    `Escherichia coli_Ciprofloxacin_norms_`,
    `Klebsiella pneumoniae_Ciprofloxacin_norms_`,
    `Morganella morganii_Ciprofloxacin_norms_`,
    `Proteus mirabilis_Ciprofloxacin_norms_`,
    `Pseudomonas aeruginosa_Ciprofloxacin_norms_`,
    `Serratia marcescens_Ciprofloxacin_norms_`
  ))
  
  CIP_norms_df <- data.frame(rbind(
    `Citrobacter freundii_Ciprofloxacin_norms_df`,
    `Enterobacter cloacae_Ciprofloxacin_norms_df`,
    `Escherichia coli_Ciprofloxacin_norms_df`,
    `Klebsiella pneumoniae_Ciprofloxacin_norms_df`,
    `Morganella morganii_Ciprofloxacin_norms_df`,
    `Proteus mirabilis_Ciprofloxacin_norms_df`,
    `Pseudomonas aeruginosa_Ciprofloxacin_norms_df`,
    `Serratia marcescens_Ciprofloxacin_norms_df`
  ))
  
  CIP_norms_df$org_name <- factor(CIP_norms_df$org_name, levels = CIP_norms_df %>%
                                    distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  assign("CIP_norms_summary",CIP_norms_summary,envir = .GlobalEnv)
  
  df
}
LVX_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Streptococcus pneumoniae",org_fullname,"",LVX,1,3,"Levofloxacin") %>% 
    res_sim(org_fullname,"Staphylococcus aureus",org_fullname,"",LVX,3,7,"Levofloxacin")
  
  
  LVX_summary <- data.frame(rbind(
    `Streptococcus pneumoniae_Levofloxacin`,
    `Staphylococcus aureus_Levofloxacin`
  ))
  
  LVXdf <- data.frame(rbind(
    `Streptococcus pneumoniae_Levofloxacindf`,
    `Staphylococcus aureus_Levofloxacindf`
  ))
  
  LVXdf$org_name <- factor(LVXdf$org_name, levels = LVXdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  LVX_plot <- ggplot(LVXdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = LVXdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Levofloxacin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(LVX_plot)
  
  assign("LVX_summary",LVX_summary,envir = .GlobalEnv)
  
  df
}
LVX_normal <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    norm_sim(org_fullname,"Staphylococcus aureus",org_fullname,"",LVX,1,1,"Levofloxacin","_norms_") %>%
    norm_sim(org_fullname,"Streptococcus pneumoniae",org_fullname,"",LVX,1,1,"Levofloxacin","_norms_")
  
  LVX_norms_summary <- data.frame(rbind(
    `Streptococcus pneumoniae_Levofloxacin_norms_`,
    `Staphylococcus aureus_Levofloxacin_norms_`
  ))
  
  LVX_norms_df <- data.frame(rbind(
    `Streptococcus pneumoniae_Levofloxacin_norms_df`,
    `Staphylococcus aureus_Levofloxacin_norms_df`
  ))
  
  LVX_norms_df$org_name <- factor(LVX_norms_df$org_name, levels = LVX_norms_df %>%
                                    distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  assign("LVX_norms_summary",LVX_norms_summary,envir = .GlobalEnv)
  
  df
}
ERY_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Streptococcus pneumoniae",org_fullname,"",ERY,3,7,"Erythromycin") %>% 
    res_sim(org_fullname,"Streptococcus Group B",org_fullname,"",ERY,3,3,"Erythromycin") %>%
    res_sim(org_fullname,"Staphylococcus aureus",org_fullname,"",ERY,3,3,"Erythromycin")
  
  
  ERY_summary <- data.frame(rbind(
    `Streptococcus pneumoniae_Erythromycin`,
    `Streptococcus Group B_Erythromycin`,
    `Staphylococcus aureus_Erythromycin`
  ))
  
  ERYdf <- data.frame(rbind(
    `Streptococcus pneumoniae_Erythromycindf`,
    `Streptococcus Group B_Erythromycindf`,
    `Staphylococcus aureus_Erythromycindf`
  ))
  
  ERYdf$org_name <- factor(ERYdf$org_name, levels = ERYdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  ERY_plot <- ggplot(ERYdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = ERYdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Erythromycin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(ERY_plot)
  
  assign("ERY_summary",ERY_summary,envir = .GlobalEnv)
  
  df
}
ERY_normal <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    norm_sim(org_fullname,"Streptococcus pneumoniae",org_fullname,"",ERY,1,1,"Erythromycin","_norms_") %>% 
    norm_sim(org_fullname,"Streptococcus Group B",org_fullname,"",ERY,1,1,"Erythromycin","_norms_") %>% 
    norm_sim(org_fullname,"Staphylococcus aureus",org_fullname,"",ERY,1,1,"Erythromycin","_norms_")
  
  ERY_norms_summary <- data.frame(rbind(
    `Streptococcus pneumoniae_Erythromycin_norms_`,
    `Streptococcus Group B_Erythromycin_norms_`,
    `Staphylococcus aureus_Erythromycin_norms_`
  ))
  
  ERY_norms_df <- data.frame(rbind(
    `Streptococcus pneumoniae_Erythromycin_norms_df`,
    `Streptococcus Group B_Erythromycin_norms_df`,
    `Staphylococcus aureus_Erythromycin_norms_df`
  ))
  
  ERY_norms_df$org_name <- factor(ERY_norms_df$org_name, levels = ERY_norms_df %>%
                                    distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  assign("ERY_norms_summary",ERY_norms_summary,envir = .GlobalEnv)
  
  df
}
CLI_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Streptococcus Group B",org_fullname,"",CLI,3,7,"Clindamycin") %>%
    res_sim(org_fullname,"Staphylococcus aureus",org_fullname,"",CLI,3,7,"Clindamycin")
  
  CLI_summary <- data.frame(rbind(
    `Streptococcus Group B_Clindamycin`,
    `Staphylococcus aureus_Clindamycin`
  ))
  
  CLIdf <- data.frame(rbind(
    `Streptococcus Group B_Clindamycindf`,
    `Staphylococcus aureus_Clindamycindf`
  ))
  
  CLIdf$org_name <- factor(CLIdf$org_name, levels = CLIdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  CLI_plot <- ggplot(CLIdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = CLIdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Clindamycin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(CLI_plot)
  
  assign("CLI_summary",CLI_summary,envir = .GlobalEnv)
  
  df
}
CLI_normal <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    norm_sim(org_fullname,"Staphylococcus aureus",org_fullname,"",CLI,1,1,"Clindamycin","_norms_") %>%
    norm_sim(org_fullname,"Streptococcus Group B",org_fullname,"",CLI,1,1,"Clindamycin","_norms_")
  
  CLI_norms_summary <- data.frame(rbind(
    `Streptococcus Group B_Clindamycin_norms_`,
    `Staphylococcus aureus_Clindamycin_norms_`
  ))
  
  CLI_norms_df <- data.frame(rbind(
    `Streptococcus Group B_Clindamycin_norms_df`,
    `Staphylococcus aureus_Clindamycin_norms_df`
  ))
  
  CLI_norms_df$org_name <- factor(CLI_norms_df$org_name, levels = CLI_norms_df %>%
                                    distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  assign("CLI_norms_summary",CLI_norms_summary,envir = .GlobalEnv)
  
  df
}
TCY_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Streptococcus pneumoniae",org_fullname,"",TCY,3,7,"Tetracycline") %>% 
    res_sim(org_fullname,"Staphylococcus aureus",org_fullname,"",TCY,1,3,"Tetracycline")
  
  
  TCY_summary <- data.frame(rbind(
    `Streptococcus pneumoniae_Tetracycline`,
    `Staphylococcus aureus_Tetracycline`
  ))
  
  TCYdf <- data.frame(rbind(
    `Streptococcus pneumoniae_Tetracyclinedf`,
    `Staphylococcus aureus_Tetracyclinedf`
  ))
  
  TCYdf$org_name <- factor(TCYdf$org_name, levels = TCYdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  TCY_plot <- ggplot(TCYdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = TCYdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Tetracycline"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(TCY_plot)
  
  assign("TCY_summary",TCY_summary,envir = .GlobalEnv)
  
  df
}
TCY_normal <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    norm_sim(org_fullname,"Staphylococcus aureus",org_fullname,"",TCY,1,1,"Tetracycline","_norms_") %>%
    norm_sim(org_fullname,"Streptococcus pneumoniae",org_fullname,"",TCY,1,1,"Tetracycline","_norms_")
  
  
  TCY_norms_summary <- data.frame(rbind(
    `Streptococcus pneumoniae_Tetracycline_norms_`,
    `Staphylococcus aureus_Tetracycline_norms_`
  ))
  
  TCY_norms_df <- data.frame(rbind(
    `Streptococcus pneumoniae_Tetracycline_norms_df`,
    `Staphylococcus aureus_Tetracycline_norms_df`
  ))
  
  TCY_norms_df$org_name <- factor(TCY_norms_df$org_name, levels = TCY_norms_df %>%
                                    distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  assign("TCY_norms_summary",TCY_norms_summary,envir = .GlobalEnv)
  
  df
}
VAN_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"^Enterococcus$",org_fullname,"",VAN,3,7,"Vancomycin")
  
  VAN_summary <- data.frame(rbind(
    `^Enterococcus$_Vancomycin`
  ))
  
  VANdf <- data.frame(rbind(
    `^Enterococcus$_Vancomycindf`
  ))
  
  VANdf$org_name <- factor(VANdf$org_name, levels = VANdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  VAN_plot <- ggplot(VANdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = VANdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Vancomycin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(VAN_plot)
  
  assign("VAN_summary",VAN_summary,envir = .GlobalEnv)
  
  df
}
VAN_normal <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    norm_sim(org_fullname,"^Enterococcus$",org_fullname,"",VAN,1,1,"Vancomycin","_norms_")
  
  
  VAN_norms_summary <- data.frame(rbind(
    `^Enterococcus$_Vancomycin_norms_`
  ))
  
  VAN_norms_df <- data.frame(rbind(
    `^Enterococcus$_Vancomycin_norms_df`
  ))
  
  VAN_norms_df$org_name <- factor(VAN_norms_df$org_name, levels = VAN_norms_df %>%
                                    distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  assign("VAN_norms_summary",VAN_norms_summary,envir = .GlobalEnv)
  
  df
}
RIF_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Staphylococcus aureus",org_fullname,"",RIF,1,3,"Rifampicin")
  
  RIF_summary <- data.frame(rbind(
    `Staphylococcus aureus_Rifampicin`
  ))
  
  RIFdf <- data.frame(rbind(
    `Staphylococcus aureus_Rifampicindf`
  ))
  
  RIFdf$org_name <- factor(RIFdf$org_name, levels = RIFdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  RIF_plot <- ggplot(RIFdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = RIFdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Rifampicin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(RIF_plot)
  
  assign("RIF_summary",RIF_summary,envir = .GlobalEnv)
  
  df
}
RIF_normal <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    norm_sim(org_fullname,"Staphylococcus aureus",org_fullname,"",RIF,1,1,"Rifampicin","_norms_")
  
  RIF_norms_summary <- data.frame(rbind(
    `Staphylococcus aureus_Rifampicin_norms_`
  ))
  
  RIF_norms_df <- data.frame(rbind(
    `Staphylococcus aureus_Rifampicin_norms_df`
  ))
  
  RIF_norms_df$org_name <- factor(RIF_norms_df$org_name, levels = RIF_norms_df %>%
                                    distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  assign("RIF_norms_summary",RIF_norms_summary,envir = .GlobalEnv)
  
  df
}
GEN_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Citrobacter freundii",org_fullname,"",GEN,1,3,"Gentamicin") %>%
    res_sim(org_fullname,"Enterobacter cloacae",org_fullname,"",GEN,1,3,"Gentamicin") %>%
    res_sim(org_fullname,"Escherichia coli",org_fullname,"",GEN,1,3,"Gentamicin") %>%
    res_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",GEN,1,3,"Gentamicin") %>%
    res_sim(org_fullname,"Morganella morganii",org_fullname,"",GEN,1,3,"Gentamicin") %>%
    res_sim(org_fullname,"Proteus mirabilis",org_fullname,"",GEN,1,3,"Gentamicin") %>%
    res_sim(org_fullname,"Pseudomonas aeruginosa",org_fullname,"",GEN,1,3,"Gentamicin") %>%
    res_sim(org_fullname,"Serratia marcescens",org_fullname,"",GEN,1,3,"Gentamicin") %>%
    res_sim(org_fullname,"Staphylococcus aureus",org_fullname,"",GEN,1,3,"Gentamicin")
  
  GEN_summary <- data.frame(rbind(
    `Citrobacter freundii_Gentamicin`,
    `Enterobacter cloacae_Gentamicin`,
    `Escherichia coli_Gentamicin`,
    `Klebsiella pneumoniae_Gentamicin`,
    `Morganella morganii_Gentamicin`,
    `Proteus mirabilis_Gentamicin`,
    `Pseudomonas aeruginosa_Gentamicin`,
    `Serratia marcescens_Gentamicin`,
    `Staphylococcus aureus_Gentamicin`
  ))
  
  GENdf <- data.frame(rbind(
    `Citrobacter freundii_Gentamicindf`,
    `Enterobacter cloacae_Gentamicindf`,
    `Escherichia coli_Gentamicindf`,
    `Klebsiella pneumoniae_Gentamicindf`,
    `Morganella morganii_Gentamicindf`,
    `Proteus mirabilis_Gentamicindf`,
    `Pseudomonas aeruginosa_Gentamicindf`,
    `Serratia marcescens_Gentamicindf`,
    `Staphylococcus aureus_Gentamicindf`
  ))
  
  GENdf$org_name <- factor(GENdf$org_name, levels = GENdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  GEN_plot <- ggplot(GENdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = GENdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Gentamicin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(GEN_plot)
  
  assign("GEN_summary",GEN_summary,envir = .GlobalEnv)
  
  df
}
GEN_normal <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    norm_sim(org_fullname,"Citrobacter freundii",org_fullname,"",GEN,1,1,"Gentamicin","_norms_") %>%
    norm_sim(org_fullname,"Enterobacter cloacae",org_fullname,"",GEN,1,1,"Gentamicin","_norms_") %>%
    norm_sim(org_fullname,"Escherichia coli",org_fullname,"",GEN,1,1,"Gentamicin","_norms_") %>%
    norm_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",GEN,1,1,"Gentamicin","_norms_") %>%
    norm_sim(org_fullname,"Morganella morganii",org_fullname,"",GEN,1,1,"Gentamicin","_norms_") %>%
    norm_sim(org_fullname,"Proteus mirabilis",org_fullname,"",GEN,1,1,"Gentamicin","_norms_") %>%
    norm_sim(org_fullname,"Pseudomonas aeruginosa",org_fullname,"",GEN,1,1,"Gentamicin","_norms_") %>%
    norm_sim(org_fullname,"Serratia marcescens",org_fullname,"",GEN,1,1,"Gentamicin","_norms_") %>%
    norm_sim(org_fullname,"Staphylococcus aureus",org_fullname,"",GEN,1,1,"Gentamicin","_norms_")
  
  GEN_norms_summary <- data.frame(rbind(
    `Citrobacter freundii_Gentamicin_norms_`,
    `Enterobacter cloacae_Gentamicin_norms_`,
    `Escherichia coli_Gentamicin_norms_`,
    `Klebsiella pneumoniae_Gentamicin_norms_`,
    `Morganella morganii_Gentamicin_norms_`,
    `Proteus mirabilis_Gentamicin_norms_`,
    `Pseudomonas aeruginosa_Gentamicin_norms_`,
    `Serratia marcescens_Gentamicin_norms_`,
    `Staphylococcus aureus_Gentamicin_norms_`
  ))
  
  GEN_norms_df <- data.frame(rbind(
    `Citrobacter freundii_Gentamicin_norms_df`,
    `Enterobacter cloacae_Gentamicin_norms_df`,
    `Escherichia coli_Gentamicin_norms_df`,
    `Klebsiella pneumoniae_Gentamicin_norms_df`,
    `Morganella morganii_Gentamicin_norms_df`,
    `Proteus mirabilis_Gentamicin_norms_df`,
    `Pseudomonas aeruginosa_Gentamicin_norms_df`,
    `Serratia marcescens_Gentamicin_norms_df`,
    `Staphylococcus aureus_Gentamicin_norms_df`
  ))
  
  GEN_norms_df$org_name <- factor(GEN_norms_df$org_name, levels = GEN_norms_df %>%
                                    distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  assign("GEN_norms_summary",GEN_norms_summary,envir = .GlobalEnv)
  
  df
}
AMK_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Enterobacter cloacae",org_fullname,"",AMK,1,3,"Amikacin") %>%
    res_sim(org_fullname,"Escherichia coli",org_fullname,"",AMK,1,3,"Amikacin") %>%
    res_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",AMK,1,3,"Amikacin") %>%
    res_sim(org_fullname,"Proteus mirabilis",org_fullname,"",AMK,1,3,"Amikacin") %>%
    res_sim(org_fullname,"Pseudomonas aeruginosa",org_fullname,"",AMK,1,3,"Amikacin")
  
  AMK_summary <- data.frame(rbind(
    `Enterobacter cloacae_Amikacin`,
    `Escherichia coli_Amikacin`,
    `Klebsiella pneumoniae_Amikacin`,
    `Proteus mirabilis_Amikacin`,
    `Pseudomonas aeruginosa_Amikacin`
  ))
  
  AMKdf <- data.frame(rbind(
    `Enterobacter cloacae_Amikacindf`,
    `Escherichia coli_Amikacindf`,
    `Klebsiella pneumoniae_Amikacindf`,
    `Proteus mirabilis_Amikacindf`,
    `Pseudomonas aeruginosa_Amikacindf`
  ))
  
  AMKdf$org_name <- factor(AMKdf$org_name, levels = AMKdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  AMK_plot <- ggplot(AMKdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = AMKdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Amikacin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(AMK_plot)
  
  assign("AMK_summary",AMK_summary,envir = .GlobalEnv)
  
  df
}
AMK_normal <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    norm_sim(org_fullname,"Enterobacter cloacae",org_fullname,"",AMK,1,1,"Amikacin","_norms_") %>%
    norm_sim(org_fullname,"Escherichia coli",org_fullname,"",AMK,1,1,"Amikacin","_norms_") %>%
    norm_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",AMK,1,1,"Amikacin","_norms_") %>%
    norm_sim(org_fullname,"Proteus mirabilis",org_fullname,"",AMK,1,1,"Amikacin","_norms_") %>%
    norm_sim(org_fullname,"Pseudomonas aeruginosa",org_fullname,"",AMK,1,1,"Amikacin","_norms_")
  
  AMK_norms_summary <- data.frame(rbind(
    `Enterobacter cloacae_Amikacin_norms_`,
    `Escherichia coli_Amikacin_norms_`,
    `Klebsiella pneumoniae_Amikacin_norms_`,
    `Proteus mirabilis_Amikacin_norms_`,
    `Pseudomonas aeruginosa_Amikacin_norms_`
  ))
  
  AMK_norms_df <- data.frame(rbind(
    `Enterobacter cloacae_Amikacin_norms_df`,
    `Escherichia coli_Amikacin_norms_df`,
    `Klebsiella pneumoniae_Amikacin_norms_df`,
    `Proteus mirabilis_Amikacin_norms_df`,
    `Pseudomonas aeruginosa_Amikacin_norms_df`
  ))
  
  AMK_norms_df$org_name <- factor(AMK_norms_df$org_name, levels = AMK_norms_df %>%
                                    distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  assign("AMK_norms_summary",AMK_norms_summary,envir = .GlobalEnv)
  
  df
}
TOB_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Citrobacter freundii",org_fullname,"",TOB,1,3,"Tobramycin") %>%
    res_sim(org_fullname,"Enterobacter cloacae",org_fullname,"",TOB,1,3,"Tobramycin") %>%
    res_sim(org_fullname,"Escherichia coli",org_fullname,"",TOB,1,3,"Tobramycin") %>%
    res_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",TOB,1,3,"Tobramycin") %>%
    res_sim(org_fullname,"Morganella morganii",org_fullname,"",TOB,1,3,"Tobramycin") %>%
    res_sim(org_fullname,"Proteus mirabilis",org_fullname,"",TOB,1,3,"Tobramycin") %>%
    res_sim(org_fullname,"Pseudomonas aeruginosa",org_fullname,"",TOB,1,3,"Tobramycin") %>%
    res_sim(org_fullname,"Serratia marcescens",org_fullname,"",TOB,1,3,"Tobramycin")
  
  TOB_summary <- data.frame(rbind(
    `Citrobacter freundii_Tobramycin`,
    `Enterobacter cloacae_Tobramycin`,
    `Escherichia coli_Tobramycin`,
    `Klebsiella pneumoniae_Tobramycin`,
    `Morganella morganii_Tobramycin`,
    `Proteus mirabilis_Tobramycin`,
    `Pseudomonas aeruginosa_Tobramycin`,
    `Serratia marcescens_Tobramycin`
  ))
  
  TOBdf <- data.frame(rbind(
    `Citrobacter freundii_Tobramycindf`,
    `Enterobacter cloacae_Tobramycindf`,
    `Escherichia coli_Tobramycindf`,
    `Klebsiella pneumoniae_Tobramycindf`,
    `Morganella morganii_Tobramycindf`,
    `Proteus mirabilis_Tobramycindf`,
    `Pseudomonas aeruginosa_Tobramycindf`,
    `Serratia marcescens_Tobramycindf`
  ))
  
  TOBdf$org_name <- factor(TOBdf$org_name, levels = TOBdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  TOB_plot <- ggplot(TOBdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = TOBdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Tobramycin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(TOB_plot)
  
  assign("TOB_summary",TOB_summary,envir = .GlobalEnv)
  
  df
}
TOB_normal <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    norm_sim(org_fullname,"Citrobacter freundii",org_fullname,"",TOB,1,1,"Tobramycin","_norms_") %>%
    norm_sim(org_fullname,"Enterobacter cloacae",org_fullname,"",TOB,1,1,"Tobramycin","_norms_") %>%
    norm_sim(org_fullname,"Escherichia coli",org_fullname,"",TOB,1,1,"Tobramycin","_norms_") %>%
    norm_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",TOB,1,1,"Tobramycin","_norms_") %>%
    norm_sim(org_fullname,"Morganella morganii",org_fullname,"",TOB,1,1,"Tobramycin","_norms_") %>%
    norm_sim(org_fullname,"Proteus mirabilis",org_fullname,"",TOB,1,1,"Tobramycin","_norms_") %>%
    norm_sim(org_fullname,"Pseudomonas aeruginosa",org_fullname,"",TOB,1,1,"Tobramycin","_norms_") %>%
    norm_sim(org_fullname,"Serratia marcescens",org_fullname,"",TOB,1,1,"Tobramycin","_norms_")
  
  TOB_norms_summary <- data.frame(rbind(
    `Citrobacter freundii_Tobramycin_norms_`,
    `Enterobacter cloacae_Tobramycin_norms_`,
    `Escherichia coli_Tobramycin_norms_`,
    `Klebsiella pneumoniae_Tobramycin_norms_`,
    `Morganella morganii_Tobramycin_norms_`,
    `Proteus mirabilis_Tobramycin_norms_`,
    `Pseudomonas aeruginosa_Tobramycin_norms_`,
    `Serratia marcescens_Tobramycin_norms_`
  ))
  
  TOB_norms_df <- data.frame(rbind(
    `Citrobacter freundii_Tobramycin_norms_df`,
    `Enterobacter cloacae_Tobramycin_norms_df`,
    `Escherichia coli_Tobramycin_norms_df`,
    `Klebsiella pneumoniae_Tobramycin_norms_df`,
    `Morganella morganii_Tobramycin_norms_df`,
    `Proteus mirabilis_Tobramycin_norms_df`,
    `Pseudomonas aeruginosa_Tobramycin_norms_df`,
    `Serratia marcescens_Tobramycin_norms_df`
  ))
  
  TOB_norms_df$org_name <- factor(TOB_norms_df$org_name, levels = TOB_norms_df %>%
                                    distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  assign("TOB_norms_summary",TOB_norms_summary,envir = .GlobalEnv)
  
  df
}
NIT_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Citrobacter freundii",org_fullname,"",NIT,1,3,"Nitrofurantoin") %>%
    res_sim(org_fullname,"Enterobacter cloacae",org_fullname,"",NIT,1,3,"Nitrofurantoin") %>%
    res_sim(org_fullname,"^Enterococcus$",org_fullname,"",NIT,1,3,"Nitrofurantoin") %>%
    res_sim(org_fullname,"Escherichia coli",org_fullname,"",NIT,1,3,"Nitrofurantoin") %>%
    res_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",NIT,1,3,"Nitrofurantoin") %>%
    res_sim(org_fullname,"Staphylococcus aureus",org_fullname,"",NIT,1,3,"Nitrofurantoin")
  
  NIT_summary <- data.frame(rbind(
    `Citrobacter freundii_Nitrofurantoin`,
    `Enterobacter cloacae_Nitrofurantoin`,
    `^Enterococcus$_Nitrofurantoin`,
    `Escherichia coli_Nitrofurantoin`,
    `Klebsiella pneumoniae_Nitrofurantoin`,
    `Staphylococcus aureus_Nitrofurantoin`
  ))
  
  NITdf <- data.frame(rbind(
    `Citrobacter freundii_Nitrofurantoindf`,
    `Enterobacter cloacae_Nitrofurantoindf`,
    `^Enterococcus$_Nitrofurantoindf`,
    `Escherichia coli_Nitrofurantoindf`,
    `Klebsiella pneumoniae_Nitrofurantoindf`,
    `Staphylococcus aureus_Nitrofurantoindf`
  ))
  
  NITdf$org_name <- factor(NITdf$org_name, levels = NITdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  NIT_plot <- ggplot(NITdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = NITdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Nitrofurantoin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(NIT_plot)
  
  assign("NIT_summary",NIT_summary,envir = .GlobalEnv)
  
  df
}
NIT_normal <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    norm_sim(org_fullname,"Citrobacter freundii",org_fullname,"",NIT,1,1,"Nitrofurantoin","_norms_") %>%
    norm_sim(org_fullname,"Enterobacter cloacae",org_fullname,"",NIT,1,1,"Nitrofurantoin","_norms_") %>%
    norm_sim(org_fullname,"^Enterococcus$",org_fullname,"",NIT,1,1,"Nitrofurantoin","_norms_") %>%
    norm_sim(org_fullname,"Escherichia coli",org_fullname,"",NIT,1,1,"Nitrofurantoin","_norms_") %>%
    norm_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",NIT,1,1,"Nitrofurantoin","_norms_") %>%
    norm_sim(org_fullname,"Staphylococcus aureus",org_fullname,"",NIT,1,1,"Nitrofurantoin","_norms_")
  
  NIT_norms_summary <- data.frame(rbind(
    `Citrobacter freundii_Nitrofurantoin_norms_`,
    `Enterobacter cloacae_Nitrofurantoin_norms_`,
    `^Enterococcus$_Nitrofurantoin_norms_`,
    `Escherichia coli_Nitrofurantoin_norms_`,
    `Klebsiella pneumoniae_Nitrofurantoin_norms_`,
    `Staphylococcus aureus_Nitrofurantoin_norms_`
  ))
  
  NIT_norms_df <- data.frame(rbind(
    `Citrobacter freundii_Nitrofurantoin_norms_df`,
    `Enterobacter cloacae_Nitrofurantoin_norms_df`,
    `^Enterococcus$_Nitrofurantoin_norms_df`,
    `Escherichia coli_Nitrofurantoin_norms_df`,
    `Klebsiella pneumoniae_Nitrofurantoin_norms_df`,
    `Staphylococcus aureus_Nitrofurantoin_norms_df`
  ))
  
  NIT_norms_df$org_name <- factor(NIT_norms_df$org_name, levels = NIT_norms_df %>%
                                    distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  assign("NIT_norms_summary",NIT_norms_summary,envir = .GlobalEnv)
  
  df
}
SXT_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Citrobacter freundii",org_fullname,"",SXT,3,7,"Trimethoprim-sulfamethoxazole") %>%
    res_sim(org_fullname,"Enterobacter cloacae",org_fullname,"",SXT,3,7,"Trimethoprim-sulfamethoxazole") %>%
    res_sim(org_fullname,"Escherichia coli",org_fullname,"",SXT,3,7,"Trimethoprim-sulfamethoxazole") %>%
    res_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",SXT,3,7,"Trimethoprim-sulfamethoxazole") %>%
    res_sim(org_fullname,"Morganella morganii",org_fullname,"",SXT,3,7,"Trimethoprim-sulfamethoxazole") %>%
    res_sim(org_fullname,"Proteus mirabilis",org_fullname,"",SXT,3,7,"Trimethoprim-sulfamethoxazole") %>%
    res_sim(org_fullname,"Serratia marcescens",org_fullname,"",SXT,3,7,"Trimethoprim-sulfamethoxazole") %>%
    res_sim(org_fullname,"Staphylococcus aureus",org_fullname,"",SXT,1,3,"Trimethoprim-sulfamethoxazole") %>%
    res_sim(org_fullname,"Stenotrophomonas maltophilia",org_fullname,"",SXT,1,3,"Trimethoprim-sulfamethoxazole") %>%
    res_sim(org_fullname,"Streptococcus pneumoniae",org_fullname,"",SXT,3,7,"Trimethoprim-sulfamethoxazole")
  
  SXT_summary <- data.frame(rbind(
    `Citrobacter freundii_Trimethoprim-sulfamethoxazole`,
    `Enterobacter cloacae_Trimethoprim-sulfamethoxazole`,
    `Escherichia coli_Trimethoprim-sulfamethoxazole`,
    `Klebsiella pneumoniae_Trimethoprim-sulfamethoxazole`,
    `Morganella morganii_Trimethoprim-sulfamethoxazole`,
    `Proteus mirabilis_Trimethoprim-sulfamethoxazole`,
    `Serratia marcescens_Trimethoprim-sulfamethoxazole`,
    `Staphylococcus aureus_Trimethoprim-sulfamethoxazole`,
    `Stenotrophomonas maltophilia_Trimethoprim-sulfamethoxazole`,
    `Streptococcus pneumoniae_Trimethoprim-sulfamethoxazole`
  ))
  
  SXTdf <- data.frame(rbind(
    `Citrobacter freundii_Trimethoprim-sulfamethoxazoledf`,
    `Enterobacter cloacae_Trimethoprim-sulfamethoxazoledf`,
    `Escherichia coli_Trimethoprim-sulfamethoxazoledf`,
    `Klebsiella pneumoniae_Trimethoprim-sulfamethoxazoledf`,
    `Morganella morganii_Trimethoprim-sulfamethoxazoledf`,
    `Proteus mirabilis_Trimethoprim-sulfamethoxazoledf`,
    `Serratia marcescens_Trimethoprim-sulfamethoxazoledf`,
    `Staphylococcus aureus_Trimethoprim-sulfamethoxazoledf`,
    `Stenotrophomonas maltophilia_Trimethoprim-sulfamethoxazoledf`,
    `Streptococcus pneumoniae_Trimethoprim-sulfamethoxazoledf`
  ))
  
  SXTdf$org_name <- factor(SXTdf$org_name, levels = SXTdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  SXT_plot <- ggplot(SXTdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = SXTdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Trimethoprim-sulfamethoxazole"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(SXT_plot)
  
  assign("SXT_summary",SXT_summary,envir = .GlobalEnv)
  
  df
}
SXT_normal <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    norm_sim(org_fullname,"Citrobacter freundii",org_fullname,"",SXT,1,1,"Trimethoprim-sulfamethoxazole","_norms_") %>%
    norm_sim(org_fullname,"Enterobacter cloacae",org_fullname,"",SXT,1,1,"Trimethoprim-sulfamethoxazole","_norms_") %>%
    norm_sim(org_fullname,"Escherichia coli",org_fullname,"",SXT,1,1,"Trimethoprim-sulfamethoxazole","_norms_") %>%
    norm_sim(org_fullname,"Klebsiella pneumoniae",org_fullname,"",SXT,1,1,"Trimethoprim-sulfamethoxazole","_norms_") %>%
    norm_sim(org_fullname,"Morganella morganii",org_fullname,"",SXT,1,1,"Trimethoprim-sulfamethoxazole","_norms_") %>%
    norm_sim(org_fullname,"Proteus mirabilis",org_fullname,"",SXT,1,1,"Trimethoprim-sulfamethoxazole","_norms_") %>%
    norm_sim(org_fullname,"Serratia marcescens",org_fullname,"",SXT,1,1,"Trimethoprim-sulfamethoxazole","_norms_") %>%
    norm_sim(org_fullname,"Staphylococcus aureus",org_fullname,"",SXT,1,1,"Trimethoprim-sulfamethoxazole","_norms_") %>%
    norm_sim(org_fullname,"Stenotrophomonas maltophilia",org_fullname,"",SXT,1,1,"Trimethoprim-sulfamethoxazole","_norms_") %>%
    norm_sim(org_fullname,"Streptococcus pneumoniae",org_fullname,"",SXT,1,1,"Trimethoprim-sulfamethoxazole","_norms_")
  
  SXT_norms_summary <- data.frame(rbind(
    `Citrobacter freundii_Trimethoprim-sulfamethoxazole_norms_`,
    `Enterobacter cloacae_Trimethoprim-sulfamethoxazole_norms_`,
    `Escherichia coli_Trimethoprim-sulfamethoxazole_norms_`,
    `Klebsiella pneumoniae_Trimethoprim-sulfamethoxazole_norms_`,
    `Morganella morganii_Trimethoprim-sulfamethoxazole_norms_`,
    `Proteus mirabilis_Trimethoprim-sulfamethoxazole_norms_`,
    `Serratia marcescens_Trimethoprim-sulfamethoxazole_norms_`,
    `Staphylococcus aureus_Trimethoprim-sulfamethoxazole_norms_`,
    `Stenotrophomonas maltophilia_Trimethoprim-sulfamethoxazole_norms_`,
    `Streptococcus pneumoniae_Trimethoprim-sulfamethoxazole_norms_`
  ))
  
  SXT_norms_df <- data.frame(rbind(
    `Citrobacter freundii_Trimethoprim-sulfamethoxazole_norms_df`,
    `Enterobacter cloacae_Trimethoprim-sulfamethoxazole_norms_df`,
    `Escherichia coli_Trimethoprim-sulfamethoxazole_norms_df`,
    `Klebsiella pneumoniae_Trimethoprim-sulfamethoxazole_norms_df`,
    `Morganella morganii_Trimethoprim-sulfamethoxazole_norms_df`,
    `Proteus mirabilis_Trimethoprim-sulfamethoxazole_norms_df`,
    `Serratia marcescens_Trimethoprim-sulfamethoxazole_norms_df`,
    `Staphylococcus aureus_Trimethoprim-sulfamethoxazole_norms_df`,
    `Stenotrophomonas maltophilia_Trimethoprim-sulfamethoxazole_norms_df`,
    `Streptococcus pneumoniae_Trimethoprim-sulfamethoxazole_norms_df`
  ))
  
  SXT_norms_df$org_name <- factor(SXT_norms_df$org_name, levels = SXT_norms_df %>%
                                    distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  assign("SXT_norms_summary",SXT_norms_summary,envir = .GlobalEnv)
  
  df
}

#Precision cross-validation
sim_validate <- function(df,colms,conditions,antibiotic,quoted_ab,n_tested,n_reps,func,abx_name,model="",n_orgs,tol_error=0.1) {
  
  antibiotic <- enquo(antibiotic)
  colms <- enquo(colms)
  
  ref_df <- df %>% filter(!is.na(!!antibiotic))#-------------------------------------Filter for available AST results for target antimicrobial
  
  ref_df$isolate_id2 <- as.character(ref_df$org_name)
  ref_df$isolate_id2[!is.na(ref_df$isolate_id2)] <- 1:sum(!is.na(ref_df$org_name))
  
  acc_df <- data.frame(matrix(ncol=10,nrow=n_reps))
  
  colnames(acc_df) <- c("acc","missed","prop1","prop2",
                        "rate_diff","p_val","ci1","ci2",
                        "n_likelihood","n_simulated")
  
  #Iterative stage
  
  for (i in 1:n_reps) {
    
    print(i)
    
    #Result deletion
    
    test_df1 <- ref_df %>%
      filter(grepl(conditions,!!colms)) %>% 
      group_by(org_fullname) %>% 
      slice_sample(n=-n_tested) %>% 
      mutate(!!antibiotic:=NA)
    
    test_df2 <- ref_df %>% anti_join(test_df1, by="isolate_id2")
    
    te_df <- data.frame(rbind(test_df1,test_df2))
    
    key <- te_df %>% filter(is.na(!!antibiotic)) %>% select(isolate_id2)
    
    #Result simulation
    
    te_df2 <- te_df %>% func()
    
    ref_df_key <- ref_df %>% semi_join(key,by="isolate_id2") %>%
      select(isolate_id2,!!antibiotic) %>% rename(ab = !!antibiotic)
    
    te_df_key <- te_df2 %>% semi_join(key,by="isolate_id2") %>% 
      select(isolate_id2,!!antibiotic) %>% rename(ab2 = !!antibiotic) %>% 
      filter(!is.na(ab2))
    
    ref_df_key <- ref_df_key %>% 
      right_join(te_df_key,by="isolate_id2") %>% 
      select(-ab2)
    
    test_key <- merge(ref_df_key,te_df_key)
    
    #Summary statistics
    
    TP <- nrow(test_key %>% filter(ab=="R" & ab2=="R"))
    FP <- nrow(test_key %>% filter(ab=="R" & ab2=="S"))
    TN <- nrow(test_key %>% filter(ab=="S" & ab2=="S"))
    FN <- nrow(test_key %>% filter(ab=="S" & ab2=="R"))
    
    
    acc <- ( TP + TN ) /( TP + FP + TN + FN )
    
    rate_diff <- nrow(test_key %>% filter(ab2=="R"))/nrow(test_key) - 
      nrow(test_key %>% filter(ab=="R"))/nrow(test_key)
    
    missed <- sum(is.na(test_key$ab2))/nrow(test_key)
    
    prop_1 <- nrow(test_key %>% filter(ab=="R"))  
    prop_2 <- nrow(test_key %>% filter(ab2=="R"))
    row_1 <- sum(!is.na(test_key$ab))
    row_2 <- sum(!is.na(test_key$ab2))
    
    
    
    propt <- prop.test(x=c(prop_1,prop_2),
                       n=c(row_1,row_2))
    
    
    
    prop <- data.frame(t(propt$estimate))
    
    
    
    p <- data.frame(propt$p.value)
    
    
    
    ci <- data.frame(t((propt$conf.int)))
    
    n_likelihood <- nrow(te_df %>% filter( grepl(conditions,!!colms) & 
                                             !is.na(!!antibiotic)))
    
    n_simulated <- nrow(te_df_key)
    
    
    acc_row <- data.frame(cbind(acc,missed,prop,rate_diff,p,ci,
                                n_likelihood,n_simulated))
    
    colnames(acc_row) <- c("acc","missed","prop1","prop2",
                           "rate_diff","p_val","ci1","ci2",
                           "n_likelihood","n_simulated")
    
    acc_df[i,] <- acc_row
    
    
  }
  
  #Summary statistics across validation
  
  m_acc <- mean(acc_df[,1],na.rm=T)
  sd_acc <- sd(acc_df[,1],na.rm=T)
  mean_miss <- nrow(acc_df %>% filter(rate_diff>=0.2|rate_diff<=-0.2))/n_reps
  sd_miss <- sd(acc_df[,2],na.rm=T)
  m_prop1 <- mean(acc_df[,3],na.rm=T)
  sd_prop1 <- sd(acc_df[,3],na.rm=T)
  m_prop2 <- mean(acc_df[,4],na.rm=T)
  sd_prop2 <- sd(acc_df[,4],na.rm=T)
  m_ratediff <- mean(acc_df[,5],na.rm=T)
  sd_ratediff <- sd(acc_df[,5],na.rm=T)
  prop_high_ratediff <- nrow(acc_df %>% filter(rate_diff>=tol_error|rate_diff<=-tol_error))/n_reps
  m_ci_low <- mean(acc_df[,7],na.rm=T)
  sd_ci_low <- sd(acc_df[,7],na.rm=T)
  m_ci_upp <- mean(acc_df[,8],na.rm=T)
  sd_ci_upp <- sd(acc_df[,8],na.rm=T)
  m_n_likelihood <- mean(acc_df[,9],na.rm=T)
  sd_n_likelihood <- sd(acc_df[,9],na.rm=T)
  m_n_simulated <- mean(acc_df[,10],na.rm=T)
  sd_n_simulated <- sd(acc_df[,10],na.rm=T)
  
  valid_frame <- data.frame(ratediff=acc_df[,5])
  valid_frame <- cbind(data.frame(Antimicrobial=abx_name),
                       data.frame(ratediff=valid_frame))
  
  valid_frame <- valid_frame %>%
    group_by(Antimicrobial) %>%
    mutate(outlier = ratediff < quantile(ratediff, .25) - 1.5*IQR(ratediff) | ratediff > quantile(ratediff, .75) + 1.5*IQR(ratediff)) %>%
    ungroup
  
  assign(glue("{quoted_ab}_validation{model}"),valid_frame,envir = .GlobalEnv)
  
  summarydf <- data.frame(cbind(m_acc,sd_acc,mean_miss,sd_miss,
                                m_prop1,sd_prop1,m_prop2,sd_prop2,m_ratediff,sd_ratediff,
                                prop_high_ratediff,m_ci_low,sd_ci_low,m_ci_upp,sd_ci_upp,
                                m_n_likelihood,sd_n_likelihood,
                                m_n_simulated,sd_n_simulated))
  
  colnames(summarydf) <- c("m_acc","sd_acc","mean_miss","sd_miss",
                           "m_prop1","sd_prop1","m_prop2","sd_prop2",
                           "m_ratediff","sd_ratediff",
                           "prop_high_ratediff",
                           "m_ci_low","sd_ci_low","m_ci_upp","sd_ci_upp",
                           "m_n_likelihood","sd_n_likelihood",
                           "m_n_simulated","sd_n_simulated")
  
  rownames(summarydf) <- quoted_ab
  
  summarydf
  
}

#Resource efficiency cross-validation
number_validate <- function(df,conds,ab,quot_ab,iters,simul_func,ab_name,model="",error_marg=0.1,starting_num=2) {
  
  ab <- enquo(ab)
  
  i <- starting_num
  
  #Cross-validation on starting number of specimens
  
  sim <- df %>% sim_validate(org_fullname,conds,!!ab,quot_ab,i,iters,simul_func,ab_name,model,tol_error=error_marg)
  
  prop <- sim$prop_high_ratediff
  
  if (prop == 0) {
    
    print(glue("Minimum testing size to avoid estimation error > {error_marg}
  in {iters} simulation iterations for {ab_name} resistance
  across {conds}
             ≤ {i} specimen(s)"))
    
    data.frame(cbind(ab_name,i))
    
  } else {
    
    #Iterative stage with increasing sample size
    
    while (prop > 0) {
      
      i <- i + 1
      
      sim <- df %>% sim_validate(org_fullname,conds,!!ab,quot_ab,i,iters,simul_func,ab_name,model,tol_error=error_marg)
      
      prop <- sim$prop_high_ratediff
      
      if (i==101) {
        
        break
        
      }
      
    }
    
    print(glue("Minimum testing size to avoid estimation error > {error_marg}
  in {iters} simulation iterations for {ab_name} resistance
  across {conds}
             = {i-1} specimen(s)"))
    
    i <- i-1
    
    #Summary statistic dataframe
    
    df1 <- data.frame(cbind(ab_name,i))
    
  }
  
}

#I to R reassignment for sensitivity analysis

sensitivity_func <- function(df) {
  
  df2 <- df %>% select(PEN:MTR)
  df2[df2=="I"] <- "R"
  df[,17:81] <- df2
  
  df
  
}