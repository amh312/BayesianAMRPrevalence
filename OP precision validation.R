#ESTIMATION BASED ON OBSERVED PREVALENCE (EOP) - PRECISION CROSS-VALIDATION

##Set dataframe to record aborted runs
missings <- data.frame(matrix(ncol=2,nrow=0))                                      #Missing iterations empty data frame
colnames(missings) <- c("Antimicrobial","Organism")

##EOP Precision CV check for each agent

###Penicillin (EOP)
PENnorm <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Streptococcus pneumoniae|Staphylococcus aureus)",
               PEN,"PEN",20,50,PEN_normal,"Benzylpenicillin","_norms")

###Ampicillin (EOP)
AMPnorm <- m_actual_intr %>%
  sim_validate(org_fullname,
               "(^Enterococcus$|Streptococcus pneumoniae|Escherichia coli|Proteus mirabilis)",
               AMP,"AMP",20,50,AMP_normal,"Ampicillin","_norms")

###Oxacillin (EOP)
OXAnorm <- m_actual_intr %>%
  sim_validate(org_fullname,
               "Staphylococcus aureus",
               OXA,"OXA",20,50,OXA_normal,"Oxacillin","_norms")

###Ampicillin-sulbactam (EOP)
SAMnorm <- m_actual_intr %>%
  sim_validate(org_fullname,
               "(Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis)",
               SAM,"SAM",20,50,SAM_normal,"Ampicillin-sulbactam","_norms")

###Piperacillin-tazobactam (EOP)
TZPnorm <- m_actual_intr %>%
  sim_validate(org_fullname,
               "(Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis|Pseudomonas aeruginosa)",
               TZP,"TZP",20,50,TZP_normal,"Piperacillin-tazobactam","_norms")

###Cefazolin (EOP)
CZOnorm <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis)",
               CZO,"CZO",20,50,CZO_normal,"Cefazolin","_norms")

###Cefuroxime (EOP)
CXMnorm <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis)",
               CXM,"CXM",20,50,CXM_normal,"Cefuroxime","_norms")

###Ceftriaxone (EOP)
CROnorm <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis)",
               CRO,"CRO",20,50,CRO_normal,"Ceftriaxone","_norms")

###Ceftazidime (EOP)
CAZnorm <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis|Pseudomonas aeruginosa)",
               CAZ,"CAZ",20,50,CAZ_normal,"Ceftazidime","_norms")

###Cefepime (EOP)
FEPnorm <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Citrobacter freundii|Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Morganella morganii|Proteus mirabilis|Pseudomonas aeruginosa|Serratia marcescens)",
               FEP,"FEP",20,50,FEP_normal,"Cefepime","_norms")

###Meropenem (EOP)
MEMnorm <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Citrobacter freundii|Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Morganella morganii|Proteus mirabilis|Pseudomonas aeruginosa|Serratia marcescens)",
               MEM,"MEM",20,50,MEM_normal,"Meropenem","_norms")

###Ciprofloxacin (EOP)
CIPnorm <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Citrobacter freundii|Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Morganella morganii|Proteus mirabilis|Pseudomonas aeruginosa|Serratia marcescens)",
               CIP,"CIP",20,50,CIP_normal,"Ciprofloxacin","_norms")

###Levofloxacin (EOP)
LVXnorm <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Streptococcus pneumoniae|Staphylococcus aureus)",
               LVX,"LVX",20,50,LVX_normal,"Levofloxacin","_norms")

###Erythromycin (EOP)
ERYnorm <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Streptococcus pneumoniae|Streptococcus Group B|Staphylococcus aureus)",
               ERY,"ERY",20,50,ERY_normal,"Erythromycin","_norms")

###Clindamycin (EOP)
CLInorm <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Streptococcus Group B|Staphylococcus aureus)",
               CLI,"CLI",20,50,CLI_normal,"Clindamycin","_norms")

###Tetracycline (EOP)
TCYnorm <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Streptococcus pneumoniae|Staphylococcus aureus)",
               TCY,"TCY",20,50,TCY_normal,"Tetracycline","_norms")

###Vancomycin (EOP)
VANnorm <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "^Enterococcus$",
               VAN,"VAN",20,50,VAN_normal,"Vancomycin","_norms")

###Rifampicin (EOP)
RIFnorm <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "Staphylococcus aureus",
               RIF,"RIF",20,50,RIF_normal,"Rifampicin","_norms")

###Gentamicin (EOP)
GENnorm <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Citrobacter freundii|Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Morganella morganii|Proteus mirabilis|Pseudomonas aeruginosa|Serratia marcescens|Staphylococcus aureus)",
               GEN,"GEN",20,50,GEN_normal,"Gentamicin","_norms")

###Amikacin (EOP)
AMKnorm <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis|Pseudomonas aeruginosa)",
               AMK,"AMK",20,50,AMK_normal,"Amikacin","_norms")

###Tobramycin (EOP)
TOBnorm <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Citrobacter freundii|Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Morganella morganii|Proteus mirabilis|Pseudomonas aeruginosa|Serratia marcescens)",
               TOB,"TOB",20,50,TOB_normal,"Tobramycin","_norms")

###Nitrofurantoin (EOP)
NITnorm <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Citrobacter freundii|Enterobacter cloacae|^Enterococcus$|Escherichia coli|Klebsiella pneumoniae|Staphylococcus aureus)",
               NIT,"NIT",20,50,NIT_normal,"Nitrofurantoin","_norms")

###Co-trimoxazole (EOP)
SXTnorm <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Citrobacter freundii|Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Morganella morganii|Proteus mirabilis|Serratia marcescens|Stenotrophomonas maltophilia|Streptococcus pneumoniae|Staphylococcus aureus)",
               SXT,"SXT",20,50,SXT_normal,"Co-trimoxazole","_norms")

##Compile summary dataframes

###Binding agent dataframes
norms <- data.frame(rbind(PENnorm,AMPnorm,OXAnorm,SAMnorm,TZPnorm,CZOnorm,CXMnorm,            #Assemble precision analysis data frame
                          CROnorm,CAZnorm,FEPnorm,MEMnorm,CIPnorm,LVXnorm,ERYnorm,CLInorm,
                          TCYnorm,VANnorm,RIFnorm,GENnorm,AMKnorm,TOBnorm,NITnorm,SXTnorm))

###Rounding for presentation
norms <- norms %>% round(3)

###Labelling dataframe
colnames(norms) <- c("Accuracy (mean)","Accuracy(sd)",
                     "Missing data rate (mean)","Missing data rate (SD)",
                     "Actual resistance rate (mean)","Actual resistance rate (SD)",
                     "Simulated resistance rate (mean)","Simulated resistance rate (SD)",
                     "Resistance rate difference (mean)","Resistance rate difference (SD)",
                     "Prop rate difference > 0.1","Lower CI bound (mean)","Lower CI bound (SD)",
                     "Upper CI bound (mean)","Upper CI bound (SD)", "n likelihood (mean)",
                     "n likelihood (SD)",
                     "n simulated (mean)","n simulated (SD)")
norms$Antimicrobial <- ab_name(rownames(norms))
norms <- norms %>% relocate(Antimicrobial, .before="Accuracy (mean)")

###Record aborted runs for report
norm_missings <- missings %>%                                                                       #Complete check for failed iterations
  count(Antimicrobial,Organism) %>% 
  rename(n_missing_iterations = n)

write_csv(norms,"BEAR_norms.csv")
