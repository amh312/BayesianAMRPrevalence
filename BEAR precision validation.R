#BAYESIAN ESTIMATION OF ANTIMICROBIAL RESISTANCE (BEAR) - PRECISION RUN

##Set up df to record aborted runs

missings <- data.frame(matrix(ncol=2,nrow=0))
colnames(missings) <- c("Antimicrobial","Organism")

##Precision cross-validation for each agent

###Penicillin (BEAR)
PENsim <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Streptococcus pneumoniae|Staphylococcus aureus)",
               PEN,"PEN",20,50,PEN_simul,"Benzylpenicillin")

###Ampicillin (BEAR)
AMPsim <- m_actual_intr %>%
  sim_validate(org_fullname,
               "(Streptococcus pneumoniae|Escherichia coli|^Enterococcus$|Proteus mirabilis)",
               AMP,"AMP",20,50,AMP_simul,"Ampicillin")

###Oxacillin (BEAR)
OXAsim <- m_actual_intr %>%
  sim_validate(org_fullname,
               "Staphylococcus aureus",
               OXA,"OXA",20,50,OXA_simul,"Oxacillin")

###Ampicillin-sulbactam (BEAR)
SAMsim <- m_actual_intr %>%
  sim_validate(org_fullname,
               "(Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis)",
               SAM,"SAM",20,50,SAM_simul,"Ampicillin-sulbactam")

###Piperacillin-tazobactam (BEAR)
TZPsim <- m_actual_intr %>%
  sim_validate(org_fullname,
               "(Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis|Pseudomonas aeruginosa)",
               TZP,"TZP",20,50,TZP_simul,"Piperacillin-tazobactam")

###Cefazolin (BEAR)
CZOsim <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis)",
               CZO,"CZO",20,50,CZO_simul,"Cefazolin")

###Cefuroxime (BEAR)
CXMsim <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis)",
               CXM,"CXM",20,50,CXM_simul,"Cefuroxime")

###Ceftriaxone (BEAR)
CROsim <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis)",
               CRO,"CRO",20,50,CRO_simul,"Ceftriaxone")

###Ceftazidime (BEAR)
CAZsim <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis|Pseudomonas aeruginosa)",
               CAZ,"CAZ",20,50,CAZ_simul,"Ceftazidime")

###Cefepime (BEAR)
FEPsim <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Citrobacter freundii|Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Morganella morganii|Proteus mirabilis|Pseudomonas aeruginosa|Serratia marcescens)",
               FEP,"FEP",20,50,FEP_simul,"Cefepime")

###Meropenem (BEAR)
MEMsim <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Citrobacter freundii|Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Morganella morganii|Proteus mirabilis|Pseudomonas aeruginosa|Serratia marcescens)",
               MEM,"MEM",20,50,MEM_simul,"Meropenem")

###Ciprofloxacin (BEAR)
CIPsim <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Citrobacter freundii|Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Morganella morganii|Proteus mirabilis|Pseudomonas aeruginosa|Serratia marcescens)",
               CIP,"CIP",20,50,CIP_simul,"Ciprofloxacin")

###Levofloxacin (BEAR)
LVXsim <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Streptococcus pneumoniae|Staphylococcus aureus)",
               LVX,"LVX",20,50,LVX_simul,"Levofloxacin")

###Erythromycin (BEAR)
ERYsim <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Streptococcus pneumoniae|Streptococcus Group B|Staphylococcus aureus)",
               ERY,"ERY",20,50,ERY_simul,"Erythromycin")

###Clindamycin (BEAR)
CLIsim <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Streptococcus Group B|Staphylococcus aureus)",
               CLI,"CLI",20,50,CLI_simul,"Clindamycin")

###Tetracycline (BEAR)
TCYsim <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Streptococcus pneumoniae|Staphylococcus aureus)",
               TCY,"TCY",20,50,TCY_simul,"Tetracycline")

###Vancomycin (BEAR)
VANsim <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "^Enterococcus$",
               VAN,"VAN",20,50,VAN_simul,"Vancomycin")

###Rifampicin (BEAR)
RIFsim <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "Staphylococcus aureus",
               RIF,"RIF",20,50,RIF_simul,"Rifampicin")

###Gentamicin (BEAR)
GENsim <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Citrobacter freundii|Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Morganella morganii|Proteus mirabilis|Pseudomonas aeruginosa|Serratia marcescens|Staphylococcus aureus)",
               GEN,"GEN",20,50,GEN_simul,"Gentamicin")

###Amikacin (BEAR)
AMKsim <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis|Pseudomonas aeruginosa)",
               AMK,"AMK",20,50,AMK_simul,"Amikacin")

###Tobramycin (BEAR)
TOBsim <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Citrobacter freundii|Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Morganella morganii|Proteus mirabilis|Pseudomonas aeruginosa|Serratia marcescens)",
               TOB,"TOB",20,50,TOB_simul,"Tobramycin")

###Nitrofurantoin (BEAR)
NITsim <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Citrobacter freundii|Enterobacter cloacae|^Enterococcus$|Escherichia coli|Klebsiella pneumoniae|Staphylococcus aureus)",
               NIT,"NIT",20,50,NIT_simul,"Nitrofurantoin")

###Co-trimoxazole (BEAR)
SXTsim <- m_actual_intr %>% 
  sim_validate(org_fullname,
               "(Citrobacter freundii|Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Morganella morganii|Proteus mirabilis|Serratia marcescens|Stenotrophomonas maltophilia|Streptococcus pneumoniae|Staphylococcus aureus)",
               SXT,"SXT",20,50,SXT_simul,"Co-trimoxazole")

##Summary data frames (BEAR)

###Bind together agent summary dfs
sims <- data.frame(rbind(PENsim,AMPsim,OXAsim,SAMsim,TZPsim,CZOsim,CXMsim,            
                         CROsim,CAZsim,FEPsim,MEMsim,CIPsim,LVXsim,ERYsim,CLIsim,
                         TCYsim,VANsim,RIFsim,GENsim,AMKsim,TOBsim,NITsim,SXTsim))

#Round for tabulation
sims <- sims %>% round(3)

#Label summary table and tidy up
colnames(sims) <- c("Accuracy (mean)","Accuracy(sd)",
                    "Missing data rate (mean)","Missing data rate (SD)",
                    "Actual resistance rate (mean)","Actual resistance rate (SD)",
                    "Simulated resistance rate (mean)","Simulated resistance rate (SD)",
                    "Resistance rate difference (mean)","Resistance rate difference (SD)",
                    "Prop rate difference > 0.1","Lower CI bound (mean)","Lower CI bound (SD)",
                    "Upper CI bound (mean)","Upper CI bound (SD)", "n likelihood (mean)",
                    "n likelihood (SD)",
                    "n simulated (mean)","n simulated (SD)")
sims$Antimicrobial <- ab_name(rownames(sims))
sims <- sims %>% relocate(Antimicrobial, .before="Accuracy (mean)")

#Tabulate counts of missing runs
sim_missings <- missings %>%                                                                       #Complete check for failed iterations
  count(Antimicrobial,Organism) %>% 
  rename(n_missing_iterations = n)
