#####BEAR RESOURCE EFFICIENCY CROSS-VALIDATION################################################

#Penicillin
numbsim_PEN <- m_actual_intr %>% number_validate(
  "(Streptococcus pneumoniae|Staphylococcus aureus)",
  PEN,"PEN",50,PEN_simul,"Benzylpenicillin",error_marg=0.1,starting_num = 2)

#Ampicillin
numbsim_AMP <- m_actual_intr %>% number_validate(
  "(Streptococcus pneumoniae|Escherichia coli|^Enterococcus$|Proteus mirabilis)",
  AMP,"AMP",50,AMP_simul,"Ampicillin",error_marg=0.1,starting_num = 2)

#Oxacillin
numbsim_OXA <- m_actual_intr %>% number_validate(
  "Staphylococcus aureus",
  OXA,"OXA",50,OXA_simul,"Oxacillin",error_marg=0.1,starting_num = 2)

#Ampicillin-sulbactam
numbsim_SAM <- m_actual_intr %>% number_validate(
  "(Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis)",
  SAM,"SAM",50,SAM_simul,"Ampicillin-sulbactam",error_marg=0.1,starting_num = 2)

#Piperacillin-tazobactam
numbsim_TZP <- m_actual_intr %>% number_validate(
  "(Escherichia coli|Klebsiella|Proteus|Pseudomonas aeruginosa)",
  TZP,"TZP",50,TZP_simul,"Piperacillin-tazobactam",error_marg=0.1,starting_num = 2)

#Cefazolin
numbsim_CZO <- m_actual_intr %>% number_validate(
  "(Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis)",
  CZO,"CZO",50,CZO_simul,"Cefazolin",error_marg=0.1,starting_num = 2)

#Cefuroxime
numbsim_CXM <- m_actual_intr %>% number_validate(
  "(Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis)",
  CXM,"CXM",50,CXM_simul,"Cefuroxime",error_marg=0.1,starting_num = 2)

#Ceftriaxone
numbsim_CRO <- m_actual_intr %>% number_validate(
  "(Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis)",
  CRO,"CRO",50,CRO_simul,"Ceftriaxone",error_marg=0.1,starting_num = 2)

#Ceftazidime
numbsim_CAZ <- m_actual_intr %>% number_validate(
  "(Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis|Pseudomonas aeruginosa)",
  CAZ,"CAZ",50,CAZ_simul,"Ceftazidime",error_marg=0.1,starting_num = 2)

#Cefepime
numbsim_FEP <- m_actual_intr %>% number_validate(
  "(Citrobacter freundii|Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Morganella morganii|Proteus mirabilis|Pseudomonas aeruginosa|Serratia marcescens)",
  FEP,"FEP",50,FEP_simul,"Cefepime",error_marg=0.1,starting_num = 2)

#Meropenem
numbsim_MEM <- m_actual_intr %>% number_validate(
  "(Citrobacter freundii|Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Morganella morganii|Proteus mirabilis|Pseudomonas aeruginosa|Serratia marcescens)",
  MEM,"MEM",50,MEM_simul,"Meropenem",error_marg=0.1,starting_num = 2)

#Ciprofloxacin
numbsim_CIP <- m_actual_intr %>% number_validate(
  "(Citrobacter freundii|Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Morganella morganii|Proteus mirabilis|Pseudomonas aeruginosa|Serratia marcescens)",
  CIP,"CIP",50,CIP_simul,"Ciprofloxacin",error_marg=0.1,starting_num = 2)

#Levofloxacin
numbsim_LVX <- m_actual_intr %>% number_validate(
  "(Streptococcus pneumoniae|Staphylococcus aureus)",
  LVX,"LVX",50,LVX_simul,"Levofloxacin",error_marg=0.1,starting_num = 2)

#Erythromycin
numbsim_ERY <- m_actual_intr %>% number_validate(
  "(Streptococcus pneumoniae|Streptococcus Group B|Staphylococcus aureus)",
  ERY,"ERY",50,ERY_simul,"Erythromycin",error_marg=0.1,starting_num = 2)

#Clindamycin
numbsim_CLI <- m_actual_intr %>% number_validate(
  "(Streptococcus Group B|Staphylococcus aureus)",
  CLI,"CLI",50,CLI_simul,"Clindamycin",error_marg=0.1,starting_num = 2)

#Tetracycline
numbsim_TCY <- m_actual_intr %>% number_validate(
  "(Streptococcus pneumoniae|Staphylococcus aureus)",
  TCY,"TCY",50,TCY_simul,"Tetracycline",error_marg=0.1,starting_num = 2)

#Vancomycin
numbsim_VAN <- m_actual_intr %>% number_validate(
  "^Enterococcus$",
  VAN,"VAN",50,VAN_simul,"Vancomycin",error_marg=0.1,starting_num = 2)

#Rifampicin
numbsim_RIF <- m_actual_intr %>% number_validate(
  "Staphylococcus aureus",
  RIF,"RIF",50,RIF_simul,"Rifampicin",error_marg=0.1,starting_num = 2)

#Gentamicin
numbsim_GEN <- m_actual_intr %>% number_validate(
  "(Citrobacter freundii|Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Morganella morganii|Proteus mirabilis|Pseudomonas aeruginosa|Serratia marcescens|Staphylococcus aureus)",
  GEN,"GEN",50,GEN_simul,"Gentamicin",error_marg=0.1,starting_num = 2)

#Amikacin
numbsim_AMK <- m_actual_intr %>% number_validate(
  "(Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis|Pseudomonas aeruginosa)",
  AMK,"AMK",50,AMK_simul,"Amikacin",error_marg=0.1,starting_num = 2)

#Tobramycin
numbsim_TOB <- m_actual_intr %>% number_validate(
  "(Citrobacter freundii|Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Morganella morganii|Proteus mirabilis|Pseudomonas aeruginosa|Serratia marcescens)",
  TOB,"TOB",50,TOB_simul,"Tobramycin",error_marg=0.1,starting_num = 2)

#Nitrofurantoin
numbsim_NIT <- m_actual_intr %>% number_validate(
  "(Citrobacter freundii|Enterobacter cloacae|^Enterococcus$|Escherichia coli|Klebsiella pneumoniae|Staphylococcus aureus)",
  NIT,"NIT",50,NIT_simul,"Nitrofurantoin",error_marg=0.1,starting_num = 2)

#Co-trimoxazole
numbsim_SXT <- m_actual_intr %>% number_validate(
  "(Citrobacter freundii|Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Morganella morganii|Proteus mirabilis|Serratia marcescens|Stenotrophomonas maltophilia|Streptococcus pneumoniae|Staphylococcus aureus)",
  SXT,"SXT",50,SXT_simul,"Trimethoprim-sulfamethoxazole",error_marg=0.1,starting_num = 2)

numbsims <- data.frame(rbind(numbsim_PEN,numbsim_AMP,numbsim_OXA,
                             numbsim_SAM,numbsim_TZP,numbsim_CZO,
                             numbsim_CXM,numbsim_CRO,numbsim_CAZ,
                             numbsim_FEP,numbsim_MEM,numbsim_CIP,
                             numbsim_LVX,numbsim_ERY,numbsim_CLI,
                             numbsim_TCY,numbsim_VAN,numbsim_RIF,
                             numbsim_GEN,numbsim_AMK,numbsim_TOB,
                             numbsim_NIT,numbsim_SXT))