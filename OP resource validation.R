######OBSERVED PREVALENCE RESOURCE EFFICIENCE CROSS-VALIDATION

numbnorm_PEN <- m_actual_intr %>% number_validate(
  "(Streptococcus pneumoniae|Staphylococcus aureus)",
  PEN,"PEN",50,PEN_normal,"Benzylpenicillin",error_marg=0.1,starting_num = 2)

numbnorm_AMP <- m_actual_intr %>% number_validate(
  "(Streptococcus pneumoniae|Escherichia coli|^Enterococcus$|Proteus mirabilis)",
  AMP,"AMP",50,AMP_normal,"Ampicillin",error_marg=0.1,starting_num = 2)

numbnorm_OXA <- m_actual_intr %>% number_validate(
  "Staphylococcus aureus",
  OXA,"OXA",50,OXA_normal,"Oxacillin",error_marg=0.1,starting_num = 2)

numbnorm_SAM <- m_actual_intr %>% number_validate(
  "(Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis)",
  SAM,"SAM",50,SAM_normal,"Ampicillin-sulbactam",error_marg=0.1,starting_num = 2)

numbnorm_TZP <- m_actual_intr %>% number_validate(
  "(Escherichia coli|Klebsiella|Proteus|Pseudomonas aeruginosa)",
  TZP,"TZP",50,TZP_normal,"Piperacillin-tazobactam",error_marg=0.1,starting_num = 2)

numbnorm_CZO <- m_actual_intr %>% number_validate(
  "(Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis)",
  CZO,"CZO",50,CZO_normal,"Cefazolin",error_marg=0.1,starting_num = 2)

numbnorm_CXM <- m_actual_intr %>% number_validate(
  "(Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis)",
  CXM,"CXM",50,CXM_normal,"Cefuroxime",error_marg=0.1,starting_num = 2)

numbnorm_CRO <- m_actual_intr %>% number_validate(
  "(Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis)",
  CRO,"CRO",50,CRO_normal,"Ceftriaxone",error_marg=0.1,starting_num = 2)

numbnorm_CAZ <- m_actual_intr %>% number_validate(
  "(Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis|Pseudomonas aeruginosa)",
  CAZ,"CAZ",50,CAZ_normal,"Ceftazidime",error_marg=0.1,starting_num = 2)

numbnorm_FEP <- m_actual_intr %>% number_validate(
  "(Citrobacter freundii|Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Morganella morganii|Proteus mirabilis|Pseudomonas aeruginosa|Serratia marcescens)",
  FEP,"FEP",50,FEP_normal,"Cefepime",error_marg=0.1,starting_num = 2)

numbnorm_MEM <- m_actual_intr %>% number_validate(
  "(Citrobacter freundii|Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Morganella morganii|Proteus mirabilis|Pseudomonas aeruginosa|Serratia marcescens)",
  MEM,"MEM",50,MEM_normal,"Meropenem",error_marg=0.1,starting_num = 2)

numbnorm_CIP <- m_actual_intr %>% number_validate(
  "(Citrobacter freundii|Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Morganella morganii|Proteus mirabilis|Pseudomonas aeruginosa|Serratia marcescens)",
  CIP,"CIP",50,CIP_normal,"Ciprofloxacin",error_marg=0.1,starting_num = 2)

numbnorm_LVX <- m_actual_intr %>% number_validate(
  "(Streptococcus pneumoniae|Staphylococcus aureus)",
  LVX,"LVX",50,LVX_normal,"Levofloxacin",error_marg=0.1,starting_num = 2)

numbnorm_ERY <- m_actual_intr %>% number_validate(
  "(Streptococcus pneumoniae|Streptococcus Group B|Staphylococcus aureus)",
  ERY,"ERY",50,ERY_normal,"Erythromycin",error_marg=0.1,starting_num = 2)

numbnorm_CLI <- m_actual_intr %>% number_validate(
  "(Streptococcus Group B|Staphylococcus aureus)",
  CLI,"CLI",50,CLI_normal,"Clindamycin",error_marg=0.1,starting_num = 2)

numbnorm_TCY <- m_actual_intr %>% number_validate(
  "(Streptococcus pneumoniae|Staphylococcus aureus)",
  TCY,"TCY",50,TCY_normal,"Tetracycline",error_marg=0.1,starting_num = 2)

numbnorm_VAN <- m_actual_intr %>% number_validate(
  "^Enterococcus$",
  VAN,"VAN",50,VAN_normal,"Vancomycin",error_marg=0.1,starting_num = 2)

numbnorm_RIF <- m_actual_intr %>% number_validate(
  "Staphylococcus aureus",
  RIF,"RIF",50,RIF_normal,"Rifampicin",error_marg=0.1,starting_num = 2)

numbnorm_GEN <- m_actual_intr %>% number_validate(
  "(Citrobacter freundii|Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Morganella morganii|Proteus mirabilis|Pseudomonas aeruginosa|Serratia marcescens|Staphylococcus aureus)",
  GEN,"GEN",50,GEN_normal,"Gentamicin",error_marg=0.1,starting_num = 2)

numbnorm_AMK <- m_actual_intr %>% number_validate(
  "(Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Proteus mirabilis|Pseudomonas aeruginosa)",
  AMK,"AMK",50,AMK_normal,"Amikacin",error_marg=0.1,starting_num = 2)

numbnorm_TOB <- m_actual_intr %>% number_validate(
  "(Citrobacter freundii|Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Morganella morganii|Proteus mirabilis|Pseudomonas aeruginosa|Serratia marcescens)",
  TOB,"TOB",50,TOB_normal,"Tobramycin",error_marg=0.1,starting_num = 2)

numbnorm_NIT <- m_actual_intr %>% number_validate(
  "(Citrobacter freundii|Enterobacter cloacae|^Enterococcus$|Escherichia coli|Klebsiella pneumoniae|Staphylococcus aureus)",
  NIT,"NIT",50,NIT_normal,"Nitrofurantoin",error_marg=0.1,starting_num = 2)

numbnorm_SXT <- m_actual_intr %>% number_validate(
  "(Citrobacter freundii|Enterobacter cloacae|Escherichia coli|Klebsiella pneumoniae|Morganella morganii|Proteus mirabilis|Serratia marcescens|Stenotrophomonas maltophilia|Streptococcus pneumoniae|Staphylococcus aureus)",
  SXT,"SXT",50,SXT_normal,"Trimethoprim-sulfamethoxazole",error_marg=0.1,starting_num = 2)

numbnorms <- data.frame(rbind(numbnorm_PEN,numbnorm_AMP,numbnorm_OXA,
                              numbnorm_SAM,numbnorm_TZP,numbnorm_CZO,
                              numbnorm_CXM,numbnorm_CRO,numbnorm_CAZ,
                              numbnorm_FEP,numbnorm_MEM,numbnorm_CIP,
                              numbnorm_LVX,numbnorm_ERY,numbnorm_CLI,
                              numbnorm_TCY,numbnorm_VAN,numbnorm_RIF,
                              numbnorm_GEN,numbnorm_AMK,numbnorm_TOB,
                              numbnorm_NIT,numbnorm_SXT))