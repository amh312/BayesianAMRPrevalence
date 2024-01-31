######PLOTS

#Precision boxplot (observed prevalence)

valid_norms <- data.frame(rbind(PEN_validation_norms,AMP_validation_norms,
                                OXA_validation_norms,SAM_validation_norms,
                                TZP_validation_norms,CZO_validation_norms,
                                CXM_validation_norms,CRO_validation_norms,
                                CAZ_validation_norms,FEP_validation_norms,
                                MEM_validation_norms,CIP_validation_norms,
                                LVX_validation_norms,ERY_validation_norms,
                                CLI_validation_norms,TCY_validation_norms,
                                VAN_validation_norms,RIF_validation_norms,
                                GEN_validation_norms,AMK_validation_norms,
                                TOB_validation_norms,NIT_validation_norms,
                                SXT_validation_norms))

valid_norms <- valid_norms %>% mutate(ratediff=ratediff*100)

levels_order <- valid_norms %>% group_by(Antimicrobial) %>% 
  mutate(iqr = IQR(ratediff)) %>% arrange(desc(iqr)) %>% 
  distinct(Antimicrobial) %>% unlist()

valid_norms$Antimicrobial <- factor(valid_norms$Antimicrobial, levels = levels_order) %>% 
  fct_rev()

qual_col_pals <-  brewer.pal.info[brewer.pal.info$category == 'qual',]

col_vector <-  unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

valid_norms_plot <- ggplot(valid_norms, aes(x=Antimicrobial,y=ratediff,fill=Antimicrobial))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(data = valid_norms %>% dplyr::filter(outlier),alpha=0.25) +
  coord_flip() +
  ggtitle(glue("AMR prevalence differences between observed prevalence estimates and actual data"))+
  xlab("Antimicrobial agent") +
  ylab("Prevalence estimation error (%)") +
  ylim(-30,30) +
  scale_fill_manual(values=col_vector) +
  theme(legend.position = "none") +
  geom_hline(yintercept=0,linetype=2) +
  geom_hline(yintercept=-10,linetype=3,colour="red") +
  geom_hline(yintercept=10,linetype=3,colour="red")

print(valid_norms_plot)





#Precision validation boxplot (BEAR)

valid <- data.frame(rbind(PEN_validation,AMP_validation,OXA_validation,
                          SAM_validation,TZP_validation,CZO_validation,
                          CXM_validation,CRO_validation,CAZ_validation,
                          FEP_validation,MEM_validation,CIP_validation,
                          LVX_validation,ERY_validation,CLI_validation,
                          TCY_validation,VAN_validation,RIF_validation,
                          GEN_validation,AMK_validation,TOB_validation,
                          NIT_validation,SXT_validation))
valid <- valid %>% mutate(ratediff=ratediff*100)
valid$Antimicrobial <- factor(valid$Antimicrobial, levels = levels_order) %>% 
  fct_rev()
qual_col_pals <-  brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <-  unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
valid_plot <- ggplot(valid, aes(x=Antimicrobial,y=ratediff,fill=Antimicrobial))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(data = valid %>% dplyr::filter(outlier),alpha=0.25) +
  coord_flip() +
  ggtitle(glue("AMR prevalence differences between BEAR estimates and actual data"))+
  xlab("Antimicrobial agent") +
  ylab("Prevalence estimation error (%)") +
  ylim(-30,30) +
  scale_fill_manual(values=col_vector) +
  theme(legend.position = "none") +
  geom_hline(yintercept=0,linetype=2) +
  geom_hline(yintercept=-10,linetype=3,colour="red") +
  geom_hline(yintercept=10,linetype=3,colour="red")
print(valid_plot)





#Resource dotplot

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