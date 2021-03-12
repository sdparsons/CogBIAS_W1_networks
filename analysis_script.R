# .R version of code for  "Emotional information-processing correlates of positive mental health in adolescence: A network analysis approach" (Parsons, Songco, Booth, & Fox, 2020, submitted)

# contains all code, without the manuscript text (as can be found in the Rmarkdown file "CCBH in adolescence - Network analysis.Rmd").


###########################################################################
###########################################################################

# include code that will install packages if not installed.

if(!"devtools" %in% rownames(installed.packages())) install.packages("devtools")

if(!"papaja" %in% rownames(installed.packages())) devtools::install_github("crsh/papaja")
if(!"tidyverse" %in% rownames(installed.packages())) install.packages("tidyverse")
if(!"foreign" %in% rownames(installed.packages())) install.packages("foreign")
if(!"bootnet" %in% rownames(installed.packages())) install.packages("bootnet")
if(!"igraph" %in% rownames(installed.packages())) install.packages("igraph")
if(!"mgm" %in% rownames(installed.packages())) install.packages("mgm")
if(!"qgraph" %in% rownames(installed.packages())) install.packages("qgraph")
if(!"NetworkComparisonTest" %in% rownames(installed.packages())) install.packages("NetworkComparisonTest")
if(!"parallel" %in% rownames(installed.packages())) install.packages("parallel")
if(!"gridExtra" %in% rownames(installed.packages())) install.packages("gridExtra")
if(!"Cairo" %in% rownames(installed.packages())) install.packages("Cairo")
if(!"splithalf" %in% rownames(installed.packages())) devtools::install_github("sdparsons/splithalf")
if(!"psych" %in% rownames(installed.packages())) install.packages("psych")
if(!"scales" %in% rownames(installed.packages())) install.packages("scales")
if(!"RColorBrewer" %in% rownames(installed.packages())) install.packages("RColorBrewer")
if(!"corrr" %in% rownames(installed.packages())) install.packages("corrr")


library("papaja")    # for APA formatting awesome-ness in Rmarkdown
library("tidyverse") # for restructuring data
library("foreign")   # for using SPSS data
library("bootnet")   # for network analysis
library("igraph")    # used for spinglass
library("mgm")       # for network analysis
library("qgraph")    # for network analysis - specifically averageLayout
library("NetworkComparisonTest") # for network model comparisons
library("parallel") # used to get number of cores for bootnet
library("gridExtra") # used for combining plots
library("Cairo")
library("splithalf") # for internal consistency of measures
library("psych")
library("scales") # for plotting the interaction figures
library("RColorBrewer") # for plotting the interaction figures
library("corrr") # for fig 1
library("glasso") # to cite


###########################################################################
###########################################################################


# This chunk loads the raw data. Note that the data provided with the paper is the data required to generate this manuscript and is not the entire CogBIAS raw dataset.
# note that the data included in teh questionnaire files is not the full dataset

# demographics, questionnaire responses, and task data all stored in separate SPSS data files.
demographics <- read.csv("Data/W1_demographics.csv", stringsAsFactors =  FALSE)
scores2 <- read.spss("Data/W1_Q2.sav", to.data.frame = TRUE)
AIBQ <- read.spss("Data/W1_AIBQ.sav", to.data.frame = TRUE)
Memory <- read.spss("Data/W1_SRET.sav", to.data.frame = TRUE)
MEM_raw <- read.spss("Data/Start_SRET_Masterfile.sav", to.data.frame = TRUE)
AB_raw <- read.spss("Data/Dot-Probe SPSS Start.sav", to.data.frame = TRUE)
AIBQ_raw <- read.spss("Data/Start_AIBQ_rawdata.sav", to.data.frame = TRUE)

###########################################################################
###########################################################################

# note: this chunk demonstrated how the variables of interest were extracted from the datasets. 

# first extract the variables of interest
Qdat <- data.frame(subject = scores2$Master_subject,
                   MHC_Tot = scores2$MCSHC_Total)

AIBQdat <- data.frame(subjectAIBQ = AIBQ$Master_subject,
                      Pos_Soc    = AIBQ$Interpretation_Pos_Social,
                      Pos_nonsoc = AIBQ$Interpretation_Pos_Nonsocial,
                      Neg_soc    = AIBQ$Interpretation_Neg_Social,
                      Neg_nonsoc = AIBQ$Interpretation_Neg_Nonsocial)

Memdat <- data.frame(subjectMEM = Memory$Master_subject,
                     MEM_pos = Memory$PosEndorsedAndRecalled,
                     MEM_neg = Memory$NegEndorsedAndRecalled)

# combine into one dataframe
CCBH <- cbind(Qdat,AIBQdat,Memdat)

# quick checks that all IDs match
sum(isFALSE(CCBH$subject == CCBH$subjectAIBQ)) 
sum(isFALSE(CCBH$subjectQ == CCBH$subjectMEM))
sum(isFALSE(CCBH$subjectAIBQ == CCBH$subjectMEM))

# remove previous data frames to save space
remove(list = c("AIBQ", "Memory", "Qdat","AIBQdat", "Memdat"))

# remove participants with incomplete data # n=443 (with that extra participant, and DPT acc < 70% removed)
CCBH <- na.omit(CCBH)

# remove one copy participant
CCBH <- CCBH %>%
  filter(subject != 301041)

#mean(CCBH$MHC_Tot)
#sd(CCBH$MHC_Tot)
#quantile(CCBH$MHC_Tot, c(.33, .66)) # 37 and 47

CCBH2 <- CCBH %>%
  select(2,4,5,6,7,9,10) %>%
  rename(MH = MHC_Tot,
         IB_S_Pos = Pos_Soc,
         IB_N_Pos = Pos_nonsoc,
         IB_S_Neg = Neg_soc,
         IB_N_Neg = Neg_nonsoc,
         MB_Pos = MEM_pos,
         MB_Neg = MEM_neg)

# save correlation and covariance matrices
full_cor <- cor(CCBH2)
full_cov <- cov(CCBH2)

#write.csv(full_cor, "Apendices/S8_FullSample_cor.csv")
#write.csv(full_cov, "Apendices/S9_FullSample_cov.csv")

###

demographics2 <- demographics %>%
  filter(Master_subject %in% unique(CCBH$subject))

median(demographics2$SES, na.rm = TRUE)

ethincity_ns <- demographics2 %>%
  group_by(ethnicity_genetic) %>%
  summarise(n = n())


###########################################################################
###########################################################################

# attention bias

# specific recodings
AB_raw$subject <- ifelse(AB_raw$subject == "315177", "351177", AB_raw$subject)
AB_raw$subject <- ifelse(AB_raw$subject == "316384", "361384", AB_raw$subject)
AB_raw$subject <- ifelse(AB_raw$subject == "551153", "331153", AB_raw$subject)
AB_raw$subject <- ifelse(AB_raw$subject == "319519", "391519", AB_raw$subject)

# this for the subjects that were in group 371 not 361

AB_raw2 <- AB_raw %>%
  separate(subject, sep = 3, into = c("first", "second"))
AB_raw2$first <- ifelse(AB_raw2$first == "361" & as.numeric(AB_raw2$second) >= 334 & as.numeric(AB_raw2$second) <= 440, "371", AB_raw2$first)
AB_raw2$subject <- paste(AB_raw2$first, AB_raw2$second, sep = "")  

# finally, remove data from participants that had <70 accuracy and other reasons

AB_raw2 <-   AB_raw2 %>% 
  filter(subject != 311002,
         subject != 311019,
         subject != 311028,
         subject != 301036,
         subject != 301041,
         subject != 331126,
         subject != 381491,
         subject != 381441,
         subject != 381442,
         subject != 381443,
         subject != 381444,
         subject != 381478 )

AB_raw2 <- filter(AB_raw2, subject %in% unique(CCBH$subject))

# now, our participant list in the CCBH data will have the same subjects as the AB_raw2

# quick check
data.frame(name = unique(AB_raw2$subject) %in% unique(CCBH$subject),
           number = unique(AB_raw2$subject))

# need to remove those without compete dataframes so that reliability is calculated fro only the final sample


# ensuring that the blockcode only contains the actual condition
AB_raw2$blockcode <- ifelse(grepl("pain", AB_raw2$blockcode), "pain", 
                            ifelse(grepl("angry", AB_raw2$blockcode), "angry",
                                   ifelse(grepl("happy", AB_raw2$blockcode), "happy", "")))

# ensuring that the congruent and incongruent labels are of the correct name fir the script
AB_raw2$congruency <- ifelse(grepl("incongruent", AB_raw2$trialcode), "Incongruent", "Congruent")

# trimming as per protocol
AB_raw2 <- AB_raw2 %>%
  group_by(subject) %>%
  filter(blockcode != "") %>%
  filter(correct == 1) %>%
  filter(latency >= 200, latency <= 3000) %>%
  group_by(subject, blockcode, congruency) %>%
  mutate(high = mean(latency)+(3*sd(latency)),
         low = mean(latency)-(3*sd(latency))) %>%
  filter(latency >= low, latency <= high) %>%
  ungroup() %>%
  as.data.frame()

AB_reliability <- splithalf(data = AB_raw2, 
                            outcome = "RT",
                            score = "difference",
                            conditionlist = c("happy", "angry", "pain"),
                            halftype = "random",
                            permutations = 5000,
                            var.RT = "latency",
                            var.condition = "blockcode",
                            var.participant = "subject", 
                            var.trialnum = "trialnum",
                            var.compare = "congruency",
                            compare1 = "Congruent",
                            compare2 = "Incongruent")


AB_reliability$final_estimates[,-2:-5]


####  AIBQ

# some renaming first
# specific recodings
AIBQ_raw$subject <- ifelse(AIBQ_raw$subject == "315177", "351177", AIBQ_raw$subject)
AIBQ_raw$subject <- ifelse(AIBQ_raw$subject == "316384", "361384", AIBQ_raw$subject)
AIBQ_raw$subject <- ifelse(AIBQ_raw$subject == "551153", "331153", AIBQ_raw$subject)
AIBQ_raw$subject <- ifelse(AIBQ_raw$subject == "319519", "391519", AIBQ_raw$subject)

# this for the subjects that were in group 371 not 361

AIBQ_raw2 <- AIBQ_raw %>%
  separate(subject, sep = 3, into = c("first", "second"))
AIBQ_raw2$first <- ifelse(AIBQ_raw2$first == "361" & as.numeric(AIBQ_raw2$second) >= 334 & as.numeric(AIBQ_raw2$second) <= 440, "371", AIBQ_raw2$first)
AIBQ_raw2$subject <- paste(AIBQ_raw2$first, AIBQ_raw2$second, sep = "")  

# finally, remove data from participants that had <70 accuracy and other reasons

# AIBQ_raw2 <-   AIBQ_raw2 %>% 
#                 filter(subject != 311002,
#                        subject != 311019,
#                        subject != 311028,
#                        subject != 301036,
#                        subject != 301041,
#                        subject != 331126,
#                        subject != 381491,
#                        subject != 381441,
#                        subject != 381442,
#                        subject != 381443,
#                        subject != 381444,
#                        subject != 381478 )

AIBQ_raw2 <- filter(AIBQ_raw2, subject %in% unique(CCBH$subject))


# checking

unique(AIBQ_raw2$subject) %in% unique(CCBH$subject)

# social positive

soc_pos <- AIBQ_raw2 %>%
  select("sit2pos_response", "sit4pos_response", "sit7pos_response", "sit9pos_response", "sit10pos_response") %>%
  psych::alpha()

soc_pos_ci <- alpha.ci(soc_pos$total$raw_alpha,n.obs=nrow(AIBQ_raw2),n.var=5,p.val=.05,digits=2)

soc_pos_o <- AIBQ_raw2 %>%
  select("sit2pos_response", "sit4pos_response", "sit7pos_response", "sit9pos_response", "sit10pos_response") %>%
  psych::omega(plot = FALSE)

# social negative

soc_neg <- AIBQ_raw2 %>%
  select("sit2neg_response", "sit4neg_response", "sit7neg_response", "sit9neg_response", "sit10neg_response") %>%
  psych::alpha()

soc_neg_ci <- alpha.ci(soc_neg$total$raw_alpha,n.obs=nrow(AIBQ_raw2),n.var=5,p.val=.05,digits=2)

soc_neg_o <- AIBQ_raw2 %>%
  select("sit2neg_response", "sit4neg_response", "sit7neg_response", "sit9neg_response", "sit10neg_response") %>%
  psych::omega(plot = FALSE)

# nonsocial positive
# 1, 3, 5, 6, 8, 

nonsoc_pos <- AIBQ_raw2 %>%
  select("sit1pos_response", "sit3pos_response", "sit5pos_response", "sit6pos_response", "sit8pos_response") %>%
  psych::alpha()

nonsoc_pos_ci <- alpha.ci(nonsoc_pos$total$raw_alpha,n.obs=nrow(AIBQ_raw2),n.var=5,p.val=.05,digits=2)

nonsoc_pos_o <- AIBQ_raw2 %>%
  select("sit1pos_response", "sit3pos_response", "sit5pos_response", "sit6pos_response", "sit8pos_response") %>%
  psych::omega(plot = FALSE)

# nonsocial negative

nonsoc_neg <- AIBQ_raw2 %>%
  select("sit1neg_response", "sit3neg_response", "sit5neg_response", "sit6neg_response", "sit8neg_response") %>%
  psych::alpha()

nonsoc_neg_ci <- alpha.ci(nonsoc_neg$total$raw_alpha,n.obs=nrow(AIBQ_raw2),n.var=5,p.val=.05,digits=2)

nonsoc_neg_o <- AIBQ_raw2 %>%
  select("sit1neg_response", "sit3neg_response", "sit5neg_response", "sit6neg_response", "sit8neg_response") %>%
  psych::omega(plot = FALSE)

AIBQ_reliability_table <- data.frame(social = c("social","social","nonsocial","nonsocial"),
                                     valence = c("pos", "neg", "pos", "neg"),
                                     lower = c(round(soc_pos_ci$lower.ci,2),
                                               round(soc_neg_ci$lower.ci,2),
                                               round(nonsoc_pos_ci$lower.ci,2),
                                               round(nonsoc_neg_ci$lower.ci,2)),
                                     alpha = c(round(soc_pos$total$raw_alpha, 2),
                                               round(soc_neg$total$raw_alpha, 2),
                                               round(nonsoc_pos$total$raw_alpha, 2),
                                               round(nonsoc_neg$total$raw_alpha, 2)),
                                     upper = c(round(soc_pos_ci$upper.ci,2),
                                               round(soc_neg_ci$upper.ci,2),
                                               round(nonsoc_pos_ci$upper.ci,2),
                                               round(nonsoc_neg_ci$upper.ci,2)),
                                     omega = c(round(soc_pos_o$omega.tot,2),
                                               round(soc_neg_o$omega.tot,2),
                                               round(nonsoc_pos_o$omega.tot,2),
                                               round(nonsoc_neg_o$omega.tot,2)))

#### reliability for the SRET

MEM2 <- MEM_raw %>%
  select(-contains("recallresponse")) %>%
  select(Master_subject,
         contains("_endorsed"),
         contains("_recalled"))

for ( col in 1:ncol(MEM2)){
  colnames(MEM2)[col] <-  sub("response.", "", colnames(MEM2)[col])
}

colnames(MEM2)

MEM3 <- MEM2 %>%
  gather(key = "trial", value = "response", -Master_subject) %>%
  separate(trial, sep = "_", into = c("word", "EndorcedRecalled")) %>%
  spread(key = EndorcedRecalled, value = response) %>%
  mutate(recalled = replace_na(recalled, 0)) %>%
  mutate(EndorcedAndRecalled = ifelse(endorsed == 1 & recalled == 1, 1, 0)) %>%
  mutate(valence = ifelse(word == "afraid", "negative",
                   ifelse(word == "alone", "negative",
                   ifelse(word == "angry", "negative",
                   ifelse(word == "anxious", "negative",
                   ifelse(word == "attractive", "positive",
                   ifelse(word == "awful", "negative",
                   ifelse(word == "bad", "negative",
                   ifelse(word == "bold", "positive",
                   ifelse(word == "boring", "negative",
                   ifelse(word == "brave", "positive",
                   ifelse(word == "cheerful", "positive",
                   ifelse(word == "clever", "positive",
                   ifelse(word == "curious", "positive",
                   ifelse(word == "exciting", "positive",
                   ifelse(word == "foolish", "negative",
                   ifelse(word == "free", "positive",
                   ifelse(word == "friendly", "positive",
                   ifelse(word == "funny", "positive",
                   ifelse(word == "happy", "positive",
                   ifelse(word == "healthy", "positive",
                   ifelse(word == "helpful", "positive",
                   ifelse(word == "hurt", "negative",
                   ifelse(word == "interesting", "positive",
                   ifelse(word == "leader", "positive",
                   ifelse(word == "lively", "positive",
                   ifelse(word == "lonely", "negative",
                   ifelse(word == "loser", "negative",
                   ifelse(word == "lucky", "positive",
                   ifelse(word == "nervous", "negative",
                   ifelse(word == "nice", "positive",
                   ifelse(word == "nobody", "negative",
                   ifelse(word == "popular", "positive",
                   ifelse(word == "proud", "positive",
                   ifelse(word == "quiet", "negative",
                   ifelse(word == "sad", "negative",
                   ifelse(word == "scared", "negative",
                   ifelse(word == "sleepy", "negative",
                   ifelse(word == "smart", "positive",
                   ifelse(word == "strange", "negative",
                   ifelse(word == "terrible", "negative",
                   ifelse(word == "tired", "negative",
                   ifelse(word == "ugly", "negative",
                   ifelse(word == "unhappy", "negative",
                   ifelse(word == "winner", "positive", 
       NA)))))))))))))))))))))))))))))))))))))))))))))

MEM_reliability <- splithalf(data = MEM3, 
                             outcome = "accuracy",
                             score = "average",
                             average = "sum",
                             conditionlist = c("positive", "negative"),
                             var.condition = "valence",
                             var.ACC = "EndorcedAndRecalled",
                             var.participant = "Master_subject")



#####   reliability for the MHC

MHC_only <- select(scores2, c("Master_subject",
                              "MCSHCSQ001",
                              "MCSHCSQ002",
                              "MCSHCSQ003",
                              "MCSHCSQ004",
                              "MCSHCSQ005",
                              "MCSHCSQ006",
                              "MCSHCSQ007",
                              "MCSHCSQ008",
                              "MCSHCSQ009",
                              "MCSHCSQ010",
                              "MCSHCSQ011",
                              "MCSHCSQ012",
                              "MCSHCSQ013",
                              "MCSHCSQ014"))

MHC_only <- filter(MHC_only, Master_subject %in% unique(CCBH$subject))

# alphas for scales for each group
MHC_alpha <- MHC_only %>%
  filter(Master_subject %in% unique(CCBH$subject)) %>%
  select(-Master_subject) %>%
  psych::alpha()

MHC_alpha_ci <- alpha.ci(MHC_alpha$total$raw_alpha,n.obs=nrow(MHC_only),n.var=14,p.val=.05,digits=2)


MHC_omega <- MHC_only %>%
  filter(Master_subject %in% unique(CCBH$subject)) %>%
  select(-Master_subject) %>%
  psych::omega(plot = FALSE)



###########################################################################
###########################################################################

### create high and low groups


# checking sample size for the tertile split
low  <- sum(CCBH$MHC_Tot < 37) # 145
# mid  <- sum(CCBH$MHC_Tot >= 37 & CCBH$MHC_Tot <= 47) # 153
high <- sum(CCBH$MHC_Tot > 47) # 150

# subsetting into three samples
low2 <- subset(CCBH, MHC_Tot < 37) # 
# mid2 <- subset(CCBH, MHC_Tot >= 37 & MHC_Tot <= 47)
high2 <- subset(CCBH, MHC_Tot > 47)

low_ppt <- unique(low2$subject)
# mid_ppt <- unique(mid2$subject)
high_ppt <- unique(high2$subject)

# selecting only the variables for inclusion in the networks
low3  <- low2[,c(4,5,6,7,9,10)]
# mid3  <- mid2[,c(4,5,6,7,9,10)]
high3 <- high2[,c(4,5,6,7,9,10)]

# changing variable names for interpretability
variable_names  <- c("IB_S_Pos", "IB_N_Pos", "IB_S_Neg", "IB_N_Neg", "MB_Pos", "MB_Neg")
colnames(low3)  <- variable_names
# colnames(mid3)  <- variable_names
colnames(high3) <- variable_names


# saving correlation and covariance matrices for the apendices. 

# low
low_cor <- cor(low3)
low_cov <- cov(low3)
high_cor <- cor(high3)
high_cov <- cov(high3)

#write.csv(low_cor, "Apendices/S2_lowMH_cor_matrix.csv")
#write.csv(low_cov, "Apendices/S3_lowMH_cov_matrix.csv")
#write.csv(high_cor, "Apendices/S5_highMH_cor_matrix.csv")
#write.csv(high_cov, "Apendices/S6_highMH_cov_matrix.csv")


###########################################################################
###########################################################################

### table1, results = 'asis', include = TRUE

m_sd <- psych::describe(CCBH2) %>%
  as.data.frame() %>%
  dplyr::select(mean, sd) %>%
  round(2)

rownames(m_sd) <- c("MH (1)",
                    "IB_S_Pos (2)",
                    "IB_N_Pos (3)",
                    "IB_S_Neg (4)",
                    "IB_N_Neg (5)",
                    "MB_Pos (6)",
                    "MB_Neg (7)")

fash1 <- fashion(cor(CCBH2), decimals = 2, leading_zeros = FALSE, na_print = "")

class(fash1) <- "data.frame"

fash1 <- fash1 %>%
  rename('(1)' = MH,
         '(2)' = IB_S_Pos,
         '(3)' = IB_N_Pos,
         '(4)' = IB_S_Neg,
         '(5)' = IB_N_Neg,
         '(6)' = MB_Pos,
         '(7)' = MB_Neg) %>%
  select(-'(7)')

fash1[upper.tri(fash1, diag = TRUE)] <- ""


table1 <- as.data.frame(cbind(m_sd, fash1))


papaja::apa_table(table1,
                  align = c('l', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r'),
                  caption = "Descriptive statistics and correlation matrix for the full sample (n = 450)",
                  note = "MH = Positive mental health; IB_S_Pos = Social Positive Interpretation Bias; IB_S_Neg = Social Negative Interpretation Bias; IB_N_Pos = Non-Social Negative Interpretation Bias; IB_N_Neg = Non-Social; MB_Pos = positive memory bias; MB_Neg = negative memory bias.")

# write.csv(table1, "Apendices/Table_1.csv")


###########################################################################
###########################################################################

### glasso networks


# estimate the networks
glasso_low  <- estimateNetwork(low3, default = "EBICglasso", tuning = .5)
#glasso_mid  <- estimateNetwork(mid3, default = "EBICglasso", tuning = .5)
glasso_high  <- estimateNetwork(high3, default = "EBICglasso", tuning = .5)

low_glas_weights <- glasso_low$graph
colnames(low_glas_weights) <- variable_names
rownames(low_glas_weights) <- variable_names

#mid_glas_weights <- glasso_mid$graph
#colnames(mid_glas_weights) <- variable_names
#rownames(mid_glas_weights) <- variable_names

high_glas_weights <- glasso_high$graph
colnames(high_glas_weights) <- variable_names
rownames(high_glas_weights) <- variable_names

max_glass <- max(low_glas_weights[which(low_glas_weights != 1)], 
                 #mid_glas_weights[which(mid_glas_weights != 1)],
                 high_glas_weights[which(high_glas_weights != 1)]) # returns highest weight in the three samples' correlation matrices - to be used


#write.csv(low_glas_weights, "./Apendices/S4_low_glasso_weights.csv")
#write.csv(mid_glas_weights, "./Apendices/mid_glasso_weights.csv")
#write.csv(high_glas_weights, "./Apendices/S7_high_glasso_weights.csv")

##### plot the three networks

layout_common <- averageLayout(glasso_low,
                               #glasso_mid,
                               glasso_high) # this is now the common layout for each figure

# plot the low and high sample networks
# save this plot

# Cairo::Cairo(file="Figures/twoglasso.png", 
#              type="png",
#              bg = "white",
#              units="in", 
#              width=16, 
#              height=8, 
#              pointsize=12, 
#              dpi=144)

layout(t(1:2))

plot(glasso_low, 
     layout = layout_common, 
     title = "low - MH plot", 
     theme = "colorblind", 
     labels = variable_names, 
     vsize = 13, 
     maximum = max_glass,
     title.cex = 2,
     label.cex = 1.1,
     border.width = 1.8)

plot(glasso_high, 
     layout = layout_common, 
     title = "high - MH plot", 
     theme = "colorblind", 
     labels = variable_names, 
     vsize = 13, 
     maximum = max_glass,
     title.cex = 2,
     label.cex = 1.1,
     border.width = 1.8)

#dev.off()

layout(1)

###########################################################################
###########################################################################

### networkcomparisontest


lh <- NCT(low3[,1:6], high3[,1:6], it=1000, gamma = .5, binary.data=FALSE)

# could also run 
#lm <- NCT(low3[,1:6], mid3[,1:6], it=1000, gamma = .5, binary.data=FALSE)
#mh <- NCT(mid3[,1:6], high3[,1:6], it=1000, gamma = .5, binary.data=FALSE)

lh$glstrinv.pval
lh$nwinv.pval

#lm$glstrinv.pval
#lm$nwinv.pval
#mh$glstrinv.pval
#mh$nwinv.pval


###########################################################################
###########################################################################

### check variance differences


low4 <- low3
low4$group <- "low"

high4 <- high3
high4$group <- "high"

mix <- rbind(low4, high4)
mix$group <- as.factor(mix$group)

bftest_IBSPos <- onewaytests::bf.test(data = mix,
                                      IB_S_Pos ~ group)
bftest_IBNPos <- onewaytests::bf.test(data = mix,
                                      IB_N_Pos ~ group)
bftest_IBSNeg <- onewaytests::bf.test(data = mix,
                                      IB_S_Neg ~ group)
bftest_IBNNeg <- onewaytests::bf.test(data = mix,
                                      IB_N_Neg ~ group)

bftest_MBPos <- onewaytests::bf.test(data = mix,
                                     MB_Pos ~ group)
bftest_MBNeg <- onewaytests::bf.test(data = mix,
                                     MB_Neg ~ group)

if(all(bftest_IBSPos$p.value < .05,
       bftest_IBNPos$p.value < .05,
       bftest_IBSNeg$p.value < .05,
       bftest_IBNNeg$p.value < .05,
       bftest_MBPos$p.value < .05,
       bftest_MBNeg$p.value < .05)) {
  bftest_min <- "< .05"
}



var(low4$IB_S_Pos)
var(high4$IB_S_Pos)

var(low4$IB_N_Pos)
var(high4$IB_N_Pos)

var(low4$IB_S_Neg)
var(high4$IB_S_Neg)

var(low4$IB_S_Neg)
var(high4$IB_S_Neg)

var(low4$MB_Pos)
var(high4$MB_Pos)

var(low4$MB_Neg)
var(high4$MB_Neg)

low_m_sd <- psych::describe(low2[,c(2,4,5,6,7,9,10)]) %>%
  as.data.frame() %>%
  select(mean, sd) %>%
  round(2)

low_m_sd$var <- round(sapply(low2[,c(2,4,5,6,7,9,10)], var), 2)

low_m_sd <- low_m_sd %>%
  rename(`low M` = mean,
         `low SD` = sd,
         `low Var` = var)

high_m_sd <- psych::describe(high2[,c(2,4,5,6,7,9,10)]) %>%
  as.data.frame() %>%
  select(mean, sd) %>%
  round(2)

high_m_sd$var <- round(sapply(high2[,c(2,4,5,6,7,9,10)], var), 2)

high_m_sd <- high_m_sd %>%
  rename(`high M` = mean,
         `high SD` = sd,
         `high Var` = var)

table2 <- cbind(low_m_sd, high_m_sd)




###########################################################################
###########################################################################


### r moderated_network

fit_obj3 <- mgm(data = CCBH2,
                type = c(rep("g", 7)),
                level = c(rep(1, 7)),
                ruleReg = "OR",
                k = 2,
                binarySign = TRUE,
                moderators = c(1)
)

# n interactions
n_int3 <- length(fit_obj3$interactions$weightsAgg[[2]])

# for the predictability aspect
p_obj3 <- predict(fit_obj3, CCBH2,
                  errorCat = c("CC","nCC","CCmarg"),
                  errorCon = c("R2"))

error_list3 <- list() # List for ring-segments
for(i in 1:7) error_list3[[i]] <- p_obj3$errors[i,2]
#beyondmarg <- p_obj2$errors[10,3]-p_obj2$errors[10,5]
#error_list[[10]] <- c(p_obj2$errors[10,5],beyondmarg)

color_list3 <- list() # List for Colors
for(i in 1:7) color_list3[[i]] <- "#90B4D4"
#color_list[[10]] <- c("#ffa500", "#ff4300")

# saving these for use in the text
pred_list <- data.frame(var = colnames(CCBH2),
                        pred = unlist(error_list3))

prediction_min <- pred_list[which.min(pred_list$pred),"pred"]
prediction_max <- pred_list[which.max(pred_list$pred),"pred"]
prediction_MH <- pred_list[pred_list$var == "MH","pred"]

mgm_figure <- qgraph(fit_obj3$pairwise$wadj * fit_obj3$pairwise$signs,
                     pie = error_list3,
                     layout="spring",
                     labels = colnames(CCBH2),
                     pieColor = color_list3,
                     curveAll = TRUE, 
                     curveDefault = .6,
                     cut = 0, 
                     labels = colnames(CCBH2),
                     theme = "colorblind",
                     vsize = 12,
                     label.cex = 1,
                     title = "")

###########################################################################
###########################################################################

### spinglass

mgm_igraph <- as.igraph(mgm_figure, attributes = TRUE)

spinglass <- data.frame(MH = NULL,
                        IB_S_Pos = NULL,
                        IB_N_Pos = NULL,
                        IB_S_Neg = NULL,
                        IB_N_Neg = NULL,
                        MB_Pos = NULL,
                        MB_Neg = NULL)
for(i in 1:1000) {
  sgc <- spinglass.community(mgm_igraph)
  spinglass[i,1:7] <- sgc$membership[1:7]
}

spinglass <- rename(spinglass,
                    MH = V1,
                    IB_S_Pos = V2,
                    IB_N_Pos = V3,
                    IB_S_Neg = V4,
                    IB_N_Neg = V5,
                    MB_Pos = V6,
                    MB_Neg = V7)
spinglass$test <- ifelse(spinglass$IB_S_Pos == 1 & spinglass$IB_N_Pos == 1 & 
                           spinglass$IB_S_Neg == 2 & spinglass$IB_N_Neg == 2 & 
                           spinglass$MB_Pos == 1 & spinglass$MB_Neg == 2, 1,
                         ifelse(spinglass$IB_S_Pos == 2 & spinglass$IB_N_Pos == 2 &
                                  spinglass$IB_S_Neg == 1 & spinglass$IB_N_Neg == 1 & 
                                  spinglass$MB_Pos == 2 & spinglass$MB_Neg == 1, 1, 0))

n_pos_neg <- sum(spinglass$test)/nrow(spinglass) # proportion of runs splitting into positive and negative

spinglass2 <- subset(spinglass, test == 1)
spinglass2$MHPos <- ifelse(spinglass2$MH == spinglass2$IB_S_Pos, 1, 0)
n_MH_pos <- sum(spinglass2$MHPos)/nrow(spinglass2) # number of runs that mental health falls into the positive group.


qgraph(fit_obj3$pairwise$wadj * fit_obj3$pairwise$signs,
       pie = error_list3,
       layout="spring",
       labels = colnames(CCBH2),
       pieColor = color_list3,
       curveAll = TRUE, 
       curveDefault = .6,
       cut = 0, 
       labels = colnames(CCBH2),
       theme = "colorblind",
       vsize = 12,
       label.cex = 1
)
dev.off()


### resampling_mgm
n_resamples <- 1000

res_obj <- resample(object = fit_obj3,
                    data = CCBH2,
                    nB = n_resamples)


###########################################################################
###########################################################################

### code to visualise the moderated networks

# Adaptation of the plotRes code (from the mgm package) to extract the label order and to adapt the plotting 
# this chunk extracts the edge weights

quantiles <- c(.05, .95)
labels = NULL
decreasing = TRUE
cut = NULL
cex.label = .75
lwd.qtl = 2
cex.mean = .5 
cex.bg = 3.5 
axis.ticks = c(-.5, -.25, 0, .25, .5, .75, 1)
labels = colnames(CCBH2)


# Get basic info
dims <- dim(res_obj$bootParameters)
p <- dims[1]
nB <- dims[3]
n_pars <- p*(p-1) / 2

# Collapse into edge x property matrix
tar_mat <- matrix(NA, nrow=n_pars, ncol = 6)
colnames(tar_mat) <- c("Variable A", "Variable B", "Mean", "qtl_low", "qtl_high", "propLtZ")

counter <- 1
for(row in 1:p) {
  for(col in row:p) {
    if(row!=col){
      
      # Variable ids
      tar_mat[counter, 1] <- row
      tar_mat[counter, 2] <- col
      
      # Quantiles
      qtls <- quantile(res_obj$bootParameters[row, col, ], probs = quantiles)
      tar_mat[counter, 3] <- mean(res_obj$bootParameters[row, col, ])
      tar_mat[counter, 4] <- qtls[1]
      tar_mat[counter, 5] <- qtls[2]
      tar_mat[counter, 6] <- mean(abs(res_obj$bootParameters[row, col, ]) > 0) # proportion estimates > 0
      
      # update counter
      counter <- counter + 1
    }
  }
}


# Order
tar_mat <- tar_mat[order(tar_mat[,3], decreasing = decreasing), ]

# Subset (cut)
if(is.null(cut)) {
  TM <- tar_mat
} else {
  TM <- tar_mat[cut, ]
}


# Generate label vector
if(is.null(labels)) {
  label_vec <- paste0(TM[, 1], " - ", TM[, 2])
} else {
  tar_mat_label <- TM[ ,1:2]
  tar_mat_label <- apply(tar_mat_label, 1:2, as.character)
  for(i in 1:p) tar_mat_label[tar_mat_label == i] <- labels[i]
  label_vec <- paste0(tar_mat_label[, 1], " - ", tar_mat_label[, 2])
}

## this chunk extracts the interaction effects

call <- list('object' = "object",
             'data' = "data",
             'nB' = "nB",
             'blocks' = "blocks",
             'pbar' = "pbar")

outlist <- list('call' = call,
                'bootParameters' = NULL,
                'bootQuantiles' = NULL,
                'models' = NULL,
                "Times" = rep(NA, nB),
                "totalTime" = NULL)

outlist$models <- res_obj$models

p <- length(fit_obj3$call$type)
# nquantiles <- length(quantiles)
nB <- n_resamples    

## Collect all estimates
collect_array <- collect_array_sign <- array(0, dim = c(p, p, nB))

#####
quantiles <- c(0.05, .95)
nquantiles <- length(quantiles)


for(b in 1:nB) {
  
  for(i in 1:length(outlist$models[[b]]$interactions$indicator[[2]][,2])){
    
    collect_array[as.numeric(outlist$models[[b]]$interactions$indicator[[2]][i,2]),
                  as.numeric(outlist$models[[b]]$interactions$indicator[[2]][i,3]),
                  b] <- as.numeric(outlist$models[[b]]$interactions$weightsAgg[[2]][i])
    
  }
  
}

for(b in 1:nB) {
  
  for(i in 1:length(outlist$models[[b]]$interactions$indicator[[2]][,2])){
    
    collect_array_sign[as.numeric(outlist$models[[b]]$interactions$indicator[[2]][i,2]),
                       as.numeric(outlist$models[[b]]$interactions$indicator[[2]][i,3]),
                       b] <- as.numeric(outlist$models[[b]]$interactions$signs[[2]][i])
    
  }
  
}     

# add sign
collect_array_wS <- collect_array
ind_negative <- which(collect_array_sign == -1, arr.ind = TRUE)
collect_array_wS[ind_negative] <- collect_array_wS[ind_negative] * -1



# Compute quantiles
quantile_array <- apply(collect_array_wS, 1:2, function(x) quantile(x, probs = quantiles))
quantile_array_res <- array(dim = c(p, p, nquantiles))
for(qu in 1:nquantiles) quantile_array_res[, , qu] <- quantile_array[qu, , ]

outlist$bootParameters <- collect_array_wS
outlist$bootQuantiles <- quantile_array_res

# this chunk extracts the interaction effects and reorders to the same as the edge weights

quantiles2 <- c(.05, .95)
labels2 = NULL
decreasing2 = TRUE
cut2 = NULL
cex.label2 = .75
lwd.qtl2 = 2
cex.mean2 = .5 
cex.bg2 = 3.5 
axis.ticks2 = c(-.2, 0, .2)
labels2 = colnames(CCBH2)

# Get basic info
dims2 <- dim(outlist$bootParameters)
p2 <- dims2[1]
nB2 <- dims2[3]
n_pars2 <- p2*(p2-1) / 2

# Collapse into edge x property matrix
tar_mat2 <- matrix(NA, nrow=n_pars2, ncol = 6)
colnames(tar_mat2) <- c("Variable A", "Variable B", "Mean", "qtl_low", "qtl_high", "propLtZ")

counter2 <- 1
for(row2 in 1:p2) {
  for(col2 in row2:p2) {
    if(row2!=col2){
      
      # Variable ids
      tar_mat2[counter2, 1] <- row2
      tar_mat2[counter2, 2] <- col2
      
      # Quantiles
      qtls2 <- quantile(outlist$bootParameters[row2, col2, ], probs = quantiles2)
      tar_mat2[counter2, 3] <- mean(outlist$bootParameters[row2, col2, ])
      tar_mat2[counter2, 4] <- qtls2[1]
      tar_mat2[counter2, 5] <- qtls2[2]
      tar_mat2[counter2, 6] <- mean(abs(outlist$bootParameters[row2, col2, ]) > 0) # proportion estimates > 0
      
      # update counter
      counter2 <- counter2 + 1
    }
  }
}

# this is where the change needs to be, the aim is to have the same order as in the first figure.

# Order
tar_mat2 <- tar_mat2[order(match(paste(tar_mat2[,1],tar_mat2[,2]), 
                                 paste(tar_mat[,1],tar_mat[,2]))
),]




# Subset (cut)
if(is.null(cut2)) {
  TM2 <- tar_mat2
} else {
  TM2 <- tar_mat2[cut2, ]
}

# extract data from the previous 2 figures

out_table <- data.frame(edge = label_vec,
                        V1 = tar_mat[,1],
                        V2 = tar_mat[,2],
                        edge_mean = tar_mat[,3],
                        edge_low95 = tar_mat[,4],
                        edge_high95 = tar_mat[,5],
                        edge_prop = tar_mat[,6],
                        int_mean = tar_mat2[,3],
                        int_low95 = tar_mat2[,4],
                        int_high95 = tar_mat2[,5],
                        int_prop = tar_mat2[,6])

levelorder <- out_table$edge

#write.csv(out_table, "Apendices/S10_FullSample_mgm_EdgeAndInteraction.csv")

over80 <- "edges"

for(i in 1:nrow(out_table)) {
  if(out_table[i,"int_prop"] > .8) {
    over80 <- paste(over80, "; and " , out_table[i,"edge"], sep = "")
  }
}

out_table %>%
  select(edge, int_prop) %>%
  subset(int_prop > .8)


###########################################################################
###########################################################################

cols <- c("#BF0000", "#E7A0A0FF", "#FFFFFF", "#ACACF1FF", "#0000D5")

# edges
edgesplot <- ggplot(data = out_table,
                    aes(
                      x = reorder(edge, edge_mean),
                      y = edge_mean,
                      group = 1,
                      level = levelorder
                    )) +
  geom_ribbon(
    data = out_table,
    aes(ymin = edge_low95, ymax = edge_high95, group = 1),
    fill = "grey80",
    alpha = .5
  ) +
  labs(y = "edge strength") +
  coord_flip() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    size = .5
  ) +  geom_line(aes(colour = edge_mean), size = 1.5) +
  geom_point(aes(colour = edge_mean), size = 3) +
  scale_colour_gradientn(
    colours = cols,
    values = rescale(c(-.2,-.03, 0, .03, .2)),
    limits = c(-.5, .5),
    guide = "none"
  ) +
  
  theme(axis.title.y = element_blank()) +
  geom_point(
    colour = "white",
    size = 7,
    shape = "circle",
    y = .55
  ) +
  geom_text(aes(label = sub(
    "^(-?)0.", "\\1.", sprintf("%.2f", edge_prop)
  ), y = .55), size = 3)

# interactions

intplot <- ggplot(data = out_table,
                  aes(
                    x = reorder(edge, edge_mean),
                    y = int_mean,
                    group = 1,
                    level = levelorder,
                    label = int_prop
                  )) +
  geom_errorbar(data = out_table, aes(ymin = int_low95, ymax = int_high95, group = 1)) +
  labs(y = "interaction strength") +
  ylim(min = -.16, max = .15) +
  coord_flip() +
  scale_colour_gradientn(
    colours = muted(cols, l = 50, c = 70),
    values = rescale(c(-.1,-.01, 0, .01, .1)),
    limits = c(-.2, .2),
    guide = "none"
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    size = .5
  ) +
  geom_point() +
  geom_point(colour = "white",
             size = 7,
             shape = "circle", y = .15) +
  geom_text(aes(label = sub(
    "^(-?)0.", "\\1.", sprintf("%.2f", int_prop)
  )),
  size = 3, y = .15)

grid.arrange(edgesplot, intplot, ncol = 2, widths = c(.6, .4))


###########################################################################
###########################################################################

### parametric bootstrap

# first, extract the pairwise model

pair_model <-  fit_obj3$pairwise$wadj * fit_obj3$pairwise$signs
pair_model[is.na(pair_model)] <- 0

# save the ggmgenerator into another name, I'm still not sure why this is necessary
gendat <- ggmGenerator()

genList <- NULL
edgeList <- NULL
int_List <- NULL

sims <- 1000

# this will simualte 1000 datasets using the null model (only pairwise associations), then attempt to fit the moderated model to them
for(i in 1:sims) {
  gen <- gendat(nrow(CCBH2), pair_model)
  
  fitMGM <- mgm(data = gen,
                type = c(rep("g", 7)),
                level = c(rep(1, 7)),
                ruleReg = "OR",
                k = 2,
                binarySign = TRUE,
                moderators = c(1))
  
  edgeList_temp <- fitMGM$pairwise$wadj * fitMGM$pairwise$signs
  
  edgeList_temp[is.na(edgeList_temp)] <- 0
  
  #genList[[i]] <- genList_temp      
  edgeList[[i]] <- edgeList_temp
  
  # interactions
  
  int_List_temp <- matrix(rep(0,7*7),nrow = 7)
  
  if(length(fitMGM$interactions$weightsAgg[[2]]) > 0) {
    
    jlist <- 1:length(fitMGM$interactions$weightsAgg[[2]])
    
    for(j in jlist){
      
      indicator <- data.frame(t(fitMGM$interactions$indicator[[2]]))
      
      row <- indicator[j,2]
      col <- indicator[j,3]
      
      int_List_temp[row, col] <- as.numeric(fitMGM$interactions$weightsAgg[[2]][j]) * fitMGM$interactions$signs[[2]][j] 
      
    }
  }
  
  int_List[[i]] <- int_List_temp
  
  print(i)
}


quantiles <- c(.05, .95)


testing <- data.frame(row = rep(NA, 21),
                      col = rep(NA, 21),
                      edge_mean = rep(NA, 21),
                      edge_qtl_low = rep(NA, 21), 
                      edge_qtl_high = rep(NA, 21),
                      edge_proportion = rep(NA, 21),
                      int_mean = rep(NA, 21),
                      int_qtl_low = rep(NA, 21), 
                      int_qtl_high = rep(NA, 21),
                      int_proportion = rep(NA, 21))
count <- 1

for(row in 1:p) {
  for(col in row:p) {
    if(row!=col){
      
      edge_temp <- sapply(edgeList, function(x) x[row, col])
      
      testing[count, "row"] <- row
      testing[count, "col"] <- col
      testing[count, "edge_mean"] <- mean(edge_temp)
      
      testing[count,"edge_proportion"] <- sum(edge_temp != 0) / sims
      
      
      edge_qtls <- quantile(edge_temp, prob = quantiles)
      
      testing[count, "edge_qtl_low"] <- edge_qtls[1]
      testing[count, "edge_qtl_high"] <- edge_qtls[2]
      
      
      int_temp <- sapply(int_List, function(x) x[row, col])
      
      testing[count, "int_mean"] <- mean(int_temp)
      
      int_qtls <- quantile(int_temp, prob = quantiles)
      
      testing[count, "int_qtl_low"] <- int_qtls[1]
      testing[count, "int_qtl_high"] <- int_qtls[2]
      
      testing[count,"int_proportion"] <- sum(int_temp != 0) / sims
      
      count <- count + 1
      
    }
  }
}

testing <- testing[order(testing[,"edge_mean"], decreasing = decreasing), ]

testing$row <- recode(testing$row, 
                      `1` = "MH",
                      `2` = "IB_S_Pos",
                      `3` = "IB_N_Pos",
                      `4` = "IB_S_Neg",
                      `5` = "IB_N_Neg",
                      `6` = "MB_Pos",
                      `7` = "MB_Neg")
testing$col <- recode(testing$col, 
                      `1` = "MH",
                      `2` = "IB_S_Pos",
                      `3` = "IB_N_Pos",
                      `4` = "IB_S_Neg",
                      `5` = "IB_N_Neg",
                      `6` = "MB_Pos",
                      `7` = "MB_Neg")


testing$edge <- paste0(testing$row, "-", testing$col)

testing <- testing %>%
  mutate(edge = fct_reorder(edge, edge_mean))

cols <- c("#BF0000", "#E7A0A0FF", "#FFFFFF", "#ACACF1FF", "#0000D5")


eff_plot <- testing %>%
  ggplot(aes(x = edge,
             y = edge_mean,
             group = 1)) +
  geom_ribbon(aes(ymin = edge_qtl_low, ymax = edge_qtl_high),
              fill = "grey80",
              alpha = .5,
              group = 1) +
  coord_flip() +
  labs(y = "edge strength") +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    size = .5
  ) +  
  geom_line(aes(colour = edge_mean), size = 1.5) +
  geom_point(aes(colour = edge_mean), size = 3) +
  theme(axis.title.y = element_blank()) +
  scale_colour_gradientn(
    colours = cols,
    values = rescale(c(-.2,-.03, 0, .03, .2)),
    limits = c(-.5, .5),
    guide = "none"
  ) +
  theme(axis.title.y = element_blank()) +
  geom_point(
    colour = "white",
    size = 7,
    shape = "circle",
    y = .55
  ) +
  geom_text(aes(label = sub(
    "^(-?)0.", "\\1.", sprintf("%.2f", edge_proportion)
  ), y = .55), size = 3)



int_plot <- testing %>%
  ggplot(aes(x = edge,
             y = int_mean)) +
  labs(y = "interaction strength") +
  ylim(min = -.16, max = .15) +
  geom_errorbar(aes(ymin = int_qtl_low, ymax = int_qtl_high)) +
  coord_flip() +
  scale_colour_gradientn(
    colours = muted(cols, l = 50, c = 70),
    values = rescale(c(-.1,-.01, 0, .01, .1)),
    limits = c(-.2, .2),
    guide = "none"
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    size = .5
  ) +
  geom_point() +
  geom_point(colour = "white",
             size = 7,
             shape = "circle",
             y = .15) +
  geom_text(aes(label = sub(
    "^(-?)0.", "\\1.", sprintf("%.2f", int_proportion)
  )),
  y = .15, size = 3)


grid.arrange(eff_plot, int_plot, ncol = 2, widths = c(.6, .4))

sum((testing[,"int_qtl_low"] < out_table[,"int_mean"]) &
      (out_table[,"int_mean"] < testing[,"int_qtl_high"])) # 0

#write.csv(testing, "Apendices/S11_ParametricBoostrap_mgm_EdgeAndInteraction.csv")

###########################################################################
###########################################################################


### Conditioning the network on Mental Health

# hist(scale(CCBH2$MH))

mgm_minus1 <- mgm::condition(object = fit_obj3,
                             values = list("1" = -1))

mgm_ave <- mgm::condition(object = fit_obj3,
                          values = list("1" = 0))

mgm_plus1 <- mgm::condition(object = fit_obj3,
                            values = list("1" = 1))


max_compare = max(c(
  mgm_minus1$pairwise$wadj,
  mgm_ave$pairwise$wadj,
  mgm_plus1$pairwise$wadj
))

layout(t(1:3))

# note: the edgelist must be multiplied by the signs list, else all edges are positive

qgraph(mgm_minus1$pairwise$wadj * mgm_minus1$pairwise$signs,
       directed = FALSE, 
       labels = colnames(CCBH2), 
       theme = "colorblind", 
       maximum = max_compare,
       vsize = 17,
       title = "-1SD MH",
       title.cex = 1.5,
       label.cex = 1.1,
       border.width = 2)

qgraph(mgm_ave$pairwise$wadj * mgm_ave$pairwise$signs,
       directed = FALSE, 
       labels = colnames(CCBH2), 
       theme = "colorblind", 
       maximum = max_compare,
       vsize = 17,
       title = "mean MH",
       title.cex = 1.5,
       label.cex = 1.1,
       border.width = 2)

qgraph(mgm_plus1$pairwise$wadj * mgm_plus1$pairwise$signs,
       directed = FALSE, 
       labels = colnames(CCBH2), 
       theme = "colorblind", 
       maximum = max_compare,
       vsize = 17,
       title = "+1SD MH",
       title.cex = 1.5,
       label.cex = 1.1,
       border.width = 2)
