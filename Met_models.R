###### Read in methylation data

mydata <- read.table("combined_27_450_methylation", header = T)
rownames(mydata) <- mydata[, 1]
mydata <- mydata[, -1]

# Calculate the number of missing observations for each site
M = apply(is.na(mydata),2,sum)

## Remove columns and rows with more than 5% NAs
rmNaCol <- function (x) {
  flag = FALSE
  num = 0
  for (i in x) {
    if (is.na(i)) {
      num = num + 1
    }
  }
  if (num > length(x) * 0.05) {
    flag = TRUE
  }
  return (flag)
}

flag_col <- apply(mydata, 2, rmNaCol)
df <- mydata[,flag_col==F]

# Impute missing values with the median values of the site
f=function(x){
	x<-as.numeric(as.character(x))
    x[is.na(x)] = median(x, na.rm=TRUE) #convert the item with NA to median value from the column
    x
}
df = data.frame(apply(df,2,f))
df = data.frame(rownames(mydata), df)
############################
# Read in the survival data
############################
clinical <- read.table("nationwidechildrens.org_clinical_patient_gbm.txt", 
                       header=T, sep="\t", na.strings=c("[Not Applicable]", "[Not Available]"), stringsAsFactors=F)
clinical <- clinical[c(-1,-2),]
# Extract the barcode, days to last followup and days to death
survival <- data.frame(clinical$bcr_patient_barcode, 
                       as.numeric(as.character(clinical$last_contact_days_to)),
                       as.numeric(as.character(clinical$death_days_to)),
                       as.character(clinical$vital_status),
                       as.numeric(as.character(clinical$birth_days_to)),
                       as.factor(clinical$gender))
colnames(survival) <- c("barcode", "days_to_last_followup", "days_to_death", "vital_status", "age", "gender")

# Trim off the patient barcode
patient = NULL
for (i in survival[,1]) {
  x = unlist(strsplit(i, "-"))
  x = as.character(x[3])
  patient = c(patient, x)
}
survival <- cbind(patient, survival)
survival_time <- apply(survival[,3:4], 1, max, na.rm=T) # choose the larger of days to last follow up and days to death
survival <- data.frame(survival[,1:4], survival_time, survival[,5:7])

# Intersect expression data with survival data 
survival_met <- survival[survival[,1] %in% rownames(mydata), ]

# Include both Dead and Alive LTS patients
cutoff = 1095
#cutoff = 1615
met_dead <- survival_met[survival_met[,6]=="Dead",]
met_alive_lts <- survival_met[survival_met[,6]=="Alive" & survival_met[,5] > cutoff,]
met_lts <- rbind(met_dead, met_alive_lts)

# Classify LTS patients
label <- ifelse(met_lts[,5] > cutoff, 1, 0)
met_lts <- na.omit(data.frame(met_lts[,-2:-6], label))
datatable <- merge(met_lts, df, by.x = 1, by.y = 1)

######### Univariate Logistic Regression #################################################

log_regression <- function(data, start_num) {
  coeff = NULL
  genes = NULL
  ngenes = dim(data)[2] - start_num + 1
  for (i in start_num:dim(data)[2]) {
    testdata <- data[, c(2,i)]
    # leave out the mutations that does not have any calls in any patient
    if (length(testdata[testdata[,2]==0,2]) == length(testdata[,2])) {
      next
    }
    colnames(testdata) <- c("LTS", "Pheno")
    model <- glm(LTS ~ Pheno, family = binomial, data = testdata)
    genes <- c(genes, colnames(data)[i])
    coeff <- rbind(coeff, as.numeric(coef(summary(model))[2,]))
  }
  # Bonferroni Correction
  pbonf <- p.adjust(coeff[,4], method = "bonferroni", n = length(coeff[,4]))

  # Benjamini-Hochberg correction
  pbh <- p.adjust(coeff[,4], method = "BH", n = length(coeff[,4]))
  coeff <- data.frame(genes, coeff, pbonf, pbh) 
  colnames(coeff) <- c("Gene", "Estimate", "Std_Error", "z_value", "p_value",
                       "pbonf", "pBH")
  return (coeff)
}

y = datatable[, -2:-3]
coeff <- log_regression(y, 3)
# write.table(coeff, "Association_LTS_methylation.txt", 
            # sep = "\t", quote  = F, col.names = T, row.names = F)
write.table(coeff, "Association_LTS_methylation_1615.txt", 
            sep = "\t", quote  = F, col.names = T, row.names = F)
sigBonf <- coeff[coeff[,6] < 0.05,]  # Bonferroni correction
sigBH <- coeff[coeff[,7] < 0.05,] # Benjamini-Hochberg correction

# write.table(sigBonf, "SigAssociation_LTS_methylation_Bonf.05.txt", 
            # sep = "\t", quote  = F, col.names = T, row.names = F)
# write.table(sigBH, "SigAssociation_LTS_methylation_BH.05.txt", 
            # sep = "\t", quote  = F, col.names = T, row.names = F)
			
write.table(sigBonf, "SigAssociation_LTS_methylation_Bonf.05_1615d.txt", 
            sep = "\t", quote  = F, col.names = T, row.names = F)
write.table(sigBH, "SigAssociation_LTS_methylation_BH.05_1615d.txt", 
            sep = "\t", quote  = F, col.names = T, row.names = F)

#########################LASSO Logistic Regression########################################
library(glmnet)
# LLR on whole set
x = as.matrix(datatable[, -1:-4])
y = datatable[, 4]
fit=glmnet(x,y,family="binomial")
auc = NULL
for (i in 1:100) {
cvfit=cv.glmnet(x,y,family="binomial",type.measure="auc", nfold = 10)
cvm = cvfit$cvm[ which(cvfit$lambda == cvfit$lambda.min) ]
#cvfit$cvsd[ which(cvfit$lambda == cvfit$lambda.min) ]
auc = c(auc, cvm)
}
mean(auc)
sd(auc)

small.lambda.index <- which(cvfit$lambda == cvfit$lambda.min)
small.lambda.betas <- cvfit$glmnet.fit$beta[,small.lambda.index]
betas <- small.lambda.betas[small.lambda.betas != 0]
length(betas)
betas <- data.frame(betas)
# write.table(betas, "LASSO_regression_methylation_GBM.txt", sep = "\t", quote = F, row.names=T, col.names=F)

write.table(betas, "LASSO_regression_methylation_GBM_1615d.txt", sep = "\t", quote = F, row.names=T, col.names=F)
############### LASSO cox model ####################################
## no demographic and pharmaceutical data
library(survival)
coxdata <- merge(survival[,c(1,5:6)], df, by.x = 1, by.y = 1)
coxdata <- na.omit(coxdata)
coxdata <- coxdata[coxdata[,2]>0,]
x <- as.matrix(coxdata[, -1:-3])
time <- coxdata[,2]
event <- coxdata[,3]
event <- ifelse(event=="Dead", 1, 0)
y <- Surv(time, event)
fit=glmnet(x,y,family="cox")
cvfit=cv.glmnet(x,y,family="cox", nfolds = 5)
small.lambda.index <- which(cvfit$lambda == cvfit$lambda.min)
small.lambda.betas <- cvfit$glmnet.fit$beta[,small.lambda.index]
betas <- small.lambda.betas[small.lambda.betas != 0]
length(betas)
betas <- data.frame(betas)
write.table(betas, "LASSO_cox_methylation_GBM.txt", sep = "\t", quote = F, row.names=T, col.names=F)




