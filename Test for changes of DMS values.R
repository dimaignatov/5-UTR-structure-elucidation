D <- read.table('Experiment.txt', sep='\t', header=TRUE)

# Exclude the positions with coverage less than 0 in any of the samples and >1% mismatches in control
D <- D[D$cvg_total>0 & D$cvg_invitro>0 &
      D$cvg_37A>0 &  D$cvg_37B>0 & D$cvg_26>0 & D$cvg_26shift>0
       & D$cvg_37shift>0 & D$cvg_hfq>0 & D$cvg_lhrA>0 & D$cvg_prfA>0,]
D <- D[D$cvg_K>0 & D$mis_K/D$cvg_K<0.01,]

# Unite samples 37A and 37B
D <- within(D, {
  cvg_37 <- cvg_37A + cvg_37B
  mis_37 <- mis_37A + mis_37B
})

# For each sample, add DMS reactivity values separately for A and C
react_A_total = 0.058; react_C_total = 0.047
react_A_invitro = 0.041; react_C_invitro = 0.026
react_A_37A = 0.048; react_C_37A = 0.042
react_A_37B = 0.050; react_C_37B = 0.045
react_A_37 = 0.049; react_C_37 = 0.043
react_A_26 = 0.035; react_C_26 = 0.030
react_A_26shift = 0.037; react_C_26shift = 0.029
react_A_37shift = 0.051; react_C_37shift = 0.041
react_A_hfq = 0.050; react_C_hfq = 0.043
react_A_lhrA = 0.047; react_C_lhrA = 0.041
react_A_prfA = 0.047; react_C_prfA = 0.041

D <- within(D, {
  react_total <- ifelse(D$Nucleotide=='A', react_A_total, react_C_total)
  react_invitro <- ifelse(D$Nucleotide=='A', react_A_invitro, react_C_invitro)
  react_37A <- ifelse(D$Nucleotide=='A', react_A_37A, react_C_37A)
  react_37B <- ifelse(D$Nucleotide=='A', react_A_37B, react_C_37B)
  react_37 <- ifelse(D$Nucleotide=='A', react_A_37, react_C_37)
  react_26 <- ifelse(D$Nucleotide=='A', react_A_26, react_C_26)
  react_26shift <- ifelse(D$Nucleotide=='A', react_A_26shift, react_C_26shift)
  react_37shift <- ifelse(D$Nucleotide=='A', react_A_37shift, react_C_37shift)
  react_hfq <- ifelse(D$Nucleotide=='A', react_A_hfq, react_C_hfq)
  react_lhrA <- ifelse(D$Nucleotide=='A', react_A_lhrA, react_C_lhrA)
  react_prfA <- ifelse(D$Nucleotide=='A', react_A_prfA, react_C_prfA)
})

# Calculate DMS values and normalize to reactivity
D <- within(D, {
 dms_total <- mis_total/(cvg_total * react_total)
 dms_invitro <- mis_invitro/(cvg_invitro * react_invitro)
 dms_37A <- mis_37A/(cvg_37A * react_37A)
 dms_37B <- mis_37B/(cvg_37B * react_37B)
 dms_37 <- mis_37/(cvg_37 * react_37)
 dms_26 <- mis_26/(cvg_26 * react_26)
 dms_26shift <- mis_26shift/(cvg_26shift * react_26shift)
 dms_37shift <- mis_37shift/(cvg_37shift * react_37shift)
 dms_hfq <- mis_hfq/(cvg_hfq * react_hfq)
 dms_lhrA <- mis_lhrA/(cvg_lhrA * react_lhrA)
 dms_prfA <- mis_prfA/(cvg_prfA * react_prfA)
})

# Remove positions with DMS values > 10
D <- D[D$dms_26 < 10 | D$dms_37 < 10 | D$dms_invitro < 10,]

# Modify and write the data to the file
D$Feature <- gsub("\\(utr\\d+\\)", "", D$Feature)
D$Feature <- gsub("utr", "lmo", D$Feature)
#write.table(D, file='Cvg_mis_dms.txt', quote=F, sep = "\t", row.names = FALSE)

# Perform Fisher exact test for differential DMS values for each sample relative to dms_37A
ft <- function(react_A, react_B, cvg_A, cvg_B, mis_A, mis_B){
  cvg_A <- round(cvg_A * react_A/react_B)
  M <- matrix(c(cvg_A, cvg_B, mis_A, mis_B), nrow = 2)
  res <- fisher.test(M)
  return(res)
}

pval_37B <- c(); odds_37B <- c()
pval_invitro <- c(); odds_invitro <- c()
pval_26 <- c(); odds_26 <- c()
pval_26shift37 <- c(); odds_26shift37 <- c()
pval_37shift26 <- c(); odds_37shift26 <- c()
pval_hfq <- c(); odds_hfq <- c()
pval_lhrA <- c(); odds_lhrA <- c()
pval_prfA <- c(); odds_prfA <- c()

for(i in 1:nrow(D)) {
  r <- D[i,]
  
  res_37B <- ft(r$react_37A, r$react_37B, r$cvg_37A, r$cvg_37B, r$mis_37A, r$mis_37B)
  pval_37B <- c(pval_37B, res_37B$p.value)
  odds_37B <- c(odds_37B, as.numeric(res_37B$estimate))
  
  res_invitro <- ft(r$react_37, r$react_invitro, r$cvg_37, r$cvg_invitro, r$mis_37, r$mis_invitro)
  pval_invitro <- c(pval_invitro, res_invitro$p.value)
  odds_invitro <- c(odds_invitro, as.numeric(res_invitro$estimate))
  
  res_26 <- ft(r$react_26, r$react_37, r$cvg_26, r$cvg_37, r$mis_26, r$mis_37)
  pval_26 <- c(pval_26, res_26$p.value)
  odds_26 <- c(odds_26, as.numeric(res_26$estimate))
  
  res_26shift37 <- ft(r$react_26, r$react_37shift, r$cvg_26, r$cvg_37shift, r$mis_26, r$mis_37shift)
  pval_26shift37 <- c(pval_26shift37, res_26shift37$p.value)
  odds_26shift37 <- c(odds_26shift37, as.numeric(res_26shift37$estimate))
  
  res_37shift26 <- ft(r$react_26shift, r$react_37, r$cvg_26shift, r$cvg_37, r$mis_26shift, r$mis_37)
  pval_37shift26 <- c(pval_37shift26, res_37shift26$p.value)
  odds_37shift26 <- c(odds_37shift26, as.numeric(res_37shift26$estimate))
  
  res_hfq <- ft(r$react_37, r$react_hfq, r$cvg_37, r$cvg_hfq, r$mis_37, r$mis_hfq)
  pval_hfq <- c(pval_hfq, res_hfq$p.value)
  odds_hfq <- c(odds_hfq, as.numeric(res_hfq$estimate))
  
  res_lhrA <- ft(r$react_37, r$react_lhrA, r$cvg_37, r$cvg_lhrA, r$mis_37, r$mis_lhrA)
  pval_lhrA <- c(pval_lhrA, res_lhrA$p.value)
  odds_lhrA <- c(odds_lhrA, as.numeric(res_lhrA$estimate))
  
  res_prfA <- ft(r$react_37, r$react_prfA, r$cvg_37, r$cvg_prfA, r$mis_37, r$mis_prfA)
  pval_prfA <- c(pval_prfA, res_prfA$p.value)
  odds_prfA <- c(odds_prfA, as.numeric(res_prfA$estimate))
}

# Define nucleotides with significant changes of DMS value for each comparison
fdr_37B <- p.adjust(pval_37B, method = 'fdr')
diff_37B <- ifelse(fdr_37B < 0.05 & 
                     (odds_37B < 0.8 | odds_37B > 1.20) &
                     abs(D$dms_37B-D$dms_37A) > 0.1, TRUE, FALSE)
D.37B <- cbind(D[,1:7], 
               cvg_37A=D$cvg_37A, mis_37A=D$mis_37A, cvg_37B=D$cvg_37B, mis_37B=D$mis_37B,
               dms_37A=D$dms_37A, dms_37B=D$dms_37B, 
               pval_37B, fdr_37B, odds_37B, diff_37B)
write.table(D.37B, file='Diff_37B.txt', quote=F, sep = "\t", row.names = FALSE)

fdr_26 <- p.adjust(pval_26, method = 'fdr')
diff_26 <- ifelse(fdr_26 < 0.05 & 
                    (odds_26 < 0.8 | odds_26 > 1.2) &
                    abs(D$dms_26-D$dms_37) > 0.1, TRUE, FALSE)
D.26 <- cbind(D[,1:7], 
                   cvg_37=D$cvg_37, mis_37=D$mis_37, cvg_26=D$cvg_26, mis_26=D$mis_26,
                   dms_37=D$dms_37, dms_26=D$dms_26, 
                   pval_26, fdr_26, odds_26, diff_26)
write.table(D.26, file='Diff_26.txt', quote=F, sep = "\t", row.names = FALSE)

fdr_26shift37 <- p.adjust(pval_26shift37, method = 'fdr')
diff_26shift37 <- ifelse(fdr_26shift37 < 0.05 & 
                    (odds_26shift37 < 0.8 | odds_26shift37 > 1.2) &
                    abs(D$dms_37shift-D$dms_26) > 0.1, TRUE, FALSE)
D.26shift37 <- cbind(D[,1:7], 
              cvg_37shift=D$cvg_37shift, mis_37shift=D$mis_37shift, cvg_26=D$cvg_26, mis_26=D$mis_26,
              dms_37shift=D$dms_37shift, dms_26=D$dms_26, 
              pval_26shift37, fdr_26shift37, odds_26shift37, diff_26shift37)
write.table(D.26shift37, file='Diff_26shift37.txt', quote=F, sep = "\t", row.names = FALSE)

fdr_hfq <- p.adjust(pval_hfq, method = 'fdr')
diff_hfq <- ifelse(fdr_hfq < 0.05 & 
                     (odds_hfq < 0.8 | odds_hfq > 1.2) &
                     abs(D$dms_hfq-D$dms_37) > 0.1, TRUE, FALSE)
D.hfq <- cbind(D[,1:7], 
                   cvg_hfq=D$cvg_hfq, mis_hfq=D$mis_hfq, cvg_37=D$cvg_37, mis_37=D$mis_37,
                   dms_37=D$dms_37, dms_hfq=D$dms_hfq, 
                   pval_hfq, fdr_hfq, odds_hfq, diff_hfq)
write.table(D.hfq, file='Diff_hfq.txt', quote=F, sep = "\t", row.names = FALSE)

fdr_lhrA <- p.adjust(pval_lhrA, method = 'fdr')
diff_lhrA <- ifelse(fdr_lhrA < 0.05 & 
                      (odds_lhrA < 0.8 | odds_lhrA > 1.2) &
                      abs(D$dms_lhrA-D$dms_37) > 0.1, TRUE, FALSE)
D.lhrA <- cbind(D[,1:7], 
               cvg_lhrA=D$cvg_lhrA, mis_lhrA=D$mis_lhrA, cvg_37=D$cvg_37, mis_37=D$mis_37,
               dms_37=D$dms_37, dms_lhrA=D$dms_lhrA, 
               pval_lhrA, fdr_lhrA, odds_lhrA, diff_lhrA)
write.table(D.lhrA, file='Diff_lhrA.txt', quote=F, sep = "\t", row.names = FALSE)

fdr_prfA <- p.adjust(pval_prfA, method = 'fdr')
diff_prfA <- ifelse(fdr_prfA < 0.05 & 
                      (odds_prfA < 0.8 | odds_prfA > 1.2) &
                      abs(D$dms_prfA-D$dms_37) > 0.1, TRUE, FALSE)
D.prfA <- cbind(D[,1:7], 
                cvg_prfA=D$cvg_prfA, mis_prfA=D$mis_prfA, cvg_37=D$cvg_37, mis_37=D$mis_37,
                dms_37=D$dms_37, dms_prfA=D$dms_prfA, 
                pval_prfA, fdr_prfA, odds_prfA, diff_prfA)
write.table(D.prfA, file='Diff_prfA.txt', quote=F, sep = "\t", row.names = FALSE)
