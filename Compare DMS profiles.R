library(Biostrings)
library(ggplot2)
library(reshape)
setwd('E:/NGS/___Materials for publication___/Statistical analysis')
ft <- readDNAStringSet("Features.fa")
D <- read.table('Experiment.txt', sep='\t', header=TRUE)

# Exclude the positions with coverage less than 250 in any of the samples and >1% mismatches in control
D <- D[D$cvg_total>0 & D$cvg_invitro>0 &
         D$cvg_37A>0 &  D$cvg_37B>0 & D$cvg_26>0 & D$cvg_26shift>0
       & D$cvg_37shift>0 & D$cvg_hfq>0 & D$cvg_lhrA>0 & D$cvg_prfA>0,]

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

# If necessary, Adjust coordinates for SRP RNA
#D <- D[grepl('SRP', D$Feature),]
#D$Position_TSS <- 2784533 - D$Position_genome

# Calculate differences of DMS values for two samples
feature_name = 'utr1364'
D.sel <- D[grepl(feature_name, D$Feature),]
D.sel <- within(D.sel, {
  minus <- dms_37shift - dms_26
  div <- dms_37shift / dms_26
})

# Generate COLORMAP profiles for visualization with VARNA
string <- ft[[feature_name]]
seq_string <- toString(string)
seq_list <- strsplit(seq_string, '')[[1]]

minus_list <- c()
div_list <- c()
nt_list <- c()
pos_list <- c()
lab_list <- c()

for(i in 1:length(string)){
  nt <- seq_list[i]
  if(is.element(i, D.sel$Position_TSS)){
    r <- D.sel[D.sel$Position_TSS==i,]
    minus <- r$minus                 # Change the sample name here!
    div <- r$div
  }else{
    minus <- 50
    div <- 50
  }
  minus_list <- c(minus_list, minus)
  div_list <- c(div_list, div)
  nt_list <- c(nt_list, nt)
  pos_list <- c(pos_list, i)
  lab_list <- c(lab_list, paste(nt, i, sep='\n'))
}
df_colormap <- cbind.data.frame(pos_list, nt_list, minus_list)
write.table(df_colormap, paste(feature_name, 'minus.colormap', sep='_'), 
            sep='\t', quote = FALSE, row.names = FALSE, col.names=FALSE)

# Plot the differences to PDF 
string <- ft[[feature_name]]
seq_string <- toString(string)
seq_list <- strsplit(seq_string, '')[[1]]

minus_list <- c()
div_list <- c()
nt_list <- c()
pos_list <- c()
lab_list <- c()

for(i in 1:length(string)){
  nt <- seq_list[i]
  if(is.element(i, D.sel$Position_TSS)){
    r <- D.sel[D.sel$Position_TSS==i,]
    minus <- r$minus                 
    div <- log(r$div, base=2)
  }else{
    minus <- 0
    div <- 0
  }
  minus_list <- c(minus_list, minus)
  div_list <- c(div_list, div)
  nt_list <- c(nt_list, nt)
  pos_list <- c(pos_list, i)
  lab_list <- c(lab_list, paste(nt, i, sep='\n'))
}
df_plot <- cbind.data.frame(pos_list, nt_list, lab_list, minus_list, div_list)
df_plot$change_list <- df_plot$minus * abs(df_plot$div_list)

file_name <- paste(feature_name, 'minus_profile.pdf', sep='_')
pdf(file_name, width=length(string)/5, height=10)
p <- ggplot(df_plot, aes(x=factor(pos_list), y=minus_list))+
  geom_col(fill='blue', col='black', width=1) +
  scale_x_discrete(labels=lab_list)  +
  coord_cartesian(ylim = c(-5, 5), expand = FALSE) +
  theme_bw()
print(p)
dev.off()

file_name <- paste(feature_name, 'div_profile.pdf', sep='_')
pdf(file_name, width=length(string)/5, height=10)
p <- ggplot(df_plot, aes(x=factor(pos_list), y=div_list))+
  geom_col(fill='blue', col='black', width=1) +
  scale_x_discrete(labels=lab_list)  +
  coord_cartesian(ylim = c(-5, 5), expand = FALSE) +
  theme_bw()
print(p)
dev.off()

file_name <- paste(feature_name, 'change_profile.pdf', sep='_')
pdf(file_name, width=length(string)/5, height=10)
p <- ggplot(df_plot, aes(x=factor(pos_list), y=change_list))+
  geom_col(fill='blue', col='black', width=1) +
  scale_x_discrete(labels=lab_list)  +
  coord_cartesian(ylim = c(-5, 5), expand = FALSE) +
  theme_bw()
print(p)
dev.off()