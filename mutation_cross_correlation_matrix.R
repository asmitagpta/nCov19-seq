library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(ggrepel)
library(stringr)
library(gridExtra)
library(ggpubr)
library(rowr)
library(factoextra)
library(corrplot)
library(ggdendro)

all.data <- read.csv('location_summary_till_date.tsv', sep='\t')
seq.data <- read.csv('sequenced_samples_till_date.txt', header = FALSE)
ind.seq <- which(all.data$INTERNAL.ID %in% seq.data$V1)
all.data.seq <- all.data[ind.seq,]
all.data.seq$DATE.OF.SAMPLE.COLLECTION <- parse_date_time(all.data.seq$DATE.OF.SAMPLE.COLLECTION, 
                                                          orders = c('m/d/y','Y-m-d','m/d/Y','d.m.Y','d-m-Y','d.m.y','m.d.y','m.d.Y'))
all.data.seq$DATE.OF.SAMPLE.COLLECTION <- as.Date(all.data.seq$DATE.OF.SAMPLE.COLLECTION)

all.data2 <- read.csv('location_summary_till_date_2.tsv', sep = '\t')
seq.data2 <- read.csv('samples_sequenced_15Sep_30Oct', header = FALSE)
ind.seq2 <- which(all.data2$INTERNAL.ID %in% seq.data2$V1)
all.data.seq2 <- all.data2[ind.seq2,]
all.data.seq2$DATE.OF.SAMPLE.COLLECTION <- parse_date_time(all.data.seq2$DATE.OF.SAMPLE.COLLECTION, 
                                                           orders = c('m/d/y','Y-m-d','m/d/Y','d.m.Y','d-m-Y','d.m.y'))
all.data.seq2$DATE.OF.SAMPLE.COLLECTION <- as.Date(all.data.seq2$DATE.OF.SAMPLE.COLLECTION)

ivarjun20 <- read.csv('all_mutations_from_ivar_set210_april_jun2020.csv', sep = '\t', stringsAsFactors = FALSE)
ivardec20 <- read.csv('all_mutations_from_iVar_set40_aug_dec2020.csv', sep = '\t', stringsAsFactors = FALSE)
ivarjul21 <- read.csv('all_mutations_from_iVar_from_jan2021.csv', sep = '\t', stringsAsFactors = FALSE)
ivaraug21 <- read.csv('all_mutations_from_iVar_from_aug2021.csv', sep = '\t', stringsAsFactors = F)
ivaraug1_16 <- read.csv('all_mutations_from_iVar_from_aug1_16_2021.csv', sep = '\t', stringsAsFactors = F)
ivar_sep <- read.csv('all_mutations_from_iVar_from_aug_sep.csv', sep = '\t', stringsAsFactors = F)
ivar_sep2 <- read.csv('all_mutations_from_iVar_sep16_30.csv', sep = '\t', stringsAsFactors = F)
ivar_oct <- read.csv('all_mutations_from_iVar_oct1_15.csv', sep = '\t', stringsAsFactors = F)
ivar_oct2 <- read.csv('all_mutations_from_iVar_oct15_30.csv', sep = '\t', stringsAsFactors = F)

#selecting dataframes by common columns
ivarjun20.sub <- ivarjun20 %>% select('POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','sample_id') 
ivardec20.sub <- ivardec20 %>% select('POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','sample_id')
ivarjul21.sub <- ivarjul21 %>% select('POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','sample_id')                
ivaraug21.sub <- ivaraug21 %>% select('POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','sample_id')                
ivaraug1_16.sub <- ivaraug1_16 %>% select('POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','sample_id')                
ivarsep.sub <- ivar_sep %>% select('POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','sample_id')                
ivarsep2.sub <- ivar_sep2 %>% select('POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','sample_id')                
ivaroct.sub <- ivar_oct %>% select('POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','sample_id')                
ivaroct2.sub <- ivar_oct2 %>% select('POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','sample_id')

ivardecjuly21 <- rbind(ivardec20.sub, ivarjul21.sub, ivaraug21.sub, ivaraug1_16.sub, ivarsep.sub, ivarsep2.sub, ivaroct.sub, ivaroct2.sub)
ivardecjuly21[,c('dep','type','codon','aa','cols5','gene','gene_type','coding','gp01','cols10','cols11','cols12')] <- str_split_fixed(ivardecjuly21$INFO,'[|]',12)

#Replacind odd portions from sample IDs
ivardecjuly21$sample_id <- str_replace(ivardecjuly21$sample_id,'../Data_12Jun/igib_upload/','')
ivardecjuly21$sample_id <- str_replace(ivardecjuly21$sample_id,'../Data_20Jul/','')
ivardecjuly21$sample_id <- str_replace(ivardecjuly21$sample_id,'../Data_28Jul/','')
ivardecjuly21$sample_id <- str_replace(ivardecjuly21$sample_id,'../Data_02Aug/','')
ivardecjuly21$sample_id <- str_replace(ivardecjuly21$sample_id,'../Data_03Jul/','')
ivardecjuly21$sample_id <- str_replace(ivardecjuly21$sample_id,'../Data_05Aug/','')
ivardecjuly21$sample_id <- str_replace(ivardecjuly21$sample_id,'../Data_30Sep/','')


ivardecjuly21[,c('samp_id_base','samp_id_rest')] <- str_split_fixed(ivardecjuly21$sample_id,'/',2)
ivardecjuly21[,c('base','rest')] <- str_split_fixed(ivardecjuly21$samp_id_base,'_S[0-9]+',2)

ivarjun20.sub[,c('dep','type','effect','gene','gp01','gene_type','gp02','coding','cols9','nt','aa','cols12')] <- str_split_fixed(ivarjun20$INFO,'[|]',12)

#fix nanopore sample names issues
ont.ind <- grep('artic',ivarjun20.sub$sample_id)
ont.df <- ivarjun20.sub[ont.ind,]
ivarjun20.tmp <- anti_join(ivarjun20.sub, ont.df)
ont.df[,c('base','rest')] <- str_split_fixed(ont.df$sample_id,'_',2)
ivarjun20.tmp[,c('base','rest')] <- str_split_fixed(ivarjun20.tmp$sample_id,'_S[0-9]+',2)
ivarjun20.sub <- rbind(ivarjun20.tmp, ont.df)

ivardecjuly21 <- ivardecjuly21 %>% select('POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','sample_id','dep','aa','gene','base','type')
ivarjun20.sub <- ivarjun20.sub %>% select('POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','sample_id','dep','aa','gene','base','type')
ivar.all <- rbind(ivarjun20.sub,ivardecjuly21)

#merge the two ont and illumina dataframes again
ivar.all$base <- str_replace(ivar.all$base,'XI','Xi')
ivar.all$base <- str_replace(ivar.all$base,'Xi758','XI758') # temporary change

samplenames <- as.data.frame(ivar.all$base)

ivar.all$aa <- str_replace(ivar.all$aa,'p.','')
ivar.all$aa <- str_replace_all(ivar.all$aa,c('Asp'='D', 'Arg'='R','Gly'='G','Pro'='P','Leu'='L','Lys'='K','Val'='V',
                                               'Ser'='S','Cys'='C','Ala'='A','Phe'='F','Glu'='E','Thr'='T','Tyr'='Y',
                                               'Ile'='I','His'='H','Asn'='N','Gln'='Q','Met'='M','Trp'='W'))

#mutation frequencies globally
data.plot <- ivar.all %>% 
  mutate(tot_samp=n_distinct(base)) %>% 
  group_by(POS, REF, ALT) %>% 
  mutate(total_num=n()) %>% 
  mutate(freq=total_num/tot_samp) %>% 
  arrange(POS)


df.split <- data.plot %>% 
  filter(!type %in% c('synonymous_variant','SILENT')) %>% 
  filter(freq >= 0.03) %>% 
  group_by(POS, REF, ALT) %>% 
  group_split()

emp.df <- data.frame()
df <- data.frame(val=character())

cbind.fill <- function(df, df2){
  nm <- list(df, df2) 
  nm<-lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
  rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

for(i in 1:length(df.split)){
  emp.df <- df.split[[i]] #to reduce type error
  colname <- unique(paste0(emp.df$gene,'_',emp.df$aa))[1]
  #colname <- unique(paste0(emp.df$gene,'_',emp.df$aa))[1]
  print(colname)
  df2 <- data.frame(charac = as.character(), stringsAsFactors = FALSE)
  df2 <- emp.df$base
  df2 <- as.data.frame(df2)
  colnames(df2) <- colname
  #df <- cbind.fill(df, df2, fill = NA\)
  df= cbind.fill(df, df2)
}
df <- df[,-1]
#data.plot$aa.num <- data.plot$aa
#data.plot$aa.num <- str_replace(data.plot$aa.num, '[A-Z]+','')
#data.plot$aa.num <- str_replace(data.plot$aa.num, '[A-Z]+','')
#data.plot$aa.num <- str_replace(data.plot$aa.num,'fs','')
#data.plot$aa.num <- str_replace(data.plot$aa.num,'del','')
#data.plot$aa.num <- as.integer(data.plot$aa.num)

# start plotting corrplot
samp.data <- df
cols <- colnames(samp.data)
#cols[grep("*_$",cols)] <- str_replace(cols[grep("*_$",cols)],'_$','_')
c = cols[grep("*_$",cols, invert = T)]
samp.data <- samp.data[,c]

#get pangolin info on samples with delta lineage
pango <- read.csv('pangolin_analysis_all_till_date.csv')
#pango.16aug <- read.csv('pangolin_16aug_2021.csv')
#pango.27aug <- read.csv('pangolin_27aug_2021.csv')
#pango.27sep <- read.csv('pangolin_27sep.csv')
#pango.30sep <- read.csv('pangolin_sep15.csv')
pango_aug_sep <- read.csv('pangolin_aug_sep.csv')
pango.30oct <- read.csv('lineage_results_till_30Oct.csv')
pango <- rbind(pango,pango_aug_sep, pango.30oct)

pango$sampID <- pango$Sequence.name
ont.ind <- grep('ONT_',pango$sampID)
ont.df <- pango[ont.ind,]
pango.sub <- anti_join(pango, ont.df)

pango.sub$sampID <- str_replace(pango.sub$sampID,'Consensus_','')
ont.df$sampID <- str_replace(ont.df$sampID,'ONT_','')
pango.sub[,c('base','rest')] <- str_split_fixed(pango.sub$sampID,'_S[0-9]+',2)
pango.sub$base <- str_replace(pango.sub$base,'XI','Xi')
pango.sub$base <- str_replace(pango.sub$base,'48001','ncov-48001')
ont.df[,c('base','rest')] <- str_split_fixed(ont.df$sampID,'_',2)
pango.all <- rbind(pango.sub, ont.df)

delta.samp.id <- pango.all %>% 
  filter(Lineage %in% c('B.1.617.2','AY.12','AY.4','B.1.617.1','AY.20','AY.26')) %>% 
  select(base, Lineage)

#select samples sets here 1. All Vaccinated /complete etc.
all.data.tot <- rbind(all.data.seq, all.data.seq2)
all.data.seq.copy <- all.data.tot %>% 
  #filter(VACCINATED != 'YES') %>% 
  select(INTERNAL.ID, DATE.OF.SAMPLE.COLLECTION)

#select delta samples from marh onwards from all data files, 87, 150, 200 samples size dataset
all.data.seq.copy <- all.data.seq.copy %>% filter(INTERNAL.ID %in% delta.samp.id$base) %>% filter(DATE.OF.SAMPLE.COLLECTION > as.Date('2021-02-28'))
all.data.seq.copy <- all.data.seq.copy[sample(nrow(all.data.seq.copy), 350),]

# make correlation matrix
options(error=recover)

for(i in colnames(samp.data)){
 # tmp <- str_split_fixed(samp_data.c.sub[[i]], '_', 2)
  print(i)
  tmp <- samp.data[,i]
  tmp <- as.data.frame(tmp)
  colnames(tmp) <- c('base')
  
  ind <- unlist(sapply(tmp$base, function(y){
    if(y %in% all.data.seq.copy$INTERNAL.ID){
      which(as.character(all.data.seq.copy$INTERNAL.ID) == as.character(y))
    }
  }))
  
  all.data.seq.copy[[i]] <- with(all.data.seq.copy, 0)  
  all.data.seq.copy[[i]][ind] <- 1 
}

all.data.seq.copy <- all.data.seq.copy %>% arrange(DATE.OF.SAMPLE.COLLECTION)
c2 <- all.data.seq.copy[,1:2]
t2 <- all.data.seq.copy[,c(-1,-2)]
t2.filter <- t2[, !sapply(t2, function(x) {sd(x) == 0} )]
all.data.seq.copy <- cbind(c2, t2.filter)

# cor test matrix
cor.mtest <- function(mat,...) {
  mat <- all.data.seq.copy[,c(-1,-2)]
  #mat <- mat[, !sapply(mat, function(x) {sd(x) == 0} )]
  #mat <- Filter(function(x) sd(x) != 0, mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat 
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(all.data.seq.copy[,c(-1,-2)])
head(p.mat[, 1:5])
#corrplot(cor(meta_copy[,3:19]), method = 'circle', order = 'hclust', type='lower', addrect = 1,
#         p.mat = p.mat, sig.level = 0.01)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
setEPS()
postscript(file='paper-II/figures/mutation_correlation_matrix_non_silent_all_samples_till_date.eps', width = 17, height = 16)
corrplot(cor(all.data.seq.copy[,c(-1,-2)]), method="square", mar = c(0.5,7,1.5,1),
         col=brewer.pal(8,'RdYlBu'),  
         type="lower", order="hclust", 
         hclust.method = 'ward.D2',
         #addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=35, tl.cex = 1.1, cl.cex = 2, cl.offset = 0.1, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = 'blank',
         # hide correlation coefficient on the principal diagonal
         diag=FALSE, addgrid.col = NA, shade.col = 'black' )
dev.off()




