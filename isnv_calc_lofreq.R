library(tidyverse)
library(ggplot2)
library(lubridate)
library(ggrepel)
library(RColorBrewer)
library(viridis)
library(rowr)

all.data <- read.csv('location_summary_till_date.tsv', sep='\t')
seq.data <- read.csv('sequenced_samples_till_date.txt', header = FALSE)

ind.seq <- which(all.data$INTERNAL.ID %in% seq.data$V1)
all.data.seq <- all.data[ind.seq,]
all.data.seq$DATE.OF.SAMPLE.COLLECTION <- parse_date_time(all.data.seq$DATE.OF.SAMPLE.COLLECTION, 
                                        orders = c('m/d/y','Y-m-d','m/d/Y','d.m.Y','d-m-Y','d.m.y','m.d.y','m.d.Y'))
all.data.seq$DATE.OF.SAMPLE.COLLECTION <- as.Date(all.data.seq$DATE.OF.SAMPLE.COLLECTION)

# read lofreq mutation file from here
lo.all <- read.csv('lofreq_mutations.csv', sep = '\t', stringsAsFactors = FALSE)

lo.all$sample_id <- str_replace(lo.all$sample_id, '.lofreq_snpeff.vcf','')
lo.all$sample_id <- str_replace(lo.all$sample_id, '.lofreq_ann.vcf','')
lo.all$sample_id <- str_replace(lo.all$sample_id,'../Data_30Sep/','')

lo.all[,c('base','rest')] <- str_split_fixed(lo.all$sample_id,'/', 2)
lo.all$base <- str_replace(lo.all$base, '_L001', '')
lo.all$base <- str_replace(lo.all$base, '_R1_001', '')
lo.all$base <- str_replace(lo.all$base, '_S[0-9]+', '')

lo.all$aa <- str_replace(lo.all$aa,'p.','')
lo.all$aa <- str_replace_all(lo.all$aa,c('Asp'='D', 'Arg'='R','Gly'='G','Pro'='P','Leu'='L','Lys'='K','Val'='V',
                                               'Ser'='S','Cys'='C','Ala'='A','Phe'='F','Glu'='E','Thr'='T','Tyr'='Y',
                                               'Ile'='I','His'='H','Asn'='N','Gln'='Q','Met'='M','Trp'='W'))

lo.all$aa <- str_replace(lo.all$aa, 'R203K','RG203KR')
lo.all$aa <- str_replace(lo.all$aa, 'G204R','RG203KR')

#identify isnv by grouping on the basis of sample id followed by position; and filtering positions mutated more than once in a sample
isnv <- lo.all %>% group_by(base, POS) %>% mutate(num_alterations=n()) %>% filter(num_alterations>1)

#now calculate the sample frequency of isnv positions
isnv.samp.freq <- isnv %>% group_by(POS) %>% mutate(num_samples=n_distinct(base)) 
isnv.samp.freq[,paste0('info',1:5)] <- str_split_fixed(isnv.samp.freq$INFO,';',5)
isnv.samp.freq[,paste0('info',6:10)] <- str_split_fixed(isnv.samp.freq$info5,'[|]',5)
isnv.samp.freq <- isnv.samp.freq %>% mutate(lab=paste0(REF,POS), lab.nt=paste0(REF,POS,ALT,'(',aa,')')) # create labels
isnv.samp.freq$info1 <- str_replace(isnv.samp.freq$info1,'DP=','') 
isnv.samp.freq$info2 <- str_replace(isnv.samp.freq$info2,'AF=','') 
isnv.samp.freq$info1 <- as.integer(isnv.samp.freq$info1)
isnv.samp.freq$info2 <- as.numeric(isnv.samp.freq$info2)

#filter allele allele with AF > 1% and sample frequency > 5

isnv.samp.freq.filtered <- isnv.samp.freq %>% filter(info2 > 0.02 & num_samples > 5)

data.labels <- isnv.samp.freq.filtered %>% 
   select(POS, num_samples, lab) %>% 
   mutate(sel_label = ifelse(num_samples>10, lab, '')) %>% 
   distinct(sel_label, .keep_all=TRUE)

#plot iSNV frequency 
p1 <- isnv.samp.freq.filtered %>% 
   ggplot(aes(x=POS, y=num_samples)) +
   geom_point(fill='firebrick3', shape=21, size=6, color='white') +
   #geom_point(shape=21, size=5, fill='firebrick', color='white') +
   geom_segment(aes(x=POS, xend=POS, y=0, yend=num_samples), color='grey76')+ 
   theme(panel.background = element_blank(), axis.line = element_line(),
         axis.text = element_text(size = 12, face = 'bold'), 
         axis.title = element_text(size=14, face = 'bold'), 
         title = element_text(size = 14, face = 'bold'), legend.position = 'bottom', 
         legend.direction = 'horizontal', legend.text = element_text(size=10)) +
   geom_text_repel(data = data.labels, label=data.labels$sel_label, size=5,
                   vjust=0.5, ylim = c(-Inf, Inf), xlim = c(-Inf, Inf),
                   segment.size = 0.3, min.segment.length = 0, max.overlaps = 20,
                   nudge_y = 6) +
  coord_cartesian(xlim=c(-7,29900), expand = T)+
  labs(x='Genomic position', y='Num. of samples') +
annotate('rect',xmin = 1, ymin = 0, xmax = 265, ymax = 380, alpha=0.2, fill='grey76')+
annotate('rect',xmin = 266, ymin = 0, xmax = 13468, ymax = 380, alpha=0.2,fill='orange')+
annotate('rect',xmin = 13488, ymin = 0, xmax = 21555, ymax = 380, alpha=0.2, fill='steelblue')+
annotate('rect',xmin = 21563, ymin = 0, xmax = 25384, ymax = 380, alpha=0.2,fill='chartreuse')+
annotate('rect',xmin = 25393 , ymin = 0, xmax = 26220, ymax = 380, alpha=0.2, fill='red') +
annotate('rect',xmin = 26245 , ymin = 0, xmax = 26472, ymax = 380, alpha=0.2, fill='olivedrab2')+
annotate('rect',xmin = 26523 , ymin = 0, xmax = 27191, ymax = 380, alpha=0.2, fill='brown')+
annotate('rect',xmin = 27202 , ymin = 0, xmax = 27387, ymax = 380, alpha=0.2, fill='skyblue')+
annotate('rect',xmin = 27394 , ymin = 0, xmax = 27759, ymax = 380, alpha=0.2, fill='yellow2')+
annotate('rect',xmin  = 27756 , ymin = 0, xmax = 27887, ymax = 380, alpha=0.2, fill='goldenrod2')+
annotate('rect',xmin = 27894 , ymin = 0, xmax = 28259, ymax = 380, alpha=0.2, fill='magenta4')+
annotate('rect',xmin = 28274 , ymin = 0, xmax = 29533, ymax = 380, alpha=0.2, fill='maroon2')+
annotate('rect',xmin = 29558 , ymin = 0, xmax = 29674, ymax = 380, alpha=0.2, fill='gray8')+
 
#add nsp lines
annotate('segment',x = 266, y = 380, xend = 266, yend = 390)+
annotate('segment',x = 805, y = 380, xend = 805, yend = 390)+
annotate('segment',x = 2719, y = 380, xend = 2719, yend = 390)+
annotate('segment',x = 8554, y = 380, xend = 8554, yend = 390)+
annotate('segment',x = 10054, y = 380, xend = 10054, yend = 390) +
annotate('segment',x = 10972, y = 380, xend = 10972, yend = 390)+
   annotate('segment',x = 11842, y = 380, xend = 11842, yend = 390)+
   annotate('segment',x = 12091, y = 380, xend = 12091, yend = 390)+
   annotate('segment',x = 12685, y = 380, xend = 12685, yend = 390)+
   annotate('segment',x  = 13024, y = 380, xend = 13024, yend = 390)+
   annotate('segment',x = 13441, y = 380, xend = 13441, yend = 390)+
   annotate('segment',x = 16236, y = 380, xend = 16236, yend = 390)+
   annotate('segment',x = 18039 , y = 380, xend = 18039, yend = 390)+
   annotate('segment',x = 19620, y = 380, xend = 19620, yend = 390)+
   annotate('segment',x = 20658 , y = 380, xend = 20658, yend = 390)+
   annotate('segment',x = 21552, y = 380, xend = 21552, yend = 390)

p1 

ggsave(p1, filename = "paper-II/figures/isnv_lofreq.png", width = 14, height = 8, dpi = 300)
ggsave(p1, filename = "paper-II/figures/isnv_lofreq.eps", device = cairo_ps, width = 14, height = 8, dpi = 300)

#analysis on iSNV sites
isnv.samp.freq.filtered <- mutate(isnv.samp.freq.filtered, ntd.change = paste0(REF,'>',ALT))
isnv.samp.freq.filtered <- isnv.samp.freq.filtered %>% mutate(lab.nt2=paste0(REF,POS,ALT)) # create alternative labels


#plot depth vs allele frequency
isnv.samp.freq.filtered %>% 
  ggplot(aes(x=info1, y=info2)) +
  geom_point() +
  geom_vline(xintercept = 100)

#plot depth AF boxplot and type of iSNVs
isnv.mut <- isnv.samp.freq.filtered %>% 
   arrange(POS) %>% 
   ggplot(aes(x=factor(lab.nt2, levels = as.ordered(unique(lab.nt2))), y=info1, fill=effect)) +
   geom_boxplot() +
   #geom_jitter(size=0.1)+
   geom_hline(yintercept = 100, color='red', linetype='dashed') +
   theme(panel.background = element_blank(), axis.line = element_line(),
         axis.text.x = element_text(size = 12, angle = 90, hjust = 1), 
         axis.text.y = element_text(size=12),
         axis.title = element_text(size=14, face = 'bold'), 
         title = element_text(size = 14, face = 'bold'), legend.position = 'bottom', 
         legend.direction = 'horizontal', legend.text = element_text(size=10)) +
   scale_fill_manual(values  = brewer.pal(5,'Set1'), name='')+
   labs(x='',y='Depth in each sample')
isnv.mut  
ggsave(isnv.mut, filename = "paper-II/figures/isnv_mutation_type.png", width = 14, height = 7, dpi = 300)
ggsave(isnv.mut, filename = "paper-II/figures/isnv_mutation_type.eps", device = cairo_ps, width = 14, height = 7, dpi = 300)

#plot isnv sites per sample   
isnv_sample <- isnv.samp.freq.filtered %>% group_by(base) %>% summarise(n=n_distinct(POS))
isnv_sample.plot <- isnv_sample %>% 
   ggplot(aes(x=n)) +
   geom_histogram(binwidth = 1, fill='black', color='black', lwd=0.2) +
   scale_x_continuous(breaks = seq(0,17,3), limits = c(0,17))+
   theme(panel.background = element_blank(),
         panel.grid.major  = element_line(linetype = 'dashed', color = 'gray87'),
         axis.text= element_text(size = 12), 
         axis.title = element_text(size=14))+
   labs(x='Number of iSNV sites',y='Number of samples')
isnv_sample.plot
ggsave(isnv_sample.plot, filename = "paper-II/figures/isnv_per_sample.png", width = 5, height = 5, dpi = 300)
ggsave(isnv_sample.plot, filename = "paper-II/figures/isnv_per_sample.eps", device = cairo_ps, width = 5, height = 5, dpi = 300)

#get sample features and ct values
isnv_sample_id <- isnv.samp.freq.filtered %>% select(base) 
samp.features <- all.data.seq[all.data.seq$INTERNAL.ID %in% isnv_sample_id$base, ]
samp.features <- read.csv('paper-II/figures/isnv_sample_features.csv', stringsAsFactors = F)
samp.features$DATE.OF.SAMPLE.COLLECTION <- parse_date_time(samp.features$DATE.OF.SAMPLE.COLLECTION, 
                                                          orders = c('m/d/y','Y-m-d','m/d/Y','d.m.Y','d-m-Y','d.m.y','m.d.y','m.d.Y'))
samp.features$DATE.OF.SAMPLE.COLLECTION <- as.Date(samp.features$DATE.OF.SAMPLE.COLLECTION)

#get ct value for each base which shows iSNV in isnv_sample dataframe
ind.ct <- unlist(sapply(isnv_sample$base, function(i){
   which(samp.features$INTERNAL.ID == i)}))

isnv_sample <- isnv_sample  %>% 
   mutate(ct=samp.features$Rdrp.Gene..Updated.[ind.ct])
ct_vx_num_isnv <- isnv_sample %>% ggplot(aes(x=n, y=ct)) +
   geom_point(shape=21, size = 2, alpha=0.7,) +
   stat_smooth(method  = 'lm')+
   theme(panel.background = element_blank(), panel.grid = element_line(linetype = 'dashed', color='gray87'),
         axis.text = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'))+
   labs(x='Number of iSNV sites',y='Ct (RdRp gene) value per sample')
ct_vx_num_isnv
ggsave(ct_vx_num_isnv, filename = "paper-II/figures/ct_vs_num_sites.png", width = 6, height = 6, dpi = 300)
ggsave(ct_vx_num_isnv, filename = "paper-II/figures/ct_vs_num_sites.eps", device = cairo_ps, width = 6, height = 6, dpi = 300)

fit <- glm(n ~ ct, data = isnv_sample, family = poisson())
summary(fit)

# KHAMMAM samples
khammam <- samp.features %>% filter(RECEIVED.FROM == 'KHAMMAM') %>% select(INTERNAL.ID)
khammam.isnv <- isnv.samp.freq.filtered %>% filter(base %in% khammam$INTERNAL.ID)

isnv.samp.freq.filtered <- isnv.samp.freq.filtered %>% 
   mutate(sv=ifelse(info2 > 0.90,'snp','isnv'))

aa_isnv <- isnv.samp.freq.filtered %>% filter(sv == 'isnv') %>% select(lab.nt2)

df.split <- isnv.samp.freq.filtered %>% 
            #filter(sv != 'snp') %>% 
            group_by(lab.nt2) %>% 
            group_split()
emp.df <- data.frame()
df <- data.frame(val=character())

cbind.fill <- function(df, df2, fill=NA){
  nm <- list(df, df2) 
  nm<-lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

for(i in 1:length(df.split)){
   emp.df <- df.split[[i]] #to reduce type error
   colname <- unique(emp.df$lab.nt2)
   #colname <- unique(paste0(emp.df$gene,'_',emp.df$aa))[1]
   print(colname)
   df2 <- data.frame(charac = as.character(), stringsAsFactors = FALSE)
   df2 <- emp.df$base
   df2 <- as.data.frame(df2)
   colnames(df2) <- colname
   df <- cbind.fill(df, df2)
}
df <- df[,-1]
samp.data <- df
cols <- colnames(samp.data)

#khammam.copy <- filter(INTERNAL.ID %in% delta.samp.id$base) 

for(i in colnames(samp.data)){
   # tmp <- str_split_fixed(samp_data.c.sub[[i]], '_', 2)
   tmp <- samp.data[,i]
   tmp <- as.data.frame(tmp)
   colnames(tmp) <- c('base')
   
   ind <- unlist(sapply(tmp$base, function(y){
      if(y %in% khammam$INTERNAL.ID){
         which(as.character(khammam$INTERNAL.ID) == as.character(y))
      }
   }))
   
   khammam[[i]] <- with(khammam, 0)  
   khammam[[i]][ind] <- 1 
}

khammam.sub <- khammam[,colSums(khammam != 0) > 0]
khamma.plot <- khammam.sub %>%    
   tidyr::gather(key = 'key', value='vals',-INTERNAL.ID) %>% 
   ggplot(aes(x=INTERNAL.ID,y=key,fill=as.character(vals))) +
   geom_tile(colour="gray89",size=0.9,show.legend = F)+
   scale_y_discrete(expand = c(0,0)) +
   scale_x_discrete(expand = c(0,0)) +
   theme(plot.background = element_blank(), panel.border = element_blank(),
         axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
         axis.text.y = element_text(size = 12), 
         axis.title = element_text(size = 14, face = 'bold'))+
   scale_fill_manual(values = c('white','black')) +
   labs(x='',y='') 
khamma.plot   

######pLOT OTHER ESTIMATES+
ggsave(khamma.plot, filename = "paper-II/figures/khammam.png", width = 7, height = 4, dpi = 300)
ggsave(khamma.plot, filename = "paper-II/figures/khammam.eps", device = cairo_ps, width = 7, height = 4, dpi = 300)

data.plot <- lo.all %>% 
   group_by(POS, REF, ALT) %>% 
   mutate(mut=n()) %>% 
   mutate(tot=n_distinct(base)) %>% 
   mutate(freq=mut/tot, ntd.change=paste0(REF,'>',ALT)) %>% 
   arrange(POS)
snp.ntd.chage <- data.plot %>% filter(freq > 0.05) %>% group_by(ntd.change) %>% summarise(n=n()) %>% arrange(desc(n)) %>% mutate(n=n/100)
isnv.ntd.change <- isnv.samp.freq.filtered %>% group_by(ntd.change) %>% summarise(n=n()) %>% arrange(desc(n)) %>% mutate(n=n*(-1))
   
ntd.changes <- ggplot() +
geom_bar(data=snp.ntd.chage, aes(x=factor(ntd.change, levels = as.ordered(ntd.change)), y=n), stat = 'identity', fill=brewer.pal(3,'Set1')[1]) +
geom_bar(data=isnv.ntd.change, aes(x=factor(ntd.change, levels = as.ordered(ntd.change)), y=n), stat = 'identity', fill=brewer.pal(3,'Set1')[2]) +
scale_y_continuous(breaks = c(-400,-300,-200,-100,0,100,200,300,400,500,600,700,800),
                   labels = c('400','300','200','100','0','10000','20000','30000','40000','50000','60000','70000','80000'))+
   coord_flip() +
   theme(panel.background = element_blank(), axis.line.y = element_line(),
         axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust = 1),
         axis.text.y = element_text(size=12)) +
   labs(y='Nucleotide changes',x='')
 ntd.changes  
ggsave(ntd.changes, filename = "paper-II/figures/isnv_ntd_changes.png", width = 8, height = 8, dpi = 300)
ggsave(ntd.changes, filename = "paper-II/figures/isnv_ntd_changes.eps", device = cairo_ps, width = 8, height = 8, dpi = 300)

# sample feature analysis showing iSNV
pango.all <- read.csv('pangolin_all_data.csv')

isnv.samp <- isnv.samp.freq.filtered %>% select(base) 
isnv.samp.lin <- pango.all %>% filter(base %in% isnv.samp.freq$base)
isnv.samp.lin <- isnv.samp.lin %>% group_by(Lineage) %>% summarise(n=n()) 
isnv.samp.lin <- isnv.samp.lin %>% 
   mutate(Variant = Lineage) %>% 
   mutate(Variant = ifelse(Variant == 'B.1.1.7','B.1.1.7',
                  ifelse(Variant == 'B.1.617.1','B.1.617.1',
                  ifelse(Variant == 'B.1.617.2','B.1.617.2',
                  ifelse(Variant == 'B.1.617.3','B.1.617.3',
                  ifelse(Variant == 'B.1.351','B.1.351',
                  ifelse(Variant == 'P.1','P.1',
                  ifelse(Variant == 'AY.1','AY.1',
                  ifelse(Variant == 'AY.3','AY.3',
                  ifelse(Variant == 'AY.12','AY.12',
                  ifelse(Variant == 'AY.4','AY.4',
                  ifelse(Variant == 'AY.5','AY.5','others')))))))))))) %>% 
   arrange(desc(n))

tot_lins <- pango.all %>% group_by(Lineage) %>% summarise(n=n())
ind <- unlist(sapply(isnv.samp.lin$Lineage, function(i){
   which(tot_lins$Lineage == i)}))

isnv.samp.lin <- isnv.samp.lin  %>% 
   mutate(tot=tot_lins$n[ind])

isnv.samp.lin <- isnv.samp.lin %>% filter(tot > 10 & !n==tot) %>% mutate(prop=n/tot) %>% arrange(desc(prop))

p2 <- isnv.samp.lin %>% 
   ggplot(aes(x=factor(Lineage, levels = as.ordered(Lineage)),y=prop, fill=factor(prop)))+
   geom_bar(stat = 'identity', width=1, color='black', size=0.2) +
   theme(panel.background = element_blank(), 
         axis.line = element_line(), 
         axis.text.x = element_text(size = 12,angle = 90, hjust = 1, vjust = 0.5), 
         axis.text.y = element_text(size = 12),
         axis.title = element_text(size=14, face = 'bold')) +
   guides(fill= 'none') +
  coord_fixed(ratio = 4) +
   scale_fill_manual(values = rev(colorspace::heat_hcl(15)), name='') +
   labs(x='',y='Proportion of samples in each lineage')
p2
ggsave(p2, filename = 'paper-II/figures/isnv_sample_lineages.png',width = 6, height = 6, units = 'in', dpi = 300)   
ggsave(p2, filename = 'paper-II/figures/isnv_sample_lineages.eps', device = cairo_ps, width = 6, height = 6, units = 'in', dpi = 300)   

isnv.samp.features <- all.data.seq %>% filter(INTERNAL.ID %in% isnv.samp$base)


isnv.s <- isnv.samp.freq.filtered %>% 
  filter(gene == 'N') 
  #filter(POS > 2719 & POS < 8555) 
   #%>% 
   #filter(lab %in% c('A21650','A23403','C21618','C23604','C24863','G22081','G22870','G24914','T24703'))

s1 <- isnv.s %>% 
    ggplot(aes(x=lab.nt, y=as.double(info2))) +
    geom_boxplot(show.legend = F, fill='white') +
    facet_wrap(. ~ lab, scales = 'free', nrow = 4, ncol = 4)

s1 <- s1 + geom_jitter(shape=16, position = position_jitter(0.2)) +
   theme(panel.background = element_blank(),
         axis.line = element_line(), axis.text.y = element_text(size = 12),
         axis.title.y = element_text(size = 14),
         axis.text.x = element_text(angle = 90, size = 14),
         legend.position = 'none',
         strip.text = element_text(size = 12)) +
   #scale_fill_manual(values = c('royalblue','darkred')) +
   labs(x='',y='Allele Frequency') 

s1
ggsave(s1, filename = 'paper-II/figures/N_isnv_samples_morethan4.png', width = 15, height = 14, dpi = 300)
ggsave(s1, filename = 'paper-II/figures/N_isnv_samples_morethan4.eps', device = cairo_ps, width = 15, height = 14, dpi = 300)


#H1101Y
nsp3.isnv <- isnv.samp.freq %>% filter(POS == 29033)
ind.date <- unlist(sapply(nsp3.isnv$base, function(i){
  which(all.data.seq$INTERNAL.ID == i)}))

nsp3.isnv <- nsp3.isnv  %>% 
  mutate(date=all.data.seq$DATE.OF.SAMPLE.COLLECTION[ind.date])

nsp3.plot <- nsp3.isnv %>% ggplot(aes(x=info9, y=info2, fill=as.Date(date)))+
  geom_point(shape=21, size=6)+
  geom_line(aes(group=base)) +
  scale_fill_date(date_breaks = '1 month', date_labels = '%b-%y',expand = c(0,0)) +
  theme(panel.background = element_blank(), axis.line = element_line())
nsp3.plot
ggsave(nsp3.plot, filename = 'paper-II/figures/nsp3_6188.png', width = 5, height = 4, dpi = 300)
ggsave(nsp3.plot, filename = 'paper-II/figures/nsp3_6188.eps', device = cairo_ps, width = 5, height = 4, dpi = 300)
