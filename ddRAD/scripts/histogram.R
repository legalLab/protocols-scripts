#!/usr/bin/env Rscript

library(ggplot2)

rm(list=ls())

re1 <- 'SdaI'
re2 <- 'Csp6I'
re <- paste(re1, re2, sep='_')
sp <- 'aotus_nancymaae'

#RAD <- read.table('aotus_nancymaae_RADtag.extracted_length', header = FALSE)
ddRAD <- read.table(paste('../cut_results/', re, '/', sp, '_ddRADtag.extracted_length', sep=''), header = FALSE)

ddRAD_log <- log10(ddRAD)

svg(paste(sp, '_', re, '.svg', sep=''))
ggplot(data=ddRAD_log, aes(ddRAD_log$V1)) + 
  geom_histogram(breaks=seq(0, 5, by=.1),
                 #col="red",
                 aes(fill=..count..)) +
  scale_fill_gradient("count", low = "red", high = "green") +
  #geom_vline(xintercept=log10(350)) +
  #geom_vline(xintercept=log10(450)) +
  ggtitle(paste(re1, re2, 'cut', sep=' ')) +
  labs(x="log10 fragment length", y="fragment count") +
  theme(plot.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=32, hjust=.5)) +
  theme(axis.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=16))
dev.off()

