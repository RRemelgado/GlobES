setwd('C:/Users/rr70wedu/Dropbox (iDiv)/AoH/01_documents/BaseData_info/')
cinfo = read.csv('IUCN_classInfo.csv', stringsAsFactors=FALSE)
cinfo$code = sprintf('%04d', cinfo$code)
gcolour = read.csv('classGroup_colors.csv', stringsAsFactors=F)

setwd('C:/Users/rr70wedu/Dropbox (iDiv)/AoH/01_documents/papers/AoH_map/validation/iucnVal/')

require(ggplot2)
require(extrafont)
loadfonts()

years = 1992:2018
ind = 1:27

#-----------------------------------------------------------------------------#

ids = read.csv('./IUCN_ecoMap_val-iucn_classValidation.csv', stringsAsFactors=F)
ids$code = sprintf('%04d', ids$code)
ids$val = ids$val * 100
ids = ids[ids$val > 0,]
ids$group = cinfo$group[match(ids$code, cinfo$code)]
ids$group = factor(ids$group, levels=gcolour$group)
ids$ghex = gcolour$hex[match(ids$group, gcolour$group)]
ids$class = cinfo$class_2[match(ids$code, cinfo$code)]

ids = ids[!ids$code %in% c('1201','1202','1301'),]
ids$class = factor(ids$class, levels=unique(ids$class[order(unique(ids$code))]))

p = ggplot(ids, aes(x=class, y=val, fill=group)) + theme_bw() + geom_bar(stat='identity') + 
  labs(fill='Ecosystem Group', y='Overall Acuraccy (%)', x='Ecosystem type') + 
  facet_grid(~group, scale='free', space='free') + 
  scale_y_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100), expand=c(0,0)) + 
  scale_fill_manual(values=gcolour$hex[gcolour$group %in% ids$group]) + 
  annotate("segment",x=Inf,xend=-Inf,y=Inf,yend=Inf,color="black",lwd=0.3) + 
  theme(panel.grid=element_blank(),
    strip.background=element_blank(),
    panel.border=element_blank(),
    panel.background=element_blank(),
    legend.position='none',
    strip.text.x =element_text(size=5, face='bold', family='Helvetica', margin=margin(0.1,0.1,0.1,0.1, "cm")),
    strip.text.y=element_text(size=5, face='bold', family='Helvetica', margin=margin(0.1,0.1,0.1,0.1, "cm")),
    axis.text.x=element_text(size=5, angle=45, hjust=1, vjust=1, family="Helvetica"), 
    axis.text.y=element_text(size=5, family='Helvetica'), 
    axis.line = element_line(colour = 'black', size=0.1),
    axis.ticks=element_line(size=0.1),
    axis.title=element_text(size=6, family='Helvetica'))

ggsave('IUCN_classValidation.png', p, width=183, height=60, units='mm', dpi=600)
ggsave('IUCN_classValidation.pdf', p, width=183, height=60, units='mm', device=cairo_pdf, dpi=600)
