setwd('C:/Users/rr70wedu/Dropbox (iDiv)/AoH/01_documents/')
cinfo = read.csv('./BaseData_info/IUCN_classInfo.csv', stringsAsFactors=FALSE)
cinfo$code = sprintf('%04d', cinfo$code)

setwd('C:/Users/rr70wedu/Dropbox (iDiv)/AoH/01_documents/papers/AoH_map/validation/gbifVal/')

require(ggplot2)
require(extrafont)
loadfonts()


ids = read.csv('IUCN_ecoMap_val-gbif_globalValidation.csv')
ids$accuracy = ids$accuracy*100

p = ggplot(ids, aes(x=resolution, y=accuracy)) + theme_bw() + geom_bar(stat='identity', fill='lightgrey') + 
  theme(axis.text=element_text(size=8, family='Garamond'), axis.title=element_text(size=8, family='Garamond'), 
        panel.background=element_blank(), panel.grid=element_blank()) + 
  scale_x_continuous(breaks=ids$resolution, labels=as.character(ids$resolution)) + 
  scale_y_continuous(limits=c(0,100), breaks=c(0,20,40,60,80,100), expand=c(0,0)) + 
  labs(x='\nResolution (km)', y='Overall Accuracy (%)\n')

ggsave('GlobalAccuracy_gbif.png', p, width=14, height=8, units='cm', dpi=600)


ids = read.csv('IUCN_ecoMap_val-gbif_regionalValidation.csv', stringsAsFactors=F)
ids$resolution = ids$buffer+1
ids$res = as.character(ids$resolution)
ids$res = factor(ids$res, levels=c(1,2,3,4,5,6,7,8,9,10))

ref = read.csv('C:/Users/rr70wedu/Dropbox (iDiv)/AoH/00_data/boundaries/UN-countryRegions.csv', stringsAsFactors=F)
ids$continent = ref$Region.Name[match(ids$region, ref$combSubReg)]

ids = ids[order(ids$region),]
rn0 = c('AUS & NZL', 'Caribbean', 'Central', 'Central', 'Eastern', 'Eastern', 'Eastern', 
        'Melanesia', 'Micronesia', 'Middle', 'Northern', 'Northern', 'Northern', 'Polynesia', 'South-eastern', 
        'South', 'Southern', 'Southern', 'Southern', 'Western', 'Western', 'Western')


ids$region2 = ''
ur = unique(ids$region)
for (r in 1:length(rn0)) ids$region2[ids$region==ur[r]] = rn0[r]

ids$val = ids$val*100
ids$oa = ''
ids$oa[ids$val < 50] = '< 50'
ids$oa[ids$val >=50 & ids$val <60] = '50 - 59'
ids$oa[ids$val >=60 & ids$val <70] = '60 - 69'
ids$oa[ids$val >=70 & ids$val <80] = '70 - 79'
ids$oa[ids$val >=80 & ids$val <90] = '80 - 89'
ids$oa[ids$val >=90] = '90 - 100'
an = rev(c('< 50', '50 - 59', '60 - 69', '70 - 79', '80 - 89', '90 - 100'))
ids$oa = factor(ids$oa, levels=an)

ids$region2 = factor(ids$region2, levels=unique(ids$region2))
ids$continent = factor(ids$continent, levels=sort(as.character(unique(ids$continent))))

cr <- colorRampPalette(c("firebrick2", "khaki2", "forestgreen"))
p = ggplot(ids, aes(x=res, y=region2, fill=oa)) + 
  geom_tile(colour='white', size=0.1) + scale_x_discrete(expand=c(0,0)) + 
  scale_y_discrete(expand=c(0,0)) + 
  scale_fill_manual(values=rev(cr(6)), guide=guide_legend(nrow=1, reverse=TRUE)) + 
  facet_grid(rows=vars(continent), scale='free', space='free', switch='y') + 
  labs(x='Pixel resolution (km)', y='Region', fill='Overall accuracy (%)', family="Helvetica") + 
  theme(aspect.ratio=1,
    strip.placement="outside",
        strip.background=element_blank(),
        legend.key.size=unit(0.3,"cm"),
        axis.text.y=element_text(size=5, family="Helvetica"),
        axis.title.x=element_text(size=6, family="Helvetica"), 
        axis.title.y=element_text(size=6, family="Helvetica"), 
        legend.title=element_text(size=6, family="Helvetica"), 
        legend.text=element_text(size=5, family="Helvetica"), 
        strip.text=element_text(size=6, family='Helvetica'),
        legend.text.align=0.5,
        legend.position='bottom',
        axis.text.x=element_text(size=6, family="Helvetica"), 
        plot.margin = margin(0,0,0,0, unit='mm'),
        strip.text.y=element_text(size=7,family='Helvetica', margin = margin(0,0,0,0, "cm")),
        axis.ticks=element_line(size=0.4),
        panel.grid=element_blank(),
        panel.border=element_blank())

library(grid)
q <- ggplotGrob(p)
lg <- linesGrob(x=unit(c(1,1),"npc"), y=unit(c(1,0),"npc"),gp=gpar(col="black", lwd=1))

for (k in grep("strip-l",q$layout$name)) {
  q$grobs[[k]]$grobs[[1]]$children[[1]] <- lg
}

q$heights[[3]] = unit(0,"in")

#grid.draw(q)
p = ggplotify::as.ggplot(q)

ggsave('regionalAccuracy.png', width=80, units='mm', dpi=600)
ggsave('regionalAccuracy.pdf', width=80, units='mm', dpi=600, device=cairo_pdf)

#-----------------------------------------------------------------------------------#

ids = read.csv('IUCN_ecoMap_val-gbif_yearCount.csv', stringsAsFactors=FALSE)
ids$count2 = ids$count/10000

p = ggplot(ids, aes(x=region, y=count2)) + theme_bw() + geom_bar(stat='identity', fill='lightgrey') + 
  theme(axis.text.x=element_text(size=6, angle=45, hjust=1, family='Garamond'), 
        axis.text.y=element_text(size=6, family='Garamond'), 
        axis.title=element_text(size=8, family='Garamond'), 
        legend.position = 'None',
        panel.grid=element_blank()) + 
  geom_text(aes(label=sprintf("%0.2f", round(count2, digits=2)), vjust=-0.5), 
            position=position_dodge(width=0.9), vjust=-0.5, 
            family='Garamond', size=2) + 
  scale_y_continuous(limits=c(0,400), expand=c(0,0)) + 
  labs(x='\nRegion', y=bquote('Sample Frequency (10'^3*')'))

ggsave('RegionalSampleFreq.png', p, width=15, height=8, units='cm', dpi=600)

#----------------------------------------------------------------------------------#

ids = read.csv('IUCN_ecoMap_val-gbif_yearCount.csv', stringsAsFactors=FALSE)

ids$count2 = ids$count/10000
ids$year2 = factor(as.character(ids$year), levels=as.character(1992:2018))

p = ggplot(ids, aes(x=year2, y=count2)) + theme_bw() + geom_bar(stat='identity', fill='lightgrey') + 
  theme(axis.text.x=element_text(size=5, angle=45, hjust=1, family='Helvetica'), 
        axis.text.y=element_text(size=5, family='Helvetica'), 
        axis.title=element_text(size=6, family='Helvetica'), 
        legend.position = 'None',
        panel.grid=element_blank(), 
        panel.background=element_blank(),
        panel.border=element_blank(), 
        axis.line=element_line(size=0.2),
        axis.ticks=element_line(size=0.2)) + 
  scale_y_continuous(limits=c(0,500), expand=c(0,0)) + 
  labs(x='Year', y=bquote('Sample frequency (10'^3*')'))

ggsave('yearlySampleFreq.png', p, height=41.500, width=86, units='mm', dpi=600)
ggsave('yearlySampleFreq.pdf', p, height=41.500, width=86, units='mm', dpi=600, device=cairo_pdf)
