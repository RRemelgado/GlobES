#-----------------------------------------------------------------------------------------------------------------#
# read validation data
# assuming there are four columns:
#   - "nr_votes" number of TRUE votes for each species per resolution
#   - "nr_eco": number of ecosystems associated to each voting species per resolution
#   - "nr_pts": size of the range map for each voting species per resolution
#   - resolution: target resolution
#-----------------------------------------------------------------------------------------------------------------#

ids = read.csv('file containing validation data for one region')

#-----------------------------------------------------------------------------------------------------------------#
# derive overal accuracy for target region, for each resolution
#-----------------------------------------------------------------------------------------------------------------#

unique.resolution = unique(ids$resolution)

ods = do.call(rbind, lapply(unique.resolution, function(r) {
  
  i = which(ids$resolution == r) # find target observations
  w1 = 1-(ids$nr_eco[i]/(max(ids$nr_eco[i])+1)) # weigth 1
  w2 = ids$nr_pts[i]/max(ids$nr_pts[i]) # weight 2
  w = w1 + w2 # final weight
  v = sum(ids$nr_votes[i]*w) / sum(ids$nr_pts[i]*w) # weighted overal average
  df = data.frame(resolution=b, oa=v, stringsAsFactors=F)
}))

