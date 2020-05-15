#-----------------------------------------------------------------------------------------------------------------#
# read validation data
# assuming there are three columns:
#   - "vote" 1 (TRUE) or 0 (FALSE), notifying if the ecosystem was found in the respective range map
#   - "nr_eco": number of ecosystems associated to each voting species
#   - "area": size of the range map for each voting species
#-----------------------------------------------------------------------------------------------------------------#

ids = read.sv('path to file with validation data for target ecosystem')

#-----------------------------------------------------------------------------------------------------------------#
# derive overal accuracy
#-----------------------------------------------------------------------------------------------------------------#

w1 = 1-(ids$nr_eco/(max(ids$nr_eco)+1)) # weight 1
w2 = 1-(ids$area/(max(ids$area[i])+1)) # weight 2
w = w1 + w2 # final weight
v = sum(ids$vote*w) / sum(w) # weighted overal accuracy

# print accuracy
return(v)
