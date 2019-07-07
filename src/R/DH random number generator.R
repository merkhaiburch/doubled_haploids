#-----------------------------------------------------#
# Title:  Consulation - Randomize Data (Merritt)      #
# Author: Brandon Monier (brandon.monier@sdstate.edu) #
# Date:   05.25.17                                    #
#-----------------------------------------------------#

# Make vector
v <- seq(751,782, 0.5)

# Randomize it
random <- sample(v, length(v), replace = FALSE)

# Assign to data frame

df.rand <- data.frame(plant.id = random)

# Write to .csv (if you want..., UNCOMMENT TO RUN)
# write.csv(df.rand, 'tmp.csv', row.names = FALSE)

write.csv(df.rand, "random.3.csv")
