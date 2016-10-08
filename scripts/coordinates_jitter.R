
coord <- read.table("coordinates", header=T, sep="\t")

ind <- coord[1]
lat <- coord[2]
long <- coord[3]

n.loc <- nrow(coord)

jitter <- 0.1

jit <- runif(n.loc, min=-jitter, max=jitter)
lat.jit <- lat + jit

jit <- runif(n.loc, min=-jitter, max=jitter)
long.jit <- long + jit

coord.jit <- data.frame(ind, lat.jit, long.jit)

write.table(coord.jit, file="coordinates.jit", col.names = T, row.names = F, quote = F)
