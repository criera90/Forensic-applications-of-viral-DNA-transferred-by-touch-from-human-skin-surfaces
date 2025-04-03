# import data
path="/Users/crieraruiz/Library/CloudStorage/OneDrive-UniversityofNebraska-Lincoln/01_SF_lab/Shared_with_committee/230411_request_from_SF_for_NIJ3_proposal/chr_counts_TO/TO_all2.txt"
data=as.matrix(read.table(path,sep='\t',header=T, row.names=1))
data

library(ggplot2)
heatmap(data,scale="column",name="Percent")
as.numeric(data)

install.packages("ComplexHeatmap")
library(ComplexHeatmap)
Heatmap(data, 
        name = "Percent", #title of legend
        column_title = "Chromosome", row_title = "Person",
        row_names_gp = gpar(fontsize = 7) # Text size for row names
)

library(tidyr)
library(dplyr)
library(ggplot2)
library(hrbrthemes)
library(plotly)
library(tibble)

x <- LETTERS[1:20]
y <- paste0("var", seq(1,20))
data <- expand.grid(X=x, Y=y)
data$Z <- runif(400, 0, 5)
data <- data %>% mutate(text = paste0("x: ", x, "\n", "y: ", y, "\n", "Value: ",round(Z,2), "\n", "What else?"))

path="/Users/crieraruiz/Library/CloudStorage/OneDrive-UniversityofNebraska-Lincoln/01_SF_lab/Shared_with_committee/230411_request_from_SF_for_NIJ3_proposal/chr_counts_TO/TO_all_ggplot_scale1-10.txt"
data=as.data.frame(read.table(path,sep='\t',header=T))
ggplot(data, aes(x = person, y = chr, fill = percent)) +
  geom_tile() +
  scale_fill_gradientn(colors = hcl.colors(30, "Oslo")) +
  coord_fixed()
data
