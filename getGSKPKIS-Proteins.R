# test to get ChEMBL stuff from R

library(RMySQL)
library(igraph)
library(dplyr)
library(ggplot2)

# Connect to the DB (local)
chembldb <- dbConnect(MySQL(), user="root", dbname="chembl_20", host="localhost")

on.exit(dbDisconnect(chembldb))

query50 = "SELECT DISTINCT md.CHEMBL_ID, cs.accession, td.chembl_id, act.published_value
        FROM molecule_dictionary md
	JOIN compound_records cr ON cr.molregno = md.molregno
	JOIN activities act ON act.molregno = cr.molregno
	JOIN assays ass ON ass.assay_id = act.assay_id
	JOIN target_dictionary td ON td.tid = ass.tid
	JOIN target_components tc ON tc.tid = td.tid
	JOIN component_sequences cs ON cs.component_id = tc.component_id
        JOIN source src ON src.src_id=ass.src_id
        WHERE src_short_name = 'GSK_PKIS'
        AND ass.description LIKE '%0.1 uM%'
        AND act.published_value >= 50"

query80 = "SELECT DISTINCT md.CHEMBL_ID, cs.accession, td.chembl_id, act.published_value
        FROM molecule_dictionary md
        JOIN compound_records cr ON cr.molregno = md.molregno
	JOIN activities act ON act.molregno = cr.molregno
	JOIN assays ass ON ass.assay_id = act.assay_id
	JOIN target_dictionary td ON td.tid = ass.tid
	JOIN target_components tc ON tc.tid = td.tid
	JOIN component_sequences cs ON cs.component_id = tc.component_id
        JOIN source src ON src.src_id=ass.src_id
        WHERE src_short_name = 'GSK_PKIS'
        AND ass.description LIKE '%0.1 uM%'
        AND act.published_value >= 80"




df50 <- dbGetQuery(chembldb, query50)
names(df50) <- c('compound', 'target_UP', 'target_ChEMBLID', 'activity')

df80 <- dbGetQuery(chembldb, query80)
names(df80) <- c('compound', 'target_UP', 'target_ChEMBLID', 'activity')

######### Build the graph before writing out and use simplify() to remove multiple edges
# There are several cases where the activity assay is repeated for the same pair of pkis-target, and they have different activity values
# filter(df50, target_UP == 'P00519' & compound=='CHEMBL1909383')
# compound target_UP target_ChEMBLID activity
# 1 CHEMBL1909383    P00519      CHEMBL1862    92.14
# 2 CHEMBL1909383    P00519      CHEMBL1862    93.90
# 3 CHEMBL1909383    P00519      CHEMBL1862    92.26
# 4 CHEMBL1909383    P00519      CHEMBL1862    94.17
# 5 CHEMBL1909383    P00519      CHEMBL1862    84.57
# 6 CHEMBL1909383    P00519      CHEMBL1862    93.30
# 7 CHEMBL1909383    P00519      CHEMBL1862    81.67

# Simplify the graph to get the mean activity in  cases of multiple activity

temp <- select(df50, target_UP, compound, activity)
names(temp) <- c('target_UP', 'compound', 'weight') #using the name weight creates a weighted graph
g50 <- simplify(graph.data.frame(temp, directed=F),remove.multiple = T, edge.attr.comb = 'mean')

temp <- select(df80, target_UP, compound, activity)
names(temp) <- c('target_UP', 'compound', 'weight') #using the name weight creates a weighted graph
g80 <- simplify(graph.data.frame(temp, directed=F),remove.multiple = T, edge.attr.comb = 'mean')

rm(temp)


# Write out the gskpkis-protein table

tmp <- get.data.frame(g50)
names(tmp) <- c('target_UP', 'compound', 'activity')
write.table(tmp, "gskpkis-protein-in50.txt", sep = '\t', quote = F,row.names = F)
rm(tmp)

tmp <- get.data.frame(g80)
names(tmp) <- c('target_UP', 'compound', 'activity')
write.table(tmp, "gskpkis-protein-in80.txt", sep = '\t', quote = F,row.names = F)
rm(tmp)

####### Inhibition 50% #######

# graph partition according to node types
drugs50 <- V(g50)[134:length(V(g50))]
targets50 <- V(g50)[1:133]

# degree distributions
drugdegree50<-degree(g, v=drugs50)
targetdegree50 <- degree(g, v=targets50)

dd50<-data.frame(table(drugdegree50))
#ggplot(dd, aes(x=drugdegree, y=Freq)) + geom_line() #doesn't work. Fix
td50 <- data.frame(table(targetdegree50))
par(mfrow=c(1,2))

######## Inhibition 80% ###########


# graph partition according to node types
drugs80 <- V(g80)[83:length(V(g80))]
targets80 <- V(g80)[1:82]

# degree distributions
drugdd80<-degree.distribution(g80, v=drugs80)
targetdd80 <- degree.distribution(g80, v=targets80)

drugdegree80<-degree(g80, v=drugs80)
targetdegree80 <- degree(g80, v=targets80)

dd80<-data.frame(table(drugdegree80))
td80 <- data.frame(table(targetdegree80))
#ggplot(dd, aes(x=drugdegree, y=Freq)) + geom_line() #doesn't work. Fix
par(mfrow=c(1,2))
plot(dd80)
plot(td80)

cat("Fraction of specific drugs")
round((length(drugdegree80[drugdegree80 ==1])/length(drugdegree80))*100, digits = 2)

cat("Fraction of promiscuous drugs")
round((length(drugdegree80[drugdegree80 >1])/length(drugdegree80))*100, digits = 2)


cat("Fraction of promiscuous targets")
round((length(targetdegree80[targetdegree80 >1])/length(targetdegree80))*100, digits = 2)


