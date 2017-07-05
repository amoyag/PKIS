# overrepresentation test of the targets of each drug on each of the FunFams


## First used in exp2016-08-01. Implements BH correction for multiple testing and it uses the total number of targets instead of the number of human proteins as N, to calculate the expected probability of a protein to be a relative of a FunFam (see below).

## Needs local access to ChEMBL and access to CATH


library(dplyr)
library(splitstackshape)
library(XLConnect)
library(RMySQL)
library(stats)





################ Functions #######################


##### COde to connect to CATH and perform a query
library(RJDBC)
library(data.table)

### Need to open a ssh tunnel to sinatra in the terminal for this to work properly
# 
### Also need to have the host sec1ssh defined in ~/.ssh/config
#  
#  
#  
#  
#  

ORACLE_JDBC_PATH = "" # path to your ojdbc6.jar file
drv <- JDBC("oracle.jdbc.OracleDriver",classPath = ORACLE_JDBC_PATH)
con <- dbConnect(drv, "") # Ask Ash for the details of the connection


## Example of query. Implemented in the function below
# target = "'I9QCD8'"
# 
# sql <- paste("select ff.superfamily_id, ff.funfam_number from CATH_DATA.funfam_member ff, CATH_DATA.CATHHMM_TO_MDA mda, CATH_DATA.UNIPROT_DESCRIPTION u where u.sequence_md5 = mda.id and u.sequence_md5 = ff.sequence_md5 and u.uniprot_acc =",target, sep="")
# 
# result <- data.table(dbGetQuery(con, sql))


### ChEMBL queries

# Connect to the DB (local)
chembldb <- dbConnect(MySQL(), user="YOUR_USER", password='YOUR_PASWORD', dbname="YOUR_CHEMBL_DB", host="localhost")
on.exit(dbDisconnect(chembldb))



get_uniprot <- function(superfamily_idfunfam_id){
        require(dplyr)
        superfamily_id <- strsplit(superfamily_idfunfam_id,split = "-")[[1]][1] # split for the FunFHMMER FunFams
        funfam_id <- strsplit(superfamily_idfunfam_id,split = "-")[[1]][2] # split for the FunFHMMER FunFams
        superfamily_id <- paste("'",superfamily_id,"'",sep = "")
        require(data.table)
        as.character(superfamily_id)
        sql <- paste("select u.uniprot_acc, ff.superfamily_id, ff.funfam_number from cath_v4_1_0.funfam_member ff, cath_v4_1_0.CATHHMM_TO_MDA mda, cath_v4_1_0.UNIPROT_DESCRIPTION u where u.sequence_md5 = mda.id and u.sequence_md5 = ff.sequence_md5 and u.taxon_id = 9606 and ff.superfamily_id =",superfamily_id, sep="")
        #return(data.table(dbGetQuery(con, sql)))
        kk <-data.table(dbGetQuery(con, sql))
        result <- filter(kk, FUNFAM_NUMBER == funfam_id)
        return(result$UNIPROT_ACC)
}

### Note: query to get data from FunFHMMER FunFams. To use the old DFX FunFams, change cath_v4_1_0 to cath_data
### call the function: get_uniprot("3.90.1640.10", 2660)




################# End of Functions ######################



# Establish the hypothesised probability of success i.e. the probability of observing a protein as a relative of a FunFam Pff = n(ff)/N ; n(ff)= relatives of the FunFam ff N=total number of relatives of all FunFams
# Update. N is the number of targets

sql1 <- "select u.uniprot_acc, ff.superfamily_id, ff.funfam_number from cath_v4_1_0.funfam_member ff, cath_v4_1_0.CATHHMM_TO_MDA mda, cath_v4_1_0.UNIPROT_DESCRIPTION u where u.sequence_md5 = mda.id and u.sequence_md5 = ff.sequence_md5 and u.taxon_id = 9606"

### Note: query to get data from FunFHMMER FunFams. To use the old DFX FunFams, change cath_v4_1_0 to cath_data

all_mda <- data.table(dbGetQuery(con, sql1))

# all_mda <- all_mda %>% mutate(FunFam_ID = paste(SUPERFAMILY_ID, FUNFAM_NUMBER, sep = ".FF")) ## This is used with the DFX FunFams
all_mda <- all_mda %>% mutate(FunFam_ID = paste(SUPERFAMILY_ID, FUNFAM_NUMBER, sep = "_")) ## This is used with the FunFHMMER FunFams
all_mda <- all_mda %>% select(FunFam_ID, UNIPROT_ACC)



# Get the drug-target data from ChEMBL


# query to get all human targets from ChEMBL
mysql.targets = "SELECT DISTINCT cs.accession FROM component_sequences cs JOIN target_components tc  ON tc.component_id = cs.component_id WHERE db_source IN('SWISS-PROT', 'TREMBL') AND ORGANISM = 'Homo sapiens'"

human.targets <- dbGetQuery(chembldb, mysql.targets)


## Get drug-target data

chembl_drug.for.target <- function(protein){
        queryprotein <- paste("'",protein,"'", sep="")
        mysql.drugs = "SELECT DISTINCT md.CHEMBL_ID, cs.accession
                        FROM component_sequences cs 
                        JOIN target_components tc ON cs.component_id=tc.component_id 
                        JOIN target_dictionary td ON tc.tid=td.tid 
                        JOIN assays ass ON ass.tid = tc.tid
                        JOIN activities act  ON ass.assay_id=act.assay_id
                        JOIN molecule_dictionary md ON act.molregno=md.molregno
                        JOIN compound_structures ms ON act.molregno=ms.molregno 
                        WHERE act.pchembl_value >= 6
                        AND(data_validity_comment IS NULL OR data_validity_comment = 'manually validated')
                        AND md.max_phase = 4
                        AND assay_type = 'B'
                        AND relationship_type = 'D'
                        AND target_type = 'SINGLE PROTEIN'
                        AND standard_relation = '='
                        AND accession = "
        sql <- paste(mysql.drugs, queryprotein, sep = "")
        return(dbGetQuery(chembldb, sql))
}  

#chembl_drug.for.target("P35968")

drug_target <- as.data.frame(do.call(rbind, lapply(human.targets$accession, chembl_drug.for.target)))

names(drug_target) <- c("CHEMBL_ID", "UNIPROT_ACC")
#drug_target$target.in.ff <- sapply(as.character(drug_target$UNIPROT_ACC), get_mda)




#targets_mda <- filter(all_mda, UNIPROT_ACC %in% unique(drug_target$UNIPROT_ACC))

targets_mda <- unique(merge(drug_target, all_mda, by = "UNIPROT_ACC") %>% dplyr::select(UNIPROT_ACC, FunFam_ID))

### Compute the expected probability that a target is a relative of a FunFam by Pff = nff/N where nff is the number of targets in the FunFam and N is the total number of targets.


N = length(unique(targets_mda$UNIPROT_ACC))
kk <- targets_mda %>% group_by(FunFam_ID) %>% summarise(Pff = length(UNIPROT_ACC)/N)

targets_mda <- merge(targets_mda, kk, by= "FunFam_ID")









drug_ff <- merge(drug_target, targets_mda, by="UNIPROT_ACC")

names(drug_ff) <- c("UNIPROT_ACC", "CHEMBL_ID",   "target.in.ff",   "Pff" )
write.table(drug_ff, file = "drug-ff.tsv", quote = F,row.names = F,sep = "\t")

# Load Toy model in results/exp2016-06-24
#drug_ff <- read.table("drug-target-ff.tsv", header = T, sep = "\t")



drug_ff <- drug_ff %>% group_by(CHEMBL_ID) %>% mutate(no.targets = length(unique(UNIPROT_ACC))) #no. targets per drug (trials)

# function to calculate the number of successes, i.e. the number of targets of a drug that are in each funfam


compute.success <- function(drug){
        list.of.targets <- unique((filter(drug_target, CHEMBL_ID==drug))$UNIPROT_ACC)
        list.of.funfams <-  unique((filter(drug_ff, CHEMBL_ID==drug))$target.in.ff)
        return(filter(drug_ff, target.in.ff %in% list.of.funfams) %>% group_by(target.in.ff) %>% mutate("success" = length(list.of.targets[list.of.targets %in% UNIPROT_ACC])) %>% filter(CHEMBL_ID==drug))
}

# Add column with the number of targets in each funfam, per drug

#drug_ff <- drug_ff %>% group_by(CHEMBL_ID) %>% mutate(no.targets.in.ff = length(unique(target.in.ff))) # number of FF in which the targets are distributed

drug_ff <- unique(as.data.frame(do.call(rbind,lapply(drug_ff$CHEMBL_ID, compute.success))))



#drug_ff <- drug_ff %>% mutate(over.under = Pff * no.targets) %>% filter(no.targets.in.ff > over.under)
drug_ff <- drug_ff %>% mutate(over.under = Pff * no.targets)

overrep <- function(success,over.under) {
        if(success > over.under) {return("Over")}
        if(success < over.under) {return("Under")}
}

drug_ff$Overrepres <- apply(drug_ff[,c("success", "over.under")], 1, function(x) overrep(x["success"], x["over.under"]))


binomial.pval <- function(succ,trials,prob){return(round((binom.test(succ,trials,prob))$p.value,digits = 6))}
drug_ff$p.val <- apply(drug_ff[,c("success", "no.targets", "Pff")], 1, function(x) binomial.pval(x["success"], x["no.targets"], x["Pff"]))


##BH corrected pval

drug_ff$BH.pVal <- p.adjust(drug_ff$p.val, method = "BH")


write.table(drug_ff, file = "drug-ff_overreptest.tsv", quote = F,row.names = F,sep = "\t")



drug_ff.over <- filter(drug_ff, Overrepres == "Over" & BH.pVal < 0.001)

write.table(drug_ff.over, file = "drug-ff_overreptest_over.tsv", quote = F,row.names = F,sep = "\t")

drug_ff.assoc <- unique(select(drug_ff.over, CHEMBL_ID, target.in.ff, p.val,BH.pVal))

write.table(drug_ff.assoc, file = "drug-ff.tsv", quote = F,row.names = F,sep = "\t")


write.table(targets_mda, file = "hstargets_ffmda.tsv", quote = F,row.names = F,sep = "\t")



ff.per.drug <- drug_ff.assoc %>% group_by(CHEMBL_ID) %>% summarise(ff.per.drug = length(unique(target.in.ff)))
drugs.per.ff <- drug_ff.assoc %>% group_by(target.in.ff) %>% summarise(drugs.per.FF = length(unique(CHEMBL_ID)))

targets.per.drug <- drug_target %>% group_by(CHEMBL_ID) %>% summarise(targets = length(unique(UNIPROT_ACC)))
drugs.per.target <- drug_target %>% group_by(UNIPROT_ACC) %>% summarise(drugs = length(unique(CHEMBL_ID)))

targets.per.ff <- targets_mda %>% group_by(FunFam_ID) %>% summarise(targets = length(unique(UNIPROT_ACC)))
targets.per.ff <- merge(targets.per.ff, targets_mda, by= "FunFam_ID") %>% dplyr::select(FunFam_ID,targets,Pff)


ff.per.target <- targets_mda %>% group_by(UNIPROT_ACC) %>% summarise(targets = length(unique(FunFam_ID)))

write.table(ff.per.drug, file = "ff_per_drug.tsv", quote = F,row.names = F,sep = "\t")
write.table(drugs.per.ff, file = "drugs_per_ff.tsv", quote = F,row.names = F,sep = "\t")

write.table(targets.per.drug, file = "targets.per.drug.tsv", quote = F,row.names = F,sep = "\t")
write.table(drugs.per.target, file = "drugs_per_target.tsv", quote = F,row.names = F,sep = "\t")

write.table(ff.per.target, file = "ff.per.target.tsv", quote = F,row.names = F,sep = "\t")


# Extra stuff. Brief data analysis

# kk <- (all_mda %>% group_by(FunFam_ID) %>% mutate(prots.in.FF = length(unique(UNIPROT_ACC))))$prots.in.FF # number of proteins per FunFam in the general mda

# ll <- (targets_mda %>% group_by(FunFam_ID) %>% mutate(prots.in.FF = length(unique(UNIPROT_ACC))))$prots.in.FF # number of proteins per FunFam in the drugs targets mda

# useful for statistical (median, range, iqr) calculations






