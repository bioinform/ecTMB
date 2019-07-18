

library(data.table)
options(scipen = 999)

args <- commandArgs(trailingOnly=TRUE)


if (length(args)==0){
	stop("Usage: Rscript create.exome.Rdata.R input_file mutation_context_file (e.g. mutation_context_96.txt) output_file_1(e.g. exome_hg19_vep.Rdata)", call.=FALSE)
}


#' Generate output_file_1 (e.g. exome_hg19_vep.Rdata)
data <- read.table(pipe(sprintf("sort -k 1,1 -k2,2n %s",args[1])),sep="\t", stringsAsFactors = F)

sel = data[data[,1] %in% paste0("chr", c(1:22,"X","Y")),]
# sel = sel[order(data[, 2]), ]
# sel = sel[order(data[, 1]), ]
sel = unique(sel[,c(1,2,9,3:8)])
colnames(sel) = paste0("V",1:9)
class(sel$V3) = 'character'


# res = vector("list",24)

# for (i in 1:24){
#   cat("Chromosome",i, "\n")

# 	if ( i ==23){
# 		sel = data[data[,1]=="chrX",]
# 	}else if (i==24){
# 		sel =data[data[,1]=="chrY",]
# 	}else{
# 		sel = data[data[,1]==paste0("chr",i),]
#   }

# 	cat("Number Positions:", nrow(sel), "\n")
# 	cat("Number Genes:", length(unique(sel[,7])), "\n")
# 	sel = sel[order(as.numeric(sel[,2])),]
# 	res[[i]] = unique(sel[,c(1,2,9,3:8)])
# 	colnames(res[[i]]) = paste0("V",1:9)
# 	class(res[[i]]$V3) = 'character'
# }

# exome = do.call(rbind, res)
# exome$V1 = as.character(exome$V1)
# exome$V3 = as.character(exome$V3)
# exome$V8 = as.character(exome$V8)
# exome$V9 = as.character(exome$V9)
exome = sel
save(exome, file=args[3])
write.table(exome, file=sub("Rdata","tsv", args[2]), sep="\t", quote=F, row.names=F, col.names=T)



# #' Generate output_file_2 (e.g. exome_gene_hg19_vep.Rdata)
# mut.context = read.table(args[2])

# preprocess.BM<-function(X,mut.context)
# {

# 	final = NULL
# 	for(i in 1:23)
#   	{
# 		if(nrow(X[[i]]) > 0)
# 		{
#      			cat("Chromosome", i, "\n")
#       			m =  X[[i]]

#       			tmp = data.table(data.frame(gene = m[,7],tag = paste0(m[,2],1),ind = m[,3], pos = m[,1]))
#       			a = data.frame(tmp[, length(pos), by = 'gene,tag,ind'])

#       			tmp = data.table(data.frame(gene = m[,7],tag = paste0(m[,2],4),ind = m[,4], pos = m[,1]))
#       			b = data.frame(tmp[, length(pos), by = 'gene,tag,ind'])

#       			tmp = data.table(data.frame(gene = m[,7],tag = paste0(m[,2],3),ind = m[,5], pos = m[,1]))
#       			c = data.frame(tmp[, length(pos), by = 'gene,tag,ind'])

#       			tmp = data.table(data.frame(gene = m[,7],tag = paste0(m[,2],2),ind = m[,6], pos = m[,1]))
#       			d = data.frame(tmp[, length(pos), by = 'gene,tag,ind'])

#       			final = rbind(final,a[a[,2] %in% mut.context[,3],])
#       			final = rbind(final,b[b[,2] %in% mut.context[,3],])
#       			final = rbind(final,c[c[,2] %in% mut.context[,3],])
#       			final = rbind(final,d[d[,2] %in% mut.context[,3],])

#      		}
#   	}
# 	return(final)
# }

# cat("Start generating output_file_2\n")

# exomgene = preprocess.BM(res,mut.context)
# exomgene = unique(exomgene)
# save(exomgene,file = args[4])

