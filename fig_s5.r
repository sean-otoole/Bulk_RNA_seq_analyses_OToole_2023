library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)
library(stringr)
library(ggbeeswarm)
library(sva)
library(gridExtra)
library(ggpubr)
source('get_path.r')

sampleinfo <- read.delim('data_repository/bulk_rna_seq_counts_and_sample_info/sample_info',header = TRUE, sep = ",")

#get the counts

starting_directory <- getwd()

setwd('data_repository/bulk_rna_seq_counts_and_sample_info/counts')

count_files <- list.files(getwd())
first_count <- count_files[1]

all_samples <- data.frame()

all_samples <- read.table(first_count, header = TRUE, sep = "\t", row.names=1)
colnames(all_samples) <- first_count

count_files <- count_files[1:length(count_files)]

for (file in count_files) {
    seqdata=read.table(file, header = TRUE, sep = "\t", row.names=1)
    colnames(seqdata) <- file
    all_samples <- all_samples[rownames(seqdata),]
    all_samples <- cbind(all_samples, seqdata)
}

all_samples <- all_samples[,c(2:ncol(all_samples))]

colnames(all_samples) <- str_sub(colnames(all_samples),1,nchar(colnames(all_samples))-6)

all_samples <- all_samples[,sampleinfo$SampleName]

all_samples <- na.omit(all_samples)

countdata <- all_samples[order(rowMeans(all_samples), decreasing = TRUE),]
countdata_filtered <- countdata[c(4:(nrow(countdata))),]
y <- DGEList(countdata_filtered)  #exclude the non gene information

group <- sampleinfo$condition
#group <- factor(group)

# Add the group information into the DGEList
y$samples$group <- group

myCPM <- cpm(countdata_filtered)
thresh <- myCPM > 5
keep <- rowSums(thresh) >= 5
summary(keep)
y <- y[keep, ]

logcounts <- cpm(y,log=TRUE)

y <- calcNormFactors(y)

design <- model.matrix(~ 0 + group)

v <- voom(y,design,plot = FALSE)

group_list <- c('1','2','3')

baz <- 2**v$E['Baz1a',]
adam <- 2**v$E['Adamts2',]
agmat <- 2**v$E['Agmat',]

conditions <- sampleinfo$condition
batch <- sampleinfo$Batch
gene.df <- as.data.frame(x = baz)
gene.df <- cbind(adam,agmat,gene.df,conditions,batch)
gene.df <- gene.df[gene.df$conditions %in% c('BazL','BazH','AgmatL','AgmatH','AdamH','AdamL'),]
font.size <- 20

first_plot <- vector()
second_plot <- vector()
third_plot <- vector()

first_plot <- NULL
second_plot < NULL
third_plot <- NULL

plot_variables <- c(first_plot,second_plot,third_plot)


# plotting function
get_the_plots <- function(current_group,input.df,upper_lim){
    plot_list <- list()
    column_names <- colnames(input.df)
    count_list <- c(adam,agmat,baz)
    groups = c('adam','agmat','baz')
    current.df <- input.df[input.df$batch == as.character(current_group),]
    for (i in 1:3){
        plot.df <- current.df[,c(i,4)]
        plot.df$batch <- NULL
        current_gene <- colnames(plot.df)[1]
        colnames(plot.df)[1] <- 'gene'
        plot.df$gene <- plot.df$gene/mean(plot.df[grep('L',plot.df$conditions),]$gene) # get fold changes
        
        n_txt_1 <- table(plot.df$conditions)[1]
        n_name_1 <- names(table(plot.df$conditions)[1])
        n_txt_2 <- table(plot.df$conditions)[2]
        n_name_2 <- names(table(plot.df$conditions)[2])
        sample_sizes <- paste(n_name_1, n_txt_1,n_name_2,n_txt_2, sep = " ")
        
        first_cond <- unique(plot.df$conditions)[1]
        second_cond <- unique(plot.df$conditions)[2]
        comp_list <- list(c(first_cond, second_cond))
        
        p <- ggplot(plot.df,aes(x = conditions, y = gene), color = 'black', size = 1) +
        geom_boxplot(width = 0.4) +
#        stat_summary(geom = "boxplot", width = 0.5, 
#             fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
#             position = "dodge") +
        geom_jitter(color="black", size=2, range = 0.05, height = 0) +
        stat_compare_means(comparisons = comp_list, label = "p.format", method="t.test", paired = FALSE, size = 5) +
        ggtitle(sample_sizes) +
        theme_classic() +
        ylim(0,upper_lim) +
        xlab('') +
        ylab(paste(current_gene,'FC', sep = " ")) +
        scale_x_discrete(limits = rev(levels(as.factor(plot.df$conditions)))) +
        theme(text = element_text(size = font.size)) +
        theme(axis.text = element_text(size = font.size)) +
        theme(legend.position = 'none')
        plot_list[[i]] <- p
        }
return(plot_list)
}

adam_plots <- get_the_plots('1',gene.df, 10)
agmat_plots <- get_the_plots('2',gene.df, 10)
baz_plots <- get_the_plots('3',gene.df, 15)

panels = list(adam_plots[[1]],adam_plots[[2]],adam_plots[[3]],
              agmat_plots[[1]], agmat_plots[[2]], agmat_plots[[3]],
              baz_plots[[1]],baz_plots[[2]],baz_plots[[3]])
    
lay <- rbind(c(1,1,1,2,2,2,3,3,3),
             c(1,1,1,2,2,2,3,3,3),
             c(4,4,4,5,5,5,6,6,6),
             c(4,4,4,5,5,5,6,6,6),
             c(7,7,7,8,8,8,9,9,9),
             c(7,7,7,8,8,8,9,9,9))


setwd(starting_directory)

fig_path = get_path(variables$figures,'fig_s5.pdf')
pdf(fig_path, width = 20, height = 20)
grid.arrange(grobs = panels, layout_matrix = lay)
dev.off()