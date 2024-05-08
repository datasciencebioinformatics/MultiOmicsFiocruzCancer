###############################################################################################################################
Network_stages_value_file               <- "/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Network_stages_value.txt" #
###############################################################################################################################
Network_stages_value_file_data     <-read.table(file = Network_stages_value_file, sep = '\t', header = TRUE,fill=TRUE)    #
###########################################################################################################################
# Relevel factor
Network_stages_value_file_data$Variable<-factor(Network_stages_value_file_data$Variable,levels=c("Number of nodes","Number of edges","Avg. number of neighbors","Network diameter","Network radius","Characteristic path length","Clustering coefficient","Network density","Network heterogeneity","Network centralization","Connected components"))

# Standard deviation of the mean as error bar
p <- ggplot(Network_stages_value_file_data, aes(x=Stage, y=Value)) +geom_bar(stat="identity")+ facet_wrap(~Variable, ncol = 5, nrow = 3, scales = "free_y") +theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# FindClusters_resolution
png(filename=paste(output_folder,"Network_panel_","summary.png",sep=""), width = 26, height = 16, res=600, units = "cm")
	p
dev.off()
