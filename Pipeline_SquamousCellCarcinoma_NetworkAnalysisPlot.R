###############################################################################################################################
Network_stages_value_file               <- "/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Network_stages_value.txt" #
###############################################################################################################################
Network_stages_value_file_data     <-read.table(file = Network_stages_value_file, sep = '\t', header = TRUE,fill=TRUE)    #
###########################################################################################################################

# Standard deviation of the mean as error bar
p <- ggplot(Network_stages_value_file_data, aes(x=Stage, y=Value)) +geom_bar(stat="identity")+ facet_wrap(~Variable, ncol = 5, nrow = 2, scales = "free_y") +theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
