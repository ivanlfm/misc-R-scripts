vector <- as.vector(RES_clado_events_tables[[b]][1])
for(i in 2:length(RES_clado_events_tables[[b]])){
    if(names(RES_clado_events_tables[[b]][i])!="daughter_nds" & names(RES_clado_events_tables[[b]][i])!="SUBdaughter_nds"){
		vector <- cbind(vector, as.vector(RES_clado_events_tables[[b]][i]))
	}
}
write.csv(vector, paste("BSM_table_rep_",r,"map_",b,".csv",sep=""))