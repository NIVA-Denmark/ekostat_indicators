
flist<-c("SE592600-181600","SE631500-190270","SE575095-164630","SE574750-164500","SE574500-164500","SE575335-165000","SE572810-164500","SE584520-172495","SE584520-172495","SE624800-181030","SE592500-191750","SE584333-172895","SE572000-163835","SE590635-182120","SE580950-170601","SE581740-170260","SE583370-165290","SE582820-165920","SE582600-165680")
df <- df %>% filter(WB %in% flist)

dflist <- df %>% split(.$WB)
cat(paste0("Time now: ",Sys.time(),"\n"))

db <- dbConnect(SQLite(), dbname="output/ekostat_coast_map.db")

dflist2 <- dflist %>% map(~ Assessment(.,nsim = 10, IndList,df.bounds,df.bounds.hypox,df.bathy,df.indicators,df.variances,db))

dbDisconnect(db)

cat(paste0("Time now: ",Sys.time(),"\n"))
save(dflist2,file="output/output_map.Rda")
cat(paste0("           time elapsed: ",Sys.time() - start_time,"\n"))


