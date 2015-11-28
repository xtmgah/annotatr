filename <- "~/Desktop/000091.tsv"

#Read in the data
ST000091 <- read.table(filename, sep = "\t", stringsAsFactors = F)

#Transpose, making variable names in the columns
ST000091 <- t(ST000091)

#name each column
colnames(ST000091) <- ST000091[1,]

#eliminate first row
ST000091 <- ST000091[-1,]

#eliminate row names
rownames(ST000091) <- NULL

#Convert to a data frame
ST000091 <- as.data.frame(ST000091)

#omitting columns with missing data
NAcols <- vector()

for (i in 3:ncol(ST000091)){
    ST000091[,i] <- as.numeric(as.character(ST000091[,i]))

    if(is.na(mean(ST000091[,i]))){
        NAcols <- append(NAcols, i)
    }
}

full_data <- ST000091[-NAcols]

#separated data into three categories: controls, deprived, and treated states
control <- full_data[1:9,]
deprived <- full_data[10:17,]
treated <- full_data[18:26,]

#initialize
p_values_control_deprived <- vector()
p_values_control_treated <- vector()
p_values_deprived_treated <- vector()

for (i in 3:ncol(full_data)){
    p_values_control_deprived[i-2] = t.test(control[,i], deprived[,i])$p.value
    p_values_control_treated[i-2] = t.test(control[,i], treated[,i])$p.value
    p_values_deprived_treated[i-2] = t.test(deprived[,i], treated[,i])$p.value
}

metabs_c_d <- names(full_data[which(p_values_control_deprived < 0.05)])
metabs_c_t <- names(full_data[which(p_values_control_treated < 0.05)])
metabs_d_t <- names(full_data[which(p_values_deprived_treated < 0.05)])




#adjusted... yields no postiv
BHY_p_values_control_deprived <- p.adjust(p_values_control_deprived, method = "fdr")
BHY_p_values_control_treated <- p.adjust(p_values_control_treated, method = "fdr")
BHY_p_values_deprived_treated <- p.adjust(p_values_deprived_treated, method = "fdr")













