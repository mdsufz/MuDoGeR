predict_memory <- function(x) {
    dataset_name <- x[1]
    kmer_33_unique <- as.numeric(x[2])
    kmer_55_unique <- as.numeric(x[3])
    
    if (kmer_55_unique > 0){
        predict_arg <- kmer_55_unique
        chosen_model <- 55
    } else if (kmer_33_unique > 0){
        predict_arg <- kmer_33_unique
        chosen_model <- 33
    } else {
        predict_arg <- NA
        chosen_model <- NA
    }
    x[4] <- chosen_model
    x[5] <- predict_arg
    
    if(!is.na(predict_arg)){
        if(chosen_model == 33){
            if(predict_arg <= max_33_under224){
                used_model <- "model33_mem_under224"
                kmer_33_unique <- predict_arg
                prediction <- predict(model33_mem_under224, data.frame(kmer_33_unique), interval = "prediction")
            }
            else if(predict_arg > max_33_under224 & predict_arg <= max_33_224to360){
                used_model <- "model33_mem_224to360"
                kmer_33_unique <- predict_arg
                prediction <- predict(model33_mem_224to360, data.frame(kmer_33_unique), interval = "prediction")
            }
            else{
                used_model <- "model55_mem_over360"
                kmer_33_unique <- predict_arg
                prediction <- predict(model33_mem_over360, data.frame(kmer_33_unique), interval = "prediction")
            }
        }
        else if(chosen_model == 55){
            if(predict_arg <= max_55_under224){
                used_model <- "model55_mem_under224"
                kmer_55_unique <- predict_arg
                prediction <- predict(model55_mem_under224, data.frame(kmer_55_unique), interval = "prediction")
            }
            else if(predict_arg > max_55_under224 & predict_arg <= max_55_224to360){
                used_model <- "model55_mem_224to360"
                kmer_55_unique <- predict_arg
                prediction <- predict(model55_mem_224to360, data.frame(kmer_55_unique), interval = "prediction")
            }else{
                used_model <- "model55_mem_over360"
                kmer_55_unique <- predict_arg
                prediction <- predict(model55_mem_over360, data.frame(kmer_55_unique), interval = "prediction")
            }
        }
    } else {
        used_model <- NA
        prediction <- data.frame("fit" = NA, "lwr" = NA, "upr" = NA)
    }
    
    x[6] <- used_model
    
    x <- t(as.matrix(x))
    prediction <- as.matrix(prediction)
    predicted_maxvmem <- ceiling(prediction[3])*1000
    results <- cbind(x,prediction,predicted_maxvmem)
    names(results) <- c("dataset","kmer_33_unique","kmer_55_unique","chosen_size","chosen_value","model","fit","lwr","upr","maxvmem_M")
    return(results)
}
