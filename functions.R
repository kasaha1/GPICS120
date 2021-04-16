

kasa.duplicationRemovalBySD <- function (x)
{
  matrix_data <- as.matrix(x[,-c(1)])
  sd <- apply(matrix_data, 1, sd)
  order_num <- seq(1:nrow(x))
  transformed <- cbind(order_num, sd, x)
  name_list <- colnames(transformed)
  colnames(transformed) <- paste0("var_", seq(1:ncol(transformed)))
  colnames(transformed)[1:3] <- c("order_num", "sd",
                                  "grouped")
  res <-transformed %>% arrange(desc(sd)) %>% group_by(grouped) %>%
    filter(row_number() == 1) %>% ungroup() %>% arrange(order_num)
  colnames(res) <- name_list
  return(res[c(-1,-2)])
}


kasa.transposeMatrix <- function (x, firstColumnName = "sample")
{
  col_names_1 <- t(x[1])
  raw_data <- t(x[-1])
  colnames(raw_data) <- col_names_1
  raw_data <- as.data.frame(raw_data)
  row_name_1 <- row.names(raw_data)
  raw_data <- cbind(row_name_1, raw_data)
  row.names(raw_data) <- NULL
  colnames(raw_data)[1] <- firstColumnName
  raw_data[, 1] <- as.character(raw_data[, 1])
  return(raw_data)
}

kasa.geneMedianCentering <- function (x) 
{
  raw.data <- as.matrix(x[-1])
  median.table <- apply(raw.data, c(1), median, na.rm = T)
  median_centered <- raw.data - median.table
  return(cbind(x[1], median_centered))
}

kasa.geneStandardization <- function (x) 
{
  raw.data <- as.matrix(x[-1])
  sd.table <- apply(raw.data, 1, sd, na.rm = T)
  res.table_1 <- raw.data/sd.table
  res <- cbind(x[1], res.table_1)
  res[is.nan(res)] <- 0
  return(res)
}
kasa.dataCleaning <- function(x){
  res <- list()
  res$classes <- sapply(x,function(y) class(y))
  res$na<- sapply(x,function(y) sum(is.na(y)))
  res$unique <- sapply(x, function(y) length(unique(y)))
  res$dulplicated <- sapply(x, function(y) sum(duplicated(y)))
  res$map <- Amelia::missmap(x, main = "Missing values vs observed")
  return(res)
}

kasaBasicFunctions::kasa.geneStandardization
