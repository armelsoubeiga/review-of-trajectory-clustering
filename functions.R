# Function to transform a table into long format with a numeric ID
transform_table <- function(data, value_name) {
  df <- as.data.frame(data)
  df$id <- 1:nrow(df)
  df %>%
    pivot_longer(cols = -id, names_to = "key", values_to = value_name)
}

# Round vector of numerics to integer while preserving their sum
sround <- function(x, digits = 0) {
  up <- 10 ^ digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y / up
}

# Synthetic longitudinal data
simul_longitudinal <- function(n, i, n_id) {
  def <- defData(varname = "grpID", formula = paste0("seq(", n_id, ", length.out =", n, ")"))
  def <- defData(def, varname = "nCount", dist = "noZeroPoisson", formula = 6)
  def <- defData(def, varname = "mInterval", dist = "gamma", formula = 30, variance = 0.01)
  def <- defData(def, varname = "vInterval", dist = "nonrandom", formula = 0.07)
  
  df <- genData(n, def)
  df <- addPeriods(df)[,c("id","grpID","period")]
  
  def_ <- defDataAdd(varname = "var_c_3m", dist = "categorical", 
                     formula = paste0(sround(rdirichlet(1, alpha = rep(10, 3)),2), collapse = ";"),
                     variance="Ibuprofen;Paracetamol;Aspirin")
  def_ <- defDataAdd(def_, varname = "var_ord_4m", dist = "categorical", 
                     formula = ".3;.2;.2;.3",
                     variance="NotClinic;ClinicNoTx;TxUncontrolled;TxControlled")
  def_ <- defDataAdd(def_, varname = "var_c_4m", dist = "categorical",
                     formula = paste0(sround(rdirichlet(1, alpha = rep(10, 4)),2), collapse = ";"),
                     variance="Rhinitis;Hypertension;Migraine;Arthritis")
  
  df <- addColumns(def_, df)
  def <- rbind(def,def_)
  def$Class <- paste0('cl_',i)
  df$Class <- paste0('cl_',i)
  rm(def_)
  return(list(def,df))
}

# Create a function to generate ordered categorical values based on probabilities
genOrderedCat <- function(df, def, var) {
  df$new <- NA
  prob_ord <- as.numeric(strsplit(def[def$varname == var, ]$formula, ";")[[1]])
  labels <- as.character(strsplit(def[def$varname == var, ]$variance, ";")[[1]])
  
  for (i in unique(df$id)) {
    df_i <- df[df$id == i, ]
    simax <- max(df_i$period) + 1
    repeat {
      idex_sim <- sample.int(n = length(prob_ord), size = simax, prob = prob_ord, replace = TRUE)
      idex_sim <- sort(idex_sim)
      if(idex_sim[1] == 1 & all(diff(idex_sim) <= 1)){
        break
      }
    }
    df$new[df$id == i] <- idex_sim
  }
  
  df$new <- factor(df$new, levels = 1:length(prob_ord), labels = labels)
  df <- df %>%  select_if(!names(.) %in% var)
  names(df)[names(df) == "new"] <- var
  return(df)
}

# Global synthetic function
synthetic <- function(n) {
  set.seed(43211234)
  n <- n
  n_id <- c(1, 51, 101, 136)
  def <- lapply(seq_along(n), function(i) {return(simul_longitudinal(n[i],i,n_id[i])[[1]])})
  def <- unique(plyr::rbind.fill(def))
  # Categorical var
  df <- lapply(seq_along(n), function(i) {return(simul_longitudinal(n[i],i,n_id[i])[[2]])})
  df <- plyr::rbind.fill(df)
  # Ordinal Categorical var
  df <- lapply(seq_along(n), function(i) {
    df_cl <- df[df$Class == paste0("cl_", i), ]
    df_cl <- genOrderedCat(df = df_cl, def = def, var = "var_ord_4m")
    return(df_cl)})
  df <- plyr::rbind.fill(df)
  # Ordinal Categorical cor
  df$var_ord_3m <- ifelse(df$var_ord_4m %in% c("NotClinic","xControlled"),"Mild",
                          ifelse(df$var_ord_4m == "ClinicNoTx", "Moderate","Severe"))
  # Save
  save(df, def, file = "Data/synthetic/synthetic.Rdata")
}

# Construct sequences
seqlist <- function(data, list_vars) {
  sequences_list <- lapply(list_vars, function(var_name) {
    seqformat(data = data,
              from = "SPELL", to = "STS",
              id = "id", 
              begin = "time", end = "time", 
              status = var_name,
              process = FALSE) %>% 
      seqdef(with.missing = TRUE, right = NA, gaps = NA, left = NA)})
  return(sequences_list)
}

# Construct distance matrices by sequence
seqdistlist <- function(sequence, method) {
  nu <- NULL
  if(method == "TWED") {nu = 0.5}
  if (method %in% c("OM", "TWED", "DHD", "OMslen")) {
    seqdist(sequence, method = method, sm = "INDELSLOG", nu = nu, with.missing = TRUE)
  } else if(method == "NMS") {
    seqdist(sequence, method = method, with.missing = TRUE)
  } else if(method == "CHI2") {
    seqdist(sequence, method = method, with.missing = TRUE,
            step = max(seqlength(sequence)))
  } else {
    seqdist(sequence, method = method, with.missing = TRUE)
  }
}

# Cost matrix function
smatrix <- function(seq_list) {
  smatrix <- lapply(seq_list, function(seq) {
    seqsubm(seq, method = "INDELSLOG", with.missing = TRUE)
  })
  return(smatrix)
}

# Distance function
diff2 <- function(X) return(as.matrix(dist(X, upper = T, diag = T)^2, nrow = nrow(X)))

# Fuzzy to Hard
fuzzy_to_hard <- function(clus) {
  return(clus$clustering)
}

# Evidential to Hard
ev_to_hard <- function(clus, C) {
  mas <- as.data.frame(clus$mass)
  cnames <- function(F) {
    F_ <- F[-c(1,nrow(F)),]
    colnames(F_) <- 1:ncol(F_)
    col_indices <-c("Cl_atypique",apply(F_, 1, function(row) {
      paste0("Cl_", paste0(which(row != 0), collapse = '_'))
    }), "Cl_incertains")
    return(col_indices)
  }
  colnames(mas) <- cnames(clus$F)
  Cluster <- colnames(mas)[apply(mas,1,which.max)]
  return(Cluster)
}


#1 - Combination of states
# This function combines states based on the specified sequence variables and distance function.
cstat <- function(data, seqvars, seqdisfunc){
  # Checking if TWED method is used, and setting parameter nu accordingly
  nu <- NULL
  if(seqdisfunc=="TWED"){nu=0.5}
  # Creating sequence list
  listseq <- seqlist(data, seqvars)
  # Computing multidimensional scaling for sequence data
  MDseq <- seqMD(listseq, what="MDseq", with.missing = TRUE)
  # Computing dissimilarity matrix based on specified distance function
  if(seqdisfunc %in% c("OM","TWED","DHD","OMslen")){
    diss <- seqdist(MDseq, method=seqdisfunc, sm="INDELSLOG", nu=nu, with.missing=TRUE)
  }else if (seqdisfunc=="LCS"){
    diss <- seqdist(MDseq, method="LCS", with.missing=TRUE)
  }else if (seqdisfunc=="CHI2"){
    diss <- seqdist(MDseq, "CHI2",
                    step = max(seqlength(MDseq)), with.missing=TRUE)
  }
  return(list(diss,MDseq))
}

#2 - Combination of costs
# This function combines costs based on the specified sequence variables and distance function.
ccost <- function(data, seqvars, seqdisfunc){
  # Similar to previous function, setting nu parameter if TWED method is used
  nu <- NULL
  if(seqdisfunc=="TWED"){nu=0.5}
  # Creating sequence list
  listseq <- seqlist(data, seqvars)
  # Creating cost matrices for sequence data
  sms <- smatrix(listseq)
  MDcost <- seqMD(listseq, sm=sms, what="cost", with.missing = TRUE)
  MDseq <- seqMD(listseq)
  # Computing dissimilarity matrix based on specified distance function and cost matrix
  diss <- seqdist(MDseq, method=seqdisfunc, sm=MDcost, nu=nu,
                  indel = attr(MDcost,"indel"), with.missing = TRUE)
  return(list(diss,MDcost,MDseq))
}

#3 - Combination of distance matrices
# This function combines distance matrices based on the specified sequence variables and distance function.
cdist <- function(data, seqvars, seqdisfunc){
  # Creating sequence list
  listseq <- seqlist(data, seqvars)
  # Computing distance matrices for each sequence variable and distance function
  listdist <- lapply(1:length(listseq), function(i) {
    seqdistlist(listseq[[i]], seqdisfunc[i])
  })
  # Combining distance matrices
  diss <- Reduce("+", listdist)
  return(list(diss,listseq,listdist))
}

#4 - Combination of clusters 
# This function combines clusters based on the specified sequence variables and distance function.
cclust<- function(data, seqvars, seqdisfunc){
  # Creating sequence list
  listseq <- seqlist(data, seqvars)
  # Computing distance matrices for each sequence variable and distance function
  listdist <- lapply(1:length(listseq), function(i) {
    seqdistlist(listseq[[i]], seqdisfunc[i])
  })
  return(listdist)
}
cclustplus <- function(listdist,listindex,k,lab=NULL,method=NULL,df=NULL){
  # Clustering each distance matrix and combining the results
  list_clustab <- lapply(1:length(listdist), function(i){
    clustering(diss=listdist[[i]],listindex,k,lab,method,df)
  })
  combined_table <- bind_rows(list_clustab)
  # Summarizing the clustering results
  result <- combined_table %>%
    group_by(method, df, index) %>%
    summarise(across(c(cah, pam, fcmdd, ecmdd), ~mean(., na.rm = TRUE)), .groups = "keep")
  result <- result[,c("index","cah","pam","fcmdd","ecmdd","method","df")]
  return(result)
}

# Multiple symPLS by MFA
# This function performs multiple symmetric Partial Least Squares (symPLS) analysis by Multifactorial Analysis (MFA).
multiplesymPLS <- function(matrices_list) {
  # Concatenating matrices
  mat_concat <- do.call(cbind, matrices_list)
  # Extracting number of groups and dimensions
  n_groupe <- sapply(matrices_list, function(x){if(is.vector(x)){return(1)}else{return(ncol(x))}})
  n_ <- length(matrices_list)
  ncp_ <- sum(n_groupe)
  # Performing MFA
  res_ <- FactoMineR::MFA(mat_concat, 
                          group = n_groupe, 
                          type = rep("s",n_),
                          name.group = paste0("matrix",1:n_),
                          ncp=ncp_,
                          graph = FALSE)
  # Computing squared differences between coordinates for each matrix
  result_sum <- 0
  for (i in 1:n_) {
    matrix_name <- paste0("matrix", i)
    result_sum <- result_sum + diff2(res_$separate.analyses[[matrix_name]]$ind$coord)
  }
  result_sum <- (result_sum)^0.5
  return(result_sum)
}

#5 - GIMSA : https://nicolas-robette.github.io/seqhandbook/reference/index.html
# This function performs Group-Informed Multi-Dimensional Scaling with Averaging (GIMSA).
gimsa <- function(data, seqvars, seqdisfunc, keep=NULL){
  # Creating sequence list
  listseq <- seqlist(data, seqvars)
  # Computing distance matrices for each sequence variable and distance function
  listdist <- lapply(1:length(listseq), function(i) {
    seqdistlist(listseq[[i]], seqdisfunc[i])
  })
  # Performing Multi-Dimensional Scaling (MDS) for each distance matrix
  listmds <- lapply(1:length(listdist), function(i) {
    if(is.null(keep)){
      k <- round(nrow(listdist[[i]])*5/100,0)
    }else{
      k <- keep
    }
    cmds <- cmdscale(listdist[[i]], k=k, eig=TRUE)
    cmdsindex <- which.min(seqmds.stress(listdist[[i]], cmds))
    return(cmds$points[,1:cmdsindex])
  })
  if(length(listseq)==2){
    # Performing symmetric PLS if there are only two sequences
    pls <- seqhandbook::symPLS(listmds[[1]],listmds[[2]])
    F <- apply(pls$F,2,scale,center=FALSE)
    G <- apply(pls$G,2,scale,center=FALSE)
    # Computing dissimilarity matrix based on symmetric PLS components
    diss <- (diff2(F)+diff2(G))^0.5
  }else{
    # Performing multiple symPLS for more than two sequences
    diss <- multiplesymPLS(listmds)
  }
  return(diss)
}

# Distance eucludienne
# This function computes Euclidean distances from a dissimilarity matrix.
eucludiss <- function(diss){
  n <- nrow(diss)
  distances <- matrix(0, n , n)
  for(i in 1:n){
    for(j in 1:n){
      distances[i, j]  <- (diss[i] - diss[j])^2
    }
  }
  return(distances)
}

# Clustering with all algorithms
# This function performs clustering using hierarchical clustering (CAH), Partitioning Around Medoids (PAM),
# Fuzzy c-means (FCMdd), and Evidential c-means (ECMdd) algorithms.
clustering <- function(diss,listindex,k,lab=NULL,method=NULL,df=NULL){
  # Hierarchical clustering (CAH)
  cat("cah \n")
  dissclust <- hclust(as.dist(diss), method = "ward.D2")
  cah_clust <- cutree(dissclust, k)
  if(listindex=="cr"){
    cah_cr <- fossil::rand.index(as.numeric(cah_clust), as.numeric(lab))
    cah_asw <- cah_hc <- NA
  }else{
    quality <- wcClusterQuality(as.dist(diss), as.numeric(cah_clust))
    cah_asw <- quality$stats["ASW"]
    cah_hc <- quality$stats["HC"]
    cah_cr <- NA
  }
  
  # Partitioning Around Medoids (PAM)
  cat("PAM \n")
  pam_clust <- cluster::pam(as.dist(diss), k)
  if(listindex=="cr"){
    pam_cr <- fossil::rand.index(as.numeric(pam_clust$clustering), as.numeric(lab))
    pam_asw <- pam_hc <- NA
  }else{
    quality <- wcClusterQuality(as.dist(diss), as.numeric(pam_clust$clustering))
    pam_asw <- quality$stats["ASW"]
    pam_hc <- quality$stats["HC"]
    pam_cr <- NA
  }
  
  # Fuzzy c-means (FCMdd)
  cat("FCMdd \n")
  fuzzy_clust <- fanny(as.dist(diss), k, diss=TRUE, memb.exp = 1.4, maxit=5e+2)
  fuzzy_clust <- fuzzy_to_hard(fuzzy_clust)
  if(listindex=="cr"){
    fuzzy_cr <- fossil::rand.index(as.numeric(fuzzy_clust), as.numeric(lab))
    fuzzy_asw <- fuzzy_hc <- NA
  }else{
    quality <- wcClusterQuality(as.dist(diss), as.numeric(fuzzy_clust))
    fuzzy_asw <- quality$stats["ASW"]
    fuzzy_hc <- quality$stats["HC"]
    fuzzy_cr <- NA
  }
  
  # Evidential c-means (ECMdd)
  cat("ECMdd \n")
  eudiss <- eucludiss(diss)
  ev_clust <- ecmdd(as.matrix(eudiss), c=k, type = 'full', alpha=0.1, beta=1.1, disp=FALSE,
                    gamma = 0, eta = 1)
  ev_clust <- ev_to_hard(ev_clust,k)
  if(listindex=="cr"){
    ev_cr <- fossil::rand.index(as.numeric(as.factor(ev_clust)), as.numeric(lab))
    ev_asw <- ev_hc <- NA
  }else{
    quality <- wcClusterQuality(diss, as.numeric(as.factor(ev_clust)))
    ev_asw <- quality$stats["ASW"]
    ev_hc <- quality$stats["HC"]
    ev_cr <- NA
  }
  
  # Combining and returning clustering results
  results <- cbind.data.frame(
    index =  c(" cr", " aws", " hc"),
    cah = c(cah_cr, cah_asw, cah_hc),
    pam = c(pam_cr, pam_asw, pam_hc),
    fcmdd = c(fuzzy_cr, fuzzy_asw, fuzzy_hc),
    ecmdd = c(ev_cr, ev_asw, ev_hc),
    method=method,
    df=df)
  return(results)
}
