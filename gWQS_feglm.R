##############################################################################################################################################
#### Spatial and Racial/ethnic Disparities in Cardiovascular Mortality Attributable to PM2.5 Components in the Contiguous United States   ####
#### Cleaned R code for the gWQS_feglm function                                                                                           ####
#### Ying Hu, Lingzhi Chu, Stefano Renzetti, Emma Zang, Ijeoma Opara, Yuan Lu, Erica S. Spatz, Harlan M. Krumholz, Kai Chen               ####
#### Yale University, New Haven, CT                                                                                                       ####
#### Feb. 20, 2025                                                                                                                        ####
##############################################################################################################################################

require(fixest); require(dplyr); require(gWQS); require(splines); require(parallel); require(multcomp); require(future); require(future.apply)

gWQS_feglm <- function(model.data = NULL,
                       usecores = 30, rh = 100, b = 100, nq = 10, validation = 0.3, lambda_value = 1, seed = 1000,
                       mix_name = NULL, covar_name = NULL, fix_name = NULL, cluster_name = NULL){
  
  # Parameter settings
  parameter.setting <- list(usecores = usecores, rh = rh, b = b, nq = nq, validation = validation, lambda_value = lambda_value, seed = seed,
                mix_name = mix_name, covar_name = covar_name, fix_name = fix_name, cluster_name = cluster_name)
  
  # WQSFeglm function
  wqsGaussFeglm <- function(initp, formula, kw, bXm, bY, bQ, kx, lambda, fix, cluster_name) {
    initp[1] <- abs(initp[1])
    w <- initp[(kx + 2):(kx + kw + 1)]
    pen <- sum(abs(w))
    w <- w^2
    w <- w / sum(w) 
    bXm[, "wqs"] <- as.numeric(as.matrix(bQ) %*% as.vector(w))
    data <- bXm
    data$y <- bY
    data[,c(colnames(fix))] <- fix
    start <- initp[1:(kx+1)]
    feglm_model <- feglm(formula, 
                         family = "gaussian",
                         cluster = cluster_name,
                         data = data,
                         start = start)
    residuals <- residuals(feglm_model)
    rss <- sum(residuals^2)
    return(list(rss = rss + pen * lambda, model = feglm_model, weights = w))
  }
  
  # Set up future parallel strategy
  plan(multisession, workers = usecores)
  options(future.globals.maxSize= 3*1024^3)
  
  # Define optimization logic as a standalone function
  optimize_weights <- function(training_data, covar_name, mix_name, fix_name, initp, formula, kw, kx, lambda_value, cluster_name, b, rh) {
    
    results <- future_lapply(1:b, function(j) {

      if (b == 1) {
        training_data_boot <- training_data
      } else {
        sampled_GEOID <- sample(unique(training_data$GEOID), size = length(unique(training_data$GEOID)), replace = TRUE)
        training_data_boot <- do.call(rbind, lapply(seq_along(sampled_GEOID), function(k) {
          id <- sampled_GEOID[k]
          subset_data <- training_data[training_data$GEOID == id, ]
          repeat_count <- sum(sampled_GEOID[1:k] == id)
          if (repeat_count > 1) subset_data$GEOID <- paste0(id, "_", repeat_count)
          subset_data
        }))
      }
      
      train_bXm <- training_data_boot[, covar_name]
      train_bY <- training_data_boot$y
      train_bQ <- training_data_boot[, mix_name]
      train_fix <- training_data_boot[, fix_name]
      
      train_weights <- rep(NA, length(mix_name))
      opt_result <- optim(par = initp,
                          fn = function(par) {
                            wqsGaussFeglm(par, formula, kw = kw, bXm = train_bXm, bY = train_bY, bQ = train_bQ, kx = kx, lambda = lambda_value, fix = train_fix, cluster_name = cluster_name)$rss
                          },
                          method = "BFGS",
                          control = list(reltol = 1e-12))
      optimized_result <- wqsGaussFeglm(opt_result$par, formula, kw = kw, bXm = train_bXm, bY = train_bY, bQ = train_bQ, kx = kx, lambda = lambda_value, fix = train_fix, cluster_name = cluster_name)
      
      list(train_weights = optimized_result$weights, train_model = optimized_result$model)
    }, future.seed = TRUE)
    
    all_b_weights <- as.data.frame(do.call(rbind, lapply(results, function(res) res$train_weights)))
    colnames(all_b_weights) <- mix_name
    mean_weight_rh <- colMeans(all_b_weights)
    all_b_weights$RH <- rh
    
    all_b_wqs <- as.data.frame(do.call(rbind, lapply(results, function(res) res$train_model$coefficients)))["wqs"]
    all_b_weights <- bind_cols(all_b_wqs, all_b_weights)
    all_b_weights_p <- all_b_weights %>% filter(all_b_weights[,1] >= 0)
    all_b_weights_n <- all_b_weights %>% filter(all_b_weights[,1] < 0)
    mean_weight_rh_p <- colMeans(all_b_weights_p[,mix_name])
    mean_weight_rh_n <- colMeans(all_b_weights_n[,mix_name])
    
    return(list(all_b_weights = all_b_weights, mean_weight_rh = mean_weight_rh, mean_weight_rh_p = mean_weight_rh_p, mean_weight_rh_n = mean_weight_rh_n))
  }
  
  # Set initial weights
  kx <- length(covar_name)
  kw <- length(mix_name)
  formula0 <- as.formula(paste("y ~ ", paste(covar_name, collapse = " + "), "|", paste(fix_name, collapse = " + ")))
  formula <- as.formula(paste("y ~ wqs +", paste(covar_name, collapse = " + "), "|", paste(fix_name, collapse = " + ")))
  initial_fit <- feglm(formula0, data = model.data, family = gaussian())
  initial_params <- initial_fit$coefficients
  initp <- c(wqs = 0, initial_params, c(rep(1/kw, kw)))
  
  # Processing the data for analysis
  q_f <- gwqs_rank(model.data, mix_name, nq)
  bQ <- q_f$Q
  data.fit <- model.data[ ,fix_name]
  data.fit[,mix_name] <- bQ
  data.fit[,covar_name] <- model.data[, covar_name]
  data.fit$y <- model.data$y
  unique_geoids <- unique(data.fit$GEOID)
  geoids_count <- length(unique_geoids)
  
  start.time <- Sys.time()
  mean_weights <- data.frame()
  all_weights <- data.frame()
  all_rss <- data.frame()
  all_wqs <- data.frame()
  set.seed(seed)
  
  all_results <- future_lapply(1:rh, function(i) {
    
    # Split data into training and validation sets
    validation_geoids <- sample(unique_geoids, size = validation * geoids_count)
    validation_data <- data.fit %>% filter(GEOID %in% validation_geoids)
    training_data <- data.fit %>% filter(!GEOID %in% validation_geoids)
    
    # Call the optimization function
    result_weights <- optimize_weights(
      training_data = training_data,
      covar_name = covar_name,
      mix_name = mix_name,
      fix_name = fix_name,
      initp = initp,
      formula = formula,
      kw = kw,
      kx = kx,
      lambda_value = lambda_value,
      cluster_name = cluster_name,
      b = b,
      rh = i
    )
    
    mean_weight_rh <- as.data.frame(bind_rows(result_weights$mean_weight_rh_p, result_weights$mean_weight_rh_n, result_weights$mean_weight_rh))
    rownames(mean_weight_rh) <- c("p", "n", "all")
    
    valid_bXm <- validation_data[, covar_name]
    valid_bY <- validation_data$y
    valid_bQ <- validation_data[, mix_name]
    valid_fix <- validation_data[, fix_name]
    all_rh_wqs <- data.frame()
    for (d in 1:3){
      if (is.na(sum(mean_weight_rh[d,]))){
        all_rh_wqs[d,] <- NA
      } else {
        weight_vector <- as.vector(as.numeric(mean_weight_rh[d,]))
        valid_bXm[, "wqs"] <- as.numeric(as.matrix(valid_bQ) %*% weight_vector)
        valid_data <- valid_bXm
        valid_data$y <- valid_bY
        valid_data[,c(colnames(valid_fix))] <- valid_fix
        feglm_model <- feglm(formula, 
                             family = "gaussian",
                             cluster = cluster_name,
                             data = valid_data)
        all_rh_wqs_d <- feglm_model[["coeftable"]]["wqs", ]
        all_rh_wqs <- bind_rows(all_rh_wqs, all_rh_wqs_d)
      }
    }
    all_b_weights <- result_weights$all_b_weights
    list(all_rh_wqs = all_rh_wqs, all_b_weights = all_b_weights, mean_weight_rh = mean_weight_rh)
  }, future.seed = TRUE)
  
  Output <- list()
  all_b_weights <- as.data.frame(do.call(rbind, lapply(all_results, function(res) res$all_b_weights)))
  all_rh_wqs <- list()
  mean_weight_rh <- list()
  wqs_np <- c("p", "n", "all") 
  for (d in 1:3){
    all_rh_wqs[wqs_np[d]] <- as.data.frame(do.call(rbind, lapply(all_results, function(res) res$all_rh_wqs$Estimate[d])))
    mean_weight_rh[wqs_np[d]] <- list(do.call(rbind, lapply(all_results, function(res) res$mean_weight_rh[d,])))
    if (is.na(sum(unlist(all_rh_wqs[wqs_np[d]]), na.rm = TRUE))){
      Output[wqs_np[d]] <- NA
    } else {
      mean_weights <- colMeans(as.data.frame(mean_weight_rh[wqs_np[d]]), na.rm =TRUE)
      LCI_weights <- apply(as.data.frame(mean_weight_rh[wqs_np[d]]), 2, function(x) quantile(x, probs = 0.025, na.rm =TRUE))
      UCI_weights <- apply(as.data.frame(mean_weight_rh[wqs_np[d]]), 2, function(x) quantile(x, probs = 0.975, na.rm =TRUE))
      weight_output <- as.data.frame(bind_rows(mean_weights, LCI_weights, UCI_weights))
      row.names(weight_output) <- c("Mean", "LCI", "UCI")
      wqs_output <- c(mean(all_rh_wqs[[wqs_np[d]]], na.rm =TRUE), quantile(all_rh_wqs[[wqs_np[d]]], probs = 0.025, na.rm =TRUE), quantile(all_rh_wqs[[wqs_np[d]]], probs = 0.975, na.rm =TRUE))
      weight_output$Estimate <- wqs_output
      Output[wqs_np[d]] <- list(weight_output)
      }
  }
  
  return(list(Output = Output,
              all_weights = all_b_weights, all_rh_wqs = all_rh_wqs,
              mean_weight_rh = mean_weight_rh,
              parameter.setting = parameter.setting,
              start.year = min(model.data$year), end.year = max(model.data$year)))
}