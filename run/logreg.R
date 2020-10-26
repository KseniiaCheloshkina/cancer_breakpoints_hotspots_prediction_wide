get_log_reg_metrics <- function(x_train, y_train, x_test, y_test){
  # fit log reg
  x_logreg <- x_train
  cols <- names(x_train)
  x_logreg$y <- ifelse(y_train == "X0", 0, 1)
  mod <- glm(
    as.formula(paste("y ~ ", paste(cols, collapse=" + "))), 
    data=x_logreg, 
    family="binomial")
  nullmod <- glm("y ~ 1", data=x_logreg, family="binomial")
  cont_pred <- exp(as.numeric(nullmod$coefficients['(Intercept)'])[1])
  # get prediction
  test_pred <- data.frame(X1=predict(mod, newdata = x_test, type = "response"))
  test_pred$target <- as.character(y_test)
  test_pred$target_num <- as.numeric(gsub(x = test_pred$target, pattern = "X", replacement = ""))
  # mcfadden r2
  logLik_mod <- sum(log(test_pred$X1*test_pred$target_num + (1-test_pred$X1)*(1-test_pred$target_num)))
  logLik_null <- sum(log(cont_pred*test_pred$target_num + (1-cont_pred)*(1-test_pred$target_num)))
  mcfadden_r2 <- 1 - logLik_mod/logLik_null
  # efron r2
  mu <- mean(test_pred$target_num)
  efron_r2 <- 1 - sum((test_pred$target_num - test_pred$X1)^2)/sum((test_pred$target_num - mu)^2)
  return(c(mcfadden_r2, efron_r2))
}

