# Prepare data clustering original list

# functions for scaling and removing outliers
preparedata <- function(input) {
               modeller_nonscaled <- str_replace_all(modeller,"_scld","")
                data.scld <- data.anal %>% 
                            filter_at(vars(one_of(modeller_nonscaled)),all_vars(!is.na(.))) %>%  
                            filter_at(vars(one_of(modeller_nonscaled)),all_vars(!is.infinite(.))) %>% 
                            group_by(sexe) %>% 
                            mutate_at(vars(one_of(modeller_nonscaled[!str_detect(modeller_nonscaled,"bin")])),
                                                    funs(scld = as.numeric(scale(.)))) %>%  # scale based on different sex
                            filter_at(vars(ends_with("scld")),all_vars(abs(.)<5)) %>% # remove the outlier
                            ungroup() 
               return(data.scld)
}






# Prepare data clustering proxy list

preparedata_w <- function(input) {
               modeller_nonscaled_w <- str_replace_all(modeller_w,"_scld","")
                data.scld_w <- data.anal_w %>% 
                            filter_at(vars(one_of(modeller_nonscaled_w)),all_vars(!is.na(.))) %>%  
                            filter_at(vars(one_of(modeller_nonscaled_w)),all_vars(!is.infinite(.))) %>% 
                            group_by(sexe) %>% 
                            mutate_at(vars(one_of(modeller_nonscaled_w[!str_detect(modeller_nonscaled_w,"bin")])),
                                                    funs(scld = as.numeric(scale(.)))) %>% # scale based on different sex, generate the new variables end with _scld
                            filter_at(vars(ends_with("scld")),all_vars(abs(.)<5)) %>% # remove the outlier
                            ungroup() 
               return(data.scld_w)
}



# Remove outliers metabolites
is_outlier <- function(object, nSD = nSD){
  object <- as.numeric(object)
  v_mean <- mean(object, na.rm = T)
  v_median <- median(object, na.rm = T)
  v_sd <- sd(object, na.rm = T)
  v_IQR <- IQR(object, na.rm = T)
  # r <- which(object > (v_median + nSD * v_IQR) | object < (v_median - nSD * v_IQR))
  r <- which(object > (v_mean + nSD * v_sd) | object < (v_mean - nSD * v_sd))
  return(r)
}


impute_min2 <- function(data) {
  imputed_data <- data
  for (col in names(imputed_data)) {
    min_val <- min(imputed_data[[col]], na.rm = TRUE) / 2
    imputed_data[[col]] <- ifelse(is.na(imputed_data[[col]]), min_val, imputed_data[[col]])
  }
  return(imputed_data)
}



log_reg <- function(data, response_var, predictor_vars) {
  
  formula <- as.formula(paste(response_var, "~", paste(predictor_vars, collapse = "+")))
  
  model <- glm(formula, data = data, family = binomial)
  
  return(model)
}



violin_plot <- function(y_var, y_label, ylim) {
  custom_colors <- c("#88c584","#7dadd4","#E7B800","#8f6799")
  pos <- position_dodge(0.9)
  
  
  ggplot(df_surv, aes(x = c_name_plot, y= y_var, fill = c_name_plot)) +
    geom_violin(alpha = 0.25, position = pos) +
    geom_boxplot(width = .2, 
                 #fatten = NULL, 
                 alpha = 0.75,
                 position = pos,
                 outlier.shape = NA) +
    scale_fill_manual(values = custom_colors) + 
    theme_classic() + 
    theme(panel.grid.major.y = element_line(),
          panel.grid.minor.y = element_line(), 
          legend.position = "none",
          panel.border = element_rect(color = "black", fill = NA, size = 0.5)) + 
    scale_y_continuous(limits = ylim) + 
    labs(x = "Cluster", y = y_label)
}




