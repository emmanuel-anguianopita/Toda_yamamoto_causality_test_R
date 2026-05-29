#################################################################################################################################
# TODA-YAMAMOTO CAUSALITY TEST
#
# Implementacion de la prueba de causalidad de Granger de Toda y Yamamoto (1995).
#
# Referencia:
#   Toda, H. Y., & Yamamoto, T. (1995). Statistical inference in vector
#   autoregressions with possibly integrated processes.
#   Journal of Econometrics, 66(1-2), 225-250.
#
# Argumentos:
#   var.model : objeto varest estimado con vars::VAR(Y, p, type)
#   d.max     : orden maximo de integracion (default 1). Se define de forma
#               exogena con base en pruebas de raiz unitaria previas.
#               El LA-VAR se estima con var.model$p + d.max rezagos totales.
#   verbose   : si TRUE imprime las restricciones exactas de cada prueba
#               (default FALSE)
#
# Output (invisible):
#   ty.results       : data.frame con cause, effect, chisq, pvalue
#   ty.augmented_var : objeto varest del LA-VAR(p + d.max) estimado
#   ty.wald          : lista con objetos wald.test por cada par
#   ty.bg            : objeto serial.test (prueba BG con 12 rezagos)
#   ty.regressors    : vector de nombres de regresores del LA-VAR
#   ty.r2            : vector con R2 ajustado por ecuacion
#
# Ultima actualizacion: 29 de mayo 2026
#################################################################################################################################
toda.yamamoto <- function(var.model, d.max = 1L, verbose = FALSE) {
  
  require(vars)
  require(aod)
  
  ty.df <- data.frame(var.model$y)
  ty.varnames <- colnames(ty.df)
  ty.regressors <- NULL
  
  # k + d_max de acuerdo con Toda y Yamamoto (1995)
  ty.lags <- var.model$p + d.max
  ty.augmented_var <- vars::VAR(ty.df, ty.lags, type = var.model$type)

  ty.regressors <- setdiff(colnames(ty.augmented_var$datamat),
                           colnames(ty.augmented_var$y))
  
  ty.results <- data.frame(cause = character(0), 
                           effect = character(0), 
                           chisq = numeric(0), 
                           pvalue = numeric(0))
  
  ty.wald <- list()
  
  cat(sprintf("\nPrueba de Causalidad de Granger -- Toda & Yamamoto (1995)\n"))
  cat(sprintf("%s\n", strrep("-", 57)))
  cat(sprintf("LA-VAR(%d+%d=%d) | k=%d\n", 
              var.model$p, d.max, ty.lags, length(ty.varnames)))
  cat(sprintf("H0: la variable en 'cause' NO causa en el sentido de Granger a 'effect'\n"))
  cat(sprintf("Significancia: *** p<0.01  ** p<0.05  * p<0.10\n\n"))
  
  for (k in 1:length(ty.varnames)) {
    for (j in 1:length(ty.varnames)) {
      if (k != j) {
        
        ty.coefres <- grep(ty.varnames[j], ty.regressors)[seq_len(var.model$p)]
        
        if (verbose) {
          cat(sprintf("  Restricciones (%s => %s):\n",
                      ty.varnames[j], ty.varnames[k]))
          cat(sprintf("    %s\n", ty.regressors[ty.coefres]))
        }
        
        wald.res <- wald.test(b     = coef(ty.augmented_var$varresult[[k]]), 
                              Sigma = vcov(ty.augmented_var$varresult[[k]]),
                              Terms = ty.coefres)
        
        par_name <- paste0(ty.varnames[j], "_causa_", ty.varnames[k])
        ty.wald[[par_name]] <- wald.res
        
        chisq  <- as.numeric(wald.res$result$chi2[1])
        pvalue <- as.numeric(wald.res$result$chi2[3])
        sig    <- if      (pvalue < 0.01) "***"
        else if (pvalue < 0.05) "**"
        else if (pvalue < 0.10) "*"
        else                    ""
        
        cat(sprintf("%-25s => %-25s | Chi2(%d) = %7.4f | p = %.4f %s\n",
                    ty.varnames[j], ty.varnames[k],
                    var.model$p, chisq, pvalue, sig))
        
        ty.results <- rbind(ty.results, data.frame(
          cause  = ty.varnames[j], 
          effect = ty.varnames[k], 
          chisq  = chisq,
          pvalue = pvalue,
          stringsAsFactors = FALSE)
        )
      }
    }
  }
  
  # --- R2 ajustado por ecuacion ---
  ty.r2 <- round(sapply(ty.augmented_var$varresult,
                        function(eq) summary(eq)$adj.r.squared), 4)
  
  cat(sprintf("\n%s\n", strrep("-", 57)))
  cat(sprintf("R2 ajustado por ecuacion:\n\n"))
  for (i in seq_along(ty.r2)) {
    cat(sprintf("  %-25s R2 adj = %.4f\n", names(ty.r2)[i], ty.r2[i]))
  }
  
  # --- Test BG multivariado sobre el LA-VAR ---
  bg.res  <- vars::serial.test(ty.augmented_var, lags.pt = 12, type = "BG")
  bg.stat <- round(as.numeric(bg.res$serial$statistic), 4)
  bg.gl   <- as.integer(bg.res$serial$parameter)
  bg.pval <- round(as.numeric(bg.res$serial$p.value), 4)
  bg.sig  <- if      (bg.pval < 0.01) "***"
  else if (bg.pval < 0.05) "**"
  else if (bg.pval < 0.10) "*"
  else                     ""
  
  cat(sprintf("\n%s\n", strrep("-", 57)))
  cat(sprintf("Test BG multivariado (lags = 12)\n"))
  cat(sprintf("H0: no hay autocorrelacion serial hasta el rezago 12\n"))
  cat(sprintf("BG = %.4f | gl = %d | p-valor = %.4f %s\n",
              bg.stat, bg.gl, bg.pval, bg.sig))
  cat(sprintf("Significancia: *** p<0.01  ** p<0.05  * p<0.10\n"))
  
  return(invisible(list(
    ty.results       = ty.results,
    ty.augmented_var = ty.augmented_var,
    ty.wald          = ty.wald,
    ty.bg            = bg.res,
    ty.regressors    = ty.regressors,
    ty.r2            = ty.r2
  )))
}