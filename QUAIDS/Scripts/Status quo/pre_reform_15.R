## ================== QUAIDS Delta CIs (plug-and-play) ==================
## Requisitos: numDeriv
if (!requireNamespace("numDeriv", quietly = TRUE))
  stop("Instale: install.packages('numDeriv')")

## --------- Util: detecção de nomes no seu dataset ---------
.detect_shares <- function(dat) {
  s <- grep("^w_", names(dat), value = TRUE)
  if (!length(s)) stop("Não encontrei colunas de shares que comecem com 'w_'.")
  # Ordena por dígito final, se houver (ex.: ...1,...2,...)
  ord <- order(suppressWarnings(as.integer(gsub(".*?(\\d+)$", "\\1", s))))
  s[ord]
}
.detect_lnexp <- function(dat) {
  cand <- intersect(c("lndespesa","ln_despesa","lnexp","ln_expenditure"), names(dat))
  if (length(cand)) return(cand[1])
  # fallback por regex
  pick <- grep("^ln.*(desp|exp)", names(dat), value = TRUE)
  if (!length(pick)) stop("Não encontrei a variável ln(x) (ex.: 'lndespesa').")
  pick[1]
}
.detect_sq <- function(dat) {
  cand <- intersect(c("z2","lndespesa2","lnexp2"), names(dat))
  if (length(cand)) return(cand[1])
  # Se não existir explícita, vamos criar on-the-fly quando preciso
  NULL
}
.detect_lnprices <- function(dat) {
  p <- grep("^ln_.*preco", names(dat), value = TRUE)
  if (!length(p)) stop("Não encontrei colunas de preços em log (ex.: 'ln_preco_por_kg#').")
  # Ordena por dígito final
  ord <- order(suppressWarnings(as.integer(gsub(".*?(\\d+)$", "\\1", p))))
  p[ord]
}

## --------- Empacotamento/Desempacotamento θ ----------
.build_theta_from_systemfit <- function(fit_sys, dat) {
  cf  <- stats::coef(fit_sys)
  Vcf <- stats::vcov(fit_sys)
  if (is.null(names(cf))) stop("Coeficientes sem nomes no systemfit.")
  
  eq_names <- names(fit_sys$eq)              # ex.: w_despesahat2,...,w_despesahat6
  shares   <- .detect_shares(dat)            # ex.: w_despesahat1:6
  K        <- length(shares)
  # qual share foi omitida na estimação?
  share_omit <- setdiff(shares, eq_names)
  if (length(share_omit) != 1) {
    stop("Não consegui identificar exatamente 1 equação omitida (shares no data vs. eqs no systemfit).")
  }
  share_omit <- share_omit[1]
  
  xvar <- .detect_lnexp(dat)
  z2var <- .detect_sq(dat) # pode ser NULL
  pvars <- .detect_lnprices(dat) # K preços em log no data
  # preços presentes no modelo (geralmente K-1, um drop)
  present_price_vars <- unique(gsub("^.*?_", "", grep(paste0("^(", paste(eq_names, collapse="|"), ")_ln"), names(cf), value=TRUE)))
  present_price_vars <- intersect(present_price_vars, pvars)
  
  # Constrói θ em ordem fixa:
  # θ = [ alpha(eq!=omit) , beta(eq!=omit) , lambda(eq!=omit) , vec(Gamma para eq!=omit e pvars presentes) ]
  theta <- numeric(0); tnames <- character(0); pick_idx <- integer(0)
  
  # α (intercepto)
  for (eq in eq_names) {
    nm <- paste0(eq, "_(Intercept)")
    val <- if (nm %in% names(cf)) cf[[nm]] else 0
    theta <- c(theta, val)
    tnames <- c(tnames, paste0("alpha|", eq))
    pick_idx <- c(pick_idx, match(nm, names(cf)))
  }
  # β (coef de ln despesa)
  for (eq in eq_names) {
    nm <- paste0(eq, "_", xvar)
    val <- if (nm %in% names(cf)) cf[[nm]] else 0
    theta <- c(theta, val)
    tnames <- c(tnames, paste0("beta|", eq))
    pick_idx <- c(pick_idx, match(nm, names(cf)))
  }
  # λ (coef de (ln desp)^2) — pode não existir; então vira 0
  for (eq in eq_names) {
    nm <- if (!is.null(z2var)) paste0(eq, "_", z2var) else NA_character_
    val <- if (!is.na(nm) && nm %in% names(cf)) cf[[nm]] else 0
    theta <- c(theta, val)
    tnames <- c(tnames, paste0("lambda|", eq))
    pick_idx <- c(pick_idx, if (!is.na(nm)) match(nm, names(cf)) else NA_integer_)
  }
  # Γ (coef. dos ln preços)
  for (eq in eq_names) {
    for (pj in present_price_vars) {
      nm <- paste0(eq, "_", pj)
      val <- if (nm %in% names(cf)) cf[[nm]] else 0
      theta <- c(theta, val)
      tnames <- c(tnames, paste0("gamma|", eq, "|", pj))
      pick_idx <- c(pick_idx, match(nm, names(cf)))
    }
  }
  
  # Matriz V_θ (preenche com zeros onde o parâmetro está como 0 não-estimado)
  p <- length(theta)
  Vth <- matrix(0, p, p, dimnames = list(tnames, tnames))
  have <- which(!is.na(pick_idx))
  if (length(have)) {
    idx <- pick_idx[have]
    Vth[have, have] <- Vcf[idx, idx, drop = FALSE]
  }
  
  list(theta = theta, V = Vth, tnames = tnames,
       eq_names = eq_names, shares = shares, share_omit = share_omit,
       xvar = xvar, z2var = z2var, pvars = pvars,
       present_price_vars = present_price_vars)
}

## Reconstrói α, β, λ, Γ COMpletos (K bens x K preços) a partir de θ
.rebuild_params_quaids <- function(theta, meta) {
  eq_names <- meta$eq_names; shares <- meta$shares
  present_price_vars <- meta$present_price_vars
  pvars_all <- meta$pvars
  K <- length(shares)
  
  # mapeia eq -> índice do bem pelo dígito final
  idx_good <- function(nm) {
    as.integer(sub(".*?(\\d+)$", "\\1", nm))
  }
  eq_idx <- idx_good(eq_names)
  omit_idx <- idx_good(meta$share_omit) # bem omitido na estimação
  
  # splits de θ
  L <- length(eq_names)
  a_vec <- theta[1:L]
  b_vec <- theta[(L+1):(2*L)]
  l_vec <- theta[(2*L+1):(3*L)]
  rest  <- theta[-(1:(3*L))]
  
  # Γ observada: L linhas (eqs) x (K-1 ou menos) colunas presentes
  G_obs <- matrix(rest, nrow = L, byrow = TRUE,
                  dimnames = list(eq_names, present_price_vars))
  
  # monta α, β, λ de tamanho K (preenche o omitido por adding-up)
  alpha <- beta <- lambda <- rep(NA_real_, K)
  alpha[eq_idx]  <- a_vec
  beta[eq_idx]   <- b_vec
  lambda[eq_idx] <- l_vec
  # adding-up
  alpha[omit_idx]  <- 1 - sum(alpha[-omit_idx], na.rm = TRUE)
  beta[omit_idx]   <- - sum(beta[-omit_idx],  na.rm = TRUE)
  lambda[omit_idx] <- - sum(lambda[-omit_idx],na.rm = TRUE)
  
  # Γ completa KxK:
  # 1) coloca blocos observados (linhas das eqs estimadas; colunas presentes)
  Gamma <- matrix(0, K, K, dimnames = list(shares, pvars_all))
  Gamma[eq_idx, present_price_vars] <- G_obs
  
  # 2) Preenche coluna do preço drop (as que NÃO estão em present_price_vars) via homogeneidade por linha
  missing_cols <- setdiff(pvars_all, present_price_vars)  # normalmente 1 coluna
  if (length(missing_cols)) {
    for (i in eq_idx) {
      Gamma[i, missing_cols] <- -sum(Gamma[i, present_price_vars, drop=TRUE])
    }
  }
  # 3) Preenche a linha da share omitida via adding-up por coluna
  for (pj in pvars_all) {
    Gamma[omit_idx, pj] <- -sum(Gamma[setdiff(seq_len(K), omit_idx), pj])
  }
  # 4) Suaviza simetria (opcional, ajuda em ln a(p))
  Gamma <- 0.5*(Gamma + t(Gamma))
  
  list(alpha = alpha, beta = beta, lambda = lambda, Gamma = Gamma)
}

## --------- Elasticidades QUAIDS (Banks–Blundell–Lewbel) ----------
.quaids_elasticities_from_params <- function(alpha, beta, lambda, Gamma,
                                             lnx, pvec, wvec,
                                             shares, pvars) {
  K <- length(shares)
  if (length(pvec) != K) stop("Comprimento de pvec diferente de K.")
  if (length(wvec) != K) stop("Comprimento de wvec diferente de K.")
  
  lnp <- log(pvec)
  # ln a(p) = sum_i alpha_i ln p_i + 1/2 sum_i sum_j Gamma_ij ln p_i ln p_j
  lna <- sum(alpha * lnp) + 0.5 * as.numeric(t(lnp) %*% Gamma %*% lnp)
  # b(p) = prod p_i^{beta_i}  => ln b = sum beta_i ln p_i
  lnb <- sum(beta * lnp); b <- exp(lnb)
  L <- lnx - lna
  
  # μ_i = ∂w_i/∂ ln m
  mu <- beta + (2 * lambda / b) * L
  
  # ∂ ln a / ∂ ln p_j  =  alpha_j + sum_k Gamma_jk ln p_k   (simetria usada)
  dlnA_dlnp <- alpha + as.numeric(Gamma %*% lnp)
  
  # μ_ij = ∂ w_i / ∂ ln p_j
  # = Gamma_ij − μ_i * dlnA_dlnp_j − (lambda_i * beta_j / b) * L^2
  L2 <- L^2
  mu_ij <- Gamma - outer(mu, dlnA_dlnp) - outer(lambda / b * L2, beta)
  
  # Elasticidades:
  # Gasto: e_i = mu_i / w_i + 1
  e_x <- mu / pmax(wvec, 1e-12) + 1
  # Marshall: e^u_ij = mu_ij / w_i − δ_ij
  delta <- diag(K)
  e_u <- sweep(mu_ij, 1, pmax(wvec, 1e-12), "/") - delta
  # Hicks: e^c_ij = e^u_ij + e_i w_j
  e_c <- e_u + outer(e_x, wvec)
  
  list(marshall = e_u, hicks = e_c, expenditure = e_x,
       lna = lna, b = b, L = L)
}

## --------- Função principal: Delta CI para QUAIDS ----------
delta_quaids_CIs <- function(fit_obj, dat, level = 0.95, method = c("pointwise","bonferroni")) {
  method <- match.arg(method)
  # aceita systemfit direto ou wrapper com $fit
  fit_sys <- if (inherits(fit_obj, "systemfit")) fit_obj else fit_obj$fit
  if (!inherits(fit_sys, "systemfit")) stop("Passe o objeto systemfit (ou um wrapper que tenha $fit).")
  
  meta  <- .build_theta_from_systemfit(fit_sys, dat)
  theta0<- meta$theta
  Vth   <- meta$V
  
  shares <- meta$shares
  pvars  <- meta$pvars
  K      <- length(shares)
  
  # ponto de avaliação:
  # w* = mediana observada de cada share
  w_star <- sapply(shares, function(s) stats::median(dat[[s]], na.rm = TRUE))
  # p* = exp(mean(ln p_j))
  p_star <- sapply(pvars, function(v) exp(mean(dat[[v]], na.rm = TRUE)))
  # ln x*:
  xvar <- meta$xvar
  lnx_star <- stats::median(dat[[xvar]], na.rm = TRUE)
  
  # g(θ): vetor = [vec(Hicks), vec(Marshall), eta]
  g_eval <- function(th) {
    par <- .rebuild_params_quaids(th, meta)
    E <- .quaids_elasticities_from_params(
      alpha  = par$alpha,
      beta   = par$beta,
      lambda = par$lambda,
      Gamma  = par$Gamma,
      lnx    = lnx_star,
      pvec   = p_star,
      wvec   = w_star,
      shares = shares,
      pvars  = pvars
    )
    # vetorizar por LINHA (byrow=TRUE)
    c(as.vector(t(E$hicks)), as.vector(t(E$marshall)), as.numeric(E$expenditure))
  }
  
  g0 <- g_eval(theta0)
  J  <- numDeriv::jacobian(g_eval, theta0)
  Vg <- J %*% Vth %*% t(J)
  se <- sqrt(pmax(diag(Vg), 0))
  
  m <- length(g0); alpha <- 1 - level
  z <- if (method=="bonferroni") qnorm(1 - alpha/(2*m)) else qnorm(1 - alpha/2)
  lo <- g0 - z*se; hi <- g0 + z*se
  
  # Reconstrói blocos
  par0 <- .rebuild_params_quaids(theta0, meta)
  E0   <- .quaids_elasticities_from_params(par0$alpha, par0$beta, par0$lambda, par0$Gamma,
                                           lnx_star, p_star, w_star, shares, pvars)
  rn <- shares; cn <- shares
  k2 <- K*K
  idx_H <- seq_len(k2)
  idx_M <- k2 + seq_len(k2)
  idx_E <- (2*k2) + seq_len(K)
  
  H_hat <- matrix(g0[idx_H], nrow=K, byrow=TRUE, dimnames=list(rn, cn))
  H_lo  <- matrix(lo[idx_H], nrow=K, byrow=TRUE, dimnames=list(rn, cn))
  H_hi  <- matrix(hi[idx_H], nrow=K, byrow=TRUE, dimnames=list(rn, cn))
  H_se  <- matrix(se[idx_H], nrow=K, byrow=TRUE, dimnames=list(rn, cn))
  
  M_hat <- matrix(g0[idx_M], nrow=K, byrow=TRUE, dimnames=list(rn, cn))
  M_lo  <- matrix(lo[idx_M], nrow=K, byrow=TRUE, dimnames=list(rn, cn))
  M_hi  <- matrix(hi[idx_M], nrow=K, byrow=TRUE, dimnames=list(rn, cn))
  M_se  <- matrix(se[idx_M], nrow=K, byrow=TRUE, dimnames=list(rn, cn))
  
  eta_hat <- setNames(g0[idx_E], rn)
  eta_lo  <- setNames(lo[idx_E], rn)
  eta_hi  <- setNames(hi[idx_E], rn)
  eta_se  <- setNames(se[idx_E], rn)
  
  # Tabelas tidy
  mk_tbl <- function(M, lo, hi, se) {
    out <- do.call(rbind, lapply(seq_len(K), function(i) {
      data.frame(i = rn[i], j = cn, est = M[i,], se = se[i,], lwr = lo[i,], upr = hi[i,], row.names = NULL)
    }))
    rownames(out) <- NULL; out
  }
  
  list(
    marshallian = mk_tbl(M_hat, M_lo, M_hi, M_se),
    hicksian    = mk_tbl(H_hat, H_lo, H_hi, H_se),
    expenditure = data.frame(good = rn, est = eta_hat, se = eta_se, lwr = eta_lo, upr = eta_hi, row.names = NULL),
    level       = level,
    method      = paste0("Delta (QUAIDS) — w* observada (mediana); x*=median(", meta$xvar, "); p*=exp(mean ln p)")
  )
}

## ----------------- Exemplo de uso -----------------
out_delta <- delta_quaids_CIs(fit_just$fit, df_iv, level = 0.95)
cat("\n== ELASTICIDADES — MARSHALL (não compensadas) — IC Delta ==\n")
print(out_delta$marshallian, row.names = FALSE)
cat("\n== ELASTICIDADES — HICKS (compensadas) — IC Delta ==\n")
print(out_delta$hicksian, row.names = FALSE)
cat("\n== ELASTICIDADES — RENDA/DESPESA — IC Delta ==\n")
print(out_delta$expenditure, row.names = FALSE)
