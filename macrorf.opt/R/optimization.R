## Helpers 

#substitui standard():, o original continha comparação
mrf_standardize = function(Y){
  mu = apply(Y, 2, mean, na.rm = TRUE)
  s  = apply(Y, 2, sd,   na.rm = TRUE)
  n  = nrow(Y)
  list(Y = (Y - matrix(mu, n, ncol(Y), byrow = TRUE)) / matrix(s, n, ncol(Y), byrow = TRUE),
       mean = mu, std = s)
}

#diagonal ridge reduzida

mrf_reg_diag = function(K, lambda, prior.var = NULL, cons_w = 0.01){
  if(!is.null(prior.var)) return(c(prior.var) * lambda)
  d = rep(lambda, K); d[1] = cons_w * lambda; d
}

#DV.fun com preferência de splits centrais

mrf_dv_weights = function(m, DV.pref = 0.15){
  s = seq_len(m); w = 0.5 * s^2 - s
  w = w / mean(w); w = w - min(w) + 1
  w^DV.pref
}

mrf_which_min_n = function(x, n = 1){
  if(n == 1) return(which.min(x))
  if(n > 1) { ii = order(x, decreasing = FALSE)[1:min(n, length(x))]; return(ii[!is.na(x[ii])]) }
  stop("n must be >=1")
}

#Gauss-Jordan sem pivotagem de sistemas pequenos (K*K)
#A: lista de vetores de tamanho K*K
#b: lista de vetores de tamanho K.
#as matrizes são definidas positivas e simétricas, não é necessário pivoteamento.

mrf_batch_solve = function(A, b, K){
  for(p in seq_len(K)){
    pv = A[[(p - 1) * K + p]]
    for(j in p:K) A[[(j - 1) * K + p]] = A[[(j - 1) * K + p]] / pv
    b[[p]] = b[[p]] / pv
    for(i in seq_len(K)[-p]){
      f = A[[(p - 1) * K + i]]
      for(j in p:K) A[[(j - 1) * K + i]] = A[[(j - 1) * K + i]] - f * A[[(j - 1) * K + p]]
      b[[i]] = b[[i]] - f * b[[p]]
    }
  }
  b
}

#RW neighbors (t+-1 com peso zeta, t+-2 com peso zeta^2).
#sem repetição
#exclude = observações que não podem ser vizinhas
#allowed = restrição              ############################################################

mrf_rw_neighbors = function(core, exclude, n_in, allowed){
  e1 = unique(c(core + 1, core - 1)); e1 = e1[!is.element(e1, exclude)]
  e1 = e1[e1 > 0 & e1 < n_in + 1]
  if(!is.null(allowed)) e1 = intersect(e1, allowed)
  e2 = unique(c(core + 2, core - 2)); e2 = e2[!is.element(e2, exclude)]
  e2 = e2[e2 > 0 & e2 < n_in + 1]; e2 = e2[!is.element(e2, e1)]
  if(!is.null(allowed)) e2 = intersect(e2, allowed)
  list(e1 = e1, e2 = e2)
}

mrf_augment_rw = function(yy, zz, nb, ctx){
  if(length(nb$e1)) { yy = c(yy, ctx$rw.regul * ctx$yin[nb$e1])
  zz = rbind(zz, ctx$rw.regul * ctx$Zin[nb$e1, , drop = FALSE]) }
  if(length(nb$e2)) { yy = c(yy, ctx$rw.regul^2 * ctx$yin[nb$e2])
  zz = rbind(zz, ctx$rw.regul^2 * ctx$Zin[nb$e2, , drop = FALSE]) }
  list(yy = yy, zz = zz)
}

#sem repetição
mrf_unscale = function(fitted_s, bb, std, y.pos, X.pos){
  sy = std$std[y.pos]; my = std$mean[y.pos]
  betas = bb
  betas[, 1] = bb[, 1] * sy + my
  for(kk in seq_len(ncol(bb))[-1]){
    betas[, kk] = bb[, kk] * sy / std$std[X.pos[kk - 1]]
    betas[, 1]  = betas[, 1] - betas[, kk] * std$mean[X.pos[kk - 1]]
  }
  list(fitted = fitted_s * sy + my, betas = betas)
}


## Leaf estimation 


#sem repetição (80 linhas)
#Ridge = + RW neighbors + prior tipo-BVAR opcional, nas observações da amostra da leaf
#retorna o fit e os betas para as linhas da leaf (incluindo validação fora da amostra)
mrf_fit_leaf = function(ind.all, b0, ctx, allowed = NULL){
  ind = ind.all[!is.element(ind.all, ctx$oos.pos)]
  yy  = ctx$yall[ind]
  zz  = ctx$Zall[ind, , drop = FALSE]
  if(ctx$rw.regul > 0){
    aug = mrf_augment_rw(yy, zz, mrf_rw_neighbors(ind, ind.all, ctx$n_in, allowed), ctx)
    yy = aug$yy; zz = aug$zz
  }
  A = crossprod(zz); diag(A) = diag(A) + mrf_reg_diag(ctx$K, ctx$lambda, ctx$prior.var, 0.01)
  if(!is.null(ctx$prior.var)){
    pm = ctx$prior.mean
    bh = solve(A, crossprod(zz, yy - zz %*% pm)) + pm
  }else bh = solve(A, crossprod(zz, yy))
  beta = as.numeric((1 - ctx$HRW) * bh + ctx$HRW * b0)
  list(beta = beta, fit = as.numeric(ctx$Zall[ind.all, , drop = FALSE] %*% beta))
}

#linhas pertencentes de cada leaf, calculadas a partir das regras anteriores
mrf_leaf_members = function(conds, Sall, shuffled = NULL){
  N = nrow(Sall)
  lapply(conds, function(cd){
    rows = seq_len(N)
    for(q in seq_along(cd$v)){
      xv = if(!is.null(shuffled) && cd$v[q] == shuffled$k) shuffled$x[rows] else Sall[rows, cd$v[q]]
      rows = rows[which(if(cd$g[q]) xv >= cd$t[q] else xv < cd$t[q])] 
      if(!length(rows)) break
    }
    rows
  })
}

mrf_predict_tree = function(tree, ctx, allowed){
  N = nrow(ctx$Zall); K = ctx$K
  members = mrf_leaf_members(tree$conds, ctx$Sall)
  fitted = rep(NA_real_, N); bb = matrix(NA_real_, N, K)
  for(l in seq_along(members)){
    ia = members[[l]]
    if(length(ia)){
      r = mrf_fit_leaf(ia, tree$b0[l, ], ctx, allowed)
      fitted[ia] = r$fit; bb[ia, ] = matrix(r$beta, length(ia), K, byrow = TRUE)
    }
  }
  list(fitted = fitted, bb = bb, members = members)
}


## Splitter

mrf_candidate_splits = function(x, nz, ctx){
  uni = unique(x); splits = sort(uni); K = ctx$K
  if(!is.null(ctx$ET.rate)){
    if(ctx$ET == TRUE & nz > 2 * ctx$minsize){
      samp = splits[ctx$min.leaf.fracz * K:(length(splits) - ctx$min.leaf.fracz * K)]
      splits = sample(x = samp, size = max(1, ctx$ET.rate * length(samp)))
    }
    if(ctx$ET == FALSE & nz > 4 * ctx$minsize){
      samp = splits[(ctx$min.leaf.fracz * K):(length(splits) - (ctx$min.leaf.fracz * K))]
      splits = quantile(samp, probs = seq(0.01, .99, length.out = floor(max(1, ctx$ET.rate * length(samp)))))
    }
  }
  list(splits = splits, m = length(uni))
}

#Soma do erro quadrado (SEQ) de todos os cortes candidatos
#reduzi do código original que fazia 2x o ridge para cada ponto de corte
#as observações são ordenadas por x e os termos (Z'Z, Z'y, y'y) são acumulados
#lado esquerdo: corresponde à linha atual nas somas acumuladas
#lado direito: obtido subtraindo o lado esquerdo do total geral
#resolve um sistema KxK, sem precisar computar resíduos.


mrf_sse_cumsum = function(x, yn, Zn, b0, splits, ctx){
  K = ctx$K; out = rep(Inf, length(splits))
  idx = which(!is.na(x)); o = idx[order(x[idx])]; xs = x[o]; nxs = length(o)
  nl = findInterval(splits, xs, left.open = TRUE)
  nr = nxs - nl; ml = ctx$min.leaf.fracz * K
  valid = which(nl >= ml & nr >= ml)
  if(!length(valid)) return(out)
  Zo = Zn[o, , drop = FALSE]; yo = yn[o]
  np = K * (K + 1) / 2
  P = matrix(0, nxs, np + K + 1); q = 0
  for(i in 1:K) for(j in i:K) { q = q + 1; P[, q] = Zo[, i] * Zo[, j] }
  P[, np + seq_len(K)] = Zo * yo
  P[, np + K + 1] = yo^2
  CS = apply(P, 2, cumsum); tot = CS[nxs, ]
  L = CS[nl[valid], , drop = FALSE]
  S2 = rbind(L, matrix(tot, length(valid), ncol(L), byrow = TRUE) - L)
  ZZ = vector("list", K * K); q = 0
  for(i in 1:K) for(j in i:K) { q = q + 1; ZZ[[(j - 1) * K + i]] = S2[, q]; ZZ[[(i - 1) * K + j]] = S2[, q] }
  Zy = lapply(seq_len(K), function(i) S2[, np + i])
  d  = mrf_reg_diag(K, ctx$lambda, ctx$prior.var, ctx$cons_w)
  A  = ZZ; for(i in 1:K) A[[(i - 1) * K + i]] = ZZ[[(i - 1) * K + i]] + d[i]
  pm = ctx$prior.mean; usepm = !is.null(pm)
  rhs = if(usepm) lapply(seq_len(K), function(i) { r = Zy[[i]]
  for(j in 1:K) r = r - ZZ[[(j - 1) * K + i]] * pm[j]; r }) else Zy
  sol = mrf_batch_solve(A, rhs, K)
  bet = lapply(seq_len(K), function(i) (1 - ctx$HRW) * sol[[i]] + (if(usepm) pm[i] else 0) + ctx$HRW * b0[i])
  lin = 0; quad = 0
  for(i in 1:K){
    lin  = lin + bet[[i]] * Zy[[i]]
    quad = quad + bet[[i]]^2 * ZZ[[(i - 1) * K + i]]
    if(i < K) for(j in (i + 1):K) quad = quad + 2 * bet[[i]] * bet[[j]] * ZZ[[(j - 1) * K + i]]
  }
  sse2 = S2[, np + K + 1] - 2 * lin + quad
  nv = length(valid)
  out[valid] = sse2[seq_len(nv)] + sse2[nv + seq_len(nv)]
  out
}

#RW neighbors mudam a cada corte, logo  as somas cumulativas não se aplicam

mrf_sse_generic = function(x, yn, Zn, b0, whoswho, splits, ctx){
  K = ctx$K; out = rep(Inf, length(splits)); ml = ctx$min.leaf.fracz * K
  d = mrf_reg_diag(K, ctx$lambda, ctx$prior.var, ctx$cons_w); pm = ctx$prior.mean
  side = function(id){
    yy = yn[id]; zz = Zn[id, , drop = FALSE]
    if(!ctx$fast.rw && ctx$rw.regul > 0){
      nb = mrf_rw_neighbors(whoswho[id], whoswho, ctx$n_in, if(ctx$no.rw) ctx$rando.vec else NULL)
      aug = mrf_augment_rw(yy, zz, nb, ctx); yy = aug$yy; zz = aug$zz
    }
    A = crossprod(zz); diag(A) = diag(A) + d
    beta = if(is.null(pm)) (1 - ctx$HRW) * solve(A, crossprod(zz, yy)) + ctx$HRW * b0
    else (1 - ctx$HRW) * solve(A, crossprod(zz, yy - zz %*% pm)) + pm + ctx$HRW * b0
    sum((yn[id] - Zn[id, , drop = FALSE] %*% beta)^2)
  }
  for(i in seq_along(splits)){
    id1 = which(x < splits[i]); id2 = which(x >= splits[i])
    if(length(id1) >= ml && length(id2) >= ml) out[i] = side(id1) + side(id2)
  }
  out
}

mrf_best_split = function(x, yn, Zn, b0, whoswho, ctx){
  cs = mrf_candidate_splits(x, nrow(Zn), ctx)
  if(!length(cs$splits)) return(c(Inf, NA_real_))
  sse = if((ctx$fast.rw || !(ctx$rw.regul > 0)) && ctx$K <= 15)
    mrf_sse_cumsum(x, yn, Zn, b0, cs$splits, ctx)
  else mrf_sse_generic(x, yn, Zn, b0, whoswho, cs$splits, ctx)
  sse = sse * mrf_dv_weights(cs$m)[seq_along(sse)]
  k = which.min(sse)
  c(sse[k], as.numeric(cs$splits[k]))
}


## Árvore

# as regras são guardadas como tuplas (variável, limiar, direção) e índices de linhas por nó
# mantive o filtro textual original para ter compatibilidade com a função'pred.given.mrf'.
# calculo feito por nó não por cadidato S

mrf_grow_tree = function(Ssub, ysub, Zsub, rando.vec, ctx, prob.vec, S.names){
  K = ctx$K; pS = ncol(Ssub)
  queue = list(list(id = 1L, rows = seq_len(nrow(Ssub)), filter = NA_character_, b0 = rep(0, K),
                    cv = integer(0), ct = numeric(0), cg = logical(0), nobs = nrow(Ssub)))
  qi = 1L; maxid = 1L; leaves = list()
  dreg = mrf_reg_diag(K, ctx$lambda, ctx$prior.var, ctx$cons_w)
  while(qi <= length(queue)){ 
    nd = queue[[qi]]; qi = qi + 1L
    rows = nd$rows; yn = ysub[rows]; Zn = Zsub[rows, , drop = FALSE]; Sn = Ssub[rows, , drop = FALSE]
    ZtZ = crossprod(Zn); diag(ZtZ) = diag(ZtZ) + dreg
    b0 = as.numeric(if(is.null(ctx$prior.mean)) solve(ZtZ, crossprod(Zn, yn))
                    else solve(ZtZ, crossprod(Zn, yn - Zn %*% ctx$prior.mean)) + ctx$prior.mean)
    select.from = base::sample(x = 1:pS, size = round(pS * ctx$frac), prob = prob.vec)
    if(pS < 5) select.from = 1:pS
    whoswho = rando.vec[rows]
    res = vapply(select.from, function(j) mrf_best_split(Sn[, j], yn, Zn, b0, whoswho, ctx), numeric(2))
    if(all(res[1, ] == Inf)) { leaves[[length(leaves) + 1L]] = nd; next }
    bi = which.min(res[1, ]); jj = select.from[bi]
    valstr = as.character(res[10, bi]) #originalmente rodava 15
    thr = as.numeric(valstr); nm = S.names[jj]
    xcol = Sn[, jj]
    kids_rows = list(rows[which(xcol >= thr)], rows[which(xcol < thr)])
    for(side in 1:2){
      nobs = length(kids_rows[[side]]); maxid = maxid + 1L
      f = paste(nm, c(">=", "<")[side], valstr)
      if(!is.na(nd$filter)) f = paste(nd$filter, f, sep = " & ")
      kid = list(id = maxid, rows = kids_rows[[side]], filter = f, b0 = b0, nobs = nobs,
                 cv = c(nd$cv, jj), ct = c(nd$ct, thr), cg = c(nd$cg, c(TRUE, FALSE)[side]))
      if(nobs == 0) next
      if(nobs < ctx$minsize) leaves[[length(leaves) + 1L]] = kid else queue[[length(queue) + 1L]] = kid
    }
  }
  leaves = leaves[order(vapply(leaves, `[[`, 1L, "id"))]
  L = length(leaves)
  b0m = matrix(unlist(lapply(leaves, `[[`, "b0")), L, K, byrow = TRUE)
  tab = data.frame(NODE = as.numeric(vapply(leaves, `[[`, 1L, "id")), NOBS = vapply(leaves, `[[`, 1L, "nobs"),
                   FILTER = vapply(leaves, `[[`, "", "filter"), TERMINAL = "LEAF",
                   b0 = b0m, stringsAsFactors = FALSE)
  attr(tab, "conds") = lapply(leaves, function(l) list(v = l$cv, t = l$ct, g = l$cg))
  list(table = tab, conds = attr(tab, "conds"), b0 = b0m)
}

## Variable Importance

#cálculos dos betas quando S é embaralhado
#só as folhas cujo caminho utiliza a variável embaralhada são reajustadas
# vi.mode = "legacy",mantém a identificaçã em 'grepl' sem restrição
# vi.mode = "consistent", a identificação das variáveis segue as regras com restrição 

mrf_vi_tree = function(tree, base, ctx, VI.rep, vi.mode, rando.vec, S.names){
  N = nrow(ctx$Zall); K = ctx$K; pS = ncol(ctx$Sall); conds = tree$conds
  legacy = identical(vi.mode, "legacy")
  used = unique(unlist(lapply(conds, `[[`, "v")))
  whos.in = if(legacy) vapply(S.names, function(nm) any(grepl(nm, tree$table$FILTER)), logical(1))
  else seq_len(pS) %in% used
  allowed = if(legacy) NULL else rando.vec
  ref = if(legacy) mrf_predict_tree(tree, ctx, NULL) else base
  bshu = array(0, c(N, K, pS + 1)); fshu = matrix(0, N, pS + 1)
  bshu[, , 1] = base$bb; fshu[, 1] = base$fitted
  for(k in seq_len(pS)){
    if(whos.in[k]){
      aff = which(vapply(conds, function(cd) k %in% cd$v, logical(1)))
      accB = matrix(0, N, K); accF = numeric(N)
      for(ii in seq_len(VI.rep)){
        shuf = sample(x = ctx$Sall[, k], replace = FALSE, size = N) #same RNG consumption as the original
        bb = ref$bb; ft = ref$fitted
        if(length(aff)){
          rows0 = unlist(ref$members[aff]); bb[rows0, ] = NA; ft[rows0] = NA
          mem = mrf_leaf_members(conds[aff], ctx$Sall, list(k = k, x = shuf))
          for(q in seq_along(aff)){
            ia = mem[[q]]
            if(length(ia)) { r = mrf_fit_leaf(ia, tree$b0[aff[q], ], ctx, allowed)
            ft[ia] = r$fit; bb[ia, ] = matrix(r$beta, length(ia), K, byrow = TRUE) }
          }
        }
        accB = ((ii - 1) / ii) * accB + bb / ii
        accF = ((ii - 1) / ii) * accF + ft / ii
      }
      bshu[, , k + 1] = accB; fshu[, k + 1] = accF
    }else{ bshu[, , k + 1] = base$bb; fshu[, k + 1] = base$fitted }
  }
  list(betas.shu = bshu, fitted.shu = fshu)
}

## Campos Aleatórios de Markov

# adicionei: 
# n.cores: número de núcleos para o loop sobre as árvores B
# os modos de vi.mode: "legacy" ou "consistent" 


MRF_fast = function(data, x.pos, oos.pos, y.pos = 1, S.pos = 2:ncol(data),
                    minsize = 10, mtry.frac = 1/3, min.leaf.frac.of.x = 1, VI = FALSE,
                    ERT = FALSE, quantile.rate = NULL, S.priority.vec = NULL,
                    random.x = FALSE, howmany.random.x = 1, howmany.keep.best.VI = 20,
                    cheap.look.at.GTVPs = TRUE,
                    prior.var = c(), prior.mean = c(),
                    subsampling.rate = 0.75, rw.regul = 0.75, keep.forest = FALSE,
                    block.size = 12, trend.pos = max(S.pos), trend.push = 1, fast.rw = TRUE,
                    ridge.lambda = 0.1, HRW = 0, B = 50, resampling.opt = 2, printb = TRUE,
                    n.cores = 1, vi.mode = c("legacy", "consistent")){
  vi.mode = match.arg(vi.mode)
  
  if(random.x || resampling.opt >= 3){
    if(!exists("MRF", mode = "function")) stop("resampling.opt>=3 / random.x=TRUE usam o MRF() original: carregue mrf_original.R.")
    message("MRF_fast: caminho nao acelerado -> MRF() original.")
    cl = match.call(); cl[[1]] = as.name("MRF"); cl$n.cores = NULL; cl$vi.mode = NULL
    return(eval(cl, parent.frame()))
  }
  
  #X.posé linear e S.pos está no estado S_t
  data = as.matrix(data); rownames(data) = NULL
  X.pos = x.pos; trend.col = trend.pos
  ET = ERT; ET.rate = quantile.rate
  VI.rep = 10 * VI; lambda = ridge.lambda; BS4.frac = subsampling.rate
  keep.VI = howmany.keep.best.VI; min.leaf.fracz = min.leaf.frac.of.x
  if(!is.null(prior.var)) prior.var = 1 / prior.var #those are in fact implemented in terms of heterogeneous lambdas

  
  #segurança de falhas
  #verificação original
  # adicionei possível erro de matriz
  
  if(length(X.pos) > 10 & (lambda < 0.25 | rw.regul == 0)) stop('For models of this size, consider using a higher ridge lambda (>0.25) and RW regularization.')
  if(min(S.pos) < 1) stop('S.pos cannot be non-postive.')
  if(max(S.pos) > ncol(data)) stop('S.pos specifies variables beyond the limit of your data matrix.')
  if(is.element(y.pos, S.pos)) { S.pos = S.pos[S.pos != y.pos]; print('Beware: foolproof 1 activated. S.pos and y.pos overlap.') }
  if(is.element(y.pos, X.pos)) { X.pos = X.pos[X.pos != y.pos]; print('Beware: foolproof 2 activated. x.pos and y.pos overlap.') }
  if(lambda < 0.0001) { lambda = 0.0001; print('Ridge lambda was too low or negative. Was imposed to be 0.0001') }
  if(!is.null(ET.rate)) if(ET.rate < 0) { ET.rate = 0; print('Quantile Rate/ERT rate was forced to 0, It cannot be negative.') }
  if(length(X.pos) < 1) stop('You need to specificy at least one X.')
  if(length(prior.var) != 0 | length(prior.mean) != 0){
    if(is.null(prior.var) & !is.null(prior.mean)) stop('Need to additionally specificy prior.var.')
    if(is.null(prior.mean) & !is.null(prior.var)) stop('Need to additionally specificy prior.mean.')
    if(length(prior.mean) != (length(X.pos) + 1) | length(prior.var) != (length(X.pos) + 1)) stop('Length of prior vectors do not match that of X.')
  }
  if(length(S.pos) < 5) print('Your S_t is small. Consider augmentating with noisy carbon copies it for better results.')
  if((length(S.pos) * mtry.frac) < 1) stop('Mtry.frac is too small.')
  if(min.leaf.fracz > 2) { min.leaf.fracz = 2; print('min.leaf.frac.of.x forced to 2. Let your trees run deep!') }
  if(is.null(colnames(data))) stop('Data frame must have column names.')
  if((min.leaf.fracz * (length(X.pos) + 1)) < 2){
    min.leaf.fracz = 2 / (length(X.pos) + 1)
    print(paste('Min.leaf.frac.of.x was too low. Thus, it was forced to ', 2 / (length(X.pos) + 1), ' -- a bare minimum. You should consider a higher one.', sep = ''))
  }
  trend.pos = which(S.pos == trend.col)
  if(length(oos.pos) == 0){
    oos.flag = TRUE
    data = rbind(data, matrix(0, 2, ncol(data), dimnames = list(NULL, colnames(data))))
    oos.pos = (nrow(data) - 1):nrow(data); fake.pos = oos.pos
  }else oos.flag = FALSE
  
  #tudo o que não depende da árvore é realizado uma única vez aqui
  #antigamente refazia a padronização e as matrizes de delineamento
  
  N = nrow(data); K = length(X.pos) + 1; pS = length(S.pos)
  in.pos = setdiff(seq_len(N), oos.pos); n_in = length(in.pos)
  std = mrf_standardize(data); Ystd = std$Y
  if(!is.null(prior.mean)){
    prior.mean[-1] = (1 / std$std[y.pos]) * (prior.mean[-1] * std$std[X.pos])
    prior.mean[1]  = (prior.mean[1] - std$mean[y.pos] + std$mean[X.pos] %*% prior.mean[-1]) / std$std[y.pos]
  }
  #evitar o usuário de travar o algoritmo
  
  if(minsize < 2 * min.leaf.fracz * (length(X.pos) + 2)) minsize = 2 * min.leaf.fracz * (length(X.pos) + 1) + 2
  Zall = cbind(1, Ystd[, X.pos, drop = FALSE]); yall = Ystd[, y.pos]
  S.names = colnames(data)[S.pos]
  ctx = list(K = K, lambda = lambda, prior.var = prior.var, prior.mean = prior.mean, HRW = HRW,
             rw.regul = rw.regul, fast.rw = fast.rw, ET = ET, ET.rate = ET.rate,
             min.leaf.fracz = min.leaf.fracz, minsize = minsize, frac = mtry.frac, cons_w = 0.01,
             no.rw = TRUE, oos.pos = oos.pos, n_in = n_in,
             Zall = Zall, yall = yall, Zin = Zall[-oos.pos, , drop = FALSE], yin = yall[-oos.pos],
             Sall = Ystd[, S.pos, drop = FALSE])
  prob.vec = if(is.null(S.priority.vec)) rep(1, pS) else S.priority.vec
  if(trend.push > 1) prob.vec[trend.pos] = trend.push
  cols_needed = unique(c(y.pos, X.pos, S.pos))

  #loop
  
  one_tree = function(b){
    if(printb == TRUE) print(paste('Tree ', b, ' out of ', B, sep = ''))
    if(resampling.opt == 0) rando.vec = 1:n_in #No Bootstrap/sub-sampling. Should not be used for looking at GTVPs.
    if(resampling.opt == 1) rando.vec = sort(sample(x = 1:n_in, replace = FALSE, size = BS4.frac * n_in)) #plain sub-sampling
    if(resampling.opt == 2){ #block sub-sampling. recommended when looking at TVPs.
      groups = sort(base::sample(x = c(1:(n_in / block.size)), size = n_in, replace = TRUE))
      w = rexp(rate = 1, n = n_in / block.size)[groups] + 0.1
      rando.vec = sort(which(w > quantile(w, 1 - BS4.frac)))
    }

    
    noise = 0.00000015 * rnorm(length(rando.vec)) #to avoid redundancy some bug
    sub = Ystd[rando.vec, , drop = FALSE]
    sub[, cols_needed] = sub[, cols_needed] + noise
    tctx = ctx; tctx$rando.vec = rando.vec
    tree = mrf_grow_tree(sub[, S.pos, drop = FALSE], sub[, y.pos], cbind(1, sub[, X.pos, drop = FALSE]),
                         rando.vec, tctx, prob.vec, S.names)
    tab = tree$table
    base = mrf_predict_tree(tree, ctx, allowed = rando.vec)
    un = mrf_unscale(base$fitted, base$bb, std, y.pos, X.pos)
    # contabilização quais observações foram utilizadas para a árvore
    in.out = rep(0, N); in.out[-rando.vec] = 1
    out = list(pred = un$fitted[oos.pos], betas = un$betas, in.out = in.out,
               tree = if(keep.forest) tab else NULL, rando.vec = if(keep.forest) rando.vec else NULL)
    if(VI.rep > 0) out$vi = mrf_vi_tree(tree, base, ctx, VI.rep, vi.mode, rando.vec, S.names)
    out
  }
  
  #Dado que as árvores B são independentes rodei junto, a VI acumulada acumulam ao invés de empilhar (economiza memória) 
  run_chunk = function(bs){
    trees = vector("list", length(bs))
    if(VI.rep > 0) acc = list(bshu = array(0, c(N, K, pS + 1)), bshu.non = array(0, c(N, K, pS + 1)),
                              fsum = matrix(0, N, pS + 1), fcnt = matrix(0, N, pS + 1))
    for(i in seq_along(bs)){
      o = one_tree(bs[i])
      if(VI.rep > 0){
        v = o$vi; oob = o$in.out == 1
        acc$bshu = acc$bshu + v$betas.shu
        nb = v$betas.shu
        #legacy
        nb[!oob, , ] = if(vi.mode == "legacy") NA else 0
        acc$bshu.non = acc$bshu.non + nb
        f = v$fitted.shu[oob, , drop = FALSE]; ok = !is.na(f)
        acc$fsum[oob, ] = acc$fsum[oob, ] + ifelse(ok, f, 0); acc$fcnt[oob, ] = acc$fcnt[oob, ] + ok
        o$vi = NULL
      }
      trees[[i]] = o
    }
    list(trees = trees, acc = if(VI.rep > 0) acc else NULL)
  }
  if(n.cores > 1 && .Platform$OS.type != "windows"){
    # RNG diferente do original por conta de haver vário núcleo (fluxos de L'Ecuyer)
    old = RNGkind(); RNGkind("L'Ecuyer-CMRG"); on.exit(RNGkind(old[1], old[2], old[3]), add = TRUE)
    parallel::mc.reset.stream()
    chunks = split(seq_len(B), cut(seq_len(B), n.cores, labels = FALSE))
    res = parallel::mclapply(chunks, run_chunk, mc.cores = n.cores, mc.set.seed = TRUE)
  }else{
    if(n.cores > 1) warning("n.cores>1 so em Unix/macOS; usando 1 nucleo.")
    res = list(run_chunk(seq_len(B)))
  }
  trees = unlist(lapply(res, `[[`, "trees"), recursive = FALSE)
  
  #junta as previsões, os betas e as médias das árvores
  commitee = do.call(rbind, lapply(trees, `[[`, "pred"))
  avg.pred = colMeans(commitee)
  betas.draws = array(unlist(lapply(trees, `[[`, "betas")), c(N, K, B))
  avg.beta = rowMeans(betas.draws, dims = 2)
  in.out.mat = sapply(trees, `[[`, "in.out")
  howmany.in = rowSums(in.out.mat)
  betas.draws.nonOVF = betas.draws
  for(b in seq_len(B)) betas.draws.nonOVF[in.out.mat[, b] == 0, , b] = NA
  nb0 = betas.draws; nb0[is.na(betas.draws.nonOVF)] = 0
  avg.beta.nonOVF = rowSums(nb0, dims = 2) / howmany.in 
  forest = if(keep.forest) lapply(trees, `[[`, "tree") else NULL
  random.vecs = if(keep.forest) lapply(trees, `[[`, "rando.vec") else NULL
  

  ## Variable Importance (VI) 
  
  
  VI_oob = VI_poos = VI_betas = VI_betas.nonOVF = impZ = NULL
  if(VI.rep > 0){
    acc = Reduce(function(a, b) list(bshu = a$bshu + b$bshu, bshu.non = a$bshu.non + b$bshu.non,
                                     fsum = a$fsum + b$fsum, fcnt = a$fcnt + b$fcnt), lapply(res, `[[`, "acc"))
    betas.shu = acc$bshu / B
    betas.shu.nonOVF = sweep(acc$bshu.non, 1, howmany.in, "/")
    #perdição da média
    oob.pred.shu = acc$fsum / acc$fcnt
    normalized.target = as.matrix(scale(data[, y.pos])[1:nrow(data)])
    crit = function(P, rows) sqrt(colMeans((P[rows, , drop = FALSE] - normalized.target[rows])^2, na.rm = TRUE))
    
    # verfica Raiz do Erro Quadrático Médio para observações fora da amostra
    oob.crit.shu = crit(oob.pred.shu, -oos.pos)
    VI_oob = oob.crit.shu / oob.crit.shu[1] - 1
    which.oob = mrf_which_min_n(-VI_oob[1:pS], keep.VI)
    VI_oob = VI_oob[-1]; which.all = unique(c(which.oob))
    
    # verifica aumento de REQM se tirar k-fatores
    if(!oos.flag){
      poos.crit.shu = crit(oob.pred.shu, oos.pos)
      VI_poos = poos.crit.shu / poos.crit.shu[1] - 1; VI_poos = VI_poos[-1]
      keep.VI = min(c(keep.VI, length(VI_poos)))
      which.poos = mrf_which_min_n(-VI_poos[1:pS], keep.VI)
      which.all = unique(c(which.poos, which.oob))
    }else print('Cannot have VI-OOS without an OOS.')
    
    ## VI para betas
    # com retorno e corte da série temporal de beta
    vib = function(A) { M = matrix(NA_real_, K, pS + 1)
    for(jj in 1:K) M[jj, -1] = sqrt(colMeans((A[, jj, -1, drop = FALSE] - as.vector(A[, jj, 1]))^2, na.rm = TRUE))
    M[, -1, drop = FALSE] }
    VI_betas = vib(betas.shu); VI_betas.nonOVF = vib(betas.shu.nonOVF)
    for(jj in 1:K) which.all = unique(append(which.all, mrf_which_min_n(-VI_betas.nonOVF[jj, ], max(floor(keep.VI / 2), 1))))
    impZ = data[, S.pos[which.all]]
  }
  
  if(oos.flag){
    avg.beta.nonOVF = avg.beta.nonOVF[-fake.pos, , drop = FALSE]; avg.beta = avg.beta[-fake.pos, , drop = FALSE]
    betas.draws.nonOVF = betas.draws.nonOVF[-fake.pos, , , drop = FALSE]; betas.draws = betas.draws[-fake.pos, , , drop = FALSE]
    VI_poos = NULL; avg.pred = NULL; commitee = NULL
    data = data[-fake.pos, , drop = FALSE]
  }
  
  #gráfico de Parâmetros Variáveis no Tempo Generalizados (PVTG)
  #sem loop
  if(cheap.look.at.GTVPs){
    if(B * (1 - BS4.frac) < 30) print('Warning: those bands may have missing values or be innacurate if B is low and subsampling.rate is high.')
    Nd = nrow(data); D = matrix(betas.draws.nonOVF, Nd * K, B)
    qs = if(requireNamespace("matrixStats", quietly = TRUE))
      matrixStats::rowQuantiles(D, probs = c(.16, .84), na.rm = TRUE, drop = FALSE)
    else t(apply(D, 1, quantile, probs = c(.16, .84), na.rm = TRUE))
    bands = array(c(qs[, 1], qs[, 2]), c(Nd, K, 2))
    col.vec = c("#386cb0", "#fdb462", "#7fc97f", "#7fc97f")
    par(mgp = c(2.2, 0.45, 0), tcl = -0.4, mar = c(1.4, 1.4, 1.1, 1.1), mfrow = if(K > 2) c(2, 2) else c(1, 2))
    Zd = cbind(1, data[, X.pos]); keep.OLS = solve(crossprod(Zd)) %*% crossprod(Zd, data[, y.pos])
    for(k in 1:K){
      ts.plot(cbind(avg.beta.nonOVF[, k], bands[, k, 1:2]), col = col.vec[c(1, 2, 2)], lwd = 2 * c(1.1, .75, .75))
      abline(h = keep.OLS[k], col = col.vec[4], lwd = 1)
      if(!oos.flag) { abline(v = min(oos.pos), col = 1, lwd = 0.55, lty = 2)
        legend("bottomleft", c('Posterior Mean', '16 & 84% quantiles', 'OLS', 'OOS start'), col = c(col.vec[c(1, 2, 4)], 1), lty = c(1, 1, 1, 2), cex = .65)
      }else legend("bottomleft", c('Posterior Mean', '16 & 84% quantiles', 'OLS'), col = col.vec[c(1, 2, 4)], lty = c(1, 1, 1), cex = .65)
    }
    par(mfrow = c(1, 1))
  }
  
  list(YandX = data[, c(y.pos, X.pos)], pred.ensemble = commitee, pred = avg.pred,
       betas.raw = avg.beta, important.S = impZ, S.names = colnames(data[, S.pos]),
       betas = avg.beta.nonOVF, betas.draws.raw = betas.draws, betas.draws = betas.draws.nonOVF,
       VI_betas = VI_betas.nonOVF, VI_oob = VI_oob, VI_oos = VI_poos, VI_betas.raw = VI_betas,
       model = list(forest = forest, data = data, regul.lambda = lambda, prior.var = prior.var,
                    prior.mean = prior.mean, rw.regul = rw.regul, HRW = HRW, no.rw.trespassing = TRUE,
                    B = B, random.vecs = random.vecs, y.pos = y.pos, S.pos = S.pos, x.pos = X.pos))
}

## Prediction with a saved forest

# Sem duplicação


pred.given.mrf.fast = function(mrf.output, newdata){
  m = mrf.output$model
  if(is.null(m$forest) || is.null(m$forest[[1]])) stop('Need to activate keep.forest in MRF function first.')
  if(any(colnames(m$data) == "")) stop('Put names on your variables with colnames, for both MRF function and this one.')
  y.pos = m$y.pos; S.pos = m$S.pos; X.pos = m$x.pos; K = length(X.pos) + 1
  std = mrf_standardize(as.matrix(m$data)); Y = std$Y
  newdata = as.matrix(newdata); one = ncol(newdata) == 1
  if(one) { newdata = t(newdata); newdata = rbind(newdata, newdata) }
  rownames(newdata) = NULL
  for(jj in 1:ncol(newdata)) newdata[, jj] = (newdata[, jj] - std$mean[jj + 1]) / std$std[jj + 1]
  n_old = nrow(Y); data2 = rbind(Y, cbind(rep(NA_real_, nrow(newdata)), newdata)); colnames(data2) = colnames(m$data)
  oos.pos = (n_old + 1):nrow(data2)
  pm = m$prior.mean
  if(!is.null(pm)) { pm[-1] = (1 / std$std[y.pos]) * (pm[-1] * std$std[X.pos])
  pm[1] = (pm[1] - std$mean[y.pos] + std$mean[X.pos] %*% pm[-1]) / std$std[y.pos] }
  Zall = cbind(1, data2[, X.pos, drop = FALSE]); yall = data2[, y.pos]
  ctx = list(K = K, lambda = m$regul.lambda, prior.var = m$prior.var, prior.mean = pm, HRW = m$HRW,
             rw.regul = m$rw.regul, oos.pos = oos.pos, n_in = n_old, Zall = Zall, yall = yall,
             Zin = Zall[-oos.pos, , drop = FALSE], yin = yall[-oos.pos], Sall = data2[, S.pos, drop = FALSE])
  avg.pred = 0
  for(b in seq_len(m$B)){
    tree = list(conds = attr(m$forest[[b]], "conds"), b0 = as.matrix(m$forest[[b]][, 5:(4 + K)]))
    base = mrf_predict_tree(tree, ctx, allowed = if(m$no.rw.trespassing) m$random.vecs[[b]] else NULL)
    fitted = base$fitted * std$std[y.pos] + std$mean[y.pos]
    avg.pred = ((b - 1) / b) * avg.pred + (1 / b) * fitted[oos.pos]
  }
  if(one) avg.pred = avg.pred[-2]
  avg.pred
}