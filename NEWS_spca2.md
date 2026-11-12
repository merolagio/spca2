# spca2 0.2.0
- spca() returns lists with weights and indices in the order they were selected
- new_spca is more efficient and takes a pca() output to increase speed

# spca2 0.1.0 (changes from spca 1.1.5)
- replaced the power method with Spectra truncated eigendecomposition
   (via RSpectra headers, `LinkingTo: RSpectra`); made it the default in
   `spca()` (`pm_weights`, `pm_varsel`) and `pca()` (`pm`). Dense
   eigendecomposition still available with `pm* = FALSE`
- variable selection uses dense eigendecomposition for small subsets
   (< 25 variables), Spectra for larger ones
- `eps_pm*` and `maxiter_pm*` arguments kept for compatibility but ignored
- `pca()` computes only the `n_comps` leading eigenpairs; for data matrices
   the cross product is not formed
- added function `eigvals()`: all eigenvalues, vexp and cvexp, with scree
   plot and qq-plot; new class `eigvals` with `scree_plot()` and
   `mp_qqplot()` methods
- `scree_plot()` warns and `mp_qqplot()` stops when a `pca` object holds
   only the leading eigenvalues
- data matrices: covariance computed as X'X/(n - 1); `eigenvalues` and
   `tot_var` are now on the covariance/correlation scale
