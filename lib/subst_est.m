function ps = subst_est(ps_est, est_ind, ps_give)
ps = ps_give;
ps(est_ind) = ps_est;