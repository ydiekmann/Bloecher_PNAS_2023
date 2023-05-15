#
#  *********************************
#  **   Written by Yoan Diekmann  **
#  **    ydiekman@uni-mainz.de    **
#  *********************************
#
#
#
library('ks') # for KDE
#
library('stringr')
# https://cran.r-project.org/web/packages/rcarbon/vignettes/rcarbon.html
#devtools::install_github('ahb108/rcarbon')
library('rcarbon')
## https://github.com/ropensci/c14bazAAR
##if(!require('remotes')) install.packages('remotes')
##remotes::install_github("ropensci/c14bazAAR")
#library('c14bazAAR')
# https://github.com/ISAAKiel/oxcAAR
#install.packages('oxcAAR')
#library('oxcAAR')
#https://github.com/UCL/ADMUR
#install.packages('ADMUR')
library('ADMUR')
#vignette('guide', package = 'ADMUR')
#
#
#
'%ni%' <- Negate('%in%')
#
prob_in_window <- function(cal_grid, start, end) {
    cum_sum = 0
    for (year in start:(end+1)) {
        if (year > cal_grid$calBP[1] | year <= cal_grid$calBP[length(cal_grid$calBP)] |
            year %ni% cal_grid$calBP | (year-1) %ni% cal_grid$calBP) { next }
        #
        cum_sum = cum_sum + (cal_grid[cal_grid$calBP==year, 'PrDens'] + cal_grid[cal_grid$calBP==(year-1), 'PrDens'])/2
    }
    cum_sum
}
#
max_window <- function(cals, indices, window_size) {
    #
    res = data.frame()
    for (year in 4000:(3000+window_size-1)) {
        log_prob = 0
        for (index in indices) {
            log_prob = log_prob + log(prob_in_window(cals[index]$grids[[1]], year, year-window_size+1))
        }
        res = rbind(res, c(year, log_prob))
    }
    colnames(res) = c('year', 'log_prob')
    as.numeric(res[which.max(res$log_prob), ])
}
#
sim_ref_distrib <- function(dat, cals, window_size, nb_sims, normalised) {
    res = max_window(cals, 1:nrow(dat), window_size)
    #
    d = kde(x=dat[, 'uncalBP_sd'])
    #
    log_sums = vector(mode = 'numeric', length = nb_sims)
    for (i in 1:nb_sims) {
        sim_means = uncalibrateCalendarDates(
                                            runif(n = nrow(dat), min = res[1]-window_size+1, max = res[1]), intcal20)
        sim_errs = rkde(fhat=d, n=nrow(dat))
        sim_cals_nn = calibrate(x=sim_means, errors=sim_errs, calCurves='intcal20', normalised = normalised)
        #
        log_sums[i] = sum(log(sapply(X = sim_cals_nn$grids, FUN = 
                                function(x) prob_in_window(cal_grid = x, start = res[1], end = res[1]-window_size+1))))
    }
    return(list('log_sums' = log_sums, 'res' = res))
}
#
analyse <- function(dat, cals, name, window_size, nb_sims, normalised) {
    file_base = paste0(name, '_', window_size, 'yw_nbsims', nb_sims)
    if (!file.exists(paste0(file_base, '.Robj'))) {
        res_ = sim_ref_distrib(dat, cals, window_size, nb_sims, normalised)
        save(res_, file = paste0(file_base, '.Robj'))
    } else {
        load(file = paste0(file_base, '.Robj'))
    }
    log_sums = res_[['log_sums']]
    res = res_[['res']]
    pdf(paste0(file_base, '.pdf'))
    dens = density(log_sums)
    pval = sum(log_sums < res[2]) / nb_sims
    plot(dens, xlim=range(c(dens$x, res[2])),
         main=paste0('p-value=', round(pval, digits = 4)))
    lines(rep(res[2], 2), par('usr')[3:4], col='red')
    dev.off()
}
#
#
#
dat_all = read.csv('C14.csv', header = TRUE, sep = ';', na.strings = c('-', ''))
dat_all = cbind(dat_all, apply(X = str_split(dat_all$uncalBP, '±', simplify = TRUE), MARGIN = 2, FUN = as.numeric))
colnames(dat_all) = c('grave', 'ind', 'sex', 'anthro_age_est', 'C14_MAMS_ID', 'uncalBP',
                      'sigma_1_68.2', 'sigma_2_95.4', 'data_shotgun_capture', 'MT', 'Y', 'uncalBP_mean', 'uncalBP_sd')
#
#
#
if (!file.exists('uncal_BP.pdf')) {
    pdf('uncal_BP.pdf')
    par(mar=c(5,10,4,2))
    perm = order(dat_all[, 'uncalBP_mean'])
    plot(dat_all[perm, 'uncalBP_mean'], 1:length(perm),
         pch=16,
         xlim = rev(range(c(dat_all[perm, 'uncalBP_mean']-1.96*dat_all[perm, 'uncalBP_sd'],
                            dat_all[perm, 'uncalBP_mean']+1.96*dat_all[perm, 'uncalBP_sd']))),
         xlab='year [uncal. BP]', ylab='', yaxt='n')
    segments(dat_all[perm, 'uncalBP_mean']-1.96*dat_all[perm, 'uncalBP_sd'], 1:length(perm),
             dat_all[perm, 'uncalBP_mean']+1.96*dat_all[perm, 'uncalBP_sd'], 1:length(perm))
    axis(2, at = 1:length(perm), labels = paste(perm, dat_all[perm, 'grave'], dat_all[perm, 'ind'], sep=' '),
         cex.axis=0.8, las=2, )
    dev.off()
}
#
#
#
cals_all_norm = calibrate(x=dat_all[, 'uncalBP_mean'], errors=dat_all[, 'uncalBP_sd'], calCurves='intcal20')
if (!file.exists('cal_BP.csv')) {
    summary_all_norm = summary(cals_all_norm)
    summary_all_norm = cbind(dat_all$grave, dat_all$ind, summary_all_norm[2:ncol(summary_all_norm)])
    write.csv(summary_all_norm, file = 'cal_BP.csv')
}
#
#
#
pdf('Burial_6_Ind_1.2sigma.pdf')
plot(cals_all_norm[which(dat_all$grave=='Burial 6' & dat_all$ind=='Ind 1')], HPD=TRUE,credMass=0.954)
dev.off()
summary(cals_all_norm[which(dat_all$grave=='Burial 6' & dat_all$ind=='Ind 1')])
#
pdf('Burial_8_Ind_1.2sigma.pdf')
plot(cals_all_norm[which(dat_all$grave=='Burial 8' & dat_all$ind=='Ind 1')], HPD=TRUE,credMass=0.954)
dev.off()
summary(cals_all_norm[which(dat_all$grave=='Burial 8' & dat_all$ind=='Ind 1')])
#
pdf('Burial_31_Ind_1.2sigma.pdf')
plot(cals_all_norm[which(dat_all$grave=='Burial 31' & dat_all$ind=='Ind 1')], HPD=TRUE,credMass=0.954)
dev.off()
summary(cals_all_norm[which(dat_all$grave=='Burial 31' & dat_all$ind=='Ind 1')])
#
if (!file.exists('cal_BP.pdf')) {
    pdf('cal_BP.pdf')
    multiplot(cals_all_norm, decreasing=TRUE, rescale=TRUE, HPD=TRUE, credMass=0.954)
    dev.off()
}
#
#
#
analyse(dat_all, cals_all_norm, 'all_norm', 1, 10000, normalised=TRUE)
analyse(dat_all, cals_all_norm, 'all_norm', 2, 10000, normalised=TRUE)
analyse(dat_all, cals_all_norm, 'all_norm', 3, 10000, normalised=TRUE)
analyse(dat_all, cals_all_norm, 'all_norm', 10, 10000, normalised=TRUE)
analyse(dat_all, cals_all_norm, 'all_norm', 20, 10000, normalised=TRUE)
analyse(dat_all, cals_all_norm, 'all_norm', 30, 10000, normalised=TRUE)
analyse(dat_all, cals_all_norm, 'all_norm', 40, 10000, normalised=TRUE)
analyse(dat_all, cals_all_norm, 'all_norm', 50, 10000, normalised=TRUE)

