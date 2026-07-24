`.sourceCpp_1_DLLInfo` <- dyn.load('/data/work/CapacityAND/sft_package/.Rtmp/RtmpP8dIWh/sourceCpp-x86_64-pc-linux-gnu-1.1.1.1.1/sourcecpp_f24ec633bcf66/sourceCpp_2.so')

nah_nak_eval_rcpp <- Rcpp:::sourceCppFunction(function(rt, cr, q, reverse) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_nah_nak_eval_rcpp')
ucip_eval_rcpp <- Rcpp:::sourceCppFunction(function(rt_list, cr_list, q, reverse) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_ucip_eval_rcpp')
ucip_score_rcpp <- Rcpp:::sourceCppFunction(function(rt_list, cr_list, signs_positive, reverse) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_ucip_score_rcpp')

rm(`.sourceCpp_1_DLLInfo`)
