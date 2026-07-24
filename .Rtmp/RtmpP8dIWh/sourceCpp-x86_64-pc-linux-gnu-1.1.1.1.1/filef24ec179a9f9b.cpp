
        #include <Rcpp.h>
        #include <vector>
        #include <algorithm>
        using namespace Rcpp;

        struct RTCR { double rt; bool cr; };

        static NumericMatrix eval_hazard(NumericVector rt, LogicalVector cr,
                                         NumericVector q, bool reverse) {
          NumericMatrix out(q.size(), 2);
          std::vector<RTCR> v;
          for (int i = 0; i < rt.size(); ++i) {
            if (!R_finite(rt[i])) continue;
            bool event = i < cr.size() && cr[i] != NA_LOGICAL && cr[i] == TRUE;
            v.push_back(RTCR{rt[i], event});
          }
          std::sort(v.begin(), v.end(), [](const RTCR& a, const RTCR& b) {
            return a.rt < b.rt;
          });
          if (v.empty()) return out;
          std::vector<double> event_rt, cum_h, cum_v;
          double h = 0.0, vv = 0.0;
          for (int i = 0; i < v.size(); ++i) {
            if (!v[i].cr) continue;
            double risk = reverse ? (i + 1.0) : (v.size() - i);
            double inc = 1.0 / risk;
            h += inc; vv += inc * inc;
            event_rt.push_back(v[i].rt);
            cum_h.push_back(h); cum_v.push_back(vv);
          }
          if (event_rt.empty()) return out;
          if (!reverse) {
            for (int j = 0; j < q.size(); ++j) {
              if (!R_finite(q[j])) continue;
              int k = std::upper_bound(event_rt.begin(), event_rt.end(), q[j]) - event_rt.begin() - 1;
              if (k >= 0) { out(j, 0) = cum_h[k]; out(j, 1) = cum_v[k]; }
            }
          } else {
            std::vector<double> tail_h(event_rt.size()), tail_v(event_rt.size());
            double th = 0.0, tv = 0.0;
            for (int i = event_rt.size() - 1; i >= 0; --i) {
              double inc_h = i == 0 ? cum_h[0] : cum_h[i] - cum_h[i - 1];
              double inc_v = i == 0 ? cum_v[0] : cum_v[i] - cum_v[i - 1];
              th += inc_h; tv += inc_v;
              tail_h[i] = -th; tail_v[i] = tv;
            }
            for (int j = 0; j < q.size(); ++j) {
              if (!R_finite(q[j])) continue;
              int k = std::lower_bound(event_rt.begin(), event_rt.end(), q[j]) - event_rt.begin();
              if (k < event_rt.size()) { out(j, 0) = tail_h[k]; out(j, 1) = tail_v[k]; }
            }
          }
          return out;
        }

        // [[Rcpp::export]]
        NumericMatrix nah_nak_eval_rcpp(NumericVector rt, LogicalVector cr,
                                        NumericVector q, bool reverse) {
          return eval_hazard(rt, cr, q, reverse);
        }

        // [[Rcpp::export]]
        NumericMatrix ucip_eval_rcpp(List rt_list, List cr_list,
                                     NumericVector q, bool reverse) {
          NumericMatrix out(q.size(), 2);
          for (int i = 0; i < rt_list.size(); ++i) {
            NumericVector rt = as<NumericVector>(rt_list[i]);
            LogicalVector cr = i < cr_list.size() && !Rf_isNull(cr_list[i])
              ? as<LogicalVector>(cr_list[i]) : LogicalVector(rt.size(), true);
            NumericMatrix one = eval_hazard(rt, cr, q, reverse);
            for (int j = 0; j < q.size(); ++j) {
              out(j, 0) += one(j, 0); out(j, 1) += one(j, 1);
            }
          }
          return out;
        }

        // [[Rcpp::export]]
        NumericVector ucip_score_rcpp(List rt_list, List cr_list,
                                      LogicalVector signs_positive, bool reverse) {
          int ncond = rt_list.size();
          std::vector<NumericVector> rts;
          std::vector<LogicalVector> crs;
          std::vector<double> pooled;
          for (int i = 0; i < ncond; ++i) {
            NumericVector rt = as<NumericVector>(rt_list[i]);
            LogicalVector cr = i < cr_list.size() && !Rf_isNull(cr_list[i])
              ? as<LogicalVector>(cr_list[i]) : LogicalVector(rt.size(), true);
            rts.push_back(rt); crs.push_back(cr);
            for (int j = 0; j < rt.size(); ++j)
              if (R_finite(rt[j])) pooled.push_back(rt[j]);
          }
          NumericVector out = NumericVector::create(
            _["capacity_numerator"] = 0.0, _["capacity_denominator"] = 0.0,
            _["numerator_variance"] = 0.0, _["denominator_variance"] = 0.0);
          if (pooled.empty()) return out;
          // Pooled event grid keeps duplicates so the risk-set columns and the
          // rightmost.closed tie handling match the R implementation exactly.
          std::sort(pooled.begin(), pooled.end());
          int T = pooled.size();

          std::vector< std::vector<double> > sorted_rt(ncond);
          for (int i = 0; i < ncond; ++i) {
            for (int j = 0; j < rts[i].size(); ++j)
              if (R_finite(rts[i][j])) sorted_rt[i].push_back(rts[i][j]);
            std::sort(sorted_rt[i].begin(), sorted_rt[i].end());
          }

          // risk_i(t): at-risk count (rt >= t) for OR, cumulative count
          // (rt <= t) for the reverse (AND/STST) hazard.
          std::vector< std::vector<double> > risk(ncond, std::vector<double>(T, 0.0));
          std::vector<double> total(T, 0.0), weight(T, 0.0);
          for (int i = 0; i < ncond; ++i) {
            const std::vector<double>& s = sorted_rt[i];
            int m = s.size();
            for (int k = 0; k < T; ++k) {
              double t = pooled[k];
              double r = reverse
                ? (double)(std::upper_bound(s.begin(), s.end(), t) - s.begin())
                : (double)(m - (std::lower_bound(s.begin(), s.end(), t) - s.begin()));
              risk[i][k] = r;
              total[k] += r;
            }
          }
          // Houpt-Townsend weight: redundant risk (row 0) times the summed
          // single-target risk, divided by the total risk.
          for (int k = 0; k < T; ++k) {
            double r0 = risk[0][k];
            double w = total[k] > 0 ? r0 * (total[k] - r0) / total[k] : 0.0;
            weight[k] = R_finite(w) ? w : 0.0;
          }

          for (int i = 0; i < ncond; ++i) {
            std::vector<double> ev;
            NumericVector rt = rts[i];
            LogicalVector cr = crs[i];
            for (int j = 0; j < rt.size(); ++j) {
              if (!R_finite(rt[j])) continue;
              bool event = j < cr.size() && cr[j] != NA_LOGICAL && cr[j] == TRUE;
              if (event) ev.push_back(rt[j]);
            }
            std::sort(ev.begin(), ev.end());
            double component = 0.0, component_var = 0.0;
            for (int e = 0; e < (int)ev.size(); ++e) {
              double t = ev[e];
              // findInterval(t, pooled, rightmost.closed = TRUE), 0-based.
              int k = (std::upper_bound(pooled.begin(), pooled.end(), t) - pooled.begin()) - 1;
              if (k == T - 1 && t == pooled[T - 1]) k = T - 2;
              if (k < 0) continue;
              double ri = risk[i][k];
              if (ri <= 0) continue;
              double term = weight[k] / ri;
              if (!R_finite(term)) continue;
              component += term;
              component_var += term * term;
            }
            if (signs_positive[i] == TRUE) {
              out[0] += component; out[2] += component_var;
            } else {
              out[1] += component; out[3] += component_var;
            }
          }
          return out;
        }
