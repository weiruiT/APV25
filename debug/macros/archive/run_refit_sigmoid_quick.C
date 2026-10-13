#include "analyze_refit_sigmoid_variance.C"

void run_refit_sigmoid_quick()
{
  analyze_refit_sigmoid_variance("Combine/Tracker-step2-ana.root",
                                 "refit_sigmoid_variance_quick.root",
                                 "refit_sigmoid_variance_quick_summary.txt",
                                 10);
}
