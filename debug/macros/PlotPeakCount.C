#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TH1.h>
#include <TLatex.h>
#include <TString.h>

void PlotPeakCount(
    const char *input = "batches/run_004/Combine/Tracker-step2-ana.root",
    const char *output = "batches/run_004/Combine/Peakcount.png")
{
    TFile file(input, "READ");
    if (file.IsZombie()) {
        Error("PlotPeakCount", "Cannot open %s", input);
        return;
    }

    auto *source = dynamic_cast<TH1 *>(file.Get("PeakCount_all_All"));
    if (!source) {
        Error("PlotPeakCount", "PeakCount_all_All not found in %s", input);
        return;
    }

    auto *hist = dynamic_cast<TH1 *>(source->Clone("peak_count_display"));
    hist->SetDirectory(nullptr);

    const double total = hist->Integral(1, hist->GetNbinsX());
    const double twoOrMore =
        hist->Integral(hist->FindBin(1.5), hist->GetNbinsX());
    const double fraction = total > 0 ? 100.0 * twoOrMore / total : 0.0;

    TCanvas canvas("peak_count_canvas", "Significant waveform peaks", 1000, 750);
    canvas.SetLogy();
    canvas.SetLeftMargin(.12);
    canvas.SetRightMargin(.04);
    canvas.SetBottomMargin(.13);
    canvas.SetTopMargin(.23);

    hist->SetTitle(";Significant peak count;Waveforms");
    hist->SetStats(false);
    hist->SetLineColor(kBlue + 1);
    hist->SetLineWidth(2);
    hist->SetFillColor(kAzure - 9);
    hist->GetXaxis()->SetRangeUser(-0.5, 6.5);
    hist->SetMinimum(0.5);
    hist->SetMaximum(hist->GetMaximum() * 5.0);
    hist->Draw("HIST");

    TLatex text;
    text.SetNDC();
    text.SetTextFont(42);
    text.SetTextSize(.032);
    text.DrawLatex(.12, .95, "All: significant waveform peaks");
    text.DrawLatex(.12, .90, "All mapped 30-sample waveforms");
    text.DrawLatex(.12, .85, Form("Waveforms = %.0f", total));
    text.DrawLatex(.56, .85,
                   Form("N_{peak} #geq 2: %.0f (%.2f%%)",
                        twoOrMore, fraction));

    text.SetTextSize(.022);
    text.DrawLatex(
        .12, .035,
        "Peak: above configured run-noise threshold, maximum within #pm3 samples; not a cluster count");

    canvas.SaveAs(output);
}