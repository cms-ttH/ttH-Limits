#!/usr/bin/env python
"""Plot a comparison of limits

All available configuration is contained in the two variables below.  The
list `channels` contains pairs of alias and y-axis labels, where the data
for the plot is extracted from files with the name
`limit_card_{alias}_125.7.log`, which are expected to contain regular
`combine` output.
"""
outfile = "limit_cmp_all.pdf"
channels = [
        ("combination", "Combination"),
        ("2lss", "Same-Sign 2l"),
        ("3l", "3l"),
        ("4l", "4l"),
        ("tt", "Hadronic #tau#tau"),
        ("bb","b#bar{b}"),
##         ("bb_8TeV","b#bar{b} (8 TeV)"),
##         ("bb_7TeV","b#bar{b} (7 TeV)"),
        ("photons", "#gamma#gamma"),
        ]

import ROOT as r

r.gROOT.SetBatch()
r.gROOT.SetStyle("Modern")
r.gStyle.SetOptStat(0)

def limit(s):
    return float(s.split()[-1])

canvas = r.TCanvas()
canvas.SetBottomMargin(.15)
canvas.SetLeftMargin(.2)
canvas.SetRightMargin(.05)
canvas.SetTopMargin(.12)
legend = r.TLegend(0.72, 0.68, 0.9, 0.87)
legend.SetBorderSize(0)

dummy = r.TH1F("dummy", ";95% CL limit on #sigma/#sigma_{SM} at m_{H} = 125.7 GeV;", 100, 0.8, 80.5)
dummy.GetXaxis().SetTitleSize(0.05)
dummy.GetYaxis().SetRangeUser(-0.5, len(channels) - 0.5)
dummy.GetYaxis().Set(len(channels), -0.5, len(channels) - 0.5)
dummy.GetYaxis().SetLabelSize(0.06)
for (n, (chan, label)) in enumerate(channels):
    dummy.GetYaxis().SetBinLabel(n + 1, label)

dummy.Draw("axis")

obs = r.TGraph(len(channels))
exp = r.TGraphAsymmErrors(len(channels))
expSig = r.TGraphAsymmErrors(len(channels))
box = r.TBox()
line = r.TLine()

tableLines = []
tableLines2 = []
for (n, (chan, label)) in enumerate(channels):
    with open("limit_{c}_125.7.log".format(c=chan)) as f:
        lines = f.readlines()
        obs_line = filter(lambda s: s.startswith("Observed Limit"), lines)[0]
        exp_lines = filter(lambda s: s.startswith("Expected "), lines)

        obs.SetPoint(n, limit(obs_line), n)
        exp.SetPoint(n, limit(exp_lines[2]), n)
        exp.SetPointError(n, limit(exp_lines[2]) - limit(exp_lines[1]), limit(exp_lines[3]) - limit(exp_lines[2]), 0, 0)

        o = limit(obs_line)
        xl2 = limit(exp_lines[0])
        xl1 = limit(exp_lines[1])
        x = limit(exp_lines[2])
        xh1 = limit(exp_lines[3])
        xh2 = limit(exp_lines[4])

        box.SetFillStyle(1001)
        box.SetFillColor(r.kYellow)
        box.DrawBox(xl2, n - 0.5, xh2, n + 0.5)
        box.SetFillColor(r.kGreen)
        box.DrawBox(xl1, n - 0.5, xh1, n + 0.5)

        line.SetLineColor(r.kBlack)
        line.SetLineWidth(2)
        line.SetLineStyle(r.kDashed)
        line.DrawLine(x, n - 0.5, x, n + 0.5)
        line.SetLineStyle(r.kSolid)
        line.DrawLine(o, n - 0.5, o, n + 0.5)

    with open("limit_{c}_125.7_sig.log".format(c=chan)) as f:
        lines = f.readlines()
        exp_line_sig = filter(lambda s: s.startswith("Observed Limit"), lines)[0]

        xs = limit(exp_line_sig)

        line.SetLineWidth(2)
        line.SetLineStyle(r.kDashed)
        line.SetLineColor(r.kRed)
        line.DrawLine(xs, n - 0.5, xs, n + 0.5)

    tableLines.append('%s & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f \\\\ ' % (label, o, xs, xl2, xl1, x, xh1, xh2))
    tableLines2.append('%s & %.1f & %.1f & %.1f \\\\' % (label, o, xs, x))


for l in reversed(tableLines):
    print l

print '\n-----\n'

for l in reversed(tableLines2):
    print l


obs.SetMarkerStyle(21)
obs.SetMarkerColor(r.kBlack)
obs.SetMarkerSize(1)
obs.SetLineStyle(r.kSolid)
obs.SetLineWidth(2)
obs.Draw("P same")
exp.SetMarkerStyle(8)
exp.SetMarkerSize(1)
exp.SetLineStyle(r.kDashed)
exp.SetLineWidth(2)
# exp.Draw("P same")
expSig.SetMarkerStyle(25)
expSig.SetMarkerColor(r.kMagenta)
expSig.SetMarkerSize(2)
expSig.SetLineStyle(r.kDashed)
expSig.SetLineColor(r.kRed)
expSig.SetLineWidth(2)

exp.SetFillColor(r.kGreen)
legend.AddEntry(exp, "Expected #pm 1 #sigma", "fl")
exp2 = exp.Clone()
exp2.SetFillColor(r.kYellow)
legend.AddEntry(exp2, "Expected #pm 2 #sigma", "fl")
legend.AddEntry(expSig,"Expected (sig. inj.)", "l")
legend.AddEntry(obs, "Observed", "lp")
legend.SetFillColor(0)
legend.SetTextFont(42)
legend.SetTextSize(0.035)
legend.Draw()

tex = r.TLatex()
tex.SetNDC()
tex.SetTextFont(42)
tex.SetTextSize(0.034)
tex.DrawLatex(0.2, 0.9, "CMS Preliminary")
tex.DrawLatex(0.5, 0.9, "#sqrt{s} = 7 TeV, L = 5.0 fb^{-1}; #sqrt{s} = 8 TeV, L = 19.5 fb^{-1}")

div = r.TLine()
div.SetLineWidth(2)
div.DrawLine(0.8,0.5,80.5,0.5)

canvas.GetPad(0).SetLogx()
canvas.GetPad(0).RedrawAxis()

canvas.SaveAs(outfile)
