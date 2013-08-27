#!/usr/bin/env python
"""Plot ratio comparisons

Currently, all customization and providing of data has to be done within
this script.  See the following two variables to get started (`channels`
contains lists of alias, y-axis label, ratio, error down, error up.)
"""
outfile = "7TeV_limit_mu_cmp.pdf"
channels = [
        ("TAU", "Hadronic #tau#tau", -0.733308, 5.24428, 6.13693),
        ("OSDIL", "Dilepton", 1.2266, 4.68643, 4.19684),
        # ("7TeV", "7TeV LJ+DIL", -2.77342, 4.87313, 4.13827),
        ("7TeV", "7TeV LJ+DIL", -2.82274, 4.9184, 4.15527),
        ("2photon", "#gamma#gamma", 0.206946, 1.46224, 2.18212),
        ("LJ", "Lepton + Jets", -0.103221, 2.57988, 2.53175),
#        ("COMBI", "Combination", 0.46764, 1.30513, 1.62569), #8TeV
        ("COMBI", "Combination", 0.740757, 1.30386, 1.3407), #all
#        ("COMBI", "Combination", 0.851967, 2.40582, 2.47273) #8TeV bb/tautau
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
legend = r.TLegend(0.6, 0.7, 0.9, 0.87)

ratio = r.TGraphAsymmErrors(len(channels))
ratio.SetMarkerStyle(21)
ratio.SetMarkerSize(2)
ratio.SetLineWidth(4)
ratio.SetLineColor(r.kRed)

dummy = r.TH1F("dummy", ";Best fit #sigma/#sigma_{SM} at m_{H} = 125 GeV;", 100, -7, 7)
dummy.GetXaxis().SetTitleSize(0.05)
dummy.GetYaxis().SetRangeUser(-0.5, len(channels) - 0.5)
dummy.GetYaxis().Set(len(channels), -0.5, len(channels) - 0.5)
dummy.GetYaxis().SetLabelSize(0.06)
for (n, (chan, label, c, m, p)) in enumerate(channels):
    dummy.GetYaxis().SetBinLabel(n + 1, label)
    ratio.SetPoint(n, c, n)
    ratio.SetPointError(n, m, p, 0, 0)

dummy.Draw("axis")
# box = r.TBox()
# box.SetFillColor(r.kGreen)
# box.SetFillStyle(1001)
# box.DrawBox(0.851967 - 2.40582, -0.5, 0.851967 + 2.47273, len(channels) - 0.5)
line = r.TLine()
line.SetLineWidth(2)
# line.DrawLine(0.851967, -0.5, 0.851967, len(channels) - 0.5)
line.DrawLine(1, -0.5, 1, len(channels) - 0.5)
ratio.Draw("same p")

exp = r.TGraphAsymmErrors(len(channels))
line = r.TLine()

# obs.SetMarkerStyle(25)
# obs.SetMarkerColor(r.kMagenta)
# obs.SetMarkerSize(2)
# obs.SetLineStyle(r.kSolid)
# obs.SetLineWidth(2)
# obs.Draw("P same")
# exp.SetMarkerStyle(8)
# exp.SetMarkerSize(1)
# exp.SetLineStyle(r.kDashed)
# exp.SetLineWidth(2)
# exp.Draw("P same")

# exp.SetFillColor(r.kGreen)
# legend.AddEntry(exp, "Expected #pm 1 #sigma", "fl")
# exp2 = exp.Clone()
# exp2.SetFillColor(r.kYellow)
# legend.AddEntry(exp2, "Expected #pm 2 #sigma", "fl")
# legend.AddEntry(obs, "Observed", "l")
# legend.SetFillColor(0)
# legend.SetTextFont(42)
# legend.SetTextSize(0.05)
# legend.Draw()

tex = r.TLatex()
tex.SetNDC()
tex.SetTextFont(42)
tex.SetTextSize(0.034)
tex.DrawLatex(0.2, 0.9, "CMS Preliminary")
tex.DrawLatex(0.5, 0.9, "#sqrt{s} = 7 TeV, L = 5.0 fb^{-1}; #sqrt{s} = 8 TeV, L = 19.5 fb^{-1}")

canvas.GetPad(0).RedrawAxis()

canvas.SaveAs(outfile)
