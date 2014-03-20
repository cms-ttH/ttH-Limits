#!/usr/bin/env python
"""Plot ratio comparisons

Currently, all customization and providing of data has to be done within
this script.  See the following two variables to get started (`channels`
contains lists y-axis label and filename to read the fit result from)
"""
outfile  = "limit_mu_cmp_all.pdf"
outfile2 = "limit_mu_cmp_all.png"
channels = [
    ('Combination', 'mlfit.root'),
    ('Same-Sign 2l', 'mlfit.root'),
    ('3l', 'mlfit.root'),
    ('4l', 'mlfit.root'),
    ('#tau_{h}#tau_{h}', 'mlfit.root'),
    ('b#bar{b}', 'mlfit.root'),
    ('#gamma#gamma', 'mlfit.root'),
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

dummy = r.TH1F("dummy", ";Best fit #sigma/#sigma_{SM} at m_{H} = 125.6 GeV;", 100, -100, 100)
dummy.GetXaxis().SetTitleSize(0.05)
dummy.GetYaxis().SetRangeUser(-0.5, len(channels) - 0.5)
dummy.GetYaxis().Set(len(channels), -0.5, len(channels) - 0.5)
dummy.GetYaxis().SetLabelSize(0.06)

xMin = 0
xMax = 0
tableLines = []
for (n, (label, fileName)) in enumerate(channels):

    ## Open the ROOT file and pull out the fit result
    #file = r.TFile.Open(fileName)
    #fr = file.Get('fit_s')
    #pars = r.RooArgSet(fr.floatParsFinal())
    #mu = pars.find('r')

    #muVal = mu.getVal()
    #muErrLo = mu.getErrorLo()
    #muErrHi = mu.getErrorHi()

    if label == 'Four-Lepton':
        muVal = -4.2
        muErrLo = -1.8
        muErrHi = 4.4
        print 'WARNING: Hand-setting 4l results to match PAS!!!'

    if label == 'Combination': 
        muVal = 2.79
        muErrLo = -0.92
        muErrHi = 1.02
        #print 'FAKE values for now'
        
    if label == 'Same-Sign 2l':
        muVal = 5.31
        muErrLo = -1.82
        muErrHi = 2.14
        
    if label == '3l':
        muVal = 3.13
        muErrLo = -1.96
        muErrHi = 2.43
        
    if label == '4l':
        muVal = -4.71
        muErrLo = -1.29
        muErrHi = 5.03
        
    if label == 'Hadronic #tau#tau':
        muVal = -1.34
        muErrLo = -5.47
        muErrHi = 6.33
        
    if label == 'b#bar{b}':
        muVal = 0.67
        muErrLo = -1.87
        muErrHi = 1.86

    if label == '#gamma#gamma':
        muVal = 2.69
        muErrLo = -1.83
        muErrHi = 2.60

    if muVal+muErrLo < xMin:
        xMin = muVal+muErrLo

    if muVal+muErrHi > xMax:
        xMax = muVal+muErrHi

    print label,muVal,muErrLo,muErrHi

    tableLines.append('%s & %.1f^{%+.1f}_{%+.1f} \\\\' % (label, muVal, muErrHi, muErrLo))
    
    dummy.GetYaxis().SetBinLabel(n + 1, label)
    ratio.SetPoint(n, muVal, n)
    ratio.SetPointError(n, -muErrLo, muErrHi, 0, 0)

print xMin, xMax
dummy.GetXaxis().SetRangeUser(1.2*xMin,1.2*xMax)
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


# Hack in a line for 4l where it is limited
line2 = r.TLine()
line2.SetLineStyle(r.kDashed)
line2.SetLineWidth(2)
line2.DrawLine(-6,2.75,-6,3.25)


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
tex.DrawLatex(0.2, 0.9, "CMS")
#tex.DrawLatex(0.2, 0.9, "CMS Preliminary")
tex.DrawLatex(0.5, 0.9, "#sqrt{s} = 7 TeV, L = 5.0 fb^{-1}; #sqrt{s} = 8 TeV, L = 19.7 fb^{-1}")

div = r.TLine()
div.SetLineWidth(2)
div.DrawLine(dummy.GetXaxis().GetBinLowEdge(dummy.GetXaxis().GetFirst()),
             0.5,
             dummy.GetXaxis().GetBinUpEdge(dummy.GetXaxis().GetLast()),
             0.5)

canvas.GetPad(0).RedrawAxis()
canvas.RedrawAxis()

canvas.SaveAs(outfile)
canvas.SaveAs(outfile2)

print '\n\n-----\n'

for l in tableLines:
    print l
