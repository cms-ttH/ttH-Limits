import csv
import math
import os
import re
import ROOT as r
import sys

# b-tag split mode:
B_OFF = 0
B_RATE = 1      # Use rate only
B_SHAPE = 2     # Split into shape and rate
B_CAT_SHAPE = 4 # Split into shape by category and rate

# Define some regular expressions to match samples and signal.  The last
# part of the former expression is a negative look-ahead making sure that
# systematics are not caught by the category definition.

# the sample string identifies things like this
# part1_part2_part3
# part1 = anything without an '_', may end in _obs
# part2 = anything without an '_'
# part3 = anything that is not followed by Up|Down,
#         this makes the sample extraction focus on the non-systematic
#         histograms
sample_re = re.compile(r'([^_]+(?:_obs)?)_([^_]+)_(.*)(?!(?:Up|Down)$)')
signal_re = re.compile(r'(ttH).*')

# This expression splits into 3 groups:
#   category_name:jet_multiplicity:parton_count
# The latter two are optional.
category_re = re.compile(r'(.*?)(?::(\d+)(?::(\d+))?)?$')

# This is the default logging instance (open file or stream)
log = sys.stderr

class IntegralException(Exception):
    """To catch integrals evaluating to zero.  `combine` does not like this
    for systematics.
    """
    pass

def get_ann_systematics(file, discriminant, categories, samples, data_sample="data_obs",
        signal_sample="ttH", is_8_tev=True):
    """For all samples and categories, this functions produces single bin
    systematics, shifts present in only one bin.  These shift will be
    produced when a process has significant expectation and the systematic
    is expected to contribute to the error.
    """
    new_sys = []

    for (c, j, p) in categories:
        data_hist = file.Get("{s}_{d}_{c}".format(s=data_sample, d=discriminant, c=c))
        sig_hist = file.Get("{s}_{d}_{c}".format(s=signal_sample + "125", d=discriminant, c=c))
        bkg_hist = None

        # Build background sum
        for (s, cats) in samples.items():
            s = s + "125" if s == "ttH" else s
            if s in (data_sample, signal_sample) or c not in cats:
                continue

            hist = file.Get("{s}_{d}_{c}".format(s=s, d=discriminant, c=c))
            if bkg_hist:
                bkg_hist.Add(hist)
            else:
                bkg_hist = hist.Clone()

        # Loop over samples for category and find low stats bins
        for (s, cats) in samples.items():
            file_s = s + "125" if s == "ttH" else s

            if c not in cats:
                continue

            hist = file.Get("{s}_{d}_{c}".format(s=file_s, d=discriminant, c=c))

            for b in range(1, hist.GetNbinsX() + 1):
                data = data_hist.GetBinContent(b)
                data_err = data_hist.GetBinError(b)

                sig = sig_hist.GetBinContent(b)
                sig_err = sig_hist.GetBinError(b)

                bkg = bkg_hist.GetBinContent(b)
                bkg_err = bkg_hist.GetBinError(b)

                val = hist.GetBinContent(b)
                val_err = hist.GetBinError(b)

                other_frac = math.sqrt(bkg_err**2 - val_err**2)

                if val < .01 or bkg_err < data_err / 3. or other_frac / bkg_err > .95 \
                        or sig / bkg < .02:
                    continue

                # FIXME Subtract 1 from bin name for comparability with
                # original C macro
                sys_name = "{s}_{c}_{e}_ANNbin{b:d}".format(
                        s=s, c=c, e="8TeV" if is_8_tev else "7TeV", b=b - 1)

                stub = "{s}_{d}_{c}_".format(s=file_s, d=discriminant, c=c)
                hist_up = hist.Clone(stub + sys_name + "Up")
                hist_up.SetBinContent(b, val + val_err)
                hist_down = hist.Clone(stub + sys_name + "Down")
                hist_down.SetBinContent(b, val - val_err)

                file.WriteObject(hist_up, hist_up.GetName())
                file.WriteObject(hist_down, hist_down.GetName())

                new_sys.append((sys_name, "shape",
                    dict([(sam, ("1" if sam == s else "-")) for sam in samples.keys()])))
    return new_sys

def get_integral(file, discriminant, category, sample="data_obs", uncertainty="", fmt="{n:.3f}", throw=False):
    """Get the integral of the histogram specified by the arguments.  The
    first argument has to be an open ROOT TFile.
    """
    if len(uncertainty) > 0:
        uncertainty = "_" + uncertainty
    h = file.Get("{s}_{d}_{c}{u}".format(s=sample, d=discriminant, c=category, u=uncertainty))
    i = h.Integral()
    if i == 0. and throw:
        raise IntegralException("The integral for {s}, {d}, {c}, {u} in {f} is zero!".format(
            s=sample, d=discriminant, c=category, u=uncertainty, f=file.GetName()))
    return fmt.format(n=h.Integral())

def get_samples(file, discriminant):
    """Reads the contents of `file`, an open ROOT TFile, and tries to
    extract available categories per sample.
    """
    samples = {}

    for k in file.GetListOfKeys():
        m = sample_re.match(k.GetName())
        if m:
            sample, disc, cat = m.groups()

            if disc != discriminant:
                continue

            m = signal_re.match(sample)
            if m:
                sample = m.group(1)

            if sample not in samples:
                samples[sample] = set()
            samples[sample].add(cat)
    return samples

def get_systematics(file, overrides={}, rename=lambda u: u, samples=False):
    """Reads `file` and returns a list of (uncertainty, type, {sample:
    value}).

    The parameter `overrides` allows to specify a dict of form
    {uncertainty: value}, and values of "x" n the systematics file are
    replaced by the one specified in the dict.

    A function passed as `rename` allows to alter the uncertainty name,
    _after_ all other processing has happened.

    The parameter `samples` can be set to `True` to obtain the available
    sample names defined in the systematics file.
    """
    reader = csv.DictReader(open(file))
    reader.fieldnames = map(str.strip, reader.fieldnames)
    sys_samples = reader.fieldnames[2:]

    if samples:
        return sys_samples

    sys = []
    # create a list whose elements are the tuples
    # (uncertainty, type, {sample:value})
    # row is a list of the items in the row
    for row in reader:
        unc = row.pop("Uncertainty").strip()
        type = row.pop("Type").strip()
        # `row` is a dictionary with sample names as keys.  Strip spaces
        # from the actual value
        row = dict(map(lambda (k,v): (k, v.strip()), row.items()))
        if unc in overrides:
            row = dict(
                    map(
                        lambda (k,v): (k, overrides[unc] if v == "x" else v),
                        row.items()))
        sys.append((rename(unc), type, row))
    return sys

def parse_btag_mode(s):
    if s == "off":
        return B_OFF
    elif s == "rate":
        return B_RATE
    elif s == "shape":
        return B_SHAPE
    elif s == "category":
        return B_CAT_SHAPE
    raise Exception("Invalid b-tag mode '{m}'".format(m=s))

def split_category_string(s):
    """Split a string of form "category:jets:partons" into components.
    Returns a list of a string and two integers.
    """
    def try_conv(v):
        try:
            return int(v)
        except:
            return v
    return map(try_conv, category_re.match(s).groups())

def split_q2(file, disc, categories):
    """Split Q2-systematics by parton number as specified in `categories`,
    which is a list containing tuples of `(category, jet_multiplicity,
    partons)`.
    """
    for (c, j, p) in categories:
        for d in ("Up", "Down"):
            try:
                oldname = "ttbar_{d}_{c}_Q2scale_ttH_ttbar{dir}".format(d=disc, c=c, dir=d)
                newname = "ttbar_{d}_{c}_Q2scale_ttH_ttbar{p}p{dir}".format(d=disc, c=c, p=p, dir=d)
                hold = file.Get(oldname)
                hnew = hold.Clone(newname)
                file.WriteObject(hnew, newname)
            except:
                log.write("Can't create Q^2 scale shifts for '{c}'\n".format(c=c))

def split_systematics(file, disc, samples, btag_mode=B_CAT_SHAPE):
    """Split b-tag uncertainties:  copy category histogram w/o systematics
    for rates, systematics of form "CMS_eff_bUp" to a shape uncertainty.

    The parameter `samples` should be a dictionary containing the
    categories per sample.

    Returns a list of systematics to be injected into the systematics file.
    """
    done = set()
    new_sys = []
    r.TH1.SetDefaultSumw2()

    if btag_mode == B_OFF:
        return new_sys

    for (s, cats) in samples.items():
        s = s + "125" if s == "ttH" else s
        for c in cats:
            stub = "_".join((s, disc, c))
            orig = file.Get(stub)
            sum = orig.Integral()
            for kind in ("eff", "fake"):
                for dir in ("Up", "Down"):
                    try:
                        # Get rate uncertainty from the shape uncertainty
                        shape_old = file.Get("{s}_CMS_{k}_b{dir}".format(s=stub, k=kind, dir=dir))
                        shape_sum = shape_old.Integral()
                        rate = orig.Clone("{s}_CMS_{k}_bRate{dir}".format(s=stub, k=kind, dir=dir))
                        rate.Scale((shape_sum / sum) if sum > 0 else 1)
                        file.WriteObject(rate, rate.GetName())

                        # Treat shape uncertainties, if desired
                        if btag_mode == B_SHAPE:
                            shape = shape_old.Clone("{s}_CMS_{k}_bShape{dir}".format(s=stub, k=kind, dir=dir))
                        elif btag_mode == B_CAT_SHAPE:
                            shape = shape_old.Clone("{s}_{c}_{k}_bShape{dir}".format(s=stub, c=c, k=kind, dir=dir))
                        else:
                            continue

                        shape_sum = shape.Integral()
                        shape.Scale((sum / shape_sum) if shape_sum > 0 else 1)

                        file.WriteObject(shape, shape.GetName())
                    except:
                        log.write("Can't create b-tag shape uncertainties for '{s}'"
                                "in '{c}'\n".format(s=c, c=c))
                if btag_mode == B_SHAPE and 'all' not in done:
                    new_sys.append((
                        "CMS_{k}_bShape".format(k=kind),
                        "shape",
                        dict([(sam, "1") for sam in samples.keys()])))
                elif btag_mode == B_CAT_SHAPE and c not in done:
                    new_sys.append((
                        "{c}_{k}_bShape".format(c=c, k=kind),
                        "shape",
                        dict([(sam, "1") for sam in samples.keys()])))
            done.add('all')
            done.add(c)
    return new_sys

def write_datacard(file, discriminant, categories, cats, samples, systematics, ofile=log):
    """
    """
    filename = file.GetName()
    observed = map(
            lambda c: get_integral(file, discriminant, c, fmt="{n:.6f}"),
            map(lambda (c, j, p): c, categories))

    if "data_obs" in cats:
        del cats["data_obs"]

    bins = []
    for (n, s) in samples:
        bins += cats[s]

    sprocs = "".join(map(lambda (n, s): (" " + s) * len(cats[s]), samples))
    nprocs = "".join(map(lambda (n, s): (" " + str(n)) * len(cats[s]), samples))

    rates = ["-1"] * len(cats["ttH"])
    for (n, s) in samples[1:]:
        for c in cats[s]:
            rates.append(str(get_integral(file, discriminant, c, s, fmt="{n:.6}")))

    # Print preamble
    ofile.write("""imax * # number of channels
jmax * # number of backgrounds
kmax * # number of nuisance parameters
---------------
bin {c}
observation {o}
---------------
shapes * * {f} $PROCESS_{d}_$CHANNEL $PROCESS_{d}_$CHANNEL_$SYSTEMATIC
shapes ttH * {f} $PROCESS$MASS_{d}_$CHANNEL $PROCESS$MASS_{d}_$CHANNEL_$SYSTEMATIC
---------------
bin {bs}
process {ps}
process {ns}
rate {rs}
---------------
""".format(
    c=" ".join(map(lambda (c, j, p): c, categories)),
        o=" ".join(map(str, observed)),
        d=discriminant,
        f=filename,
        bs=" ".join(bins),
        ps=sprocs,
        ns=nprocs,
        rs=" ".join(rates)))

    active_unc = []
    debugUncert = False
    for (unc, type, vals) in systematics:
        if debugUncert: log.write("-----------------------------------------------")
        if debugUncert: log.write("Considering uncert %s\n   type = %s\n   vals = %s\n" % (unc, type, vals))
        active = False

        ofile.write("{u} {t}".format(u=unc, t=type))
        for (n, s) in samples:
            if debugUncert: log.write("This is sample %s (also %s) \n" % (s,n))
            file_s = s + "125" if s == "ttH" else s
            for c in cats[s]:
                if debugUncert: log.write("This is category %s\n" %c)
                if type == "shape" and vals[s] != "-":
                    try:
                        get_integral(file, discriminant, c, file_s, unc + "Up", throw=True)
                        get_integral(file, discriminant, c, file_s, unc + "Down", throw=True)
                        ofile.write(" " + vals[s])
                        active = True
                    except IntegralException, e:
                        ofile.write(" -")
                        # Print for everything _except_ for b-tag shape or ANN
                        # uncertainties with inappropriate category
                        if not (not unc.startswith(c) and ("bShape" in unc or "ANNbin" in unc)):
                            log.write("Integral zero for {s}, {c}, {u}: disabling "
                                    "systematics\n".format(s=s, c=c, u=unc))
                    except:
                        ofile.write(" -")
                        if not (not unc.startswith(c) and ("bShape" in unc or "ANNbin" in unc)):
                            log.write("Integral not available for {s}, {c}, {u}: disabling "
                                    "systematics\n".format(s=s, c=c, u=unc))
                elif vals[s] in ["-", "1"] and not \
                        (unc == "Q2scale_ttH_V" and s in ("wjets", "zjets")):
                    ofile.write(" " + vals[s])
                    if vals[s] == "1":
                        active = True
                else:
                    # Here is how you add an uncertainty that is applied to only one
                    # category

                    # Q2 scale for wjets/zjets
                    # This uncertainty depends on the jet multiplicty
                    # or parton multiplicity you are considering
                    if unc == "Q2scale_ttH_V" and s in ("wjets", "zjets"):
                        try:
                            # If this category is in your list of categories
                            # find out the parton multiplicty
                            # This pulls out the first entry in the list returned
                            # by the filter
                            # Then pullos out the first entry in the list, then
                            # gets the first value associated with it, which is
                            # the number of jets
                            mult = filter(lambda (cat, j, p): cat == c, categories)[0][1]

                            # The uncertainty will be 10% per jet
                            vals[s] = str(1 + .1 * mult)

                        # exception will occur if this category is not in the list of catgories?
                        # not sure how we get here
                        except:
                            ofile.write(" -")
                            continue
                    # end if unc
                    # if you are doing the NPSF
                    if unc == "NPSF_4j1t":
                        # if you are not in the right category
                        if c != "SS_ge4je1t":
                            # this is not the right category
                            # print a blank
                            ofile.write(" -")
                            continue
                        # end if not in right category
                    # end if NPSF
                    if unc == "NPSF_4j2t":
                        if c != "SS_ge4jge2t":
                            ofile.write(" -")
                            continue
                        # end if
                    # end if
                    if unc == "NPSF_3j2t":
                        if c != "SS_e3jge2t":
                            ofile.write(" -")
                            continue
                    active = True
                    new_val = math.e ** (math.sqrt(math.log(1 + (float(vals[s]) - 1)**2)))
                    ofile.write(" {n:.3f}".format(n=new_val))
        ofile.write("\n")

        if active:
            active_unc.append(unc)

    ofile.write("---------------\n")
    return active_unc

def create_datacard(ifile, ofile, disc, all_categories,
        disabled_systematics=[], btag_mode=B_CAT_SHAPE,
        print_summary=False):
    """Create a datacard for `ifile` (an open ROOT file) using the
    discriminant `disc` and categories, jet multiplicities, parton counts
    defined in `all_categories`.
    """
    print os.path.dirname(__file__)
    sysfile = os.path.join(os.path.dirname(__file__), "systematics.csv")
    all_category_names = map(lambda (c, j, p): c, all_categories)

    is_8_tev = True
    def rename(unc):
        if unc == "lumi":
            return unc + "_8TeV" if is_8_tev else unc + "_7TeV"
        return unc

    # This replaces "x" in the systematics csv file with the values specified
    # for certain uncertainties
    overrides = {
            "lumi": "1.044" if is_8_tev else "1.022",
            "CMS_ttH_eff_lep": "1.04" if is_8_tev else "1.018",
            "CMS_ttH_QCDscale_ttbb": "1.5"}

    # Retrieve list of samples (ordered) from systematics file
    samples = get_systematics(sysfile, samples=True)

    # Get available categories for every sample
    all_samples = get_samples(ifile, disc)

    # Trim previous to the samples defined in the systematics file and filter
    # categories to the ones defined above
    cats = dict(map(
            lambda (k, cs): (k, filter(lambda c: c in cs, all_category_names)),
            filter(
                lambda (k, cs): k in samples,
                all_samples.items())))

    # Enumerate samples for combine
    # nums = dict([(s, n) for (n, s) in enumerate(samples)])
    samples = enumerate(samples)

    systematics = get_systematics(sysfile, overrides=overrides, rename=rename)
    all_uncertainties = map(lambda (u, t, vs): u, systematics)
    systematics = filter(lambda (u, t, vs): u not in disabled_systematics, systematics)
    systematics += split_systematics(ifile, disc, cats, btag_mode)
    systematics += get_ann_systematics(ifile, disc, all_categories, cats, is_8_tev=is_8_tev)

    new_cats = set()
    for (s, cs) in cats.items():
        for c in cs:
            new_cats.add(c)
    new_cats = list(new_cats)
    categories = filter(lambda (c, j, p): c in new_cats, all_categories)

    # keep only essential samples
    samples = filter(lambda (n, s): s in cats, samples)

    split_q2(ifile, disc, all_categories)

    active_unc = write_datacard(ifile, disc, categories, cats, samples, systematics,
            ofile=ofile)

    if not print_summary:
        return

    category_names = map(lambda (c, j, p): c, categories)

    cstring = " ".join(map(lambda c: c if c in category_names else "[" + c + "]",
        all_category_names)).replace("] [", " ")
    sstring = " ".join(map(lambda (n, s): s if s in all_samples.keys() else "[" + s + "]",
        samples)).replace("] [", " ")
    ustring = " ".join(map(lambda u: u if u in active_unc else "[" + u + "]",
        all_uncertainties)).replace("] [", " ")

    log.write("""
Executive Summary
=================

Included in Datacard
--------------------
Categories:    {cs}
Samples:       {ss}
Uncertainties: {us}
~~~
Disabled objects present within []s
""".format(cs=cstring, ss=sstring, us=ustring))
