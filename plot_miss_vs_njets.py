#!/bin/env python

# Plot the number of hdw-sim mismatches vs. N jets
#
# davide.gerbaudo@gmail.com
# Jul 2017

import itertools
import optparse
import re
import os
import inspect
import string
import sys
from collections import defaultdict
from pprint import pprint
from math import pi, sqrt

cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../L1TopoValidation/L1TopoCheck/python/")))
if cmd_subfolder not in sys.path:
     sys.path.insert(0, cmd_subfolder)


import ROOT as R
R.PyConfig.IgnoreCommandLineOptions = True # don't let root steal your cmd-line options
R.gROOT.SetBatch(1)                        # go batch!
R.gErrorIgnoreLevel = 9999                 # suppress messages about missing dict
                                           # (can't get rid of the 'duplicate' ones?)
R.gROOT.ProcessLine('#include "L1TopoCheck/TriggerBits.h"')
R.gROOT.ProcessLine('#include "L1TopoCheck/AlgorithmBits.h"')
R.gROOT.Macro('$ROOTCOREDIR/scripts/load_packages.C')
EmTob = R.L1Topo.offline.EmTOB
JetTob = R.L1Topo.offline.JetTOB
# MuonTob = R.MuonTOB
# MuonTob = R.L1Topo.EnhancedMuonTOB
MuonTob = R.L1Topo.MuonTOB
TauTob = R.L1Topo.offline.TauTOB

import utils

def main():
    usage = ("Usage : %prog [options] filename"
             "\n Examples :"
             "\n %prog  -v tmptrig.root"
             )

    parser = optparse.OptionParser(usage = usage)
    parser.add_option('-n', '--num-events', default=None, type=int, help='number of events to process (default all)')
    parser.add_option('-s', '--skip-events', default=None, type=int, help='number of events to skip (default none)')
    parser.add_option('--skip-overflow', action='store_true', help='do not use events with overflow bits')
    parser.add_option('-v', '--verbose', default=False, action='store_true')
    parser.add_option('-d', '--debug', default=False, action='store_true')
    parser.add_option('-t', '--treename', default='trig')

    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")
    verbose = options.verbose
    debug = options.debug
    if verbose:
        utils.print_running_conditions(parser, options)

    input_filenames = utils.read_filename_arguments(args[0], options)
    if verbose:
        print 'Input files:'
        print '\n'.join(input_filenames)
    chain = R.TChain(options.treename)
    for input_filename in input_filenames:
        chain.Add(input_filename)  #chain beomes an array with the various files .root introduced
    num_available = chain.GetEntries()
    num_skip = options.skip_events
    num_toprocess = number_of_entries_to_process(num_available, options)
    if verbose:
        print "About to process %s (out of %d) entries: " % (num_toprocess, num_available)
    # test: print all branches available in the tree
    # print 'branches:'
    # print 'list of branches ',chain.GetListOfBranches()
    # print '\n'.join([k.GetName() for k in chain.GetListOfBranches()])
    # print
    # made from Murrough's suggestion (email "Jet sorting and
    # (possibly) different FPGA U1/U2 jet fibre order") with
    # egrep '(.*INVM.*-AJ30s6-AJ20s6|HT.*-J.*s5.ETA.*|0DET10-Js1-Js2|10MINDPHI-J.*s2-XE50.*|.*MINDPHI-AJj20s6-XE0)' ../L1TopoValidation/L1TopoCheck/share/algorithm_labels.txt
    algorithm_names = ['900INVM9999-AJ30s6-AJ20s6',
                       '800INVM9999-AJ30s6-AJ20s6',
                       '700INVM9999-AJ30s6-AJ20s6',
                       '500INVM9999-AJ30s6-AJ20s6',
                       '400INVM9999-AJ30s6-AJ20s6',
                       '300INVM9999-AJ30s6-AJ20s6',
                       '200INVM9999-AJ30s6-AJ20s6',
                       '10MINDPHI-J20s2-XE50',
                       '100INVM9999-AJ30s6-AJ20s6',
                       'HT150-J20s5.ETA31',
                       'HT190-J15s5.ETA21',
                       ]
    jet_binning   = (13 , -0.5, 12.5)
    tprof_mismatches = {k : R.TProfile('n_mismatches_vs_n_jets_'+k, k+';N JetTOB; N hdw!=sim', *jet_binning)
                        for k in algorithm_names}
    hist_mismatches = {k : R.TH1F('mismatches_vs_n_jets_'+k, 'Mismatches '+k+';N JetTOB; N hdw!=sim1', *jet_binning)
                       for k in algorithm_names}
    hist_hdwpasses = {k : R.TH1F('hdwpass_vs_n_jets_'+k, k+';N JetTOB; N hdw==1', *jet_binning)
                      for k in algorithm_names}
    sim_bit_number = R.L1Topo.TOPOSIM # see L1TopoCheck/AlgorithmBits.h
    hdw_bit_number = R.L1Topo.FIRED
    ovf_bit_number = R.L1Topo.OVERFLOWN
    print 'Checking these algorithm bits:'
    print('sim_bit_number ',sim_bit_number)
    print('hdw_bit_number ',hdw_bit_number)
    print('ovf_bit_number ',ovf_bit_number)

    for algorithm_name in algorithm_names:
        branch_name = algorithm_name.replace('-', '_') + '_0_aBits'
        elist_name = 'elist_'+algorithm_name
        selection = ("{bn:s}.TestBitNumber({sbn:d})"
                     " != "
                     "{bn:s}.TestBitNumber({hbn:d})").format(**{'bn':branch_name,
                                                                'sbn':sim_bit_number,
                                                                'hbn':hdw_bit_number})
        selection_skip_overflow = "{bn:s}.TestBitNumber({obn:d})==0".format(**{'bn':branch_name,
                                                                               'obn':ovf_bit_number})
        selection = (selection+' && '+selection_skip_overflow) if options.skip_overflow else selection

        print 'selection: ',selection
        chain.Draw('>> '+elist_name, selection, 'entrylist')
        event_list = R.gDirectory.Get(elist_name)
        chain.SetEntryList(event_list)
        print 20*'-'
        print branch_name+' : got entry list with ',event_list.GetN(),' entries'
        print 20*'-'
        tprof = tprof_mismatches[algorithm_name]
        h_mis = hist_mismatches[algorithm_name]
        for iEntry in xrange(event_list.GetN()):
            chain.GetEntry(event_list.GetEntry(iEntry))
            aBits = getattr(chain, branch_name)
            jetTobs = [JetTob(w) for w in chain.jetTobs]
            n_jets = float(len(jetTobs))
            tprof.Fill(n_jets, 1.0)
            h_mis.Fill(n_jets)
            print ("runNumber = {0:d}  eventNumber = {1:d}".format(chain.runNumber, chain.eventNumber))
            print 'bit ',sim_bit_number,' (SIM) ', aBits.TestBitNumber(sim_bit_number),' bit ',hdw_bit_number,' (HDW)', aBits.TestBitNumber(hdw_bit_number)
            print 'jetTobs[%d]' % len(jetTobs)
            for i, jt in enumerate(jetTobs):
                print "[%d] (%f, %f, %f)"%(i, jt.energy8x8(), jt.eta(), jt.phi())

    out_filename = ('miss_vs_njets_no_ovf.root' if options.skip_overflow else
                    'miss_vs_njets.root')
    out_file = R.TFile(out_filename, 'recreate')
    out_file.cd()
    # for t in tprof_mismatches.values():
    #     t.SetDirectory(out_file)
    for h in hist_mismatches.values():
        h.SetDirectory(out_file)
    out_file.Write()
    out_file.Close()
 
def number_of_entries_to_process(available_entries, options=None):
    N = available_entries
    n = options.num_events
    s = options.skip_events
    to_process = (min([N, n, N-s]) if n and s else
                  min([N, n]) if n else
                  N-s if s else
                  N)
    to_process = to_process if to_process > 0 else 0
    return to_process

if __name__=='__main__':
    main()
