#!/bin/env python

# Example script to validate LFV topo triggers
#
# Example input file:
# /afs/cern.ch/user/g/gerbaudo/public/tmp/for_marco/user.olya.11640704._000003.tmptrig.root
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
    # # test: print all branches available in the tree
    # print 'chain ',chain
    # print 'branches:'
    # print 'list of branches ',chain.GetListOfBranches()
    # print '\n'.join([k.GetName() for k in chain.GetListOfBranches()])
    # print
    # return

    iEntry = 0
    possible_outcomes = ['pass_em_pass_hw', 'pass_em_fail_hw', 'fail_em_pass_hw', 'fail_em_fail_hw']
    algo_counters = defaultdict(int)
    item_counters = defaultdict(int)
    valid_counters = {k:0 for k in ['overflow'] + possible_outcomes}

    histos = {}
    histos2= {}
    lvl1item_name  = 'L1_LFV-MU6'
    algorithm_name = '0DR15-2MU6ab'

    l = lvl1item_name
    num_binning   = (9 , -0.5, 8.5)
    dr_binning    = (50 , 0.0 , 1.5)
    pt_binning    = (10, 3000.0 , 11000.0) 
    angle_binning = (50, -3.5, 3.5)
    for k in possible_outcomes: #initialize the histograms, they will still be empty after 
        histos[k] = {
            'n_mu'          : R.TH1F('n_mu'+'_'+k          , l+'; N input l1mus'                   , *num_binning),
            'n_mu6ab'       : R.TH1F('n_mu6ab'+'_'+k       , l+'; N mu6 muons'                     , *num_binning),
            'n_pairs_0dr15' : R.TH1F('n_pairs_0dr15'+'_'+k , l+'; N 0dr15 pairs'                   , *num_binning),
            'n_pairs_mu6ab' : R.TH1F('n_pairs_mu6ab'+'_'+k , l+'; N mu6 pairs'                     , *num_binning),
            'n_cand_pairs'  : R.TH1F('n_cand_pairs'+'_'+k  , l+'; N candidate pairs'               , *num_binning),
            'dr_min'        : R.TH1F('dr_min'+'_'+k        , l+'; min #DeltaR best candidate pair' , *dr_binning),
            'dr_any'        : R.TH1F('dr_any'+'_'+k        , l+'; #DeltaR any candidate pair'      , *dr_binning),
            'Phi_mu6'       : R.TH1F('Phi_mu6'+'_'+k       , l+'; Phi angle any mu6 muon'          , *angle_binning),
            'pt_any'        : R.TH1F('pt_any'+'_'+k        , l+'; #Pt any muon'                    , *pt_binning),
            }
        histos2[k] = R.TH2F('PhiEta_mu6'+'_'+k   , l+'; Phi angle any mu6; Eta angle any mu6' , *2*angle_binning)

    histo_names = [name for name, histo in histos[possible_outcomes[0]].items()]

    for iEvent, event in enumerate(chain):
        if num_skip and iEvent<num_skip: continue
        if iEntry > num_toprocess: break
        # see how branches are created in TopoNtuple.py
        item_bits = item_tbits2bits(getattr(event,
                                            lvl1item_name.replace('-','_')+'_iBits'))
        increment_counters(item_counters, item_bits)
        # # These are not filled, so don't bother for now # DG-2016-06-23
        # algo_bits = algo_tbits2bits(getattr(event, algorithm_name+'_0_aBits'))
        # increment_counters(algo_counters, algo_bits)
        pass_hw = item_bits['TDT_TBP']
        pass_sim = item_bits['L1SIMULATED']
        
        # overflown = algo_bits['OVERFLOWN']
        # if overflown:
        #     valid_counters['overflow'] += 1
        #     continue
        # emTobs = [EmTob(w) for w in event.emTobs]
        # jetTobs = [JetTob(w) for w in event.jetTobs]
        # tauTobs = [TauTob(w) for w in event.tauTobs]
        # if debug:
        #     print 'emTobs[%d]' % len(emTobs)
        #     for i, et in enumerate(emTobs):
        #         print "[%d] (%f, %f)"%(i, et.eta(), et.phi())
        #     print 'jetTobs[%d]' % len(jetTobs)
        #     for i, jt in enumerate(jetTobs):
        #         print "[%d] (%f, %f)"%(i, jt.eta(), jt.phi())
        #     print 'tauTobs[%d]' % len(tauTobs)
        #     for i, tt in enumerate(tauTobs):
        #         print "[%d] (%f, %f)"%(i, tt.eta(), tt.phi())

        # these are EnhancedMuonTOB objects
        muons = [Muon(tob.pt, tob.eta, tob.phi, tob.pt) for tob in event.hdwMuonTOB
                 if tob.bcn==0] # only pick the ones from bunch crossing number 0
        n_mu = len(muons) #number of events
        list_mu6ab = algo_MU6ab(muons) #mu6 list
        n_mu6ab = len(list_mu6ab)
        list_0dr15 = algo_0DR15(muons, n_mu) #0dr15 couplelist
        # list_0DR15 = algo_0DR15(list_mu6ab, n_mu) #0dr15 couplelist only for mu6ab this should be more correct

        #TODO 0dr2pi

        n_0dr15 = len(list_0dr15)
        pass_emul = n_mu6ab-2 and n_0dr15   #returns true if 2mu6ab and 0dr15
        
        outcome = ('pass_em_pass_hw' if pass_hw and pass_emul else
                   'pass_em_fail_hw' if pass_emul else
                   'fail_em_pass_hw' if pass_hw else
                   'fail_em_fail_hw')
        valid_counters[outcome] += 1
        fill_histos(histos[outcome], histos2[outcome], muons, list_mu6ab, list_0dr15) #fill histograms
        if debug and pass_hw:
            print "passed, %d muons" % len(muons)
        iEntry += 1


    print 'algo_counters:'
    pprint(dict(algo_counters))
    print 'item_counters:'
    pprint(dict(item_counters))
    print 'valid counters:'
    pprint(dict(valid_counters))

    c = R.TCanvas('c') #plot
    
    for name in histo_names:
        i = 1
        c.Clear()
        c.Divide(2,2)
        for outcome, hs in histos.items():
            h = histos[outcome][name]
            c.cd(i)
            h.Draw('h text')
            c.Update()
            i+=1
        c.SaveAs(name+'.png')
    
    i=1
    c.Clear()
    c.Divide(2,2)
    for outcome, h in histos2.items(): 
        c.cd(i)
        h.Draw('Colz')
        c.Update()
        i+=1
    c.SaveAs(name+'.png')

def algo_0DR15(muons, n_mu): #retuns ordered list with any couple of muons satisfying 0DR15
    couples_0dr15 = []
    for i in range(n_mu-1): #check all muon couples to see if Delta r is lower than 1.5
        for j in range(i+1,n_mu):
            dr = muons[i].p4.DeltaR(muons[j].p4)
            if dr<1.5:    #not necessary to impose dr>0
                couples_0dr15.append((muons[i],muons[j], dr))
    
    couples_0dr15.sort(key = lambda couple_0dr15: couple_0dr15[2]) #sort list
    return couples_0dr15


def algo_MU6ab(muons): #returns list with all muons satifying MU6 sorted by energy
    mu6ab_list = []

    for muon in muons:
        pt = muon.p4.Pt()
        if pt>6000:
            mu6ab_list.append((muon,pt))
    
    mu6ab_list.sort(key = lambda mu6ab: mu6ab[1]) #sort list
    return mu6ab_list 

    
def fill_histos(histos, histos2, muons, list_mu6ab, list_0dr15): #fills histograms
   n_mu = len(muons)
   n_mu6= len(list_mu6ab)
   n_0dr15 = len(list_0dr15)
   histos['n_mu'         ].Fill(n_mu)
   histos['n_mu6ab'      ].Fill(n_mu6)
   histos['n_pairs_0dr15'].Fill(n_mu6*(n_mu6-1.)/2) #number of mu6 pairs
   histos['n_pairs_mu6ab'].Fill(n_0dr15)            #number of 0dr15 pairs
   histos['n_cand_pairs' ].Fill(n_mu*(n_mu-1.)/2)   #number of candidate pairs
   
   if n_0dr15: #fill histograms of dr
       histos['dr_min'].Fill(list_0dr15[0][2])
       
       for couple in list_0dr15:
           histos['dr_any'].Fill(couple[2])
   
   for muon, pt in list_mu6ab: #fill histograms of angles
       Phi = muon.p4.Phi()
       histos['Phi_mu6'].Fill(Phi)
       histos2.Fill(Phi, muon.p4.Eta())
   
   for muon in muons: #fill histogram of momentums
       histos['pt_any'].Fill(muon.p4.Pt())


class Muon(object):
    def __init__(self, pt, eta, phi, energy):
         #tlv = R.TLorentzVector() # four-momentum
         #self.p4 = tlv.SetPtEtaPhiE(pt, eta, phi, energy)
         self.p4 = R.TLorentzVector() # four-momentum
         self.p4.SetPtEtaPhiE(pt, eta, phi, energy)

def algo_bit_names_and_numbers():
    "Bits stored in the TBits for each L1 algorithm"
    return (('FIRED'        , R.L1Topo.FIRED       ),
            ('OVERFLOWN'    , R.L1Topo.OVERFLOWN   ),
            ('TOB_EMULATED' , R.L1Topo.TOB_EMULATED),
            ('ROI_EMULATED' , R.L1Topo.ROI_EMULATED),
            ('CTP_TIP'      , R.L1Topo.CTP_TIP     ),
            ('TOPOSIM'      , R.L1Topo.TOPOSIM     ))

def item_bit_names_and_numbers():
    "Bits stored in the TBits for each L1 item"
    return (('TDT_TBP'    , R.L1Topo.TDT_TBP    ),
            ('TDT_TAP'    , R.L1Topo.TDT_TAP    ),
            ('TDT_TAV'    , R.L1Topo.TDT_TAV    ),
            ('L1EMULATED' , R.L1Topo.L1EMULATED ),
            ('L1SIMULATED', R.L1Topo.L1SIMULATED),
            ('CTP_TBP'    , R.L1Topo.CTP_TBP    ))

def algo_tbits2bits(bits=None):
    "convert TBits to dict"
    return dict((k, bits.TestBitNumber(b)) for k, b in algo_bit_names_and_numbers())

def item_tbits2bits(bits=None):
    "convert TBits to dict"
    return dict((k, bits.TestBitNumber(b)) for k, b in item_bit_names_and_numbers())

def increment_counters(counters={}, bits={}):
    for k, b in bits.items():
        counters[k] += 1 if b else 0


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
