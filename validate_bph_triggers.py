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
    lvl1item_name  = 'L1_BPH-2M9-MU6MU4_BPH-0DR15-MU6MU4'           #ASK DAVIDE FOR NAME
    algorithm_name = '0DR15-MU4MU6'

    l = lvl1item_name
    num_binning   = (9 , -0.5, 8.5)
    dr_binning    = (30 , -0.1 , 5.9)
    binning_2im9 = (15 , 0.0 , 1.5)
    pt_binning    = (8, 3500.0 , 11500.0) 
    angle_binning = (28, -3.5, 3.5)
    for k in possible_outcomes: #initialize the histograms, they will still be empty after 
        histos[k] = {
            'n_mu'                     : R.TH1F('n_mu'+'_'+k                     , l+'; N input l1mus'                   , *num_binning),
            'n_mu6mu4'                 : R.TH1F('n_mu6mu4'+'_'+k                 , l+'; N mu6mu4 muons'                  , *num_binning),
            'n_pairs_mu6mu4_2m9_0dr15' : R.TH1F('n_pairs_mu6mu4_2m9_0dr15'+'_'+k , l+'; N mu6mu4_2m9_0dr15 pairs'        , *num_binning),
            'n_pairs_mu6mu4'           : R.TH1F('n_pairs_mu6mu4'+'_'+k           , l+'; N mu6mu4 pairs'                  , *num_binning),
            'n_cand_pairs'             : R.TH1F('n_cand_pairs'+'_'+k             , l+'; N candidate pairs'               , *num_binning),
            'dr_min_mu6mu4'            : R.TH1F('dr_min_mu6mu4'+'_'+k            , l+'; min #DeltaR best candidate pair' , *dr_binning),
            'dr_mu6mu4'                : R.TH1F('dr_mu6mu4'+'_'+k                , l+'; #DeltaR mu6mu4 pairs'            , *dr_binning),
            'dr_mu6mu4_2m9_0dr15'      : R.TH1F('dr_mu6mu4_2m9_0dr15'+'_'+k      , l+'; #DeltaR mu6mu4_2m9_0dr15'        , *binning_2im9),
            'Phi_mu6mu4'               : R.TH1F('Phi_mu6mu4'+'_'+k               , l+'; Phi angle any mu6mu4 muon'       , *angle_binning),
            'Eta_mu6mu4'               : R.TH1F('Eta_mu6mu4'+'_'+k               , l+'; Eta angle any mu6mu4 muon'       , *angle_binning),
            'pt_any'                   : R.TH1F('pt_any'+'_'+k                   , l+'; #Pt any muon'                    , *pt_binning),
            }
        histos2[k] = R.TH2F('PhiEta_mu6mu4'+'_'+k   , l+'; Phi angle any mu6mu4; Eta angle any mu6mu4' , *2*angle_binning)

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
        muons = [Muon(tob.pt, tob.eta, tob.phi) for tob in event.hdwMuonTOB
                 if tob.bcn==0] # only pick the ones from bunch crossing number 0

        list_mu4 = algo_MU4(muons) #mu4 list
        list_mu6mu4_2m9_0dr15_pairs, list_mu6mu4_pairs = algo_2M9_0DR15(list_mu4) #2im9_0dr15 couplelist
        list_2m9_0dr15_mu6mu4 = algo_MU6MU4_pairs(list_mu6mu4_2m9_0dr15_pairs)

        pass_emul = len(list_mu6mu4_2m9_0dr15_pairs)   #returns true if mu6mu4, 2m9 and 0dr15

        outcome = ('pass_em_pass_hw' if pass_hw and pass_emul else
                   'pass_em_fail_hw' if pass_emul else
                   'fail_em_pass_hw' if pass_hw else
                   'fail_em_fail_hw')
        valid_counters[outcome] += 1
        fill_histos(histos[outcome], histos2[outcome], muons, list_mu4,
                    list_mu6mu4_2m9_0dr15_pairs, list_mu6mu4_pairs) #fill histograms
        
        if debug and pass_hw:
            print "passed, %d muons" % len(muons)
        iEntry += 1


    print 'algo_counters:'
    pprint(dict(algo_counters))
    print 'item_counters:'
    pprint(dict(item_counters))
    print 'valid counters:'
    pprint(dict(valid_counters))

    #print errors
    p_p=valid_counters['pass_em_pass_hw']
    p_f=valid_counters['pass_em_fail_hw']
    f_p=valid_counters['fail_em_pass_hw']
    f_f=valid_counters['fail_em_fail_hw']

    total_imputs = p_p+p_f+f_p+f_f
    total_pass_em = p_p+p_f
    total_pass_hw = f_p+p_p
    total_discordance = 100.*(f_p+p_f)/total_imputs
    pass_em_discordance = 100.*p_f/total_pass_em
    pass_hw_discordance = 100.*f_p/total_pass_hw
    print('  total   error {:.2f}%'.format(total_discordance))
    print('  em pass error {:.2f}%'.format(pass_em_discordance))
    print('  hw pass error {:.2f}%'.format(pass_hw_discordance))

    c = R.TCanvas('c')
    order = [2,4,3,1]
    
    for name in histo_names:
        i = 0
        c.Clear()
        c.Divide(2,2)
        for outcome, hs in histos.items():
            h = histos[outcome][name]
            c.cd(order[i])
            h.Draw('h text')
            c.Update()
            i+=1
        
        c.SaveAs(name+'.png')
        c.SaveAs(name+'.root')
    
    i=0
    c.Clear()
    c.Divide(2,2)
    for outcome, h in histos2.items(): 
        c.cd(order[i])
        h.Draw('Colz')
        c.Update()
        i+=1
        if verbose:
            h.Print("all")
    if verbose:
        print('\n')

    c.SaveAs('PhiEta_mu6mu4.png')
    c.SaveAs('PhiEta_mu6mu4.root')

def algo_MU6MU4_pairs(list_2m9_0dr15_pairs):
    list_2m9_0dr15_mu6mu4 = []
    for couple in list_2m9_0dr15_pairs:
        if couple[0].p4.Pt()<8000 and couple[1].p4.Pt()<8000:
            list_2m9_0dr15_mu6mu4.append(couple)

    return list_2m9_0dr15_mu6mu4


def algo_2M9_0DR15(list_mu4): #retuns ordered list with couples of mu6mu4 muons satisfying 2M9_0DR15
    couples_any   = []
    list_mu6 = [muon for muon in list_mu4 if muon[1]>5500]
    n_mu4 = len(list_mu4)
    n_mu6 = len(list_mu6)
    n_dif = n_mu4-n_mu6
    for i in range(n_dif): #check all mu4 mu6 couples to see if 2M9_0DR15
        for muon, pt in list_mu6:
            dr = list_mu4[i].p4.DeltaR(muon.p4)
            m = (list_mu4[i].p4+muon.p4).M()
            couples_any.append((list_mu4[i][0], muon, dr, m))

    for i in range(n_mu6-1): #check all mu6 mu6 couples to see if 2M9_0DR15
        for j in range(i+1,n_mu6):
            dr = list_mu6[i].p4.DeltaR(list_mu6[j].p4)
            m = (list_mu6[i].p4+list_mu6[j].p4).M()
            couples_any.append((list_mu6[i][0], list_mu6[j][0], dr, m))



    couples_any.sort(key = lambda couple: couple[2]) #sort list
    couples_2m9_0dr15 = [couple for couple in couples_any if couple[2]<1.505 and (couple[3]>1999 and couple[3]<9001)] #take only 2mu9_0dr15
        return (couples_2m9_0dr15, couples_any, list_2m9_0dr15)

    return (couples_2m9_0dr15, couples_any)



def algo_MU4(muons): #returns sorted list of muons satistfying MU4
    mu6_list = []
    mu4_list = []

    for muon in muons:
        pt = muon.p4.Pt()
        #some error when doing >=4000 solved with >3500
        if pt>3500:                               
            mu4_list.append((muon,pt))

    mu4_list.sort(key = lambda muon: muon[1])
    
    return mu4_list

    
def fill_histos(histos, histos2, muons, list_mu4,
                list_mu6mu4_2m9_0dr15_pairs, list_mu6mu4_pairs): #fills histograms
    n_mu = len(muons)
    n_mu4= len(list_mu4)
    n_mu6mu4_pairs = len(list_mu6mu4_pairs)
    n_mu6mu4_2m9_0dr15_pairs = len(list_mu6mu4_2m9_0dr15_pairs)
    histos['n_mu'                    ].Fill(n_mu)
    histos['n_mu6mu4'                ].Fill(n_mu4)
    histos['n_pairs_mu6mu4_2m9_0dr15'].Fill(n_mu6mu4_2m9_0dr15_pairs)  #number of mu6mu4_2m9_0dr15 pairs
    # same number obtained here than using the formula n*(n-1)/2
    histos['n_pairs_mu6mu4'          ].Fill(n_mu6mu4_pairs)            #number of mu6mu4 pairs
    histos['n_cand_pairs'            ].Fill(n_pairs)                   #number of candidate pairs

    
    if n_mu6mu4_2m9_0dr15_pairs:
        for couple in list_mu6mu4_2m9_0dr15_pairs:
            histos['dr_mu6mu4_2m9_0dr15'].Fill(couple[2])
    
    if n_mu6mu4_pairs: 
        histos['dr_min_mu6mu4'].Fill(list_pairs[0][2])
        for couple in list_mu6mu4_pairs:
            histos['dr_mu6mu4'].Fill(couple[2])

   
    for muon in muons: #fill histogram of momentums
        histos['pt_any'].Fill(muon.p4.Pt())
  
    for muon, pt in list_mu4: #fill histograms of angles
        Phi = muon.p4.Phi()
        Eta = muon.p4.Eta()
        histos['Phi_mu6mu4'].Fill(Phi)
        histos['Eta_mu6mu4'].Fill(Eta)
        histos2.Fill(Phi, Eta)


class Muon(object):
    def __init__(self, pt, eta, phi):
        muon_mass = 105.65
        #tlv = R.TLorentzVector() # four-momentum
        #self.p4 = tlv.SetPtEtaPhiE(pt, eta, phi, energy)
        self.p4 = R.TLorentzVector() # four-momentum
        self.p4.SetPtEtaPhiM(pt, eta, phi, muon_mass)

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
