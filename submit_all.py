#!/bin/env python

import os

def main():
    # generated with
    # eos ls /eos/atlas/atlastier0/rucio/data17_13TeV/physics_L1Topo/00327764/data17_13TeV.00327764.physics_L1Topo.merge.RAW
    base_dir = '/afs/cern.ch/user/m/mpraderi/public/L1TopoStudies/Version0'
    input_filelist = base_dir+'/MarcoL1Topo/files.txt'

    for iLine, line in enumerate(open(input_filelist).readlines()):
        input_file = 'file'+str(iLine)+'.txt'
        aux_file = open(input_file, 'w')
        aux_file.write(line)
        aux_file.close()
        job_script_path = write_job(iJob = iLine,
                                    input_file=input_file,
                                    base_dir = base_dir)

        # cmd = "bsub -L /bin/bash -q 1nd -o job_log_%d.oe %s" % (iJob, job_script_path)
        cmd = "bsub -L /bin/bash -q 1nh -o job_log_%d.oe %s" % (iLine, job_script_path)
        print cmd
    

def write_job(iJob=0, input_file='', base_dir=''):
    "write the script and chmod it"
    parameters = {'iJob' : iJob,
                  'input_file' : input_file,
                  'base_dir' : base_dir
                  }
    script_file_path = base_dir+'/MarcoL1Topo/'+"script%d.sh" % iJob
    script_file = open(script_file_path, 'w')
    script_file.write("""#!/bin/bash

tmp_dir=${{TMPDIR}}/LSF_${{LSB_JOBID}}
mkdir -vp ${{tmp_dir}}
cd ${{tmp_dir}}
echo "working from `pwd`"

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${{ATLAS_LOCAL_ROOT_BASE}}/user/atlasLocalSetup.sh --quiet

cd {base_dir:s}
setupATLAS --quiet
rcSetup Base,2.4.28
rc find_packages

cd ${{tmp_dir}}

out_filename="log{iJob}.txt"

python {base_dir:s}/MarcoL1Topo/validate_lfv_triggers.py -v {base_dir:s}/MarcoL1Topo/{input_file} >> ${{out_filename}}


dest_dir="{base_dir:s}/MarcoL1Topo/Histograms/file{iJob}"

mkdir ${{dest_dir}}

rsync ${{out_filename}} ${{dest_dir}}/${{out_filename}}

rsync n_mu.png ${{dest_dir}}/n_mu.png
rsync n_mu.root ${{dest_dir}}/n_mu.root
rsync n_mu6ab.png ${{dest_dir}}/n_mu6ab.png
rsync n_mu6ab.root ${{dest_dir}}/n_mu6ab.root
rsync n_pairs_mu6_0dr15.png ${{dest_dir}}/n_pairs_mu6_0dr15.png
rsync n_pairs_mu6_0dr15.root ${{dest_dir}}/n_pairs_mu6_0dr15.root
rsync n_pairs_0dr15.png ${{dest_dir}}/n_pairs_0dr15.png
rsync n_pairs_0dr15.root ${{dest_dir}}/n_pairs_0dr15.root
rsync n_pairs_mu6ab.png ${{dest_dir}}/n_pairs_mu6ab.png
rsync n_pairs_mu6ab.root ${{dest_dir}}/n_pairs_mu6ab.root
rsync n_cand_pairs.png ${{dest_dir}}/n_cand_pairs.png
rsync n_cand_pairs.root ${{dest_dir}}/n_cand_pairs.root
rsync dr_min.png ${{dest_dir}}/dr_min.png
rsync dr_min.root ${{dest_dir}}/dr_min.root
rsync dr_0dr15.png ${{dest_dir}}/dr_0dr15.png
rsync dr_0dr15.root ${{dest_dir}}/dr_0dr15.root
rsync dr_mu6.png ${{dest_dir}}/dr_mu6.png
rsync dr_mu6.root ${{dest_dir}}/dr_mu6.root
rsync dr_mu6_0dr15.png ${{dest_dir}}/dr_mu6_0dr15.png
rsync dr_mu6_0dr15.root ${{dest_dir}}/dr_mu6_0dr15.root
rsync dr_any.png ${{dest_dir}}/dr_any.png
rsync dr_any.root ${{dest_dir}}/dr_any.root
rsync Phi_mu6.png ${{dest_dir}}/Phi_mu6.png
rsync Phi_mu6.root ${{dest_dir}}/Phi_mu6.root
rsync pt_0dr15.png ${{dest_dir}}/pt_0dr15.png
rsync pt_0dr15.root ${{dest_dir}}/pt_0dr15.root
rsync pt_any.png ${{dest_dir}}/pt_any.png
rsync pt_any.root ${{dest_dir}}/pt_any.root
rsync PhiEta_mu6.png ${{dest_dir}}/PhiEta_mu6.png
rsync PhiEta_mu6.root ${{dest_dir}}/PhiEta_mu6.root

""".format(**parameters))
    script_file.close()
    os.chmod(script_file_path, 0755)
    return script_file_path

if __name__=='__main__':
    main()