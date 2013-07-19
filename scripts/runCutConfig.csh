#!/bin/tcsh

set queue = cmslong

# write cuts file only if the 1st argument is 1
set writefile = 0
if ($#argv > 0) then
  set writefile = $1
  echo "write : $writefile "
endif 

# submit the job only if the 2nd argument is 1
set run = 0
if ($#argv > 1) then
  set run = $2
  echo "run : $run "
endif 

# app to run
set app = ./tmp/analysis
if ($#argv > 2) then
  set app = $3
  echo "app : $app"
endif 

#set samples = ( VBF_HToGG_7TeV_powheg_M120_23nov DiPhotonBox_Pt25to250_TrackingParticles_7TeV-pythia6 DiPhotonJets_7TeV-madgraph GJet_Pt-20_doubleEMEnriched_TuneZ2 )
#set higgssamples = ( VBF_HToGG_7TeV_powheg_M120_23nov  )
#set samples = ( DiPhotonBox_Pt25to250_TrackingParticles_7TeV-pythia6_bis )

#set samples = ( VBF_HToGG_7TeV_powheg_M120_23nov DiPhotonBox_Pt25to250_TrackingParticles_7TeV-pythia6_bis DiPhotonJets_7TeV-madgraph )
#set samples = ( GJet_Pt-20_doubleEMEnriched_TuneZ2 )
#set samples = ( DiPhotonJets_7TeV-madgraph )

#set samples = ( VBF_HToGG_7TeV_powheg_M100_22nov VBF_HToGG_7TeV_powheg_M110_24nov_bis VBF_HToGG_7TeV_powheg_M115_22nov VBF_HToGG_7TeV_powheg_M120_23nov VBF_HToGG_7TeV_powheg_M130_22nov VBF_HToGG_7TeV_powheg_M140_22nov )

#set samples = ( QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6_1 QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6_2  )

set samples = ( QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6_a QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6_b )

set outdir = output.pt1stgam50.looseid
set logdir = log.pt1stgam50.looseid

mkdir -p $outdir  $logdir

set j2 = 15
foreach j1 ( 20 )
  foreach deta ( 2.5  )
    foreach zepp ( 2.5  )
      foreach mjj ( 300.0 )
        set cutfile = "cuts_${j1}_${j2}_${deta}_${zepp}_${mjj}.txt"
        if ( $writefile == 1 ) then
           echo $j1 >> $cutfile
           echo $j2 >> $cutfile
           echo $deta >> $cutfile
           echo $zepp >> $cutfile
           echo $mjj >> $cutfile
           echo "cutfile created: $cutfile"
        endif
        foreach sample ( $samples )
           set rootfile = "${outdir}/out_${sample}_${j1}_${j2}_${deta}_${zepp}_${mjj}.root"
           set jobname = "${sample}_${j1}_${j2}_${deta}_${zepp}_${mjj}"
           set logfile = "${logdir}/${sample}_${j1}_${j2}_${deta}_${zepp}_${mjj}.txt"
           set command = "bsub -q ${queue} -o $logfile -J ${jobname} ${app} ${sample}.list ${rootfile} 2 ${cutfile}"
           echo "submitting job ${jobname}"
           if($run == 1) then
             ${command}
             sleep 1
           endif
        end
      end
    end
  end
end 
