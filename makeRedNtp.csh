#!/bin/tcsh

set queue = cmslong

# submit the job only if the 2nd argument is 1
set run = 0
if ($#argv > 0) then
  set run = $1
  echo "run : $run "
endif 

# app to run
set app = ./tmp/redntpApp

if(! -e $app ) then
  echo "missing executable $app"
  exit 0
endif 

#set samples = ( VBF_HToGG_7TeV_powheg_M120_23nov DiPhotonBox_Pt25to250_TrackingParticles_7TeV-pythia6_bis DiPhotonJets_7TeV-madgraph )
#set samples = ( GJet_Pt-20_doubleEMEnriched_TuneZ2 )


#set bkgsamples = ( DiPhotonBox_Pt25to250_TrackingParticles_7TeV-pythia6_bis )

set higgssamples = ( VBF_HToGG_7TeV_powheg_M100_22nov VBF_HToGG_7TeV_powheg_M110_24nov_bis VBF_HToGG_7TeV_powheg_M115_22nov VBF_HToGG_7TeV_powheg_M120_23nov VBF_HToGG_7TeV_powheg_M130_22nov VBF_HToGG_7TeV_powheg_M140_22nov WH_ZH_TTH_HToGG_M-120_7TeV-pythia6_Fall10_8dec )


set qcdsamples = ( )
foreach i (`seq -w 0 58 `)
  set qcdsamples = ( $qcdsamples  QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6_a+b_${i} )
  #echo $qcdsamples[$i]
end

set gjetsamples = ( )
foreach i (`seq -w 0 10 `)
  set gjetsamples = ( $gjetsamples  GJet_Pt-20_doubleEMEnriched_TuneZ2_${i} )
end

set boxsamples = ( )
foreach i (`seq -w 0 10 `)
  set boxsamples = ( $boxsamples  DiPhotonBox_Pt25to250_TrackingParticles_7TeV-pythia6_bis_${i} )
end

set digammajetsamples = ( )
foreach i (`seq -w 0 10 `)
  set digammajetsamples = ( $digammajetsamples  DiPhotonJets_7TeV-madgraph_${i} )
end

set datasamples = ( )
foreach i (`seq -w 0 15 `)
  set datasamples = ( $datasamples  data_Cert_132440-149442_7TeV_StreamExpress_Collisions10_JSON_v3_NOVRERECO_2010B_smaller_new2_${i} )
end

set samples = ( $higgssamples $qcdsamples $gjetsamples  $boxsamples  $digammajetsamples  $datasamples )
#set samples = ( $datasamples )


set outdir = redntp.pt1stgam50.looseid.V5
set logdir = ${outdir}/log

mkdir -p $outdir  $logdir


set j1 = 20
set j2 = 15
set deta = 2.5
set zepp = 2.5
set mjj = 300.0

foreach sample ( $samples )
   set rootfile = "${outdir}/redntp_${sample}_${j1}_${j2}_${deta}_${zepp}_${mjj}.root"
   set jobname = "${sample}_${j1}_${j2}_${deta}_${zepp}_${mjj}"
   set logfile = "${logdir}/${sample}_${j1}_${j2}_${deta}_${zepp}_${mjj}.txt"
   set listfile = "list/${sample}.list"
   if(! -e $listfile ) then
      echo "skipping non-existent file $listfile"
      continue
   endif
   set command = "bsub -q ${queue} -o $logfile -J ${jobname} ${app} list/${sample}.list ${rootfile}"
   echo "submitted job: ${jobname}"
   #echo "   command: ${command}"
   if($run == 1) then
     ${command}
     #sleep 1
   endif
end
