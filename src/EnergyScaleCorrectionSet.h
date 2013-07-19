#include "EnergyScaleCorrection.h"
#include "TString.h"

void fillCorrections(TString& escaleCorrections,EnergyScaleCorrection::energyScaleParameters& scaleCorrections)
{ 
  if (escaleCorrections.CompareTo("May10PromptV4") == 0 )
    {
      std::cout << "Initializing May10PromptV4 corrections" << std::endl;
      // 160404 - 163869 (203.5 Pb)
      // 165071 - 165970 (139.6 pb)
      // 165971 - 166502 (187.9 pb)
      // 166503 - 166861 (191.6 pb)
      // 166862 - 167784 (252.8 pb) 
      scaleCorrections.n_categories=4;
      scaleCorrections.categoryType="2CatR9_EBEE";
      scaleCorrections.parameterSetName="May10PromptV4";

      EnergyScaleOffset energyCorrection_MC(0,10);
      //EB
      energyCorrection_MC.scale_offset["EBHighR9"]=1;
      energyCorrection_MC.scale_offset_error["EBHighR9"]=0.;
      energyCorrection_MC.scale_offset["EBLowR9"]=1.0;
      energyCorrection_MC.scale_offset_error["EBLowR9"]=0.0;
      //EE	       MC
      energyCorrection_MC.scale_offset["EEHighR9"]=1.0;
      energyCorrection_MC.scale_offset_error["EEHighR9"]=0.0;
      energyCorrection_MC.scale_offset["EELowR9"]=1.0;
      energyCorrection_MC.scale_offset_error["EELowR9"]=0.;

      energyCorrection_MC.smearing["EBHighR9"]=0.0141;
      energyCorrection_MC.smearing_error["EBHighR9"]=0.0007;
      energyCorrection_MC.smearing["EBLowR9"]=0.0200;
      energyCorrection_MC.smearing_error["EBLowR9"]=0.0004;
      //EE	       MC
      energyCorrection_MC.smearing["EEHighR9"]=0.0474;
      energyCorrection_MC.smearing_error["EEHighR9"]=0.0014;
      energyCorrection_MC.smearing["EELowR9"]=0.0361;
      energyCorrection_MC.smearing_error["EELowR9"]=0.0014;

      scaleCorrections.scale_offset_byrun.push_back(energyCorrection_MC);      
      // RunRange 160404 - 163869
      EnergyScaleOffset energyCorrection_0(160404,163869);
      //EB
      energyCorrection_0.scale_offset["EBHighR9"]=1-0.0046;
      energyCorrection_0.scale_offset_error["EBHighR9"]=0.0006;
      energyCorrection_0.scale_offset["EBLowR9"]=1.0+0.0025;
      energyCorrection_0.scale_offset_error["EBLowR9"]=0.0004;
      //EE
      energyCorrection_0.scale_offset["EEHighR9"]=1.0+0.0050;
      energyCorrection_0.scale_offset_error["EEHighR9"]=0.0019;
      energyCorrection_0.scale_offset["EELowR9"]=1.0-0.0023;
      energyCorrection_0.scale_offset_error["EELowR9"]=0.0016;

      energyCorrection_0.smearing["EBHighR9"]=0.0141;
      energyCorrection_0.smearing_error["EBHighR9"]=0.0007;
      energyCorrection_0.smearing["EBLowR9"]=0.0200;
      energyCorrection_0.smearing_error["EBLowR9"]=0.0004;
      //EE	       0
      energyCorrection_0.smearing["EEHighR9"]=0.0474;
      energyCorrection_0.smearing_error["EEHighR9"]=0.0014;
      energyCorrection_0.smearing["EELowR9"]=0.0361;
      energyCorrection_0.smearing_error["EELowR9"]=0.0014;

      scaleCorrections.scale_offset_byrun.push_back(energyCorrection_0);

      // RunRange 165071 - 165970
      EnergyScaleOffset energyCorrection_1(165071,165970);
      //EB
      energyCorrection_1.scale_offset["EBHighR9"]=1.0-0.0007;
      energyCorrection_1.scale_offset_error["EBHighR9"]=0.0007;
      energyCorrection_1.scale_offset["EBLowR9"]=1.0+0.0049;
      energyCorrection_1.scale_offset_error["EBLowR9"]=0.0004;
      //EE	       1
      energyCorrection_1.scale_offset["EEHighR9"]=1.0+0.0249;
      energyCorrection_1.scale_offset_error["EEHighR9"]=0.0022;
      energyCorrection_1.scale_offset["EELowR9"]=1.0+0.0062;
      energyCorrection_1.scale_offset_error["EELowR9"]=0.0018;

      energyCorrection_1.smearing["EBHighR9"]=0.0141;
      energyCorrection_1.smearing_error["EBHighR9"]=0.0007;
      energyCorrection_1.smearing["EBLowR9"]=0.0200;
      energyCorrection_1.smearing_error["EBLowR9"]=0.0004;
      //EE	       1
      energyCorrection_1.smearing["EEHighR9"]=0.0474;
      energyCorrection_1.smearing_error["EEHighR9"]=0.0014;
      energyCorrection_1.smearing["EELowR9"]=0.0361;
      energyCorrection_1.smearing_error["EELowR9"]=0.0014;

      scaleCorrections.scale_offset_byrun.push_back(energyCorrection_1);

      // RunRange 165971 - 166502
      EnergyScaleOffset energyCorrection_2(165071,166502);
      //EB
      energyCorrection_2.scale_offset["EBHighR9"]=1.0+0.0003;
      energyCorrection_2.scale_offset_error["EBHighR9"]=0.0006;
      energyCorrection_2.scale_offset["EBLowR9"]=1.0+0.0067;
      energyCorrection_2.scale_offset_error["EBLowR9"]=0.0004;
      //EE	       2
      energyCorrection_2.scale_offset["EEHighR9"]=1.0+0.0376;
      energyCorrection_2.scale_offset_error["EEHighR9"]=0.0019;
      energyCorrection_2.scale_offset["EELowR9"]=1.0+0.0133;
      energyCorrection_2.scale_offset_error["EELowR9"]=0.0015;

      energyCorrection_2.smearing["EBHighR9"]=0.0141;
      energyCorrection_2.smearing_error["EBHighR9"]=0.0007;
      energyCorrection_2.smearing["EBLowR9"]=0.0200;
      energyCorrection_2.smearing_error["EBLowR9"]=0.0004;
      //EE	       2
      energyCorrection_2.smearing["EEHighR9"]=0.0474;
      energyCorrection_2.smearing_error["EEHighR9"]=0.0014;
      energyCorrection_2.smearing["EELowR9"]=0.0361;
      energyCorrection_2.smearing_error["EELowR9"]=0.0014;

      scaleCorrections.scale_offset_byrun.push_back(energyCorrection_2);

      // RunRange 166503 - 166861
      EnergyScaleOffset energyCorrection_3(166503,166861);
      //EB
      energyCorrection_3.scale_offset["EBHighR9"]=1.0+0.0011;
      energyCorrection_3.scale_offset_error["EBHighR9"]=0.0006;
      energyCorrection_3.scale_offset["EBLowR9"]=1.0+0.0063;
      energyCorrection_3.scale_offset_error["EBLowR9"]=0.0004;
      //EE	       3
      energyCorrection_3.scale_offset["EEHighR9"]=1.0+0.045;
      energyCorrection_3.scale_offset_error["EEHighR9"]=0.002;
      energyCorrection_3.scale_offset["EELowR9"]=1.0+0.0178;
      energyCorrection_3.scale_offset_error["EELowR9"]=0.0015;

      energyCorrection_3.smearing["EBHighR9"]=0.0141;
      energyCorrection_3.smearing_error["EBHighR9"]=0.0007;
      energyCorrection_3.smearing["EBLowR9"]=0.0200;
      energyCorrection_3.smearing_error["EBLowR9"]=0.0004;
      //EE	       3
      energyCorrection_3.smearing["EEHighR9"]=0.0474;
      energyCorrection_3.smearing_error["EEHighR9"]=0.0014;
      energyCorrection_3.smearing["EELowR9"]=0.0361;
      energyCorrection_3.smearing_error["EELowR9"]=0.0014;

      scaleCorrections.scale_offset_byrun.push_back(energyCorrection_3);

      // RunRange 166862 - 167913
      EnergyScaleOffset energyCorrection_4(166862,167913);
      //EB
      energyCorrection_4.scale_offset["EBHighR9"]=1.0+0.0014;
      energyCorrection_4.scale_offset_error["EBHighR9"]=0.0005;
      energyCorrection_4.scale_offset["EBLowR9"]=1.0+0.0074;
      energyCorrection_4.scale_offset_error["EBLowR9"]=0.0003;
      //EE	       4
      energyCorrection_4.scale_offset["EEHighR9"]=1.0+0.0561;
      energyCorrection_4.scale_offset_error["EEHighR9"]=0.0018;
      energyCorrection_4.scale_offset["EELowR9"]=1.0+0.0273;
      energyCorrection_4.scale_offset_error["EELowR9"]=0.0013;

      energyCorrection_4.smearing["EBHighR9"]=0.0141;
      energyCorrection_4.smearing_error["EBHighR9"]=0.0007;
      energyCorrection_4.smearing["EBLowR9"]=0.0200;
      energyCorrection_4.smearing_error["EBLowR9"]=0.0004;
      //EE	       4
      energyCorrection_4.smearing["EEHighR9"]=0.0474;
      energyCorrection_4.smearing_error["EEHighR9"]=0.0014;
      energyCorrection_4.smearing["EELowR9"]=0.0361;
      energyCorrection_4.smearing_error["EELowR9"]=0.0014;

      scaleCorrections.scale_offset_byrun.push_back(energyCorrection_4);
      
    }
  else if (escaleCorrections.CompareTo("Jul05") == 0 )
    {
      std::cout << "Initializing Jul05 corrections" << std::endl;
      // 160404 - 163869 (203.5 Pb)
      // 165071 - 165970 (139.6 pb)
      // 165971 - 166502 (187.9 pb)
      // 166503 - 166861 (191.6 pb)
      // 166862 - 167784 (252.8 pb) 
      scaleCorrections.n_categories=4;
      scaleCorrections.categoryType="2CatR9_EBEE";
      scaleCorrections.parameterSetName="Jul05";

      EnergyScaleOffset energyCorrection_MC(0,10);
      //EB
      energyCorrection_MC.scale_offset["EBHighR9"]=1;
      energyCorrection_MC.scale_offset_error["EBHighR9"]=0.;
      energyCorrection_MC.scale_offset["EBLowR9"]=1.0;
      energyCorrection_MC.scale_offset_error["EBLowR9"]=0.0;
      //EE	       MC
      energyCorrection_MC.scale_offset["EEHighR9"]=1.0;
      energyCorrection_MC.scale_offset_error["EEHighR9"]=0.0;
      energyCorrection_MC.scale_offset["EELowR9"]=1.0;
      energyCorrection_MC.scale_offset_error["EELowR9"]=0.;

      energyCorrection_MC.smearing["EBHighR9"]=0.0092;
      energyCorrection_MC.smearing_error["EBHighR9"]=0.0020;
      energyCorrection_MC.smearing["EBLowR9"]=0.0170;
      energyCorrection_MC.smearing_error["EBLowR9"]=0.0042;
      //EE	       MC
      energyCorrection_MC.smearing["EEHighR9"]=0.0292;
      energyCorrection_MC.smearing_error["EEHighR9"]=0.0050;
      energyCorrection_MC.smearing["EELowR9"]=0.0289;
      energyCorrection_MC.smearing_error["EELowR9"]=0.0042;

      scaleCorrections.scale_offset_byrun.push_back(energyCorrection_MC);      
      // RunRange 160404 - 163869
      EnergyScaleOffset energyCorrection_0(160404,167913);
      //EB
      energyCorrection_0.scale_offset["EBHighR9"]=1-0.0003;
      energyCorrection_0.scale_offset_error["EBHighR9"]=0.0006;
      energyCorrection_0.scale_offset["EBLowR9"]=1.0+0.0054;
      energyCorrection_0.scale_offset_error["EBLowR9"]=0.0034;
      //EE
      energyCorrection_0.scale_offset["EEHighR9"]=1.0-0.0003;
      energyCorrection_0.scale_offset_error["EEHighR9"]=0.0026;
      energyCorrection_0.scale_offset["EELowR9"]=1.0+0.0024;
      energyCorrection_0.scale_offset_error["EELowR9"]=0.0026;

      energyCorrection_0.smearing["EBHighR9"]=0.0092;
      energyCorrection_0.smearing_error["EBHighR9"]=0.0020;
      energyCorrection_0.smearing["EBLowR9"]=0.0170;
      energyCorrection_0.smearing_error["EBLowR9"]=0.0042;
      //EE	       0
      energyCorrection_0.smearing["EEHighR9"]=0.0292;
      energyCorrection_0.smearing_error["EEHighR9"]=0.0050;
      energyCorrection_0.smearing["EELowR9"]=0.0289;
      energyCorrection_0.smearing_error["EELowR9"]=0.0042;

      scaleCorrections.scale_offset_byrun.push_back(energyCorrection_0);
    }
  else if (escaleCorrections.CompareTo("LP11") == 0 )
    {
      std::cout << "Initializing LP11 corrections" << std::endl;

      scaleCorrections.n_categories=4;
      scaleCorrections.categoryType="2CatR9_EBEE";
      scaleCorrections.parameterSetName="LP11";

      EnergyScaleOffset energyCorrection_MC(0,10);
      //EB
      energyCorrection_MC.scale_offset["EBHighR9"]=1;
      energyCorrection_MC.scale_offset_error["EBHighR9"]=0.;
      energyCorrection_MC.scale_offset["EBLowR9"]=1.0;
      energyCorrection_MC.scale_offset_error["EBLowR9"]=0.0;
      //EE	       MC
      energyCorrection_MC.scale_offset["EEHighR9"]=1.0;
      energyCorrection_MC.scale_offset_error["EEHighR9"]=0.0;
      energyCorrection_MC.scale_offset["EELowR9"]=1.0;
      energyCorrection_MC.scale_offset_error["EELowR9"]=0.;

      energyCorrection_MC.smearing["EBHighR9"]=0.0100;
      energyCorrection_MC.smearing_error["EBHighR9"]=0.0020;
      energyCorrection_MC.smearing["EBLowR9"]=0.0170;
      energyCorrection_MC.smearing_error["EBLowR9"]=0.0042;
      //EE	       MC
      energyCorrection_MC.smearing["EEHighR9"]=0.0303;
      energyCorrection_MC.smearing_error["EEHighR9"]=0.0050;
      energyCorrection_MC.smearing["EELowR9"]=0.0296;
      energyCorrection_MC.smearing_error["EELowR9"]=0.0042;

      scaleCorrections.scale_offset_byrun.push_back(energyCorrection_MC);      
      // RunRange 160404 - 163869
      EnergyScaleOffset energyCorrection_0(160404,167913);
      //EB
      energyCorrection_0.scale_offset["EBHighR9"]=1-0.0003;
      energyCorrection_0.scale_offset_error["EBHighR9"]=0.0006;
      energyCorrection_0.scale_offset["EBLowR9"]=1.0+0.0054;
      energyCorrection_0.scale_offset_error["EBLowR9"]=0.0034;
      //EE
      energyCorrection_0.scale_offset["EEHighR9"]=1.0-0.0003;
      energyCorrection_0.scale_offset_error["EEHighR9"]=0.0026;
      energyCorrection_0.scale_offset["EELowR9"]=1.0+0.0024;
      energyCorrection_0.scale_offset_error["EELowR9"]=0.0026;

      energyCorrection_0.smearing["EBHighR9"]=0.0092;
      energyCorrection_0.smearing_error["EBHighR9"]=0.0020;
      energyCorrection_0.smearing["EBLowR9"]=0.0170;
      energyCorrection_0.smearing_error["EBLowR9"]=0.0042;
      //EE	       0
      energyCorrection_0.smearing["EEHighR9"]=0.0292;
      energyCorrection_0.smearing_error["EEHighR9"]=0.0050;
      energyCorrection_0.smearing["EELowR9"]=0.0289;
      energyCorrection_0.smearing_error["EELowR9"]=0.0042;

      scaleCorrections.scale_offset_byrun.push_back(energyCorrection_0);

      EnergyScaleOffset energyCorrection_1(167914,172619);
      //EB
      energyCorrection_1.scale_offset["EBHighR9"]=1.0+0.0020;
      energyCorrection_1.scale_offset_error["EBHighR9"]=0.0005;
      energyCorrection_1.scale_offset["EBLowR9"]=1.0+0.0081;
      energyCorrection_1.scale_offset_error["EBLowR9"]=0.0033;
      //EE	       1
      energyCorrection_1.scale_offset["EEHighR9"]=1.0-0.0109;
      energyCorrection_1.scale_offset_error["EEHighR9"]=0.0010;
      energyCorrection_1.scale_offset["EELowR9"]=1.0-0.0043;
      energyCorrection_1.scale_offset_error["EELowR9"]=0.0026;

      energyCorrection_1.smearing["EBHighR9"]=0.0092;
      energyCorrection_1.smearing_error["EBHighR9"]=0.0020;
      energyCorrection_1.smearing["EBLowR9"]=0.0170;
      energyCorrection_1.smearing_error["EBLowR9"]=0.0042;
      //EE	       1
      energyCorrection_1.smearing["EEHighR9"]=0.0292;
      energyCorrection_1.smearing_error["EEHighR9"]=0.0050;
      energyCorrection_1.smearing["EELowR9"]=0.0289;
      energyCorrection_1.smearing_error["EELowR9"]=0.0042;

      scaleCorrections.scale_offset_byrun.push_back(energyCorrection_1);

      EnergyScaleOffset energyCorrection_2(172620,172802);
      //EB
      energyCorrection_2.scale_offset["EBHighR9"]=1.0+0.0035;
      energyCorrection_2.scale_offset_error["EBHighR9"]=0.0005;
      energyCorrection_2.scale_offset["EBLowR9"]=1.0+0.0097;
      energyCorrection_2.scale_offset_error["EBLowR9"]=0.0033;
      //EE	       2
      energyCorrection_2.scale_offset["EEHighR9"]=1.0-0.0215;
      energyCorrection_2.scale_offset_error["EEHighR9"]=0.0023;
      energyCorrection_2.scale_offset["EELowR9"]=1.0-0.0158;
      energyCorrection_2.scale_offset_error["EELowR9"]=0.0030;

      energyCorrection_2.smearing["EBHighR9"]=0.0092;
      energyCorrection_2.smearing_error["EBHighR9"]=0.0020;
      energyCorrection_2.smearing["EBLowR9"]=0.0170;
      energyCorrection_2.smearing_error["EBLowR9"]=0.0042;
      //EE	       2
      energyCorrection_2.smearing["EEHighR9"]=0.0292;
      energyCorrection_2.smearing_error["EEHighR9"]=0.0050;
      energyCorrection_2.smearing["EELowR9"]=0.0289;
      energyCorrection_2.smearing_error["EELowR9"]=0.0042;

      scaleCorrections.scale_offset_byrun.push_back(energyCorrection_2);
    }
  else if (escaleCorrections.CompareTo("noCalib") == 0 )
    {
      std::cout << "Not applying any ad-hoc scale corrections" << std::endl;
    }
  else
    {
      std::cout << "Unknowm parameter set required. Not applying any ad-hoc scale corrections" << std::endl;
    }
}
