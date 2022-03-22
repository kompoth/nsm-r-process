#include <fstream>
#include <string>
#include <vector>
#include <valarray>

#include "skynet/BuildInfo.hpp"
#include "skynet/EquationsOfState/HelmholtzEOS.hpp"
#include "skynet/EquationsOfState/SkyNetScreening.hpp"
#include "skynet/Network/ReactionNetwork.hpp"
#include "skynet/Network/NSE.hpp"
#include "skynet/Network/NSEOptions.hpp"
#include "skynet/Reactions/REACLIBReactionLibrary.hpp"
#include "skynet/DensityProfiles/ConstantFunction.hpp"
#include "skynet/DensityProfiles/ExpTMinus3.hpp"
#include "skynet/Utilities/Interpolators/PiecewiseLinearFunction.hpp"


int
main(int argc, char *argv[]) {
  // Nuclei library
  auto nuclib = NuclideLibrary::CreateFromWebnucleoXML(SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml");
  
  // Options for reaction libraries
  NetworkOptions opts;
  opts.IsSelfHeating = true;
  opts.EnableScreening = true;
   
  // Reaction libraries (read user input)
  std::string user_strong = SkyNetRoot + "data/reaclib";
  std::string user_weak = SkyNetRoot + "data/reaclib";
  char ch;
  while ((ch = getopt(argc, argv, "s:w:")) != -1) {
    switch (ch) {
      case 's': user_strong = optarg; break;
      case 'w': user_weak = optarg; break;
    }
  }
  std::cout << "Using strong reaction lib: " << user_strong << std::endl;
  std::cout << "Using weak reaction lib: " << user_weak << std::endl;

  SkyNetScreening screen(nuclib);
  HelmholtzEOS helm(SkyNetRoot + "/data/helm_table.dat");
  REACLIBReactionLibrary strongReactionLibrary(user_strong, ReactionType::Strong, false, 
    LeptonMode::TreatAllAsDecayExceptLabelEC, "Strong reactions", nuclib, opts, true);
  REACLIBReactionLibrary symmetricFission(SkyNetRoot + "data/netsu_panov_symmetric_0neut", ReactionType::Strong, false,
    LeptonMode::TreatAllAsDecayExceptLabelEC, "Symmetric neutron induced fission with 0 free neutrons", nuclib, opts, false);
  REACLIBReactionLibrary spontaneousFission(SkyNetRoot + "data/netsu_sfis_Roberts2010rates", ReactionType::Strong, false,
    LeptonMode::TreatAllAsDecayExceptLabelEC, "Spontaneous fission", nuclib, opts, false);
  REACLIBReactionLibrary weakReactionLibrary(user_weak, ReactionType::Weak, false, 
    LeptonMode::TreatAllAsDecayExceptLabelEC, "Weak reactions", nuclib, opts, true);

  // Reaction rates output
  strongReactionLibrary.SetDoOutput(false);
  symmetricFission.SetDoOutput(false);
  spontaneousFission.SetDoOutput(false);
  weakReactionLibrary.SetDoOutput(false);
  
  // Reaction network
  ReactionLibs reactionLibraries{&strongReactionLibrary, &symmetricFission, &spontaneousFission, &weakReactionLibrary};
  ReactionNetwork net(nuclib, reactionLibraries, &helm, &screen, opts);
  
  double Ye = 0.1;
  double tau = 0.01;
  double s = 10;
  double T0 = 6.0; 

  // Run NSE with the temperature and entropy to find the initial density
  NSE nse(net.GetNuclideLibrary(), &helm, &screen);
  auto nseResult = nse.CalcFromTemperatureAndEntropy(T0, s, Ye);
  auto rhofunc = ExpTMinus3(nseResult.Rho(), tau);

  // Run simulation
  auto output = net.EvolveSelfHeatingWithInitialTemperature(
    nseResult.Y(), 0.0, 1.0, T0, &rhofunc, "history");

  // Save results
  std::vector<double> finalYVsA = output.FinalYVsA();
  FILE * f = fopen("final_abundance.txt", "w");
  for (unsigned int A = 0; A < finalYVsA.size(); ++A)
    fprintf(f, "%6i %30.20E\n", A, finalYVsA[A]);
}
