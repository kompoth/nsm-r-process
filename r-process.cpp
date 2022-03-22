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

#include "LinearExpansion.hpp"

#define LN10 2.30258509299

int
main(int argc, char *argv[]) {
  // Library of nuclei
  auto nuclib = NuclideLibrary::CreateFromWebnucleoXML(SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml");
  
  // Options for reaction libraries
  NetworkOptions opts;
  opts.ConvergenceCriterion = NetworkConvergenceCriterion::Mass;
  opts.MassDeviationThreshold = 1.0E-10;
  opts.IsSelfHeating = true;
  opts.EnableScreening = true;
   
  // Reaction libraries (read user input)
  std::string user_strong = SkyNetRoot + "data/reaclib";
  if (argc > 1) {
    user_strong = argv[1];
  }
  std::cout << "using strong reaction lib: " << user_strong << std::endl;
  SkyNetScreening screen(nuclib);
  HelmholtzEOS helm(SkyNetRoot + "/data/helm_table.dat");
  REACLIBReactionLibrary strongReactionLibrary(user_strong, ReactionType::Strong, false, 
    LeptonMode::TreatAllAsDecayExceptLabelEC, "Strong reactions", nuclib, opts, true);
  REACLIBReactionLibrary symmetricFission(SkyNetRoot + "data/netsu_panov_symmetric_0neut", ReactionType::Strong, false,
    LeptonMode::TreatAllAsDecayExceptLabelEC, "Symmetric neutron induced fission with 0 free neutrons", nuclib, opts, false);
  REACLIBReactionLibrary spontaneousFission(SkyNetRoot + "data/netsu_sfis_Roberts2010rates", ReactionType::Strong, false,
    LeptonMode::TreatAllAsDecayExceptLabelEC, "Spontaneous fission", nuclib, opts, false);
  REACLIBReactionLibrary weakReactionLibrary(SkyNetRoot + "data/reaclib", ReactionType::Weak, false, 
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
  //double rho0 = 9.58e+9;
  double tau = 0.01;
  double T0 = 6.; 

  // run NSE with the temperature and entropy to find the initial density
  NSEOptions nse_opts;
  nse_opts.MaxIterations = 1000;
  nse_opts.DoScreening = true;
  NSE nse(net.GetNuclideLibrary(), &helm, &screen, nse_opts);
  std::cout << "max iterations: " << nse.GetNSEOptions().MaxIterations << std::endl;
  auto nseResult = nse.CalcFromTemperatureAndDensity(T0, rho0, Ye);

  LinearExpansion rho_profile(rho0, tau);
  auto output = net.EvolveSelfHeatingWithInitialTemperature(
    nseResult.Y(), 0.0, 1.0, T0, &rho_profile, "history");

  // Save results
  std::vector<double> finalYVsA = output.FinalYVsA();
  FILE * f = fopen("final_abundance.txt", "w");
  for (unsigned int A = 0; A < finalYVsA.size(); ++A)
    fprintf(f, "%6i %30.20E\n", A, finalYVsA[A]);
}
