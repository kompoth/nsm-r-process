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


using namespace std;

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
  
  double T0 = 1.2;
  double Rho0 = 1e+8;
  auto t9func = ConstantFunction(T0);
  auto rhofunc = ConstantFunction(Rho0);

  // Initial concentrations
  int N = net.GetNuclideLibrary().NumNuclides();
  cout << "# number of nuclei " << N << endl;
  vector<double> y0(N);
  ifstream abund("abund.dat");
  if (!abund)
    throw runtime_error("error loading initial concentrations");
  string dummy;
  getline(abund, dummy);
  getline(abund, dummy);
  double xsum = 0;
  double ysum = 0;
  map<string, int> nuclide_names = net.GetNuclideLibrary().NuclideIdsVsNames();
  while (abund) {
    string name;
    double x;
    abund >> name >> x;
    if (!abund.fail()) {
      if (name == "neut") name = "n";
      if (name == "h1") name = "p";
      if (name == "h2") name = "d";
      if (name == "h3") name = "t";
      if (nuclide_names.count(name) <= 0) {
        cout << "isotop '" << name << "' does not exist in reaction database" << endl;
        throw runtime_error("Error loading initial concentrations");
      }
      int index = nuclide_names[name];
      y0[index] = x;
      xsum += x;
      ysum += x * net.GetNuclideLibrary().Nuclides()[index].A();
    }
  }
  cout << "# xsum = " << xsum << endl;
  cout << "# ysum = " << ysum << endl;
  for (size_t i = 0; i < y0.size(); i++)
    y0[i] *= 1.0 / ysum;
  
  // Run simulation
  auto output = net.Evolve(y0, 0.0, 1.0, &t9func, &rhofunc, "history");

  // Save results
  std::vector<double> finalYVsA = output.FinalYVsA();
  FILE * f = fopen("final_abundance.txt", "w");
  for (unsigned int A = 0; A < finalYVsA.size(); ++A)
    fprintf(f, "%6i %30.20E\n", A, finalYVsA[A]);
}
