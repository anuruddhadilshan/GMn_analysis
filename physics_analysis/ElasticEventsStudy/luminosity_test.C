#include "exprconstants.h"



void luminosity_test()
{
  std::cout << "Beam current in uA: ";
  double beam_curr{0.};
  std::cin >> beam_curr;

  double luminosity = Exprconstants::lD2_lumifactor*beam_curr;

  std::cout << "Luminosity in s^-1.cm^-2: " << luminosity << '\n';

}
