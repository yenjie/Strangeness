#include "StrangenessMessenger.h"

void TestStrangeness()
{
   TFile file("../../sample/Strangeness/merged_mc_v2.root");

   StrangenessTreeMessenger M(file, "Tree");

   long long n = M.GetEntries();
   for(long long i = 0; i < n; i++)
   {
      if(!M.GetEntry(i)) continue;
      if(M.PassAll == 0) continue;  // for example

      // Now you can access, e.g.
      // M.NKShort, M.KShortPx[j], M.NPhi, M.PhiPx[k], etc.
   }
}
