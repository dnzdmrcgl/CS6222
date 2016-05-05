#include <Rcpp.h>
using namespace Rcpp;

int binOffsets[] = {512+64+8+1, 64+8+1, 8+1, 1, 0};

// [[Rcpp::export]]
int binFromRangeStandard(int start, int end) {
  int startBin = start, endBin = end-1, i;
  startBin >>= 17;
  endBin >>= 17;
  for (i=0; i<5; ++i)
  {
    if (startBin == endBin)
      return binOffsets[i] + startBin;
    startBin >>= 3;
    endBin >>= 3;
  }
  return 0;
}
