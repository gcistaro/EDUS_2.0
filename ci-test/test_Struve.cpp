#include "RytovaKeldysh/RytovaKeldysh.hpp"

int main()
{
  printf("\nStruve\n  x      v=0      v=1      v=2      v=3      v=4      v=5");
  for(double x=-5; x<6; x++)
  {
    printf("\nx=%.1lf",x);
    for(int v=0;v<=5;v++)
      printf(" %8.6lf", struve(x,v));
  }
}