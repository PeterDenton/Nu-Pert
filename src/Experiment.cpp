/*
Nominal values for various long baseline experiments.

exp_num:
	0 - DUNE
	1 - NOvA
	2 - T2K/T2HK
*/
#include <cassert>

#include "Experiment.h"

void set_experimental_parameters(int exp_num, double *Yrho, double *L)
{
	assert (exp_num >= 0 && exp_num <= 2); // 0: DUNE, 1: NOvA, 2: T2K/T2HK
	switch (exp_num)
	{
		case 0: // DUNE
			*Yrho = 1.4; // YrhoE over E
			*L = 1300;
			break;
		case 1: // NOvA
			*Yrho = 1.4; // YrhoE over E
			*L = 810;
			break;
		case 2: // T2K/T2HK
			*Yrho = 1.15; // YrhoE over E
			*L = 295;
			break;
	}
}
