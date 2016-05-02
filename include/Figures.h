// Writes data out to make figures
#ifndef Figures_H
#define Figures_H

namespace Figures
{
void Eigenvalues();			// outputs Ye*rho*E, lambda1,2,3 - all to zeroth (first) order
void Phases();				// outputs Ye*rho*E, phi, psi
void Expansion_Parameter();	// outputs Ye*rho*E, |eps|, |eps'|
void Pmu2e_Precision();		// outputs E, P_exact; and the absolute value of the fractional errors of: P_GF0, P_GF1, P_Check2
void Eigenvalues_Bases();	// outputs Ye*rho*E, lambda{a,b,c,-,0,+,1,2,3}
} // namespace Figures

#endif
