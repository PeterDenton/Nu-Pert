#ifndef Verify_H
#define Verify_H

namespace Verify
{
bool is_accurate(double P1, double P2, double Precision = 1e-2);
void Precision_Scan(double (*P_To_Verify)(int, int, double, double, double, int), int order = 0, double Precision = 1e-2);
} // namespace Verify

#endif
