#ifndef Warning_H
#define Warning_H

#include <string>

extern int N_Warning;
extern int N_IndexWarning;
extern int N_EdgeCaseWarning;
extern int N_ProbabilityWarning;

class Warning
{
	protected:
		std::string name;
		std::string message;
		bool warn;
		int Warning_number;
	public:
		Warning();
		~Warning();
};
class IndexWarning : public Warning
{
	private:
		int min, max;
		double val;
	public:
		IndexWarning(int min, int max, double val);
};
class EdgeCaseWarning : public Warning
{
	private:
		double min, max, val;
	public:
		EdgeCaseWarning(double min, double max, double val);
};
class ProbabilityWarning : public Warning
{
	public:
		ProbabilityWarning(double *P);
};

#endif
