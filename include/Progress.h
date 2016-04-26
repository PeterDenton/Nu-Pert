#ifndef Progress_H
#define Progress_H

class Progress_Bar
{
	public:
		Progress_Bar();
		virtual ~Progress_Bar();
		void update(double progress); // progress should be in 0-1

		// if linear then assume linear scale, if not linear, assume log scale
		// that is, x += 0.1 or x *= 1.1
		void update(double min, double max, double val, bool linear);
	private:
		int bar_width;
		double progress, start_time, elapsed_time, remaining_time;

		// for stats
		double times[2][100];
		int N_stats;
};

void set_Progress_Bar_visibility(bool b);
void set_Progress_Bar_statistics(bool b);

extern bool Progress_Bar_visible;
extern bool Progress_Bar_statistics;

#endif
