#include <iostream>
#include <cmath>
#include <ctime>

#include "Progress.h"

Progress_Bar::Progress_Bar()
{
	progress = 0;
	bar_width = 55;
	start_time = time(NULL);
	elapsed_time = 0;
	N_stats = 0;
	update(0);
}

Progress_Bar::~Progress_Bar()
{
	if (not Progress_Bar_visible)
		return;
	update(1);
	elapsed_time = time(NULL) - start_time;
	std::cerr << std::endl;
	std::cerr << "Elapsed time: " << elapsed_time << " ";
	if (Progress_Bar_statistics)
	{
		double sum = 0;
		double sumsq = 0;
		double res;
		for (int i = 0; i < N_stats; i++)
		{
			res = times[0][i] + times[1][i] - elapsed_time;
			sum += res;
			sumsq += pow(res, 2);
		}
		std::cerr << "Mean: " << sum / N_stats << " STD: " << sqrt(sumsq / N_stats - pow(sum / N_stats, 2)) << " N: " << N_stats;
	}
	std::cerr << std::endl;
}

void Progress_Bar::update(double _progress)
{
	if (not Progress_Bar_visible)
		return;
	if (std::abs(_progress - progress) < 0.009 and _progress < 1)
		return;

	progress = _progress;
	elapsed_time = time(NULL) - start_time;
	remaining_time = elapsed_time * (1. - progress) / progress;

	if (Progress_Bar_statistics and fmod(elapsed_time, 5) == 0)
	{
		if (N_stats < 100)
		{
			times[0][N_stats] = elapsed_time;
			times[1][N_stats] = remaining_time;
			N_stats++;
		}
	}

	std::cerr << "  " << int(progress * 100.) << " % [";
	int position = bar_width * progress;
	for (int i = 0; i < bar_width; ++i)
	{
		if (i < position)
			std::cerr << "=";
		else if (i == position)
			std::cerr << ">";
		else
			std::cerr << " ";
	}
	std::cerr << "] ";
	if (progress > 0.15 or elapsed_time > 10)
		std::cerr << int(remaining_time) << "   ";
	std::cerr << "\r";
	std::cerr.flush();
}

void Progress_Bar::update(double min, double max, double val, bool linear = true)
{
	if (linear)
		update((val - min) / (max - min));
	else
		update(log(val / min) / log(max / min));
}

void set_Progress_Bar_visibility(bool b)
{
	Progress_Bar_visible = b;
}

void set_Progress_Bar_statistics(bool b)
{
	Progress_Bar_statistics = b;
}

bool Progress_Bar_visible = false;
bool Progress_Bar_statistics = false;
