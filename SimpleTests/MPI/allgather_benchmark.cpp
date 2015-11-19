/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <chrono>
#include <vector>

#include <tclap/CmdLine.h>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/variance.hpp>
namespace ba = boost::accumulators;

#include "BaseLib/BuildInfo.h"
#include "BaseLib/RunTime.h"
#include "BaseLib/CPUTime.h"

#include "Applications/ApplicationsLib/LinearSolverLibrarySetup.h"
#include "Applications/ApplicationsLib/LogogSetup.h"

#include "ProcessLib/NumericsConfig.h"

void benchmark0(std::size_t const, GlobalSetupType)
{
}

void benchmark1(std::size_t const global_vector_size, GlobalSetupType _global_setup)
{
	auto _x = _global_setup.createVector(global_vector_size);

	//int rank;
	//MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//INFO("[%d] %u global, %u local", rank, _x->size(), _x->getLocalSize());
	//INFO("[%d] from %u to %u", rank, _x->getRangeBegin(), _x->getRangeEnd());
	delete _x;
}

void benchmark2(std::size_t const global_vector_size,
                GlobalSetupType _global_setup)
{
	auto _x = _global_setup.createVector(global_vector_size);
	std::vector<double> u(_x->size());
	delete _x;
}

void benchmark3(std::size_t const global_vector_size,
                GlobalSetupType _global_setup)
{
	auto _x = _global_setup.createVector(global_vector_size);
	std::vector<double> u(_x->size());
	_x->getValues(&u[0]);
	delete _x;
}

template <typename F, typename... Args>
std::chrono::duration<double> timeit(F const& f, Args&&... args)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	BaseLib::CPUTime cpu_time;
	BaseLib::RunTime run_time;

	MPI_Barrier(MPI_COMM_WORLD);
	auto start = std::chrono::high_resolution_clock::now();
	cpu_time.start();
	run_time.start();
	f(std::forward<Args>(args)...);
	auto end = std::chrono::high_resolution_clock::now();
	INFO("[%d] chrono %ld, cputime %.16lf, runtime %.16lf", rank, (end - start).count(),
	     cpu_time.elapsed(), run_time.elapsed());
	return end - start;
}

template <typename F, typename... Args>
std::vector<double> collect_stats(int n_tries, F const& f, Args&&... args)
{
	std::vector<double> times;
	times.reserve(n_tries);
	while (n_tries-- > 0)
	{
		auto const t = timeit(f, std::forward<Args>(args)...).count();
		times.push_back(t);
	}
	return times;
}

// Returns pair of statistic names and corresponding values.
std::pair<std::vector<std::string>, std::vector<double>> evaluate_stats(
    std::vector<double> const& data)
{
	ba::accumulator_set<double,
	                    ba::stats<ba::tag::min,
	                              ba::tag::mean,
	                              ba::tag::median,
	                              ba::tag::max,
	                              ba::tag::variance>> acc;

	for (auto d : data)
		acc(d);

	return {{"min", "mean", "median", "max", "variance"},
	        {ba::extract::min(acc), ba::extract::mean(acc),
	         ba::extract::median(acc), ba::extract::max(acc),
	         ba::extract::variance(acc)}};
}

void print_statistics(std::vector<double> const& times)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int n_processes;
	MPI_Comm_size(MPI_COMM_WORLD, &n_processes);

	/*
	{
		MPI_Barrier(MPI_COMM_WORLD);
		// This are per rank statistics.
		auto stats = evaluate_stats(times);
		// Print names.
		if (rank == 0)
		{
			std::string print_names;
			for (auto n : stats.first)
				print_names += n + " ";
			INFO("[ ] : %s", print_names.c_str());
		}
		MPI_Barrier(MPI_COMM_WORLD);

		// Print values for each rank.
		std::string print_values;
		for (auto s : stats.second)
			print_values.append(std::to_string(s) + " ");
		INFO("[%d] : %s", rank, print_values.c_str());
	}
	*/

	// Gather times from all processes and evaluate total statistics.
	std::vector<double> all_times;
	if (rank == 0) all_times.resize(times.size() * n_processes);

	MPI_Gather(times.data(), times.size(), MPI_DOUBLE, all_times.data(),
	           times.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// Print all times' statistics.
	if (rank == 0)
	{
		auto stats = evaluate_stats(all_times);
		// Print values for each rank.
		std::ostringstream print_values;
		print_values << std::setprecision(
		    std::numeric_limits<long double>::digits10 + 1);
		for (auto s : stats.second)
			print_values << s << " ";
		INFO("[ ] : %s", print_values.str().c_str());
	}
}

int main(int argc, char* argv[])
{
	// Parse CLI arguments.
	TCLAP::CmdLine cmd(
	    "OpenGeoSys-6 software.\n"
	    "Copyright (c) 2012-2015, OpenGeoSys Community "
	    "(http://www.opengeosys.org) "
	    "Distributed under a Modified BSD License. "
	    "See accompanying file LICENSE.txt or "
	    "http://www.opengeosys.org/project/license",
	    ' ',
	    BaseLib::BuildInfo::git_describe);

	TCLAP::UnlabeledValueArg<std::size_t> global_size_arg(
	    "global_size", "Size of the global vector", true, 10, "INT");
	cmd.add(global_size_arg);

	TCLAP::UnlabeledValueArg<std::size_t> n_tries_arg(
	    "n_tries", "Number of repetitions.", false, 10, "INT");
	cmd.add(n_tries_arg);

	cmd.parse(argc, argv);

	ApplicationsLib::LogogSetup logog_setup;
	ApplicationsLib::LinearSolverLibrarySetup linear_solver_library_setup(argc,
	                                                                      argv);
	GlobalSetupType global_setup;

	int const n_tries = n_tries_arg.getValue();
	print_statistics(collect_stats(n_tries, benchmark0,
	                               global_size_arg.getValue(), global_setup));

	print_statistics(collect_stats(n_tries, benchmark1,
	                               global_size_arg.getValue(), global_setup));

	print_statistics(collect_stats(n_tries, benchmark2,
	                               global_size_arg.getValue(), global_setup));

	print_statistics(collect_stats(n_tries, benchmark3,
	                               global_size_arg.getValue(), global_setup));
	return 0;
}
