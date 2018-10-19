import skopt
import argparse
import numpy as np
import subprocess
import os
import csv


setting_names = [
    "populationCount",
    "mutationRate",
    "crossOverRate",
    "parentTournamentSize",
    "migrationTournamentSize",
    "offspringCount",
    "islandCount",
    "migrationTournamentSize",
    "islandCount",
    "migrationCount",
    "migrationInterval",
    "tau",
    "tauPrime",
    "adaptationBoundary"
]

def format_args(args):
    return ["-D{0!s}={1!s}".format(name, value) for name, value in args.items()]


ITERATION = 0

def get_func(num_runs, num_jobs, function_name, no_islands):
    if no_islands:
        score_file = open("{0!s}-no-islands.scores.csv".format(function_name), "w")
        settings_file = open("{0!s}-no-islands.settings.csv".format(function_name), "w")
    else:
        score_file = open("{0!s}.scores.csv".format(function_name), "w")
        settings_file = open("{0!s}.settings.csv".format(function_name), "w")

    score_writer = csv.DictWriter(score_file, fieldnames=["iteration", "run", "score"])
    settings_writer = csv.DictWriter(settings_file, fieldnames=["iteration"] + setting_names)

    score_writer.writeheader()
    settings_writer.writeheader()

    def func(parameters):
        global ITERATION

        population_count = parameters[5]
        mutation_rate = parameters[0]
        cross_over_rate = parameters[1]
        parent_tournament_size = int(np.maximum(1, parameters[2] * population_count))
        migration_tournament_size = int(np.maximum(1, parameters[3] * population_count))
        offspring_count = int(np.ceil(parameters[4] * population_count / 2)) * 2
        if no_islands:
            island_count = 1
            migration_interval = 2*20
        else:
            island_count = parameters[6]
            migration_interval = parameters[8]
        migration_count = int(np.round(parameters[7] * population_count))

        tau = parameters[9]
        tau_prime = parameters[10]
        adaptation_boundary = parameters[11]

        parameter_settings = {
            "populationCount": population_count,
            "mutationRate": mutation_rate,
            "crossOverRate": cross_over_rate,
            "parentTournamentSize": parent_tournament_size,
            "migrationTournamentSize": migration_tournament_size,
            "offspringCount": offspring_count,
            "islandCount": island_count,
            "migrationCount": migration_count,
            "migrationInterval": migration_interval,
            "tau": tau,
            "tauPrime": tau_prime,
            "adaptationBoundary": adaptation_boundary
        }

        settings_writer.writerow(
            {
                "iteration": ITERATION,
                **parameter_settings
            }
        )
        settings_file.flush()

        results = []
        for i in range(0, num_runs, num_jobs):
            processes = []

            for j in range(np.minimum(num_jobs, num_runs - i)):
                args = [
                    "docker", "run", "--rm", "-e", "LD_LIBRARY_PATH=/code", "ec",
                    "java",
                ] + format_args(parameter_settings) + [
                    "-jar", "testrun.jar", "-submission=player66",
                    "-evaluation={0!s}".format(function_name),
                    "-seed={0!s}".format(np.random.randint(2**31)),
                ]

                processes.append(
                    subprocess.Popen(
                        args,
                        stdout=subprocess.PIPE
                    )
                )

            for process in processes:
                process.wait()

                output = process.stdout.readlines()

                print(output)

                if len(output) < 2:
                    continue

                results.append(float(output[0][7:]))

        for i, score in enumerate(results):
            score_writer.writerow({
                "iteration": ITERATION,
                "run": i,
                "score": score
            })

        score_file.flush()

        ITERATION += 1

        return 10 - np.mean(results)

    return func

        

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--function", "-f", type=str, required=True,
        choices=["KatsuuraEvaluation", "SphereEvaluation", "BentCigarFunction", "SchaffersEvaluation"],
        help="The function to optimize"
    )
    parser.add_argument(
        "--num-trials", type=int, default=1000, required=False,
        help="The amount of trials of gp minimize to perform."
    )
    parser.add_argument(
        "--num-random-trials", type=int, default=50, required=False,
        help="The amount of trials that are random used for initialization."
    )
    parser.add_argument(
        "--num-runs", type=int, default=100, required=False,
        help="The amount of runs with different random seeds to perform for each trial."
    )
    parser.add_argument(
        "--n-jobs", type=int, default=4, required=False,
        help="The amount of concurrent runs to perform."
    )
    parser.add_argument(
        "--no-islands", default=False, action="store_true",
        help="Disable the island model."
    )

    args = parser.parse_args()

    np.random.seed(42)

    optimization_space = [
        skopt.space.Real(0, 10, "uniform"),  # Mutation rate.
        skopt.space.Real(0, 1, "uniform"),  # Cross over rate.
        skopt.space.Real(0, 1, "uniform"),  # Fraction of population used parent selection.
        skopt.space.Real(0, 1, "uniform"),  # Fraction of population used for migration selection.
        skopt.space.Real(1, 100, "log-uniform"),  # Offspring factor.
        skopt.space.Integer(1, 1000),  # Population size.
        skopt.space.Integer(1, 100),  # Island count.
        skopt.space.Real(0, 1, "uniform"),  # Fraction of population that migrates.
        skopt.space.Integer(1, 10000),  # Migration interval.
        skopt.space.Real(1.0e-15, 10, "log-uniform"),  # The self adaptive tau parameter.
        skopt.space.Real(1.0e-15, 10, "log-uniform"),  # The self adaptive tau prime parameter.
        skopt.space.Real(1.0e-15, 10, "log-uniform")  # The mutation boundary parameter.
    ]

    results = skopt.gp_minimize(
        get_func(args.num_runs, args.n_jobs, args.function, args.no_islands),
        optimization_space,
        n_random_starts=args.num_random_trials,
        n_calls=args.num_trials,
        verbose=True
    )

    print(results)
