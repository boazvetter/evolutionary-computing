#
# Run this script by first running
#  > make player66.class
#  > make submission
#  > docker build -t ec_diversity:v1 .
# and then
#  > python3 diversity_results.py -f SphereEvaluation
# or the general form:
#  > python3 diversity_results.py -f functionname
#

import skopt
import argparse
import numpy as np
import subprocess
import os
import csv
import json


def format_args(args):
    return ["-D{0!s}={1!s}".format(name, value) for name, value in args.items()]


def run_island_model(island_count, function, n_jobs, num_runs):


    parameter_settings = {
        "islandCount": island_count
    }

    results = []
    for i in range(0, num_runs, n_jobs):
        processes = []

        for j in range(np.minimum(n_jobs, num_runs - i)):
            args = [
                "docker", "run", "--rm", "-e", "LD_LIBRARY_PATH=/code", "ec_diversity:v1",
                "java",
            ] + format_args(parameter_settings) + [
                "-jar", "testrun.jar", "-submission=player66",
                "-evaluation={0!s}".format(function),
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

            if len(output) < 2:
                continue

            results.append(json.loads(output[0]))

    mean_results = np.mean(results, axis=0)

    for t in range(len(results[0])):
        score_writer.writerow({
            "island_count": island_count,
            "time": t,
            "score": mean_results[t]
        })

    score_file.flush()






if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--function", "-f", type=str, required=True,
        choices=["KatsuuraEvaluation", "SphereEvaluation", "BentCigarFunction", "SchaffersEvaluation"],
        help="The function to optimize"
    )
    parser.add_argument(
        "--n-jobs", type=int, default=4, required=False,
        help="The amount of concurrent runs to perform."
    )
    parser.add_argument(
        "--num-runs", type=int, default=8, required=False,
        help="The amount of runs with different random seeds to perform for each trial."
    )

    args = parser.parse_args()

    np.random.seed(42)


    score_file = open("{0!s}-variation.csv".format(args.function), "w")
    score_writer = csv.DictWriter(score_file, fieldnames=["island_count", "time", "score"])
    score_writer.writeheader()

    for island_count in range(1,100):
        print("Island count: ", island_count)
        run_island_model(island_count, args.function, args.n_jobs, args.num_runs)
