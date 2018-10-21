import subprocess
import csv
import re
import numpy as np

def format_args(args):
    return ["-D{0!s}={1!s}".format(name, value) for name, value in args.items()]

if __name__ == "__main__":
    runs = 50
    n_jobs = 4

    # Function to run.
    function_name = "SchaffersEvaluation"

    base_population_count = 690
    base_offspring_count = 3298
    migration_fraction = 0.2

    # Change this to actual counts
    island_counts = [1, 2, 5, 10, 20, 30, 40, 50, 75, 100, 200, 500]

    migration_intervals = [25, 50, 75, 100]

    results_file = open("results-{0!s}.csv".format(function_name), "w")
    results_writer = csv.DictWriter(results_file, fieldnames=["islandCount", "migrationInterval", "run", "generation", "diversity",])
    results_writer.writeheader()

    for island_count in island_counts:
        for migration_interval in migration_intervals:
            parameter_settings = {
                "islandCount": island_count,
                "populationCount": int(np.maximum(1, base_population_count // island_count)),
                "offspringCount": int(np.maximum(1, base_offspring_count // island_count)),
                "migrationCount": int(np.maximum(1, (base_population_count // island_count) * migration_fraction)),
                "migrationInterval": migration_interval
            }

            for i in range(0, runs, n_jobs):

                results = []

                for _ in range(0, min(runs - i, n_jobs)):
                    processes = []

                    for j in range(np.minimum(n_jobs, runs - i)):
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

                        r = []

                        # Parse the output to get the required results.

                        for line in output:
                            m = re.match(r"Diversity for generation (?P<generation>\d+) is (?P<diversity>.*)", line.decode())

                            if m:
                                generation = int(m.group("generation"))
                                diversity = float(m.group("diversity"))

                                r.append({
                                    "generation": generation,
                                    "diversity": diversity
                                })

                        results.append(r)

                for i, result in enumerate(results):
                    # Here you write out the results.

                    for r in result:
                        results_writer.writerow({
                            "islandCount": island_count,
                            "migrationInterval": migration_interval,
                            "run": i,
                            **r
                        })

                    results_file.flush()