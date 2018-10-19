import subprocess
import csv

if __name__ == "__main__":
    runs = 50
    n_jobs = 4

    # Function to run.
    function_name = "SchaffersEvaluation"

    # Change this to actual counts
    island_counts = [10, 20, 30, 40]

    results_file = open("results-{0!s}.csv".format(function_name))
    results_writer = csv.DictWriter(results_file, fieldnames=["islandCount", "run", "score"])
    results_writer.writeheader()    

    for island_count in island_counts:
        results = []

        parameter_settings = {
            "islandCount": island_count
        }

        for i in range(0, runs, n_jobs):

            processes = []

            for _ in range(0, min(runs - i, n_jobs)):
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

                    # Parse the output to get the required results.

                    results.append(float(output[0][7:]))

                for i, result in enumerate(results):
                    # Here you write out the results.
                    results_writer.writerow({
                        "islandCount": island_count,
                        "run": i,
                        "score": result
                    })