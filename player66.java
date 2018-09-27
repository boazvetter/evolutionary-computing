import org.vu.contest.ContestSubmission;
import org.vu.contest.ContestEvaluation;

import java.util.Random;
import java.util.Properties;
import java.util.Arrays;

public class player66 implements ContestSubmission
{
	Random rnd_;
	ContestEvaluation evaluation_;
    private int evaluations_limit_;

	public player66()
	{
		rnd_ = new Random();
	}

	public void setSeed(long seed)
	{
		// Set seed of algortihms random process
		rnd_.setSeed(seed);
	}

	public void setEvaluation(ContestEvaluation evaluation)
	{
		// Set evaluation problem used in the run
		evaluation_ = evaluation;

		// Get evaluation properties
		Properties props = evaluation.getProperties();
        // Get evaluation limit
        evaluations_limit_ = Integer.parseInt(props.getProperty("Evaluations"));
		// Property keys depend on specific evaluation
		// E.g. double param = Double.parseDouble(props.getProperty("property_name"));
        boolean isMultimodal = Boolean.parseBoolean(props.getProperty("Multimodal"));
        boolean hasStructure = Boolean.parseBoolean(props.getProperty("Regular"));
        boolean isSeparable = Boolean.parseBoolean(props.getProperty("Separable"));

		// Do sth with property values, e.g. specify relevant settings of your algorithm
        if(isMultimodal){
            // Do sth
        }else{
            // Do sth else
        }
    }

	public void run()
	{
		// Run your algorithm here


        int evals = 0;
        // init population
        // calculate fitness
        Instance[] population = init_population(10);

        while (evals < this.evaluations_limit_) {
            for (int i = 0; i < population.length && evals < this.evaluations_limit_; i += 1) {
                population[i].calculate_fitness(this.evaluation_);
                evals += 1;
            }

            // Sort the population according to their fitness.
            // They are sorted from worst to best.
            Arrays.sort(population);

            Instance[] new_population = new Instance[population.length];
            System.arraycopy(population, 0, new_population, 0, population.length);

            // Now we make a new population through parent selection and cross over.

            int[] selections = rank_based_selection(population, population.length);

            for (int i = 0; i < selections.length; i += 2) {
                Instance[] children = Instance.one_point_crossover(
                    population[selections[i]],
                    population[selections[i + 1]],
                    this.rnd_
                );

                new_population[i] = children[0];
                new_population[i + 1] = children[1];
            }

            population = new_population;

            for (int i = 0; i < population.length; i += 1) {
                population[i].mutate(this.rnd_);
            }
        }
    }
    
    public int[] rank_based_selection(Instance[] parents, int amount) {
        double[] probabilities = new double[parents.length];

        double s = 1.5;

        for (int i = 0; i < parents.length; i += 1) {
            probabilities[i] = (2 - s) / parents.length + (2 * i * (s - 1)) / (parents.length * (parents.length - 1));
        }

        return select_parents(probabilities, amount);
    }

    public int[] select_parents(double[] probabilities, int amount) {
        int[] selections = new int[amount];

        for (int p=0; p < amount; p += 1){
            double probability_sum = 0.0;
            double prob = rnd_.nextDouble();

            int chosen_parent = 0;
            for (int i = 0 ; i < probabilities.length; i += 1){
                probability_sum += probabilities[i];

                if (prob < probability_sum){
                    chosen_parent = i;
                    break;
                }
            }

            selections[p] = chosen_parent;
        }

        return selections;
    }

	public Instance[] init_population(int n){
		Instance population[] = new Instance[n];

		for(int x=0; x < n; x++){

			double child[] = new double[10];
			for (int j=0; j < 10; j++){
				child[j] = rnd_.nextDouble() * 10.0 - 5.0;
			}

			population[x] = new Instance(child);
		}

		return population;
	}
}

class Instance implements Comparable<Instance>
{
    private double[] _genes;
    private Double _fitness;
    private Instance[] _parents;

    public Instance(double[] genes) {
        this._genes = genes;
    }

    public Instance(double[] genes, Instance[] parents) {
        this._genes = genes;
        this._parents = parents;
    }

    public double calculate_fitness(ContestEvaluation evaluation) {
        if (this._fitness == null) {
            this._fitness = (Double)evaluation.evaluate(this._genes);
        }

        return this._fitness;
    }

    public int compareTo(Instance other) {
        return Double.compare(this._fitness, other._fitness);
    }

    public void mutate(Random rnd) {
        if (this._fitness != null) {
            return;
        }

		//init
		double lowerlim = 0.99;
		double upperlim = 1.01;
		//System.out.print("Gene to be mutated: ");
		//System.out.println(Arrays.toString(gene));
		// Probability of mutating a gene
		boolean mutateprob = rnd.nextInt(5)==0;
		// Uniform mutation of chromosomes
		if (mutateprob) {
			//System.out.print("Mutate: ");
			//System.out.println(mutateprob);
			int chromosome = rnd.nextInt(this._genes.length);
			this._genes[chromosome] = this._genes[chromosome] * lowerlim + rnd.nextDouble() * (upperlim - lowerlim);
			//System.out.print("Same gene after mutation: ");
			//System.out.println(Arrays.toString(gene));
        }
        
        this._fitness = null;
    }

    public static Instance[] one_point_crossover(Instance p1, Instance p2, Random rnd) {
        int cross_over_point = rnd.nextInt(p1._genes.length);

        double[] c1 = new double[p1._genes.length];
        double[] c2 = new double[p1._genes.length];

        for (int i = 0; i < cross_over_point; i += 1) {
            c1[i] = p1._genes[i];
            c2[i] = p2._genes[i];
        }

        for (int i = cross_over_point; i < p1._genes.length; i += 1) {
            c1[i] = p2._genes[i];
            c2[i] = p1._genes[i];
        }
        
        Instance[] children = new Instance[2];

        children[0] = new Instance(c1);
        children[1] = new Instance(c2);

        return children;
    }
}
