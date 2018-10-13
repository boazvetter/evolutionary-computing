import org.vu.contest.ContestSubmission;
import org.vu.contest.ContestEvaluation;

import java.util.Random;
import java.util.Properties;
import java.util.Arrays;
import java.lang.Math;

public class player66 implements ContestSubmission
{
	Random rnd_;
    private ContestWrapper _contest;

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

		// Get evaluation properties
		Properties props = evaluation.getProperties();
        // Get evaluation limit
        int evaluations_limit_ = Integer.parseInt(props.getProperty("Evaluations"));
		// Property keys depend on specific evaluation
		// E.g. double param = Double.parseDouble(props.getProperty("property_name"));
        boolean isMultimodal = Boolean.parseBoolean(props.getProperty("Multimodal"));
        boolean hasStructure = Boolean.parseBoolean(props.getProperty("Regular"));
        boolean isSeparable = Boolean.parseBoolean(props.getProperty("Separable"));

        this._contest = new ContestWrapper(evaluation, evaluations_limit_);

		// Do sth with property values, e.g. specify relevant settings of your algorithm
        if(isMultimodal){
            // Do sth
        }else{
            // Do sth else
        }
    }

	public void run()
	{
		// Initialization
        IMutationOperator mutationOperator = new SelfAdaptiveMutation(0.07, 0.22);
        ICrossOverOperator crossOverOperator = new OnePointCrossOver();
		Instance[] population = init_population(100);

        
        // calculate fitness       
        while (!this._contest.isDone()) {
            for (int i = 0; i < population.length && !this._contest.isDone(); i += 1) {
                this._contest.evaluate(population[i]);
            }

            // Sort the population according to their fitness.
            // They are sorted from worst to best.                   
            Arrays.sort(population);
            // System.out.println("ranked: ");
            // for (int i = 0; i < population.length; i ++){
            // 	population[i].getInfo();
            // }

            Instance[] new_population = new Instance[population.length];
            System.arraycopy(population, 0, new_population, 0, population.length);

            // Now we make a new population through parent selection and cross over.
            int[] selections = rank_based_selection(population, 50);

            for (int i = 0; i < selections.length; i += 2) {
                Instance[] children = crossOverOperator.crossOver(
                    new Instance[] {
                        population[selections[i]],
                        population[selections[i + 1]]
                    },
                    this.rnd_
                );

                new_population[i] = children[0];
                new_population[i + 1] = children[1];
            }

            population = new_population;

            // Perform mutation
            for (int i = 0; i < population.length; i += 1) {
                population[i].mutate(this.rnd_, mutationOperator);
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
            double[] child = new double[10];
            double[] mutationRates = new double[10];
			for (int j=0; j < 10; j++){
                child[j] = rnd_.nextDouble() * 10.0 - 5.0;
                mutationRates[j] = rnd_.nextGaussian();
			}
			population[x] = new Instance(child, mutationRates);
		}
		return population;
	}
}

class ContestWrapper
{
    private ContestEvaluation _contest;
    private int _evaluationLimit;
    private int _evauationCount;

    public ContestWrapper(ContestEvaluation contest, int evaluationLimit) {
        this._contest = contest;
        this._evaluationLimit = evaluationLimit;
        this._evauationCount = 0;
    }

    public boolean isDone() {
        return this._evaluationLimit == this._evauationCount;
    }

    public double evaluate(Instance instance) {
        if (!instance.hasFitness()) {
            this._evauationCount += 1;
            return instance.calculate_fitness(this._contest);
        }
        else {
            return instance.getFitness();
        }
    }
}

interface ICrossOverOperator
{
    public Instance[] crossOver(Instance[] parents, Random rnd);
}

class OnePointCrossOver implements ICrossOverOperator
{
    public Instance[] crossOver(Instance[] parents, Random rnd) {
        int geneCount = parents[0].getGenes().length;
        int cross_over_point = rnd.nextInt(geneCount);       

        double[] c1 = new double[geneCount];
        double[] c2 = new double[geneCount];
        double[] m1 = new double[geneCount];
        double[] m2 = new double[geneCount];

        for (int i = 0; i < cross_over_point; i += 1) {
            c1[i] = parents[0].getGenes()[i];
            c2[i] = parents[1].getGenes()[i];
            m1[i] = parents[0].getMutationRates()[i];
            m2[i] = parents[1].getMutationRates()[i];
        }

        for (int i = cross_over_point; i < geneCount; i += 1) {
            c1[i] = parents[1].getGenes()[i];
            c2[i] = parents[0].getGenes()[i];
            m1[i] = parents[1].getMutationRates()[i];
            m2[i] = parents[0].getMutationRates()[i];
        }
        
        Instance[] children = new Instance[2];

        children[0] = new Instance(c1, m1);
        children[1] = new Instance(c2, m2);


        return children;   
    }
}

interface IMutationOperator
{
    public void mutate(double[] genes, double[] mutationRates, Random rnd);
}

class IdentityMutation implements IMutationOperator {
    public IdentityMutation() {

    }

    public void mutate(double[] genes, double[] mutationRates, Random rnd) {
        // Here we just do nothing.
    }
}

class SelfAdaptiveMutation implements IMutationOperator
{
    private double _tau;
    private double _tauPrime;

    public SelfAdaptiveMutation(double tau, double tauPrime) {
        this._tau = _tau;
        this._tauPrime = tauPrime;
    }

    public void mutate(double[] genes, double[] mutationRates, Random rnd) {
        // Perform self adaptive mutation.

        // First we mutate the array of mutation rates.
        double globalMutationRate = this._tauPrime * rnd.nextGaussian();
        for (int i = 0; i < mutationRates.length; i += 1) {
            double individualMutationRate = this._tau * rnd.nextGaussian();
            mutationRates[i] = mutationRates[i] * Math.exp(globalMutationRate + individualMutationRate);
        }

        // Now that we have adjusted the mutation rates for each gene we can apply those.
        for (int i = 0; i < genes.length; i += 1) {
            genes[i] = Utils.clamp(
                Constants.MIN_VALUE,
                Constants.MAX_VALUE,
                genes[i] + mutationRates[i] * rnd.nextGaussian()
            );
        }
    }
}

class Utils {
    public static double clamp(double lowerBound, double upperBound, double value) {
        if (value > upperBound) {
            return upperBound;
        }
        else if (value < lowerBound) {
            return lowerBound;
        }
        else {
            return value;
        }
    }
}

class Constants {
    public static double MAX_VALUE = 5.0;
    public static double MIN_VALUE = -5.0;
}

class Instance implements Comparable<Instance>
{
    private double[] _genes;
    private Double _fitness;
    private Instance[] _parents;
    private double[] _mutationRates;

    public Instance(double[] genes, double[] mutationRates) {
        this._genes = genes;
        this._mutationRates = mutationRates;
    }

    public Instance(double[] genes, Instance[] parents, double[] mutationRates) {
        this._genes = genes;
        this._parents = parents;
        this._mutationRates = _mutationRates;
    }

    public boolean hasFitness() {
        return this._fitness != null;
    }

    public double getFitness() {
        return this._fitness;
    }

    public double[] getGenes() {
        return this._genes;
    }

    public double[] getMutationRates() {
        return this._mutationRates;
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

    public void mutate(Random rnd, IMutationOperator mutationOperator) {
        if (this._fitness != null) {
            return;
        }

        mutationOperator.mutate(this._genes, this._mutationRates, rnd);
    }
}
