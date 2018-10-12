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
		// Initialization 
        int evals = 0;
        String mutationForm = "nonuniform"; // Select the kind of mutation
		double learningRate = 0.05;
		Instance[] population = init_population(100);

        
        // calculate fitness       
        while (evals < this.evaluations_limit_) {
            for (int i = 0; i < population.length && evals < this.evaluations_limit_; i += 1) {
                population[i].calculate_fitness(this.evaluation_);
                evals += 1;
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

            // Perform mutation
            for (int i = 0; i < population.length; i += 1) {
                population[i].mutate(this.rnd_, mutationForm, learningRate);
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
		double stepSize;

		for(int x=0; x < n; x++){
			double child[] = new double[10];
			for (int j=0; j < 10; j++){
				child[j] = rnd_.nextDouble() * 10.0 - 5.0;
			}
			stepSize = rnd_.nextDouble();
			population[x] = new Instance(child, stepSize);
		}
		return population;
	}
}

class Instance implements Comparable<Instance>
{
    private double[] _genes;
    private Double _fitness;
    private Instance[] _parents;
    private double _stepSize;

    public Instance(double[] genes, double stepSize) {
        this._genes = genes;
        this._stepSize = stepSize;
    }

    public Instance(double[] genes, Instance[] parents, double stepSize) {
        this._genes = genes;
        this._parents = parents;
        this._stepSize = stepSize;
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

    public void mutate(Random rnd, String mutationForm, double learningRate) {
        if (this._fitness != null) {
            return;
        }

		// Initialize
        boolean mutateprob = true; // p = 1
        //boolean mutateprob = rnd.nextInt(10)==0; // p = 1/10th
		
		// Uniform mutation of alleles <x1,...,xn>
		if (mutateprob) {			
            int allele = rnd.nextInt(this._genes.length); // Select one allele at random
            switch (mutationForm.toLowerCase()) {
            	case "uniform":
            		double lowerlim = 0.75;
					double upperlim = 1.25;
            		this._genes[allele] = this._genes[allele] * (lowerlim + (rnd.nextDouble() * (upperlim - lowerlim)));
            		break;
            	case "nonuniform":
            		// Sigma is mutated every time step by multiplying it by a term e ^ (random variable from normal distribution with mean 0 and standard deviation t)
            		double sigmaPrime = this._stepSize * Math.exp(learningRate * rnd.nextGaussian()); // Equation 4.2 in book
            		//System.out.println(sigmaPrime);
            		double mutationFactor = (sigmaPrime * rnd.nextGaussian()) + 1;  // Normal distribution multiplied with step size and a mean of '1'
            		this._genes[allele] = this._genes[allele] * mutationFactor; // Equation 4.3 in book
            		this._stepSize = sigmaPrime;
            		break;
            }
			
            // Curtail if the mutation is outside of bounds
            if(this._genes[allele] > 5){
                    this._genes[allele] = 5;
                }
                else if(this._genes[allele] < -5){
                    this._genes[allele] = -5;
                }       
        }
        
        this._fitness = null;
    }

    public void getInfo() {
    	System.out.print("Current objects' allele values: ");
    	System.out.print(this);
    	System.out.print(" || Step size: ");
    	System.out.print(this._stepSize);
    	System.out.print(" || Fitness: ");
    	System.out.println(this._fitness);
    	for (int i = 0; i < this._genes.length ; i++){
    		System.out.print(" ");
    		System.out.println(this._genes[i]);
    	}
    	System.out.println("");
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

        children[0] = new Instance(c1, p1._stepSize); // TO BE CHANGED - WHOSE STEP SIZE DOES THE CHILD INHERIT?
        children[1] = new Instance(c2, p2._stepSize); // TO BE CHANGED - WHOSE STEP SIZE DOES THE CHILD INHERIT?


        return children;
    }
}
