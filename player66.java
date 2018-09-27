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
        double population[][] = init_population(10);
		//System.out.println(Arrays.toString(population[0]));
		// onepointcross(population[1], population[2]);

		//Mutate population
		for (int n = 0; n < population.length; n++){
			mutate(population[n]);
		}

        evaluations_limit_ = 20;
        while(evals<evaluations_limit_){
            // System.out.println(Arrays.deepToString(population));
            // System.out.println(population.length);

            // Select parents (stochastic component, can result in bigger population)
            population = select_parents(population, 10);

            // System.out.println(Arrays.toString(population[0]));


            // Apply crossover / mutation operators
            double child[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
            double childs[][] = new double[2][10];
            childs = onepointcross(population[0], population[1]);
            population[0] = childs[0];
            population[1] = childs[1];

            // System.out.println(Arrays.toString(population[0]));

            // Check fitness of unknown fuction
            // Double fitness = (double) evaluation_.evaluate(childs[0]);

            evals++;
            // Select survivors (deterministic, select n best individuals from population)
        }
        // System.out.println(fitness);

	}


    public double[][] select_parents(double population[][], int nr_parents){

        //
        // CALCULATING FITNESS
        //
        double fitnesses[] = new double[population.length];
        double sum = 0.0;
        for(int i=0; i < population.length; i++){
            fitnesses[i] = -1.0 / (double) evaluation_.evaluate(population[i]);
            sum += fitnesses[i];
        };
        // System.out.print("Fitness of individuals: ");
        // System.out.println(Arrays.toString(fitnesses));
        // System.out.print("Fitness sum of population: ");
        // System.out.println(sum);

        double relative_fitnesses[] = fitnesses;
        for(int i=0; i < population.length; i++){
            relative_fitnesses[i] = fitnesses[i] / sum;
        }
        System.out.println(Arrays.toString(relative_fitnesses));

        //
        // ROULETTE WHEEL PARENT SELECTION
        //
        double parents[][] = new double[nr_parents][10];

        for(int p=0; p<nr_parents; p++){
            double probability_sum = 0.0;
            double prob = rnd_.nextDouble();

            int chosen_parent = 0;
            for(int i=0; i<population.length; i++){
                probability_sum += fitnesses[i];

                if(prob < probability_sum){
                    chosen_parent = i;
                    break;
                }
            }
            relative_fitnesses = [ 0.1 0.05 0.05 0.1 0.7]

            parents[p] = population[chosen_parent];

        }



        // double bestParents[][] = init_population(10);
        // bestParents = Arrays.copyOfRange(population, 0, 5);
        return parents;
    }

	public double[][] init_population(int n){
		double population[][] = new double[n][10];

		for(int x=0; x < n; x++){

			double child[] = new double[10];
			for (int j=0; j < 10; j++){
				child[j] = rnd_.nextDouble() * 10.0 - 5.0;
			}

			population[x] = child;
		}

		return population;
	}


	public double[][] onepointcross (double parent1[], double parent2[]){
		// Set a random point for crossover
		int point = new Random().nextInt(10);
		// System.out.print("Crossover point: ");
  //       System.out.println(point);
		// Initialize two children
		double child1[] = new double[10];
		double child2[] = new double[10];
		//Loop the parents before the point to make
		// first part of the children
		for (int i=0; i < 10; i++){
			if(i < point){
				child1[i] = parent1[i];
				child2[i] = parent2[i];
			}
			else{
				child1[i] = parent2[i];
				child2[i] = parent1[i];
			}
		}

		// System.out.print("parent 1: ");
		// System.out.println(Arrays.toString(parent1));

  //       System.out.print("parent 2: ");
  //       System.out.println(Arrays.toString(parent2));

		// System.out.print("child 1:  ");
		// System.out.println(Arrays.toString(child1));

  //       System.out.print("child 2:  ");
  //       System.out.println(Arrays.toString(child2));

        double childs[][] = new double[2][10];
        childs[0] = child1;
        childs[1] = child2;
        return childs;
	}

	public void mutate(double gene[]){
		//init
		double lowerlim = 0.99;
		double upperlim = 1.01;
		//System.out.print("Gene to be mutated: ");
		//System.out.println(Arrays.toString(gene));
		// Probability of mutating a gene
		boolean mutateprob = rnd_.nextInt(5)==0;
		// Uniform mutation of chromosomes
		if(mutateprob){
			//System.out.print("Mutate: ");
			//System.out.println(mutateprob);
			int chromosome = rnd_.nextInt(10);
			gene[chromosome] = gene[chromosome] * lowerlim + rnd_.nextDouble() * (upperlim - lowerlim);
			//System.out.print("Same gene after mutation: ");
			//System.out.println(Arrays.toString(gene));
		}
	}


}
