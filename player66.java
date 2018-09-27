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
        double population[][] = init_population(5);
		//System.out.println(Arrays.toString(population[0]));
		onepointcross(population[1], population[2]);
		mutate(population[0]);
        while(evals<evaluations_limit_){
            // Select parents
            // Apply crossover / mutation operators
            double child[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
            // Check fitness of unknown fuction
            Double fitness = (double) evaluation_.evaluate(child);

            evals++;	
            // Select survivors
        }
        // System.out.println(fitness);

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


	public void onepointcross (double parent1[], double parent2[]){
		// Set a random point for crossover 
		int point = new Random().nextInt(10);
		System.out.println(point);
		// Initialize two children
		double child1[] = new double[10];
		double child2[] = new double[10];
		//Loop te the parents before the point to make 
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
		//Loop te the parents after the point to make 
		// second part of the children
		//for (int i = point; i < 10; i++){
			
	//	}
		//System.out.println("parent 1:");
		//System.out.println(Arrays.toString(parent1));

		//System.out.println("child 1:");
		//System.out.println(Arrays.toString(child1));
	}

	public void mutate(double gene[]){
		//System.out.print("Gene to be mutated: ");
		//System.out.println(Arrays.toString(gene));
		// Probability of mutating a gene
		boolean mutateprob = rnd_.nextInt(5)==0;
		// Uniform mutation of chromosomes
		if(mutateprob){
			//System.out.print("Mutate: ");
			//System.out.println(mutateprob);
			int chromosome = rnd_.nextInt(10);
			gene[chromosome] = rnd_.nextDouble() * 10.0 - 5.0;
			//System.out.print("Same gene after mutation: ");
			//System.out.println(Arrays.toString(gene));
		}
	}


}
