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
    IMutationOperator mutationOperator;

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
        if (isMultimodal && !hasStructure && !isSeparable) {
            // Set the settings for katsuura function.
            mutationOperator = new UniformMutation(Double.parseDouble(System.getProperty("mutationRate")));
        }
        else if (!isMultimodal && hasStructure && isSeparable) {
            // Set the settings for the sphere function.
            mutationOperator = new SelfAdaptiveMutation(0.07, 0.22, 1.0e-9);
        }
        else if (!isMultimodal && !hasStructure && !isSeparable) {
            // Set the settings for the bent cigar function.
            mutationOperator = new SelfAdaptiveMutation(0.07, 0.22, 1.0e-9);
        }
        else {
            // Set the setting for the schaffers function.
            mutationOperator = new SelfAdaptiveMutation(0.07, 0.22, 1.0e-9);
        }
    }

	public void run()
	{
		// Initialization

        ICrossOverOperator crossOverOperator = new UniformCrossOver(Double.parseDouble(System.getProperty("crossOverRate")));
        IParentSelectionOperator parentSelectionOperator = new TournamentSelection(Integer.parseInt(System.getProperty("parentTournamentSize")));
        ISurvivorSelectionMethod survivorSelectionMethod = new CommaSelection();
        IParentSelectionOperator migrationSelector = new TournamentSelection(Integer.parseInt(System.getProperty("migrationTournamentSize")));
        int offspringCount = Integer.parseInt(System.getProperty("offspringCount"));
        int populationCount = Integer.parseInt(System.getProperty("populationCount"));
        int islandCount = Integer.parseInt(System.getProperty("islandCount"));
        int migrationCount = Integer.parseInt(System.getProperty("migrationCount"));
        int migrationInterval = Integer.parseInt(System.getProperty("migrationInterval"));
        Instance[][] islands = new Instance[islandCount][];
        
        for (int i = 0; i < islandCount; i += 1) {
            islands[i] = init_population(populationCount);
        }

        // calculate fitness
        int generationCount = 0;
        
        while (!this._contest.isDone()) {
            generationCount += 1;

            for (int k = 0; k < islandCount; k += 1) {
                Instance[] population = islands[k];

                // Evaluate fitness of all population members
                for (int i = 0; i < population.length && !this._contest.isDone(); i += 1) {
                    this._contest.evaluate(population[i]);
                }

                if (this._contest.isDone()){
                    break;
                }

                Instance[] parents;

                if (generationCount % migrationInterval == 0) {
                    // Perform migration.
                    int candidateSize = migrationCount * (islandCount - 1) + population.length;
                    Instance[] candidates = new Instance[candidateSize];

                    for (int i = 0; i < population.length; i += 1) {
                        candidates[i] = population[i];
                    }

                    for (int i = 0; i < islandCount; i += 1) {
                        if (i < k) {
                            Instance[] selection = migrationSelector.selectParents(islands[i], migrationCount, this.rnd_);

                            for (int j = 0; j < selection.length; j += 1) {
                                candidates[population.length + (i * migrationCount) + j] = selection[j];
                            }
                        }
                        else if (i > k) {
                            Instance[] selection = migrationSelector.selectParents(islands[i], migrationCount, this.rnd_);

                            for (int j = 0; j < selection.length; j += 1) {
                                candidates[population.length + ((i - 1) * migrationCount) + j] = selection[j];
                            }
                        }
                    }

                    parents = parentSelectionOperator.selectParents(candidates, offspringCount, this.rnd_);
                }
                else {
                    // Parent selection
                    parents = parentSelectionOperator.selectParents(population, offspringCount, this.rnd_);
                }

                // Crossover
                Instance[] offspring = new Instance[parents.length];
                for (int i = 0; i < parents.length; i += 2) {
                    Instance[] new_offspring = crossOverOperator.crossOver(
                        new Instance[] {
                            parents[i],
                            parents[i + 1]
                        },
                        this.rnd_
                    );

                    offspring[i] = new_offspring[0];
                    offspring[i + 1] = new_offspring[1];
                }

                // Mutation
                for (int i = 0; i < offspring.length; i += 1) {
                    offspring[i].mutate(this.rnd_, mutationOperator);
                }

                // Evaluate the new offspring
                for (int i = 0; i < offspring.length && !this._contest.isDone(); i += 1) {
                    this._contest.evaluate(offspring[i]);
                }
                if (this._contest.isDone()){
                    break;
                }

                // Survivor selection
                population = survivorSelectionMethod.selectSurvivors(population, offspring, this.rnd_);

                islands[k] = population;
            }
        }
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

// Handy class to keep track of the amount of evaluations that have been performed.
class ContestWrapper
{
    private ContestEvaluation _contest;
    private int _evaluationLimit;
    private int _evaluationCount;

    public ContestWrapper(ContestEvaluation contest, int evaluationLimit) {
        this._contest = contest;
        this._evaluationLimit = evaluationLimit;
        this._evaluationCount = 0;
    }

    public boolean isDone() {
        return this._evaluationLimit == this._evaluationCount;
    }

    public double evaluate(Instance instance) {
        if (!instance.hasFitness()) {
            this._evaluationCount += 1;
            return instance.calculate_fitness(this._contest);
        }
        else {
            return instance.getFitness();
        }
    }
}


// Common interface for the parent selection operators.
interface IParentSelectionOperator
{
    public Instance[] selectParents(Instance[] parents, int parentCount, Random rnd);
}

// Implements linear rank based parent selection.
class RankBasedSelection implements IParentSelectionOperator
{
    private double _selectionWeight;

    public RankBasedSelection(double selectionWeight) {
        this._selectionWeight = selectionWeight;
    }

    public Instance[] selectParents(Instance[] population, int parentCount, Random rnd) {

        if (this._selectionWeight < 1.0 || this._selectionWeight > 2.0){
            throw new IllegalArgumentException("Selection pressure should be between 1.0 and 2.0");
        }

        // Sort them base on their performance (from worst to best)
        Arrays.sort(population);

        Instance[] selections = new Instance[parentCount];

        double[] probabilities = new double[population.length];

        for (int i = 0; i < population.length; i += 1) {
            double p = (2 - this._selectionWeight) / population.length;
            probabilities[i] = p + (2 * i * (this._selectionWeight - 1)) / (population.length * (population.length - 1));
        }


        for (int i = 0; i < parentCount; i += 1) {
            selections[i] = population[Utils.rouletteWheelSelection(probabilities, rnd)];
        }

        return selections;
    }
}

// Implements tournament selection
class TournamentSelection implements IParentSelectionOperator
{
    private int _tournamentSize;

    public TournamentSelection(int tournamentSize) {
        this._tournamentSize = tournamentSize;
    }

    public Instance[] selectParents(Instance[] population, int parentCount, Random rnd) {
        Instance[] selections = new Instance[parentCount];

        for (int i = 0; i < parentCount; i += 1) {
            Instance[] participants = new Instance[this._tournamentSize];

            for (int j = 0; j < participants.length; j += 1) {
                participants[j] = population[rnd.nextInt(population.length)];
            }

            Arrays.sort(participants);

            selections[i] = participants[participants.length - 1];
        }

        return selections;
    }
}

// Implements uniform parent selection.
class UniformSelection implements IParentSelectionOperator
{
    public Instance[] selectParents(Instance[] population, int parentCount, Random rnd) {
        Instance[] selections = new Instance[parentCount];

        for (int i = 0;i < parentCount; i += 1) {
            selections[i] = population[rnd.nextInt(population.length)];
        }

        return selections;
    }
}




// Common interface for the cross over operator.
interface ICrossOverOperator
{
    public Instance[] crossOver(Instance[] parents, Random rnd);
}

// Implements the identity cross over operator.
// Basically clones the parents into children.
class IdentityCrossOver implements ICrossOverOperator
{
    public Instance[] crossOver(Instance[] parents, Random rnd) {
        Instance[] children = new Instance[parents.length];

        for (int i = 0; i < parents.length; i += 1) {
            children[i] = new Instance(
                parents[i].getGene(),
                parents[i].getMutationRates()
            );
        }

        return children;
    }
}

// Implements simple arithmetic recombintation
class SingleArithmeticRecombination implements ICrossOverOperator
{
    private double _weight;

    public SingleArithmeticRecombination(double weight) {
        this._weight = weight;
    }

    public Instance[] crossOver(Instance[] parents, Random rnd) {
        int geneCount = parents[0].getGene().length;
        int crossOverPoint = rnd.nextInt(geneCount);

        double[] c1 = new double[geneCount];
        double[] c2 = new double[geneCount];
        double[] m1 = new double[geneCount];
        double[] m2 = new double[geneCount];

        for (int i = 0; i < geneCount; i += 1) {
            if (i == crossOverPoint) {
                c1[i] = parents[0].getGene()[i] * this._weight + parents[1].getGene()[i] * (1 - this._weight);
                c2[i] = parents[1].getGene()[i] * this._weight + parents[0].getGene()[i] * (1 - this._weight);
                m1[i] = parents[0].getGene()[i] * this._weight + parents[1].getGene()[i] * (1 - this._weight);
                m2[i] = parents[1].getGene()[i] * this._weight + parents[0].getGene()[i] * (1 - this._weight);
            }
            else {
                c1[i] = parents[0].getGene()[i];
                c2[i] = parents[1].getGene()[i];
                m1[i] = parents[0].getMutationRates()[i];
                m2[i] = parents[1].getMutationRates()[i];
            }
        }

        return new Instance[] {
            new Instance(c1, m1),
            new Instance(c2, m2)
        };
    }
}

// Implements simple arithmetic recombintation
class SimpleArithmeticRecombination implements ICrossOverOperator
{
    private double _weight;

    public SimpleArithmeticRecombination(double weight) {
        this._weight = weight;
    }

    public Instance[] crossOver(Instance[] parents, Random rnd) {
        int geneCount = parents[0].getGene().length;
        int crossOverPoint = rnd.nextInt(geneCount);

        double[] c1 = new double[geneCount];
        double[] c2 = new double[geneCount];
        double[] m1 = new double[geneCount];
        double[] m2 = new double[geneCount];

        for (int i = 0; i < crossOverPoint; i += 1) {
            c1[i] = parents[0].getGene()[i];
            c2[i] = parents[1].getGene()[i];
            m1[i] = parents[0].getMutationRates()[i];
            m2[i] = parents[1].getMutationRates()[i];
        }

        for (int i = crossOverPoint; i < geneCount; i += 1) {
            c1[i] = parents[0].getGene()[i] * this._weight + parents[1].getGene()[i] * (1 - this._weight);
            c2[i] = parents[1].getGene()[i] * this._weight + parents[0].getGene()[i] * (1 - this._weight);
            m1[i] = parents[0].getGene()[i] * this._weight + parents[1].getGene()[i] * (1 - this._weight);
            m2[i] = parents[1].getGene()[i] * this._weight + parents[0].getGene()[i] * (1 - this._weight);
        }

        return new Instance[] {
            new Instance(c1, m1),
            new Instance(c2, m2)
        };
    }
}

// Implements whole arithmetic recombintation
class WholeArithmeticRecombination implements ICrossOverOperator
{
    private double _weight;

    public WholeArithmeticRecombination(double weight) {
        this._weight = weight;
    }

    public Instance[] crossOver(Instance[] parents, Random rnd) {
        int geneCount = parents[0].getGene().length;

        double[] c1 = new double[geneCount];
        double[] c2 = new double[geneCount];
        double[] m1 = new double[geneCount];
        double[] m2 = new double[geneCount];

        for (int i = 0; i < geneCount; i += 1) {
            c1[i] = parents[0].getGene()[i] * this._weight + parents[1].getGene()[i] * (1 - this._weight);
            c2[i] = parents[1].getGene()[i] * this._weight + parents[0].getGene()[i] * (1 - this._weight);
            m1[i] = parents[0].getGene()[i] * this._weight + parents[1].getGene()[i] * (1 - this._weight);
            m2[i] = parents[1].getGene()[i] * this._weight + parents[0].getGene()[i] * (1 - this._weight);
        }

        return new Instance[] {
            new Instance(c1, m1),
            new Instance(c2, m2)
        };
    }
}

// Implements the one point cross over operator.
class OnePointCrossOver implements ICrossOverOperator
{
    public Instance[] crossOver(Instance[] parents, Random rnd) {
        int geneCount = parents[0].getGene().length;
        int cross_over_point = rnd.nextInt(geneCount);

        double[] c1 = new double[geneCount];
        double[] c2 = new double[geneCount];
        double[] m1 = new double[geneCount];
        double[] m2 = new double[geneCount];

        for (int i = 0; i < cross_over_point; i += 1) {
            c1[i] = parents[0].getGene()[i];
            c2[i] = parents[1].getGene()[i];
            m1[i] = parents[0].getMutationRates()[i];
            m2[i] = parents[1].getMutationRates()[i];
        }

        for (int i = cross_over_point; i < geneCount; i += 1) {
            c1[i] = parents[1].getGene()[i];
            c2[i] = parents[0].getGene()[i];
            m1[i] = parents[1].getMutationRates()[i];
            m2[i] = parents[0].getMutationRates()[i];
        }

        return new Instance[] {
            new Instance(c1, m1),
            new Instance(c2, m2)
        };
    }
}

// Implement the uniform cross over operation.
class UniformCrossOver implements ICrossOverOperator
{
    private double _crossOverProbabilty;

    public UniformCrossOver(double crossOverProbability) {
        this._crossOverProbabilty = crossOverProbability;
    }

    public Instance[] crossOver(Instance[] parents, Random rnd) {
        int geneCount = parents[0].getGene().length;

        double[] c1 = new double[geneCount];
        double[] c2 = new double[geneCount];
        double[] m1 = new double[geneCount];
        double[] m2 = new double[geneCount];

        for (int i = 0; i < geneCount; i += 1) {
            if (rnd.nextDouble() < this._crossOverProbabilty) {
                c1[i] = parents[0].getGene()[i];
                c2[i] = parents[1].getGene()[i];
                m1[i] = parents[0].getMutationRates()[i];
                m2[i] = parents[1].getMutationRates()[i];
            }
            else {
                c1[i] = parents[1].getGene()[i];
                c2[i] = parents[0].getGene()[i];
                m1[i] = parents[1].getMutationRates()[i];
                m2[i] = parents[0].getMutationRates()[i];
            }
        }

        return new Instance[] {
            new Instance(c1, m1),
            new Instance(c2, m2)
        };
    }
}



// The common interface for the mutation operator.
interface IMutationOperator
{
    public void mutate(double[] genes, double[] mutationRates, Random rnd);
}

// Implements the identity mutation.
class IdentityMutation implements IMutationOperator {
    public void mutate(double[] genes, double[] mutationRates, Random rnd) {
        // Here we just do nothing.
    }
}

// Implements self adaptive mutation with per gene mutation rates.
class SelfAdaptiveMutation implements IMutationOperator {
    private double _tau;
    private double _tauPrime;
    private double _mutationBoundary;

    public SelfAdaptiveMutation(double tau, double tauPrime, double mutationBoundary) {
        this._tau = _tau;
        this._tauPrime = tauPrime;
        this._mutationBoundary = mutationBoundary;
    }

    public void mutate(double[] genes, double[] mutationRates, Random rnd) {
        // Perform self adaptive mutation.

        // First we mutate the array of mutation rates.
        double globalMutationRate = this._tauPrime * rnd.nextGaussian();
        for (int i = 0; i < mutationRates.length; i += 1) {
            double individualMutationRate = this._tau * rnd.nextGaussian();
            mutationRates[i] = mutationRates[i] * Math.exp(globalMutationRate + individualMutationRate);

            if (mutationRates[i] < this._mutationBoundary) {
                mutationRates[i] = this._mutationBoundary;
            }
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

//Implements the Uniform Mutation operator
class UniformMutation implements IMutationOperator {
    private double _mutationRate;

    public UniformMutation(double mutationRate) {
        this._mutationRate = mutationRate;
    }
    // Give the mutationrate, the chance of mutation.
    public void mutate(double[] genes, double[] mutationRates, Random rnd){
        for (int i = 0; i < genes.length; i += 1) {
            double randomNum = rnd.nextDouble(); // rnd.nextInt(1/mutationRate);
                if(randomNum < this._mutationRate ){
                    genes[i] = rnd.nextDouble() * 10.0 - 5.0;
                }
        }
    }
}



// Common interface for survivor selection methods.
interface ISurvivorSelectionMethod
{
    public Instance[] selectSurvivors(Instance[] parents, Instance[] children, Random rnd);
}

// Implements genetor survivor selection method I.E. replace worst (plus selection (μ+λ))
class GenetorSelection implements ISurvivorSelectionMethod
{
    public Instance[] selectSurvivors(Instance[] parents, Instance[] offspring, Random rnd)
    {
        Instance[] final_population = new Instance[parents.length];
        Instance[] temp_population = new Instance[parents.length + offspring.length];

        // Fill temporary population array with both parents and offspring
        System.arraycopy(parents, 0, temp_population, 0, parents.length);
        System.arraycopy(offspring, 0, temp_population, parents.length, offspring.length);

        // Sort the temporary population by fitness (worst to best)
        Arrays.sort(temp_population);

        // Reverse order so it becomes best to worst
        for(int i = 0; i < temp_population.length / 2; i++){
            Instance temp = temp_population[i];
            temp_population[i] = temp_population[temp_population.length - i - 1];
            temp_population[temp_population.length - i - 1] = temp;
        }

        // Only keep first μ best individuals
        // (analogous to throwing away λ worst individuals as discussed in the book)
        System.arraycopy(temp_population, 0, final_population, 0, parents.length);

        return final_population;
    }
}

// Implements comma (μ, λ) survivor selection method
class CommaSelection implements ISurvivorSelectionMethod
{
    public Instance[] selectSurvivors(Instance[] parents, Instance[] offspring, Random rnd)
    {
        Instance[] final_population = new Instance[parents.length];

        // Sort the offspring by fitness (worst to best)
        Arrays.sort(offspring);

        // Reverse order so it becomes best to worst
        offspring = Utils.reversePopulationOrder(offspring);

        // Only keep first μ best individuals of offspring
        System.arraycopy(offspring, 0, final_population, 0, parents.length);

        return final_population;
    }
}


// Utility class to hold usefull functions.
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

    public static int rouletteWheelSelection(double[] probabilities, Random rnd) {
        double prob = rnd.nextDouble();
        double probability_sum = 0.0;

        for (int i = 0 ; i < probabilities.length; i += 1) {
            probability_sum += probabilities[i];

            if (prob < probability_sum) {
                return i;
            }
        }
        return probabilities.length - 1;
    }

    public static Instance[] reversePopulationOrder(Instance[] myPopulation){
        for(int i = 0; i < myPopulation.length / 2; i++){
            Instance temp = myPopulation[i];
            myPopulation[i] = myPopulation[myPopulation.length - i - 1];
            myPopulation[myPopulation.length - i - 1] = temp;
        }
        return myPopulation;
    }
}

// Constants used.
class Constants {
    // Maximum value of a gene.
    public static double MAX_VALUE = 5.0;
    // Minimum value for a gene.
    public static double MIN_VALUE = -5.0;
}

// Class that keeps track of all the information pertaining a single individual.
class Instance implements Comparable<Instance>
{
    private double[] _gene;
    private Double _fitness;
    private Instance[] _parents;
    private double[] _mutationRates;

    public Instance(double[] genes, double[] mutationRates) {
        this._gene = genes;
        this._mutationRates = mutationRates;
    }

    public Instance(double[] genes, Instance[] parents, double[] mutationRates) {
        this._gene = genes;
        this._parents = parents;
        this._mutationRates = _mutationRates;
    }

    public boolean hasFitness() {
        return this._fitness != null;
    }

    public double getFitness() {
        return this._fitness;
    }

    public double[] getGene() {
        return this._gene;
    }

    public double[] getMutationRates() {
        return this._mutationRates;
    }

    public double calculate_fitness(ContestEvaluation evaluation) {
        if (this._fitness == null) {
            this._fitness = (Double)evaluation.evaluate(this._gene);
        }

        return this._fitness;
    }

    public int compareTo(Instance other) {
        return Double.compare(this._fitness, other._fitness);
    }

    public void mutate(Random rnd, IMutationOperator mutationOperator) {
        mutationOperator.mutate(this._gene, this._mutationRates, rnd);
    }
}
