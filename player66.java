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
        if (isMultimodal && !hasStructure && !isSeparable) {
            // Set the settings for katsuura function.
        }
        else if (!isMultimodal && hasStructure && isSeparable) {
            // Set the settings for the sphere function.
        }
        else if (!isMultimodal && !hasStructure && !isSeparable) {
            // Set the settings for the bent cigar function.
        }
        else {
            // Set the setting for the schaffers function.
        }
    }

	public void run()
	{
		// Initialization
        IMutationOperator mutationOperator = new UniformMutation(0.1); //0.07, 0.22
        ICrossOverOperator crossOverOperator = new SingleArithmeticRecombination(0.3);
        IParentSelectionOperator parentSelectionOperator = new TournamentSelection(5);
        ISurvivorSelectionMethod survivorSelectionMethod = new Genetor();
        Instance[] population = init_population(100);
        int offSpringCount = 50;


        // calculate fitness
        while (!this._contest.isDone()) {
            for (int i = 0; i < population.length && !this._contest.isDone(); i += 1) {
                this._contest.evaluate(population[i]);
            }

            // Sort the population according to their fitness.
            // They are sorted from worst to best.
            Arrays.sort(population);

            Instance[] new_population = new Instance[population.length];
            System.arraycopy(population, 0, new_population, 0, population.length);

            // Parent selection
            Instance[] parents = parentSelectionOperator.selectParents(population, offSpringCount, this.rnd_);

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

            // Parent selection
            population = survivorSelectionMethod.selectSurvivors(population, offspring, this.rnd_);

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
        // Sort them base on their performance.
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
                parents[i].getGenes(),
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
        int geneCount = parents[0].getGenes().length;
        int crossOverPoint = rnd.nextInt(geneCount);

        double[] c1 = new double[geneCount];
        double[] c2 = new double[geneCount];
        double[] m1 = new double[geneCount];
        double[] m2 = new double[geneCount];

        for (int i = 0; i < geneCount; i += 1) {
            if (i == crossOverPoint) {
                c1[i] = parents[0].getGenes()[i] * this._weight + parents[1].getGenes()[i] * (1 - this._weight);
                c2[i] = parents[1].getGenes()[i] * this._weight + parents[0].getGenes()[i] * (1 - this._weight);
                m1[i] = parents[0].getGenes()[i] * this._weight + parents[1].getGenes()[i] * (1 - this._weight);
                m2[i] = parents[1].getGenes()[i] * this._weight + parents[0].getGenes()[i] * (1 - this._weight);
            }
            else {
                c1[i] = parents[0].getGenes()[i];
                c2[i] = parents[1].getGenes()[i];
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
        int geneCount = parents[0].getGenes().length;
        int crossOverPoint = rnd.nextInt(geneCount);

        double[] c1 = new double[geneCount];
        double[] c2 = new double[geneCount];
        double[] m1 = new double[geneCount];
        double[] m2 = new double[geneCount];

        for (int i = 0; i < crossOverPoint; i += 1) {
            c1[i] = parents[0].getGenes()[i];
            c2[i] = parents[1].getGenes()[i];
            m1[i] = parents[0].getMutationRates()[i];
            m2[i] = parents[1].getMutationRates()[i];
        }

        for (int i = crossOverPoint; i < geneCount; i += 1) {
            c1[i] = parents[0].getGenes()[i] * this._weight + parents[1].getGenes()[i] * (1 - this._weight);
            c2[i] = parents[1].getGenes()[i] * this._weight + parents[0].getGenes()[i] * (1 - this._weight);
            m1[i] = parents[0].getGenes()[i] * this._weight + parents[1].getGenes()[i] * (1 - this._weight);
            m2[i] = parents[1].getGenes()[i] * this._weight + parents[0].getGenes()[i] * (1 - this._weight);
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
        int geneCount = parents[0].getGenes().length;

        double[] c1 = new double[geneCount];
        double[] c2 = new double[geneCount];
        double[] m1 = new double[geneCount];
        double[] m2 = new double[geneCount];

        for (int i = 0; i < geneCount; i += 1) {
            c1[i] = parents[0].getGenes()[i] * this._weight + parents[1].getGenes()[i] * (1 - this._weight);
            c2[i] = parents[1].getGenes()[i] * this._weight + parents[0].getGenes()[i] * (1 - this._weight);
            m1[i] = parents[0].getGenes()[i] * this._weight + parents[1].getGenes()[i] * (1 - this._weight);
            m2[i] = parents[1].getGenes()[i] * this._weight + parents[0].getGenes()[i] * (1 - this._weight);
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
        int geneCount = parents[0].getGenes().length;

        double[] c1 = new double[geneCount];
        double[] c2 = new double[geneCount];
        double[] m1 = new double[geneCount];
        double[] m2 = new double[geneCount];

        for (int i = 0; i < geneCount; i += 1) {
            if (rnd.nextDouble() < this._crossOverProbabilty) {
                c1[i] = parents[0].getGenes()[i];
                c2[i] = parents[1].getGenes()[i];
                m1[i] = parents[0].getMutationRates()[i];
                m2[i] = parents[1].getMutationRates()[i];
            }
            else {
                c1[i] = parents[1].getGenes()[i];
                c2[i] = parents[0].getGenes()[i];
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

// Implements genetor survivor selection method I.E. replace worst.
class Genetor implements ISurvivorSelectionMethod
{
    public Instance[] selectSurvivors(Instance[] parents, Instance[] children, Random rnd)
    {
        Instance[] new_population = new Instance[parents.length];

        Arrays.sort(parents);

        // Copy over the children.
        for (int i = 0; i < children.length; i += 1) {
            new_population[i] = children[i];
        }

        // Copy over the best from the parents.
        for (int i = children.length; i < parents.length; i += 1) {
            new_population[i] = parents[i];
        }

        return new_population;
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
        mutationOperator.mutate(this._genes, this._mutationRates, rnd);
    }
}
