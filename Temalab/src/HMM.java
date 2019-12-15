
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

public class HMM {
	private List<Kmer> kmers;
	private double k = 0.9;
	private Integrator integrator;
	private double [][] overlapMatrix;
	private double [][] transitionMatrix;
	private int[] ref;
	private Kmer start;
	private Kmer end;
	
	public HMM(List<Kmer> init) {
		kmers = init;
		for(int i = 0; i < kmers.size(); i++) {
			kmers.get(i).setUpTransitions();
		}
		integrator = new Integrator();
	}
	
	public int[] getRef() { return ref; }
	public double[][] getOverlapMatrix() { return overlapMatrix; }
	public void setOverlapMatrix(double[][] om) { overlapMatrix = om; }
	public double[][] getTransitionMatrix() { return transitionMatrix; }
	public void setTransitionMatrix(double[][] tm) { transitionMatrix = tm; }
	public int getStart() { return getPosByBaseChain(start.getBaseChain()); }
	public int getEnd() { return getPosByBaseChain(end.getBaseChain()); }
	
	public void setUpOverlapMatrix() {
		overlapMatrix = new double[kmers.size()][kmers.size()];
		for(int i = 0; i < kmers.size(); i++) {
			for(int j = 0; j < kmers.size(); j++) {
				if(j > i)
					overlapMatrix[i][j] = overlap(kmers.get(i).getMu(), kmers.get(i).getSigma(), kmers.get(j).getMu(), kmers.get(j).getSigma());
				else if(j == i)
					overlapMatrix[i][j] = 1.0;
				else
					overlapMatrix[i][j] = overlapMatrix[j][i];
			}
			if(i % 16 == 0)
				System.out.println(kmers.get(i).getBaseChain());
		}
	}
	
	public void setUpTransitionMatrix() {
		transitionMatrix = new double[kmers.size()][kmers.size()];
		for(int i = 0; i < kmers.size(); i++) {
			for(int j = 0; j < kmers.size(); j++) {
				transitionMatrix[i][j] = transition(kmers.get(i), kmers.get(j));
			}
		}
	}
	
	public double overlap(double mu1, double sg1, double mu2, double sg2) {
		TwoGauss ol = new TwoGauss(mu1, mu2, sg1, sg2);
		double overlapRange = 4.0;
		int overlapPrecision = 5000;
		double a = mu1 < mu2 ? mu1 - overlapRange * sg1 : mu2 - overlapRange * sg2;
		double b = mu1 > mu2 ? mu1 + overlapRange * sg1 : mu2 + overlapRange * sg2;
		return integrator.trapezium(a, b, overlapPrecision, ol);
	}
	
	public double overlapFromMatrix(int pos1, int pos2) {
		return overlapMatrix[pos1][pos2];
	}
	
	public double emission(Kmer kmer, double mu, double sigma) {
		return overlap(kmer.getMu(), kmer.getSigma(), mu, sigma);
	}
	
	public double transition(Kmer start, Kmer end) {
		if(end.getBaseChain().equals(start.getBaseChain()))
			return k;
		else {
			if(start.hasTransitionTo(end))
				return (1.0 - k) / (double)(start.getTransitions().size() - 1);
		}
		return 0.0000000000000001;
	}
	
	public double transitionFromMatrix(int pos1, int pos2) {
		return transitionMatrix[pos1][pos2];
	}
	
	public void findSimilarPaths(int[] emissions){
		//System.out.println("Calculating probabilities...\nTick (total: " + emissions.length + "): ");
		State[][] result = new State[emissions.length][kmers.size()];
		for(int x = 0; x < emissions.length; x++) {
			//System.out.print(x + ",");
			//if(x % 50 == 0)
				//System.out.print("\n");
			for(int y = 0; y < kmers.size(); y++) {
				if(x == 0) {
					if(y == emissions[0])
						result[x][y] = new State(1.0, -1);
					else
						result[x][y] = new State(0.0, -1);
				}
				else {
					if(x == 1) {
						double Phere = transitionFromMatrix(getPosByBaseChain(start.getBaseChain()), y) * overlapFromMatrix(emissions[1], y) / 2.0;
						Phere = Math.log(Phere);
						result[x][y] = new State(Phere, emissions[0]);
					}
					else{
						State s = new State(0.0, -1);
						for(int i = 0; i < kmers.size(); i++) {
							double Phere = result[x-1][i].P() + Math.log(transitionFromMatrix(i, y) * overlapFromMatrix(emissions[x], y) / 2.0);
							if(i == 0 || Phere > s.P()) {
								s.setP(Phere);
								s.setPrev(i);
							}
						}
						result[x][y] = s;
					}
				}
			}
		}
		//System.out.print("\ngetting the best matching paths...\n");
		
		List<State> best = new ArrayList<>();
		for(int i = 0; i < kmers.size(); i++) {
			double Phere = result[emissions.length-2][i].P() + Math.log(transitionFromMatrix(i, getPosByBaseChain(end.getBaseChain()))) / 2.0;
			best.add(new State(Phere, i));
			//best.add(result[emissions.length-1][i]);
			if(i >= 3) {
				Collections.sort(best, new StateSorter());
				best.remove(best.size()-1);
			}
		}
		
		for(int i = 0; i < 3; i++) {
			int[] path = new int[emissions.length];
			State tmp = best.get(i);
			path[0] = emissions[0];
			path[emissions.length-1] = emissions[emissions.length-1];
			for(int x = emissions.length - 2; x > 0; x--) {
				path[x] = tmp.prev();
				tmp = result[x][tmp.prev()];
			}
			//System.out.print("path number " + i + " (log probability: " + best.get(i).P() + "):\n");
			//printPath(path);
			printPathToFile(path, i+1, best.get(i).P());
		}
		
		//System.out.println("finished");
	}
	
	public double similarity(Kmer start1, Kmer start2, Kmer end1, Kmer end2) {
		if(start1.hasTransitionTo(end1) && start2.hasTransitionTo(end2)) {
			return overlapFromMatrix(getPosByBaseChain(end1.getBaseChain()), getPosByBaseChain(end2.getBaseChain())) / 2.0;
		}
		return -1;
	}
	
	public double max(double[] numbers) {
		double res = 0;
		for(int i = 0; i < 5; i++) {
			if(numbers[i] > res)
				res = numbers[i];
		}
		return res;
	}
	
	public Kmer getKmerByBaseChain(String bc) {
		return kmers.get(getPosByBaseChain(bc));
	}
	
	public int getPosByBaseChain(String bc) {
		int pos = 0;
		for(int i = 0; i < bc.length(); i++) {
			if(bc.charAt(i) == 'C') {
				pos += Math.pow(4, 5-i);
			}
			else if(bc.charAt(i) == 'G') {
				pos += Math.pow(4, 5-i) * 2;
			}
			else if(bc.charAt(i) == 'T') {
				pos += Math.pow(4, 5-i) * 3;
			}
		}
		return pos;
	}
	
	public void generateRandomPath() {
		Random rand = new Random();
		int tick = rand.nextInt(20) + 40;
		int startpos = rand.nextInt(4096);
		
		start = kmers.get(startpos);
		
		ref = new int[tick];
		ref[0] = startpos;
		int hereSince = 1;
		for(int i = 1; i < tick; i++) {
			ref[i] = randomNextKmer(kmers.get(ref[i-1]), hereSince);
			if(ref[i] == ref[i-1])
				hereSince++;
			else
				hereSince = 1;
		}
		end = kmers.get(ref[tick-1]);
	}
	
	public void printPath(int[] path) {
		Kmer tmp = start;
		int counter = 1;
		System.out.print(tmp.getBaseChain() + "(" + tmp.getMu() + ", " + tmp.getSigma() + ")");
		for(int i = 1; i < path.length; i++) {
			if(!kmers.get(path[i]).getBaseChain().equals(tmp.getBaseChain())) {
				tmp = kmers.get(path[i]);
				System.out.print(" -> " + tmp.getBaseChain() + "(" + tmp.getMu() + ", " + tmp.getSigma() + ")");
				counter++;
				if(counter == 6) {
					counter = 0;
					System.out.print("\n\n");
				}
			}
		}
		System.out.print("\n\n");
	}
	
	public void printPathToFile(int[] path, int pathNum, double logProb) {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter("C:\\Users\\Vili\\Documents\\Egyetem\\Témalabor\\results.txt", true));
			
			Kmer tmp = start;
			int counter = 1;
			if(pathNum == 0) {
				writer.append("Reference sequence:\n");
			}
			else {
				writer.append("Path number " + pathNum + " (log probability: " + logProb + "):\n");
			}
			writer.append(tmp.getBaseChain() + "(" + tmp.getMu() + ", " + tmp.getSigma() + ")");
			for(int i = 1; i < path.length; i++) {
				if(!kmers.get(path[i]).getBaseChain().equals(tmp.getBaseChain())) {
					tmp = kmers.get(path[i]);
					writer.append(" -> " + tmp.getBaseChain() + "(" + tmp.getMu() + ", " + tmp.getSigma() + ")");
					counter++;
					if(counter == 6) {
						counter = 0;
						writer.append("\n\n");
					}
				}
			}
			writer.append("\n\n");
			if(pathNum == 3)
				writer.append("------------------------------------\n\n");
			
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private int randomNextKmer(Kmer kmer, int hereSince) {
		Random rand = new Random();
		double x1 = rand.nextDouble();
		String nextbc;
		if(x1 > Erlang(hereSince)) {
			nextbc = kmer.getBaseChain();
		}
		else {
			int x2 = rand.nextInt(4);
			nextbc = kmer.getTransitions().get(x2);
		}
		
		return getPosByBaseChain(nextbc);
	}

	private double Erlang(int hereSince) {
		double x = (double) hereSince;
		return 1.0 - Math.exp(-x/2.0) - Math.exp(-x/2.0) * x / 2.0;
	}
}
