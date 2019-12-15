import java.util.ArrayList;
import java.util.List;

public class Kmer {
	private String basechain;
	private double sigma;
	private double mu;
	private List<String> transitions;
	
	public Kmer(String bc, double sg, double m) {
		basechain = bc;
		sigma = sg;
		mu = m;
		transitions = new ArrayList<>();
	}
	
	public String getBaseChain() {
		return basechain;
	}
	
	public double getMu() {
		return mu;
	}
	
	public double getSigma() {
		return sigma;
	}
	
	public List<String> getTransitions(){
		return transitions;
	}
	
	public boolean hasTransitionTo(Kmer end) {
		for(int i = 0; i < transitions.size(); i++) {
			if(transitions.get(i).equals(end.getBaseChain()))
				return true;
		}
		return false;
	}
	
	public void setUpTransitions() {
		transitions.add(basechain.substring(1, 6) + "A");
		transitions.add(basechain.substring(1, 6) + "G");
		transitions.add(basechain.substring(1, 6) + "T");
		transitions.add(basechain.substring(1, 6) + "C");
		if(!(basechain.equals("AAAAAA")
				|| basechain.equals("GGGGGG")
				|| basechain.equals("CCCCCC")
				|| basechain.equals("TTTTTT")))
			transitions.add(basechain);
	}
	
	public void print() {
		System.out.println("Kmer: " + basechain);
		System.out.println("mu = " + mu + ", sigma = " + sigma);
		System.out.println("possible next Kmers:");
		for(int i = 0; i < getTransitions().size(); i++) {
			System.out.println(getTransitions().get(i));
		}
	}
}
