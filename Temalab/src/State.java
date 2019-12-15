
public class State {
	private double value;
	private int posBefore;
	
	public State(double p, int prev) {
		value = p;
		posBefore = prev;
	}
	
	public void setP(double p) {
		value = p;
	}
	public void setPrev(int prev) {
		posBefore = prev;
	}
	
	public double P() {
		return value;
	}
	
	public int prev() {
		return posBefore;
	}
}
