
public class Integrator {
	
	public double trapezium(double a, double b, int n, FPFunction f) {
		double range = checkParamsGetRange(a, b, n);
		double nFloat = (double)n;
		double sum = 0.0;
		
		for(int i = 0; i < n; i++) {
			double x = a + range * (double)i / nFloat;
			sum += f.eval(x);
		}
		sum += (f.eval(a) + f.eval(b)) / 2.0;
		
		return sum * range / nFloat;
	}
	
	public double checkParamsGetRange(double a, double b, int n) {
		if(n <= 0)
			throw new IllegalArgumentException("Invalid value of n");
		double range = b-a;
		if(range <= 0)
			throw new IllegalArgumentException("Invalid range");
		
		return range;
	}
}
