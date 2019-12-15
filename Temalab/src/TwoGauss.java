
public class TwoGauss implements FPFunction {
	private double mu1;
	private double mu2;
	private double sigma1;
	private double sigma2;
	
	public TwoGauss(double m1, double m2, double sg1, double sg2) {
		mu1 = m1;
		mu2 = m2;
		sigma1 = sg1;
		sigma2 = sg2;
	}

	@Override
	public double eval(double n) {
		double a1 = 1.0 / (sigma1 * Math.sqrt(2.0 * Math.PI));
		double b1 = -1.0 * (n - mu1) * (n - mu1) / (2 * sigma1 * sigma1);
		double res1 = a1 * Math.exp(b1);
		
		double a2 = 1.0 / (sigma2 * Math.sqrt(2.0 * Math.PI));
		double b2 = -1.0 * (n - mu2) * (n - mu2) / (2 * sigma2 * sigma2);
		double res2 = a2 * Math.exp(b2);
		
		if(res1 < res2)
			return res1;
		else
			return res2;
	}
	
}
