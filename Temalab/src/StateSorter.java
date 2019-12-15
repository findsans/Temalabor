import java.util.Comparator;

public class StateSorter implements Comparator<State> {

	@Override
	public int compare(State arg0, State arg1) {
		if(arg1.P()-arg0.P() > 0)
			return 1;
		else if(arg1.P()-arg0.P() < 0)
			return -1;
		else
			return 0;
	}

}
