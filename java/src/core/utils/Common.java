package core.utils;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import core.ParallelSubpopulations;
import core.Subpopulation;

public final class Common {

	public static double rESOLUTION = 0.01;


	public static int MAX_PM = 6;
	public static Double[] allowed_FREQUENCIES;

	static{
		allowed_FREQUENCIES=seq(0,0.99,rESOLUTION);
	}




	public static Double[] seq(double min, double max, double step) {
		Double[] d = new Double[(int) ((step + max - min) / step)];
		for (int i = 0; i < d.length; i++) {
			d[i] = min + i * step;
		}
		return d;
	}



	public static <C extends Object> List<C> toList(C[] value) {
		List<C> l = new ArrayList<C>();
		for (C d : value) {
			l.add(d);
		}
		return l;
	}

//	public static Double[] intersect(Double[] a1, Double[] a2) {
//		HashSet<Double> intersection = new HashSet<Double>();
//		for (Double i : a1) {
//			for (Double j : a2) {
//				if (Math.abs(i-j)<=rESOLUTION) {
//					intersection.add(i);
//					break;
//				}
//			}
//		}
//		Double[] i = intersection.toArray(new Double[intersection.size()]);
//		Arrays.sort(i);
//		return i;
//	}

	public static void setALLOWED_FREQUENCIES(double[] aLLOWED_FREQUENCIES) {
		allowed_FREQUENCIES = new Double[aLLOWED_FREQUENCIES.length];
		for(int i =0; i<allowed_FREQUENCIES.length; i++){
			allowed_FREQUENCIES[i]=aLLOWED_FREQUENCIES[i];
		}

	}
	
	public static void setMAX_PM(int max_PM) {
		MAX_PM=max_PM;
	}

//	public static void setRESOLUTION(double r) {
//		rESOLUTION = r;
//	}

	public static Double[] below(double f, Double[] a) {
		List<Double> l =new ArrayList<Double>();
		for(Double d: a){
			if(d<=f){
				l.add(d);
			}
		}
		return l.toArray(new Double[l.size()]);
	}

	public static Double[] above(double f, Double[] a) {
		List<Double> l =new ArrayList<Double>();
		for(Double d: a){
			if(d>=f){
				l.add(d);
			}
		}
		return l.toArray(new Double[l.size()]);
	}

	
	public static void main(String[] args) {
		try {
//			Common.allowed_FREQUENCIES=new Double[]{0.2063427, 0.3126855, 0.6529823};
//			Common.rESOLUTION=0.02126855;
			ParallelSubpopulations sps = new ParallelSubpopulations(1.278829, 0.8969697, 1,true);
			Map<Double, Double> d = sps.getCellFreq2ProbabilityMap();
			Map<Subpopulation, Double> ds = sps.getSubpopulation2ProbabilityMap();
			System.out.println(sps.getBestFittingSubpopulation());
		} catch (NotAValidCompositionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}


	}



	
}
