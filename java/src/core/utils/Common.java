package core.utils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import core.ParallelSubpopulations;
import core.Subpopulation;

public final class Common {

	public static double rESOLUTION = 0.01;


	public static int MAX_PM = 6;
	
	/**
	 * Holds the allowed frequency spectrum for subpopulations which evolved through point mutations
	 */
	public static Set<Double> allowed_SP_FREQUENCIES;
	
	/**
	 * Holds the allowed frequency spectrum for subpopulations which evolved through copy number variations, if different than for point mutations
	 */
	private static Set<Double> allowed_SP_CNV_FREQUENCIES=null;
	
	/**
	 * Holds the allowed frequency spectrum of sibling subpopulations mapped to the siblings subpopulation size.
	 */
	private static Map<Double,Set<Double>> allowedSiblings_SP_FREQUENCIES=null;

	static{
		allowed_SP_FREQUENCIES=seq(0,0.99,rESOLUTION);
	}




	public static Set<Double> seq(double min, double max, double step) {
		Set<Double> d = new HashSet<Double>();
		double v = min-step;
		while(v<max){
			v+=step;
			d.add(v);
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

	/**
	 * Resetting general SP frequency spectrum also resets @allowedSiblings_SP_FREQUENCIES
	 * @param aLLOWED_FREQUENCIES
	 */
	public static void setALLOWED_SP_FREQUENCIES(double[] aLLOWED_FREQUENCIES) {
		allowed_SP_FREQUENCIES = new HashSet<Double>();
		for(int i =0; i<aLLOWED_FREQUENCIES.length; i++){
			allowed_SP_FREQUENCIES.add(aLLOWED_FREQUENCIES[i]);
		}
		allowedSiblings_SP_FREQUENCIES=null;
		allowed_SP_CNV_FREQUENCIES=null;

	}
	
	public static void setALLOWED_SP_CNV_FREQUENCIES(double[] aLLOWED_FREQUENCIES) {
		allowed_SP_CNV_FREQUENCIES = new HashSet<Double>();
		for(int i =0; i<aLLOWED_FREQUENCIES.length; i++){
			allowed_SP_CNV_FREQUENCIES.add(aLLOWED_FREQUENCIES[i]);
		}

	}
	
	public static void setAllowedSiblings_SP_FREQUENCIES(double sp, double sibling) {
		if(allowedSiblings_SP_FREQUENCIES==null){
		allowedSiblings_SP_FREQUENCIES= new HashMap<Double,Set<Double>>();
		}
		if(!allowedSiblings_SP_FREQUENCIES.containsKey(sp)){
			allowedSiblings_SP_FREQUENCIES.put(sp, new HashSet<Double>());
		}
		allowedSiblings_SP_FREQUENCIES.get(sp).add(sibling);
	}
	
	
	/**
	 * Changes default in favor of parental over sibling relations.
	 * This will allow only parental relations among SPs; 
	 * unless siblings are explicitly noted for a given subpopulation. 
	 */
	public static void setAllowedSiblings_SP_FREQUENCIES() {
		if(allowedSiblings_SP_FREQUENCIES==null){
			allowedSiblings_SP_FREQUENCIES= new HashMap<Double,Set<Double>>();
		}
	}
	
	
	public static Set<Double> getallowed_SP_FREQUENCIES() {
		return allowed_SP_FREQUENCIES;
	}
	
	/**
	 * Per default the allowed subpopulation frequency spectrum is identical for point mutations and copy number variations
	 * @return the allowed frequency spectrum for subpopulations which evolved through copy number variations
	 */
	public static Set<Double> getAllowed_SP_CNV_FREQUENCIES() {
		if(allowed_SP_CNV_FREQUENCIES==null || allowed_SP_CNV_FREQUENCIES.size()==0){
			return allowed_SP_FREQUENCIES;
		}
		return allowed_SP_CNV_FREQUENCIES;
	}
	
	/**
	 * Per default the allowed subpopulation frequency spectrum is identical for all phylogenies.
	 * This function allows the frequency spectrum of siblings to be modified for every subpopulation:
	 * @param sp
	 * from which a sibling can diverge.
	 * >> Note: as soon as function has been called for at least one subpopulation, << 
	 * >> default behavior will no longer apply for any subpopulation. << 
	 * @return frequency spectrum of all siblings allowed to diverge from given subpopulation.
	 */
	public static Set<Double> getAllowedSiblings_SP_FREQUENCIES(Double sp) {
		if(allowedSiblings_SP_FREQUENCIES==null){
			return allowed_SP_FREQUENCIES;
		}else if(!allowedSiblings_SP_FREQUENCIES.containsKey(sp)){
			return new HashSet<Double>();
		}
		return allowedSiblings_SP_FREQUENCIES.get(sp);
	}
	
	/**
	 * 
	 * @param sp
	 * @return all subpopulations that don't have a fraternal relationship with this sp
	 */
	public static Set<Double> getAllowedPaternal_SP_FREQUENCIES(Double sp) {
		if(allowedSiblings_SP_FREQUENCIES==null){
			return allowed_SP_FREQUENCIES;
		}
		Set<Double> v=new HashSet<Double>();
		v.addAll(allowed_SP_FREQUENCIES);
		v.removeAll(getAllowedSiblings_SP_FREQUENCIES(sp));
		return(v);
	}
	
	
	public static void setMAX_PM(int max_PM) {
		MAX_PM=max_PM;
	}

//	public static void setRESOLUTION(double r) {
//		rESOLUTION = r;
//	}

	public static Double[] below(double f, Set<Double> a) {
		List<Double> l =new ArrayList<Double>();
		for(Double d: a){
			if(d<=f){
				l.add(d);
			}
		}
		return l.toArray(new Double[l.size()]);
	}

	public static Double[] above(double f, Set<Double> a) {
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
			Common.setMAX_PM(6);
			Common.setALLOWED_SP_FREQUENCIES(new double[]{   0.343, 0.494 ,0.737, 0.646, 0.859 ,0.919});
			Common.setAllowedSiblings_SP_FREQUENCIES();
			ParallelSubpopulations sps = new ParallelSubpopulations(  2.317312, 0.852459, 1,false,2);
			//PM=2, PM_cnv=2, PM_B=1, SP=0.23
			Map<Double, Double> d = sps.getCellFreq2ProbabilityMap();
			Map<Subpopulation, Double> ds = sps.getSubpopulation2ProbabilityMap();
			Entry<Subpopulation, Double> sp = sps.getBestFittingSubpopulation();
			System.out.println(sp);
			if(sp.getKey().getChild()!=null){
				System.out.println("Child--> "+sp.getKey().getChild());
			}else {
				System.out.println("Parent--> "+sp.getKey().getParent());
			}
		} catch (NotAValidCompositionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}


	}






	
}
