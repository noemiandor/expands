package core;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import core.utils.AFtimesCopyNumber;
import core.utils.Common;
import core.utils.CopyNumber;
import core.utils.NotAValidCompositionException;

/**
 * Class modeling all possible subpopulation scenarios that can exist for a given mutated locus. Only one of these scenarios can account for the allele frequency and copy number measured at the locus.
 * For each scenario, model starts with a germline population that will be the root of all other modeled subpopulations. 
 * First subpopulation (SP1) modeled to evolve from the germline population is always the one carrying a CNV. 
 * Subsequent subpopulations are always defined by an SNV and are modeled in relation to SP1, either as:
 * - its parent
 * - its child (including twin - i.e. the two SPs are de-facto the same)
 * - its sibling
 * @author noemi
 *
 */
public class ParallelSubpopulations {
	/**
	 * Subpopulation scenarios where the subpopulation carrying the point mutation gave rise to the subpopulation carrying the CNV.
	 * Child acquiring the CNV is of interest too.
	 */
	private LinkedHashSet<Subpopulation> snvBeforeCNV;

	/**
	 * Subpopulation scenarios where the subpopulation carrying the point mutation evolved from the subpopulation carrying the CNV.
	 * Parent acquiring the CNV is of interest too.
	 */
	private LinkedHashSet<Subpopulation> snvAfterCNV;

	/**
	 * Subpopulation scenarios where the subpopulation carrying the point mutation evolved along with the subpopulation carrying the CNV from a common ancestor.
	 * Other sibling acquiring the CNV is of interest too.
	 */
	private LinkedHashSet<Subpopulation> snvSiblingOfCNV;


	/**
	 * 
	 * @param cn
	 * @param af
	 * @param pnb
	 * @param enforceCoocurrence - if true, the only scenarios modeled will be the ones where the point mutation and the CNV were propagated during the same clonal expansion.If set to false, all scenarios will be modeled.
	 * @throws NotAValidCompositionException
	 */
	public ParallelSubpopulations(double cn, double af, int pnb, boolean enforceCoocurrence) throws NotAValidCompositionException {
		snvBeforeCNV = new LinkedHashSet<Subpopulation>();
		snvAfterCNV = new LinkedHashSet<Subpopulation>();
		snvSiblingOfCNV = new LinkedHashSet<Subpopulation>();

		Subpopulation root = new Subpopulation(1, 2, pnb, null);
		HashSet<Subpopulation> sps = root.getPotentialChildren(new CopyNumber(
				cn)); //First subpopulation modeled is always the one carrying a CNV
		Iterator<Subpopulation> iter = sps.iterator();
		while (iter.hasNext()) {
			Subpopulation sp = iter.next();
			AFtimesCopyNumber afxcn = new AFtimesCopyNumber(af * cn);
			if(enforceCoocurrence || root.getPmb()==1){ //Co-occurrence always enforced when we're dealing with LOH <=> no more than 1 event per aberration-type
				Common.setALLOWED_FREQUENCIES(new double[]{sp.getCellularFrequency()});
			} else{
				snvBeforeCNV.addAll(sp.getPotentialParents(afxcn)); //there is a parental SP with a point mutation that gave rise to this SP
				snvSiblingOfCNV.addAll(sp.getPotentialSiblings(afxcn));//the two SPs are evolved from a common ancestor rather than from each other  <=> SNV and CNV are not tight
			}
			snvAfterCNV.addAll(sp.getPotentialChildren(afxcn));//this SP gave rise to a SP with a point mutation
		}

	}


	//	/**
	//	 * @return - only one frequency per SP family (parent/child/sibling) instead of all frequencies, while using pre-calculated overall-probability; 
	//	 */
	//	public Map<Double,Double> getCellFreq2ProbabilityMap() {
	//		Map<Subpopulation, Double> m=getSubpopulation2ProbabilityMap();
	//		Map<Double,Double> d = new HashMap<Double,Double>();
	//		for(Entry<Subpopulation, Double> e: m.entrySet()){
	//			Subpopulation sp = e.getKey();
	//			double value=e.getValue();
	//			if(d.containsKey(sp.getCellularFrequency())){
	//				value+=d.get(sp.getCellularFrequency());
	//			}
	//			d.put(sp.getCellularFrequency(), value);
	//		}
	//
	//		return (d);
	//	}

	/**
	 * The better a cellular frequency can account for the measured allele frequency and copy number, the higher the likelihood that it actually exists.
	 * The fewer somatic events are necessary to explain a subpopulation's evolution from the germline population (root), the higher the likelihood that it actually exists.
	 * @return a map of each cellular frequency to the likelihood that it manifests as reality. 
	 */
	public Map<Double,Double> getCellFreq2ProbabilityMap() {
		Map<Subpopulation, Double> m=getSubpopulation2ProbabilityMap();
		Map<Double,Double> d = new HashMap<Double,Double>();
		for(Subpopulation tsp: m.keySet()){
			List<Subpopulation> sps=new ArrayList<>();
			sps.add(tsp);
			if(tsp.getChild()!=null){
				sps.add(tsp.getChild());
			}
			if(tsp.getSibling()!=null){
				sps.add(tsp.getSibling());
			}
			if(!tsp.getParent().equals(tsp.getRoot())){
				sps.add(tsp.getParent());
			}

			for(Subpopulation sp:sps){
				//				double value=1/(sp.getDeviation()*sp.getRoot().eventsTo(sp));
				//increasing weight of somatic events count minimization: 
				double value=1/(sp.getDeviation()*Math.exp(1.0*sp.getRoot().eventsTo(sp)));
				if(d.containsKey(sp.getCellularFrequency())){
					value+=d.get(sp.getCellularFrequency());
				}
				d.put(sp.getCellularFrequency(), value);
			}
		}

		return (d);
	}

	/**
	 * 
	 * @return the subpopulation that best explains the measured allele frequency and copy number.
	 */
	public Entry<Subpopulation, Double> getBestFittingSubpopulation(){
		Map<Subpopulation,Double> all=getSubpopulation2ProbabilityMap();
		double pmax=Double.MIN_VALUE;
		Entry<Subpopulation, Double> sp=null;
		for(Entry<Subpopulation, Double> e: all.entrySet()){
			if(pmax<=e.getValue()){
				pmax=e.getValue();
				sp=e;
				//				System.out.println(sp.toString()+":"+pmax);
			}
		}
		return(sp);
	}


	/**
	 * The better a subpopulation can account for the measured allele frequency and copy number, the higher the likelihood that it actually exists.
	 * The fewer somatic events are necessary to explain a subpopulation's evolution from the germline population (root), the higher the likelihood that it actually exists.
	 * @TODO: Can we simplify to count # of events from root rather than parent/child, regardless of evolutionary scenario?
	 * @return a map of each subpopulation scenario to the likelihood that it manifests as reality. 
	 */
	public Map<Subpopulation,Double> getSubpopulation2ProbabilityMap() {
		Map<Subpopulation,Double> d = new HashMap<Subpopulation,Double>();
		Iterator<Subpopulation> iter = snvBeforeCNV.iterator();
		// Child acquiring the CNV is of interest too
		while (iter.hasNext()) {
			Subpopulation sp = iter.next();
			double value=1/(sp.getDeviation()  	* sp.getRoot().eventsTo(sp)+
					sp.getChild().getDeviation()* sp.eventsTo(sp.getChild()));
			if(d.containsKey(sp)){
				value+=d.get(sp);
			}
			d.put(sp, value);
		}

		// Parent acquiring the CNV is of interest too
		iter = snvAfterCNV.iterator();
		while (iter.hasNext()) {
			Subpopulation sp = iter.next();
			double x=sp.getParent().getDeviation() * sp.getRoot().eventsTo(sp.getParent());
			double value;
			if(sp.getParent().eventsTo(sp)==0){  //the two SPs are de-facto the same
				value=1/(x+ sp.getDeviation() * sp.getRoot().eventsTo(sp));
			}else{ //parent gave rise to child - not the same SPs
				value=1/(x+ sp.getDeviation() * sp.getParent().eventsTo(sp));
			}
			if(d.containsKey(sp)){
				value+=d.get(sp);
			}
			d.put(sp, value);
		}

		// Sibling acquiring the CNV is of interest too
		iter = snvSiblingOfCNV.iterator();
		while (iter.hasNext()) {
			Subpopulation sp = iter.next();
			double value=1/(sp.getDeviation()  	   * sp.getRoot().eventsTo(sp)+
					sp.getSibling().getDeviation() * sp.getRoot().eventsTo(sp.getSibling()));
			if(d.containsKey(sp)){
				value+=d.get(sp);
			}
			d.put(sp, value);
		}
		return (d);
	}


}
