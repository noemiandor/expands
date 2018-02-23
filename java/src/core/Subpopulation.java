package core; 
import java.util.HashSet;

import core.utils.AFtimesCopyNumber;
import core.utils.Common;
import core.utils.CopyNumber;
import core.utils.NotAValidCompositionException;

/**
 * Class modeling a virtual subpopulation at a single mutated genomic locus. The locus is characterized by: 
 * - a somatic point mutation (has to be present) 
 * - a copy number aberration (may or may not be present).
 * These two characteristics of a mutated locus enforce every somatic subpopulation to have exactly one relative of either:
 * - a child
 * - a sibling
 * - a parent that is not the root,
 * whereby the subpopulation itself represents one characteristic, while its relative represents the other.
 * Any subpopulation that does not have any of above mentioned relatives is the root itself, and is thus germline, not somatic.
 * 
 * @author noemi
 *
 */
public class Subpopulation {
	/**
	 * The total ploidy of the mutated locus
	 */
	private int pm;

	/**
	 * The ploidy of the mutated allele
	 */
	private int pmb;

	/**
	 * The cellular frequency of the subpopulation
	 */
	private double f;

	/**
	 * The subpopulation from which this subpopulation evolved
	 */
	private Subpopulation parent;

	/**
	 * The subpopulation to which this subpopulation gave rise 
	 */
	private Subpopulation child;

	/**
	 * The subpopulation that evolved along with this subpopulation from a common ancestor.
	 */
	private Subpopulation sibling;


	/**
	 * Difference between the actual bulk population measurements and the ploidy constellation of this subpopulation (including all its relatives). 
	 * The higher the value of this variable, the higher the likelihood that the virtual subpopulation does not exist in reality. 
	 */
	private double deviation;



	public Subpopulation(double f, int pm, int pmb, Subpopulation parent) throws NotAValidCompositionException {
		if (f < 0 || f > 1) {
			throw new IllegalArgumentException(
					"A cellular frequency must be between 0 and 1.");
		}
		if (pm < 0) {
			throw new IllegalArgumentException(
					"Subpopulation's ploidy must be >= 0.");
		}
		if (pmb < 0) {
			throw new IllegalArgumentException(
					"Subpopulation's mutated allele ploidy must be >= 0.");
		}
		if (pmb > pm) {
			throw new IllegalArgumentException(
					"Subpopulation's mutated allele ploidy can't exceed its total ploidy.");
		}
		this.parent=parent;
		this.f = f;
		this.pm = pm;
		this.pmb = pmb;
		if(parent!=null && !isValidAsSomatic()){ //The root population is the only population that doesn't have a parent and that is allowed to be germline
			throw new NotAValidCompositionException();
			//Not a valid somatic subpopulation - this looks exactly like germline.
		}
	}

	/**
	 * Verifies the constraints that define a somatic subpopulation: either over-representation of mutated allele or presence of copy number aberration
	 * @return true if this is a valid somatic subpopulation, false otherwise 
	 */
	private boolean isValidAsSomatic() {
		boolean c2=this.pm!=this.getRoot().pm;
		if (c2){
			return true;
		}
		boolean c1=this.pmb/(double)this.pm > this.getRoot().pmb/(double)this.getRoot().pm;
		return c1;
	}

	/**
	 * 
	 * @return the germline population from which this subpopulation descends.
	 */
	public Subpopulation getRoot() {
		Subpopulation root = this;
		while (root.parent != null) {
			root = root.parent;
		}
		return (root);
	}



	//	/**
	//	 * >>>For LOH only<<<
	//	 * Modeling subpopulations carrying the germline SNP and a CNV that can have given rise to this subpopulation
	//	 * through a 2nd copy number variation. Any parent must have a cellular frequency above
	//	 * that of this subpopulation: @param f. 
	//	 * Parents' pm constrains its pmb and must be different from the root's ploidy (normally 2) and from this subpopulation's copy number (i.e. parent does NOT have this subpopulation's cnv). 
	//	 * Parents' pmb is >=1 and depends on @param afxcn.
	//	 * This subpopulation cannot have a parent if it has the same ploidy as the root (normally pm=2), because then a CNV cannot explain its evolution;
	//	 * @param af
	//	 * @return
	//	 * @throws NotAValidCompositionException 
	//	 */
	//	public HashSet<Subpopulation> getPotentialParentsWithCNV(AFtimesCopyNumber afxcn) {
	//
	//		int pnb=getRoot().getPmb();
	//		HashSet<Subpopulation> sps = new HashSet<Subpopulation>();
	//		if(this.pm==this.getRoot().pm || pnb==0){
	//			return(sps);
	//		}
	//		Double[] freqs=Common.above(f, Common.getAllowedPaternal_SP_FREQUENCIES(this.f));
	//		for(double pf: freqs) {
	//			for (double pm : Common.seq(pnb	, Common.MAX_PM, 1)) {
	//				if(pm==this.getRoot().getPm() || pm==this.pm){
	//					continue;
	//				}
	//				for (double pmb : Common.seq(pnb, pm, 1)) {
	//					try{
	//						Subpopulation s = new Subpopulation(pf, (int) pm, (int) pmb,this.getRoot());
	//						s.setChild(this);
	//						s.deviation = afxcn.getValue()
	//								- (s.pmb * (s.f-s.child.f) + s.child.pmb*s.child.f + this.getRoot().pmb * (1 - s.f));
	//						sps.add(s);
	//					}catch(NotAValidCompositionException e){
	//						continue;
	//					}
	//
	//				}
	//
	//			}
	//		}
	//		return sps;
	//	}

	/**
	 * Modeling subpopulations carrying a point mutation that can have given rise to this subpopulation
	 * through a copy number variation. Any parent must have a cellular frequency above
	 * that of this subpopulation: @param f. 
	 * Parents' pm must be equal to root's pm (i.e. parent does NOT have this subpopulation's cnv). 
	 * Parents' pmb depends on @param afxcn and is constrained by its pm.
	 * This subpopulation cannot have a parent if it has the same ploidy as the root (normally pm=2), because then a CNV cannot explain its evolution;
	 * @param af
	 * @return
	 * @throws NotAValidCompositionException 
	 */
	public HashSet<Subpopulation> getPotentialParents(AFtimesCopyNumber afxcn) {

		HashSet<Subpopulation> sps = new HashSet<Subpopulation>();
		if(this.pm==this.getRoot().pm){
			return(sps);
		}
		Double[] freqs=Common.above(f, Common.getAllowedPaternal_SP_FREQUENCIES(this.f));
		int PM=this.getRoot().getPm();
		for(double pf: freqs) {
			for (double pmb : Common.seq(1, PM, 1)) {
				try{
					Subpopulation s = new Subpopulation(pf, PM, (int) pmb,this.getRoot());
					s.setChild(this);
					s.deviation = afxcn.getValue()
							- (s.pmb * (s.f-s.child.f) + s.child.pmb*s.child.f + this.getRoot().pmb * (1 - s.f));
					sps.add(s);
				}catch(NotAValidCompositionException e){
					continue;
				}

			}

		}
		return sps;
	}

	/**
	 * Modeling subpopulation that evolved along with this subpopulation from a common ancestor. 
	 * Sibling diverged from this subpopulation through a point mutation. 
	 * Sibling may have any cellular frequency for as long as it doesn't exceed 1 when added to this subpopulation's frequency.
	 * Sibling's pm must be equal to root's pm.
	 * Sibling's pmb depends on @param afxcn and is constrained by its pm.
	 * This subpopulation cannot have a sibling if it already has the point mutation itself, because then the point mutation cannot explain the divergence.
	 * @param afxcn
	 * @return
	 */
	public HashSet<Subpopulation> getPotentialSiblings(
			AFtimesCopyNumber afxcn) {
		HashSet<Subpopulation> sps = new HashSet<Subpopulation>();
		if(this.pmb>0){
			return(sps);
		}

		Double[] freqs=Common.below(1-f, Common.getAllowedSiblings_SP_FREQUENCIES(this.f));
		int PM=this.getRoot().getPm();
		for(double pf: freqs) {
			for (double pmb : Common.seq(1, PM, 1)) {

				try{
					Subpopulation s = new Subpopulation(pf, PM, (int) pmb,this.getRoot());
					s.setSibling(this);
					s.deviation = afxcn.getValue()
							- (s.pmb * s.f + this.getRoot().pmb * (1 - s.f));
					sps.add(s);
				}catch(NotAValidCompositionException e){
					continue;
				}

			}

		}
		return sps;
	}



	/**
	 * Modeling the subpopulations that can evolve from this subpopulation
	 * through a point mutation. Any child must have a cellular frequency below
	 * that of this subpopulation: @param f.
	 * Child's pm must be equal to this subpopluation's pm. 
	 * Child's pmb depends on @param afxcn and is constraint by its pm.
	 * Child's pmb must be at least as high as this subpopluation's pmb.
	 * - This subpopulation cannot be a parent if it has already lost the germline variant (if any), 
	 * because then it cannot give rise to a subpopulation that again harbors that same germline variant.  
	 * - This subpopulation also cannot be a parent if it has already lost all copies of the locus, 
	 * because then it cannot give rise to a subpopulation that again harbors a somatic variant at that same locus.  
	 * - Finally, this subpopulation cannot be a parent if it has already acquired the somatic variant,
	 * because then the variant cannot explain the child's evolution (exception if child and parent are de-facto the same).
	 * @param afxcn
	 * @return
	 */
	public HashSet<Subpopulation> getPotentialChildren(AFtimesCopyNumber afxcn) {
		HashSet<Subpopulation> sps = new HashSet<Subpopulation>();
		if(this.pmb<this.getRoot().pmb || this.pm==0){ 
			return(sps);
		}

		Double[] freqs=Common.below(f, Common.getAllowedPaternal_SP_FREQUENCIES(this.f));
		if(this.getRoot().pmb==0 && this.pmb>0){
			freqs=new Double[]{f};//Set to only identical parent-child relation: de-facto same subpopulation
		}

		for(double pf: freqs) {
			for (double pmb : Common.seq(Math.max(1, this.pmb), pm, 1)) {
				try {
					Subpopulation s = new Subpopulation(pf, pm, (int) pmb, this.getRoot());
					s.setParent(this);
					s.deviation = afxcn.getValue()
							- (s.pmb * s.f + this.getRoot().pmb * (1 - s.f));
					sps.add(s);
				} catch (NotAValidCompositionException e) {
					continue;
				}

			}
		}
		return sps;
	}

	/**
	 * Modeling the subpopulations that can evolve from this subpopulation
	 * through a copy number variation. Any child must have a cellular frequency
	 * below that of this subpopulation: @param f.
	 * Child's pm depends on @param cn and constrains its pmb. Child's pm is not constraint by this subpopulation's pm.
	 * @TODO: Guarantee that this method will be called exclusively by a root population.
	 * @param cn
	 * @return
	 */
	public HashSet<Subpopulation> getPotentialChildren(CopyNumber cn) {

		HashSet<Subpopulation> sps = new HashSet<Subpopulation>();
		Double[] freqs=Common.below(f, Common.getAllowed_SP_CNV_FREQUENCIES());
		for(double pf: freqs) {
			for (double pm : Common.seq(0, Common.MAX_PM, 1)) {
				double dev = cn.getValue()
						- (pm * pf + this.getRoot().pm * (1 - pf));

				for(double pmb: Common.seq(0, pm, 1)){
					try {
						Subpopulation s1= new Subpopulation(pf, (int) pm, (int)pmb, this);
						s1.setParent(this);
//						s1.deviation = dev;
						s1.deviation = Math.signum(dev)*Math.pow(1+Math.abs(dev),2); //CN deviation has stronger weight than AF deviation
						sps.add(s1);
					} catch (NotAValidCompositionException e) {
					}
				}
			}

		}

		return sps;
	}

	@Override
	/**
	 * Two subpoulations are equal if and only if they have:
	 * - the same cellular frequency
	 * - the same ploidy constellations
	 * - the same exact relative (either one of parent, child or sibling)
	 */
	public boolean equals(Object o) {
		if(o instanceof Subpopulation){
			Subpopulation os = (Subpopulation)o;
			if(os.f==f && os.pm==pm && os.pmb==pmb){
				//parent
				if(this.parent!=null){
					if(!this.parent.equals(os.parent)){
						return false;
					}
				}else if(os.parent!=null){
					return false;
				}
				//child
				if(this.child!=null){
					if(!this.child.equals(os.child)){
						return false;
					}
				}else if(os.child!=null){
					return false;
				}
				return true;
			}
		}
		return false;
	}

	/**
	 * 
	 * @param os
	 * @return true if the two subpopulations have the same ploidy constellations, false otherwise
	 */
	public boolean shallowEquals(Subpopulation os) {
		if(os.pm==pm && os.pmb==pmb){
			return true;
		}
		return false;
	}

	/**
	 * Sets the parent of this subpopulation. Parent's cellular frequency has to be >= to this subpopulation's frequency.
	 * The parent's cellular frequency may be equal to this subpopulation's frequency if and only if, the parent also has the exact same ploidy constellation (i.e. is de-facto the same population).
	 * @TODO: NotAValidCompositionException if subpopulation already has a child or sibling 
	 * @param subpopulation
	 * @throws NotAValidCompositionException
	 */
	private void setParent(Subpopulation subpopulation) throws NotAValidCompositionException {
		if(subpopulation.f<this.f){
			throw new IllegalArgumentException("Parent can't be smaller than its child's subpopulation size");
		}else if((subpopulation.f==this.f && !this.shallowEquals(subpopulation)) || 
				(subpopulation.f!=this.f &&  this.shallowEquals(subpopulation)) ||
				(this.pm-this.pmb>0 && subpopulation.pm-subpopulation.pmb==0) ){//child, but not parent, has WT allele
			throw new NotAValidCompositionException();
			//"Identical subpopulations can't have divergent PM and PM_B"
		}
		this.parent = subpopulation;
		//@TODO: wouldn't subpopulation.setChild(this); work?
	}

	/**
	 * Sets the sibling of this subpopulation. 
	 * The sibling's cellular frequency may be equal to this subpopulation's frequency if and only if, the sibling also has the exact same ploidy constellation (i.e. is de-facto the same population).
	 * @TODO: NotAValidCompositionException if subpopulation already has a child or non-root parent
	 * @param subpopulation
	 * @throws NotAValidCompositionException
	 */
	private void setSibling(Subpopulation subpopulation) throws NotAValidCompositionException {
		if(subpopulation.f+this.f>1){
			throw new IllegalArgumentException("Two siblings' cummulative cellular frequency cannot exceed 1");
		}else if((subpopulation.f==this.f && !this.shallowEquals(subpopulation)) || 
				(subpopulation.f!=this.f &&  this.shallowEquals(subpopulation)) ){
			throw new NotAValidCompositionException();
			//"Identical subpopulations can't have divergent PM and PM_B"
		}
		this.sibling = subpopulation;
	}

	/**
	 * Sets the child of this subpopulation. Child's cellular frequency has to be <= to this subpopulation's frequency.
	 * The child's cellular frequency may be equal to this subpopulation's frequency if and only if, the child also has the exact same ploidy constellation (i.e. is de-facto the same population).
	 * @TODO: NotAValidCompositionException if subpopulation already has a non-root parent or sibling
	 * @param subpopulation
	 * @throws NotAValidCompositionException
	 */
	private void setChild(Subpopulation subpopulation) throws NotAValidCompositionException {
		if(subpopulation.f>this.f){
			throw new IllegalArgumentException("Child can't be larger than its parent subpopulation size");
		}else if((subpopulation.f==this.f && !this.shallowEquals(subpopulation)) || 
				(subpopulation.f!=this.f &&  this.shallowEquals(subpopulation)) ||
				(subpopulation.pm-subpopulation.pmb>0 && this.pm-this.pmb==0)){ //child, but not parent, has WT allele
			throw new NotAValidCompositionException();
			//"Identical subpopulations can't have divergent PM and PM_B"
		}
		this.child = subpopulation;
	}

	public double getCellularFrequency() {
		return f;
	}

	/**
	 * Guaranteed not to return zero, even for perfect agreement between subpopulation's constellation and measurement.
	 * @return the discrepancy between this subpopulation and the measured allele frequency/copy number. 
	 */
	public double getAbsoluteDeviation() {
		return Math.max(Double.MIN_NORMAL, Math.abs(deviation));
	}

	public double getDeviation() {
		return (deviation);
	}

	public Subpopulation getChild() {
		return child;
	}

	public Subpopulation getSibling() {
		return sibling;
	}

	public Subpopulation getParent() {
		return parent;
	}

	public int getPm() {
		return pm;
	}

	public int getPmb() {
		return pmb;
	}

	public String toString(){
		return f+":pm="+pm+";pmb="+pmb;
	}

	/**
	 * 
	 * @param sp - another subpopulation
	 * @return the number of somatic events that can explain the evolution of @param sp from this subpopulation
	 */
	public double eventsTo(Subpopulation sp) {
		double rootOriginSlowDown=1.0; //If one of the populations is germline (root), a mutation is less likely to occur --> polyclonal origin unlikely
		if(!this.isValidAsSomatic()	|| !sp.isValidAsSomatic()){
			rootOriginSlowDown=1.0; //One somatic event in a mutated cell is as likely as 10 somatic events in the root	
		}

		int mut_d=sp.pmb-pmb; //Difference in mutated allele ploidy
		int d=sp.pm-pm; //Difference in overall ploidy
		//Can difference in mutated allele ploidy be explained by gain/amplification?
		double v;
		//@TODO: Consider being slightly less stringent with deviation from background ploidy state, e.g. use "Math.pow(1.5,x)" instead of "Math.exp(x)"
		if(this.getPmb()>0 && mut_d==d){
			v= rootOriginSlowDown * Math.exp( Math.abs(mut_d) );
		}else{
			v=rootOriginSlowDown * Math.exp( Math.abs(d) + Math.abs(mut_d) );
		}
		return(v);
	}


}
