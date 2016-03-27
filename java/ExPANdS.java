import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;

/**
 * Class used by the R-package EXPANDS (Expanding Ploidy and Allele Frequency on
 * Nested Subpopulations), that characterizes coexisting SPs in a tumor using
 * copy number and B-allele frequencies derived from exome- or whole genome
 * sequencing input data. Functionality provides exhaustive iteration through
 * scenarios that explain the B-allele frequency and copy number measured for a
 * given locus, while considering cell type specific ploidies, cellular
 * frequencies and noise in copy number measurements . Functionalities concern
 * individual point mutations only, rather than mutation clusters.
 * 
 * @author Noemi
 * 
 */
public class ExPANdS {

	/**
	 * Stores default value of maximal ploidy in mutated cells, accepted as
	 * solution.
	 */
	public static final int DEFAULT_MAX_PM = 6;

	private static final int SNVONLY_FLAG = 1;
	private static final int CNVONLY_FLAG = 2;
	private static final int SNVANDCNV_FLAG = 3;
	private static final int SNVBEFORECNV_FLAG = 4;

//	/**
//	 * Minimal ploidy in mutated cells (PM=0 is not accepted, since B-allele
//	 * ploidy equation not applicable)
//	 */
//	private static int min_PM = 1;;

	/**
	 * Total ploidy in normal cells
	 */
	private static int PN = 2;

	/**
	 * Alternative solutions, consisting of (PM, PM_B, f). Each entry in the
	 * list consists of fixed values for PM and PM_B and variable f (frequency
	 * of mutated cell)
	 */
	private List<FitAlternatives> solutions;

	/**
	 * Observed/measured copy number
	 */
	private double CN;

	/**
	 * Observed/measured B-allele frequency
	 */
	private double AF;

	/**
	 * Maximal value for estimated copy-number error (default 0.1)
	 */
	private static double maxCN_Error = 0.1;

	/**
	 * Observed/measured ploidy of B-allele in normal cells
	 */
	private int PN_B;

	/**
	 * The error associated with the currently best solution
	 */
	private double minErr;

	/**
	 * Copy number alternatives obtained by adding different levels of noise to
	 * the measured copy number: alternativeCN:=(CN - maxCN_Error):(CN +
	 * maxCN_Error), where CN is the measured copy number and maxCN_error is the
	 * maximal level of noise.
	 */
	private double[] alternativeCN;

	/**
	 * Cellular frequency alternatives for fixed PM and PM_B
	 */
	private FitAlternatives bestSolution;

	/**
	 * Maximal ploidy in mutated cells, accepted as solution.
	 */
	private int max_PM;

	/**
	 * The ploidy of the subpopulation harboring the CNV (set only if SNV
	 * occured before CNV)
	 */
	private int PM = Integer.MIN_VALUE;

	/**
	 * The size of the subpopulation harboring the CNV (set only if SNV occured
	 * before CNV)
	 */
	private double f_CNV = Double.NaN;

	/**
	 * Deviation between observed and expected values for B-allele frequency and
	 * copy number. The value of this function is 0, if and only if the
	 * components of the best solution yield exactly the same copy number and
	 * allele frequency as measured.
	 * 
	 * @return the deviation of the observed from the expected copy number and
	 *         allele frequency, for the inferred solution
	 */
	public double getDeviation() {
		return minErr;
	}

	/**
	 * Calculates alternative solutions, consisting of (PM, PM_B, f) =(total
	 * ploidy in mutated cells, B-allele ploidy in mutated cells, frequency of
	 * mutated cells), that can explain B-allele frequency and copy number
	 * measured for a single mutated locus.
	 * 
	 * @param af
	 *            measured B-allele frequency
	 * @param cn
	 *            measured copy number
	 * @param pnb
	 *            ploidy of the B-allel in normal cells
	 * @param max_PM
	 *            maximal total ploidy in mutated cells
	 */
	public ExPANdS(double af, double cn, int pnb, int max_PM) {
		this.max_PM = max_PM;
		this.PN_B = pnb;
		init(af, cn);
	}

	/**
	 * This constructor is to be invoked if the SNV happened in a clonal
	 * expansion that preceded the clonal expansion of a CNV affecting the same
	 * locus. Consequence is that PM_B can have >2 values in a subpopulation
	 * dependent manner.
	 * 
	 * Calculates alternative solutions, consisting of (PM, PM_B, f) =(total
	 * ploidy in mutated cells, B-allele ploidy in mutated cells, frequency of
	 * mutated cells), that can explain B-allele frequency and copy number
	 * measured for a single mutated locus.
	 * 
	 * @param af
	 *            measured B-allele frequency
	 * @param cn
	 *            measured copy number
	 * @param f_cnv
	 *            size of the subpopulation harboring the CNV
	 * @param PM
	 *            ploidy of the subpopulation harboring the CNV
	 */
	public ExPANdS(double af, double cn, double f_cnv, int PM, int pnb) {
		this.f_CNV = f_cnv;
		this.PM = PM;
		this.PN_B=pnb;
		init(af, cn);
	}

	private void init(double af, double cn) {
		bestSolution = new FitAlternatives(1, -1, -1);
		bestSolution.set(0, Double.NaN, Double.NaN);
		this.AF = af;
		this.CN = cn;
		double step = 0.005;
		int alt = (int) Math.floor(((CN + maxCN_Error) - Math.max(0, CN
				- maxCN_Error))
				/ step);
		this.alternativeCN = new double[alt];
		alternativeCN[0] = Math.max(0, CN - maxCN_Error);
		for (int i = 1; i < alternativeCN.length; i++) {
			alternativeCN[i] = alternativeCN[i - 1] + step;
		}
	}

	public void run(int snv_cnv_flag) {
		if (snv_cnv_flag == SNVONLY_FLAG) {
			runSNVonly();
		} else if (snv_cnv_flag == CNVONLY_FLAG) {
			runCNVonly();
		} else if (snv_cnv_flag == SNVANDCNV_FLAG) {
			runSNVandCNV();
		} else if (snv_cnv_flag == SNVBEFORECNV_FLAG) {
			if (Double.isNaN(f_CNV) || f_CNV <= 0) {
				System.err
						.println("The size of the subpopulation harboring the CNV is set to "
								+ f_CNV
								+ ". Invalid value for selected evolutionary scenario. Aborting.");
				System.exit(1);
			}
			runSNVbeforeCNV();
		} else {
			System.err
					.print("Invalid value found for snv_cnv_flag: "
							+ snv_cnv_flag
							+ ". Valid values are: 1 (SNV only), 2 (CNV only) or 3 (SNV and CNV). ");
			System.exit(1);
		}
	}

	/**
	 * Exhaustive iteration through all possible combinations of (PM,PM_B,f),
	 * while considering copy number noise
	 */
	private void runSNVandCNV() {
		int min_PM=1;//Minimal ploidy in mutated cells (PM=0 is not accepted, since B-allele ploidy equation not applicable)
		minErr = Double.POSITIVE_INFINITY;
		int nIter = max_PM - min_PM + 1;
		int[][] PM = new int[nIter][nIter];
		int[][] PM_B = new int[nIter][nIter];
		solutions = new ArrayList<ExPANdS.FitAlternatives>(PM.length
				* PM_B.length);
		for (int pm = min_PM; pm <= max_PM; pm++) {

			for (int pmb = min_PM; pmb <= pm; pmb++) {
				if (pmb / (double) pm <= PN_B / (double) PN) {
					continue;
				}
				
				// Affects the distribution of PM_B: increasing copyPenality shifts PM_B
				// to 0.5*PM - more mutations are heterozygous
//				double copyPenality = Math.pow(Math.abs(pmb-(PN_B+Math.ceil((pm-PN_B)*0.5)))+ 1,3);
				double copyPenality = 1; // No copy penality
				
				RealMatrix coefficients = new Array2DRowRealMatrix(
						new double[][] { { pm - PN }, { pmb - PN_B } }, false);

				// Assign
				int i = pm - min_PM;
				int j = pmb - min_PM;
				PM[i][j] = pm;
				PM_B[i][j] = pmb;

				// FitAlternatives for same PM and PM_B considering error in CN
				// measurement
				FitAlternatives alternatives = new FitAlternatives(
						alternativeCN.length, pm, pmb);
				for (int k = 0; k < alternativeCN.length; k++) {
					double cn = alternativeCN[k];
					RealVector constants = new ArrayRealVector(new double[] {
							cn - PN, cn * AF - PN_B }, false);
					solve(coefficients, constants, copyPenality, cn,
							alternatives, k);
				}

				solutions.add(alternatives);

			}
		}
	}

	/**
	 * Finds alternative solutions to the size of subpopulations harboring an
	 * SNV, under the assumption that a CNV overlapping the mutated locus
	 * happened in a subsequent clonal expansion.
	 */
	private void runSNVbeforeCNV() {
		minErr = Double.POSITIVE_INFINITY;
		solutions = new ArrayList<ExPANdS.FitAlternatives>(PM * 2);

		for (int pmb =1; pmb <= 2; pmb++) {
			double copyPenality = 1; // No copy penality
//			double copyPenality = Math.pow(Math.abs(pmb-(PN_B+Math.ceil((PM-PN_B)*0.5)))+ 1,3);

			// case deletion or copy neutral loss of the mutated allele
			if (pmb / 2.0 > PN_B / (double) PN && PM <= 2) {
				// # of mutated copies lost due to deletion
				for (int pmb_f = 1; pmb_f <= Math.max(1, Math.min(pmb, 2 - PM)); pmb_f++) {
					// pmb*x-pmb_f*f=AF*cn;
					RealMatrix coefficients = new Array2DRowRealMatrix(
							new double[][] { { pmb + pmb_f - PN_B} }, false);
				
					// FitAlternatives for same PMB and PMB_f considering error
					// in CN measurement
					FitAlternatives alternatives = new FitAlternatives(
							alternativeCN.length, -1, pmb);
					for (int k = 0; k < alternativeCN.length; k++) {
						double cn = alternativeCN[k];
						RealVector constants = new ArrayRealVector(
								new double[] { cn * AF + pmb_f * f_CNV - PN_B}, false);
						
						solve(coefficients, constants, copyPenality, cn,
								alternatives, k);
					}
					solutions.add(alternatives);
				}
			}

			// case amplification or copy neutral loss of the WT allele
			if (PM >= 2) {
				for (int pmb_f = pmb + 1; pmb_f <= PM; pmb_f++) {
					if (Math.max(pmb / 2.0,pmb_f/PM) <= PN_B / (double) PN) {
						continue;
					}
					// pmb*x+pmb_f*f=AF*cn;
					RealMatrix coefficients = new Array2DRowRealMatrix(
							new double[][] { { pmb - pmb_f - PN_B} }, false);
					// FitAlternatives for same PMB and PMB_f considering error
					// in CN measurement
					FitAlternatives alternatives = new FitAlternatives(
							alternativeCN.length, -1, pmb);
					for (int k = 0; k < alternativeCN.length; k++) {
						double cn = alternativeCN[k];
						RealVector constants = new ArrayRealVector(
								new double[] { cn * AF - pmb_f * f_CNV - PN_B}, false);
						solve(coefficients, constants, copyPenality, cn,
								alternatives, k);
					}
					solutions.add(alternatives);
				}
			}
		}
	}

	private void solve(RealMatrix coefficients, RealVector constants,
			double copyPenality, double cn, FitAlternatives alternatives, int k) {
		DecompositionSolver solver = new SingularValueDecomposition(
				coefficients).getSolver();
		RealVector solution = solver.solve(constants);
		RealVector errorterm = coefficients.operate(solution).subtract(
				constants);

		double sumErr = Math.abs(CN - cn);
		for (int i = 0; i < errorterm.getDimension(); i++) {
			sumErr += Math.abs(errorterm.getEntry(i));
		}
		sumErr = sumErr * copyPenality;

		alternatives.set(k, solution.getEntry(0), sumErr);
		if (solution.getEntry(0) > -0.1 && solution.getEntry(0) <= 1.1
				&& minErr > sumErr) {
			minErr = sumErr;
			this.bestSolution = alternatives;
			// System.out.println(solution.getEntry(0)+"\t"+errorterm.getEntry(0)
			// +
			// "\t"+errorterm.getEntry(1));
		}
	}

	private void runSNVonly() {
		int min_PM=1;//Minimal ploidy in mutated cells (PM=0 is not accepted, since B-allele ploidy equation not applicable)
		minErr = Double.POSITIVE_INFINITY;
		int nIter = max_PM - min_PM + 1;
		int[] PM_B = new int[nIter];
		solutions = new ArrayList<ExPANdS.FitAlternatives>(PM_B.length);

		for (int pmb = min_PM; pmb <= max_PM; pmb++) {
//			double copyPenality = Math.pow(Math.abs(pmb-(PN_B+Math.ceil((CN-PN_B)*0.5)))+ 1,3);//a-priori maximum likelihood: when half of the existing loci are mutated
			double copyPenality = 1; // No copy penality

			RealMatrix coefficients = new Array2DRowRealMatrix(
					new double[][] { { pmb - PN_B } }, false);
			// Assign
			int j = pmb - min_PM;
			PM_B[j] = pmb;

			// FitAlternatives for same PM and PM_B considering error in CN
			// measurement
			FitAlternatives alternatives = new FitAlternatives(
					alternativeCN.length, -1, pmb);
			for (int k = 0; k < alternativeCN.length; k++) {
				double cn = alternativeCN[k];
				RealVector constants = new ArrayRealVector(new double[] { cn
						* AF - PN_B }, false);
				solve(coefficients, constants, copyPenality, cn, alternatives,
						k);
			}
			solutions.add(alternatives);
		}
	}

	private void runCNVonly() {
		int min_PM=0;//Minimal ploidy in mutated cells (PM=0 accepted, since B-allele is not in same SP)
		minErr = Double.POSITIVE_INFINITY;
		int nIter = max_PM - min_PM + 1;
		int[] PM = new int[nIter];
		solutions = new ArrayList<ExPANdS.FitAlternatives>(PM.length);
		for (int pm = min_PM; pm <= max_PM; pm++) {

			// double copyPenality = Math
			// .pow(Math.abs(pm - Math.round(CN)) + 1, 3);
			double copyPenality = 1; // No copy penality

			RealMatrix coefficients = new Array2DRowRealMatrix(
					new double[][] { { pm - PN } }, false);

			// Assign
			int i = pm - min_PM;
			PM[i] = pm;

			// FitAlternatives for same PM and PM_B considering error in CN
			// measurement
			FitAlternatives alternatives = new FitAlternatives(
					alternativeCN.length, pm, -1);
			for (int k = 0; k < alternativeCN.length; k++) {
				double cn = alternativeCN[k];
				RealVector constants = new ArrayRealVector(new double[] { cn
						- PN }, false);
				solve(coefficients, constants, copyPenality, cn, alternatives,
						k);
			}

			solutions.add(alternatives);

		}
	}

	/**
	 * B-allele ploidy of this locus, obtained by minimizing the deviation
	 * between observed and expected B-allele frequency and copy number.
	 * 
	 * @return ploidy of the B-allele
	 */
	public int getPM_B() {
		return bestSolution.getPM_B();
	}

	/**
	 * Total ploidy of this locus, obtained by minimizing the deviation between
	 * observed and expected B-allele frequency and copy number.
	 * 
	 * @return total ploidy
	 */
	public int getPM() {
		return bestSolution.getPM();
	}

	/**
	 * Frequency of cells harboring the mutation at this locus, obtained by
	 * minimizing the deviation between observed and expected B-allele frequency
	 * and copy number.
	 * 
	 * @return cellular frequency
	 */
	public double getF() {
		try {
			return bestSolution.getBest().getF();
		} catch (Exception e) {
			return Double.NaN;
		}
	}

	public void printFitsToFile(File f) throws IOException {
		BufferedWriter w = Common.getWriter(f.getAbsolutePath());
		w.write(Common.toString(SOLUTION_ENTITIES, "\t"));
		w.newLine();
		for (FitAlternatives fits : this.solutions) {
			Iterator<Solution> iter = fits.iterator();
			while (iter.hasNext()) {

				w.write(iter.next().toString());
				w.newLine();

			}
		}
		w.flush();
		w.close();
	}

	/**
	 * Alternative solutions (PM_B, PM, f)=(B-allele ploidy, total ploidy,
	 * cellular frequency) computed to explain NGS derived allele frequency and
	 * copy number measurements.
	 * 
	 * @return List of alternative solutions
	 */
	public Collection<FitAlternatives> solutions() {
		return solutions;
	}

	public static void main(String[] args) {
		ExPANdS expandsa = new ExPANdS(0.39, 1.33, 0.66, 1,0);
		expandsa.runSNVbeforeCNV();
		
		File[] files = (new File(System.getProperty("user.dir"))
				.listFiles(new FilenameFilter() {

					@Override
					public boolean accept(File arg0, String arg1) {
						return arg1.endsWith(".snv");
					}
				}));
		Arrays.sort(files);
		for (File f : files) {
			System.out.println(f.getName());

			File pair = new File(f.getAbsolutePath() + ".expands");
			if (pair.exists()) {
				System.out.println("ExPANdS file for " + f.getName()
						+ " already exists. Skipped.");
				continue;
			}

			int count = 0;
			try {
				BufferedReader r = Common.getReader(f.getAbsolutePath());
				BufferedWriter w = Common.getWriter(pair.getAbsolutePath());
				String[] header = r.readLine().trim().split("\\s+");
				w.write(Common.toString(header, "\t"));
				w.newLine();
				int countI = Common.firstIndexOf("Count", header);
				int afI = Common.firstIndexOf("AF_Tumor", header);
				int cnI = Common.firstIndexOf("CN_Estimate", header);
				int pnbI = Common.firstIndexOf("PN_B", header);
				int pABBBI = Common.firstIndexOf("pABBB", header);
				int pAABBI = Common.firstIndexOf("pAABB", header);
				int pAAABI = Common.firstIndexOf("pAAAB", header);
				int fI = Common.firstIndexOf("f", header);
				int pmI = Common.firstIndexOf("PM", header);
				int pmbI = Common.firstIndexOf("PM_B", header);
				int devI = Common.firstIndexOf("dev", header);
				for (String l = r.readLine(); l != null; l = r.readLine()) {
					String[] features = l.trim().split("\t");
					if (features.length > header.length) {
						features = Arrays.copyOfRange(features, 1,
								features.length);
					}
					int id = (int) Double.parseDouble(features[countI]);
					int pnb = (int) Math.round(Double
							.parseDouble(features[pABBBI]));
					try {
						ExPANdS expands = new ExPANdS(
								Double.parseDouble(features[afI]),
								Double.parseDouble(features[cnI]), pnb,
								DEFAULT_MAX_PM);
						expands.run(3);
						double pAABB = Double.parseDouble(features[pAABBI]);
						double pABBB = Double.parseDouble(features[pABBBI]);
						double pAAAB = Double.parseDouble(features[pAAABI]);
						if (pABBB >= 0.5 || pAABB + pAAAB >= 0.9) {
							expands.printFitsToFile(new File(f
									.getAbsolutePath() + "." + id + ".fit"));
							count++;
						}
						features[pmI] = "" + expands.getPM();
						features[pmbI] = "" + expands.getPM_B();
						features[fI] = "" + expands.getF();
						features[devI] = "" + expands.getDeviation();
						features[pnbI] = "" + pnb;
					} catch (Exception e) {
						e.printStackTrace();
					}

					w.write(Common.toString(features, "\t").replace("Infinity",
							"Inf"));
					w.newLine();

				}
				w.flush();
				w.close();
			} catch (FileNotFoundException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		}
	}

	/**
	 * Stores alternative cellular frequencies f, that fit the equations:
	 * PM*f+PN*(1-f)=CN and PM_B*f+PN_B*(1-f)=CN*AF, where CN and AF are copy
	 * number and B-allele frequency measurements derived from NGS, while (PM,
	 * PN, PM_B, PN_B)=(total ploidy in mutated cells, total ploidy in normal
	 * cells, B-allele ploidy in mutated cells, B-allele ploidy in normal cells)
	 * are fixed.
	 * 
	 * @author Noemi
	 * 
	 */
	public class FitAlternatives implements Iterator<Solution>,
			Iterable<Solution> {

		/**
		 * Alternative cellular frequencies
		 */
		private double[] f;
		/**
		 * Error associated with each solution.
		 */
		private double[] dev;

		/**
		 * Total ploidy
		 */
		private int pm;

		/**
		 * B-allele ploidy
		 */
		private int pmb;

		/**
		 * Pointer
		 */
		int idx = -2;

		/**
		 * Constructs exhaustive iterator for
		 * 
		 * @param numberAlternatives
		 *            alternative solutions of cellular frequencies for fixed
		 * @param pm
		 *            - total ploidy and
		 * @param pmb
		 *            - B-allele ploidy,
		 */
		private FitAlternatives(int numberAlternatives, int pm, int pmb) {
			f = new double[numberAlternatives];
			dev = new double[numberAlternatives];
			this.pm = pm;
			this.pmb = pmb;
		}

		public int getPM_B() {
			return pmb;
		}

		public int getPM() {
			return pm;
		}

		public Iterator<Solution> iterator() {
			idx = -1;
			return this;
		}

		/**
		 * 
		 * @param k
		 * @return the cellular frequency associated with the k-th solution
		 */
		public double getF(int k) {
			return f[k];
		}

		/**
		 * 
		 * @param k
		 * @return the error term associated with the k-th solution
		 */
		public double getDev(int k) {
			return dev[k];
		}

		/**
		 * 
		 * @return the number of alternative solutions
		 */
		public int size() {
			return dev.length;
		}

		/**
		 * 
		 * @return the solution with the smallest error
		 */
		public Solution getBest() {
			int bestIdx = Common.argmin(dev);
			return new Solution(PN, PN_B, pm, pmb, CN, AF, f[bestIdx],
					dev[bestIdx]);
		}

		/**
		 * Sets a cellular frequency as the i-th solution
		 * 
		 * @param i
		 *            - index of the solution
		 * @param f
		 *            - cellular frequency
		 * @param dev
		 *            - deviation associated with this solution
		 */
		void set(int i, double f, double dev) {
			this.f[i] = f;
			this.dev[i] = dev;
		}

		@Override
		public boolean hasNext() {
			return idx > -2 && idx < f.length - 1;
		}

		@Override
		public Solution next() {
			idx++;
			if (idx >= dev.length) {
				return null;
			}
			return new Solution(PN, PN_B, pm, pmb, CN, AF, f[idx], dev[idx]);
		}

		@Override
		public void remove() {
			// TODO Auto-generated method stub

		}

	}

	/**
	 * To be used as column headers - components constituting the inputs and
	 * outputs of a solution
	 */
	public static final String[] SOLUTION_ENTITIES = { "PN", "PN_B", "PM",
			"PM_B", "CN_Estimate", "AF_Tumor", "f", "dev" };

	/**
	 * Stores alternative solutions (PM_B, PM, f)=(B-allele ploidy, total
	 * ploidy, cellular frequency) computed for the NGS derived measures: AF -
	 * B-allele frequency, CN - copy number and PN_B - ploidy of the B-allele in
	 * normal cells, that are associated with a mutated locus.
	 * 
	 * @author Noemi
	 * 
	 */
	public final class Solution {

		/**
		 * Cellular frequency (output)
		 */
		private double f;
		/**
		 * B-allele frequency among all cells (input)
		 */
		private double AF;
		/**
		 * Average copy number among all cells (input)
		 */
		private double CN;

		/**
		 * B-allele ploidy in mutated cells (output)
		 */
		private int PM_B;

		/**
		 * Total B-allele ploidy in mutated cells (output)
		 */
		private int PM;

		/**
		 * B-allele ploidy in normal cells (input)
		 */
		private int PN_B;

		/**
		 * Total B-allele ploidy in normal cells (input, should equal 2)
		 */
		private int PN;

		/**
		 * Error term assessing how well outputs can explain inputs
		 */
		private double dev;

		public Solution(int PN, int PN_B, int PM, int PM_B, double CN,
				double AF, double f, double dev) {
			this.PN = PN;
			this.PN_B = PN_B;
			this.PM = PM;
			this.PM_B = PM_B;
			this.CN = CN;
			this.AF = AF;
			this.f = f;
			this.dev = dev;

		}

		public double getF() {
			return f;
		}

		public double[] toDouble() {
			return new double[] { PN, PN_B, PM, PM_B, CN, AF, f, dev };
		}

		public String toString() {
			return Common.toString(toDouble(), "\t");
		}
	}
}
