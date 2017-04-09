package core.utils;


public class AFtimesCopyNumber extends Measurement {

	public AFtimesCopyNumber(double value) {
		super(value);
		if (value < 0) {
			throw new IllegalArgumentException(
					"Product of allele frequency and copy number estimate must be above 0.");
		}
	}

}
