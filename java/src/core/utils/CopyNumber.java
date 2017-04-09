package core.utils;


/**
 * @(#) CopyNumber.java
 */

public class CopyNumber extends Measurement {

	public CopyNumber(double value) {
		super(value);
		if (value < 0) {
			throw new IllegalArgumentException(
					"An absolute copy number estimate must be above 0.");
		}
	}

}
