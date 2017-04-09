package core.utils;
/**
 * @(#) Measurement.java
 */

public abstract class Measurement {
	private double value;

	public double getValue() {
		return value;
	}

	public Measurement(double value) {
		this.value = value;
	}

}
