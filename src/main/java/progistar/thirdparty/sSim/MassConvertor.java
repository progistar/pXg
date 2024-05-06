package progistar.thirdparty.sSim;

public class MassConvertor {

	public static double mzToMass (double mz, double charge) {
		double mass = charge * (mz);
		return mass;
	}

	public static double massToMz (double mass, double charge) {
		double mz = (mass)/charge;
		return mz;
	}
}