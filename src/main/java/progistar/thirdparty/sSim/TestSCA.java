package progistar.thirdparty.sSim;

import java.text.DecimalFormat;
import java.util.ArrayList;

public class TestSCA {

	public static void main(String[] args) {
		ArrayList<double[]> experimentalPeaks = new ArrayList<>();
		ArrayList<double[]> syntheticPeaks = new ArrayList<>();





		double[] value1 = {147.1122589,1778.018920898400}; double[] value4 = {305.1812744,4569.834960937500};
		double[] value2 = {175.0709076,2376.356445312500}; double[] value5 = {404.2494812,2525.513183593800};
		double[] value3 = {234.1444092,3922.728515625000}; double[] value6 = {532.3082886,636.505737304700};

		double[] value7 = {589.3294678,6332.264648437500};
		double[] value8 = {752.3916626,1908.116577148400};
		double[] value9 = {839.4284058,226.291168212900};

		double[] value1_ = {147.1002502,18800.380859375000}; double[] value4_ = {234.1450348,267222.218750000000};
		double[] value2_ = {147.1128082,330881.875000000000}; double[] value5_ = {305.1817322,287241.968750000000};
		double[] value3_ = {175.0712585,130659.507812500000}; double[] value6_ = {404.2502441,154387.562500000000};

		double[] value7_ = {532.3081055,37581.449218750000};
		double[] value8_ = {589.3311768,371156.750000000000};
		double[] value9_ = {752.3934937,115102.343750000000};
		double[] value10_ = {839.4251709,15495.421875000000};

		experimentalPeaks.add(value1);
		experimentalPeaks.add(value2);
		experimentalPeaks.add(value3);
		experimentalPeaks.add(value4);
		experimentalPeaks.add(value5);
		experimentalPeaks.add(value6);
		experimentalPeaks.add(value7);
		experimentalPeaks.add(value8);
		experimentalPeaks.add(value9);

		syntheticPeaks.add(value1_);
		//syntheticPeaks.add(value2_);
		syntheticPeaks.add(value3_);
		syntheticPeaks.add(value4_);
		syntheticPeaks.add(value5_);
		syntheticPeaks.add(value6_);
		syntheticPeaks.add(value7_);
		syntheticPeaks.add(value8_);
		syntheticPeaks.add(value9_);
		syntheticPeaks.add(value10_);

		double maxThr = 42;
		double maxExp = 9;

		int sizeOfExp = experimentalPeaks.size();
		int sizeOfThr = syntheticPeaks.size();
		double[][] expPeaks = new double[sizeOfExp][2];
		double[][] thrPeaks = new double[sizeOfThr][2];

		System.out.println("Normalized by maxInt");

		System.out.print("{");
		for(int i=0; i<sizeOfThr; i++) {
			thrPeaks[i][0] = syntheticPeaks.get(i)[0];
			thrPeaks[i][1] = syntheticPeaks.get(i)[1];
			thrPeaks[i][1] = thrPeaks[i][1] / maxThr;

			if(i!=0) {
				System.out.print(",");
			}
			System.out.print(thrPeaks[i][1]);
		}
		System.out.println("}");

		for(int i=0; i<sizeOfExp; i++) {
			expPeaks[i][0] = experimentalPeaks.get(i)[0];
			expPeaks[i][1] = experimentalPeaks.get(i)[1];
			expPeaks[i][1] = expPeaks[i][1] / maxExp;
		}



		// normalization theoretical peaks.
		double norThr = 0;
		for(int i=0; i<sizeOfThr; i++) {
			norThr += Math.pow(thrPeaks[i][1],2);
			thrPeaks[i][1] = thrPeaks[i][1];
		}
		norThr = Math.sqrt(norThr);

		System.out.println("Normalized by sum of sqrt");

		System.out.print("{");
		for(int i=0; i<sizeOfThr; i++) {
			thrPeaks[i][1] = thrPeaks[i][1] / norThr;
			if(i!=0) {
				System.out.print(",");
			}
			System.out.print(thrPeaks[i][1]);
		}
		System.out.println("}");

		// normalization experimental peaks.
		double norExp = 0;
		for(int i=0; i<sizeOfExp; i++) {
			norExp += Math.pow(expPeaks[i][1], 2);
			expPeaks[i][1] = expPeaks[i][1];
		}
		norExp = Math.sqrt(norExp);
		for(int i=0; i<sizeOfExp; i++) {
			expPeaks[i][1] = expPeaks[i][1] / norExp;
		}

		//
		double innerproduct = 0;
		// inner product
		int thrIndex = 0;
		int expIndex = 0;
		boolean doMatchFurther = true;
		while(doMatchFurther) {
			double delta = thrPeaks[thrIndex][0] - expPeaks[expIndex][0];
			if(Math.abs(delta) <= 0.02) {
				innerproduct += thrPeaks[thrIndex][1] * expPeaks[expIndex][1];
				expIndex++;
				thrIndex++;
			}
			else if(delta > 0) {
				expIndex++;
			} else if(delta < 0) {
				thrIndex++;
			} else if(delta == 0) {
				expIndex++;
				thrIndex++;
			}

			if(thrIndex >= thrPeaks.length || expIndex >= expPeaks.length) {
				doMatchFurther = false;
			}
		}

		// calculate experimental Norm.
		DecimalFormat decimalFormat = new DecimalFormat("#.################");
		innerproduct = Double.parseDouble(decimalFormat.format(innerproduct));
		double arccos = Double.parseDouble(decimalFormat.format(Math.acos(innerproduct)));
		double midval = Double.parseDouble(decimalFormat.format(arccos / Math.PI));
		double sca = 1 - 2 * midval;
		System.out.println(sca);
	}
}
