package progistar.peaks;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class RatioOfncPSMsWithBA {

	public static void main(String[] args) throws IOException {
		String fileName = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/7.Unmodified/pXg/S2.PEAKS.UNMOD.rank10.pXg.BA.xCorr";
		BufferedReader BR = new BufferedReader(new FileReader(fileName));

		String line = null;
		int alcScoreIdx = 7;
		int mhcBAIdx = 43;
		int isCanonicalIdx = 36;
		
		BR.readLine(); // skip header
		
		double[] accBAs = new double[101];
		double[] accNonBAs = new double[101];
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			
			if(fields[isCanonicalIdx].equalsIgnoreCase("TRUE")) {
				continue;
			}
			String isBA = fields[mhcBAIdx];
			Integer alcScore = Integer.parseInt(fields[alcScoreIdx]);
			if(isBA.equalsIgnoreCase("NB")) {
				accNonBAs[alcScore] += 1;
			} else {
				accBAs[alcScore] += 1;
			}
		}
		
		for(int i=accBAs.length-2; i>=0; i--) {
			accBAs[i] += accBAs[i+1];
			accNonBAs[i] += accNonBAs[i+1];
		}
		
		for(int i=accBAs.length-1; i>=0; i--) {
			if(accBAs[i] != 0 || accNonBAs[i] != 0) {
				System.out.println(i+"\t"+(accBAs[i]/(accBAs[i]+accNonBAs[i])));
			}
		}
		
		BR.close();
	}
}
