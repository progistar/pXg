package progistar.thirdparty.netMHCpan;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class App {

	public static void main(String[] args) {
		int tPeptideIndex = 20;
		
		String pXgOutputFileName = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/5.withCalibrationAddScanNumWithoutDeami/pXg/PeptideAnnotationS4_5ppm_002_recal.scanNum.rep1.rank10.pXg";
		NetMHCpanResult netMHCpanResult = NetMHCpanParser.parseNetMHCpan("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/5.withCalibrationAddScanNumWithoutDeami/pXg/S4.NetMHCpan");
		
		try {
			BufferedReader BR = new BufferedReader(new FileReader(pXgOutputFileName));
			BufferedWriter BW = new BufferedWriter(new FileWriter(pXgOutputFileName+".BA"));
			String line = BR.readLine() +"\t"+netMHCpanResult.getHeader(); // header;
			
			// append header
			BW.append(line);
			BW.newLine();
			
			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				String peptide = fields[tPeptideIndex];
				
				BW.append(line).append("\t").append(netMHCpanResult.getHLATyping(peptide));
				BW.newLine();
			}
			
			BW.close();
			BR.close();
		}catch(IOException ioe) {
			
		}
	}
}
