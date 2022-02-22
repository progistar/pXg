package progistar.thirdparty.netMHCpan;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class App {

	public static void main(String[] args) {
		int tPeptideIndex = 20;
		
		String pXgOutputFileName = "C:\\Users\\progi\\Desktop\\Projects\\pXg\\Laumont_NatCommun2016\\Results\\2.withCalibration\\2.R_pXg\\PeptideAnnotationS2_5ppm_002_recal.rep1.rank10.pXg";
		NetMHCpanResult netMHCpanResult = NetMHCpanParser.parseNetMHCpan("C:\\Users\\progi\\Desktop\\Projects\\pXg\\Laumont_NatCommun2016\\Results\\2.withCalibration\\2.R_pXg\\PeptideAnnotationS2_5ppm_002_recal.rep1.rank10.netMHCpan_4.1");
		
		try {
			BufferedReader BR = new BufferedReader(new FileReader(pXgOutputFileName));
			BufferedWriter BW = new BufferedWriter(new FileWriter(pXgOutputFileName+".netMHCpan"));
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
