package progistar.thirdparty.netMHCpan;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class App_pXg {

	public static void main(String[] args) {
		int tPeptideIndex = 3;
		
		String pXgOutputFileName = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/9.Unmodified_10ppm/pXg_Unmatched/S4.RAW.PEAKS.csv.top1.unided";
		NetMHCpanResult netMHCpanResult = NetMHCpanParser.parseNetMHCpan("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/9.Unmodified_10ppm/pXg_Unmatched/S4_NetMHCpan.xls");
		
		try {
			BufferedReader BR = new BufferedReader(new FileReader(pXgOutputFileName));
			BufferedWriter BW = new BufferedWriter(new FileWriter(pXgOutputFileName+".BA"));
			String line = BR.readLine() +"\t"+netMHCpanResult.getHeader(); // header;
			
			// append header
			BW.append(line.replaceAll("\\,", "\t"));
			BW.newLine();
			
			while((line = BR.readLine()) != null) {
				line = line.replaceAll("\\,", "\t");
				String[] fields = line.split("\t");
				String peptide = fields[tPeptideIndex];
				
				if(peptide.length() >7 && peptide.length() < 16) {
					BW.append(line).append("\t").append(netMHCpanResult.getHLATyping(peptide));
					BW.newLine();
				}
				
			}
			
			BW.close();
			BR.close();
		}catch(IOException ioe) {
			
		}
	}
}
