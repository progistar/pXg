package progistar.thirdparty.netMHCpan;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class App_pXg {

	public static void main(String[] args) {
		int tPeptideIndex = 20;
		
		String pXgOutputFileName = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/7.Unmodified/pXg/S1.PEAKS.UNMOD.rank10.pXg";
		NetMHCpanResult netMHCpanResult = NetMHCpanParser.parseNetMHCpan("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/7.Unmodified/NetMHCpan/S1.NetMHCpan.xls");
		
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
