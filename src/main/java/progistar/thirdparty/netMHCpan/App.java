package progistar.thirdparty.netMHCpan;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class App {

	public static void main(String[] args) {
		int tPeptideIndex = 19;
		
		String pXgOutputFileName = "C:\\Users\\progi\\Desktop\\Projects\\pXg\\Qi_MCP2021\\Results\\h1975wc.pxg";
		NetMHCpanResult netMHCpanResult = NetMHCpanParser.parseNetMHCpan("C:\\Users\\progi\\Desktop\\Projects\\pXg\\Qi_MCP2021\\NetMHCpan\\H1975_WC.txt");
		
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
