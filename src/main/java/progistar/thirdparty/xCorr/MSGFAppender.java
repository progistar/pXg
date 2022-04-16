package progistar.thirdparty.xCorr;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

public class MSGFAppender {

	public static final int MSGF_PEPTIDE_IDX = 9;
	public static final int MSGF_TITLE_IDX = 3;
	
	public static void main(String[] args) throws IOException {
		File xCorrTSV = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/MSGF/unmodified/xCorr/msgfXCorrOutput.tsv");
		File[] msgfIDList = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/MSGF/unmodified/fdr_q5").listFiles();
		
		BufferedReader BR_ = new BufferedReader(new FileReader(xCorrTSV));
		String line = null;
		
		BR_.readLine();// skip header
		Hashtable<String, String> titleToXcorr = new Hashtable<String, String>();
		while((line = BR_.readLine()) != null) {
			String[] fields = line.split("\t");
			String title = fields[0];
			String xCorr = fields[5];
			
			titleToXcorr.put(title, xCorr);
		}
		
		BR_.close();
		
		
		Hashtable<String, String> msgfPSMs = new Hashtable<String, String>();
		// aggregate PSMs
		for(File file : msgfIDList) {
			if(file.getName().startsWith(".")) {
				continue;
			}
			
			if(file.getName().endsWith(".fdr.BA")) {
				BufferedReader BR = new BufferedReader(new FileReader(file));
				BufferedWriter BW = new BufferedWriter(new FileWriter(file.getAbsolutePath()+".xCorr"));
				String header = BR.readLine(); // skip header
				header += "\txCorr";
				
				BW.append(header);
				BW.newLine();
				
				while((line = BR.readLine()) != null) {
					String[] fields = line.split("\t");
					String key = fields[MSGF_TITLE_IDX];
					String peptide = fields[MSGF_PEPTIDE_IDX];
					if(msgfPSMs.get(key) == null) {
						peptide = peptide.replaceAll("[+-.\\(\\)0123456789*]", "");
						msgfPSMs.put(key, peptide);
						
						String xCorr = titleToXcorr.get(key);
						BW.append(line+"\t"+xCorr);
						BW.newLine();
					}
				}
				
				BR.close();
				BW.close();
			}
		}
	}
}
