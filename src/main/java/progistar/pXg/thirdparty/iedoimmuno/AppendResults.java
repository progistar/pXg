package progistar.pXg.thirdparty.iedoimmuno;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

public class AppendResults {

	public static void main(String[] args) throws IOException {
		File pXgFile = new File("/Users/gistar/projects/pXg/immunogenicity/ids.tsv");
		File predictionFile = new File("/Users/gistar/projects/pXg/immunogenicity/three_tools_prediction.tsv");
		
		// key = peptide_hla
		Hashtable<String, String> iedbTool = new Hashtable<String, String>();
		Hashtable<String, String> deepHLApan = new Hashtable<String, String>();
		Hashtable<String, String> deepNetBim = new Hashtable<String, String>();
		
		BufferedReader BR = new BufferedReader(new FileReader(predictionFile));
		String line = null;
		BR.readLine();
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[3];
			String hla = fields[2];
			String immunogenicity = fields[4];
			String key = peptide+"_"+hla;
			if(fields[0].equalsIgnoreCase("DeepHLApan")) {
				deepHLApan.put(key, immunogenicity);
			} else if(fields[0].equalsIgnoreCase("DeepNetBim")) {
				deepNetBim.put(key, immunogenicity);
			} else {
				iedbTool.put(key, immunogenicity);
			}
		}
		
		BR.close();
		
		
		BR = new BufferedReader(new FileReader(pXgFile));
		BufferedWriter BW = new BufferedWriter(new FileWriter(pXgFile.getAbsolutePath().replace(".tsv", ".immunogenicity.tsv")));
		String header = BR.readLine();
		header += "\tImmunogenicity3.0\tDeepNetBim\tDeepHLApan";
		BW.append(header);
		BW.newLine();
		
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[24];
			String hla = fields[50];
			String key = peptide+"_"+hla;
			String iedb = iedbTool.get(key);
			String dnb = deepNetBim.get(key);
			String dhp = deepHLApan.get(key);
			
			if(iedb == null) {
				iedb = "";
			}
			if(dnb == null) {
				dnb = "";
			}
			if(dhp == null) {
				dhp = "";
			}
			
			BW.append(line).append("\t").append(iedb).append("\t").append(dnb).append("\t").append(dhp);
			BW.newLine();
		}
		
		BW.close();
		BR.close();
	}
}
