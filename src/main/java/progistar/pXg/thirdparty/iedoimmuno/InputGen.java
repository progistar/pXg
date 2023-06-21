package progistar.pXg.thirdparty.iedoimmuno;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;

public class InputGen {

	public static void main(String[] args) throws IOException {
		StringBuilder batchScript = new StringBuilder();
		
		File targetFile = new File("/Users/gistar/tools/iedb_immunogenecity3.0/10samples/canonical.txt");
		
		BufferedReader BR = new BufferedReader(new FileReader(targetFile));
		BufferedWriter BW = new BufferedWriter(new FileWriter(targetFile.getAbsolutePath().replace(".txt", ".sh")));
		
		String line = null;
		
		Hashtable<String, ArrayList<String>> hlaPepts = new Hashtable<String, ArrayList<String>>();
		BR.readLine(); // skip header
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String hla = fields[0];
			String peptide = fields[1];
			
			ArrayList<String> peptides = hlaPepts.get(hla);
			if(peptides == null) {
				peptides = new ArrayList<String>();
				hlaPepts.put(hla, peptides);
			}
			peptides.add(peptide);
		}
		
		BR.close();
		
		Iterator<String> hlas = (Iterator<String>) hlaPepts.keys();
		int idx = 0;
		while(hlas.hasNext()) {
			String hla = hlas.next();
			idx ++;
			String outputFileName = targetFile.getAbsolutePath().replace(".txt", "."+idx+".txt");
			BufferedWriter hlaBW = new BufferedWriter(new FileWriter(outputFileName));
			
			BW.append("python predict_immunogenicity.py -a ").append(hla).append(" ").append(outputFileName);
			BW.newLine();
			
			for(String pept : hlaPepts.get(hla)) {
				hlaBW.append(pept);
				hlaBW.newLine();
			}
			hlaBW.close();
		}
		
		BW.close();
	}
}
