package progistar.thirdparty.examples;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class GTFSelection {

	
	public static void main(String[] args) throws IOException {
		String gtfFileName = "/Users/gistar/resources/Models/gencode.v38.chr_patch_hapl_scaff.annotation.gtf";
		
		String[] geneNames = {"MBNL1", "LAPTM5", "TNFRSF13C"};
		
		
		BufferedReader BR = new BufferedReader(new FileReader(gtfFileName));
		BufferedWriter BW = new BufferedWriter(new FileWriter(gtfFileName+".examples"));
		String line = null;
		
		while((line = BR.readLine()) != null) {
			if(line.startsWith("#")) {
				BW.append(line);
				BW.newLine();
			} else {
				String[] fields = line.split("\t");
				// pass gene
				if(fields[2].equalsIgnoreCase("gene")) {
					continue;
				}
				
				String attrs = fields[8];
				boolean isSelected = false;
				for(String geneName : geneNames) {
					if(attrs.contains(geneName)) {
						isSelected = true;
					}
				}
				
				if(isSelected) {
					BW.append(line);
					BW.newLine();
				}
			}
			
		}
		
		BR.close();
		BW.close();
	}
}
