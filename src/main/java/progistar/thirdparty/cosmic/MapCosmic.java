package progistar.thirdparty.cosmic;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Pattern;

public class MapCosmic {

	public static void main(String[] args) throws IOException {
		Pattern gLociPattern = Pattern.compile("(chr[0-9]+:[0-9]+)");
		
		String[] cosmicFilePath = {"/Users/gistar/projects/pXg/COSMICMutations/CosmicNCV.tsv", "/Users/gistar/projects/pXg/COSMICMutations/CosmicMutantExport.tsv"};
		String[] formats = {"format1","format2"};
		String pXgFilePath = "/Users/gistar/projects/pXg/COSMICMutations/pXgMutations.pXg.log";
		
		Cosmic cosmic = new Cosmic(cosmicFilePath, formats);
		
		BufferedReader BR = new BufferedReader(new FileReader(pXgFilePath));
		String line = null;
		
		
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String[] mutations = fields[24].split("\\|");
			
			for(String mutation : mutations) {
				ArrayList<CosmicMutation> cMutations = cosmic.mutations.get(mutation);
				if(cMutations != null) {
					for(CosmicMutation cMutation : cMutations) {
						System.out.println(mutation+","+cMutation.cosmicID+","+cMutation.somaticInfo);
					}
				}
			}
		}
		
		BR.close();
		
	}
}
