package progistar.thirdparty.sSim;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class MergeSpectra {

	
	public static void main(String[] args) throws IOException {
		File[] fileList = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/ProteomeToolsHLA/selectedMGF_pXg").listFiles();
		File outputFile = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/ProteomeToolsHLA/selectedMGF_pXg.mgf");
		
		BufferedWriter BW = new BufferedWriter(new FileWriter(outputFile));
		
		for(File file : fileList) {
			if(file.getName().startsWith(".")) continue;
			if(!file.getName().endsWith(".mgf")) continue;
			
			BufferedReader BR = new BufferedReader(new FileReader(file));
			
			String line = null;
			
			while((line = BR.readLine()) != null) {
				BW.append(line);
				BW.newLine();
			}
			
			BR.close();
		}
		
		BW.close();
	}
}
