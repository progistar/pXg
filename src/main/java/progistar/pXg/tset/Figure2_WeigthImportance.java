package progistar.pXg.tset;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class Figure2_WeigthImportance {

	public static void main(String[] args) throws IOException {
		File[] files = new File("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_features").listFiles();
		String[] featureSet = {"MainScore", "DeltaScore", "Log2Reads", "absppm", "Charge1", "Charge2", "Charge3", "Charge4", "SA",
				"BestDeltaRT", "Log2MeanQScore"};
		for(int i=0; i<featureSet.length; i++) {
			System.out.print(featureSet[i]+"\t");
		}
		System.out.println("Sample");
		for(File file : files) {
			if(!file.getName().endsWith("feat2.weights")) {
				continue;
			}
			
			String sampleName = file.getName().split("\\.")[0];
			
			BufferedReader BR = new BufferedReader(new FileReader(file));
			String line = null;
			while((line = BR.readLine()) != null) {
				if(line.startsWith("#")) continue;
				String[] fields = line.split("\t");
				String[] features = BR.readLine().split("\t");
				
				for(int i=0; i<featureSet.length; i++) {
					boolean hasFeature = false;
					for(int j=0; j<fields.length; j++) {
						if(fields[j].equalsIgnoreCase(featureSet[i])) {
							System.out.print(features[j]+"\t");
							hasFeature = true;
						}
					}
					if(!hasFeature) {
						System.out.print("0\t");
					}
				}
				System.out.println(sampleName);
				BR.readLine();
			}
			BR.close();
		}
	}
}
