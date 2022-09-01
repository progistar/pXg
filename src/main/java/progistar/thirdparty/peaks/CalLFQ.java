package progistar.thirdparty.peaks;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

public class CalLFQ {

	public static void main(String[] args) throws IOException {
		File pXgResult = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/9.Unmodified_10ppm/pXg/S4.RAW.PEAKS.pxg.BA");
		File peaksFeature = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/9.Unmodified_10ppm/PEAKS/S4.RAW.PEAKS.FEAT.csv");
		BufferedWriter BW = new BufferedWriter(new FileWriter(pXgResult.getAbsolutePath()+".LFQ"));
		BufferedReader BR = new BufferedReader(new FileReader(pXgResult));
		String line = null;
		ArrayList<String> pXgRecords = new ArrayList<String>();
		Hashtable<String, ArrayList<String>> peptToFeature = new Hashtable<String, ArrayList<String>>();
		
		String header = BR.readLine() +"\tAbundance";
		while((line = BR.readLine()) != null) {
			pXgRecords.add(line);
			String[] fields = line.split("\t");
			String pept = fields[3];
			String featureID = fields[2];
			
			ArrayList<String> features = peptToFeature.get(pept);
			if(features == null) {
				features = new ArrayList<String>();
				peptToFeature.put(pept, features);
			}
			features.add(featureID);
		}
		
		BR.close();
		
		BR = new BufferedReader(new FileReader(peaksFeature));
		Hashtable<String, Double> lfqTable = new Hashtable<String, Double> ();
		
		BR.readLine(); // skip header
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\\,");
			
			String key = "F"+fields[0]+":"+fields[1];
			Double quant = Double.parseDouble(fields[7]);
			lfqTable.put(key, quant);
		}
		
		BR.close();
		
		BW.append(header);
		BW.newLine();
		
		for(int i=0; i<pXgRecords.size(); i++) {
			BW.append(pXgRecords.get(i));
			String[] fields = pXgRecords.get(i).split("\t");
			String key = fields[3];
			ArrayList<String> features = peptToFeature.get(key);
			ArrayList<Double> lfqValues = new ArrayList<Double>();
			
			for(String feature : features) {
				Double quant = lfqTable.get(feature);
				if(quant != null) {
					lfqValues.add(quant);
				}
			}
			
			if(lfqValues.size() == 0) {
				BW.append("\t0");
			} else {
				Collections.sort(lfqValues);
				int leftIdx = (lfqValues.size()-1)/2;
				int rightIdx = (lfqValues.size())/2;
				double value = (lfqValues.get(leftIdx)+lfqValues.get(rightIdx)) / 2;
				BW.append("\t"+value);
			}
			BW.newLine();
		}
		
		BW.close();
	}
}
