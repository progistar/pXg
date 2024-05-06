package progistar.thirdparty.netMHCpan;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

public class NetMHCpanRatio {

	public static void main(String[] args) throws IOException {
		File file = new File("C:\\Users\\progi\\Desktop\\Projects\\pXg\\Laumont_NatCommun2016\\Results\\3.S3_pXg\\PeptideAnnotationS3_5ppm_002.rep1.rank10.pXg.netMHCpan_4.1");

		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;

		BR.readLine(); // skip header
		Hashtable<String, String> isDuplicated = new Hashtable<>();
		ArrayList<String> records = new ArrayList<>();

		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String key = fields[0]+"_"+fields[1]+"_"+fields[4];

			if(isDuplicated.get(key) == null) {
				records.add(line);
				isDuplicated.put(key, "");
			}
		}

		BR.close();

		Hashtable<String, Integer> countIDs = new Hashtable<>();
		Hashtable<String, Integer> countMAPs = new Hashtable<>();
		for (String record : records) {
			String[] fields = record.split("\t");
			String isCanonical = fields[36];
			String class_ = fields[43];
			String score = fields[6];
			if(isCanonical.equalsIgnoreCase("TRUE")) {
				Integer count = countIDs.get(score);
				if(count == null) {
					count = 0;
				}

				count++;
				countIDs.put(score, count);

				if(class_.equalsIgnoreCase("NB")) {

				} else {
					count = countMAPs.get(score);
					if(count == null) {
						count = 0;
					}
					count++;
					countMAPs.put(score, count);
				}
			}
		}

		for(int i=99; i>0; i--) {
			String score = (i+"");
			Integer countID = countIDs.get(score);
			Integer countMAP = countMAPs.get(score);

			if(countID == null) {
				countID = 0;
			}
			if(countMAP == null) {
				countMAP = 0;
			}

			System.out.println(score+"\t"+countID+"\t"+countMAP);
		}
	}
}
