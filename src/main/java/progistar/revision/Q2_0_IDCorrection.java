package progistar.revision;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

public class Q2_0_IDCorrection {

	public static void main(String[] args) throws IOException {
		File mgfFile = new File("/Users/gistar/projects/pXg/MSMS/B_LCL1.mgf");
		File cometFile = new File("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/4.ReadTranslation/B_LCL1.txt");

		BufferedReader BR = new BufferedReader(new FileReader(mgfFile));
		BufferedWriter BW = new BufferedWriter(new FileWriter(cometFile.getAbsolutePath().replace(".txt", ".title.txt")));
		String line = null;
		Hashtable<String, String> map = new Hashtable<>();
		int msms=0;

		String title = "";
		String scanNr = "";
		double pepMass = 0;
		double charge = 0;
		String rt = "";
		while((line = BR.readLine()) != null) {
			if(line.startsWith("TITLE=")) {
				msms++;
				title = line.split("\\=")[1];
				scanNr = title.split("\\.")[1];
				charge = Double.parseDouble(title.split("\\.")[3]);

			} else if(line.startsWith("RT")) {
				rt = line.split("\\=")[1].split("\\.")[0];
			} else if(line.startsWith("PEPMASS")) {
				pepMass = Double.parseDouble(line.split("\\=")[1].split("\\s")[0]);
				pepMass = pepMass*charge - (charge)*1.0072765;


				String key = scanNr+"_"+pepMass;
				key = key.split("\\.")[0]+"."+key.split("\\.")[1].substring(0, 2)+"_"+rt;


				if(map.get(key) != null) {
					System.out.println("Dup:"+scanNr+"\t"+key);
				}
				map.put(key, title);
			}
		}

		BR.close();

		BR = new BufferedReader(new FileReader(cometFile));
		BW.append(BR.readLine());
		BW.newLine();
		BW.append(BR.readLine());
		BW.newLine();

		int cnt = 0;
		int cCnt = 0;
		while((line = BR.readLine()) != null) {
			cCnt++;
			String[] fields = line.split("\t");
			String thisScanNr = fields[0];
			String mass = fields[3]+"0000";
			rt = fields[18].split("\\.")[0];
			String key = thisScanNr+"_"+mass.split("\\.")[0]+"."+mass.split("\\.")[1].substring(0,2)+"_"+rt;

			title = map.get(key);
			if(title == null) {
				int end = Integer.parseInt(key.split("\\_")[2]);
				key = thisScanNr+"_"+mass.split("\\.")[0]+"."+mass.split("\\.")[1].substring(0,2)+"_"+(end-1);
				title = map.get(key);
				if(title == null) {
					System.out.println(key);
				} else {
					cnt++;
				}
			} else {
				cnt++;
			}

			BW.append(title);
			for(int i=1; i<fields.length; i++) {
				BW.append("\t"+fields[i]);
			}
			BW.newLine();
		}

		System.out.println(cCnt);
		System.out.println(cnt+"/"+msms);
		BR.close();
		BW.close();
	}
}
